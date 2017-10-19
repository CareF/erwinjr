#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include "complex.h"
#include <omp.h>

#define ANG 1E-10 /*angstrom in m */
#define KVpCM 1E5 /*kV/cm in V/m (SI unit) */
#define sq(X) (X)*(X)

const double hbar = 6.626e-34 / 2 / 3.141592653589793;
const double m0 = 9.109e-31;
const double e0 = 1.602e-19;
const double pi = 3.1415926535897932385;
/*
 * TODO: using only selected MP functions
 * */

void psiFn(double Eq, int startpoint, int xPsiSize, double xres, 
		const double *xVc, const double *xEg, const double *xF, 
		const double *xEp, const double *xESO, const double *xMc, 
		double *xMcE, double *xPsi)
{ /*A ODE solver for wave function, using Euler method.  
	(TB improve: Numerov method)
	The formular is Eq.(2.6) in the Thesis for ODE solver
TODO: try remove Eq
	
INPUT:
	Eq is the eigen value in eV (same for xEg, xEp, xESO)
	startpoint controls the area/startpoint of the wavefunction, 
		so that the result is the states confined in the first well 
		after startpoint (if xPsi[end] is 0)
	xPsiSize is the range for the following arrays
	xres for the resolution (length per pixal, in Angstrom)
	xVc[x] is the electrical potential at x, including:
		- material band offset (for Gamma point and conduction band)
		- external field
		- electron distribution (self-consistency calculation required)
	xEg[x] is the band gap at Gamma point. Correction in Eq.(2.15) should 
		already been included. 
		(following from Eq.(2.20) and Eq.(2.17) in the Thesis)
	xF[x] is the Kane parameter representing the second-order k.p perturbation 
	xEp[x] is is the energy-unit representatiation of the momentum matrix 
		element between the s-like conduction bands and p-like valence bands
	xESO[x] is spin-orbit splitting (noted as \Delta_{SO} in the thesis)
		--- Q_\epsilon is negeleted for approximation
	xMc[x] is effective mass calculate for 2nd order perturbation (no kinetic 
		energy dependence), and is actually not used

RESUTL:
	xPsi[x] is the wave function
	xMcE[x] is the effective mass at position x, 
		and is calculated according to Eq.(2.20)
	*/
	for(int q=0; q<xPsiSize; q++)
	{ 
		if(1) //?
			/*Calculate effective mass by Eq.(2.20)
			 * but subtract xVc for non-material effect */
			xMcE[q] = m0 / ( 1 + 2*xF[q] + xEp[q]/3 * (2 / (Eq - xVc[q] + xEg[q]) 
						+ 1 / (Eq - xVc[q] + xEg[q] + xESO[q]) ));
		else
			xMcE[q] = m0 * xMc[q] * (1 - (xVc[q] - Eq) / xEg[q]);
		if(q>1)
			/*Linear interpolation to move half step, for the recursion of xPsi*/
			xMcE[q-1] = 0.5 * (xMcE[q] + xMcE[q-1]);
	}
	for(int q=0; q<startpoint; q++) 
		xPsi[q] = 0; 
	xPsi[startpoint] = 1;
	//*xPsi=0;
	//*(xPsi+1)=1;
	for(int q=startpoint; q<xPsiSize-1; q++)
	{
		/* (xVc[q] - Eq) is in unit eV, need a e0 factor to transform unit */
		xPsi[q+1] = (
				(2*sq(xres*ANG /hbar) * (xVc[q] - Eq)*e0 
				 + 1 / xMcE[q] + 1 / xMcE[q-1]) * xPsi[q] 
				- xPsi[q-1] / xMcE[q-1]
				) * xMcE[q];
	}
	return;
}

int psiFnEnd(const double *eEq, int eEqSize, int xPsiSize, double xres, 
		double EField, const double *xVc, const double *xEg, const double *xF, 
		const double *xEp, const double *xESO, const double *xMc, 
		double *xMcE, double *xPsi, double *xPsiEnd)
{ /*To get boundary dependence of energy Eq (Fig.2.1 left in the thesis) 
	to help decide the eigenvalue for zero boundary condition 
 
INPUT:
	eEq[n] is the required energy series
	eEqSize is the size of eEq
	xPsiSize, xres, xVc, ..., xPsi are parameters for psiFn
	EField is static external electrical field, in unit kV/cm

OUTPUT:
	xPsiEnd is wavefunction (non-normalized) at the end, 
		supposed to be 0 for eigenenergy 
  */
	const double extLength=200; /*nm, the extend length for start point*/ 
#pragma omp parallel for
	for(int q=0; q<eEqSize; q++)
	{
		double * xMcE_para = (double *) malloc(xPsiSize * sizeof(double));
		double * xPsi_para = (double *) malloc(xPsiSize * sizeof(double));
		double Eq = eEq[q];
		/* set start point, according to energy offset and external field */
		int startpoint = xPsiSize - ceil(
				(Eq - eEq[0])/(EField * KVpCM * ANG * xres) + extLength/xres);
		if(startpoint<1) 
			startpoint = 1;

		psiFn(Eq, startpoint, xPsiSize, xres, 
				xVc, xEg, xF, xEp, xESO, xMc, xMcE_para, xPsi_para);
		xPsiEnd[q] = xPsi_para[xPsiSize-1];
		free(xMcE_para);
		free(xPsi_para);
		//printf("%d: %g %d        ", q, eEq[q], startpoint);
		//printf("%d  ", startpoint);
	}

	return 1;
}

int inv_quadratic_interp(const double *x, const double *y, 
		const int64_t *idxs, int idxLength, double *root)
{ /* Using inverse quadratic interpolation to get x so that y(x)=0, 
	and thus x is the eigen-energy satisfies zero boundary condition

INPUT:
	y[n] is the boundary value for E=x[n]
	x[ idxs[n] ] is the approximate n-th eigen-energy, 
	idxLength is size of idxs

OUTPUT:
	root[n] is the interpolation result for y(x)=0, 
		with length idxLength */
	int n;
	int q;
	for(q=0; q<idxLength; q++)
	{
		n = idxs[q];
		/* printf("No. %d at %d\n",q, n); */
		/* printf("xnear: %f, %f, %f\n",x[n-1]-x[n],x[n],x[n+1]-x[n]); */
		root[q] = x[n-1]*y[n]*y[n+1] / ( (y[n-1]-y[n]  )*(y[n-1]-y[n+1]) ) 
			      + x[n]*y[n-1]*y[n+1] / ( (y[n]  -y[n-1])*(y[n]  -y[n+1]) )
			      + x[n+1]*y[n-1]*y[n] / ( (y[n+1]-y[n-1])*(y[n+1]-y[n]  ) );
	}
	return 1;
}


int returnme()
{return 42;}

int psiFill(int xPsiSize, double xres, int EigenESize, double *EigenE, 
		double *xVc, double *xEg, double *xF, double *xEp, double *xESO, 
		double *xMc, double *xMcE, double *xyPsi)
{ /* To calculate a series of wave function according to given eigen energy
INPUT:
	EigenE[n] is the n-th eigen-energy, with length EigenESize
	others see psiFn
OUTPUT:
	xyPsi[n] is the wave function corresponding to EigenE[n] */
	for(int col=0; col<EigenESize; col++) // loop on column
	{
		double Eq = EigenE[col];
		psiFn(Eq, 1, xPsiSize, xres, xVc, xEg, xF, xEp, xESO, xMc, xMcE, 
				xyPsi+col*xPsiSize);
		/* Normalization */
		/*double PsiInt = 1; //one for xPsi[1]
		for(int q=1; q<xPsiSize-1; q++)
		{ 
			PsiInt += sq(xyPsi[q+col*xPsiSize]) * ( 
					1 + ( Eq-xVc[q] ) / ( Eq-xVc[q]+xEg[q] ) );
		}*/
		/* TODO: positive charge holes? */
		double PsiInt = 0; 
		for(int q=0; q<xPsiSize; q++) 
			/*Usual normalization*/
			/* PsiInt += sq(xyPsi[col*xPsiSize + q]);  */
			/* Eq.(2.55) in the thesis*/
			PsiInt += sq(xyPsi[col*xPsiSize + q]) * (
						1 + ( Eq - xVc[q] ) / ( Eq - xVc[q] + xEg[q] ) );
		double NormFactor = 1 / sqrt(xres * ANG * PsiInt);

		for(int q=0; q<xPsiSize; q++)
			xyPsi[col*xPsiSize+q] *= NormFactor;
	}
	return 1;
}


typedef struct {complex aa; complex ab; complex ba; complex bb;} matrix;

matrix identity()
{
    matrix c;
    c.aa = cmplx(1,0);
    c.ab = cmplx(0,0);
    c.ba = cmplx(0,0);
    c.bb = cmplx(1,0);
    return c;
}

matrix mmult(matrix m1, matrix m2)
{
    matrix c;
    c.aa = cxadd(cxmul(m1.aa,m2.aa) , cxmul(m1.ab,m2.ba));
    c.ab = cxadd(cxmul(m1.aa,m2.ab) , cxmul(m1.ab,m2.bb));
    c.ba = cxadd(cxmul(m1.ba,m2.aa) , cxmul(m1.bb,m2.ba));
    c.bb = cxadd(cxmul(m1.ba,m2.ab) , cxmul(m1.bb,m2.bb));
    return c;
}

#define MAXLENGTH 30
void chiImag_array(double wavelength, const double *thicknesses, 
		const double *indexesReal, const double *indexesImag, int numLayers, 
		const double *betaInReal, const double *betaInImag, 
		int numBetas, double *chiImag)
{

	double k = 2 * pi / wavelength;
	double z0 = 0.003768;

	complex index[MAXLENGTH];

	for (int j=0; j < numLayers; j++)
		index[j] = cmplx(indexesReal[j],indexesImag[j]);

	int q=0;
	for (q=0; q<numBetas; q++)
	{
		int j = 0;

		complex beta = cmplx(betaInReal[q], betaInImag[q]);

		//alpha = sqrt(self.stratumRIndexes**2-beta**2)
		complex alpha[MAXLENGTH];
		for (j=0; j < numLayers; j++)
			alpha[j] = cxsqrt(cxsub(cxsqr(index[j]),cxsqr(beta)));
		//make sure correct sign of alphac and alphas are chosen, see Chilwell
		//TODO: why only alpha0 and alpha[n-1]
		if(imag(alpha[0]) < 0)
			alpha[1] = cxneg(alpha[0]); // Why neg? in ThePhysics.py it's conj
		if(imag(alpha[numLayers-1]) < 0)
			alpha[numLayers-1] = cxneg(alpha[numLayers-1]);

		//gamma = z0*alpha/self.stratumRIndexes**2
		complex gamma[MAXLENGTH];
		for (j=0; j < numLayers; j++)
			gamma[j] = cmuld(cxdiv(alpha[j],cxsqr(index[j])),z0);

		//phi = k*self.stratumThicknesses*alpha
		complex phi[MAXLENGTH];
		for (j=0; j < numLayers; j++)
			phi[j] = cmuld(alpha[j],k*thicknesses[j]);

		double zeta[MAXLENGTH];
		for (j=0; j < numLayers; j++)
			zeta[j] = k*thicknesses[j]/z0;

		complex chi;
		matrix m=identity(), mt=identity();
		complex ni=cmplx(0,-1);

		for (j=numLayers-1; j>-1; j--)
			//array([[cos(phi[q]),              -1j/gamma[q]*sin(phi[q])],
			//       [-1j*gamma[q]*sin(phi[q]), cos(phi[q])]])
		{
			if (real(index[j]) == real(beta) && imag(index[j]) == imag(beta))
				mt.ab = cxmul(cmuld(ni,zeta[j]),cxsqr(index[j]));  // -i*k*thickness*n^2/z0
			else
				mt.ab = cxdiv(cxmul(ni,cxsin(phi[j])),gamma[j]);  // -i*sin(phi)/gamma

			mt.aa = cxcos(phi[j]);
			mt.ba = cxmul(cxmul(ni,cxsin(phi[j])),gamma[j]);   // -i*sin(phi)*gamma
			mt.bb = mt.aa;

			m = mmult(mt,m); 
		}
		chi = cxadd4(cxmul(gamma[numLayers-1],m.aa) , 
				cxmul(gamma[0],cxmul(gamma[numLayers-1],m.ab)) , 
				m.ba , 
				cxmul(gamma[0],m.bb) );
		//chi = gammac*M[0,0] + gammac*gammas*M[0,1] + M[1,0] + gammas*M[1,1]

		chiImag[q] = imag(chi);
	}
}

double abschi_find(double wavelength, const double *thicknesses, const double *indexesReal, 
                const double *indexesImag, int numLayers, double betaInReal, double betaInImag)
{

    double pi = 3.1415926535897932385;
    double k = 2 * pi / wavelength;
    double z0 = 0.003768;
    int j = 0;

    complex beta = cmplx(betaInReal, betaInImag);
    complex index[MAXLENGTH];
    for (j=0; j < numLayers; j++)
        index[j] = cmplx(indexesReal[j],indexesImag[j]);
    
    //alpha = sqrt(self.stratumRIndexes**2-beta**2)
    complex alpha[MAXLENGTH];
    for (j=0; j < numLayers; j++)
        alpha[j] = cxsqrt(cxsub(cxsqr(index[j]),cxsqr(beta)));
    //make sure correct sign of alphac and alphas are chosen, see Chilwell
    if(imag(alpha[0]) < 0)
        alpha[0] = cxneg(alpha[0]);
    if(imag(alpha[numLayers-1]) < 0)
        alpha[numLayers-1] = cxneg(alpha[numLayers-1]);

    //gamma = z0*alpha/self.stratumRIndexes**2
    complex gamma[MAXLENGTH];
    for (j=0; j < numLayers; j++)
        gamma[j] = cmuld(cxdiv(alpha[j],cxsqr(index[j])),z0);
    
    //phi = k*self.stratumThicknesses*alpha
    complex phi[MAXLENGTH];
    for (j=0; j < numLayers; j++)
        phi[j] = cmuld(alpha[j],k*thicknesses[j]);

    double zeta[MAXLENGTH];
    for (j=0; j < numLayers; j++)
        zeta[j] = k*thicknesses[j]/z0;
    
    complex chi;
    matrix m=identity(), mt=identity();
    complex ni=cmplx(0,-1);

    for (j=numLayers-1; j>-1; j--)
    //array([[cos(phi[q]), -1j/gamma[q]*sin(phi[q])],[-1j*gamma[q]*sin(phi[q]), cos(phi[q])]])
    {
        if (real(index[j]) == real(beta) && imag(index[j]) == imag(beta))
            mt.ab = cxmul(cmuld(ni,zeta[j]),cxsqr(index[j]));  // -i*k*thickness*n^2/z0
        else
            mt.ab = cxdiv(cxmul(ni,cxsin(phi[j])),gamma[j]);  // -i*sin(phi)/gamma

        mt.aa = cxcos(phi[j]);
        mt.ba = cxmul(cxmul(ni,cxsin(phi[j])),gamma[j]);   // -i*sin(phi)*gamma
        mt.bb = mt.aa;
        
        m = mmult(mt,m);
    }
    chi = cxadd4(cxmul(gamma[numLayers-1],m.aa) , 
                       cxmul(gamma[0],cxmul(gamma[numLayers-1],m.ab)) , 
                       m.ba , 
                       cxmul(gamma[0],m.bb) );
    //chi = gammac*M[0,0] + gammac*gammas*M[0,1] + M[1,0] + gammas*M[1,1]
    return cxabs(chi);  // return a double value
}

int argmin(const double *values, int numElements)
{
    int idx=0;
    double minValue = values[0];
    int j=0;
    for (j=1; j<numElements; j++)
        idx = values[j] < minValue ? j : idx;
    return idx;
}

void beta_find(double wavelength, const double *thicknesses, const double *indexesReal, 
               const double *indexesImag, int numLayers, double betaInReal, double betaInImag,
               double beta_find_precision, double *betaOut)
{
    //initialize beta0
    complex beta0 = cmplx(betaInReal, betaInImag);

    complex rInc = cmplx(0.0001,0);
    complex iInc = cmplx(0,1.0e-6);

    int numBetas = 9;
    complex betas[9] = {0};

    double abschiOld=0.0, abschiNew=0.0, absdeltachi=0.0;
    int chiMinIdx = 0;
    int q=0;
    do{
        //set betas
        betas[0] = beta0;
        betas[1] = cxadd(beta0,rInc); 
        betas[2] = cxsub(beta0,rInc);
        betas[3] = cxadd(beta0,iInc); 
        betas[4] = cxsub(beta0,iInc);
        betas[5] = cxadd(cxadd(beta0,rInc),iInc);
        betas[6] = cxsub(cxsub(beta0,rInc),iInc);
        betas[7] = cxadd(cxsub(beta0,rInc),iInc);
        betas[8] = cxsub(cxadd(beta0,rInc),iInc);
    
        double abschi[numBetas];
        int j=0;
        for (j=0; j<numBetas; j++)
        {
            abschi[j] = abschi_find(wavelength, thicknesses, indexesReal, indexesImag, 
                                    numLayers, real(betas[j]), imag(betas[j]));
        }
        chiMinIdx = argmin(abschi, numBetas);
        abschiOld = abschiNew;
        abschiNew = abschi[chiMinIdx];
        beta0 = betas[chiMinIdx];

        //abschidelta was added in this way because abs() doesn't seem to give correct answer in Mac OS X
        absdeltachi = abschiNew > abschiOld ? abschiNew-abschiOld : abschiOld-abschiNew;

    }while(absdeltachi > beta_find_precision);

    *betaOut = real(beta0);
    *(betaOut+1) = imag(beta0);
}

