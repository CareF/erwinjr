#include <stdio.h>
#include <math.h>
#include <stdint.h>
#ifdef __MP /*openmp support*/
#include <stdlib.h>
#include <omp.h>
#endif 

#define ANG 1E-10 /*angstrom in m */
#define KVpCM 1E5 /*kV/cm in V/m (SI unit) */
#define sq(X) (X)*(X)

const double hbar = 6.626e-34 / 2 / 3.141592653589793;
const double m0 = 9.109e-31;
const double e0 = 1.602e-19;
const double pi = 3.1415926535897932385;

#ifdef _WINDLL
typedef int32_t numpyint;
#else
typedef int64_t numpyint;
#endif // _WINDLL

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WINDLL
	__declspec(dllexport)
#endif // _WINDLL
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

#ifdef _WINDLL
	__declspec(dllexport)
#endif // _WINDLL
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
	int q;
#ifdef __MP
#pragma omp parallel private(xMcE, xPsi)
	{
		xMcE = (double *)malloc(xPsiSize * sizeof(double));
		xPsi = (double *)malloc(xPsiSize * sizeof(double));
#pragma omp for
#endif
	for(q=0; q<eEqSize; q++)
	{
#ifdef __MP
#endif
		double Eq = eEq[q];
		/* set start point, according to energy offset and external field */
		int startpoint = xPsiSize - ceil(
				(Eq - eEq[0])/(EField * KVpCM * ANG * xres) + extLength/xres);
		if(startpoint<1) 
			startpoint = 1;

		psiFn(Eq, startpoint, xPsiSize, xres, 
				xVc, xEg, xF, xEp, xESO, xMc, xMcE, xPsi);
		xPsiEnd[q] = xPsi[xPsiSize-1];
		//printf("%d: %g %d        ", q, eEq[q], startpoint);
		//printf("%d  ", startpoint);
	}
#ifdef __MP
		free(xMcE);
		free(xPsi);
	}
#endif

	return 1;
}

#ifdef _WINDLL
	__declspec(dllexport)
#endif // _WINDLL
int inv_quadratic_interp(const double *x, const double *y, 
		const numpyint *idxs, int idxLength, double *root)
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

#ifdef _WINDLL
__declspec(dllexport)
#endif // _WINDLL
int psiFill(int xPsiSize, double xres, int EigenESize, const double *EigenE, 
		const double *xVc, const double *xEg, const double *xF, 
		const double *xEp, const double *xESO, const double *xMc, 
		double *xMcE, double *xyPsi)
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

#ifdef _WINDLL
__declspec(dllexport)
#endif // _WINDLL
double inv_tau_int(int xPsiSize, double xres, double kl, 
		const double * xPoints, const double *psi_i, const double *psi_j) {
	/* To calculate the LO phonon life time between two given wave functions,
	 * The integral part
	 * see Eq.(2.65) in Kale's*/
	double Iij = 0;
	int i;
#ifdef __MP
#pragma omp parallel for reduction(+:Iij)
#endif
	for(i=0; i < xPsiSize; i++){
		for(int j=0; j < xPsiSize; j++){
			Iij += psi_i[i]*psi_j[i] * exp(-kl*ANG*fabs(xPoints[i] - xPoints[j]))
					* psi_i[j] * psi_j[j]; 
		}
	}
	return Iij * sq(xres*ANG);
}


#ifdef _WINDLL
__declspec(dllexport)
#endif // _WINDLL
int inv_alpha()
{return 137;}

#ifdef __cplusplus
}
#endif
