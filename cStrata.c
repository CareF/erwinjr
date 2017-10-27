#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include "complex.h"
#ifdef __MP /*openmp support*/
#include <stdlib.h>
#include <omp.h>
#endif 

const double pi = 3.1415926535897932385;

#ifdef __cplusplus
extern "C" {
#endif

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

#ifdef _WINDLL
__declspec(dllexport)
#endif // _WINDLL
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

#ifdef _WINDLL
__declspec(dllexport)
#endif // _WINDLL
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

#ifdef _WINDLL
__declspec(dllexport)
#endif // _WINDLL
int argmin(const double *values, int numElements)
{
    int idx=0;
    double minValue = values[0];
    int j=0;
    for (j=1; j<numElements; j++)
        idx = values[j] < minValue ? j : idx;
    return idx;
}

#define numBetas 9

#ifdef _WINDLL
__declspec(dllexport)
#endif // _WINDLL
void beta_find(double wavelength, const double *thicknesses, const double *indexesReal, 
               const double *indexesImag, int numLayers, double betaInReal, double betaInImag,
               double beta_find_precision, double *betaOut)
{
    //initialize beta0
    complex beta0 = cmplx(betaInReal, betaInImag);

    complex rInc = cmplx(0.0001,0);
    complex iInc = cmplx(0,1.0e-6);

    //const int numBetas = 9;
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


#ifdef _WINDLL
__declspec(dllexport)
#endif // _WINDLL
int answer()
{return 42;}

#ifdef __cplusplus
}
#endif
