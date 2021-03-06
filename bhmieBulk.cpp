//%function bhmie calculates amplitude scatteing matrix elements and efficiencies for extinction,
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <complex>

using namespace std;

// Global Variable Prototypes
double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062;
complex <double> imaginary = sqrt(complex<double>(-1,0)); // initialize an imaginary number

// Function Prototypes
int bhmie(double x,double refrel,int nang, double* Qscat_p, double* Qext_p, double* Qback_p, complex<double>* S1_p, complex<double>* S2_p, double* diffPoint);
double trapz(double x[], double y[], int size);

// Main Function
int main(){

    // Variables //

    // User Inputs //

    // Refractive Index
    double refMed = 1.33;
    double refPart = 1.45;//(1.3, 0.008); // relative refractive index
    double refRel = refPart/refMed;

    // Wavelength
    double lambda = 0.532; // lidar wavelength in a vaccuum (um)
    double lambdaMed = lambda/refMed; // Lidar Wavelength in Medium (um)
    double kMed = 2*pi/lambdaMed; // Lidar wavenumber in Medium;

    // Angles
    int nang = 91; // number of angles between 0-90
    int nang1 = nang*2-1; // number of angles 0-180

    // Particle Size Distribution Parameters

    // Set PSD parameters
    double Dmin = 0.1; // minimum particle diameter (um);
    double Dmax = 150.0; // maximum particle diameter (um);
    int diamBin = 100; // # of diameter bins;
    double fac = pow((Dmax/Dmin),(1.0/(diamBin-1.0))); // exponential factor necessary for defining logarithmically spaced diameter bins

    // Initialize Particle Size Arrays
    double D[diamBin]; double r[diamBin]; double sizeParam[diamBin];
    double diffNumDistribution[diamBin];
    double* diffPoint = &diffNumDistribution[0];

    // Set PSD array values

    // Diameter Array
    for (int i; i<diamBin; i++){ // generate an array of particle diameters
      D[i]=Dmin*pow(fac,(i)); // define the diameter bins
    }
    // Radius Array
    for (int i; i<diamBin; i++){ // generate an array of particle diameters
      r[i]=D[i]/2; // define the diameter bins
    }
    // Size Parameter Array
    for (int i; i<diamBin; i++){ // generate an array of size parameters
      sizeParam[i] = 2*pi*r[i]*refMed / lambda; // mie theory size parameter
    }

    // Define Mie Output Variables and Pointers
    double Qscat; double Qext; double Qback;  // Mie scattering efficiencies
    double* Qscat_p = &Qscat; double* Qext_p = &Qext; double* Qback_p = &Qback; // pointers to Mie Scattering efficiencies

    complex<double> S1[2*nang]; complex<double> S2[2*nang];
    complex<double>* S1_p = S1; complex<double>* S2_p = S2;

    // Mueller Matrix elements
    double s11[2*nang-1][diamBin];
    double s12[2*nang-1][diamBin];
    double s33[2*nang-1][diamBin];
    double s34[2*nang-1][diamBin];

    double integrandS11[2*nang-1][diamBin];
    double integrandS12[2*nang-1][diamBin];
    double integrandS33[2*nang-1][diamBin];
    double integrandS34[2*nang-1][diamBin];

    double integrandArray11[diamBin];
    double integrandArray12[diamBin];
    double integrandArray33[diamBin];
    double integrandArray34[diamBin];

    double s11bar[2*nang-1];
    double s12bar[2*nang-1];
    double s33bar[2*nang-1];
    double s34bar[2*nang-1];

    // Define Distribution //
    double k = 5E18; // differential number concentration at particle size D0

    double jungeSlope = 4.0; // slope of the junge distribution

    for (int i = 0; i<diamBin; i++){
      diffNumDistribution[i] = k*pow((D[i]/D[0]),(-1*jungeSlope)); // # of particles m^-3 um^-1
    }


    // Mie Calculations for Each Size Parameter in the distribution
    //j+1 is used to convert from fortran indexing to c++indexing
    for (int i = 0; i<diamBin; i++){
      bhmie(sizeParam[i],refRel,nang,Qscat_p, Qext_p, Qback_p, S1_p, S2_p, diffPoint);

      for (int j = 0; j<2*nang-1; j++){
        s11[j][i] = 0.5 * (pow(abs(S2[j+1]),2) + pow(abs(S1[j+1]),2));
        s12[j][i] = 0.5 * (pow(abs(S2[j+1]),2) - pow(abs(S1[j+1]),2));
        s33[j][i] = real(S1[j+1]*conj(S2[j+1]));
        s34[j][i] = imag(S2[j+1]*conj(S1[j+1]));
      }
    }


    // Define integrand to calculate bulk mueller atrix properties
    for (int i=0; i<diamBin; i++){
      for (int j=0; j<2*nang-1; j++){
        integrandS11[j][i] = diffNumDistribution[i] * s11[j][i]; // good
        integrandS12[j][i] = diffNumDistribution[i] * s12[j][i]; // good
        integrandS33[j][i] = diffNumDistribution[i] * s33[j][i]; // good
        integrandS34[j][i] = diffNumDistribution[i] * s34[j][i]; // good
    }
  }
    // Checked throught here. SOmething May be wrong with the way you are doing things below
    // Integrate of Size Distribution For Each Angles
    for (int i=0; i<2*nang-1; i++){
      for (int j=0; j<diamBin; j++){
        integrandArray11[j] = integrandS11[i][j];
        integrandArray12[j] = integrandS12[i][j];
        integrandArray33[j] = integrandS33[i][j];
        integrandArray34[j] = integrandS34[i][j];

      }
      s11bar[i] = (1.0/(kMed*kMed))*trapz(sizeParam,integrandArray11,diamBin);
      s12bar[i] = (1.0/(kMed*kMed))*trapz(sizeParam,integrandArray12,diamBin);
      s33bar[i] = (1.0/(kMed*kMed))*trapz(sizeParam,integrandArray33,diamBin);
      s34bar[i] = (1.0/(kMed*kMed))*trapz(sizeParam,integrandArray34,diamBin);
    }


    // for (int i=0; i<2*nang-1; i++){
    //   cout << s12bar[i] << endl;
    // }











    //bhmie(x, refrel, nang, Qscat_p, Qext_p, Qback_p, S1_p, S2_p);

    // cout << Qscat << endl;
    // cout << Qext << endl;
    // cout << Qback << endl;
    // cout << "S1 = " << endl;
    //
    // for (int i = 1; i<=2*nang-1; i++){
    //   cout << S1[i] << endl;
    // }
    //
    // cout << "S2 = " << endl;
    //
    // for (int i = 1; i<=2*nang-1; i++){
    //   cout << S2[i] << endl;
    // }

    return 0;
}




//// Function Definitions
//
//
int bhmie(double x, double refrel, int nang, double* Qscat_p, double* Qext_p, double* Qback_p, complex<double>* S1_p, complex<double>* S2_p, double * diffPoint){

    // Variable Definitions
    complex<double> y; double dx; double nstop; double ymod; int nmx; double dang; double theta;
    int nn;
    double RN; double DN; double FN;
    complex<double> AN; complex<double> BN;
    double PSI; double PSI0; double PSI1;
    double CHI; double CHI0; double CHI1;
    double APSI; double APSI0; double APSI1;
    complex <double> XI; complex <double> XI0; complex <double> XI1;
    double P; double T;

    // Assign variables that will be exported to the value of their pointers
    double Qscat_f = *Qscat_p; double Qext_f = *Qext_p; double Qback_f = *Qback_p;

    // Array Definitions
    double PI[nang+1]; double PI0[nang+1]; double PI1[nang+1];
    complex<double> S1[2*nang+1]; complex<double> S2[2*nang+1];
    double AMU[nang+1]; double TAU[nang+1];

    dx = x;
    y = x * refrel;

    nstop = ceil(x + 4 * pow(x,0.3333) +2);
    ymod = abs(y);
    nmx = max(nstop,ceil(ymod))+15;
    complex<double> D[nmx];
    dang = pi/2/(nang-1);


    for (int i = 1; i<=nang; i++){
        theta = (double)(i-1) * dang;
        AMU[i] =  cos(theta);
    }

    //Logarithmic derivative D(j) calculated by downward recurence
    D[nmx]= complex <double>(0,0);
    nn = nmx-1;

    for (int n = 1; n<=nn; n++){
        RN = nmx-n+1;
        D[nmx-n]=(RN/y)-(1.0/(D[nmx-n+1]+RN/y));
    }

    for (int j = 1; j <= nang; j++){
        PI0[j] = 0.0;
        PI1[j] = 1.0;
    }

    nn = 2*nang-1;

    for (int j = 1; j <= nang; j++){
        S1[j] = 0.0;
        S2[j] = 0.0;
    }

    //Riccati-Bessel functins with real argument x calculated by upward recurence
    PSI0=cos(dx);
    PSI1=sin(dx);
    CHI0=-sin(x);
    CHI1=cos(x);
    APSI0=PSI0;
    APSI1=PSI1;
    XI0=APSI0-CHI0*imaginary;
    XI1=APSI1-CHI1*imaginary;
    Qscat_f=0.0;

    int n=1;
    while(n-1-nstop<0){
        DN=n;
        RN=n;
        FN=(2*RN+1)/(RN*(RN+1));
        PSI=(2*DN-1)*PSI1/dx-PSI0;
        APSI=PSI;
        CHI=(2*RN-1)*CHI1/x-CHI0;
        XI=APSI-CHI*imaginary;
        AN=((D[n]/refrel+RN/x)*APSI-APSI1)/((D[n]/refrel+RN/x)*XI-XI1);
        BN=((refrel*D[n]+RN/x)*APSI-APSI1)/((D[n]*refrel+RN/x)*XI-XI1);
        Qscat_f=Qscat_f+(2.0*RN+1.0)*(abs(AN)*abs(AN)+abs(BN)*abs(BN));


        for (int j=1; j<=nang; j++){
            int jj=2*nang-j;
            PI[j]=PI1[j];
            TAU[j]=RN*AMU[j]*PI[j]-(RN+1)*PI0[j];
            P=pow(-1.0,(n-1));
            S1[j]=S1[j]+FN*(AN*PI[j]+BN*TAU[j]);
            T=pow((-1),n);
            S2[j]=S2[j]+FN*(AN*TAU[j]+BN*PI[j]);
            if(j != jj){
                S1[jj]=S1[jj]+FN*(AN*PI[j]*P+BN*TAU[j]*T);
                S2[jj]=S2[jj]+FN*(BN*PI[j]*P+AN*TAU[j]*T);

            }
        }

        PSI0=PSI1;
        PSI1=PSI;
        APSI1=PSI1;
        CHI0=CHI1;
        CHI1=CHI;
        XI1=APSI1-CHI1*imaginary;
        n=n+1;
        RN=n;
        for (int j=1; j<=nang; j++){
            PI1[j]=((2*RN-1)/(RN-1))*AMU[j]*PI[j]-RN*PI0[j]/(RN-1);
            PI0[j]=PI[j];
        }
    }

    Qscat_f=(2.0/(x*x))*Qscat_f;
    Qext_f=(4.0/(x*x))*real(S1[1]);
    //Note: Qback is not Qbb, but the radar back scattering.
    Qback_f=(4.0/(x*x))*(abs(S1[2*nang-1])*abs(S1[2*nang-1]));
    //cout << "Qext_f = " << Qext_f << endl;


    // Move scattering efficiencies out of the function
    *Qscat_p = Qscat_f; // update the value of Qscat in main using a pointer
    *Qext_p  = Qext_f; // update the value of Qext in main using a pointer
    *Qback_p = Qback_f; // update the value of Qback in main using a pointer

    // Move S1 and S2 out of the function
    for (int i=1; i<=nang*2-1; i++){
      S1_p[i] = S1[i];
    }

    for (int i=1; i<=nang*2-1; i++){
      S2_p[i] = S2[i];
    }
    return 0;
}

double trapz(double x[], double y[], int size){
  double s;
  double sTemp = 0.0;
  for (int i = 0; i<size-1; i++){
    sTemp = sTemp+(x[i+1]-x[i])*(y[i+1]+y[i]);
  }
  s = 0.5*sTemp;
  return s;
}
