//Combined FLBM with MC coupled through shared kappa and theta.

//Need to optimize a great deal.
/***************************************************/
/*                      TO DO                      */
/* -Second moment                                  */
/* -Examine Discrepancy at different lattice Sizes */
/* -Compare with diffuse D1Q3                      */
/***************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <complex.h>
#include <fftw3.h>
#include <mygraph.h>

#define xdim 32

//LB arrays
double n[xdim];
double f0[xdim], f1[xdim], f2[xdim];
double feq0[xdim], feq1[xdim], feq2[xdim];
double M0[xdim], M1[xdim], M2[xdim];
double mu[xdim];
double u[xdim];
double Ak[xdim];
double nk[xdim], nkave[xdim];
double Fk0nid[xdim], Fk1nid[xdim], Fk0ave[xdim], Fk1ave[xdim];
double w[3];
double psi[3];

//LB Parameters
double n0 = 100;
double theta = 1./3.;
double kap = 0.;
double A = 0;
double B = 0;
double aven;
double sumn = 0;
double tau[3] = {1., 1., 1.};
double forcescale = 1.;
double noisescale = 1.;
double aven2lb;
int v[3] = {0, 1, -1};
int forcecheck = 0;
int noisecheck = 1;
int iterations = 0;

//MC Arrays
int nMC[xdim];
double NMC[xdim];
double nkMC[xdim], nkaveMC[xdim];
double Fk0nidMC[xdim], Fk1nidMC[xdim], Fk0aveMC[xdim], Fk1aveMC[xdim];
double AkMC[xdim];
double aven2MC;

//MC Parameters
int n0MC = 100;
int sumnMC = 0;
int avenMC;

//Histogram Arrays
double Poisson[200];
int MChist[200];
double MChistprint[200];
int LBhist[200];
double LBhistprint[200];
double Moment2Th;

//GUI controls
int repeat = 1;
int next = 0;
int Pause = 1;
int done = 0;

//FFTW variables
int nFFTcount = 0;
int nMCFFTcount = 0;

//FFTW declarations
fftw_complex *in, *out, *inMC, *outMC;
fftw_plan plan, planMC;

//Setup FFT
void SetupFFT() {
  in = fftw_alloc_complex(xdim);
  out = fftw_alloc_complex(xdim);
  plan = fftw_plan_dft_1d(xdim, in, out, -1, FFTW_ESTIMATE);
  inMC = fftw_alloc_complex(xdim);
  outMC = fftw_alloc_complex(xdim);
  planMC = fftw_plan_dft_1d(xdim, inMC, outMC, -1, FFTW_ESTIMATE);
}

//Function to calculate the FFT
int FourierXform(double g[xdim], int count) {
  count++;
  double sum = 0;
  for (int i = 0; i < xdim; i++) {
    sum += g[i];
    in[i] = g[i];
  }
  sum /= xdim;

  for (int i = 0; i < xdim; i++) {
    in[i] -= sum;
    in[i] = in[i]/xdim;
  }

  fftw_execute(plan);

  for (int i = 0; i < xdim; i++ ) {
    nkave[i] = (nkave[i] * count + cabs(out[i])*cabs(out[i]))/((count+1));
    nk[i] = cabs(out[i]) * cabs(out[i]);
  }
  return count;
}

//Reset Fourier counter
void ResetAveCount() {
  nFFTcount = 0;
  nMCFFTcount = 0; // Add in later
}

void Histogram() {
  for (int i = 0; i < xdim; i++) {
    MChist[nMC[i]] =  MChist[nMC[i]] + 1;
    LBhist[(int)round(n[i])] = LBhist[(int)round(n[i])] + 1;
    MChistprint[nMC[i]] = MChist[nMC[i]]/((double)iterations*xdim);
    LBhistprint[(int)round(n[i])] = LBhist[(int)round(n[i])]/((double)iterations*xdim);
  }
}

void FEAnalysis(double A[xdim], double F0out[xdim], double F1out[xdim], int m) {
  double sum0nid = 0;
  double sum1nid = 0;
  for (int i = 0; i < xdim; i++) {
    sum0nid = sum0nid + kap*A[m]*A[m]*sin(2.*M_PI*m*i/xdim)*(sin(2.*M_PI*m*i/xdim)-0.5*(sin(2.*M_PI*m*((i-1+xdim)%xdim)/xdim)+sin(2.*M_PI*m*((i+1)%xdim)/xdim)));
    sum1nid = sum1nid + kap/4.*A[m]*A[m]*(pow(sin(2.*M_PI*m*((i+1)%xdim)/xdim),2)+pow(sin(2.*M_PI*m*((i-1+xdim)%xdim)/xdim),2)-2.*sin(2.*M_PI*m*((i+1)%xdim)/xdim)*sin(2.*M_PI*m*((i-1+xdim)%xdim)/xdim));
  }
  F0out[m] = sum0nid;
  F1out[m] = sum1nid;
}

void FEAveraging(double F0[xdim], double F1[xdim], double F0ave[xdim], double F1ave[xdim]) {
  for (int k = 0; k < xdim; k++) {
    F0ave[k] = (F0ave[k]*nFFTcount+F0[k])/(nFFTcount+1);
    F1ave[k] = (F1ave[k]*nFFTcount+F1[k])/(nFFTcount+1);
  }

}

void SecondMoment() {
  double sumlb = 0;
  double summc = 0;
  Moment2Th = n0*n0+n0;
  for (int i = 0; i < xdim; i++) {
    sumlb = sumlb + n[i]*n[i];
    summc = summc + nMC[i]*nMC[i];
  }
  aven2lb = sumlb/((double)xdim);
  aven2MC = summc/((double)xdim);
  
}

void Averaging() {
  sumn = 0;
  for (int i = 0; i < xdim; i++) {
    sumn += n[i];
  }
  aven = sumn/((double)xdim);
}

void AveragingMC() {
  sumnMC=0;
  for (int i = 0; i < xdim; i++) {
    sumnMC += nMC[i];
  }
  avenMC = sumnMC/((double)xdim);
}

void Weights(){
  w[0] = 1. - theta;
  w[1] = theta/2.;
  w[2] = theta/2.;
}

void SetInteractionParameters() {
  psi[1] = A + kap;
  psi[0] = psi[2] = -kap/2.;
}

//Calculate Chcmical Potential
void MU(int x) {
  mu[x] = theta*log(n[x]) + 2.*n[x]*psi[1] + psi[0]*(n[(x-1+xdim)%xdim] + n[(x+1)%xdim]) + psi[2]*(n[(x-1+xdim)%xdim] + n[(x+1)%xdim]);
}

//Calculate velocity space equilibrium distributions
void SetEqDist(int x) {

  feq0[x] = w[0] * n[x] * (1. - u[x] * u[x]/(2.*theta));
  feq1[x] = w[1] * n[x] * (1. + u[x]/theta + u[x] * u[x]/(2.*theta*theta) - u[x] * u[x]/(2.*theta));
  feq2[x] = w[2] * n[x] * (1. - u[x]/theta + u[x] * u[x]/(2.*theta*theta) - u[x] * u[x]/(2.*theta));
}

//Calculate Gradients
double Gradient(double g[xdim], int x) {
  double grad = 0;
  for (int c = 0; c < 3; c++) {
    grad -= w[c]*v[c]*g[(x+xdim+v[c])%xdim];
  }
  return 3.*grad;
}

double CalculateAk(double g[xdim], int m) {
  double sum = 0;
  for (int i = 0; i < xdim; i++) {
    sum = sum + g[i]*sin(2*M_PI*m*i/xdim);
  }
  return sum;
}

int FourierXformMC(double g[xdim], int count) {
  count++;
  int sum = 0;
  for (int i = 0; i < xdim; i++) {
    sum += g[i];
    inMC[i] = g[i];
  }
  sum /= xdim;

  for (int i = 0; i < xdim; i++) {
    inMC[i] -= sum;
    inMC[i] = inMC[i]/xdim;
  }

  fftw_execute(planMC);

  for (int i = 0; i < xdim; i++) {
    nkaveMC[i] = (nkaveMC[i] * count + cabs(outMC[i])*cabs(outMC[i]))/((count+1));
    nkMC[i] = cabs(outMC[i]) * cabs(outMC[i]);
  }

  return count;
}

double SetFreeEnergy(int g[xdim], int x) {
  double sum = 0;
  for (int deltax = -1; deltax < 2; deltax++) {
    sum = sum + g[x]*theta + g[x] * g[(x+deltax+xdim)%xdim] *psi[deltax+1];
  }
  return sum;
}

void initMC() {
  for (int i = 0; i < xdim; i++) {
    nMC[i] = n0MC;
    NMC[i] = nMC[i];
  }
  nMCFFTcount = FourierXformMC(NMC, nMCFFTcount);
  AveragingMC();
}

void MCCheck(int h[xdim]) {
  int hprime[xdim];
  for (int i = 0; i < xdim; i++) {
    hprime[i] = h[i];
  }
  for (int i = 0; i < xdim; i++) {
    int r = rand();
    int x1 = r & (xdim -1);
    r /= xdim;
    int nc = r % (xdim * n0MC);
    int x2 = 0;
    int nn = hprime[x2];
    while (nc > nn) {
      nn += hprime[++x2];
    }

    hprime[x1] += 1;
    hprime[x2] -= 1;

    double deltaFE = SetFreeEnergy(hprime, (x1-1+xdim)%xdim) - SetFreeEnergy(h, (x1-1+xdim)%xdim);
    deltaFE += SetFreeEnergy(hprime, x1) - SetFreeEnergy(h, x1);
    deltaFE += SetFreeEnergy(hprime, (x1+1)%xdim) - SetFreeEnergy(h, (x1+1)%xdim);
    deltaFE += SetFreeEnergy(hprime, (x2-1+xdim)%xdim) - SetFreeEnergy(h, (x2-1+xdim)%xdim);
    deltaFE += SetFreeEnergy(hprime, x2) - SetFreeEnergy(h, x2);
    deltaFE += SetFreeEnergy(hprime, (x2+1)%xdim) - SetFreeEnergy(h, (x2+1)%xdim);

    if (deltaFE < 0) {
      h[x1] = hprime[x1];
      h[x2] = hprime[x2];
    } else {
      double AcceptProbability = exp(-deltaFE/theta);
      double RandProb = (1.0 * rand()/RAND_MAX);
      if (AcceptProbability >= RandProb) {
        h[x1] = hprime[x1];
        h[x2] = hprime[x2];
      } else {
        hprime[x1] -= 1;
        hprime[x2] += 1;
      }
    }
  }
}

void iterationMC() {
  MCCheck(nMC);
  nMCFFTcount = FourierXformMC(NMC, nMCFFTcount);

  for (int i = 0; i < xdim; i++) {
    NMC[i] = nMC[i];
  }

  for (int k = 0; k < xdim; k++) {
    FEAnalysis(AkMC, Fk0nidMC, Fk1nidMC, k);
  }

  FEAveraging(Fk0nidMC, Fk1nidMC, Fk0aveMC, Fk1aveMC);
 
  AveragingMC();
}

void PoissonDist() {
  Poisson[0] = exp(-n0);
  for (int i = 0; i < 200; i++) {
    Poisson[i+1] = n0/(i+1)*Poisson[i];
  }
}

void init() {
  initMC(); 
  iterations = 0;
  PoissonDist(); // Do this after simulation works
  Weights();
  SetInteractionParameters();

  for (int i = 0; i < xdim; i++) {
    n[i] = n0;
    f0[i] = w[0] * n[i];
    f1[i] = w[1] * n[i];
    f2[i] = w[2] * n[i];
    u[i] = (f1[i] - f2[i])/n[i];
    nk[i] = 0;
    nkave[i] = 0;
    Fk0nid[i] = Fk1nid[i] = Fk0ave[i] = Fk1ave[i] = 0;
    Fk0nidMC[i] = Fk1nidMC[i] = Fk0aveMC[i] = Fk1aveMC[i] = 0;
  }

  nFFTcount = FourierXform(n, nFFTcount);

  for (int i = 0; i < xdim; i++) {
    SetEqDist(i);
    MU(i);
    //nktheory[i] = CalculateRhoKTheory(i);  // Do this after the actual simulation works
  }
  for (int i = 0; i < 200; i++) {
    MChist[i] = MChistprint[i] = 0;
    LBhist[i] = LBhistprint[i] = 0;
  }
  Averaging();
  SecondMoment();
  ResetAveCount();
}

void iteration() {
  SetInteractionParameters();
  iterationMC();
  Weights();
  PoissonDist();

  double Meq0[xdim], Meq1[xdim], Meq2[xdim];

  double m[3][3] = {
    {1., 1., 1.,},
    {0., sqrt(1./theta), -sqrt(1./theta)},
    {-sqrt(theta/(1.-theta)), sqrt((1.-theta)/theta), sqrt((1.-theta)/theta)}
  };

  for (int i = 0; i < xdim; i++) {
    n[i] = f0[i] + f1[i] + f2[i];
    u[i] = (f1[i] - f2[i])/n[i];
    SetEqDist(i);
    M0[i] = M1[i] = M2[i];
    Meq0[i] = Meq1[i] = Meq2[i] = 0;

    //Forward Transform 
    M0[i] += m[0][0] * f0[i];
    M0[i] += m[0][1] * f1[i];
    M0[i] += m[0][2] * f2[i];
    M1[i] += m[1][0] * f0[i];
    M1[i] += m[1][1] * f1[i];
    M1[i] += m[1][2] * f2[i];
    M2[i] += m[2][0] * f0[i];
    M2[i] += m[2][1] * f1[i];
    M2[i] += m[2][2] * f2[i];

    Meq0[i] += m[0][0] * feq0[i];
    Meq0[i] += m[0][1] * feq1[i];
    Meq0[i] += m[0][2] * feq2[i];
    Meq1[i] += m[1][0] * feq0[i];
    Meq1[i] += m[1][1] * feq1[i];
    Meq1[i] += m[1][2] * feq2[i];
    Meq2[i] += m[2][0] * feq0[i];
    Meq2[i] += m[2][1] * feq1[i];
    Meq2[i] += m[2][2] * feq2[i];

    //Calculate and Include Noise
    double noise[3];

    //Only adding noise to non-conserved quantities...
    //Only 1 non-conserved quantity related to v_i*v_i.

    noise[0] = 0;
    noise[1] = 0;
    if (n[i] > 0) {
      noise[2] = noisescale*sqrt(n[i]*12.)*((double)rand()/RAND_MAX - 0.5);
    } else {
      noise[2] = 0;
    }
    if (noisecheck == 0) {
      noise[2] = 0;
    }

    //Collision in moment space
    M0[i] = 1./tau[0]*(Meq0[i] + (tau[0] - 1.) * M0[i] + sqrt(2.*tau[0] - 1.)*noise[0]);
    M1[i] = 1./tau[1]*(Meq1[i] + (tau[1] - 1.) * M1[i] + sqrt(2.*tau[1] - 1.)*noise[1]);
    M2[i] = 1./tau[2]*(Meq2[i] + (tau[2] - 1.) * M2[i] + sqrt(2.*tau[2] - 1.)*noise[2]);

    f0[i] = f1[i] = f2[i] = 0;

    //Back Transform
    f0[i] += m[0][0] * w[0] * M0[i];
    f0[i] += m[1][0] * w[0] * M1[i];
    f0[i] += m[2][0] * w[0] * M2[i];
    f1[i] += m[0][1] * w[1] * M0[i];
    f1[i] += m[1][1] * w[1] * M1[i];
    f1[i] += m[2][1] * w[1] * M2[i];
    f2[i] += m[0][2] * w[2] * M0[i];
    f2[i] += m[1][2] * w[2] * M1[i];
    f2[i] += m[2][2] * w[2] * M2[i];
  }

  for (int i = 0; i < xdim; i++) {
    MU(i);
    //nktheory[i] = CalculateRhoKTheory(i);
  }

  //Calculate and include forcing terms
  if (forcecheck == 1) {
    for (int i = 0; i < xdim; i++) {
      double FF = 0;
      double FFn = 0;
      double F0, F1, F2;
      //Try using (n[i+1] + n[i-1])/2 in stead of n[i]
      FF = (n[(i+1)%xdim]-n[(i-1+xdim)%xdim])/2. * (1. - 1./2.)*Gradient(mu,i);  // This is the kappa offending function...
      //FF = (1. - 1./2.)*Gradient(mu,i);      // When using this without the density prefactor, we have stability for a large range of kappa. upto kap = 0.5
      FFn = (1.-1./2.)*theta*Gradient(n,i);

      F0 = 0;
      F1 = (FF - FFn)/2.;
      F2 = (FFn - FF)/2.;

      f0[i] += F0;
      f1[i] += F1;
      f2[i] += F2;
    }
  }

  //Calculate Fourier xforms
  nFFTcount = FourierXform(n, nFFTcount);
  for (int k = 0; k < xdim; k++) {
    Ak[k] = CalculateAk(n,k);
    AkMC[k] = CalculateAk(NMC, k);
  }

  //Streaming step
  double tmp1, tmp2;
  tmp1=f1[xdim-1];
  tmp2=f2[0];
  memmove(&f1[1],&f1[0],(xdim-1)*sizeof(double));
  memmove(&f2[0],&f2[1],(xdim-1)*sizeof(double));
  f2[xdim-1]=tmp2;
  f1[0] = tmp1;

  for (int k = 0; k < xdim; k++) {
    FEAnalysis(Ak, Fk0nid, Fk1nid, k);
  }

  FEAveraging(Fk0nid, Fk1nid, Fk0ave, Fk1ave);

  Averaging();
  SecondMoment();
  Histogram();

  iterations++;
}

void GUI() {
  static int Xdim = xdim;
  static int Pois = 200;

  DefineGraphN_R("n", &n[0], &Xdim, NULL);
  DefineGraphN_R("mu", &mu[0], &Xdim, NULL);
  DefineGraphN_R("u", &u[0], &Xdim, NULL);
  DefineGraphN_R("NMC", &NMC[0], &Xdim, NULL);
  NewGraph();
  DefineGraphN_R("nk", &nk[0], &Xdim, NULL);
  DefineGraphN_R("nkave", &nkave[0], &Xdim, NULL);
  DefineGraphN_R("nkMC", &nkMC[0], &Xdim, NULL);
  DefineGraphN_R("nkMCave", &nkaveMC[0], &Xdim, NULL);
  //DefineGraphN_R("nktheory", &nktheory[0], &Xdim, NULL);
  DefineGraphN_R("Fk0nidave", &Fk0ave[0], &Xdim, NULL);
  DefineGraphN_R("Fk1nidave", &Fk1ave[0], &Xdim, NULL);
  DefineGraphN_R("Fk0nidaveMC", &Fk0aveMC[0], &Xdim, NULL);
  DefineGraphN_R("Fk1nidaveMC", &Fk1aveMC[0], &Xdim, NULL);
  NewGraph();
  DefineGraphN_R("Poisson", &Poisson[0], &Pois, NULL);
  SetDefaultColor(2);
  DefineGraphN_R("MC Dist", &MChistprint[0], &Pois, NULL);
  SetDefaultColor(4);
  DefineGraphN_R("LB Dist", &LBhistprint[0], &Pois, NULL);

  StartMenu("Fluctuating D1Q3 LB and MC", 1);
    DefineInt("Iterations", &iterations);
    StartMenu("LB Parameters", 0);
      DefineDouble("n0", &n0);
      DefineDouble("kappa", &kap);
      DefineDouble("A", &A);
      DefineDouble("B", &B);
      DefineDouble("theta", &theta);
      DefineDouble("Noise Scale", &noisescale);
      DefineDouble("Force Scale", &forcescale);
      DefineBool("LB Noise", &noisecheck);
      DefineBool("LB Forcing", &forcecheck);
      StartMenu("Relaxation Times", 0);
        DefineDouble("tau[0]", &tau[0]);
        DefineDouble("tau[1]", &tau[1]);
        DefineDouble("tau[2]", &tau[2]);
      EndMenu();
      DefineDouble("Average n", &aven);
      DefineDouble("2nd Moment LB", &aven2lb);
      DefineDouble("2nd Moment Theory", &Moment2Th);
    EndMenu();
    StartMenu("MC Parameters", 0);
      DefineInt("n0MC", &n0MC);
      DefineDouble("kappa", &kap);
      DefineDouble("theta", &theta);
      DefineInt("Average nMC", &avenMC);
      DefineDouble("2nd Moment MC", &aven2MC);
      DefineDouble("2nd Moment Theory", &Moment2Th);
    EndMenu();
    DefineFunction("Initialize", &init);
    DefineFunction("Reset FFT Averaging", &ResetAveCount);
    SetActiveGraph(0);
    DefineGraph(curve2d_, "Rho Graphs");
    SetActiveGraph(1);
    DefineGraph(curve2d_, "FFT Graphs");
    SetActiveGraph(2);
    DefineGraph(curve2d_, "Poisson");
    DefineInt("Repeat", &repeat);
    DefineBool("Next", &next);
    DefineBool("Pause", &Pause);
    DefineBool("Close", &done);
  EndMenu();
}

int main(int argc, char *argv[]) {
  int newdata = 1;

  SetupFFT();
  init();
  GUI();            //Put this in last

  while (done == 0) {
    Events(newdata);
    DrawGraphs();
    if (next || !Pause) {
      newdata = 1;
      next = 0;
      for (int i = 0; i < repeat; i++) {
        iteration();
      }
    }
    else sleep(1);
  }

  return 0;
} 
