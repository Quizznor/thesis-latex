/*
 \file util.h
 \author Ugo Giaccari (u.giaccari@astro.ru.nl)
 \author Claudio Galelli 
 \author Lorenzo Caccianiga
 \date 20 June 2021
 \last modification 4 May 2022
 */

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <limits>

#include <cmath>
#include <string>
#include <math.h>  // fmod
#include <sstream>
#include <time.h>

//---- ROOT -----//
#include "Riostream.h"
#include "TRandom.h"
#include "TH2F.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1D.h"
#include "TPad.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TRint.h"
#include "TStyle.h"
#include "TMarker.h"
#include "TSystem.h"

#include <sys/time.h>

#include "STClibrary.h"
#include "healpixmap.h"
#include "../Data/DataPath.h"

using namespace std;

//-----------------------------------//
// -------- constants ---------------//
//-----------------------------------//
const double d2r = TMath::Pi() / 180.;
const double r2d = 180 / TMath::Pi();

// Number of events
const int nEventsOfficial = 2635;
// Binning in angular scale
// dTheta = 0.25° (for 1° < theta <= 5°) | dTheta = 1.° (for theta > 5° &&  theta < 5°)
const int nbin = 42;
const int nbinbl = 30;
const double dx1 = 0.25;
const double dx2 = 1.;
const int nEbin = 49; // number of bins in energy, from 32 to 80 EeV in step of 1 EeV
const double minEn = 32.; // minimum energy
const double maxEn = 81.; // maximum energy (in the plot)

// default number of simulated skies for local p-value computation
const int nsimBl=1E7; // for Blind it is the number of simulated events rather than skies
const int nsimAC=1E5; // Autocorrelation is quite slow
const int nsimCA=1E7; // Centaurus A has a lower local p-value
const int nsimT=1E6; // all the other analyses

// default number of simulation for penalized p-values
const int nsimBlPen=1E4; // Blind 
const int nsimACPen=1E4; // Autocorrelation is quite slow
const int nsimCAPen=5E6; // centaurus A has a lower local p-value
const int nsimTPen=1E5; // all the other analyses

// declination limit
const double maxDec = 45.0;

// nSide of the Healpix map
const int nSide = 64;

// exposure function
const string fileInExpo= "../../Data/"+exponame;

//! Get Seed
void GetSeed(unsigned int &seed);

//! Get angular distance btw two points on a sphere
double GetAngDist(double alphaEv, double deltaEv, double alphaObj, double deltaObj);

//! Print on the screen some info about the analysis 
void PrintInfo(string type, int nsky);

//! Set file for storing local p-value results
string SetNameFile(string type);

//! Check existence exposure file
bool CheckExpo();

//! Set declinations of the healpix grid for dec < 45°
vector<double> SetDecAxis(int nSide);

//! Get ipix of the pixels in the Auger fov (dec < 45°)
vector<int> GetIpix(int nSide);

//! Set binning in search radius
void SetAngScaleBinning(double *xbins);

//! Get binning in search radius for blind searches
void SetAngScaleBinningBl(double *xbins);

//! Set binning in energy
void SetEbin(double *enBin);

//! Get index for number events/pairs at each Eth and search radius
int GetIndex(double angBin, double enBin);

//! Set the object coordinates
void SetObjectCoords(string type, double &alphaObj, double &deltaObj);

//! Integrate th2 histo (E >= Eth , angular scale <=\Phi)
TH2F GetNormIntegral(TH2F h2, double nsky);

//! Get number of events in each energy bin
void GetNEbin(vector <double> enEv, double *nOut);

//! Get number of events for each threshold
void GetNETh(vector <double> enEv, double *nOut);

//! Get observed number of events/pairs
TH2F GetTH2data(string type);

//! Calculate binomial probability histo for Centaurus A
void GetBinProbCen(TH2F &hData, TH2F &hSim, TH2F &h2);

//! Get binomial probability
double GetBinProb(double nSim, double nData, double nTot);

//! Get minimum p-value
void GetMinimum(TH2F &hData, TH2F &hMean, TH2F &hProb);

//! check allowed analyses
bool CheckAnalysis(string type);

//! Li & Ma estimator
double GetLiMa(double nOn, double nOff, double alpha);

//! Progress bar
void ProgressBar(double progress);
