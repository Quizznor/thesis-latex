/*
 \file util.cc
 \author Ugo Giaccari (u.giaccari@astro.ru.nl)
 \author Claudio Galelli
 \author Lorenzo Caccianiga
 \date 20 June 2021
 \last modification 4 May 2022
 */

#include "utils.h"

void GetSeed(unsigned int &seed)
{
  //----- seed changed --------//
  gRandom->SetSeed(0);
  struct timeval myTimeVal;
  struct timezone myTimeZone;
  gettimeofday(&myTimeVal, &myTimeZone);
  seed = (unsigned int) (myTimeVal.tv_usec+(myTimeVal.tv_sec % 1000)*1000000);
}

double GetAngDist(double alphaEv, double deltaEv, double alphaObj, double deltaObj)
{
  double deltaEvR = deltaEv*d2r;
  double deltaObjR = deltaObj*d2r;
  double diffAlpha = alphaEv - alphaObj;
  double da = sin (deltaEvR) * sin (deltaObjR );
  double dg = cos (deltaEvR) * cos (deltaObjR) * cos( (diffAlpha) * d2r );
  double argAcos = da + dg;
  return acos (argAcos) * r2d;
}

void PrintInfo(string type, int nsky)
{
  cout << "\n";
  cout << " --------------------------------------------- \n";
  if (!type.compare("ac"))
    cout << " Auto-correlation analysis " << "\n";
  if (!type.compare("gc"))
    cout << " Analysis around the Galactic Center " << "\n";
  else if (!type.compare("gp"))
    cout << " Analysis around the Galactic Plane " << "\n";
  else if (!type.compare("sgp"))
    cout << " Analysis around the Super Galactic Plane " << "\n";
  else if (!type.compare("ca"))
    cout << " Analysis around Centaurus A Region " << "\n";
  cout << " Number of simulated isotropic skies " << nsky << endl;
  cout << " --------------------------------------------- \n";
  cout << "\n";
}

string SetNameFile(string type)
{
  // create Results directory
  if (gSystem->AccessPathName("Results"))
    gSystem->Exec("mkdir Results");
  string file;
  if (!type.compare("ac"))
    file = "Results/Figure2AC.root";
  if (!type.compare("gc"))
    file = "Results/Figure2GC.root";
  else if (!type.compare("gp"))
    file = "Results/Figure2GP.root";
  else if (!type.compare("sgp"))
    file = "Results/Figure2SGP.root";
  else if (!type.compare("ca"))
    file = "Results/Figure5.root";
  return file;
}

bool CheckExpo()
{
  if(gSystem->AccessPathName(fileInExpo.c_str())){
    cout << "\n";
    cout << " Exposure file does not exist. Run ./exposure.exe in the Data folder." << endl;
    cout << "\n";
    return FALSE;
  } else {
    return TRUE;
  }
}

vector<double> SetDecAxis(int nSide)
{
  THealpixMap map(nSide,'Q');
  unsigned int npix = map.NPix();
  vector<double> decrange;
  double ndec = 0.;
  for (unsigned int i = (npix-1); i > 0; --i) {
    double ra = 0.;
    double dec = 0.;
    map.GiveLB(i, ra, dec);
    if (dec<=maxDec) {
      if(dec != ndec) {
        ndec = dec;
        decrange.push_back(dec);
      }
    }
  }
  decrange.push_back(maxDec);
  return decrange;
}

vector<int> GetIpix(int nSide)
{
  THealpixMap map(nSide,'Q');
  unsigned int npix = map.NPix();
  vector<int> ipix;
  for (unsigned int i = 0; i < npix; ++i) {
    double ra = 0.;
    double dec = 0.;
    map.GiveLB(i, ra, dec);
    if (dec<=maxDec) {
      ipix.push_back(i);
    }
  }
  ipix.push_back(npix + 1);
  return ipix;
}

void SetAngScaleBinning(double *xbins)
{
  for (int i=0; i<=nbin; ++i) {
  if (i<2)
    xbins[i] = i * dx2;
  else if (i>=2 && i<=16)
    xbins[i] = 1+  (i-1) * dx1;
  else
    xbins[i] = 5 + (i-17) * dx2;
  }
}

void SetAngScaleBinningBl(double *xbins)
{
  for (int i=0; i<=nbinbl; ++i)
    xbins[i] = (i * dx2);
}

void SetEbin(double *enBin)
{
  for (int z = 0; z <= nEbin; ++z)
    enBin[z] = minEn + z;
}

int GetIndex(double angBin, double enBin)
{
  int index = ( (floor(enBin) - 1) * nbin ) + floor(angBin) -1;
  return index;
}

void SetObjectCoords(string type, double &alphaObj, double &deltaObj)
{
  if (!type.compare("ca")) { // CentaurusA equatorial coordinates
    alphaObj = 201.4;
    deltaObj = -43.0;
  }
  else { // Gal Center 
    alphaObj = 0.;
    deltaObj = 0.;
  }
}

TH2F GetNormIntegral(TH2F h2, double nsky)
{
  TH2F hFinal(h2);
  double content;
  for (int i=1; i<=nEbin; ++i) { // energy loop
    for (int j=1; j<=nbin; ++j) { // angular scale loop
      content = h2.Integral(1,j,i,nEbin);
      hFinal.SetBinContent(j, i, content / double (nsky) );
    } // end angular scale loop
  } // end energy loop
  return hFinal;
}

// Get number of events in each differential energy bin
void GetNEbin(vector<double> enEv, double *nOut)
{
  int nEvents=enEv.size();
  for (int i=0; i<nEbin; ++i) {
    nOut[i] = 0;
    for ( int j =0; j<nEvents; ++j) {
      if (i!=(nEbin-1)) { //exclude the last bin above 80 EeV
        if (enEv[j]>= (minEn+i) && enEv[j]< (minEn+1+i) )
          ++nOut[i];
      }
      else {
        if (enEv[j]>= (minEn+i))
        ++nOut[i];
      }
    }
  }
}

// Get number of events for each threshold
void GetNETh(vector <double> enEv, double *nOut)
{
  int nEvents=enEv.size();
  for (int i=0; i<nEbin; ++i) { // energy loop
    nOut[i] = 0;
    double eth = i + minEn;
    for (int z=0; z<nEvents; ++z) { // event loop
      if (enEv[z] >= eth)
        ++nOut[i];
    } // end event loop
  } // end energy loop
}

TH2F GetTH2data(string type)
{
  // set object coordinates
  double alphaObj = 0.;
  double deltaObj = 0.;
  SetObjectCoords(type, alphaObj, deltaObj);

  // set angular scale bin
  double xbins[nbin+1];
  SetAngScaleBinning(xbins);

  vector <double> vId,vZen,vTheta,alphaEv,deltaEv,enEv,vexp;
  string pathFile= "../../Data/"+dataname;
  if (!type.compare("ca") || !type.compare("ac") || !type.compare("bl")) 
    LoadAugerData(pathFile , vId, vZen, alphaEv,deltaEv,enEv, vexp,"eq");
  else if (!type.compare("gc") || !type.compare("gp"))  
    LoadAugerData(pathFile , vId, vZen, alphaEv,deltaEv,enEv, vexp,"gal");
  else if (!type.compare("sgp")) 
    LoadAugerData(pathFile , vId, vZen, alphaEv,deltaEv,enEv, vexp,"sgal");
  // TH2 map histo: differential (hDataTmp) and integrated form (hData)
  TH2F hData("hData","hData",nbin,xbins,nEbin,minEn,maxEn);
  TH2F hDataTmp("hDataTmp","hDataTmp",nbin,xbins,nEbin,minEn,maxEn);
  int nEvents=enEv.size();
  
  // for 2-point correlations
  if (!type.compare("ac")) {
    for (int zz=0; zz<nEvents; ++zz) {
      for (int zj=zz+1; zj<nEvents; ++zj) {
        double gamma = 0.;
        gamma = GetAngDist(alphaEv[zz], deltaEv[zz], alphaEv[zj], deltaEv[zj]);
        if (enEv[zz]<= 80) 
          hDataTmp.Fill(gamma, enEv[zz]);
        else 
          hDataTmp.Fill(gamma, 80);   
        }
      }
  } else { // for the other analyses
    for (int z=0; z<nEvents; ++z) {
      double gamma = 0.;
      // choose the angular distance: from an object (CenA, GC) or from a plane (GP, SGP)
      if (!type.compare("ca")) {
        gamma = GetAngDist(alphaEv[z], deltaEv[z], alphaObj, deltaObj);
      } else if (!type.compare("gc")) {
        gamma = GetAngDist(alphaEv[z], deltaEv[z], alphaObj, deltaObj);
      } else if (!type.compare("gp")) {
        gamma = abs(deltaEv[z]);
      } else if (!type.compare("sgp")) {
        gamma = abs(deltaEv[z]);
      }
      // fill the histo
      if (enEv[z] <= 80)
        hDataTmp.Fill(gamma, enEv[z]);
      else
        hDataTmp.Fill(gamma, 80);
    }
  }
  // final TH2 map: observed events as a function of energy and angular scale
  hData  = GetNormIntegral(hDataTmp, 1);
  hData.SetName("hData");

  return hData;
}

void GetMinimum(TH2F &hData, TH2F &hMean, TH2F &hProb)
{  
  // search for the minimum
  int nbinx = 0;
  int nbiny = 0;
  int nbinz = 0;

  hProb.GetMinimumBin(nbinx, nbiny, nbinz); 

  double x = hProb.GetXaxis()->GetBinLowEdge(nbinx + 1); 
  double y = hProb.GetYaxis()->GetBinLowEdge(nbiny);
  double nData = hData.GetBinContent(nbinx, nbiny); 
  double nExp  = hMean.GetBinContent(nbinx, nbiny); 
  double min = hProb.GetMinimum();
    
  // if p-value = 0 simulate more isotropic skies
  if (min==0)
    cout << " More simulations needed to estimate the local p-value!" << endl;
  else {
    cout << endl;
    cout << endl;
    cout << " ---- Results ----------------------------------------------------------------- " << endl;
    cout << " Local p-value " << min << endl;
    cout << " Energy threshold / EeV " << y << " Angular scale / deg " << x << endl;
    cout << " Observed number of events " << nData
      << " Expected number of events (on average) " << nExp << endl;
    cout << " ------------------------------------------------------------------------------ " << endl;
    cout << endl;
    cout << "\n";
    cout << " Now go ahead with the penalization! " << "\n";
    cout << "\n";
  }
}

void GetBinProbCen(TH2F &hData, TH2F &hSim, TH2F &h2)
{
  // get data
  vector <double> vId,vZen,vTheta,alphaEv,deltaEv,enEv,vexp;
  string pathFile= "../../Data/"+dataname;
  LoadAugerData(pathFile , vId, vZen, alphaEv,deltaEv,enEv, vexp,"eq");

  // Get number of events for each threshold
  double nOut[nEbin];
  GetNETh(enEv, nOut);

  // binomial probability
  for (int i=1; i<=nEbin; ++i) {
    double nTot = nOut[i-1];
    for (int j=1; j<=nbin; ++j) {
      double prob = 0.;
      double nData = hData.GetBinContent(j,i);
      double nSim = hSim.GetBinContent(j,i);
      if (nData)
        prob = TMath::BetaIncomplete(nSim / nTot, nData, nTot - nData + 1);
      else
        prob = 1.;
      h2.SetBinContent(j,i,prob);
    } // end angular scale loop
  } // end energy loop

  // print the results on the screen
  GetMinimum(hData, hSim, h2);
}

bool CheckAnalysis(string type)
{
  if ( (!type.compare("ca") ) && (!type.compare("gc")) &&
       (!type.compare("gp")) && (!type.compare("sgp")) ) {
    cout << endl;
    cout << " Type of allowed analyses: ca (CenA), gc (Galactic Center), gp (GalPlane), sgp (SGalPlane) " 
      << endl;
    cout << endl;
    return false;
  }
  else
    return true;
}

double GetBinProb(double nSim, double nData, double nTot)
{
  double prob = 0.;
  if (nData)
    prob = TMath::BetaIncomplete(nSim / nTot, nData, nTot - nData + 1);
  else
    prob = 1;
  return prob;
}

double GetLiMa(double nOn, double nOff, double alpha)
{
  double excess = nOn - (nOff * alpha);
  double n1 = 0.;
  double n2 = 0.;
  double sign = 0.;
  double a = 0.;

  if (excess > 0) {
    n1 = nOn;
    n2 = nOff;
    sign = 1.0;
    a = alpha;
  }
  else {
    n2 = nOn;
    n1 = nOff;
    sign = -1.0;
    a = 1.0 / alpha;
  }
  double nTot = n1 + n2;
  double t1 = n1 * log( ((1+a)/a) * (n1/nTot) );
  double t2 = n2 * log( (1+a) * (n2/nTot) );
  double significance = sqrt(2.) * sqrt( t1 + t2 );

  return sign*significance;
}

void ProgressBar(double progress)
{
  int barWidth = 70;

  std::cerr << " [";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos)
      cerr << "=";
    else if (i == pos)
      cerr << ">";
    else
      cerr << " ";
  }
  cerr << "] " << int(progress * 100.0) << " %\r";
  cerr.flush();
}
