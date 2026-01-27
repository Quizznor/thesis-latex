/*
 \file 2_penalization.cc
 \author Ugo Giaccari (u.giaccari@astro.ru.nl)
 \author Claudio Galelli
 \author Lorenzo Caccianiga
 \date 20 June 2021
 \last modification 4 May 2022
 */

#include "../utils.h"

// Toolkit
#include "healpixmap.h"
#include "projmap.h"

void simSky(double *alphaSim, double *deltaSim, double *enSim)
{
  // Get data
  vector <double> vId,vZen,vTheta,alphaEv,deltaEv,enEv,vexp;
  string pathFile= "../../Data/"+dataname;
  LoadAugerData(pathFile , vId, vZen, alphaEv,deltaEv,enEv, vexp,"eq");

  // Get number of events in each energy bin
  double nEnbin[nEbin];
  GetNEbin(enEv, nEnbin);

  // Get seed
  unsigned int seed;
  GetSeed(seed);
  TRandom3 myR(seed);

  // read exposure function
  TFile *f = TFile::Open(fileInExpo.c_str(),"READ");
  TF1 *f1;
  f->GetObject("fcos",f1);

  int count = 0;
  // energy loop
  for (int j = 0; j<nEbin; ++j) {
    double eSim = double(j + minEn);
    // events loop
    int iCr = 0;
    while (iCr < nEnbin[j]){
      double raSim = myR.Rndm() *360.;
      double decSim = myR.Rndm() * 180. - 90.;
      double y = myR.Rndm() * f1->GetMaximum();
      double yFun = f1->Eval(decSim);
      if (y<=yFun) {
        ++iCr;
        alphaSim[count] = raSim;
        deltaSim[count] = decSim;
        enSim[count] = eSim;
        ++count;
      } //end if
    } //end while
  } // end energy loop
  f->Close();
}

int main(int argc, char* argv[])
{
  //-----------------------------//
  //------ arguments ------------//
  //-----------------------------//
  int nSky; 
  if (argc==1) 
    nSky=nsimBlPen;// use default number of isotropic skies
  else if (argc==2) 
    nSky = atoi(argv[1]); //user defined number of isotropic skies
  else {
    cout << " Synopsys :" << argv[0]
         << "[optional] <number of isotropic skies>"
      << endl;
    return 0;
  }

  // penalized pvalue
  double pvalue = 0.;

  // get local p-value for the penalization
  double minval = 0.;
  double ramin = 0.;
  double decmin = 0.;
  double nexpmin = 0.;
  double nobsmin = 0.;
  double emin = 0.;
  double radiusmin = 0.;

  ifstream file;
  string title = "blindsearch_localp";
  file.open(title.c_str());

  while (!file.eof()) {
    file >> minval >> ramin >> decmin >> nexpmin >> nobsmin >> emin >> radiusmin;
    if (!file.good())
      break;
  }

  cout << "\n";
  cout << " Get Global p-value for blind searches analysis " << "\n";
  cout << " Local p-value " << "\n";
  cout << " Binomial probability " << minval
    << " Ra/deg " << ramin << " Dec/deg " << decmin
    << " Nexp " << nexpmin << " Nobs " << nobsmin
    << " Eth/EeV " << emin << " angular scale [deg] " << radiusmin << endl;
  cout << "\n";

  // read isotropic expectation
  string pv = "phisto.root";
  TFile fv(pv.c_str(),"READ");
  TH2F *hFinal;
  fv.GetObject("hFinal",hFinal);

  // xaxis
  double scale[nbinbl+1];
  SetAngScaleBinningBl(scale);

  // get ipix
  THealpixMap map(nSide,'Q');
  vector<int> ipix = GetIpix(nSide);
  int binz = ipix.size() - 1;

  // yaxis
  double enrange[nEbin];
  SetEbin(enrange);

  // get data
  vector <double> vId,vZen,vTheta,alphaEv,deltaEv,enEv,vexp;
  string pathFile= "../../Data/"+dataname;
  LoadAugerData(pathFile , vId, vZen, alphaEv,deltaEv,enEv, vexp,"eq");
  int nEvents = enEv.size();

  // Get number of events for each energy threshold
  double nOut[nEbin];
  GetNETh(enEv, nOut);

  // loop over the number of simulations
  for (int kk = 0; kk < nSky; ++ kk) {
    ProgressBar((double)kk / (double)nSky);
    double minimum = 10e16;

    // simulate one isotropic sky
    double alphaSim[nEvents];
    double deltaSim[nEvents];
    double enSim[nEvents];
    simSky(alphaSim, deltaSim, enSim);

    // pixels loop
    for (int i = 0; i < binz; ++i) {

      TH2F h2("h2","h2",nbinbl, scale, nEbin, enrange);
      TH2F h2f("h2f","h2f",nbinbl, scale, nEbin, enrange);
      double ra = 0.;
      double dec = 0.;
      map.GiveLB(ipix.at(i), ra, dec);

      // events loop
      for (int k = 0; k<nEvents; ++k) {
        double gamma = 0.;
        gamma = GetAngDist(ra, dec, alphaSim[k], deltaSim[k]);
        if (enSim[k] <= 80)
          h2.Fill(gamma, enSim[k]);
        else
          h2.Fill(gamma, 80);
      } // end event loop

      // integrate the th2 histo
      h2f = GetNormIntegral(h2, 1);
      // angular scale loop
      for (int z = 1 ; z <= nbinbl; ++z) {
        int bin = hFinal->FindBin(scale[z-1], dec);
        double p = hFinal->GetBinContent(bin);

        // energy loop
        for (int j = 0 ; j < nEbin; ++j) {
          double nex  = p*nOut[j];
          double nobs = h2f.GetBinContent(z, j+1);
          double prob = GetBinProb(nex, nobs, nOut[j]);
          if (prob <= minimum) {
            minimum = prob;
          } // end if
        } // end energy loop
      } // end angular scale loop
      h2f.Reset();
      h2.Reset();
    } // end
    cout << "\n";
    cout << " sky " << kk << " p-value " << minimum << endl;

    if (minimum<=minval)
      ++pvalue;
  }

  if (!pvalue) {
    pvalue = 1 / nSky;
    cout <<"\n";
    cout << " Penalized p-value <= " << pvalue << endl;
    cout << " More simulations needed!" << endl;
    cout <<"\n";
  }   
  else {
    pvalue /= nSky;
    cout <<"\n";
    cout << " Penalized p-value " << pvalue << endl;
    cout <<"\n";
  }  

  return 0;
}
