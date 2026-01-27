/*
 \file 2_penalization.cc
 \author Ugo Giaccari (u.giaccari@astro.ru.nl)
 \author Claudio Galelli
 \author Lorenzo Caccianiga
 \date 20 June 2021
 \last modification 4 May 2022
 */

#include "../utils.h"

// function for penalization
double GetPenAc(string type, int nsky)
{
  double pvalue = 0.;

  // ------------------------------//
  // take info 		   //	
  // ------------------------------//
  string file = SetNameFile(type);
  // take the values
  TFile *f = new TFile(file.c_str());
  file.clear();
  TH2F *hData = (TH2F*)f->Get("hData");
  TH2F *hMean = (TH2F*)f->Get("hMean");
  TH2F *hP = (TH2F*)f->Get("hProb");
  int nHisto = nEbin*nbin;
  TH1F *histSim[nHisto];
  for (int z=0; z<nHisto; ++z) {
    char *histname = new char[10];
    sprintf(histname, "h%d",z);
    histSim[z] = (TH1F*)f->Get(histname);
    delete [] histname;
  }
  GetMinimum(*hData, *hMean, *hP);
  double min = hP->GetMinimum();
  cout << " p-value to penalize " << min << endl;
  // ------------------------------//

  // Set angular scale bin
  double xbins[nbin+1]; //42+1 bins
  SetAngScaleBinning(xbins);

  // TH2 histo(s)
  // temporary histos used during the analysis
  TH2F hSimTmp("hSimTmp","hSimTmp",nbin,xbins,nEbin,minEn,maxEn);
  TH2F hSimTmp2("hSimTmp2","hSimTmp2",nbin,xbins,nEbin,minEn,maxEn);
  TH2F hSim("hSim","hSim",nbin,xbins,nEbin,minEn,maxEn);

  // Get data
  vector <double> vId,vZen,vTheta,alphaEv,deltaEv,enEv,vexp;
  string pathFile= "../../Data/"+dataname;
  LoadAugerData(pathFile , vId, vZen, alphaEv,deltaEv,enEv, vexp,"eq");
  int nEvents = enEv.size(); 

  // Get number of events in each energy bin
  double nEnbin[nEbin]; 
  GetNEbin(enEv, nEnbin);

  // Get number of events for each threshold
  double nOut[nEbin];
  GetNETh(enEv, nOut);

  // Get seed
  unsigned int seed;
  GetSeed(seed);
  TRandom3 myR(seed);

  // read exposure function
  TFile *fexpo = TFile::Open(fileInExpo.c_str(),"READ");
  TF1 *f1;
  fexpo->GetObject("fcos",f1);

  // make the simulations //

  // loop over isotropic skies
  for (int i = 1; i <= nsky; ++i) {
    ProgressBar((double)i/(double)nsky);

    // to store simulations
    double alphaSim[nEvents];
    double deltaSim[nEvents];
    double enSim[nEvents];

    int count=0;
    // energy loop
    for (int j = 0; j<nEbin; ++j) {
      double eSim = double(j + minEn);
      // events loop
      int iCr = 0;
      while (iCr < nEnbin[j]){
        double raSim = myR.Rndm()*360;
        double decSim = myR.Rndm()*180. - 90.;
        double y = myR.Rndm()*f1->GetMaximum();
        double yFun = f1->Eval(decSim);
        if (y<=yFun) {
          ++iCr;
          alphaSim[count] = raSim;
          deltaSim[count] = decSim;
          enSim[count] = eSim;
          ++count;
          //--------------------//
        } // if
      } // while
    } // energy loop

    // fill the histos
    for (int zz=0; zz<nEvents; ++zz) {
      for (int zj=zz+1; zj<nEvents; ++zj) {
        double gamma = 0.;
        gamma = GetAngDist(alphaSim[zz], deltaSim[zz], alphaSim[zj], deltaSim[zj]);
        if (enSim[zz] <= 80) {
          hSimTmp.Fill(gamma, enSim[zz]);
        }  
        else {
          hSimTmp.Fill(gamma, 80);
        }    
      }
    }

    hSimTmp2  = GetNormIntegral(hSimTmp, 1);

    double localMin = 1.;

    // compare single sky map with event map
    for (int k=1; k<=nEbin; ++k) {
      for (int y=1; y<=nbin; ++y) {
        double nDataSim  = hSimTmp2.GetBinContent(y,k);
        double prob = 0.;
        hP = (TH2F*)f->Get("hP");
        int nhisto = GetIndex(y,k);
        double entries = histSim[nhisto]->GetEntries();
        double integral = histSim[nhisto]->Integral(nDataSim+1,nEvents);
        prob = integral / entries;
        if (!prob)
          prob = 1. / nsky;
        if (prob < localMin)
          localMin = prob;
      } // end angular scale loop
    } // end energy loop

    if (localMin <= min)
      ++pvalue;

    hSimTmp.Reset();
    hSimTmp2.Reset();
  } //  loop
  delete f1;
  delete fexpo;
  delete f;

  if (!pvalue) {
    pvalue = 1 / nsky;
    cout << " Penalized p-value <= " << pvalue << endl;
    cout << " More simulations needed!" << endl;
    cout << "\n";
  }  
  else {
    pvalue /= nsky;
    cout << " Penalized p-value " << pvalue << endl;
    cout << "\n";
  }

  return pvalue;
}

int main(int argc, char* argv[])
{

  //-----------------------------//
  //------ arguments ------------//
  //-----------------------------//
  int nsky; 
  string type = "ac";
  if (argc==1) 
    nsky=nsimACPen;
  else if (argc==2) 
    nsky = atoi(argv[1]); // number of isotropic skies
  else {
    cout << " Synopsys :" << argv[0]
         << " [optional] <number of isotropic skies>"
      << endl;
    return 0;
  }

  if (!CheckAnalysis(type))
    return 0;

  cout << "\n";
  cout << " Getting post-trial p-value for " << "\n";
  cout << "\n";
  PrintInfo(type, nsky);
  GetPenAc(type, nsky);

  return 0;
}
