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
double GetPenalization(string type, int nsky, int flag)
{
  double *pvalue = new double[nEbin]; // penalized p-value(s)

  // get data
  vector <double> vId,vZen,vTheta,alphaEv,deltaEv,enEv,vexp;
  string pathFile= "../../Data/"+dataname;
  if (!type.compare("ca")) 
    LoadAugerData(pathFile , vId, vZen, alphaEv,deltaEv,enEv, vexp,"eq");
  else if (!type.compare("gc") || !type.compare("gp"))  
    LoadAugerData(pathFile , vId, vZen, alphaEv,deltaEv,enEv, vexp,"gal");
  else if (!type.compare("sgp")) 
    LoadAugerData(pathFile , vId, vZen, alphaEv,deltaEv,enEv, vexp,"sgal");
  int nEvents = enEv.size();

  // ------------------------------//
  // take the values
  string file = SetNameFile(type); 
  TFile *f = new TFile(file.c_str());
  file.clear();
  TH2F *hData = (TH2F*)f->Get("hData");
  TH2F *hMean = (TH2F*)f->Get("hMean");  
  TH2F *hProb = (TH2F*)f->Get("hProb");
  int nHisto = nEbin*nbin;
  TH1F *histSim[nHisto];
  for (int z=0; z<nHisto; ++z) {
    char *histname = new char[10];
    sprintf(histname, "h%d",z);
    histSim[z] = (TH1F*)f->Get(histname);
    delete [] histname;
  } 
  GetMinimum(*hData, *hMean, *hProb); 
  double* min = new double[nEbin];
  if (!flag) { 
    for (int j = 1; j<=nEbin; ++j) {
      min[j-1] = 1.; 
      for (int y=1; y<=nbin; ++y) {
        double a = hProb->GetBinContent(y ,j); 
        if (a<min[j-1]) 
          min[j-1] = a;
      }
    }
  }
  else  min[0] = hProb->GetMinimum();
  // ------------------------------//

  // make the simulations //

  // set object coordinates: CenA, GC
  double alphaObj = 0.;
  double deltaObj = 0.;
  SetObjectCoords(type, alphaObj, deltaObj);

  // set angular scale bin
  double xbins[nbin+1]; 
  SetAngScaleBinning(xbins);

  // TH2 histo(s)
  // temporary histos used during the analysis
  TH2F hSimTmp("hSimTmp","hSimTmp",nbin,xbins,nEbin,minEn,maxEn);
  TH2F hSimTmp2("hSimTmp2","hSimTmp2",nbin,xbins,nEbin,minEn,maxEn);
  TH2F hSim("hSim","hSim",nbin,xbins,nEbin,minEn,maxEn);

  // Get number of events in each energy bin
  double nEnbin[nEbin]; 
  GetNEbin(enEv, nEnbin);

  // get number of events for each threshold
  double nOut[nEbin];
  GetNETh(enEv, nOut);

  // get seed
  unsigned int seed;
  GetSeed(seed);
  TRandom3 myR(seed);

  // read exposure function
  TFile *fexpo = TFile::Open(fileInExpo.c_str(),"READ");
  TF1 *f1;
  fexpo->GetObject("fcos",f1);

  // loop over isotropic skies
  for (int i = 1; i <= nsky; ++i) {
    ProgressBar((double)i/(double)nsky);  
    // energy loop
    for (int j = 0; j<nEbin; ++j) {
      double eSim = double(j + minEn);
      // events loop
      int iCr = 0;
      while (iCr < nEnbin[j]) {
        double raSim = myR.Rndm()*360;
        double decSim = myR.Rndm()*180. - 90.;
        double y = myR.Rndm()*f1->GetMaximum();
        double yFun = f1->Eval(decSim);
        if (y<=yFun) {
          ++iCr;
          // calculate the distance (gamma) from the corresponding object or plane
          double gamma = 0.;
          if (!type.compare("ca")) {
            gamma = GetAngDist(raSim, decSim, alphaObj, deltaObj);
          } else if (!type.compare("gc")) {
            double ll = 0.;
            double bb = 0.;
            radec2gal(raSim/15.0, decSim, &ll, &bb,2000);
            gamma = GetAngDist(ll, bb, alphaObj, deltaObj);
          } else if (!type.compare("gp")) {
            double ll = 0.;
            double bb = 0.;
            radec2gal(raSim/15.0, decSim, &ll, &bb,2000);
            gamma = abs(bb);
          } else if (!type.compare("sgp")) {
            double ll(0.),bb(0.);
            radec2gal(raSim/15.0, decSim, &ll, &bb,2000);
            double sll = 0.;
            double sbb = 0.;
            gal2Sgal(ll, bb, sll, sbb);
            gamma = abs(sbb);
          }
          //--------------------//
          // Fill the histo
          hSimTmp.Fill(gamma, eSim);
          //--------------------//
        } // if
      } // while
    } // energy loop
    hSimTmp2  = GetNormIntegral(hSimTmp, 1);

    double localMin = 1.;

    // compare single sky map with event map
    for (int k=1; k<=nEbin; ++k) {
      double nTot = nOut[k-1];
      for (int y=1; y<=nbin; ++y) {
        double nDataSim  = hSimTmp2.GetBinContent(y,k);
        double prob = 0.;
        // Cenaturus A analysis - binomial probability
        if (!type.compare("ca")) {
          double mean = hMean->GetBinContent(y,k);
          if (nDataSim)
            prob = TMath::BetaIncomplete(mean / nTot, nDataSim, nTot - nDataSim + 1);
          else
            prob = 1.;
        }  
        // others analyses with simulations
        else {
          int nhisto = GetIndex(y,k);
          double entries = histSim[nhisto]->GetEntries();
          double integral = histSim[nhisto]->Integral(nDataSim+1,nEvents);
          prob = integral / entries;
          if (!prob)
            prob = 1. / nsky;
        }
        if (prob < localMin)
          localMin = prob;
      } // end angular scale loop
      if (!flag && localMin <= min[k-1])
        ++pvalue[k-1];
    } // end energy loop
    hSimTmp.Reset();
    hSimTmp2.Reset();
    if(flag && localMin <= min[0])
      ++pvalue[0];
  } // end energy loop
  delete f1;
  delete f;

  if (flag) {
    if (!pvalue[0]) {
      pvalue[0] = 1 / nsky;
      cout <<"\n";
      cout << " Penalized p-value <= " << pvalue[0] << endl;
      cout << " More simulations needed!" << endl;
      cout <<"\n";
    }  
    else {
      pvalue[0] /= nsky;
      cout <<"\n";
      cout << " Penalized p-value " << pvalue[0] << endl;
      cout <<"\n";
    }
  } else {
    ofstream ofile;
    string title = "fig6tmp";
    ofile.open(title.c_str());
    for (int k=0; k<nEbin; ++k) { 
      if (!pvalue[k])
        pvalue[k] = 1 / nsky;
      else
        pvalue[k] /= nsky;
      ofile << k+32 << " " << pvalue[k] << endl;
    }
  } 

  delete [] pvalue;

  return 0;
}

int main(int argc, char* argv[])
{
  //-----------------------------//
  //------ arguments ------------//
  //-----------------------------//
  string type = argv[1]; 
  int nsky;
  int flag;

  if (argc>6 || argc<2) {
    cout << " Synopsys :" << argv[0]
         << " <Type of analysis: ca, gc, gp, sgp> "
         << " [optional] -f <0 for angle-only, any other int for normal penalization>"
         << " [optional] -n <number of isotropic skies>"
      << endl;
    return 0;
  } 
  if (!type.compare("ca")) 
    nsky = nsimCAPen; // default number of isotropic skies for CenA
  else 
    nsky = nsimTPen; // default number of isotropic skies 
  for (int i = 2; i < argc; ++i) {
    if (strstr(argv[i], "-n") == argv[i]) {
      nsky = atoi(argv[i+1]); 
      ++i;
    }
    if (strstr(argv[i], "-f") == argv[i]) {
      flag = atoi(argv[i+1]); 
      ++i;
    }
  }

  if (!CheckAnalysis(type))
    return 0;
  cout << "\n";
  cout << " Getting post-trial p-value for " << "\n";
  cout << "\n";
  PrintInfo(type, nsky);
  if (type.compare("ca") && !flag) {
    cout << "\nWARNING: you asked to penalize an analysis which is not Centaurus A " 
      << "for the scan in angle only and not for Energy." <<endl;
    cout << "This combination is not used for any analysis in the original paper "<<endl;
    cout<< "if you meant to do a full penalization use -f 1 (or any other int !=0) " 
      << "or simply launch the program without the -f option. \n "<<endl;
    return 0;  
  }
  // get penalized p-value
  GetPenalization(type, nsky, flag);

  return 0;
}
