/*
 \file 1_pvalue.cc
 \author Ugo Giaccari (u.giaccari@astro.ru.nl)
 \author Claudio Galelli
 \author Lorenzo Caccianiga
 \date 20 June 2021
 \last modification 4 May 2022
 */

#include "../utils.h"

// Calculate expected average number of events from the isotropic expectation (hMean)
// and the p-value (hProb)
void GetPval(int nsky, TH2F& hMean, TH2F& hProb, TH2F hData, TH1D *histSim[])
{
  // set angular scale bin
  double xbins[nbin+1];
  SetAngScaleBinning(xbins);

  // TH2 histo(s)
  // temporary histos used during the analysis
  TH2F hSimTmp("hSimTmp","hSimTmp",nbin,xbins,nEbin,minEn,maxEn);
  TH2F hSimTmp2("hSimTmp2","hSimTmp2",nbin,xbins,nEbin,minEn,maxEn);
  TH2F hSim("hSim","hSim",nbin,xbins,nEbin,minEn,maxEn);

  int nHisto = nEbin*nbin;
  //set the maximum number of pairs for the histograms as twice the ones find in the real data
  double nPairs = 2 * hData.GetMaximum(); 
  for (int z=0; z<nHisto; ++z) {
    char *histname = new char[100];
    sprintf(histname, "h%d",z);
    histSim[z] = new TH1D(histname, histname, nPairs, 0, nPairs);
    delete [] histname;
  }

  // get data
  vector <double> vId,vZen,vTheta,alphaEv,deltaEv,enEv,vexp;
  string pathFile= "../../Data/"+dataname;
  LoadAugerData(pathFile , vId, vZen, alphaEv,deltaEv,enEv, vexp,"eq");
  int nEvents = enEv.size();

  // get number of events in each energy bin
  double nEnbin[nEbin]; 
  GetNEbin(enEv, nEnbin);

  // get seed
  unsigned int seed;
  GetSeed(seed);
  TRandom3 myR(seed);

  // read exposure function
  TFile *f = TFile::Open(fileInExpo.c_str(),"READ");
  TF1 *f1;
  f->GetObject("fcos",f1);

  // make simulations!

  // loop over isotropic skies
  for (int i = 1; i <= nsky; ++i) {
    ProgressBar((double)i/(double)nsky);

    // isotropic sky map
    double alphaSim[nEvents] = {0.};
    double deltaSim[nEvents] = {0.};
    double enSim[nEvents] = {0.};

    int count = 0;
    for (int j = 0; j<nEbin; ++j) { // energy loop
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

    // fill the histos
    for (int zz=0; zz<nEvents; ++zz) {
      for (int zj=zz+1; zj<nEvents; ++zj) {
        double gamma = 0.;
        gamma = GetAngDist(alphaSim[zz], deltaSim[zz], alphaSim[zj], deltaSim[zj]);
        if (enSim[zz] <= 80) {
          hSimTmp.Fill(gamma, enSim[zz]);
          hSim.Fill(gamma, enSim[zz]);
        }  
        else {
          hSimTmp.Fill(gamma, 80);
          hSim.Fill(gamma, 80);
        }    
      }
    }

    // integratation
    hSimTmp2  = GetNormIntegral(hSimTmp, 1);

    // compare single sky map with event map (for the p-value calculation)
    for (int k=1; k<=nEbin; ++k) {
      for (int y=1; y<=nbin; ++y) {
        double nData = hData.GetBinContent(y,k);
        double nSim  = hSimTmp2.GetBinContent(y,k);
        int nhisto = GetIndex(y,k);
        histSim[nhisto]->Fill(nSim);
        if (nSim>=nData){
          hProb.Fill(xbins[y-1],k+31, 1 / (double)nsky);
        }
      } // end angular scale loop
    } // end energy loop
    hSimTmp.Reset();
    hSimTmp2.Reset();
  } // end events loop
  delete f1;
  delete f;
  // now integrate the TH2 histo
  hMean  = GetNormIntegral(hSim, nsky);
  hMean.SetName("hMean");
  GetMinimum(hData, hMean, hProb);
}

int main(int argc, char* argv[])
{
  //-----------------------------//
  //------ arguments ------------//
  //-----------------------------//
  int nsky;
  if (argc==1) 
    nsky = nsimAC;// default number of isotropic skies
  else if (argc==2) 
    nsky = atoi(argv[1]); // number of isotropic skies provided by user
  else {
    cout << " Synopsys :" << argv[0]
         << " [optional] <number of isotropic skies>"
      << endl;
    return 0;
  }

  string type = "ac";    // autocorrelation analysis
  PrintInfo(type, nsky);
  CheckExpo(); // check exposure file exists

  // Set angular scale bin
  double xbins[nbin+1];
  SetAngScaleBinning(xbins);

  //-------------------------------------------------------------------------------------//
  // get number of observed pairs as a function of energy and angular scale              //
  //-------------------------------------------------------------------------------------//
  TH2F hData = GetTH2data(type);

  //-------------------------------------------------------------------------------------//
  // expected number of events as a function of the energy and angular scale             //
  //-------------------------------------------------------------------------------------//
  TH2F hMean("hMean","hMean",nbin,xbins,nEbin,minEn,maxEn);

  //-------------------------------------------------------------------------------------//
  // get local p-value                                                                   //
  //-------------------------------------------------------------------------------------//
  TH2F hProb("hProb","hProb",nbin,xbins,nEbin,minEn,maxEn);

  // histograms used for the penalization
  int nhistoTot = nEbin*nbin;
  TH1D *histSim[nhistoTot];
  GetPval(nsky, hMean, hProb, hData, histSim);

  //-------------------------------------------------------------------------------------//
  // save results in a root file                                                         //
  //-------------------------------------------------------------------------------------//
  string nameFile = SetNameFile(type);
  TFile OutputFile (nameFile.c_str(),"recreate");
  OutputFile.cd();
  for (int h=0; h<nhistoTot; ++h)
    histSim[h]->Write();
  hData.Write();
  hMean.Write();
  hProb.Write();

  cout << "\n";
  cout << " Now go ahead with the penalization! " << "\n";
  cout << "\n";

  return 0;
}
