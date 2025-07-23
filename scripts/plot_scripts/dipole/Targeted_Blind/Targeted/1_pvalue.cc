/*
 \file 1_pvalue.cc
 \author Ugo Giaccari (u.giaccari@astro.ru.nl)
 \author Claudio Galelli
 \author Lorenzo Caccianiga
 \date 20 June 2021
 \last modification 4 May 2022
 */

#include "../utils.h"

// Calculate expected number of events from the isotropic expectation (hMean)
// and the p-value (hProb)
void GetPval(string type, int nsky, TH2F& hMean, TH2F& hProb, TH2F hData, TH1F *histSim[])
{
  // set object coordinates: CenA, GC, GalPlane, SGalPlane
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

  // Get data
  vector <double> vId,vZen,vTheta,alphaEv,deltaEv,enEv,vexp;
  string pathFile= "../../Data/"+dataname;
  if (!type.compare("ca"))
    LoadAugerData(pathFile , vId, vZen, alphaEv,deltaEv,enEv, vexp,"eq");
  else if (!type.compare("gc") || !type.compare("gp"))
    LoadAugerData(pathFile , vId, vZen, alphaEv,deltaEv,enEv, vexp,"gal");
  else if (!type.compare("sgp"))
    LoadAugerData(pathFile , vId, vZen, alphaEv,deltaEv,enEv, vexp,"sgal");
  int nEvents = enEv.size();

  int nHisto = nEbin*nbin;
  for (int z=0; z<nHisto; ++z) {
    char *histname = new char[10];
    sprintf(histname, "h%d",z);
    histSim[z] = new TH1F(histname,histname,nEvents,0,nEvents);
    delete [] histname;
  }

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

  // make simulations!

  // loop over isotropic skies
  for (int i = 1; i <= nsky; ++i) {
    ProgressBar((double)i/(double)nsky);
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
          hSim.Fill(gamma, eSim);
          //--------------------//
        } //end if
      } //end while
    } // end energy loop

    // for the other analyses
    hSimTmp2  = GetNormIntegral(hSimTmp, 1);

    // compare single sky map with event map (for p-value calculation)
    for (int k=1; k<=nEbin; ++k) {
      for (int y=1; y<=nbin; ++y) {
        double nData = hData.GetBinContent(y,k);
        double nSim  = hSimTmp2.GetBinContent(y,k);
        int nhisto = GetIndex(y,k);
        histSim[nhisto]->Fill(nSim);
        if (nSim>=nData)
          hProb.Fill(xbins[y-1],k+31, 1 / (double)nsky);
      } // end angular scale loop
    } // end energy loop
    hSimTmp.Reset();
    hSimTmp2.Reset();
  } // end events loop
  delete f1;

  // now integrate the TH2 histo
  hMean  = GetNormIntegral(hSim, nsky);
  hMean.SetName("hMean");

  // print results only for not CenA analyses
  if (type.compare("ca") != 0)
    GetMinimum(hData, hMean, hProb);
  else {
    GetBinProbCen(hData, hMean, hProb);
  }     
}

int main(int argc, char* argv[])
{
  //-----------------------------//
  //------ arguments ------------//
  //-----------------------------//
  int nsky;
  string type;
  if (argc==2) {
    type = argv[1]; // choose the analysis, possible options: ca, gc, gp, sgp
    if (!type.compare("ca"))
      nsky = nsimCA; // default number of isotropic skies for CenA
    else
      nsky = nsimT; // default number of isotropic skies
  }
  else if (argc==3) {
    type = argv[1]; // choose the analysis, possible options: ca, gc, gp, sgp
    nsky = atoi(argv[2]); // number of isotropic skies provided by user
  }
  else {
    cout << " Synopsys :" << argv[0]
         << " <Type of analysis: ca, gc, gp, sgp> "
         << " [optional] <number of isotropic skies>"
      << endl;
    return 0;
  }

  if (!CheckAnalysis(type))
    return 0;
  PrintInfo(type, nsky);
  CheckExpo(); // check exposure file exists

  // Set angular scale bin
  double xbins[nbin+1];
  SetAngScaleBinning(xbins);

  //-------------------------------------------------------------------------------------//
  // observed events as a function of energy and angular scale                           //
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
  int nhistoTot = nEbin*nbin;
  TH1F *histSim[nhistoTot];
  GetPval(type, nsky, hMean, hProb, hData, histSim);
  
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

  return 0;
}
