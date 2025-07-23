/*
 \file 1_pvalue.cc
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

int main(int argc, char* argv[])
{
  //-----------------------------//
  //------ arguments ------------//
  //-----------------------------//
  int nev;
  if (argc==1) 
    nev=nsimBl;//default number of isotropic events
  else if (argc==2) 
    nev = atoi(argv[1]); //number of isotropic events specified by user
  else {
    cout << " Synopsys :" << argv[0]
      << "[optional] <Number of isotropic events to simulate>"
      << endl;
    return 0;
  }

  cout << "\n";
  cout << "------------------------------ " << "\n";
  cout << " make isotropic simulations    " << "\n";
  cout << "------------------------------ " << "\n";
  cout << "\n";

  ///////////////////////////////////
  //  Expected number of events    //
  ///////////////////////////////////
  // get ipix
  THealpixMap map(nSide,'Q');
  vector<int> ipix = GetIpix(nSide);
  int binz = ipix.size() - 1;

  // decrange
  vector<double> decrange = SetDecAxis(nSide);
  int ynbin = decrange.size() - 1;

  // xaxis
  double scale[nbinbl+1];
  SetAngScaleBinningBl(scale);

  // yaxis
  double enrange[nEbin];
  SetEbin(enrange);

  CheckExpo(); // check exposure file exists
  // exposure function
  TFile *f = TFile::Open(fileInExpo.c_str(),"READ");
  TF1 *f1;
  f->GetObject("fcos",f1);

  TH2F hFinal("hFinal", "hFinal", nbinbl, scale, ynbin, decrange.data());
  TH2F hSim("hSim", "hSsim", nbinbl, scale, ynbin, decrange.data());

  // Get seed
  unsigned int seed;
  GetSeed(seed);
  TRandom3 myR(seed);

  // events loop
  int iCr = 0;
  while (iCr < nev) {
    ProgressBar((double)iCr/(double)nev);
    double raSim = myR.Rndm() *360.;
    double decSim = myR.Rndm() * 180. - 90.;
    double y = myR.Rndm() * f1->GetMaximum();
    double yfun = f1->Eval(decSim);
    if (y <= yfun) {
      ++iCr;
      // loop on sky positions
      for ( int u = 0; u < ynbin; ++u) {
        double dist = GetAngDist(raSim, decSim, 0, decrange[u]);
        hSim.Fill(dist, decrange[u], 1);
      }
    }
  }

  // integrate TH2F
  for (int i = 1; i <= ynbin; ++i) {
    double content = 0.;
    for ( int u = 1; u <= nbinbl; ++u) {
      content = hSim.Integral(1, u, i, i);
      hFinal.SetBinContent(u, i, content / (double) nev );
    }
  }

  //--------------------------------------------------//
  // save the 2d histo in a root file                 //
  //--------------------------------------------------//
  string fileOutPut = "phisto.root";
  TFile fOutPut(fileOutPut.c_str(),"recreate");
  hFinal.Write();
  fOutPut.Close();

  ///////////////////////////////////
  //  Observed events computation  //
  ///////////////////////////////////

  cout << "\n";
  cout << "------------------------------ " << "\n";
  cout << " get local p-value " << "\n";
  cout << "------------------------------ " << "\n";
  cout << "\n";

  // get data
  vector <double> vId,vZen,vTheta,alphaEv,deltaEv,enEv,vexp;
  string pathFile= "../../Data/"+dataname;
  LoadAugerData(pathFile , vId, vZen, alphaEv,deltaEv,enEv, vexp,"eq");
  int nEvents=enEv.size();

  // Get number of events for each threshold
  double nOut[nEbin];
  GetNETh(enEv, nOut);

  TH2F h2("h2","h2",nbinbl, scale, nEbin, enrange);
  TH2F h2f("h2f","h2f",nbinbl, scale, nEbin, enrange);

  ///////////////////////////////////
  // for Pvalue computation        //
  ///////////////////////////////////
  ofstream myfile;
  string title = "blindsearch_localp";
  myfile.open(title.c_str());

  // minimum p-value
  double minval = 99999.;
  double decmin = 0.;
  double ramin = 0.;
  double nobsmin = 0.;
  double nexpmin = 0.;
  double emin = 0;
  double radiusmin = 0;

  // ipix loop
  for(int i = 0; i < binz; ++i) {
    ProgressBar((double)i/(double)binz);

    // get ra and get of the pixel (centre)
    double ra = 0.;
    double dec = 0.;
    map.GiveLB(ipix.at(i), ra, dec);

    // events loop
    for (int k = 0; k<nEvents; ++k) {
      double gamma = 0.;
      gamma = GetAngDist(ra, dec, alphaEv[k], deltaEv[k]);
      if (enEv[k] <= 80)
        h2.Fill(gamma, enEv[k]);
      else
        h2.Fill(gamma, 80);
    } // end event loop

    // integrate the th2 histo
    h2f  = GetNormIntegral(h2, 1);

    // angular scale loop
    for (int z = 1; z <= nbinbl; ++z) {
      // get isotropic probabilities for each declination and search radius
      int bin = hFinal.FindBin(scale[z-1], dec);
      double p = hFinal.GetBinContent(bin);

      // energy loop
      for (int j =0; j < nEbin; ++j) {
        double nex  = p * nOut[j];
        double nobs = h2f.GetBinContent(z, j+1);
        double prob = GetBinProb(nex, nobs, nOut[j]);

        // minimum p-value
        if (prob<minval) {
          minval = prob;
          decmin = dec;
          ramin = ra;
          nexpmin = nex;
          nobsmin = nobs;
          emin = enrange[j];
          radiusmin = scale[z];
        }
      } // end for
    } // end for
    h2f.Reset();
    h2.Reset();
  } // end angular scale loop

  double ntot = nOut[(int)emin-32]; // total number of event above the Eth with the lowest p-value
  double sign = GetLiMa(nobsmin, (ntot - nobsmin), (nexpmin / ( ntot - nexpmin) ) );

  cout << "\n";
  cout << " ------------------------------------ " << "\n";
  cout << " local p-value " << "\n";
  cout << " binomial probability " << minval << "\n";
  cout << " ra " << ramin << " dec " << decmin << "\n";
  cout << " nexp " << nexpmin << " nobs " << nobsmin << "\n";
  cout << " emin " << emin << " radius " << radiusmin << "\n";
  cout << " number of events above " << emin << " EeV: " << ntot << "\n";
  cout << " Li&Ma significance " << sign << endl;
  cout << " ------------------------------------ " << "\n";
  cout << "\n";

  myfile << minval
    << " " << ramin << " " << decmin
    << " " << nexpmin << " " << nobsmin
    << " " << emin << " " << radiusmin << endl;

  cout << "\n";
  cout << " Now go ahead with the penalization! " << "\n";
  cout << "\n";

  return 0;
}
