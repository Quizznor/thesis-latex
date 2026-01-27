#include <stdio.h>
#include <math.h>
#include <cstring>
#include <string>
#include "DataPath.h"

// ROOT
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"

// HEALPIX
#include "healpixmap.h"
#include "STClibrary.h"

#define PI 3.14159265358979323846

using namespace std;

//----- exposure function * cos(declination) for a given zenith range -----//
Double_t fCos(Double_t *x, Double_t *par)
{
  // parameters: 0 latitude | 1 thetamin | 2 thetamax
  Float_t xx = x[0] * TMath::DegToRad(); //declination
  double ksiM = ( cos(par[2]) - sin(par[0]) * sin(xx) ) / ( cos(par[0]) * cos(xx) );
  double ksim = ( cos(par[1]) - sin(par[0]) * sin(xx) ) / ( cos(par[0]) * cos(xx) );
    
  double aM = 0.;
  double am = 0.;
  if (ksiM<-1)
    aM = PI;
  else if (ksiM>=-1 && ksiM<=1)
    aM = acos(ksiM);
  if (ksim<-1)
    am = PI;
  else if (ksim>=-1 && ksim<=1)
    am = acos(ksim);

  double omega = (cos(par[0]) * cos(xx) * (sin(aM) - sin(am)) + (aM-am) * sin(par[0]) * sin(xx)) * cos(xx);
  return omega;
}

//---- sum of the vertical and inclined exposure ------//
Double_t sum(Double_t *x, Double_t *par)
{
  // parameters: 0 thetamax vertical | 1 latitude | 2 normalization vertical events
  // 3 thetamin | 4 thetamax | 5 normalization inclined events

  Float_t xx = x[0] * TMath::DegToRad(); //declination

  // vertical part
  double alpha = 0.;
  Double_t f = ( cos(par[1]) - (sin(par[0]) * sin(xx)) ) / ( cos(par[0]) * cos(xx) );
  if (f>1)
    alpha = 0.;
  else if (f<-1)
    alpha = PI;
  else
    alpha = acos(f);
  double omegav = ( (cos(par[0]) * cos(xx) * sin(alpha) + alpha * sin(par[0]) * sin(xx)) ) * par[2];

  // inclined part
  double aM = 0.;
  double am = 0.;
  double ksiM = ( cos(par[4]) - sin(par[0]) * sin(xx) ) / ( cos(par[0]) * cos(xx) );
  double ksim = ( cos(par[3]) - sin(par[0]) * sin(xx) ) / ( cos(par[0]) * cos(xx) );
  if (ksiM<-1)
    aM = PI;
  else if (ksiM>=-1 && ksiM<=1)
    aM = acos(ksiM);
  if (ksim<-1)
    am = PI;
  else if(ksim>=-1 && ksim<=1)
    am = acos(ksim);
  double omegah = (cos(par[0]) * cos(xx) * (sin(aM) - sin(am)) + (aM-am) * sin(par[0]) * sin(xx)) * par[5];

  return (omegav + omegah) * cos(xx);
}

int main()
{

string pathFile= "./"+dataname; //Dataname specified in Data/DataPath.h
vector <double> vId,vZen,vTheta,alphaEv,deltaEv,enEv,vexp;
  LoadAugerData(pathFile , vId, vZen, alphaEv,deltaEv,enEv, vexp,"eq");
   
  // number of vertical & inclined events
  double nVertOfficial = 2040.;
  double nHorOfficial   = 595.;
  double nVert = 0.;
  double nHor = 0.;
  for(unsigned int i=0;i<vZen.size();i++) {
	if(vZen[i]<60.)nVert++;
	else nHor++;
  }
  // zenith angle range
  double thetaMin = 60. * TMath::DegToRad();
  double thetaMax = 80. * TMath::DegToRad();

  // Auger latitude site
  double augerLat = -35.2 * TMath::DegToRad();

  cout << "\n";
  cout << "Calculate differential exposure normalize with the number of events" << "\n";
  cout << "Read events from file "<< pathFile <<" found "<<enEv.size()<<" events" << "\n";
  cout << " Vertical events: " << nVert << " - inclined events: " << nHor << endl;
  if(nVert!=nVertOfficial || nHor!=nHorOfficial)   cout << " WARNING! Are you using a dataset different from the one released with the paper? Expected Vertical events: " << nVertOfficial << " - inclined events: " << nHorOfficial << endl;
  cout << "\n";

  //---------------------------------------------//
  // EXPOSURE OF THE AUGER HORIZONTAL EVENTS     //
  //---------------------------------------------//
  TF1 fHas("fHas",fCos,-90.,90.,3);
  fHas.SetParameters(augerLat,thetaMin,thetaMax);

  // calculate the normalization factor
  double intexpohor = fHas.Integral(-90.,90.);
  double normHor = nHor * 1. / (intexpohor*TMath::DegToRad()) * 1. / (2*PI);

  //---------------------------------------------//
  // EXPOSURE OF THE AUGER VERTICAL EVENTS       //
  //---------------------------------------------//
  TF1 f("f",fCos,-90.,90.,3);
  f.SetParameters(augerLat, 0, thetaMin);
  //normalization factor
  double intexpo = f.Integral(-90.,90.);
  double normVert = nVert * 1. / (intexpo*TMath::DegToRad()) * 1. / (2*PI);

  //---------------------------------------//
  // FOR THE STRUCTURES ANALYSIS           //
  //---------------------------------------//
  // total exposure (vertical + inclined)
  TF1 fcos("fcos",sum,-90.,90.,6);
  fcos.SetParameters(augerLat, thetaMin, normVert,thetaMin,thetaMax,normHor);

  // save declination distribution in a root file //
  string fileOutPut = "exposure.root";
  TFile fOutPut(fileOutPut.c_str(),"recreate");
  fcos.Write();
  fOutPut.Close();

  //---------------------------------------//
  // FOR THE LIKELIHOOD ANALYSIS           //
  //---------------------------------------//

  //--------------------------------------------------//
  // create HealpixMap of the sky                     //
  // verticals and horizontals			              //
  //--------------------------------------------------//
  int nside=64;
  THealpixMap covV(nside, 'G');
  THealpixMap covH(nside, 'G');
  THealpixMap covTot(nside, 'G');

  double th, ph, ra, dec; 
  for (unsigned int i = 0; i < covV.NPix(); ++i) {
    covV.GiveLB(i, th, ph);
    gal2radec(th, ph, &ra,  &dec);
    covV[i]=f.Eval(dec)/TMath::Cos(TMath::DegToRad()*dec);
    covH[i]=fHas.Eval(dec)/TMath::Cos(TMath::DegToRad()*dec);
  }
  covTot = covV/covV.Total()*nVert + covH/covH.Total()*nHor;

  //--------------------------------------------------//
  // save Healpix maps in fits files                  //
  //--------------------------------------------------//
  string fitfile = "exposure.fits";
  covTot.WriteFits((char*)fitfile.c_str());

  cout << " ....done!" << "\n";
  cout << "\n";

  return 0;
}
