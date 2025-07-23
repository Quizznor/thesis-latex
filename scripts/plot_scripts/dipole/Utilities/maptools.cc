#include <iostream>
#include <fstream>

#include "harmotools.h"
#include "maptools.h"
#include "common.h"
#include "STClibrary.h"

#ifdef gcc323
char* operator+( std::streampos&, char* );
#endif

using namespace std;



THealpixMap map_events(unsigned int nside, const vector<double> & l, const vector<double> & b)
{
  THealpixMap map(nside);
  vector<long> ipix = map.Ip(l,b);
  for(unsigned int i = 0; i < l.size(); i++) map[ipix[i]]++;
  return(map);
}

THealpixMap map_sources(unsigned int nside, const vector<double> & l, const vector<double> & b,const vector<double> &w)
{
  THealpixMap map(nside);
  vector<long> ipix = map.Ip(l,b);
  for(unsigned int i = 0; i < l.size(); i++) map[ipix[i]]+=w[i];
  return(map);
}
vector<long> GetPixFromLB(const THealpixMap & map, double l, double b, double radius)
{
  vector<long> pixKeep;
  vector<long> allPix(map.NPix());
  for(unsigned int i = 0; i < map.NPix(); i++) allPix[i] = i;
  vector<double> lAllPix, bAllPix;
  map.GiveLB(allPix,lAllPix,bAllPix);

  vector<double> uv0 = ll2uv(l,b), uvPix;
  double ang;
  for(unsigned int i = 0; i < map.NPix(); i++)
    {
      uvPix = ll2uv(lAllPix[i],bAllPix[i]);
      ang = acos(uvPix[0]*uv0[0]+uvPix[1]*uv0[1]+uvPix[2]*uv0[2]);
      if( ang <= radius*DTOR ) pixKeep.push_back(allPix[i]);
    }
  return pixKeep;
}



vector<double> ll2uv(double l, double b)
{
  double lon = l*DTOR;
  double lat = b*DTOR;
  double cos_lat = cos(lat);
  double sin_lat = sin(lat);
  double cos_lon = cos(lon);
  double sin_lon = sin(lon);
  vector<double> uv(3);
  uv[0] = cos_lat * cos_lon;
  uv[1] = cos_lat * sin_lon;
  uv[2] = sin_lat;
  return uv;
}



double AngularDistance(double lSky1, double bSky1, double lSky2, double bSky2)
{
  vector<double> uv1 = ll2uv(lSky1, bSky1);
  vector<double> uv2 = ll2uv(lSky2, bSky2);
  double cosAngle, angle;
  cosAngle = uv1[0]*uv2[0]+uv1[1]*uv2[1]+uv1[2]*uv2[2];
  angle = acos(cosAngle)*1./DTOR;

  return angle;
}



double TP_AngularDistance(double theta1, double phi1, double theta2, double phi2)
{
  return ( (1./DTOR)*acos( sin(theta1*DTOR)*cos(phi1*DTOR)*sin(theta2*DTOR)*cos(phi2*DTOR)
			   + sin(theta1*DTOR)*sin(phi1*DTOR)*sin(theta2*DTOR)*sin(phi2*DTOR)
			   + cos(theta1*DTOR)*cos(theta2*DTOR) ) );
}



bool XYtoAngLambertAzimuthal(double x, double y, double & theta, double & phi, double phi0, double cl, double sl)
{
  double rho, c, sc, cc, val;
  bool status = true;
  double thex2=x*x;
  rho=sqrt(thex2+y*y);
  c=2*asin(rho/2);
  cc=cos(c);
  sc=sin(c);
  phi=mod(phi0+atan2(-x*sc,rho*cl*cc-y*sl*sc),TwoPi);
  val=cc*sl+y*sc*cl/rho;
  if( fabs(val)<=1 ) theta=asin(val);
  else status=false;
  return status;
}



bool AngtoXYLambertAzimuthal(double theta, double phi, double & xrad, double & yrad, double lambda0, double cp1, double sp1)
{  
  double cp=cos(theta), sp=sin(theta);
  double dphi=phi-lambda0, sdphi=sin(dphi), cdphi=cos(dphi);
  double k=sqrt(2./(1+sp1*sp+cp1*cp*cdphi));  
  xrad=-k*cp*sdphi;
  yrad=k*(cp1*sp-sp1*cp*cdphi);
  return true; // the transormation is always possible
}



bool XYtoAngMollweide(double x, double y, double & theta, double & phi, double, double, double)
{
  static double prevy = y; // initialization
  static double factor = TwoPi/sin(acos(2.*y));
  if( y != prevy )// recompute the factor
    {
      factor = TwoPi/sin(acos(2.*y));
      prevy=y;
    }
  theta = (0.5-y)*M_PI;
  phi = -x*factor;
  bool status;
  if( phi <= M_PI && phi >=( -M_PI) ) status = true;
  else status = false;
  return status;
}



bool AngtoXYMollweide(double theta, double phi, double & x, double & y, double, double, double)
{
  double factor = sin(acos(theta/PiOver2));
  if( phi>M_PI ) phi -= TwoPi;
  x = -phi*factor;
  y = theta;
  return true;
}



void xyz2thetaphi(double x, double y, double z, double & theta, double & phi)
{
  theta = 0., phi = 0.;
  double norm = sqrt(x*x+y*y+z*z);
  if( norm != 0 )
    {
      theta = acos(z/norm);
      if( x != 0. ) phi = atan2(y,x);
      if( x == 0. ) phi = M_PI/2.;
    }
  if( phi < 0.) phi = phi+TwoPi;
}

