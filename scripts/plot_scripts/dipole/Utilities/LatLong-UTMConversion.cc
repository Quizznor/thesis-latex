#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "STCconstants.h"
#include "LatLong-UTMConversion.h"

/* 
   Converts lat/long to UTM coords. Equations from USGS Bulletin 1532. East Longitudes are positive, 
   West longitudes are negative.North latitudes are positive, South latitudes are negative Lat and
  Long are in decimal degrees.
*/
void LLtoUTM(int ReferenceEllipsoid, double Lat, double Long,
             double &UTMNorthing, double &UTMEasting, string& UTMZone) {

  double a = kEllipsoid[ReferenceEllipsoid].EquatorialRadius;
  double eccSquared = kEllipsoid[ReferenceEllipsoid].eccentricitySquared;
  double k0 = 0.9996;

  double LongOrigin;
  double eccPrimeSquared;
  double N, T, C, A, M;

  // Make sure the longitude is between -180.00 .. 179.9
  double LongTemp = (Long+180)-int((Long+180)/360)*360-180; // -180.00 .. 179.9;

  double LatRad = Lat*kSTC::DTOR;
  double LongRad = LongTemp*kSTC::DTOR;
  double cLatRad = cos(LatRad);
  double sLatRad = sin(LatRad);
  double tLatRad = sLatRad/cLatRad;
  double LongOriginRad;
  int    ZoneNumber;

  ZoneNumber = int((LongTemp + 180)/6) + 1;

  if( Lat >= 56.0 && Lat < 64.0 && LongTemp >= 3.0 && LongTemp < 12.0 )
    ZoneNumber = 32;

  // Special zones for Svalbard
  if( Lat >= 72.0 && Lat < 84.0 ) {
    if(      LongTemp >= 0.0  && LongTemp <  9.0 )
      ZoneNumber = 31;
    else if( LongTemp >= 9.0  && LongTemp < 21.0 )
      ZoneNumber = 33;
    else if( LongTemp >= 21.0 && LongTemp < 33.0 )
      ZoneNumber = 35;
    else if( LongTemp >= 33.0 && LongTemp < 42.0 )
      ZoneNumber = 37;
  }
  LongOrigin = (ZoneNumber - 1)*6 - 180 + 3;  // +3 puts origin in middle of zone
  LongOriginRad = LongOrigin * kSTC::DTOR;

  // compute the UTM Zone from the latitude and longitude
  char utmz[4];
  sprintf(utmz, "%d%c", ZoneNumber, UTMLetterDesignator(Lat));
  UTMZone = utmz;

  eccPrimeSquared = (eccSquared)/(1-eccSquared);

  N = a/sqrt(1-eccSquared*sLatRad*sLatRad);
  T = tLatRad*tLatRad;
  C = eccPrimeSquared*cLatRad*cLatRad;
  A = cLatRad*(LongRad-LongOriginRad);

  double e4=eccSquared*eccSquared;
  double e6=e4*eccSquared;
  double e8=e4*e4;

  double c1 = 1-eccSquared/4
              -3*e4/64
              -5*e6/256
              -175*e8/16384;

  double c2 = -(3*eccSquared/8
                +3*e4/32
                +45*e6/1024
                +105*e8/4096);

  double c3 = 15*e4/256
              +45*e6/1024
              +525*e8/16384;

  double c4 = -(35*e6/3072
                +175*e8/12288);

  double c5 = 315*e8/131072;

  M = a*( c1*LatRad + c2*sin(2*LatRad) + c3*sin(4*LatRad) + c4*sin(6*LatRad) + c5*sin(8*LatRad) );

  double A2 = A*A;
  double A3 = A2*A;
  double T2 = T*T;

  UTMEasting = (double)(k0*N*(A+(1-T+C)*A3/6
                              + (5-18*T+T2+72*C-58*eccPrimeSquared)*A2*A3/120)
                        + 500000.0);

  UTMNorthing = (double)(k0*(M+N*tLatRad*
                             (A2/2+(5-T+9*C+4*C*C)*A*A3/24
                              + (61-58*T+T2+600*C-330*eccPrimeSquared)*A3*A3/720)));

  if(Lat < 0)
    UTMNorthing += 10000000.0; // 10000000 meter offset for southern hemisphere

}

char UTMLetterDesignator(double Lat) {
  char LetterDesignator;

  if((84 >= Lat) && (Lat >= 72))
    LetterDesignator = 'X';
  else if( ( Lat >= -80 ) && ( Lat < -32 ) )
    LetterDesignator = (char)((int)floor((Lat+80)/8)+'C');
  else if( ( Lat >= -32 ) && ( Lat < 8 ) )
    LetterDesignator = (char)((int)floor((Lat+32)/8)+'J');
  else if( ( Lat >= 8 ) && ( Lat < 72 ) )
    LetterDesignator = (char)((int)floor((Lat-8)/8)+'P');
  else if( ( Lat >= 72 ) && ( Lat <= 84 ) )
    LetterDesignator = 'X';
  else
    LetterDesignator = 'Z';

  return LetterDesignator;
}


void UTMtoLL(int ReferenceEllipsoid, double UTMNorthing, double UTMEasting, const string& UTMZone,
             double& Lat,  double& Long ) {
  double k0 = 0.9996;
  double a = kEllipsoid[ReferenceEllipsoid].EquatorialRadius;
  double eccSquared = kEllipsoid[ReferenceEllipsoid].eccentricitySquared;
  double eccPrimeSquared;
  double e1 = (1-sqrt(1-eccSquared))/(1+sqrt(1-eccSquared));
  double N1, T1, C1, R1, D, M;
  double LongOrigin;
  double mu, phi1, phi1Rad;
  double x, y;
  int ZoneNumber;
  char* ZoneLetter;
  int NorthernHemisphere; // 1 for northern hemispher, 0 for southern

  x = UTMEasting - 500000.0; // remove 500,000 meter offset for longitude
  y = UTMNorthing;

  ZoneNumber = strtoul(UTMZone.c_str(), &ZoneLetter, 10);
  if((*ZoneLetter - 'N') >= 0)
    NorthernHemisphere = 1; // point is in northern hemisphere
  else {
    NorthernHemisphere = 0; // point is in southern hemisphere
    y -= 10000000.0; // remove 10,000,000 meter offset used for southern hemisphere
  }

  LongOrigin = (ZoneNumber - 1)*6 - 180 + 3; // +3 puts origin in middle of zone

  eccPrimeSquared = (eccSquared)/(1-eccSquared);
  double e4=eccSquared*eccSquared;
  double e6=e4*eccSquared;
  double e8=e4*e4;
  double e1_2 = e1*e1;

  M = y / k0;
  mu = M/(a*(1-eccSquared/4-3*e4/64-5*e6/256
             -175*e8/16384));

  phi1Rad = mu	+ (3*e1/2-27*e1_2*e1/32)*sin(2*mu)
            + (21*e1_2/16-55*e1_2*e1_2/32)*sin(4*mu)
            +(151*e1*e1_2/96)*sin(6*mu)+(1097*e1_2*e1_2/512)*sin(8*mu);
  phi1 = phi1Rad*kSTC::RTOD;

  double cphi1 = cos(phi1Rad);
  double sphi1 = sin(phi1Rad);
  double tphi1 = sphi1/cphi1;

  N1 = a/sqrt(1-eccSquared*sphi1*sphi1);
  T1 = tphi1*tphi1;
  C1 = eccPrimeSquared*cphi1*cphi1;
  R1 = a*(1-eccSquared)/pow(1-eccSquared*sphi1*sphi1, 1.5);
  D = x/(N1*k0);
  double D2 = D*D;
  double D4 = D2*D2;
  double C1_2 = C1*C1;
  double T1_2 = T1*T1;

  Lat = phi1Rad - (N1*tphi1/R1)*
        (D2/2-(5+3*T1+10*C1-4*C1_2-9*eccPrimeSquared)*D4/24
         +(61+90*T1+298*C1+45*T1_2-252*eccPrimeSquared-3*C1_2)*D2*D4/720);
  Lat = Lat * kSTC::RTOD;

  Long = (D-(1+2*T1+C1)*D*D2/6+(5-2*C1+28*T1-3*C1_2+8*eccPrimeSquared+24*T1_2)*D*D4/120)/cphi1;

  Long = LongOrigin + Long * kSTC::RTOD;

}

void ne2xy(const TObservatory& obs, double n, double e, double *x, double *y)
{
  /// cf GAP-2001-038
  *x = (1+obs.GetBeta())*(e-obs.GetRefEasting())+obs.GetAlpha()*(n-obs.GetRefNorthing());
  *y = (1+obs.GetBeta())*(n-obs.GetRefNorthing())-obs.GetAlpha()*(e-obs.GetRefEasting());
}

void nea2xyz(const TObservatory& obs, double n, double e, double a, double *x, double *y, double *z)
{
  /// cf GAP-2001-038
  ne2xy(obs,n,e,x,y);
  *z = a-obs.GetRefAltitude()-(pow(n-obs.GetRefNorthing(),2)+pow(e-obs.GetRefEasting(),2))*obs.GetGamma();
}

void xy2ne(const TObservatory& obs, double x, double y, double *n, double *e)
{
  double upb = 1 + obs.GetBeta();
  double coeff = obs.GetAlpha()*obs.GetAlpha()+upb*upb;
  *n = obs.GetRefNorthing()+(obs.GetAlpha()*x+upb*y)/coeff;
  *e = obs.GetRefEasting()+(upb*x-obs.GetAlpha()*y)/coeff;
}

void xy2ne(const TObservatory& obs, double x, double y, double dx, double dy, double *n, double *e,
	   double *dn, double *de)
{
  double upb = 1+obs.GetBeta();
  double coeff = obs.GetAlpha()*obs.GetAlpha()+upb*upb;
  *n = obs.GetRefNorthing()+(obs.GetAlpha()*x+upb*y)/coeff;
  *e = obs.GetRefEasting()+(upb*x-obs.GetAlpha()*y)/coeff;
  *dn = (obs.GetAlpha()*fabs(dx)+upb*fabs(dy))/coeff;
  *de = (upb*fabs(dx)+obs.GetAlpha()*fabs(dy))/coeff;
}

