#ifndef _STCCONSTANTS_H_
#define _STCCONSTANTS_H_

// Short description of the file
/*! \file STCconstants.h
  @brief Related to orbital elements
*/



// Detailed description of the file
/*! \file STCconstants.h
  This file enclosed the kSTC namespace, which enclose physical, mathematical and astronomical constants used 
  all along the code and the TObservatory class that describes the observation site.
*/



#include <string>
#include <ctime>

using namespace std;

//! Physical, mathematical and astronomical constants
namespace kSTC
{

#if 0
  // Constants for changes in coordinates
  // UTM - xyz - LatLong part
  const string AugerUTMZone("19H");
  const int  WGS84  = 23;       ///< UTM stuff
  const double Alpha = 0.00252;
  const double Beta  = 0.000603;
  const double Gamma = 0.00000007853;     ///< UTM conversion constants
  /// UTM -> cartesian conversion, easting and northing references (69.25,35.25 cf GAP 2001-038)
  const double E0  = 477256.0;
  const double N0  = 6099203.0;
  const double A0  =   1400.0;
  const double NorthDeclination = 3.53;
#endif

  //! Converts hours (longitude for instance) into degrees
  const double Hrtod = 15.;
  
  //! Converts degrees into hours
  const double Dtohr = 0.066666666666666;
  
  //! Converts radians into hours
  const double Rtohr = 3.8197186342;
  
  //! Converts hours into radians
  const double Hrtor = 0.261799387799;
  
  //! Julian date of 1/1/1970 UTC (origin for unix seconds)
  const double UnixStartJD = 2440587.5;
  
  // Ecliptic and galactic values at J2000. All values are in decimal degrees.
  //! Reference equinox for the computations
  const int RefEquinox = 2000;
  
  //! Ecliptic obliquity  = 23 deg 26' 21''.448 (J2000)
  const double Eps_2000 = 23.4392911111;
  
  //! Right Ascension of galactic North Pole
  const double AlphaG_2000 = 192.85948;
  
  //! Declination of galactic North Pole
  const double DeltaG_2000 = 27.12825;

  //! Galactic latitude of the equatorial South Pole 
  const double GalacticSouthPoleB  = -27.128250; // dec = -90

  //! Galactic longitude of the equatorial South Pole
  const double GalacticSouthPoleL  = 302.93192; // ra = 0

  //! Galactic latitude of the equatorial North Pole 
  const double GalacticNorthPoleB  = 27.128250; // dec = 90

  //! Galactic longitude of the equatorial North Pole
  const double GalacticNorthPoleL  = 122.93192; // ra = 0

  //! Galactic longitude of celestial equator
  const double Lomega_2000 = 32.93192;
  
  //! Ecliptic longitude of galactic North Pole
  const double AlphaE_2000 = 180.02322;
  
  //! Ecliptic latitude of galactic North Pole
  const double DeltaE_2000 = 29.811438523;
  
  //! Galactic longitude of ecliptic equator
  const double Eomega_2000 = 6.3839743;  
  
  //! Value of J2000 in Julian days. A julian year is 365.25 days
  const double J2000 = 2451545.0;    

  //! Speed of light in \f$ m.\mu s^{-1} \f$
  const double CMICRO  =  299.792458;
  
  //! Speed of light in \f$ m.ns^{-1} \f$
  const double CSPEED  =  0.299792458;
  
  //! Speed of light in \f$ m.s^{-1} \f$
  const double CSPEEDSI=  299792458;

  //! Radian to degree conversion
  const double RTOD =  57.2957795130823;

  //! Degree to radian conversion
  const double DTOR =  0.017453292519943295;
  
  //! \f$ 2\pi \f$
  const double TwoPi = 6.28318530717958623200;
  
  //! \f$ \frac{\pi}{2} \f$
  const double PiOver2 = 1.57079632679489661923;
}


//! Characteristics of the observation site. See GAP2001_038 for more explanations.
class TObservatory
{
 public:
  //! Default constructor
  TObservatory(){Init();}

  //! Constructor
  TObservatory(string utm, int ref, double z, double e, double n, double alpha, double beta, double gamma)
    : _UTMZone(utm),_RefEllipsoid(ref),_RefAltitude(z),_RefNorthing(n),_RefEasting(e),_Alpha(alpha),
    _Beta(beta),_Gamma(gamma)
    {
    }
  
  //! Set the UTM zone of the observatory
  void SetUTMZone(string utm) {_UTMZone = utm;}

  //! Set the reference ellipsoid to describe the shape of the Earth
  void SetRefEllipsoid(int id) {_RefEllipsoid = id;}

  //! Set the altitude of the observation site
  void SetRefAltitude(double z) {_RefAltitude = z;}

  //! Set the UTM northing of the observatory
  void SetRefNorthing(double n) {_RefNorthing = n;}

  //! Set the UTM easting of the observatory
  void SetRefEasting(double e) {_RefEasting = e;}

  //! Set the longitude of the observatory
  void SetLongitude(double l) {_Longitude = l;}

  //! Set the latitude of the observatory
  void SetLatitude(double l) {_Latitude = l;}

  //! Set \f$ \alpha \f$.
  void SetAlpha(double x) {_Alpha = x;}

  //! Set \f$ \beta \f$.
  void SetBeta(double x) {_Beta = x;}

  //! Set \f$ \gamma \f$.
  void SetGamma(double x) {_Gamma = x;}

  //! Get the UTM zone
  string GetUTMZone() const {return _UTMZone;}

  //! Get the reference ellipsoid used
  int GetRefEllipsoid() const {return _RefEllipsoid;}

  //! Get the Altitude of the observation site
  double GetRefAltitude() const {return _RefAltitude;}

  //! Get the UTM easting of the observation site
  double GetRefEasting() const {return _RefEasting;}
  
  //! Get the UTM northing of the observation site
  double GetRefNorthing() const {return _RefNorthing;}
  
  // Get the longitude of the observation site
  double GetLongitude() const {return _Longitude;}

  //! Get the latitude of the observation site
  double GetLatitude() const {return _Latitude;}

  //! Get \f$ \alpha \f$
  double GetAlpha() const {return _Alpha;}
  
  //! Get \f$ \beta \f$
  double GetBeta() const {return _Beta;}
  //! Get \f$ \gamma \f$
  double GetGamma() const {return _Gamma;}

  //! Initialize some of the private attributes
  void Init()
  {
    // defaults values
    _RefAltitude = 0;
    // no correction for Earth shape :
    // Coefficients de translation et d'homtétie pour corriger de la courbure de la terre pour la conversion UTM en X et Y
    _Alpha = 0;
    _Beta = 0;
    _Gamma = 0;
  }

 private:
  
  //! UTM zone
  string _UTMZone;

  //! Reference ellipsoid
  int _RefEllipsoid;

  //! Altitude of the observatory
  double _RefAltitude;

  //! UTM northing of the observatory
  double _RefNorthing;

  //! UTM easting of the observatory
  double _RefEasting;

  //! Longitude of the observatory
  double _Longitude;

  //! Latitude of the observatory
  double _Latitude;
  
  //! \f$ \alpha \f$. Used in the UTM coordinates to regional cartesian coordinates conversion.
  double _Alpha;
  
  //! \f$ \beta \f$. Used in the UTM coordinates to regional cartesian coordinates conversion.
  double _Beta;
  
  //! \f$ \gamma \f$. Used in the UTM coordinates to regional cartesian coordinates conversion.
  double _Gamma;
};

#endif
