#ifndef _LATLONGCONV_H_
#define _LATLONGCONV_H_

// Short description of the file
/*! \file LatLong-UTMConversion.h
  @brief Latitude-Longitude to UTM coordinates conversions.
*/

// Detailed description of the file
/*! \file LatLong-UTMConversion.h
  This file gathers the functions you need to convert the latitude-longitude into UTM coordinates. The class
  Ellipsoid, which described the shape of the Earth, is also defined here. 
*/

#include <string>
#include "STCconstants.h"

using namespace std;

//! Converts latitude-longitude into UTM. Lat and Long must be in degrees.
void LLtoUTM(int ReferenceEllipsoid, double Lat, double Long,
             double &UTMNorthing, double &UTMEasting, string& UTMZone);

//! Converts UTM into latitude-longitude. Returns lat and long in degrees
void UTMtoLL(int ReferenceEllipsoid, double UTMNorthing, double UTMEasting, const string& UTMZone,
             double& Lat,  double& Long );

/*! This routine determines the correct UTM letter designator for the given latitude in degrees.\n
  Returns 'Z' if latitude is outside the UTM limits of 84N to 80S written by 
  <a href="mailto:chuck.gantz@globalstar.com">Chuck Gantz</a>. The letters 'I' and 'O' are not used.
*/
char UTMLetterDesignator(double Lat);


// Detailed description of the class
/*!
  Reference ellipsoids derived from <a href="mailto:pdana@mail.utexas.edu">Peter H. Dana</a>\n 
  Department of Geography, University of Texas at Austin\n 
  \n
  Source Defense Mapping Agency. 1987b. DMA Technical Report :\n
  Supplement to Department of Defense World Geodetic System 1984\n
  Technical Report. Part I and II. Washington, DC: Defense Mapping
*/


// Short description of the class
//! Reference ellipsoid used in the description of the shape of the Earth.
class Ellipsoid {
public:
  //! Default constructor
  Ellipsoid() {}

  //! Constructor 
  Ellipsoid(int Id, string name, double radius, double ecc) {
    id = Id;
    ellipsoidName = name;
    EquatorialRadius = radius;
    eccentricitySquared = ecc;
  } 

  //! Index of the reference ellipsoid (default = 23 = WGS-84)
  int id; 
  
  //! Name of the reference ellipsoid 
  string ellipsoidName;

 //! Equatorial radius associated to the ellipsoid 
  double EquatorialRadius;
  
  //! Square of the eccentricity of the ellipsoid
  double eccentricitySquared; 

};

//! UTM Northing-Easting to X-Y 
void ne2xy(const TObservatory& obs, double n, double e, double *x, double *y);

//! UTM Northing-Easting-Altitude to X-Y-Z 
void nea2xyz(const TObservatory& obs, double n, double e, double a, double *x, double *y, double *z);

//! X-Y to UTM Northing-Easting 
void xy2ne(const TObservatory& obs, double x, double y, double *n, double *e);

//! X-Y to UTM Northing-Easting. Manages errors on position
void xy2ne(const TObservatory& obs, double x, double y, double dx, double dy, double *n, double *e,
	   double *dn, double *de);

//! Various ellipsoids of reference commonly used to describe the shape of the earth
const Ellipsoid kEllipsoid[] =
  {
    //  id, Ellipsoid name, Equatorial Radius, square of eccentricity
    Ellipsoid( -1, "Placeholder", 0, 0), //placeholder only, To allow array indices to match id numbers
    Ellipsoid( 1, "Airy", 6377563, 0.00667054),
    Ellipsoid( 2, "Australian National", 6378160, 0.006694542),
    Ellipsoid( 3, "Bessel 1841", 6377397, 0.006674372),
    Ellipsoid( 4, "Bessel 1841 (Nambia) ", 6377484, 0.006674372),
    Ellipsoid( 5, "Clarke 1866", 6378206, 0.006768658),
    Ellipsoid( 6, "Clarke 1880", 6378249, 0.006803511),
    Ellipsoid( 7, "Everest", 6377276, 0.006637847),
    Ellipsoid( 8, "Fischer 1960 (Mercury) ", 6378166, 0.006693422),
    Ellipsoid( 9, "Fischer 1968", 6378150, 0.006693422),
    Ellipsoid( 10, "GRS 1967", 6378160, 0.006694605),
    Ellipsoid( 11, "GRS 1980", 6378137, 0.00669438),
    Ellipsoid( 12, "Helmert 1906", 6378200, 0.006693422),
    Ellipsoid( 13, "Hough", 6378270, 0.00672267),
    Ellipsoid( 14, "International", 6378388, 0.00672267),
    Ellipsoid( 15, "Krassovsky", 6378245, 0.006693422),
    Ellipsoid( 16, "Modified Airy", 6377340, 0.00667054),
    Ellipsoid( 17, "Modified Everest", 6377304, 0.006637847),
    Ellipsoid( 18, "Modified Fischer 1960", 6378155, 0.006693422),
    Ellipsoid( 19, "South American 1969", 6378160, 0.006694542),
    Ellipsoid( 20, "WGS 60", 6378165, 0.006693422),
    Ellipsoid( 21, "WGS 66", 6378145, 0.006694542),
    Ellipsoid( 22, "WGS-72", 6378135, 0.006694318),
    Ellipsoid( 23, "WGS-84", 6378137, 0.00669438)
  };

#endif
