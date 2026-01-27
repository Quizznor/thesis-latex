#ifndef _HEALPIXMAP_H
#define _HEALPIXMAP_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <complex>
#include <sys/time.h>

#include "common.h"

extern "C" 
{
#include "chealpix.h"
#include "fitsio.h"
}

using namespace std;


/*! \file healpixmap.h
  One can find in this file all the functions related to map manipulation. It has been written with help from 
  Sebastien Metrot and various pieces of code have been grabbed from Alexandre Amblard and Dominique Yvon.
*/


/*!
  This class provides an interface for C++ users to the \b HEALPix package widely used in Cosmic Microwave 
  Background. It provides the main features of \b HEALPixsuch as RING and NEST pixel numbering schemes, Map2Cl, 
  Map2Alm, Alm2Map, Cl2Map and Map filtering. It also gives the relation between directions on the sky and pixel 
  numbers. The THealpixMap class inherits from a vector<double> so that pixel values may be accessed through 
  THealpixMap[i].
*/

//! Related to \b HEALPix.
class THealpixMap : public vector<double>
{
 public:
  //! Constructor using nside. It allocates \f$ 12\times nside^2 \f$ pixels.
  THealpixMap(unsigned int nside = 0, char coordsys = 'G');

  //! Constructor from a vector<double>.
  THealpixMap(const vector<double> & v, char coordsys = 'G');  

  //! Copy Constructor.
  THealpixMap(const THealpixMap & map);

  //! Constructor from a fits FITS containing a \b HEALPix map.
  THealpixMap(char * filename);
  
  //! Destructor.
  virtual ~THealpixMap();   
 
  //! Returns the sum of all pixels in the map.
  double Total() const;

  //! Returns the average of all pixels in the map.
  double Mean() const;

  //! Returns the RMS of the map.
  double RMS() const;

  //! Returns the minimum value of the map.
  double Min() const;

  //! Returns the maximum value of the map.
  double Max() const; 

  //! Returns the pixel index and the corresponding \f$ \theta \f$ and \f$ \phi \f$ of the maximum of the map.
  void GetMaxPosition(double & theta, double & phi, unsigned int & ipix) const;

  //! Writes the map into a FITS file.
  int WriteFits(char * filename) const;

  //! Set #_nside, #_npix and #fIpix.
  void SetNSide(unsigned int nside);

  //! Returns the size of the pixel according to #_npix.
  double GetPixSize() const;

  //! Calls ud_grade executable to upgrade/degrade the map.
  THealpixMap Map2Map(int nside_out);

  //! Calls anafast executable to compute the \f$ C_\ell \f$ power spectrum of the map up to lmax multipole.
  vector<double> Map2Cl(int lmax) const;

  /*! 
    Calls synfast executable to fill the map with a realisation of the power spectrum given in variable 
    \f$ C_\ell \f$ convolved with the beam given by beam. The phases are randomly drawn by anafast.
  */
  void Cl2Map(vector<double> & cl, float beam = 0.);

  //! Calls anafast executable to compute the \f$ a_{\ell m} \f$ coefficients of the map up to lmax multipole.
  vector<vector<complex<double> > > Map2Alm(int lmax); 

  //! Calls synfast executable to transform \f$ a_{\ell m} \f$ coefficients into a \b HEALPix map.
  void Alm2Map(vector<vector<complex<double> > > & alm, float beam = 0.);

  /*! 
    Filters the map in harmonic space with the beam given by variables theta and lobe such as the beam is 
    lobe(\f$ \theta \f$). The map is filtered with the legendre transform of the lobe.
  */
  THealpixMap Filter(const vector<double> & theta, const vector<double> & lobe); 

  /*! 
    Integrates the map through the beam. Uses the Filter method and then normalizes with the effective area of 
    the beam (that is renormalized to a maximum of 1).
  */
  THealpixMap IntBeam(const vector<double> & theta, const vector<double> & lobe); 

  /*! 
    Filters the map in harmonic space with the beam transfer function \f$ B_\ell \f$ which should be the legendre 
    transform of the real space beam to be applied to the map.
  */
  THealpixMap FiltBl(const vector<double>& bl);

  /*! 
    Compute the new map in a new coordinate system. The conversions G->Q and reciprocal are allowed.
  */
  THealpixMap ChCoordSys(char newcoordsys) const;

  //! Returns the pixel indices of the directions given by l and b.
  vector<long> Ip(const vector<double> & l, const vector<double> & b) const;   

  //! Returns the pixel indice of the direction given by l and b.
  long Ip(double l, double b) const; 

  //! Returns the pixels values of the directions given by l and b.
  vector<double> Values(const vector<double>& l, const vector<double>& b) const; 

  //! Returns the pixel value of the direction given by l and b.
  double Value(double l, double b) const;

  //! Returns the l and b direction of the pixel indice given in ipix.
  void GiveLB(long ipix, double& l, double& b) const;

  //! Returns the l and b directions of the pixels indices given in ipix.
  void GiveLB(const vector<long>& ipix, vector<double>& l, vector<double>& b) const;

  /*! 
    Returns a 2 dimensional map (nx,ny) of the Mollweide projection of the map. The background of the map 
    (not on the sphere) is set to the value bg and the image can be written into the FITS file filename 
    that may be read and displayed by an external viewing program.
  */
  vector<vector<double> > Map2MollProj(int sizeX, int sizeY, double background = 0., double decLimit = 90., char * filename = (char*)" ") const;

  /*! 
    Returns a 2 dimensional map (nx,ny) of the Lambert azimuthal projection of the map around the location 
    specified by lcenter, bcenter with the radius "radius". The image can be written into the FITS file 
    filename that may be read and displayed by an external viewing program. The formulae used are from 
    <a href="http://mathworld.wolfram.com/LambertAzimuthalEqual-AreaProjection.html">WolframMathWorld</a> web site
  */
  vector<vector<double> > Map2LambertAzimuthalProj(int sizeX, int sizeY, double longRef, double latRef, double radius, double background, double decLimit = 90., char * filename = (char*)" ") const;

  /*! 
    Returns the local maxima of the \b HEALPix map. Values are returned and pixel are given in ipmax. This 
    routines spawns the "hotspots_cxx" executable provided in the \b HEALPix release. The resulting vectors 
    are sorted in decreasing map value order.
  */
  vector<double> FindMaxima(vector<long> & ipmax) const;

  //! Returns #ordering.
  const char * Ordering() const;

  //! Returns #coordsys.
  char CoordSys() const;

  //! Returns #units.
  const char * Units() const;

  //! Set #ordering.
  void SetOrdering(const char* pOrdering);

  //! Set #coordsys.
  void SetCoordSys(char coord);

  //! Returns #_nside
  unsigned int NSide() const; 

  //! Returns #_npix 
  unsigned int NPix() const; 

  //! Surdefinition of the operator + between two maps. 
  THealpixMap operator+(THealpixMap map2) const;
  
  //! Surdefinition of the operator - between two maps. 
  THealpixMap operator-(THealpixMap map2) const;

  //! Surdefinition of the operator * between two maps. 
  THealpixMap operator*(THealpixMap map2) const;

  //! Surdefinition of the operator / between two maps. 
  THealpixMap operator/(THealpixMap map2) const;

  //! Surdefinition of the operator + between a map and a double. 
  THealpixMap operator+(double x) const;

  //! Surdefinition of the operator - between a map and a double. 
  THealpixMap operator-(double x) const;

  //! Surdefinition of the operator * between a map and a double. 
  THealpixMap operator*(double x) const;

  //! Surdefinition of the operator / between a map and a double. 
  THealpixMap operator/(double x) const;

  //! Redefine operator = (copy)
  THealpixMap & operator=(const THealpixMap &);

  //! Redefine =
  THealpixMap & operator = (double x);

  //! Redefine *=
  THealpixMap & operator *= (double x);

  //! Redefine /=
  THealpixMap & operator /= (double x);
  
  //! Redefine +=
  THealpixMap & operator += (double x);

  //! Redefine -=
  THealpixMap & operator -= (double x);

  //! Pixel indexes of the map
  vector<long> fIpix;

 private:
  //! Ordering scheme. Could be either RING or NEST.
  char ordering[256];

  //! Coordinate system. Could be either 'G' or 'Q'.
  char coordsys[256];

  //! Units in which the map is expressed.
  char units[256];

 protected:
  //! nside parameter of the map.
  unsigned int _nside;

  //! Number of pixels in the map.
  unsigned int _npix;
};

//! Defines x * THealpixMap = THealpixMap * x
inline THealpixMap operator*(double x, const THealpixMap& map) {return map*x;}

//! Defines x + THealpixMap = THealpixMap + x
inline THealpixMap operator+(double x, const THealpixMap& map) {return map+x;}

//! Defines x - THealpixMap = (-1)*THealpixMap + x
inline THealpixMap operator-(double x, const THealpixMap& map) {return (-1)*map+x;}

#endif
