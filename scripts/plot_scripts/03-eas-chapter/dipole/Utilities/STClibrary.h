#ifndef _STCLIBRARY_H_
#define _STCLIBRARY_H_


// Begin the main page documentation
/*! \mainpage
  This package deals with transformation of space and time coordinates. You will find usefull examples in 
  the main directory, which can be used on a day to day basis :
  \arg <B>gps2ts</B> transforms a GPS time in timestamp. The GPS reference is 00h00 UT (midnight) on 
  January 6, 1980. The GPS time is the number of seconds ellapsed since this time.
  \arg <B>gps2unix</B> goes from GPS time to UNIX time, which is the number of seconds ellapsed since 
  01/01/1970 00h00 UTC NOT corrected for leap seconds so that there is a constant shift between GPS time 
  and UNIX time.  
  \arg <B>gps2utcs</B> returns for a given GPS time the corresponding UTC seconds. The UTC seconds are the 
  number of seconds ellapsed since 01/01/1970 00h00 UTC corrected for leap seconds. 
  \arg <B>ts2gps</B> transforms a timestamp in GPS time.
  \arg <B>ts2utcs</B> transforms a timestamp in UTC seconds.
  \arg <B>unix2gps</B> is the gps2unix reciprocal.
  \arg <B>utcs2gps</B> is the gps2utcs reciprocal.
  \arg <B>utcs2ts</B> is the ts2utcs reciprocal.
  \arg <B>ExampleREADME</B> is the piece of code used by the example in the README. Basically, it computes the 
  galactic coordinates corresponding to a given elevation, azimuth and UTC date.
  
  <B>Author :</B> <a href="mailto:revenu@in2p3.fr">B. Revenu</a>\n
  \n
  <B>Date :</B> January 2010\n
  \n
  <B>To Do :</B> The computation of precession is not yet completely included here so that all the coordinates 
  are referred to equinox J2000. The effect of nutation is not taken into account for the moment.
*/
// End the main page documentation



// Short description of the file
/*! \file STClibrary.h
  @brief Defines time and space coordinates.\n
*/



// Detailed description of the file
/*! \file STClibrary.h
  Time and dates :\n
  The julian day is a continuous count of days and fractions from the beginning of the year -4712. By tradition, 
  the julian day begins at Greenwich mean noon, that is, at 12h UT (Universal Time). The variable "day" is always 
  expressed in unit and fraction of day.\n
  \n
  Coordinates:\n
  ################# ATTENTION ##################\n
  ALL the angles (longitude, latitude, lambda, beta, l, b, dec)\n 
  are given in decimal degrees EXCEPT ra (right ascension)\n 
  which is in decimal hours.\n
  #############################################\n
  Remember : the longitude is measured positively eastwards.\n
  The azimuth is measured westwards from the South.\n
*/


#include <stdio.h>
#include <stdlib.h>
#include "STCconstants.h"
#include "LatLong-UTMConversion.h"

// Time *********************************************************************
//! Converts sexagesimal values into decimal value
double sexa2dec(int hr, int min, double sec);

//! Converts a decimal value into sexagesimal values
void dec2sexa(double dec, int * hr, int * min, double *sec);

//! This function is the modulo function. Returns x mod y
double mod(double x, double y = 0);

//! Conversion of UTC hour, min, sec in fraction of day. The fraction of the day is between 0 and 1
double hms2day(int hour, int min, double sec);

//! Conversion of a fraction of day into hour, min and sec. The fraction of the day must be between 0 and 1
void day2hms(double day, int * hour, int * min, double * sec);

/*!
  @brief converts gps seconds since january 6, 1980 00h00 UTC into date (UTC)
  @param gps_seconds is the value of time given by the GPS receiver+ns
  @param year is the corresponding year
  @param month is the corresponding month
  @param day is the corresponding day
  @param time is the corresponding UTC-time
*/
void gps2utcdate(double gps_seconds, int * year, int * month, int * day, double * time);

//! Converts UTC seconds to UTC julian day
void utcs2jd(double utcseconds, double * jd);

//! Converts julian day to UTC seconds
void jd2utcs(double jd, double * utcseconds);

//! Converts UTC seconds to UTC date
void utcs2date(double utcseconds, int* year, int* month, int* day, double* utc);

//! Converts a UTC date to UTC seconds
void date2utcs(int year, int month, int day, double utc, double * utcseconds);

//! Converts a UTC date (with decimal day) into Julian Date
void date2jd(int year, int month, double day, double * jd);

//! Converts a UTC date (with decimal hour) into Julian Date
void date2jd(int year, int month, int day, double hour, double * jd);

//! Converts a timestamp second into utc seconds
void ts2utcs(const char* ts, double *utcsecond);

/*! Converts a UTC date into Modified Julian Date. The modified Julian day (MJD) is Julian day - 2400000.5. 
  Contrarily to the JD, the MJD begins at Greenwich mean midnight !
*/
void date2mjd(int year, int month, double day, double * mjd);

///! Reciprocal of date2jd
void jd2date(double jd, int * year, int * month, double * day);

//! Returns the amount of time in days between 2 julian dates
double dateinterval(double jd1, double jd2);

//! Returns the amount of time in days between 2 UTC dates
double dateinterval(int year, int month, double day,
                    int year1, int month1, double day1);

//! Given a date, this function calculates the day of the week (0=sunday, 7=saturday)
int daynamefromdate(int year, int month, double day);

//! Day number in the year
int date2daynumber(int year, int month, int day);

/*! Given a year and a the number of the day in the year, this function converts into month and day. This is 
  the reciprocal of date2daynumber.
*/
void daynumber2date(int year, int daynb, int * month, int * day);

//! Check wether the given year is leap year or not. Rreturns 1 or 0.
int leapyear(int year);

//! Computes the easter sunday of the input year
void easterdate(int year, int * month, int * day);

/*! Converts the Julian Date (UTC) into the LOCAL sidereal time (in hours)\n
  (for the Greenwich sidereal time, just set longitude to zero)\n
  \n
  The result lst is in decimal hours. The longitude is the geographical longitude of the observation site and 
  is expressed in decimal degrees. REMEMBER : the longitudes are measured POSITIVELY EASTWARDS. For instance, 
  Washington = -77, Vienna = +16 (opposite to J. Meeus book). This convention is the same than that of the 
  IAU but is the opposite than those of the other planets. The parameter refjd is useful when you need high 
  accuracy so that the full jd is jd+refjd.
*/
void ct2lst(double longitude, double jd, double * lst, double refjd=0);
// *********************************************************************************



// coordinates *********************************************************************
/*! Computes the precession coefficients between the two equinoxes. The output coefficients are in 
  arcseconds (decimal)
 */
void precess_coeffs(int equinox_input, int equinox_output,
                    double * zeta, double * z, double * theta);

/*! Computes the value of the ecliptic obliquity angle for this value of JD. The returned value is in 
  decimal degrees
*/
double ecliptic_obliquity(double jd);

//! Converts equatorial coordinates with ra in decimal hours into ecliptic coordinates.
void radec2ecl(double ra, double dec, double * lambda,
               double * beta, int eqn=kSTC::RefEquinox);

//! Converts ecliptic coordinates into equatorial coordinates. Remind that ra is returned in decimal hours
void ecl2radec(double lambda, double beta, double * ra,
               double * dec, int eqn=kSTC::RefEquinox);

//! Converts equatorial coordinates with ra in decimal hours into galactic coordinates.
void radec2gal(double ra, double dec, double * l,
               double * b, int eqn=kSTC::RefEquinox);

//! Converts galactic coordinates into equatorial coordinates. Remind that ra is returned in decimal hours
void gal2radec(double l, double b, double * ra,
               double * dec, int eqn=kSTC::RefEquinox);

/*! Converts azimuth-elevation coordinates into equatorial coordinates. If refjd is not equal to 0 then the 
  full julian day is jd+refjd
*/
void azel2radec(double latitude, double longitude, double az, double el,
                double jd, double * ra, double * dec, double refjd=0);

/*! Converts equatorial coordinates into azimuth-elevation coordinates. If refjd is not equal to 0 then the 
  full julian day is jd+refjd
*/
void radec2azel(double latitude, double longitude, double ra, double dec,
                double jd, double * az, double * el, double refjd=0);

//! Converts ecliptic coordinates into galactic coordinates
void ecl2gal(double lambda, double beta, double * l, double * b, int eqn=kSTC::RefEquinox);

//! onverts galactic coordinates into ecliptic coordinates
void gal2ecl(double l, double b, double * lambda, double * beta, int eqn=kSTC::RefEquinox);

/*! Converts ecliptic coordinates into azimuth-elevation coordinates. Iif refjd is not equal to 0 then the 
  full julian day is jd+refjd
*/
void ecl2azel(double latitude, double longitude, double lambda,
              double beta, double jd, double * az, double * el, double refjd=0);

/*! Converts azimuth-elevation coordinates into ecliptic coordinates. If refjd is not equal to 0 then the full 
  julian day is jd+refjd
*/
void azel2ecl(double latitude, double longitude, double az, double el,
              double jd, double * lambda, double * beta, double refjd=0);

/*! Converts galactic coordinates into azimuth-elevation coordinates. If refjd is not equal to 0 then the full 
  julian day is jd+refjd
*/
void gal2azel(double latitude, double longitude, double l, double b,
              double jd, double * az, double * el, double refjd=0, int eqn=kSTC::RefEquinox);

/*! Converts azimuth-elevation coordinates into galactic coordinates. If refjd is not equal to 0 then the full 
  julian day is jd+refjd
*/
void azel2gal(double latitude, double longitude, double az, double el,
              double jd, double * l, double * b, double refjd=0, int eqn=kSTC::RefEquinox);

//! Converts the latitude on an ellipsoid into isometric latitude
void Latitude2IsoLatitude(double latitude, double ecc, double& isolatitude);

//! Converts galactic coordinates to Super-galactic coordinates
void gal2Sgal(double ll, double bb, double &sll, double &sbb);


//! Converts Super-galactic coordinates to galactic coordinates
void Sgal2gal(double sll, double sbb, double &ll, double &bb);

#endif
