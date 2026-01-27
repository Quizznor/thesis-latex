#include <cmath>
#include <iostream>

#include "STClibrary.h"
#include "STCconstants.h"

using namespace std;

/**
 * @file STClibrary.cc
 * @brief this file implements the functions defined in STClibrary.h
 * some astronomical data are defined
 */

double sexa2dec (int hr, int min, double sec)
{
	int sign = 1;
	if (hr < 0)
	{
		sign = -1;
		hr *= sign;
	}
	return sign * (hr + min / 60. + sec / 3600.);
}

void dec2sexa (double dec, int *hr, int *min, double *sec)
{
	double sign = 1;
	double epsilon = 1e-10;
	if (dec < 0)
	{
		sign = -1;
		dec *= sign;
	}
	double deg = floor (dec);
	double ex_min = 60. * (dec - deg);
	double cmin = ceil (ex_min);
	if (cmin - ex_min < epsilon)
		ex_min = cmin;
	double lmin = floor (ex_min);
	double lsec = 60. * (ex_min - lmin);
	if (fabs (lsec) < epsilon)
		lsec = 0;

	*hr = (int) (sign * deg);
	*min = (int) lmin;
	*sec = lsec;
}

double mod (double d, double periode)
{
	if (d >= 0)
		return d - floor (d / periode) * periode;
	else
		return d - ceil (d / periode) * periode + periode;
}

/**
 * @var ecliptic_obliquity_coeffs
 * @brief are the coefficients used to compute the exact ecliptic obliquity.\n
 * The values are in arcseconds.
 */
static double ecliptic_obliquity_coeffs[10] = {
	-4680.93,
	-1.55,
	1999.25,
	-51.38,
	-249.67,
	-39.05,
	7.12,
	27.87,
	5.79,
	2.45
};

/**
 * @var zeta_coeff
 * @brief are the coefficients used to compute the precession coefficients.\n
 * The values are in arcseconds.
 */
static double zeta_coeff[6] = {
	2306.2181,
	1.39656,
	-0.000139,
	0.30188,
	-0.000344,
	0.017998
};

/**
 * @var z_coeff
 * @brief are the values used to compute the precession coefficients.\n
 * The values are in arcseconds.
 */
static double z_coeff[6] = {
	2306.2181,
	1.39656,
	-0.000139,
	1.09468,
	0.000066,
	0.018203
};

/**
 * @var theta_coeff
 * @brief are the values used to compute the precession coefficients.\n
 * The values are in arcseconds.
 */
static double theta_coeff[6] = {
	2004.3109,
	-0.85330,
	-0.000217,
	-0.42665,
	-0.000217,
	-0.041833
};

/**
 * @var sidereal_coeff
 * @brief are the coefficients used to compute the sidereal time.\n
 * The values are in degrees.
 */
static double sidereal_coeff[4] = {
	280.46061837,
	360.98564736629,
	0.000387933,
	38710000
};

double hms2day (int hour, int min, double sec)
{
	return (hour + min / 60. + sec / 3600.) / 24.;
}

void day2hms (double day, int *hour, int *min, double *sec)
{
	double frac = (day - floor (day)) * 24.;
	*hour = (int) floor (frac);
	*min = (int) (floor (60. * (frac - (*hour))));
	*sec = 3600. * (frac - (*min) / 60. - (*hour));
	//  cout << *sec-60 << endl;
	if (*sec == 60)
	{
		//      cout << "icisec" << endl;
		*min++;
		*sec = 0;
	}
	if (*min == 60)
	{
		*hour++;
		*min = 0;
	}
	//  cout << "day2hms : " << frac << " " << day << " " << *hour << " " << *min << " " << *sec << endl;
}

void utcs2jd (double utcseconds, double *jd)
{
	*jd = kSTC::UnixStartJD + utcseconds / 3600. / 24.;
}

void jd2utcs (double jd, double *utcseconds)
{
	*utcseconds = 3600. * 24. * (jd - kSTC::UnixStartJD);
}

void utcs2date (double utcseconds, int *year, int *month, int *day, double *utc)
{
	double jd;
	utcs2jd(utcseconds,&jd);
	double decimalday;
	jd2date(jd,year,month,&decimalday);
	*day = (int) floor(decimalday);
	int hour, min;
	double sec;
	day2hms (decimalday, &hour, &min, &sec);
	*utc = hour + min * 1. / 60 + sec * 1. / 3600;
}

void date2utcs (int year, int month, int day, double utc, double *utcseconds)
{
	double jd;
	date2jd(year,month,day,utc,&jd);
	jd2utcs(jd,utcseconds);
}

void date2jd (int year, int month, double day, double *jd)
{
	if (month <= 2)
	{
		year--;
		month += 12;
	}
	int b;
	if (year < 1582
			|| (year == 1582 && (month < 10 || (month == 10 && day < 15))))
	{
		// we are in the Julian calendar
		b = 0;
	}
	else
	{
		// we are in the Gregorian calendar
		int a = (int) (floor (year / 100.));
		b = 2 - a + a / 4;
	}
	*jd =
		floor (365.25 * (year + 4716)) + floor (30.6001 * (month + 1)) + day + b -
		1524.5;
}

void date2jd (int year, int month, int day, double hour, double *jd)
{
	int hr, min;
	double sec;
	dec2sexa (hour, &hr, &min, &sec);
	double dday = hms2day (hr, min, sec);	// fraction of day
	date2jd (year, month, day + dday, jd);
}

void date2mjd (int year, int month, double day, double *mjd)
{
	// contrary to the JD, the MJD begins at Greenwich mean midnight !
	double jd;
	date2jd (year, month, day, &jd);
	*mjd = jd - 2400000.5;
}

void jd2date (double jd, int *year, int *month, double *day)
{
	jd += 0.5;
	int z = (int) floor (jd);
	double f = jd - z;
	int a, b, c, d, e;
	int alpha;

	if (z < 2299161)
		a = z;
	else
	{
		alpha = (int) floor ((z - 1867216.25) / 36524.25);
		a = z + 1 + alpha - (int) floor (alpha / 4.);
	}

	b = a + 1524;
	c = (int) floor ((b - 122.1) / 365.25);
	d = (int) floor (365.25 * c);
	e = (int) floor ((b - d) / 30.6001);

	*day = b - d - floor (30.6001 * e) + f;
	*month = (e < 14) ? (e - 1) : (e - 13);
	*year = (*month < 3) ? (c - 4715) : (c - 4716);
}

double dateinterval (double jd1, double jd2)
{
	return jd1 - jd2;
}

double dateinterval (int year, int month, double day, int year1, int month1,
		double day1)
{
	double jd, jd1;
	date2jd (year, month, day, &jd);
	date2jd (year1, month1, day1, &jd1);
	return jd - jd1;
}

int daynamefromdate (int year, int month, double day)
{
	// 0 = sunday, ..., 6 = saturday
	double jd;
	date2jd (year, month, floor (day), &jd);	// same day but at 0h

	double tmp = jd + 1.5;
	return (int) (tmp - 7. * floor (tmp / 7.));
}

int leapyear (int year)
{
	static int sleap = 0;
	static int syear = 0;
	if (year == syear)
		return sleap;
	else
		syear = year;
	if (year - 4 * (int) floor (year / 4.) == 0)
	{
		if (year < 1582 || (year - 400 * (int) floor (year / 400.) == 0))
			sleap = 1;
		else
			sleap = 0;
	}
	else
		sleap = 0;
	return sleap;
}

int date2daynumber (int year, int month, int day)
{
	int k = 2 - leapyear (year);
	return (int) floor (275. * month / 9.) -
		k * (int) floor ((month + 9.) / 12.) + day - 30;
}

void daynumber2date (int year, int daynb, int *month, int *day)
{
	int k = 2 - leapyear (year);
	*month = (int) floor (9. * (k + daynb) / 275. + 0.98);
	if (daynb < 32)
		*month = 1;
	*day =
		daynb - (int) floor (275. * (*month) / 9.) +
		k * (int) floor ((*month + 9.) / 12.) + 30;
}

void easterdate (int year, int *month, int *day)
{
	// Easter Sunday is the first sunday after the full moon that
	// happens on or next after the March equinox

	int a, b, c, d, e, f, g, h, i, k, l, m, n, p, tmp;

	static int syear = 0;
	static int emonth = 0;
	static int eday = 0;

	if (year == syear)
	{
		*month = emonth;
		*day = eday;
	}
	else
	{
		syear = year;
		if (year > 1582)
		{
			// Gregorain calendar
			a = year - 19 * (int) floor (year / 19.);
			b = (int) floor (year / 100.);
			c = year - 100 * b;
			d = (int) floor (b / 4.);
			e = b - 4 * d;
			f = (int) floor ((b + 8.) / 25.);
			g = (int) floor ((b - f + 1.) / 3.);
			tmp = 19 * a + b - d - g + 15;
			h = tmp - 30 * (int) floor (tmp / 30.);
			i = (int) floor (c / 4.);
			k = c - 4 * i;
			tmp = 32 + 2 * e + 2 * i - h - k;
			l = tmp - 7 * (int) floor (tmp / 7.);
			m = (int) floor ((a + 11 * h + 22 * l) / 451.);
			tmp = h + l - 7 * m + 114;
			n = (int) floor (tmp / 31.);
			p = tmp - 31 * n;
			*month = n;
			*day = p + 1;		// on Sunday
		}
		else
		{
			// Julian calendar
			a = year - 4 * (int) floor (year / 4.);
			b = year - 7 * (int) floor (year / 7.);
			c = year - 19 * (int) floor (year / 19.);
			tmp = 19 * c + 15;
			d = tmp - 30 * (int) floor (tmp / 30.);
			tmp = 2 * a + 4 * b - d + 34;
			e = tmp - 7 * (int) floor (tmp / 7.);
			tmp = d + e + 114;
			f = (int) floor (tmp / 31.);
			g = tmp - 31 * f;
			*month = f;
			*day = g + 1;
		}
		// remember the results
		emonth = *month;
		eday = *day;
	}
	return;
}

void ct2lst (double longitude, double jd, double *lst, double refjd)
{
  // longitude stays in degrees here and is positive eastwards
  double dt;
  if( !refjd )
    {
      int intpjd = int(floor(jd));
      int intdt = intpjd-int(kSTC::J2000);
			double frac = jd-intpjd;
      dt = intdt+frac;
    }
  else
    {
      dt = jd + (refjd - kSTC::J2000);
    }
  double dt0 = dt / 36525.;
  double theta = sidereal_coeff[0] + dt * sidereal_coeff[1] + dt0 * dt0 * (sidereal_coeff[2] - dt0 / sidereal_coeff[3]);	// is in degrees
  // original J. Meeus: longitude is positive westwards  *lst = mod(theta-longitude,360.) * kSTC::Dtohr;
  // with the new convention for longitudes (positive eastwards) :
  *lst = mod (theta + longitude, 360.) * kSTC::Dtohr;
}

void precess_coeffs (int eqn1, int eqn2, double *zeta, double *z, double *theta)
{
	double gt = (eqn1 - kSTC::J2000) / 36525.;	// in Julian centuries since kSTC::J2000
	double dt = (eqn2 - eqn1) / 36525.;
	// zeta, z and theta are in arcseconds
	*zeta = dt * (zeta_coeff[0] +
			gt * (zeta_coeff[1] + zeta_coeff[2] * gt) +
			dt * (zeta_coeff[3] + zeta_coeff[4] * gt +
				zeta_coeff[5] * dt));

	*z = dt * (z_coeff[0] +
			gt * (z_coeff[1] + z_coeff[2] * gt) +
			dt * (z_coeff[3] + z_coeff[4] * gt + z_coeff[5] * dt));

	*theta = dt * (theta_coeff[0] +
			gt * (theta_coeff[1] + theta_coeff[2] * gt) +
			dt * (theta_coeff[3] + theta_coeff[4] * gt +
				theta_coeff[5] * dt));
}

double ecliptic_obliquity (double jd)
{
	// var is expressed in units of 10^4 Julian years from kSTC::J2000
	double var = (jd - kSTC::J2000) / (100 * 36525);
	double tmp = 0;
	double pow = 1;
	for (int i = 0; i < 10; i++)
	{
		pow *= var;
		tmp += ecliptic_obliquity_coeffs[i] * pow;
	}
	tmp /= 3600.;			// the correction is in arcseconds
	return kSTC::Eps_2000 + tmp;
}

void radec2ecl (double ra, double dec, double *lng, double *lat, int eqn)
{
	ra *= kSTC::Hrtor;		// from decimal hour to radians
	dec *= kSTC::DTOR;
	double eps = kSTC::Eps_2000;	// in degrees
	if (eqn != kSTC::RefEquinox)
	{
		double jd;
		date2jd (eqn, 1, 1.5, &jd);	// 1st january, noon
		eps = ecliptic_obliquity (jd);
	}
	eps *= kSTC::DTOR;
	double sra = sin (ra);
	double ce = cos (eps);
	double se = sin (eps);
	*lng = mod (atan2 (sra * ce + tan (dec) * se, cos (ra)) * kSTC::RTOD, 360.);	// in degrees
	*lat = asin (sin (dec) * ce - cos (dec) * se * sra) * kSTC::RTOD;	// in degrees
}

void ecl2radec (double lng, double lat, double *ra, double *dec, int eqn)
{
	lng *= kSTC::DTOR;
	lat *= kSTC::DTOR;
	double eps = kSTC::Eps_2000;
	if (eqn != kSTC::RefEquinox)
	{
		double jd;
		date2jd (eqn, 1, 1.5, &jd);	// 1.5 = 1st january, noon
		eps = ecliptic_obliquity (jd);
	}
	eps *= kSTC::DTOR;
	double sl = sin (lng);
	double ce = cos (eps);
	double se = sin (eps);
	*ra = mod (atan2 (sl * ce - tan (lat) * se, cos (lng)) * kSTC::Rtohr, 24.);	// in decimal hours
	*dec = asin (sin (lat) * ce + cos (lat) * se * sl) * kSTC::RTOD;	// in decimal degrees
}

void radec2gal (double ra, double dec, double *l, double *b, int eqn)
{
	// the equinox specifies the equinox used for ra and dec
	// lambda and beta are given with respect to this equinox
	if (eqn != 2000)
	{
		cout << "At this time, the only equinox possible is 2000" << endl;
		return;
	}
	ra *= kSTC::Hrtod;
	dec *= kSTC::DTOR;
	double dg = kSTC::DeltaG_2000 * kSTC::DTOR;
	double da = kSTC::DTOR * (kSTC::AlphaG_2000 - ra);
	double x = atan2 (sin (da), cos (da) * sin (dg) - tan (dec) * cos (dg));
	*l = mod (270. + kSTC::Lomega_2000 - x * kSTC::RTOD, 360.);
	*b =
		asin (sin (dec) * sin (dg) +
				cos (dec) * cos (dg) * cos (da)) * kSTC::RTOD;
}

void gal2radec (double l, double b, double *ra, double *dec, int eqn)
{
	if (eqn != 2000)
	{
		cout << "At this time, the only equinox possible is 2000" << endl;
		return;
	}
	l *= kSTC::DTOR;
	b *= kSTC::DTOR;
	double da = l - (90. + kSTC::Lomega_2000) * kSTC::DTOR;
	double dg = kSTC::DeltaG_2000 * kSTC::DTOR;
	double y = atan2 (sin (da), cos (da) * sin (dg) - tan (b) * cos (dg));
	*ra = mod ((y * kSTC::RTOD + kSTC::AlphaG_2000 - 180.) * kSTC::Dtohr, 24.);
	*dec =
		asin (sin (b) * sin (dg) + cos (b) * cos (dg) * cos (da)) * kSTC::RTOD;
}

void azel2radec (double latitude, double longitude, double az, double el,
		double jd, double *ra, double *dec, double refjd)
{
	double lst;
	ct2lst (longitude, jd, &lst, refjd);
	longitude *= kSTC::DTOR;	// positive eastwards
	latitude *= kSTC::DTOR;
	az *= kSTC::DTOR;
	el *= kSTC::DTOR;
	double slat = sin (latitude);
	double clat = cos (latitude);
	double saz = sin (az);
	double caz = cos (az);
	double sel = sin (el);
	double cel = cos (el);
	double ha = atan2 (saz, slat * caz + sel * clat / cel);
	ha *= kSTC::Rtohr;		// in hours
	*ra = mod (lst - ha, 24.);
	*dec = asin (slat * sel - clat * cel * caz) * kSTC::RTOD;	// in degrees
}

void radec2azel (double latitude, double longitude, double ra, double dec,
		double jd, double *az, double *el, double refjd)
{
	double lst;
	ct2lst (longitude, jd, &lst, refjd);
	latitude *= kSTC::DTOR;
	double ha = (lst - ra) * kSTC::Hrtor;	// in radians
	dec *= kSTC::DTOR;
	double slat = sin (latitude);
	double clat = cos (latitude);
	double sha = sin (ha);
	double cha = cos (ha);
	double sdec = sin (dec);
	double cdec = cos (dec);
	*el = asin (slat * sdec + clat * cdec * cha) * kSTC::RTOD;	// in degrees
	*az = mod (atan2 (sha, cha * slat - sdec * clat / cdec) * kSTC::RTOD, 360.);	// in degrees
}

void ecl2gal (double lambda, double beta, double *l, double *b, int eqn)
{
	// should use eqn
	eqn += 0;
	lambda *= kSTC::DTOR;
	beta *= kSTC::DTOR;
	double ae = kSTC::AlphaE_2000 * kSTC::DTOR;
	double de = kSTC::DeltaE_2000 * kSTC::DTOR;
	double dl = lambda - ae;
	double tmp = atan2 (-sin (de) * cos (dl) + cos (de) * tan (beta), sin (dl));
	*l = mod (kSTC::Eomega_2000 + tmp * kSTC::RTOD, 360.);
	*b =
		asin (sin (de) * sin (beta) +
				cos (de) * cos (beta) * cos (dl)) * kSTC::RTOD;
}

void gal2ecl (double l, double b, double *lambda, double *beta, int eqn)
{
	// should use eqn
	eqn += 0;
	l *= kSTC::DTOR;
	b *= kSTC::DTOR;
	double e0 = kSTC::Eomega_2000 * kSTC::DTOR;
	double de = kSTC::DeltaE_2000 * kSTC::DTOR;
	double dl = l - e0;
	double tmp = atan2 (sin (de) * sin (dl) - cos (de) * tan (b), cos (dl));
	*lambda = mod (90. + kSTC::AlphaE_2000 + tmp * kSTC::RTOD, 360.);
	*beta =
		asin (sin (de) * sin (b) + cos (de) * cos (b) * sin (dl)) * kSTC::RTOD;
}

void ecl2azel (double latitude, double longitude, double lambda,
		double beta, double jd, double *az, double *el, double refjd)
{
	double ra, dec;
	ecl2radec (lambda, beta, &ra, &dec, kSTC::RefEquinox);
	radec2azel (latitude, longitude, ra, dec, jd, az, el, refjd);
}

void azel2ecl (double latitude, double longitude, double az, double el,
		double jd, double *lambda, double *beta, double refjd)
{
	double ra, dec;
	azel2radec (latitude, longitude, az, el, jd, &ra, &dec, refjd);
	radec2ecl (ra, dec, lambda, beta, kSTC::RefEquinox);
}

void gal2azel (double latitude, double longitude, double l, double b,
		double jd, double *az, double *el, double refjd, int eqn)
{
	double ra, dec;
	gal2radec (l, b, &ra, &dec, eqn);
	radec2azel (latitude, longitude, ra, dec, jd, az, el, refjd);
}

void azel2gal (double latitude, double longitude, double az, double el,
		double jd, double *l, double *b, double refjd, int eqn)
{
	double ra, dec;
	azel2radec (latitude, longitude, az, el, jd, &ra, &dec, refjd);
	radec2gal (ra, dec, l, b, eqn);
}


void ts2utcs(const char* ts, double * utcs)
{
  int y,m,d,h,min,s;
  sscanf(ts,"%04d%02d%02d%02d%02d%02d",&y,&m,&d,&h,&min,&s);
  double jd;
  date2jd(y,m,d,h+min/60.+s/3600.,&jd);
  jd2utcs(jd,utcs);
}

void Latitude2IsoLatitude (double latitude, double ecc, double &isolatitude)
{
	double sl = sin (latitude * kSTC::DTOR);
	isolatitude =
		log (tan (M_PI / 4. + latitude * kSTC::DTOR / 2.) *
				pow ((1 - ecc * sl) / (1 + ecc * sl), ecc / 2.));
}

void gal2Sgal(double ll, double bb, double &sll, double &sbb){
	ll *= kSTC::DTOR;
	bb *= kSTC::DTOR;

	double sgl1,sgl2,sgl3,gl1,gl2,gl3; /* used for rotation */
	gl1 = cos (bb) * cos (ll);
	gl2 = cos (bb) * sin (ll);
	gl3 = sin (bb);
	sgl1 = -0.735742574804 * gl1 +  0.677261296414 * gl2 + 0.000000000000 * gl3;
	sgl2 = -0.074553778365 * gl1 + -0.080991471307 * gl2 + 0.993922590400 * gl3;
	sgl3 =  0.673145302109 * gl1 +  0.731271165817 * gl2 + 0.110081262225 * gl3;

	sbb = asin(sgl3) / kSTC::DTOR;
	sll = atan2(sgl2,sgl1) / kSTC::DTOR;

	/* Range Check */
	if (sbb < -90.0) sbb += 180.0; /* sbb is in the range -90 to 90 */
	if (sbb > +90.0) sbb -= 180.0; /* sbb is in the range -90 to 90 */
	if (sll < 0.0) sll += 360.0; /* sll is in the range 0 to 360 */
	if (sll > 360.0) sll -= 360.0; /* sll is in the range 0 to 360 */
}

void Sgal2gal(double sll, double sbb, double &ll, double &bb){
	double sgl1,sgl2,sgl3,gl1,gl2,gl3; /* used for rotation */
	sll *= kSTC::DTOR;
	sbb *= kSTC::DTOR;

	sgl1 = cos (sbb) * cos (sll);
	sgl2 = cos (sbb) * sin (sll);
	sgl3 = sin (sbb);
	gl1 = -0.735742574804 * sgl1 + -0.074553778365 * sgl2 + .673145302109 * sgl3;
	gl2 =  0.677261296414 * sgl1 + -0.080991471307 * sgl2 + .731271165817 * sgl3;
	gl3 =  0.000000000000 * sgl1 +  0.993922590400 * sgl2 + .110081262225 * sgl3;

	bb = asin(gl3) / kSTC::DTOR;
	ll = atan2(gl2,gl1) / kSTC::DTOR;

	/* Range Check */
	if (bb < -90.0) bb += 180.0; /* bb is in the range -90 to 90 */
	if (bb > +90.0) bb -= 180.0; /* bb is in the range -90 to 90 */
	if (ll < 0.0) ll += 360.0; /* ll is in the range 0 to 360 */
	if (ll > 360.0) ll -= 360.0; /* ll is in the range 0 to 360 */
}

