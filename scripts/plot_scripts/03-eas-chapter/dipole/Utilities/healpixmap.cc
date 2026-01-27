#include <fstream>
#include <algorithm>

#include "healpixmap.h"
#include "harmotools.h"
#include "common.h"
#include "STClibrary.h"
#include "maptools.h"

using namespace std;



THealpixMap::THealpixMap(unsigned int nside, char thecoordsys)
{
  coordsys[0]=thecoordsys;
  strcpy(ordering,"RING");
  _nside=nside;
  _npix=nside2npix(nside);
  resize(_npix);
  fIpix.resize(_npix);
  for(unsigned int  k=0; k<_npix; k++)
    {
      at(k)=0.;
      fIpix[k] = k;
    }
}



THealpixMap::THealpixMap(const vector<double> & v, char thecoordsys)
{
  coordsys[0]=thecoordsys;
  strcpy(ordering,"RING");
  _npix=v.size();
  _nside=(long)pow(_npix*1./12,0.5);
  resize(_npix);
  fIpix.resize(_npix);
  for(unsigned int i=0; i<_npix; i++)
    {
      at(i)=v[i];
      fIpix[i] = i;
    }
}



THealpixMap::THealpixMap(const THealpixMap & m)
{
  *this = m;
}



THealpixMap::THealpixMap(char * filename)
{
  float *themap;
  long nside;
  themap=read_healpix_map(filename,&nside,coordsys,ordering);
  long npix=nside2npix(nside);
  _npix=npix;
  _nside=nside;
  if(!coordsys[0]) coordsys[0]='G';
  resize(_npix);
  fIpix.resize(_npix);
  for(unsigned long i=0;i<_npix;i++)
    {
      at(i) = themap[i];
      fIpix[i] = i;
    }
  delete themap;
}



THealpixMap::~THealpixMap()
{
  clear();
}



double THealpixMap::Total() const
{
  double bla = 0.;
  for(unsigned int i = 0; i < _npix; i++) bla += at(i);
  return bla;
}



double THealpixMap::Mean() const
{
  return Total()/_npix;
}



double THealpixMap::RMS() const
{
  double bla = 0.;
  for(unsigned int i=0; i<_npix; i++) bla += at(i)*at(i);
  return sqrt(bla/_npix-pow(Mean(),2));
}



double THealpixMap::Min() const
{
  return *min_element(begin(),end());
}



double THealpixMap::Max() const
{
  return *max_element(begin(),end());
}



void THealpixMap::GetMaxPosition(double & theta, double & phi, unsigned int & ipix) const
{
  vector<double>::const_iterator iter=find(begin(),end(),Max());
  unsigned int shift=iter-begin();
  ipix=shift;
  if( strncmp(ordering,"RING",4)==0 ) pix2ang_ring(_nside,ipix,&theta,&phi);
  else if( strcmp(ordering,"NEST")==0 ) pix2ang_nest(_nside,ipix,&theta,&phi);
  if( coordsys[0] == 'G' || coordsys[0] == 'Q' )
    {
      theta=(PiOver2-theta)*RTOD;
      phi*=RTOD;
    }
}



int THealpixMap::WriteFits(char * filename) const
{
  float *signal = new float[_npix];
  string fn(filename);
  fn="!"+fn;
  for(unsigned int i=0; i<_npix; i++) signal[i]=(float)at(i);
  char nest = strncmp(Ordering(),"RING",4);
  char* pTemp = strdup(fn.c_str());
  char * tmp=new char[strlen(coordsys)+1];
  strcpy(tmp,coordsys);
	//Note: edit by J. Biteau due to healpix upgrade//////
	//TODO: status is now meaningless
	//int status= write_healpix_map(signal,NSide(),pTemp,nest,tmp);
  int status = 1;
	///////////////////////////////////////////////////////
	write_healpix_map(signal,NSide(),pTemp,nest,tmp);
  delete [] tmp;
  free(pTemp);
  delete [] signal;
  return status;
}



void THealpixMap::SetNSide(unsigned int nside)
{
  coordsys[0] = 'G';
  strcpy(ordering,"RING");
  _nside = nside;
  _npix = nside2npix(nside);
  resize(_npix);
  fIpix.clear();
  fIpix.resize(_npix);
  for(unsigned int k = 0; k < _npix; k++ )
    {
      at(k) = 0.;
      fIpix[k] = k;
    }
}



double THealpixMap::GetPixSize() const {return sqrt(4.*M_PI*RTOD*RTOD/_npix);}



const char * THealpixMap::Ordering() const {return ordering;}



char THealpixMap::CoordSys() const {return* coordsys;}



const char * THealpixMap::Units() const {return units;}



void THealpixMap::SetOrdering(const char* pOrdering) {strcpy(ordering, pOrdering);}



void THealpixMap::SetCoordSys(char coord) {coordsys[0] = coord;}



unsigned int THealpixMap::NSide() const {return _nside;}



unsigned int THealpixMap::NPix() const {return _npix;}



THealpixMap THealpixMap::Map2Map(int nside_out)
{
  return ud_grade_healpix(*this, nside_out);
}



vector<double> THealpixMap::Map2Cl(int lmax) const
{
  return anafast_healpix(*this, lmax);
}



vector<vector<complex<double> > > THealpixMap::Map2Alm(int lmax)
{
  return anafast_alm_healpix(*this, lmax);
}



void THealpixMap::Cl2Map(vector<double> & cl, float beam)
{
  THealpixMap tmp(synfast_healpix(NSide(),cl.size()-1,cl,beam));
  coordsys[0] = tmp.CoordSys();
  strcpy(ordering, tmp.Ordering());
  strcpy(units, tmp.Units());
  _nside=tmp.NSide();
  _npix=tmp.NPix();
  resize(_npix);
  for(unsigned int i = 0; i < _npix; i++) at(i) = tmp[i];
}



void THealpixMap::Alm2Map(vector<vector<complex<double> > > & alm, float beam)
{
  THealpixMap tmp(synfast_alm_healpix(NSide(),alm.size()-1,alm,beam));
  coordsys[0] = tmp.CoordSys();
  strcpy(ordering, tmp.Ordering());
  strcpy(units, tmp.Units());
  _nside=tmp.NSide();
  _npix=tmp.NPix();
  resize(_npix);
  for(unsigned int i = 0; i < _npix; i++) at(i) = tmp[i];
}



THealpixMap THealpixMap::Filter(const vector<double> & theta, const vector<double> & lobe)
{
  long lmax = NSide()*3-1;
  vector<double> bl = legendrelobe(theta,lobe,lmax);
  THealpixMap themap(FiltBl(bl));
  return themap;
}



THealpixMap THealpixMap::IntBeam(const vector<double> & theta, const vector<double> & lobe)
{
  // Normalisation of the beam : maximum of one
  double max=*max_element(lobe.begin(),lobe.end());
  long sz = lobe.size();
  vector<double> lobenorm(sz);
  for(long i=0;i<sz;i++) lobenorm[i]=lobe[i]*1./max*sin(theta[i]*DTOR);
  double area=integrate_nc5(theta,lobenorm)*2*M_PI/DTOR;
  THealpixMap mapfilt;
  mapfilt = Filter(theta,lobe);
  mapfilt = mapfilt*(double)(mapfilt.size()*area/(4.*M_PI/(DTOR*DTOR)));
  return mapfilt;
}



THealpixMap THealpixMap::FiltBl(const vector<double> & bl)
{
  unsigned long lmax=3*NSide()-1;
  if(bl.size()<lmax+1) lmax = bl.size()+1;
  vector<vector<complex<double> > > alm;
  alm=Map2Alm(lmax);
  for(unsigned long l = 0;l < lmax+1; l++)
    {
      for(unsigned long m = 0; m < l+1; m++) alm[l][m] = alm[l][m]*bl[l];
    }
  
  THealpixMap mapfilt(NSide());
  mapfilt.Alm2Map(alm);
  return mapfilt;
}



THealpixMap THealpixMap::ChCoordSys(char newcoordsys) const
{
  if( newcoordsys == coordsys[0] || ((newcoordsys != 'Q') && (newcoordsys != 'G')) ){
	  THealpixMap tmp = (*this);
		tmp.coordsys[0] = newcoordsys;
		return tmp;
	}
  THealpixMap tmp = (*this)*0.;
  if( !strcmp(ordering,"RING") )// RING scheme
    { 
      double * theta = new double[_npix];
      double * phi = new double[_npix];
      for(unsigned int i = 0; i < _npix; i++) pix2ang_ring(_nside,i,&theta[i],&phi[i]);
      tmp.coordsys[0] = newcoordsys;
      if( coordsys[0] == 'G' && newcoordsys == 'Q' )
        {
          long ipix;
          double l, b;
          for(unsigned int i = 0; i< _npix; i++)
            {
              radec2gal(phi[i]/(15.*DTOR),90.-theta[i]/DTOR,&l,&b);
              ang2pix_ring(_nside,(90.-b)*DTOR,l*DTOR,&ipix);
              tmp[i] = at(ipix);
            }
        }
      else if( coordsys[0] == 'Q' && newcoordsys == 'G' )
        {
          long ipix;
          double ra, dec;
          for( unsigned int i=0;i<_npix;i++ )
            {
              gal2radec(phi[i]/DTOR,90.-theta[i]/DTOR,&ra,&dec);
              ra *= 15; // hr to deg
              ang2pix_ring(_nside,(90.-dec)*DTOR,ra*DTOR,&ipix);
              tmp[i] = at(ipix);
            }
        }
      delete [] theta;
      delete [] phi;
    }
  return tmp;
}



vector<long> THealpixMap::Ip(const vector<double> & l, const vector<double> & b) const
{
  long sz = l.size();
  vector<long> ip(sz);
  long pix;
  if(!strcmp(ordering,"RING"))
    { 
      for(long i = 0; i < sz; i++)
        {
          ang2pix_ring(NSide(),(90.-b[i])*DTOR,l[i]*DTOR,&pix);
          ip[i] = pix;
        }
    }
  else
    {
      for(long i = 0; i < sz; i++)
        {
          ang2pix_nest(NSide(),(90.-b[i])*DTOR,l[i]*DTOR,&pix);
          ip[i] = pix;
        }
    }
  return ip;
}



long THealpixMap::Ip(double l, double b) const
{
  long ip;
  if(!strcmp(ordering,"RING")) ang2pix_ring(NSide(),(90.-b)*DTOR,l*DTOR,&ip);
  else ang2pix_nest(NSide(),(90.-b)*DTOR,l*DTOR,&ip);

  return ip;
}



vector<double> THealpixMap::Values(const vector<double> & l, const vector<double> & b) const
{
  vector<long> ip;
  ip = Ip(l,b);
  vector<double> val;
  long sz = ip.size();
  val.resize(sz);
  for(long i=0;i<sz;i++) val[i]=at(ip[i]);
  return val;
}



double THealpixMap::Value(double l, double b) const
{
  long ip;
  ip = Ip(l,b);
  double val = at(ip);
  return val;
}



void THealpixMap::GiveLB(const vector<long> & ipix, vector<double> & l,  vector<double> & b) const
{
  long sz = ipix.size();
  l.resize(sz);
  b.resize(sz);
  double theta,phi;
  if( !strcmp(ordering,"RING") )
    { 
      for(long i = 0; i < sz; i++)
        {
          pix2ang_ring(NSide(),ipix[i],&theta,&phi);
          l[i] = phi/DTOR;
          b[i] = 90.-theta/DTOR;
        }
    }
  else
    {
      for(long i = 0; i < sz; i++)
        {
          pix2ang_nest(NSide(),ipix[i],&theta,&phi);
          l[i] = phi/DTOR;
          b[i] = 90.-theta/DTOR;
        }
    }  
}



void THealpixMap::GiveLB(long ipix, double & l, double & b) const
{
  double theta,phi;
  if( !strcmp(ordering,"RING") )
    { 
      pix2ang_ring(NSide(),ipix,&theta,&phi);
      l = phi/DTOR;
      b = 90.-theta/DTOR;
    }
  else
    {
      pix2ang_nest(NSide(),ipix,&theta,&phi);
      l = phi/DTOR;
      b = 90.-theta/DTOR;
    }
}  



vector< vector<double> >  THealpixMap::Map2MollProj(int sizeX, int sizeY, double background, double decLimit, char* filename) const
{
  int nj = sizeX;
  int ni = sizeY;
  vector<vector<double> > proj(nj);
  for(int i = 0; i < sizeX; i++) proj[i].resize(ni);
  double yd, theta, xa, phi, ra, dec;
  long pix;
  for(int i = 0; i < ni; i++)// latitudes
    {
      yd = (i+0.5)/ni-0.5;
      for(int j = 0; j < nj; j++)// longitudes
        {
          xa = (j+0.5)/nj-0.5;
          if( XYtoAngMollweide(xa,yd,theta,phi) )
            {
              gal2radec(phi*RTOD, 90-theta*RTOD, &ra, &dec);
              if(!strcmp(ordering,"RING")) ang2pix_ring(_nside,theta,phi,&pix);
              else ang2pix_nest(_nside,theta,phi,&pix);
              if (dec <= decLimit) proj[j][i] = at(pix);
              else proj[j][i] = background;
            }
          else proj[j][i] = background;
        }
    }
  
  if( strcmp(filename," ") )
    {
      cout << "Writing Mollweide projection to file : " << filename << endl;
      write_fits_image(proj,filename);
    }
  cout << "Returning " << endl;
  return proj;
}



vector<vector<double> > THealpixMap::Map2LambertAzimuthalProj(int sizeX, int sizeY, double longRef, double latRef, double radius, double background, double decLimit, char *filename) const
{
  long pix;
  vector< vector<double> > proj(sizeX);
  for(int i = 0; i < sizeX; i++) proj[i].resize(sizeY);
  double minL = -radius*DTOR;
  double maxL = radius*DTOR;
  double minB = -radius*DTOR;
  double maxB = radius*DTOR;
  double thex, they, theta, phi, ra, dec;
  double stepX = (maxL-minL)/(sizeX-1);
  double stepY = (maxB-minB)/(sizeY-1);
  bool status;
  double cosLatRef = cos(latRef*DTOR), sinLatRef = sin(latRef*DTOR);
  double longRefRad = longRef*DTOR;
  for(int i = 0; i < sizeX; i++)
    {
      thex = minL+stepX*i;
      for(int j = 0; j < sizeY; j++)
        {
          they = minB+stepY*j;
          if( status = XYtoAngLambertAzimuthal(thex,they,theta,phi,longRefRad,cosLatRef,sinLatRef) )
            {
              gal2radec(phi*RTOD, theta*RTOD, &ra, &dec);
              if( !strcmp(ordering,"RING") ) ang2pix_ring(_nside,PiOver2-theta,phi,&pix);
              else ang2pix_nest(_nside,PiOver2-theta,phi,&pix);
              if (dec <= decLimit) proj[i][j] = at(pix);
              else proj[i][j] = background;
            }
          else proj[i][j] = background;
        }
    }
  if( strcmp(filename," ") )
    {
      cout << "Writing Lambert azimuthal projection to file : " << filename << endl;
      write_fits_image(proj,filename);
    }
  return proj;
}



vector<double> THealpixMap::FindMaxima(vector<long> & ipmax) const
{
  string healpixdir(getenv("HEALPIX_DIR"));
  unsigned int seed;
  struct timeval mytimeval;
  struct timezone mytimezone;
  gettimeofday(&mytimeval,&mytimezone);
  seed = (unsigned int) (mytimeval.tv_usec + (mytimeval.tv_sec % 1000)*1000000);
  srandom(seed);
  char random_str[256];
  sprintf(random_str,"%d",(int)random());
  // Definition of the name of the files
  char fitsmapin[256], fitsmapout[256], maxifile[256];
  char minifile[256],txtfile[256],maxifile_ord[256];
  char * tmpbatch;
  tmpbatch = getenv("TMPBATCH");
  if ( tmpbatch == NULL )
    {
      tmpbatch = new char[256];
      sprintf(tmpbatch,"/tmp");
    }
  sprintf(fitsmapin,"%s/tmpmapin%s.fits",tmpbatch,random_str);
  sprintf(fitsmapout,"%s/tmpmapout%s.fits",tmpbatch,random_str);
  sprintf(maxifile,"%s/maxifile%s.txt",tmpbatch,random_str);
  sprintf(maxifile_ord,"%s/maxifile_ord%s.txt",tmpbatch,random_str);
  sprintf(minifile,"%s/minifile%s.txt",tmpbatch,random_str);
  sprintf(txtfile,"%s/tmptxt%s.txt",tmpbatch,random_str);
  // Temporary FITS file
  WriteFits(fitsmapin);
  // Temporay txt file
  ofstream oftxt(txtfile);
  oftxt << "infile = " << fitsmapin << endl;
  oftxt << "outmap = " << fitsmapout << endl;
  oftxt << "maxfile = " << maxifile << endl;
  oftxt << "minfile = " << minifile << endl;
  oftxt.close();
  // system command
  char commande[1000];
  sprintf(commande,"%s/bin/hotspots_cxx %s > /dev/null",healpixdir.c_str(),txtfile);
  system(commande);
  // Sort the file
  sprintf(commande,"cat %s | awk '{print $2,$1}' | sort -nr > %s",maxifile,maxifile_ord);
  system(commande);
  // Read the file
  ifstream ifs(maxifile_ord);
  vector<double> maxivals;
  long ip;
  double theval;
  ipmax.clear();
  while( ifs >> theval )
    {
      ifs >> ip;
      maxivals.push_back(theval);
      ipmax.push_back(ip);
    }
  // effacement des fichiers temporaires
  sprintf(commande,"rm -f %s %s %s %s %s %s",fitsmapin,fitsmapout,txtfile,minifile,maxifile,maxifile_ord);
  system(commande);
 
  delete [] tmpbatch;
 
  return maxivals;
}



THealpixMap THealpixMap::operator+(THealpixMap map2) const
{
  THealpixMap temp(_nside);
  if(map2.NSide()==_nside) for(unsigned int i=0;i<_npix;i++) temp[i]=at(i)+map2[i]; 
  
  return temp;
}



THealpixMap THealpixMap::operator-(THealpixMap map2) const
{
  THealpixMap temp(_nside);
  if( map2.NSide() == _nside ) for(unsigned int i = 0; i < _npix; i++) temp[i] = at(i)-map2[i]; 
  
  return temp;
}



THealpixMap THealpixMap::operator*(THealpixMap map2) const
{
  THealpixMap temp(_nside);
  if( map2.NSide() == _nside ) for(unsigned int i = 0; i < _npix; i++) temp[i] = at(i)*map2[i]; 
  
  return temp;
}



THealpixMap THealpixMap::operator/(THealpixMap map2) const
{
  THealpixMap temp(_nside);
  if( map2.NSide() == _nside ) for(unsigned int i = 0; i < _npix; i++) if(map2[i] != 0) temp[i] = at(i)/map2[i]; 
  
  return temp;
}



THealpixMap THealpixMap::operator+(double x) const
{
  THealpixMap temp(_nside);
  for(unsigned int i = 0; i < _npix; i++) temp[i] = at(i)+x; 

  return temp;
}



THealpixMap THealpixMap::operator-(double x) const
{
  THealpixMap temp(_nside);
  for(unsigned int i = 0; i < _npix; i++) temp[i] = at(i)-x; 

  return temp;
}


THealpixMap THealpixMap::operator*(double x) const
{
  THealpixMap temp(_nside);
  for(unsigned int i = 0; i < _npix; i++) temp[i] = at(i)*x; 

  return temp;
}



THealpixMap THealpixMap::operator/(double x) const
{
  THealpixMap temp(_nside);
  for(unsigned int i = 0; i < _npix; i++) temp[i] = at(i)/x; 

  return temp;
}



THealpixMap & THealpixMap::operator=(const THealpixMap& m)
{
  coordsys[0] = m.CoordSys();
  strcpy(ordering, m.Ordering());
  strcpy(units, m.Units());
  _nside=m.NSide();
  _npix=m.NPix();
  resize(_npix);
  for(unsigned int i = 0; i < _npix; i++) at(i) = m[i];

  return *this;
}



THealpixMap & THealpixMap::operator = (double x)
{
  for(unsigned int i = 0; i < _npix; i++) at(i) = x; return *this;
}



THealpixMap & THealpixMap::operator *= (double x) 
{
  for(unsigned int i = 0; i < _npix; i++) at(i) *= x; return *this;
}



THealpixMap & THealpixMap::operator /= (double x)
{
  for(unsigned int i = 0; i < _npix; i++) if( x!=0 ) at(i) /= x; return *this;
}



THealpixMap & THealpixMap::operator += (double x)
{
  for( unsigned int i = 0; i < _npix; i++) at(i) += x; return *this;
}



THealpixMap & THealpixMap::operator -= (double x)
{
  for( unsigned int i=0; i<_npix; i++) at(i) -= x; return *this;
}

