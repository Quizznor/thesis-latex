#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <tuple>
#include <sys/stat.h>
#include "common.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
// Toolkit
#include "maptools.h"
#include "projmap.h"
#include "common.h"
#include "STClibrary.h"
#include "healpixmap.h"

bool CheckFile(string fileName)
{
  char fileNameCopy[100];
  struct stat fileStat;
  strcpy(fileNameCopy,fileName.c_str());
  int status = lstat(fileNameCopy,&fileStat);

  if( status == -1 )
    {
      cout << "-----------------------------------------------------" << endl;
      cout << "File " << fileName << " does not exist." << endl;
      cout << "-----------------------------------------------------" << endl;
      return false;
    } 
  cout << "-----------------------------------------------------" << endl;
  cout << "Checking " << fileName << " (" << fileStat.st_size << " bytes)." << endl;
  if( fileStat.st_size == 0 )
    {
      cout << "file size = 0 !!" << endl << endl;
      return false;
    }
  cout << "-----------------------------------------------------" << endl;
  return true;
}


//Load the data
void LoadAugerData(string filename,  vector< double >& vId, vector< double >& vTheta, vector< double >& vl,  vector< double >& vb,  vector< double >& vE, vector< double >& vexp, string coord, bool verbose)
{
	//Following numbers depend on the format of the file!!
	unsigned int nentries_per_line = 9;
	//Yr	JD	UTC	Th	Ph	RA	Dec	E	Expo
	//0		1	2	3	4	5	6	7	8
	
	//coord can be gal, eq or sgal 	

	//Clear the vectors to file
	vl.clear();
	vb.clear();
	vTheta.clear();
	vId.clear();
	vE.clear();
	vexp.clear();

	//Load data from file
	vector< tuple<double, double, double, double, double, double> > vdata;
	std::ifstream myfile;	
	myfile.open(filename.c_str());
	if(myfile.is_open()){
		//loops on the lines
		while(myfile.good()){
			string line;
			getline (myfile,line);
			std::stringstream ss(line);
			vector< double > vbuf;
			double buf;;
			while(ss>>buf) vbuf.push_back(buf);
			if((vbuf.size()==nentries_per_line) and (vbuf[0]>2000) and (vbuf[1]<2021) ){
				double ra = vbuf[5];
				double dec = vbuf[6];
				if(!coord.compare("eq"))vdata.push_back(make_tuple(vbuf[7], vbuf[2], vbuf[3], ra, dec, vbuf[8]));//E, UTC, Theta, Ra, Dec, exp
				else if(!coord.compare("sgal")){
				double l=0, b=0;
				double sl = 0.;
    				double sb = 0.;
				radec2gal(ra/360*24, dec, &l, &b);
				gal2Sgal(l, b, sl, sb);
				vdata.push_back(make_tuple(vbuf[7], vbuf[2], vbuf[3], sl, sb, vbuf[8]));//E, UTC, Theta, supergal l, supergal b, exp
				}
				else {
				double l=0, b=0;
				radec2gal(ra/360*24, dec, &l, &b);
				vdata.push_back(make_tuple(vbuf[7], vbuf[2], vbuf[3], l, b, vbuf[8]));//E, UTC, Theta,  l,  b, exp
				}			}
			else if(verbose)  cout<<"Line unread in event file: "<<line<< endl;
		}
		myfile.close();
	}
	else cout<<"Could not open file: "<<filename<< endl;
	
	//Sort by increasing energy
	sort(vdata.begin(), vdata.end());
	for(unsigned int i=0; i<vdata.size(); i++){
		vE.push_back(get<0>(vdata[i]));
		vId.push_back(get<1>(vdata[i]));
		vTheta.push_back(get<2>(vdata[i]));
		vl.push_back(get<3>(vdata[i]));
		vb.push_back(get<4>(vdata[i])); 	
		vexp.push_back(get<5>(vdata[i])); 	
	}
}


// Code from Xavier Bertou - Temp fix by JB stderr -> stdout (because of ugly fix to remove message from healpix)
TProgressBar::TProgressBar() {Zero();}



void TProgressBar::Zero()
{
  fCurrentPercentage = 0;
  fReached = 0;
  fBegin = 0;
  fEnd = 0;
}



void TProgressBar::InitPercent()
{
  fprintf(stdout,"[                                                  ]\r[");
  fCurrentPercentage = 0;
  fReached = 0;
}



void TProgressBar::EndPercent() const
{
  fprintf(stdout,"\r[##################################################]\n");
}



void TProgressBar::PrintPercent(unsigned int value)
{
  fReached=value;
  unsigned int newper=fReached*100/(fEnd-fBegin);
  // step of 2%
  while( newper>fCurrentPercentage+2 )
    {
      fprintf(stdout,"#");
      fflush(stdout);
      fCurrentPercentage+=2;
    }
}



template<typename T> vector< vector<unsigned int> > GetContiguousPoints(unsigned int size, const T *xx, T ignore)
{
  unsigned int xi(0), xstart;
  vector< vector<unsigned int> > good;
  // The idea is to find in xx contiguous portions of data.These portions are separated by the value ignore. 
  // For instance, for : 0 1 2 3 -1 4 5 6 -1 7 8 9
  // the function return 3 vectors of indices corresponding to the elements 0 1 2 3 then 4 5 6 and finally 7 8 9.
  // This function is used to correctly draw the isolongitude and isolatitude lines. Indeed, if you try to connect 
  // points lying inside your map with points lying outside then you will obtain large staright lines. This is what 
  // you will if you try to connect the -90 latitude with the 90 latitude
  // find first == ignore  
  while( xx[xi] != ignore && xi < size ) xi++;
  if( xi == size )
    {
      // only good data
      vector<unsigned int> tmp(size);
      for( unsigned int i = 0; i < size; i++ ) tmp[i] = i;
      // then return all indices as they are in the original vector
      good.push_back(tmp);
      return good;
    }
  xstart = xi;
  xi++; // skip detected ignore point
  while( xi != xstart ) // have to scan the complete vector until we are back to first ignore
    {
      // find first != ignore
      while( xx[xi] == ignore && xi != xstart ) xi=(xi+1+size)%size;
      // find next == ignore
      vector<unsigned int> contiguous;
      while( xx[xi] != ignore && xi != xstart )
        {
          contiguous.push_back(xi);
          xi = (xi+1+size)%size;
        }
      if( contiguous.size() ) good.push_back(contiguous);
    }
  return good;
}


int * GetPalette(int MaxPretty)
{
  TColor *color;
  TColor *col1 = gROOT->GetColor(1);
  float saturation = 1;
  float lightness = 0.5;
  float MaxHue = 280;
  float MinHue = 0;
  float hue;
  float r, g, b;
  Int_t *palette = new Int_t[MaxPretty];
  for(Int_t i=0; i<MaxPretty; i++)
    {
      hue = MaxHue-(i+1)*((MaxHue-MinHue)/MaxPretty);
      color = gROOT->GetColor(51+i);
      if ( color )
        {
          color->HLStoRGB(hue,lightness,saturation,r,g,b);
          color->SetRGB(r,g,b);
        }
      else
        {
          col1->HLStoRGB(hue,lightness,saturation,r,g,b);
          color = new TColor(51+i,r,g,b);
        }
      palette[i] = 51+i;
    }
  return palette;
}



//! instanciation of templates
template vector< vector<unsigned int> > GetContiguousPoints(unsigned int size, const double * xx, double ignore_value);
template vector< vector<unsigned int> > GetContiguousPoints(unsigned int size, const float * xx, float ignore_value);
template vector< vector<unsigned int> > GetContiguousPoints(unsigned int size, const int * xx, int ignore_value);
