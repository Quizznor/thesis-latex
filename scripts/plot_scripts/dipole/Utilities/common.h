#ifndef _COMMON_H
#define _COMMON_H

#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <TH1F.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TColor.h>

#include "STCconstants.h"

using namespace std;


/*
  Uses the kSTC namespace defined in the STCoordinates package. It contains constants such as degrees to radians 
  conversion, usefull sky-transformations constants and so on.
*/
using namespace kSTC;



//! Progress bar in the terminal. 
class TProgressBar
{
public:
  //! Constructor.
  TProgressBar();

  //! Initializes everything.
  void Zero();

  //! First value in the loop (0\%).
  unsigned int fBegin;

  //! Last value in the loop (100\%).
  unsigned int fEnd;

  //! Current value of the loop.
  unsigned int fReached;

  //! Current corresponding percentage value.
  unsigned int fCurrentPercentage;

  //! Initializes stderr.
  void InitPercent();

  //! Dump last stderr (full bar).
  void EndPercent() const;

  //! Print current percentage to stderr.
  void PrintPercent(unsigned int current);
};

/*!
  Given a input array of data containing good data and bad data (with 
  value ignore), detect contiguous points (not separated by the ignore 
  value) and return the list of the indexes of contiguous points. Example : 
  data = [10,11,12,13,-1,-1,14,15,16,-1,-1] with ignore = -1 returns [6,7,8] 
  and [0,1,2,3], index lists of good data.
*/
template<typename T> vector< vector<unsigned int> > GetContiguousPoints(unsigned int size, const T *input, T ignore_value);



//! Get the ROOT color table.
int * GetPalette(int maxpretty);

//Reads the data file.
void LoadAugerData(string filename,  vector< double >& vId, vector< double >& vTheta, vector< double >& vl,  vector< double >& vb,  vector< double >& vE, vector< double >& vexp, string coord ="gal", bool verbose=false); // coord can be gal, eq or sgal

#endif


