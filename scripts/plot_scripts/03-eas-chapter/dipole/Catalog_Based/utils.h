///////////////////////////////////////////////////////
//                                                   //
//   Dev. by Jonathan Biteau (jbiteau.pro@gmail.com) //
//                       last edit: 2021-11-25       //
///////////////////////////////////////////////////////

// C/C++ classics
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <sys/stat.h>
#include <numeric>
#include <time.h>
#include <tuple>
#include <utility>

// Toolkit
#include "projmap.h"
#include "healpixmap.h"

// ROOT
#include "TROOT.h"
#include "TH2D.h"
#include "TGraph.h"

using namespace std;


//Global variable for the maps
extern int sizeX;
extern int sizeY;
extern double decLimit;//max to plot all the map
extern double longStep;
extern double latStep;
//Cosmetics
extern double max_marker_size;
extern int marker_style;
extern int marker_color, coord_color;
extern double radius_write_source;

//Global variables needed by Minuit
extern vector<double> viso_GLOB, vsigma_theta_GLOB;
extern vector< vector<double> > vmodel_GLOB;

/// Data Loading //////////////////////////////////////////////////////////////////////////////////////

//Remove first elements - ASSUMES events are energy ordered!!!!!
void Trim(double Eth, vector< double >& vE, vector< double >& vZen, vector< double >& vl, vector< double >& vb, vector< double >& vexp);

//Load the source model
void LoadModel(string filename, vector< string >& vSrcName, vector< double >& vl, vector< double >& vb, vector< double >& vw, bool verbose=false);


/// Likelihood ////////////////////////////////////////////////////////////////////////////////////////

//Smoothing function
double smoothing_fun(double theta, double sig_theta);

//Smoothing vectors
vector< vector < double > > LoadSmoothingVectors(double theta);

//Load a source map for a given smoothing radius
THealpixMap LoadSrcMap(int nSide, double theta, vector<double> vl, vector<double> vb, vector<double> vw);
	
vector< THealpixMap > LoadSmoothedModelMaps(int nSide, vector<double> vl_model, vector<double> vb_model, vector<double> vw_model, bool low_resolution=false, bool mute=true);

//Used in interpolations to find the element just before x
int utils_GetIndex(double x, vector< double > vx, bool verbose = false);

//Interpolate a function = vector< vector < > >[i_x][i_y] at x
vector< double > utils_Interpolate(vector< vector<double> > vData, vector<double> vx, double x);

//Compute the likelihood of a model map for a given event map 
double LogLikelihood(double theta, double alpha);

void Minuit_MLogL_All(int &, double *, double &f, double *par, int );

//Find best likelihood vs theta
double FindMaxLogL_All(double& theta, double& dtheta, double& alpha, double& dalpha, bool fix_theta=false, bool fix_alpha=false);

//Returns the max TS and best-fit parameters for a given UHECR dataset
double TS_fixed_E(vector< double > vl, vector< double > vb, THealpixMap exposureMap, vector< THealpixMap > vsmoothedMap, double& sig_theta_bf, double& alpha_bf, double& logL_0, bool mute=true);

//Compute the likelihood vector of a model map for each event 
vector<double> VectorTSevents(double theta, double alpha);

//Returns the TS for fixed parameters and a given UHECR dataset
double TS_fixed_E_Par(vector< double > vl, vector< double > vb, THealpixMap exposureMap, vector< THealpixMap > vsmoothedMap, double sig_theta_bf, double alpha_bf, bool mute=true);

//Utilities////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Load a 2D histogram
TH2D *LoadHisto(string htitle, string ztitle, string xtitle, string ytitle, int nx, double xmin, double xmax, int ny, double ymin, double ymax, vector<double> vx,vector<double> vy,vector<double> vz, double default0 = 1E-5);

//Load a contour from a graph
TGraph* MyGraphContour(TH2D* hist, double contour_level);

//Load the best-fit value & 1/2sigma contours based on TS map
vector< TGraph* > LoadContoursAndBestFit(TH2D *h, vector<double > vcont_level, double x0, double y0, bool verbose=true);

//Load the best-fit value & 1sigma uncertainties
vector< double > LoadBestFitParam(TH2D *h, double cont_level, double x0, double y0, bool verbose=true);
