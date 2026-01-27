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
#include "maptools.h"
#include "projmap.h"
#include "common.h"
#include "STClibrary.h"
#include "healpixmap.h"

// ROOT
#include "TROOT.h"
#include "TH2D.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TObjArray.h"

// UTILITIES CATALOG ANALYSIS
#include "utils.h"

using namespace std;

//Global variables for the likelihood
vector<double> viso_GLOB, vsigma_theta_GLOB;
vector< vector<double> > vmodel_GLOB;
int nEvents;
/// Data Loading //////////////////////////////////////////////////////////////////////////////////////

void Trim(double Eth, vector< double >& vE, vector< double >& vZen, vector< double >& vl, vector< double >& vb, vector< double >& vexp){
	while(vE.front()<Eth){
		if(vE.front()>vE[1]) cout<<"EVENTS ARE NOT ORDERED IN ENERGY!!!! WRONG SCAN"<< endl;
		vE.erase(vE.begin());
		vl.erase(vl.begin());
		vb.erase(vb.begin());
		vZen.erase(vZen.begin());
		vexp.erase(vexp.begin());
	}
	if(vE.size()!=vl.size() || vE.size()!=vb.size() || vE.size()!=vZen.size() || vE.size()!=vexp.size())  cout<<"Problem in trimming!!!!"<< endl;
}

//Load the source model
void LoadModel(string filename, vector< string >& vSrcName, vector< double >& vl, vector< double >& vb, vector< double >& vw, bool verbose)
{
	//Following numbers depend on the format of the file!!
	unsigned int nentries_per_line = 4;

	//Clear the vectors to file
	vSrcName.clear();
	vl.clear();
	vb.clear();
	vw.clear();

	//Load data from file
	ifstream myfile;	
	myfile.open(filename.c_str());
	if(myfile.is_open()){
		//loops on the lines
		while(myfile.good()){
			string line;
			getline (myfile,line);
			stringstream ss(line);
			//get source name
			string srcName;
			ss>>srcName;
			//get source numbers
			vector< double > vbuf;
			double buf;;
			while(ss>>buf) vbuf.push_back(buf);
			///fill
			if(vbuf.size()==nentries_per_line-1){
				double ra = vbuf[0], dec = vbuf[1];
				double l=0, b=0;
				radec2gal(ra/360*24, dec, &l, &b);
				vSrcName.push_back(srcName);
				vl.push_back(l);
				vb.push_back(b);
				vw.push_back(vbuf[2]);
			}
			else if(verbose) cout<<"Line unread in event file: "<<line<<endl;
		}
		myfile.close();
	}
	else cout<<"Could not open file: "<<filename<<endl;
}



/// Likelihood ////////////////////////////////////////////////////////////////////////////////////////
//Smoothing function
double smoothing_fun(double theta, double sig_theta)
{
	//Fisher - Von Mises, e.g. arXiv:0703168v1
	double theta_rad = theta*DTOR;
	double K = 0.5/(1-TMath::Cos(sig_theta*DTOR));
	return exp(K*cos(theta_rad));//the normalization does not matter for the smoothing
}

//Smoothing vectors
vector< vector < double > > LoadSmoothingVectors(double theta){
	double theta_max_smooth = 180.;
	unsigned int ntheta_smooth = int(theta_max_smooth/0.2);
	vector< double > vtheta, vlobe;
	for(unsigned int j=0; j<ntheta_smooth; j++){
		vtheta.push_back(j*theta_max_smooth/(ntheta_smooth-1.));
		vlobe.push_back(smoothing_fun(vtheta[j],theta));
	}

	vector< vector < double > > vsmoothing;
	vsmoothing.push_back(vtheta);
	vsmoothing.push_back(vlobe);
	
	return vsmoothing;
}


//Load a source map for a given smoothing radius
THealpixMap LoadSrcMap(int nSide, double theta, vector<double> vl, vector<double> vb, vector<double> vw)
{
	THealpixMap inputMap(nSide,'G');
	for(unsigned int i=0; i<vl.size(); i++){
		double ip = inputMap.Ip(vl[i], vb[i]);
		inputMap[ip]+=vw[i];
	}

	vector< vector< double > > vsmoothing = LoadSmoothingVectors(theta);

	inputMap = inputMap.Filter(vsmoothing[0], vsmoothing[1]);

	inputMap /= inputMap.Max();
	
	return inputMap;
}
	
vector< THealpixMap > LoadSmoothedModelMaps(int nSide, vector<double> vl_model, vector<double> vb_model, vector<double> vw_model, bool low_resolution, bool mute)
{
	if(mute) freopen("/dev/null", "w", stderr);

	vector< THealpixMap > vsmoothedMap;
	
	//steps in sigma_theta
	double sig_theta_min = 3., sig_theta_max = 40., dsig_theta = 1.;
	//to speed up the results obtained with simulations
	//low-resolution results degraded by at most ~0.5 TS unit
	if(low_resolution){
		sig_theta_min = 3.;
		sig_theta_max = 39.;
		dsig_theta = 9.;
	}
	vsigma_theta_GLOB.clear();
	vsigma_theta_GLOB.push_back(sig_theta_min);
	while(vsigma_theta_GLOB.back()<sig_theta_max) vsigma_theta_GLOB.push_back(vsigma_theta_GLOB.back()+dsig_theta);
	
	for(unsigned int i=0; i<vsigma_theta_GLOB.size(); i++) vsmoothedMap.push_back(LoadSrcMap(nSide, vsigma_theta_GLOB[i], vl_model, vb_model, vw_model));

	return vsmoothedMap;
}

//Used in interpolations to find the element just before x
int utils_GetIndex(double x, vector< double > vx, bool verbose)
{
	vector< double >::iterator it_lower_bound = lower_bound (vx.begin(), vx.end(), x);//first element larger than x
	if(it_lower_bound==vx.begin()){
		it_lower_bound++;
		double relative_distance = (vx.front()-x)/(vx.back()-vx.front())*(vx.size()-1.);
		if(relative_distance>0.1 && verbose) cout<<"Value out of range (>xmax) in interpolation - extrapolation performed, distance to xmin / step = "<<relative_distance*100<<"%"<< endl;
	}
	else if(it_lower_bound==vx.end()){
		double relative_distance = (x-vx.back())/(vx.back()-vx.front())*(vx.size()-1.);
		if(relative_distance>0.1 && verbose) cout<<"Value out of range (>xmax) in interpolation - extrapolation performed, distance to xmax / step = "<<relative_distance*100<<"%"<< endl;
	}
	int iup = it_lower_bound-vx.begin();
	if(iup>=(int)vx.size()) iup=vx.size()-1.;//In case of off-bound events

	return iup;
}

//Interpolate a function = vector< vector < > >[i_x][i_y] at x
vector< double > utils_Interpolate(vector< vector<double> > vData, vector<double> vx, double x)
{
	vector<double> vInterpol;
	int iup = utils_GetIndex(x,vx);
	for(unsigned int j=0; j<vData[iup].size(); j++) vInterpol.push_back(vData[iup-1][j] + (vData[iup][j]-vData[iup-1][j])*(x-vx[iup-1])/(vx[iup]-vx[iup-1]));
	return vInterpol;
}	

//Compute the likelihood of a model map for a given event map 
double LogLikelihood(double theta, double alpha)
{
	vector < double > vmodel = utils_Interpolate(vmodel_GLOB, vsigma_theta_GLOB, theta);
	double logL=0.;
	for(unsigned int i=0; i<viso_GLOB.size(); i++){
		double nmodel =  alpha*vmodel[i] + (1.-alpha)*viso_GLOB[i];
		if(nmodel<=0)	return -1.E99;//default value in case of troubles
		else logL+=log(nmodel);
	}
	return logL;
}

//Compute the likelihood vector of a model map for each event 
vector<double> VectorTSevents(double theta, double alpha)
{
	vector < double > vmodel = utils_Interpolate(vmodel_GLOB, vsigma_theta_GLOB, theta);
	vector< double > vTS;
	for(unsigned int i=0; i<viso_GLOB.size(); i++){
		double nmodel =  alpha*vmodel[i] + (1.-alpha)*viso_GLOB[i];
		double nmodel0 = viso_GLOB[i];
		double TS = -1.E99;
		if(nmodel>0 && nmodel0>0) TS = 2.*log(nmodel/nmodel0);
		vTS.push_back(TS);
	}
	return vTS;
}

void Minuit_MLogL_All(int &, double *, double &f, double *par, int )
{
	f = -LogLikelihood(par[0], par[1]);
	return; 
}

//Find best likelihood vs theta
double FindMaxLogL_All(double& theta, double& dtheta, double& alpha, double& dalpha, bool fix_theta, bool fix_alpha)
{
	
	//Minuit Initialization
	Int_t npar = 2;
	TMinuit* ptMinuit = new TMinuit(npar);  //initialize TMinuit with a maximum of npar
	ptMinuit->SetPrintLevel(-1);//-1 no output - 1 std output
	double arglist[10];  int ierflg = 0;
	ptMinuit->mnexcm("SET NOW",arglist,1,ierflg);//No warnings
	ptMinuit->SetErrorDef(0.5);//1 for chi2, 0.5 for -logL
	ptMinuit->mnparm(0,"sig_theta",theta,dtheta,3.,40.,ierflg);//range for theta (pixels too small below 3 degrees)
	ptMinuit->mnparm(1,"alpha",alpha,dalpha,0.00,1.00,ierflg);
	if(fix_theta) ptMinuit->FixParameter(0);
	if(fix_alpha) ptMinuit->FixParameter(1);

	ptMinuit->SetFCN(Minuit_MLogL_All);
	int nmax_iterations = 100;
	ptMinuit->SetMaxIterations(nmax_iterations);
	arglist[0]=1;//1 std, 2 try to improve (slower)
	ptMinuit->mnexcm("SET STR",arglist,1,ierflg);	
	ptMinuit->mnexcm("SIMPLEX",0,0,ierflg);	

	double MLogLmin = 0.;
	double fedm, errdef;
	int npari, nparx, istat;
	ptMinuit->mnstat(MLogLmin,fedm, errdef, npari, nparx, istat);
	ptMinuit->GetParameter(0,theta,dtheta);
	ptMinuit->GetParameter(1,alpha,dalpha);
	delete ptMinuit;\
	return -MLogLmin;
}

//Returns the max TS and best-fit parameters for a given UHECR dataset
double TS_fixed_E(vector< double > vl, vector< double > vb, THealpixMap exposureMap, vector< THealpixMap > vsmoothedMap, double& sig_theta_bf, double& alpha_bf, double& logL_0, bool mute){
	
	if(mute) freopen("/dev/null", "w", stderr);
	
	unsigned int nevents = vl.size();
	unsigned int nSide = exposureMap.NSide();

	// Load the list of pixels with non-zero content and the event map for final plotting purpose	
	THealpixMap eventMap(nSide,'G');
	vector< unsigned int > vIpix;
	for(unsigned int i=0; i<nevents; i++){
		int pix = exposureMap.Ip(vl[i],vb[i]);
		vIpix.push_back(pix);
		eventMap[pix]+=1;
	}

	//Load the isotropic map
	THealpixMap modelMapIsotropic(nSide,'G');
	modelMapIsotropic=exposureMap*(nevents/exposureMap.Total());
	viso_GLOB.clear();
	viso_GLOB.resize(vIpix.size());
	for(unsigned int i=0; i<nevents; i++) viso_GLOB[i] = modelMapIsotropic[vIpix[i]]; 

	//Load the normalized smoothed maps
	vmodel_GLOB.clear();
	vmodel_GLOB.resize(vsigma_theta_GLOB.size());
	for(unsigned int i=0; i<vsmoothedMap.size(); i++){
		THealpixMap modelMap(nSide,'G');
		modelMap = vsmoothedMap[i]*exposureMap;
		modelMap=modelMap*(nevents/modelMap.Total());

		//Load the vector of event counts vs iPix
		vmodel_GLOB[i].clear();
		vmodel_GLOB[i].resize(vIpix.size());
		for(unsigned int j=0; j<vIpix.size(); j++) vmodel_GLOB[i][j] = modelMap[vIpix[j]]; 
	}

	//Likelihood for the data
	sig_theta_bf=10.;
	double dsig_theta_bf=5;
	bool fix_theta=false; 
	alpha_bf=0.1;
	double dalpha_bf=0.05;
	bool fix_alpha=false;
	double logL = FindMaxLogL_All(sig_theta_bf, dsig_theta_bf, alpha_bf, dalpha_bf, fix_theta, fix_alpha);

	if(logL<-1.E50 || isnan(logL)) cout<<"!!! Danger: computed chi2 is not right:"<<logL<<endl;
	logL_0 = LogLikelihood(sig_theta_bf, 0.);
	if(logL_0<-1.E50 || isnan(logL_0)) cout<<"!!! Danger: reference chi2 is not right:"<<logL_0<<endl;

	//Compute significance	
	double TS_data = 2.*(logL-logL_0);

	return TS_data;
}

//Returns the max TS at given parameters for a given UHECR dataset
double TS_fixed_E_Par(vector< double > vl, vector< double > vb, THealpixMap exposureMap, vector< THealpixMap > vsmoothedMap, double sig_theta_bf, double alpha_bf, bool mute){
	
	if(mute) freopen("/dev/null", "w", stderr);
	
	unsigned int nevents = vl.size();
	unsigned int nSide = exposureMap.NSide();

	// Load the list of pixels with non-zero content and the event map for final plotting purpose	
	THealpixMap eventMap(nSide,'G');
	vector< unsigned int > vIpix;
	for(unsigned int i=0; i<nevents; i++){
		int pix = exposureMap.Ip(vl[i],vb[i]);
		vIpix.push_back(pix);
		eventMap[pix]+=1;
	}

	//Load the isotropic map
	THealpixMap modelMapIsotropic(nSide,'G');
	modelMapIsotropic=exposureMap*(nevents/exposureMap.Total());
	viso_GLOB.clear();
	viso_GLOB.resize(vIpix.size());
	for(unsigned int i=0; i<nevents; i++) viso_GLOB[i] = modelMapIsotropic[vIpix[i]]; 

	//Load the normalized smoothed maps
	vmodel_GLOB.clear();
	vmodel_GLOB.resize(vsigma_theta_GLOB.size());
	for(unsigned int i=0; i<vsmoothedMap.size(); i++){
		THealpixMap modelMap(nSide,'G');
		modelMap = vsmoothedMap[i]*exposureMap;
		modelMap=modelMap*(nevents/modelMap.Total());

		//Load the vector of event counts vs iPix
		vmodel_GLOB[i].clear();
		vmodel_GLOB[i].resize(vIpix.size());
		for(unsigned int j=0; j<vIpix.size(); j++) vmodel_GLOB[i][j] = modelMap[vIpix[j]]; 
	}

	//Likelihood for the data
	double logL = LogLikelihood(sig_theta_bf, alpha_bf);
	if(logL<-1.E50 || isnan(logL)) cout<<"!!! Danger: computed chi2 is not right:"<<logL<<endl;
	double logL_0 = LogLikelihood(sig_theta_bf, 0.);
	if(logL_0<-1.E50 || isnan(logL_0)) cout<<"!!! Danger: reference chi2 is not right:"<<logL_0<<endl;

	//Compute significance	
	double TS_data = 2.*(logL-logL_0);

	return TS_data;
}

//Utilities////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Load a 2D histogram
TH2D *LoadHisto(string htitle, string ztitle, string xtitle, string ytitle, int nx, double xmin, double xmax, int ny, double ymin, double ymax, vector<double> vx,vector<double> vy,vector<double> vz, double default0)
{

	TH2D *h = new TH2D(htitle.c_str(),"",nx,xmin,xmax,ny,ymin,ymax);
	h->SetStats(0);
	h->SetDirectory(0);
	h->SetMinimum(0.);
	h->SetMaximum(*max_element(vz.begin(),vz.end()));
	h->GetXaxis()->SetTitle(xtitle.c_str());
	h->GetYaxis()->SetTitle(ytitle.c_str());
	h->GetZaxis()->SetTitle(ztitle.c_str());
	h->GetXaxis()->SetTitleOffset(1.2);
	h->GetYaxis()->SetTitleOffset(1.25);
	h->GetZaxis()->SetTitleOffset(1.0);
	h->GetXaxis()->CenterTitle();
	h->GetYaxis()->CenterTitle();
	h->GetZaxis()->CenterTitle();
	for(unsigned int i=0; i<vz.size(); i++){
		h->Fill(vx[i],vy[i],vz[i]);
		if(vx[i]==0) h->Fill(vx[i],vy[i],default0);
	}

	return h;
}

//Load a contour from a graph
TGraph* MyGraphContour(TH2D* hist, double contour_level)
{
	TH2D* mapProj = (TH2D*)hist->Clone();
	double list_contours[1] = {contour_level};	
	mapProj->SetContour(1,list_contours);

	TCanvas can_tmp;
	can_tmp.cd();
	mapProj->Draw("CONT Z LIST");		
	can_tmp.Update();

	TObjArray *plah = (TObjArray *)gROOT->GetListOfSpecials()->FindObject("contours");
	TList *list = (TList*)plah->At(0);
	TGraph* Gcont = (TGraph*)list->First();
	Gcont->SetLineWidth(2);

	return  (TGraph*)Gcont->Clone();
}

//Load the best-fit value & 1/2sigma contours based on TS map
vector< TGraph* > LoadContoursAndBestFit(TH2D *h, vector<double > vcont_level, double x0, double y0, bool verbose)
{
	vector< TGraph* > vG;
	vector< double > vxmin, vxmax, vymin, vymax;
	for(unsigned int i=0; i<vcont_level.size(); i++){
		vG.push_back(MyGraphContour(h,vcont_level[i]));
		double *x = vG.back()->GetX();
		double *y = vG.back()->GetY();
		double xmin=1, xmax=0, ymin=90, ymax=0;
		for(int j=0; j<vG.back()->GetN(); j++){
			if(x[j]<xmin) xmin=x[j];
			if(x[j]>xmax) xmax=x[j];
			if(y[j]<ymin) ymin=y[j];
			if(y[j]>ymax) ymax=y[j];
		}
		vxmin.push_back(xmin);
		vxmax.push_back(xmax);
		vymin.push_back(ymin);
		vymax.push_back(ymax);
	}

	double x[1] = {x0}, y[1] = {y0};
	TGraph *G_bf = new TGraph(1,x,y);
	G_bf->SetMarkerStyle(5);
	G_bf->SetMarkerSize(3.);
	vG.push_back(G_bf);

	if(verbose){
		cout<<"Best fit at: "<<endl;
		cout<<"alpha +/-1sigma  = "<<setprecision(3)<<100.*x0<<" -"<<100.*(x0-vxmin.front())<<" +"<<100.*(vxmax.front()-x0)<<" %"<<endl;
		cout<<"theta +/-1sigma = "<<setprecision(3)<<y0<<" -"<<y0-vymin.front()<<" +"<<vymax.front()-y0<<" deg"<<endl;
	}

	return vG;
}

//Load the best-fit value & 1sigma uncertainties
vector< double > LoadBestFitParam(TH2D *h, double cont_level, double x0, double y0, bool verbose)
{
	TGraph* G = MyGraphContour(h,cont_level);
	double *x = G->GetX();
	double *y = G->GetY();
	double xmin=1, xmax=0, ymin=90, ymax=0;
	for(int j=0; j<G->GetN(); j++){
		if(x[j]<xmin) xmin=x[j];
		if(x[j]>xmax) xmax=x[j];
		if(y[j]<ymin) ymin=y[j];
		if(y[j]>ymax) ymax=y[j];
	}
	double dxm = x0-xmin;
	if(dxm<0) dxm = x0; 
	vector< double > vres = {x0, dxm, xmax-x0, y0, y0-ymin, ymax-y0};

	return vres;
}
