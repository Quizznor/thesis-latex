///////////////////////////////////////////////////////
//                                                   //
//   Dev. by Jonathan Biteau (biteau@in2p3.fr)       //
//                       last edit: 2022-03-14       //
///////////////////////////////////////////////////////

// C/C++ classics
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <tuple>

// Toolkit
#include "healpixmap.h"

// ROOT
#include "TRint.h"
#include "TROOT.h"
#include "TH2D.h"
#include "TMath.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TF1.h"

// UTILITIES CATALOG ANALYSIS
#include "utils.h"

#ifdef gcc323
char* operator+( streampos&, char* );
#endif

using namespace std;

//Progress bar
static TProgressBar gProgB;

//Global variables defining the data, catalogs and composition scenarios

//Data///////////////////
string InputData_seventfile = "../Data/AugerApJS2022_Yr_JD_UTC_Th_Ph_RA_Dec_E_Expo";
string InputData_sexpofitsfile = "../Data/exposure.fits";	

//Catalogs///////////////
vector< string > InputData_vssources = {"Starburst galaxies (radio)", "Jetted AGN (#gamma-rays)", "All AGN (hard X-rays)", "Galaxies > 1 Mpc (IR)"};//Associated display name
vector< string > InputData_vsfile = {"sbg", "fermi_lat_agn", "swift_bat_agn", "2mass_hyperleda"} ;//Associated folder, either: sbg, fermi_lat_agn, swift_bat_agn, 2mass_hyperleda

//Composition////////////
vector< string > InputData_vscenarios = {"EPO1st"};
vector< string > InputData_vname_scenarios = {"EPOS-LHC 1st minimum"};

//Cosmetics//////////////
Color_t col[4] = {kAzure+6, kMagenta+1, kOrange+8, kSpring-8};
int lstyle[4] = {1, 7, 2, 9};
int lwidth[4] = {2, 3, 3, 2};
int mstyle[4] = {21, 25, 20, 27};


void Usage(string sinput)
{
	cout << endl;
	cout << " Synopsis : " << endl;
	cout << sinput << " -h, --help to obtain this message" << endl;
	cout << " Description :" << endl;
	cout << " Simulate isotropic events and scan as a function of threshold energy"<<endl
		<< " -- that's all folks! -- "
		<< endl<< endl;

	exit(0);
}

/// Isotropic datasets: Eth, l, b//////////////////////////////////////////////////////////////////////
void IsotropicDataset(vector< unsigned int>& vEsim, vector< double >& vlsim, vector< double >& vbsim, vector< unsigned int > vnevts,  vector< unsigned int > vEth, THealpixMap exposureMap, tuple< TRandom3, TRandom3, TRandom3 >& vr){
	
	//Clean
	vEsim.clear();
	vlsim.clear();
	vbsim.clear();

	//Random variables for generation on the sphere
	int nmaxsimu=1E6;
	TRandom3 ru(get<0>(vr).Integer(nmaxsimu)), rv(get<1>(vr).Integer(nmaxsimu));//for long, lat, picking
	TRandom3 rp(get<2>(vr).Integer(nmaxsimu));//for rejection/acceptance

	//Load the simulated evts: same number of events above decreasing Eth
	double expo_max = exposureMap.Max();
	unsigned int count = 0;
	vector< tuple<unsigned int, double, double> > vevts;
	for(int k=vEth.size()-1; k>=0; k--){
		while(count<vnevts[k]){
			//Uniform picking on the sphere
			double phi = 2*M_PI*ru.Uniform();// in [0;2pi]
			double theta = TMath::ACos(2.*rv.Uniform()-1.);// in [0;pi]
			double l = (phi-M_PI)*RTOD;// in [-180;180]
			double b = (theta-0.5*M_PI)*RTOD;// in [-90;90]

			//Acceptance/rejection criterion
			long pix = exposureMap.Ip(l,b);
			if(rp.Uniform()<=exposureMap[pix]/expo_max){
				count++;
				vevts.push_back(make_tuple(vEth[k], l, b));
			}
		}
	}

	//Sort the simulated events by increasing energy
	sort(vevts.begin(), vevts.end());		
	for(unsigned int i=0; i<vevts.size(); i++){
		vEsim.push_back(get<0>(vevts[i]));
		vlsim.push_back(get<1>(vevts[i]));
		vbsim.push_back(get<2>(vevts[i])); 	
	}
}
	
void Trim(unsigned int Eth, vector< unsigned int >& vE, vector< double >& vl, vector< double >& vb){
	while(vE.front()<Eth){
		if(vE.front()>vE[1]) cout<<"EVENTS ARE NOT ORDERED IN ENERGY!!!! WRONG SCAN"<< endl;
		vE.erase(vE.begin());
		vl.erase(vl.begin());
		vb.erase(vb.begin());
	}
	if(vE.size()!=vl.size() || vE.size()!=vb.size())  cout<<"Problem in trimming!!!!"<< endl;
}
	
/// Max TS above all thresholds///////////////////////////////////////////////////////////////////////
double LoadTS(vector< unsigned int > vEsim, vector< double > vlsim, vector< double > vbsim, vector< unsigned int > vEth, THealpixMap exposureMap, vector< THealpixMap > vsmoothedMap){
	double TS_max = 0;
	for(unsigned int i=0; i<vEth.size(); i++){
		Trim(vEth[i],vEsim, vlsim, vbsim);
		double theta_bf, alpha_bf, logL_0;
		double TS_simplex = TS_fixed_E(vlsim,vbsim, exposureMap, vsmoothedMap, theta_bf, alpha_bf, logL_0);
		if(TS_simplex>TS_max) TS_max = TS_simplex;	
	}
	return TS_max;
}


int count_above(TH1D* hTS, int i_first_count_bin){
	int n=0;
	for(int i=i_first_count_bin; i<=hTS->GetNbinsX(); i++) n+=hTS->GetBinContent(i);
	return n;
}

//Downward curved parabolic fit in lin-log
double FitNorm(TGraphErrors *G, double xmin = 0, double xmax = 10){

	TF1* f = new TF1("f","[0]*TMath::Exp([1]*x+[2]*x*x)",xmin, xmax);
	f->SetParameters(1,-0.5,0);
	f->SetParLimits(2,-1,0);//prevents upward going curve
	G->Fit(f,"","",xmin, xmax);
	double p0 = f->GetParameter(0);
	
	delete f;

	return p0;
}

//p-value vs TS
TGraphErrors* LoadGraphTS(TH1D* hTS){
	int nTS_tot = count_above(hTS, 1);
	vector< double > vTS, vpval, vdpval;
	for(int i=1; i<=hTS->GetNbinsX(); i++){
		double TS = hTS->GetBinCenter(i);
		int ncounts_above = count_above(hTS, i);
		
		if(TS>=0 && ncounts_above>0){
			double pval = ncounts_above*1./nTS_tot;
			vTS.push_back(TS);
			vpval.push_back(pval);
			vdpval.push_back(pval/sqrt(1.*ncounts_above));//Note that uncertainites are correlated
		}
	}
	
	//Load the graphs - normalized through a parabolic fit in lin-log
	TGraphErrors *Gpval_no_norm = new TGraphErrors(vTS.size(), &vTS[0], &vpval[0], 0, &vdpval[0]);
	double norm = FitNorm(Gpval_no_norm);
	for(unsigned int i=0; i<vTS.size(); i++){
		vpval[i]/=norm;
		vdpval[i]/=norm;
	}
	TGraphErrors *Gpval = new TGraphErrors(vTS.size(), &vTS[0], &vpval[0], 0, &vdpval[0]);
	
	return Gpval;
}

//Compute penalty factor = post-scan p-value / expected p-value from chi2 with 2 d.o.f.
TGraphErrors* LoadPenalty(TGraphErrors* Gpval){
	TF1 *f = new TF1("f","TMath::Prob(x,2)",0.,50);	//expected distribution w/o E-scan

	vector< double > vTSpen, vPen, vdPen;
	for(int i=0; i<Gpval->GetN(); i++){
		double TS = Gpval->GetX()[i];
		double norm = f->Eval(TS);
		vTSpen.push_back(TS);
		vPen.push_back(Gpval->GetY()[i]/norm);
		vdPen.push_back(Gpval->GetErrorY(i)/norm);		
	}
	TGraphErrors *G = new TGraphErrors(vTSpen.size(), &vTSpen[0], &vPen[0], 0, &vdPen[0]);
		G->SetMarkerStyle(Gpval->GetMarkerStyle());
		G->SetMarkerSize(Gpval->GetMarkerSize());
		G->SetMarkerColor(Gpval->GetMarkerColor());
		G->SetLineColor(Gpval->GetLineColor());

	return G;
}

/// Main //////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
	////////////////////////////////////////////////////////////////////////////
	//                                                                        //
	//                            Initialization                              //
	//                                                                        //
	////////////////////////////////////////////////////////////////////////////

	// Command line
	if(argc != 1) Usage(argv[0]);

	// Initialize root
	int fargc = 1;
	TRint *rint = NULL;
	rint = new TRint("Test plot", &fargc, argv);
	gROOT->SetStyle("Plain");
	gStyle->SetTitleFont(30,"TITLE");

	//Load the data to determine nevts
	vector< double > vId, vZen, vl, vb, vE, vexp;
	LoadAugerData(InputData_seventfile , vId, vZen, vl, vb, vE, vexp);

	//Range over which the fit is ran: E > 32...80 EeV
	vector< unsigned int > vEth, vnevts;
	for(unsigned int Eth=32; Eth<81; Eth++){
		//Trim the data below threshold
		Trim(Eth,vE, vZen, vl,vb, vexp);
		
		vEth.push_back(Eth);
		vnevts.push_back(vE.size());
	}	
	
	//Exposure
	THealpixMap exposureMap((char*)InputData_sexpofitsfile.c_str());
	exposureMap.SetCoordSys('G');
	unsigned int nSide = exposureMap.NSide();
		
	//Composition model
	int compo_model = 0;//by default EPOS-LHC composition model
	
	//Degraded binning for model computation
	bool low_resolution = true;//9deg steps instead of 1deg
	
	//Load the model - take that above 40 EeV for fast computation (note: expected result indep. from model)
	vector < vector< THealpixMap > > vsmoothedMap;
	for(unsigned cat = 0; cat<InputData_vssources.size(); cat++){
		string ssrc_file = InputData_vsfile[cat];		
		stringstream ssfile_compo;
		ssfile_compo<<"Catalogs/ModelsUHECR/"<<ssrc_file<<"/"<<ssrc_file;
		if(InputData_vscenarios[compo_model].length()>1) ssfile_compo<<"_"<<InputData_vscenarios[compo_model];
		ssfile_compo<<"_threshold40";	
				
		vector< string > vSrcName;
		vector< double > vl_model, vb_model, vw_model;
		LoadModel(ssfile_compo.str(),vSrcName, vl_model, vb_model, vw_model);	
		
		vsmoothedMap.push_back(LoadSmoothedModelMaps(nSide, vl_model, vb_model, vw_model, low_resolution));
	}
	
	////////////////////////////////////////////////////////////////////////////
	//                                                                        //
	//                            Simulate and fit                            //
	//                                                                        //
	////////////////////////////////////////////////////////////////////////////

	// Storage
	vector< TH1D* > vhTS;
	double xmin = 2., xmax = 50.;
	double bin_size = 1;
	int nbins = (int)(xmax-xmin)/bin_size;
	for(unsigned cat = 0; cat<InputData_vssources.size(); cat++){
		TH1D* hTS = new TH1D(("hTS"+to_string(cat)).c_str(),"",nbins,xmin,xmax);
		vhTS.push_back(hTS);
	}

	// Random seeds
	unsigned int nseed = 123;	
	TRandom3 r0(nseed+0), r1(nseed+1), r2(nseed+2);
	tuple< TRandom3, TRandom3, TRandom3 > vr = make_tuple(r0, r1, r2);

	// Simulations
	unsigned int nsims = 10000;//10,000 ~ 20h on a single standard CPU
	
	cout<<"########### Starting "<<nsims<<" simulations ###########"<<endl;
	gProgB.Zero();
	gProgB.fBegin = 0;
	gProgB.fEnd = nsims;
	gProgB.InitPercent();

	vector< unsigned int> vEsim(vnevts[0]);
	vector< double > vlsim(vnevts[0]), vbsim(vnevts[0]);
	for(unsigned int i=0; i<nsims; i++){
		gProgB.PrintPercent(i);
		IsotropicDataset(vEsim, vlsim, vbsim, vnevts, vEth, exposureMap, vr);
		for(unsigned cat = 0; cat<InputData_vssources.size(); cat++){
			vhTS[cat]->Fill(LoadTS(vEsim, vlsim, vbsim, vEth, exposureMap, vsmoothedMap[cat]));	
		}
	}
	gProgB.EndPercent();	

	////////////////////////////////////////////////////////////////////////////
	//                                                                        //
	//                            p-val vs TS                                 //
	//                                                                        //
	////////////////////////////////////////////////////////////////////////////
	vector< TGraphErrors* > vG_pval, vGpenalty;
	for(unsigned cat = 0; cat<vhTS.size(); cat++){
		vG_pval.push_back(LoadGraphTS(vhTS[cat]));
		vG_pval.back()->SetMarkerStyle(mstyle[cat]);
		vG_pval.back()->SetMarkerSize(1.);
		vG_pval.back()->SetMarkerColor(col[cat]+1);
		vG_pval.back()->SetLineColor(col[cat]+1);
		
		vGpenalty.push_back(LoadPenalty(vG_pval.back()));
	}
	
	//Plot pval vs TS
	TH1D *hGraph0 = new TH1D("hGraph0","",100,0.,30.);
		hGraph0->SetStats(0);
		hGraph0->SetDirectory(0);
		hGraph0->SetMaximum(1.5);
		hGraph0->SetMinimum(1E-6);
		hGraph0->GetXaxis()->SetTitle("Test Statistics");
		hGraph0->GetYaxis()->SetTitle("p-value");
		hGraph0->GetXaxis()->SetTitleOffset(1.2);
		hGraph0->GetYaxis()->SetTitleOffset(1.3);
		hGraph0->GetXaxis()->CenterTitle();
		hGraph0->GetYaxis()->CenterTitle();	

	TLegend *leg_TS = new TLegend(0.50,0.85-0.05*(InputData_vssources.size()+1),0.88,0.85);
		leg_TS->SetLineColor(kWhite); 
		leg_TS->SetFillColor(kWhite);
		leg_TS->SetMargin(0.25); 
		leg_TS->SetTextSize(0.04); 
		for(unsigned int i=0; i<InputData_vssources.size(); i++){
			stringstream catname;
			catname<<"#color["<<col[i]+1<<"]{"<<InputData_vssources[i]<<"}";
			leg_TS->AddEntry(vG_pval[i],("#bf{"+catname.str()+"}").c_str(),"lp");
		}
		TF1 *fchi2 = new TF1("fchi2","TMath::Prob(x,2)",0.,50);	//expected distribution w/o E-scan
		fchi2->SetLineWidth(2);
		leg_TS->AddEntry(fchi2,"#bf{SDF(#chi^{2}_{2})}","l");		
				
	TCanvas *c0 = new TCanvas("p-value","p-value");	
		c0->SetLeftMargin(0.1);
		c0->SetRightMargin(0.1);
		c0->SetTopMargin(0.1);
		c0->SetBottomMargin(0.1);	
		c0->SetTickx();
		c0->SetTicky();
		c0->SetLogy();
		
		hGraph0->Draw();
		leg_TS->Draw();
		fchi2->Draw("samel");
		for(unsigned int i=0; i<vG_pval.size(); i++) vG_pval[i]->Draw("same pz");
		
	////////////////////////////////////////////////////////////////////////////
	//                                                                        //
	//                     Plot penalty vs TS                                 //
	//                                                                        //
	////////////////////////////////////////////////////////////////////////////

	//Plot penalty vs TS
	TH1D *hGraph1 = new TH1D("hGraph1","",100,0.,30.);
		hGraph1->SetStats(0);
		hGraph1->SetDirectory(0);
		hGraph1->SetMaximum(15.);
		hGraph1->SetMinimum(0);
		hGraph1->GetXaxis()->SetTitle("Test Statistics");
		hGraph1->GetYaxis()->SetTitle("Penalty");
		hGraph1->GetXaxis()->SetTitleOffset(1.2);
		hGraph1->GetYaxis()->SetTitleOffset(1.3);
		hGraph1->GetXaxis()->CenterTitle();
		hGraph1->GetYaxis()->CenterTitle();

	TLegend *leg1 = new TLegend(0.12,0.85-0.06*InputData_vssources.size(),0.50,0.85);
		leg1->SetLineColor(kWhite); 
		leg1->SetFillColor(kWhite);
		leg1->SetMargin(0.3); 
		leg1->SetTextSize(hGraph1->GetXaxis()->GetTitleSize());
		for(unsigned int j=0; j<vGpenalty.size(); j++){
			stringstream catname;
			catname<<"#color["<<col[j]+1<<"]{"<<InputData_vssources[j]<<"}";
			leg1->AddEntry(vGpenalty[j],("#bf{"+catname.str()+"}").c_str(),"lp");
		}

	TCanvas *cgraph1 = new TCanvas("penalty","penalty");
		cgraph1->SetLeftMargin(0.10);
		cgraph1->SetTopMargin(0.10);
		cgraph1->SetRightMargin(0.10);
		c0->SetBottomMargin(0.1);			
		cgraph1->SetTickx();
		cgraph1->SetTicky();
		
		hGraph1->Draw();
		leg1->Draw();
		vector< double > vcoeff_pen;
		for(unsigned j=0; j<vGpenalty.size(); j++){
			vGpenalty[j]->Draw("same pz");
			TF1 *f = new TF1(("f"+to_string(j)).c_str(),"1+[0]*x",2,30);
				f->SetParameter(0,0.27);
				f->SetLineColor(vGpenalty[j]->GetLineColor());
				f->SetLineWidth(2);
			vGpenalty[j]->Fit(f);
			vcoeff_pen.push_back(f->GetParameter(0));
		}
		
		double sum = std::accumulate(vcoeff_pen.begin(), vcoeff_pen.end(), 0.0);
		double mean = sum / vcoeff_pen.size();
		std::vector<double> diff(vcoeff_pen.size());
		std::transform(vcoeff_pen.begin(), vcoeff_pen.end(), diff.begin(), bind2nd(std::minus<double>(), mean));
		double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
		double stdev = std::sqrt(sq_sum / (vcoeff_pen.size()-1));
		cout<<"-------------------------------------------"<<endl;
		cout<<"Penalty coefficient: "<<mean<<" +/- "<<stdev<<endl;
		cout<<"-------------------------------------------"<<endl;
		
	////////////////////////////////////////////////////////////////////////////
	//                                                                        //
	//                         End of the code                                //
	//                                                                        //
	////////////////////////////////////////////////////////////////////////////

	cout << endl <<"------ Program Finished Normally: Good Job! ------" << endl;
//	gSystem->Exit(0);
	rint->Run(kTRUE);
}
