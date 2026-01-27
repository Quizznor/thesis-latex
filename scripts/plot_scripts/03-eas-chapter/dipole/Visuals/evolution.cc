//////////////////////////////////////////////////////
//                                                   //
//   Dev. by Jonathan Biteau (biteau@in2p3.fr)       //
//                       last edit: 2022-04-14       //
///////////////////////////////////////////////////////

// C/C++ classics
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

// ROOT
#include "TRint.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TSystem.h"
#include "Math/QuantFuncMathCore.h"

// UTILITIES
#include "../Catalog_Based/utils.h"


#ifdef gcc323
char* operator+( streampos&, char* );
#endif

using namespace std;

//Global variables defining the data, catalogs and composition scenarios

//Data///////////////////
string InputData_seventfile = "../Data/AugerApJS2022_Yr_JD_UTC_Th_Ph_RA_Dec_E_Expo";
string InputData_sexpofitsfile = "../Data/exposure.fits";	

//Catalogs///////////////
vector< string > InputData_vssources = {"Starburst galaxies (radio)", "Jetted AGN (#gamma-rays)", "All AGN (hard X-rays)", "Galaxies > 1 Mpc (IR)"};//Associated display name
vector< string > InputData_vsfile = {"sbg", "fermi_lat_agn", "swift_bat_agn", "2mass_hyperleda"} ;//Associated folder, either: sbg, fermi_lat_agn, swift_bat_agn, 2mass_hyperleda

//Composition////////////
vector< string > InputData_vscenarios = {"EPO1st","EPO2nd","Sib1st",""};
vector< string > InputData_vname_scenarios = {"EPOS-LHC 1st minimum","EPOS-LHC 2nd minimum","Sibyll 1st minimum","no attenuation"};

//Table 3 in arXiv:1612.07155
double Rcut_min = TMath::Power(10, 18.69-18), Rcut_best = TMath::Power(10, 18.72-18), Rcut_max = TMath::Power(10, 18.77-18);

//Exposure vs time
TGraph *Gyr2exp = new TGraph("../Data/time_exposure.dat");
double expo_year(double *x, double *par){ return Gyr2exp->Eval(x[0])/1E3;}

//Cosmetics
Color_t col[5] = {kAzure+7, kMagenta+1, kOrange+8, kSpring-8, kOrange-1};
int lstyle[5] = {1, 7, 2, 9, 9};
int lwidth[5] = {2, 3, 3, 2, 3};
int mstyle[5] = {21, 25, 20, 27, 24};

//Usage
void Usage(string sinput)
{
	cout << endl;
	cout << " Synopsis : " << endl;
	cout << sinput << " -h, --help to obtain this message" << endl;
	cout << sinput << " " << endl << endl;

	cout << " Description :" << endl;
	cout << " Return a comparison of TS profile for a given catalog and in the Cen region"<<endl
			 << " -- that's all folks! -- "
			 << endl<< endl;

	exit(0);
}

//Load Centaurus region results
TGraph* LoadCenProfile(){
	TGraph* Gcen = new TGraph("../Targeted_Blind/Targeted/fig6tmp");		
	return Gcen;
}

//Load the results vs Eth obtained with a given composition model and catalog, defined in InputData
void LoadScan_interpretation(int cat=0){

	////////////////////////////////////////////////////////////////////////////
	//                                                                        //
	//                    Load the results                                    //
	//                                                                        //
	////////////////////////////////////////////////////////////////////////////
	
	//Load catalog results
	stringstream catname;
	catname<<"#color["<<col[cat]<<"]{"<<InputData_vssources[cat]<<"}";

	string filename = "../Catalog_Based/fig3.root";
	TFile *f = new TFile(filename.c_str());
	TCanvas *c_res = (TCanvas *)f->Get("c_res");
	stringstream TSname;
	TSname<<"G_TS"<<cat;
	TGraph* G_TS = (TGraph*)c_res->FindObject(TSname.str().c_str());
		G_TS->SetLineStyle(lstyle[cat]);
		G_TS->SetMarkerStyle(mstyle[cat]);
		G_TS->SetLineColor(col[cat]);
		G_TS->SetLineWidth(lwidth[cat]);
		G_TS->SetMarkerColor(col[cat]);
		G_TS->SetMarkerSize(0.85);	

	//Load Cen results
	TGraph* Gcen = LoadCenProfile();
	for (int i=0;i<Gcen->GetN();i++) Gcen->GetY()[i]=ROOT::Math::chisquared_quantile_c(Gcen->GetY()[i],2); //convert into format to superimpose with TS
		Gcen->SetLineStyle(lstyle[4]);
		Gcen->SetMarkerStyle(mstyle[4]);
		Gcen->SetLineColor(col[4]);
		Gcen->SetLineWidth(lwidth[4]);
		Gcen->SetMarkerColor(col[4]);
		Gcen->SetMarkerSize(1.3);
	
	////////////////////////////////////////////////////////////////////////////
	//                                                                        //
	//                   Plot the results                                     //
	//                                                                        //
	////////////////////////////////////////////////////////////////////////////

	double x_offset = 1.05, y_offset = 1.0;
	double top_margin = 0.12, bottom_margin = 0.11, side_margin = 0.11;
	double text_size = 0.05;

	//Range of the TS double
	double TS_min = 0, TS_max = 31.5;
	double Emin = 31, Emax = 81;
	
	TH1D *hGraph_TS = new TH1D("hGraph_TS","",100,Emin, Emax);
		hGraph_TS->SetStats(0);
		hGraph_TS->SetDirectory(0);
		hGraph_TS->SetMinimum(TS_min);
		hGraph_TS->SetMaximum(TS_max);
		hGraph_TS->GetYaxis()->SetTitle("Test statistic, TS");
		hGraph_TS->GetXaxis()->SetTitle("Threshold energy, #it{E}_{th}  [EeV]");
		hGraph_TS->GetYaxis()->SetTitleSize(text_size);
		hGraph_TS->GetYaxis()->SetLabelSize(text_size);
		hGraph_TS->GetXaxis()->SetTitleSize(text_size);
		hGraph_TS->GetXaxis()->SetLabelSize(text_size);
		hGraph_TS->GetXaxis()->SetTitleOffset(x_offset);
		hGraph_TS->GetYaxis()->SetTitleOffset(y_offset);
		hGraph_TS->GetXaxis()->CenterTitle();
		hGraph_TS->GetYaxis()->CenterTitle();
	
	TGaxis *Glocal_pval = new TGaxis(Emax,TS_max,Emax+0.001,0,TMath::Prob(TS_max,2.),1.,510,"G-");
		Glocal_pval->SetTitle("Pre-trial #it{p}-value, #it{P}_{#chi^{2}}(TS,2)");
		Glocal_pval->CenterTitle();
		Glocal_pval->SetTitleSize(text_size);
		Glocal_pval->SetLabelSize(text_size);
		Glocal_pval->SetTitleOffset(1.4);
		Glocal_pval->SetLabelOffset(0.06);

	//Upper axis
	TGaxis *GZ = new TGaxis(Emin,TS_max,Emax,TS_max,Emin/Rcut_best,Emax/Rcut_best,510,"-");
		GZ->SetTitle("Minimum bulk charge, #it{Z}_{min} = #it{E}_{th}/#it{R}_{cut, best}(#Phi, X_{max})");
		GZ->CenterTitle();
		GZ->SetTitleSize(text_size);
		GZ->SetLabelSize(text_size);
		GZ->SetTitleSize(text_size);
		GZ->SetLabelSize(text_size);			
		GZ->SetTitleOffset(x_offset);
		GZ->SetLabelOffset(0.00);	
	
	//Load points corresponding to different Z	
	vector< TGraph *> vGZ;
	vector< TGraphAsymmErrors *> vGZ_arrows;

	for(unsigned int z=6; z<16; z++){
		double Eth = z*Rcut_best;
		double y = TS_max-1-1.*( (z-6)%2 )- 0.2*( (z-6)%4 );
		
		vector< double > vZ = {Eth, Eth};
		vector< double > vy = {TS_min, TS_max};
		TGraph* G = new TGraph(vZ.size(), &vZ[0], &vy[0]);
			G->SetLineColor(kGray+2-(z%2));
			G->SetLineWidth(1);			
			G->SetLineStyle(2);					
		vGZ.push_back(G);
		
		vector< double > vZa = {Eth};
		vector< double > vya = {y};	
		vector< double > vem_Za = {z*(Rcut_best-Rcut_min)};		
		vector< double > vep_Za = {z*(Rcut_max-Rcut_best)};						
		
		TGraphAsymmErrors* Ga = new TGraphAsymmErrors(vZa.size(), &vZa[0], &vya[0], &vem_Za[0], &vep_Za[0], 0, 0);	
			Ga->SetLineColor(kGray+2-(z%2));
			Ga->SetFillColor(kGray+2-(z%2));		
			Ga->SetMarkerColor(kGray+2-(z%2));			
			Ga->SetLineWidth(1.);	
			Ga->SetMarkerSize(1.);
			Ga->SetMarkerStyle(20);
		vGZ_arrows.push_back(Ga);	
	}
		
	//Legend
	TLegend *leg_TS = new TLegend(0.13,0.13,0.47,0.13+0.042*3);
		leg_TS->SetLineColor(kWhite); 
		leg_TS->SetFillColor(kWhite);
		leg_TS->SetMargin(0.2); 
		leg_TS->SetTextSize(hGraph_TS->GetXaxis()->GetTitleSize()); 
		stringstream ss0;
		ss0<<"#bf{#color["<<col[4]<<"]{Centaurus region}}";
		leg_TS->AddEntry(Gcen,ss0.str().c_str(),"lp");
		leg_TS->AddEntry(G_TS,("#bf{"+catname.str()+"}").c_str(),"lp");
		
	TLatex tl;
		tl.SetNDC();
		tl.SetTextAlign(12);
		tl.SetTextSize(hGraph_TS->GetXaxis()->GetTitleSize());
		tl.SetTextFont(hGraph_TS->GetXaxis()->GetTitleFont());	
		
	TLatex tl_odd;
		tl_odd.SetTextAlign(12);
		tl_odd.SetTextSize(hGraph_TS->GetXaxis()->GetTitleSize());
		tl_odd.SetTextFont(hGraph_TS->GetXaxis()->GetTitleFont());		
		tl_odd.SetTextColor(kGray+1);				
	TLatex tl_even;
		tl_even.SetTextAlign(12);
		tl_even.SetTextSize(hGraph_TS->GetXaxis()->GetTitleSize());
		tl_even.SetTextFont(hGraph_TS->GetXaxis()->GetTitleFont());		
		tl_even.SetTextColor(kGray+2);						

	//Graph of TS
	TCanvas *c_TS_E = new TCanvas("cTS","cTS");
		c_TS_E->SetLeftMargin(side_margin);
		c_TS_E->SetRightMargin(side_margin);
		c_TS_E->SetTopMargin(top_margin);
		c_TS_E->SetBottomMargin(bottom_margin);		
		
		

		
		hGraph_TS->Draw();
		Glocal_pval->Draw();
		GZ->Draw();
		for(unsigned int i=0; i<vGZ.size(); i++) vGZ[i]->Draw("same l");
		for(unsigned int i=0; i<vGZ_arrows.size(); i++) vGZ_arrows[i]->Draw("same lpz");		
		leg_TS->Draw();

		G_TS->Draw("same lp");
		Gcen->Draw("same lp");
		c_TS_E->SaveAs("fig6.pdf");
}

//Load the growth of the signal TS
TGraph* LoadGrowthTS(double theta_bf, double alpha_bf, vector< double > vexp){

	////////////////////////////////////////////////////////////////////////////
	//                                                                        //
	//     	            Event-by-event TS sorted by exposure                  //
	//                                                                        //
	////////////////////////////////////////////////////////////////////////////

	vector<double> vTS = VectorTSevents(theta_bf,alpha_bf);

	//Sort the events by increasing exposure
	vector< pair< double, double > > vdata;
	for(unsigned int i=0; i<vTS.size(); i++) vdata.push_back(make_pair(vexp[i], vTS[i]));
	sort(vdata.begin(), vdata.end());//, compare);

	vector<double> vexp_sorted(vTS.size());
	for(unsigned int i=0; i<vTS.size(); i++) vexp_sorted[i] = vdata[i].first/1000;

	//Accumulate TS on data
	vector< double > vTS_acc(vTS.size());
	for(unsigned int i=0; i<vTS.size(); i++){		
		vTS_acc [i] = vdata[i].second;
		if(i>0) vTS_acc[i]+=vTS_acc[i-1];

	}

	//Graph of the observed distribution
	TGraph *G = new TGraph(vexp_sorted.size(), &vexp_sorted[0], &vTS_acc[0]);

	return G;
}

double GetAngDist(double alphaEv, double deltaEv, double alphaObj, double deltaObj)
{

  double d2r = TMath::Pi() / 180.;
  double r2d = 180 / TMath::Pi();

  double deltaEvR = deltaEv*d2r;
  double deltaObjR = deltaObj*d2r;
  double diffAlpha = alphaEv - alphaObj;
  double da = sin (deltaEvR) * sin (deltaObjR );
  double dg = cos (deltaEvR) * cos (deltaObjR) * cos( (diffAlpha) * d2r );
  double argAcos = da + dg;
  return acos (argAcos) * r2d;
}

TGraph* LoadGrowthXS(double Eth, double angle, double alpha){
	double RACA = 201.4;
	double decCA = -43.0;
	
	// 9 columns in the event file
  	double theta, phi, UTC, yr, jd, en, RA, dec, expo;

	int n_in = 0;
	int n_out = 0;

	vector<double> xs;
	vector<double> xp;

	ifstream in;	
	in.open(InputData_seventfile .c_str());	
	in.ignore(10000, '\n')	;

	while (!in.eof()) {
		in >> yr >> jd >> UTC >> theta >> phi >> RA >> dec >> en >> expo;
		if(!in.good()) break;
		if(en>Eth){
		xp.push_back(expo/1000.);
		if(GetAngDist(RA, dec, RACA, decCA) < angle) n_in +=1;	
		else n_out +=1;
		xs.push_back((double)n_in-((double)(n_in+n_out))*alpha);
		}
	}
	TGraph* G = new TGraph(xp.size(), &xp[0], &xs[0]);		
	return G;
}

//Return 2.3, 15.9, 50, 84.1 and 97.7 percentiles
vector< double > Quantiles(vector<double> vx){
	vector< double > vres;

	vector< double > vq = {0.023, 0.159, 0.50, 0.841, 0.977};//median, 1sig, 2sig
	sort(vx.begin(), vx.end());
	for(unsigned int i=0; i<vq.size(); i++){
		int ix = round(vq[i]*vx.size());
		if(ix>=vx.size()) ix = vx.size()-1;
		vres.push_back(vx[ix]);
	}	
	
	return vres;
}


//Return contour to be filled
TGraph* LoadContour(vector< double > vx0, vector< double > vlower, vector< double > vupper){

	vector< double >  vx, vy;
	for(unsigned int i=0; i<vx0.size(); i++){
		vx.push_back(vx0[i]);
		vy.push_back(vlower[i]);
	}
	for(int i=vx0.size()-1; i>-1; i--){
		vx.push_back(vx0[i]);
		vy.push_back(vupper[i]);
	}
	vx.push_back(vx.front());
	vy.push_back(vy.front());
	TGraph *G = new TGraph(vx.size(), &vx[0], &vy[0]); 
	
	return G;
}

//Return 1 and 2 sigma confidence region of simulated XS
vector< TGraph* > LoadConfRegionsCen(vector< double > vexp_sorted, double Ncen, double Nout, double alpha, double exp_tot, int nsimu, int nseed=123){

	//Initialize random number
	TRandom3 r0(nseed);
		
	//Generate sets of Poisson fluctuations in region
	unsigned int nevts_rec = vexp_sorted.size();
	vector< vector < double > > vvXSsim(nevts_rec, std::vector<double>(nsimu));	
	for(int ievts=0; ievts<nevts_rec; ievts++){
		double Ncen_sim = Ncen*vexp_sorted[ievts]/exp_tot;
		double Nout_sim = Nout*vexp_sorted[ievts]/exp_tot;	
		for(int isim=0; isim<nsimu; isim++){
			int ncen_sim = r0.Poisson(Ncen_sim);
			int nout_sim = r0.Poisson(Nout_sim);
			vvXSsim[ievts][isim] = (1-alpha)*ncen_sim - alpha*nout_sim;
		}
	}
	
	//Reconstruct appropriate quantiles
	vector< double > vlow_2sig, vlow_1sig, vup_1sig, vup_2sig;
	for(unsigned int ievts=0; ievts<nevts_rec; ievts++){
		vector< double > vXSq = Quantiles(vvXSsim[ievts]);
		vlow_2sig.push_back(vXSq[0]);
		vlow_1sig.push_back(vXSq[1]);
		//median here
		vup_1sig.push_back(vXSq[3]);
		vup_2sig.push_back(vXSq[4]);
	}

	//Export the graphs
	vector< TGraph* > vG;
	vG.push_back(LoadContour(vexp_sorted, vlow_1sig, vup_1sig));
	vG.push_back(LoadContour(vexp_sorted, vlow_2sig, vup_2sig));		
	
	return vG;
}



//Return 1 and 2 sigma confidence region of simulated TS
vector< TGraph* > LoadConfRegionsTS(double theta_bf, double alpha_bf, THealpixMap modelMap_norm, THealpixMap exposureMap, vector< THealpixMap > vsmoothedMap, vector< double > vexp_sorted, int nsimu, int nseed=123){

	//Initialize random numbers
	vector< TRandom3 > vr0;
	for(unsigned int i=0; i<3; i++){
		TRandom3 r0(nseed+i);
		vr0.push_back(r0);
	}
	
	//Generate sets of events
	unsigned int nevts_rec = vexp_sorted.size();
	vector< vector < double > > vvTSsim(nevts_rec, std::vector<double>(nsimu));
	for(int isim=0; isim<nsimu; isim++){
		//Random variables for generation on the sphere
		int nmaxsimu=1E6;
		TRandom3 ru(int(vr0[0].Uniform(nmaxsimu))), rv(int(vr0[1].Uniform(nmaxsimu)));//for long, lat, picking
		TRandom3 rp(int(vr0[2].Uniform(nmaxsimu)));//for rejection/acceptance

		//Load the simulated evts: same number of events above decreasing Eth
		vector< double > vlsim, vbsim;	
		while(vlsim.size()<nevts_rec){
			//Uniform picking on the sphere
			double phi = 2*M_PI*ru.Uniform();// in [0;2pi]
			double theta = TMath::ACos(2.*rv.Uniform()-1.);// in [0;pi]
			double l = (phi-M_PI)*RTOD;// in [-180;180]
			double b = (theta-0.5*M_PI)*RTOD;// in [-90;90]

			//Acceptance/rejection criterion
			long pix = modelMap_norm.Ip(l,b);
			if(rp.Uniform()<=modelMap_norm[pix]){
				vlsim.push_back(l);
				vbsim.push_back(b);
			}
		}

		//Reconstruction
		double TSsim = TS_fixed_E_Par(vlsim, vbsim, exposureMap, vsmoothedMap, theta_bf, alpha_bf);
		vector< double > vTSsim = VectorTSevents(theta_bf, alpha_bf);
		vvTSsim[0][isim] = vTSsim[0];
		for(unsigned int ievts=1; ievts<nevts_rec; ievts++) vvTSsim[ievts][isim] = vvTSsim[ievts-1][isim] + vTSsim[ievts];
	}
	
	//Reconstruct appropriate quantiles
	vector< double > vlow_2sig, vlow_1sig, vup_1sig, vup_2sig;
	for(unsigned int ievts=0; ievts<nevts_rec; ievts++){
		vector< double > vTSq = Quantiles(vvTSsim[ievts]);
		vlow_2sig.push_back(vTSq[0]);
		vlow_1sig.push_back(vTSq[1]);
		//median here
		vup_1sig.push_back(vTSq[3]);
		vup_2sig.push_back(vTSq[4]);
	}

	//Export the graphs
	vector< TGraph* > vG;
	vG.push_back(LoadContour(vexp_sorted, vlow_1sig, vup_1sig));
	vG.push_back(LoadContour(vexp_sorted, vlow_2sig, vup_2sig));		
	
	return vG;
}

//Load the growth of TS and XS as a function of accumulated exposure
void LoadScan_evolution(int nsimu){

	////////////////////////////////////////////////////////////////////////////
	//                                                                        //
	//     	                 FIT THE DATA ABOVE Eth_TS                         //
	//                                                                        //
	////////////////////////////////////////////////////////////////////////////

	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
	//Hard-coded entries - only valid for Pierre Auger Collab. 2022//
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//	
	//Catalog and associated best-fit threshold
	int cat=0;
	double Eth_TS=38;

	//Load the data to determine nevts
	vector< double > vId, vZen, vl, vb, vE, vexp;
	LoadAugerData(InputData_seventfile , vId, vZen, vl, vb, vE, vexp);
	Trim(Eth_TS,vE, vZen, vl,vb, vexp);//Trim the data below threshold
	
	//Exposure
	THealpixMap exposureMap((char*)InputData_sexpofitsfile.c_str());
	exposureMap.SetCoordSys('G');
	unsigned int nSide = exposureMap.NSide();	

	//Composition model
	int compo_model=0;
	string ssrc_file = InputData_vsfile[cat];
	stringstream ssfile_compo;
	ssfile_compo<<"../Catalog_Based/Catalogs/ModelsUHECR/"<<ssrc_file<<"/"<<ssrc_file;
	if(InputData_vscenarios[compo_model].length()>1) ssfile_compo<<"_"<<InputData_vscenarios[compo_model];
	ssfile_compo<<"_threshold"<<int(Eth_TS);	


	//Load the model
	vector< string > vSrcName;
	vector< double > vl_model, vb_model, vw_model;
	LoadModel(ssfile_compo.str(),vSrcName, vl_model, vb_model, vw_model);
	vector< THealpixMap > vsmoothedMap = LoadSmoothedModelMaps(nSide, vl_model, vb_model, vw_model);

	//Fit the data
	double theta_bf, alpha_bf, logL_0;
	double TS_simplex = TS_fixed_E(vl,vb, exposureMap, vsmoothedMap, theta_bf, alpha_bf, logL_0);
	
	//Load the generation model map
	THealpixMap srcMap = exposureMap*LoadSrcMap(nSide, theta_bf, vl_model, vb_model, vw_model);
	srcMap /=srcMap.Total();
	THealpixMap modelMap = exposureMap/exposureMap.Total();//isotropic fraction
	modelMap = alpha_bf*srcMap + (1.-alpha_bf)*modelMap;
	modelMap /= modelMap.Max();//normalize max to 1 for acceptance/rejection		

	////////////////////////////////////////////////////////////////////////////
	//                                                                        //
	//                    Load the results                                    //
	//                                                                        //
	////////////////////////////////////////////////////////////////////////////

	//Load TS
	TGraph * G_TS = LoadGrowthTS(theta_bf, alpha_bf, vexp);
		G_TS->SetLineColor(col[cat]);
		G_TS->SetLineWidth(2);
		G_TS->SetMarkerSize(0);
		
	TGraph* G_TS_side = new TGraph(G_TS->GetN(), G_TS->GetX(), G_TS->GetY());
		G_TS_side->SetLineColor(col[cat]);		
		G_TS_side->SetLineWidth(2);			

	//Simulated TS growth
	vector< int > vCL = {68, 95};
	vector< int > vcolTS = {kAzure+6, kAzure+2};
	sort(vexp.begin(), vexp.end());
	vector< double > vexp_sorted;
	for(unsigned int i=0; i<vexp.size(); i++) vexp_sorted.push_back(vexp[i]/1000);//for display purpose
	vector< TGraph* > vG_TS_1_2sigma = LoadConfRegionsTS(theta_bf, alpha_bf, modelMap, exposureMap, vsmoothedMap, vexp_sorted, nsimu);
	for(unsigned int i=0; i<vG_TS_1_2sigma.size(); i++){
		vG_TS_1_2sigma[i]->SetFillColorAlpha(vcolTS[i],1);
		vG_TS_1_2sigma[i]->SetLineWidth(0);
	}

	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
	//Hard-coded entries - only valid for Pierre Auger Collab. 2022//
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//	
	//Best-fit results in Centaurus region
	double Eth_Cen=38.;			//best threshold
	double ang_Cen=27.;			//best angle
	double Ncen = 215.;			//number of observed events in Cen region of Psi = 27deg
	double Nout = 1621.-215.;			//number of events outside of that region 
	double rel_expo = 9.37/100.;	//relative exposure in Cen region of Psi = 27deg
	double exp_tot = 122.;		//total exposure in 1E3 km^2 sr yr
	double TS_ref = 25., XS_ref = 67.;//2 additionnal events in Cen to reach equiv of TS = 25
	double XS_2_TS = TS_ref/XS_ref;

	//Load Cen results
	TGraph* Gcen_side = LoadGrowthXS(Eth_Cen, ang_Cen, rel_expo);
		Gcen_side->SetLineColor(col[4]);
		Gcen_side->SetLineWidth(2);	
		
	vector< double > cenTS;
	for(int i=0; i<Gcen_side->GetN(); i++) cenTS.push_back(Gcen_side->GetY()[i]*XS_2_TS);
	TGraph* Gcen = new TGraph(Gcen_side->GetN(), Gcen_side->GetX(), &cenTS[0]);
		Gcen->SetLineColor(col[4]);
		Gcen->SetLineWidth(3);
		Gcen->SetMarkerSize(0);

	

	//Simulated Cen growth
	vector< TGraph* > vGcen_1_2sigma = LoadConfRegionsCen(vexp_sorted, Ncen, Nout, rel_expo, exp_tot, nsimu);
	vector< int > vcolCen = {kOrange-2,kOrange+3};
	for(unsigned int i=0; i<vGcen_1_2sigma.size(); i++){
		vGcen_1_2sigma[i]->SetFillColorAlpha(vcolCen[i],1);
		vGcen_1_2sigma[i]->SetLineWidth(0);
	}


	////////////////////////////////////////////////////////////////////////////
	//                                                                        //
	//                   Plot the results                                     //
	//                                                                        //
	////////////////////////////////////////////////////////////////////////////
	
	double expo_max = 121.4, TS_min = 0., TS_max = 32;
	double x_offset = 1.05, y_offset = 0.95;
	double top_margin = 0.12, bottom_margin = 0.12, side_margin = 0.10;
	double text_size = 0.05;
	
	//Right-hand side axis
	TGaxis *G = new TGaxis(expo_max,TS_min,expo_max,TS_max,TS_min/XS_2_TS,TS_max/XS_2_TS,510,"+");
		G->SetTitle("Centaurus region excess");
		G->CenterTitle();
		G->SetTitleSize(text_size);
		G->SetLabelSize(text_size);
		G->SetTitleOffset(y_offset);
		G->SetLabelOffset(0.04);
		G->SetLineColor(col[4]);
		G->SetTextColor(col[4]);
		G->SetLabelColor(col[4]); 
		
	//Upper year axis
	TF1 *fup = new TF1("fup",expo_year,2004.001,2021,0);
	TGaxis *Gup = new TGaxis(0,TS_max,expo_max,TS_max,"fup",210,"-");
		Gup->SetTitleSize(text_size);
		Gup->SetLabelSize(0.0);
		Gup->SetTitle("Year");
		Gup->CenterTitle();
		Gup->SetTitleOffset(1.2);
		Gup->SetLabelOffset(0.0);

	TGaxis *Gup_minor = new TGaxis(0,TS_max,expo_max,TS_max,"fup",220,"L-");
		Gup_minor->SetLabelSize(0.0);
		Gup_minor->SetLineColor(kGray);	
	
	//Plot anciliaries
	TH1D *hGraph1 = new TH1D("hGraph1","",100,0,expo_max);
		hGraph1->GetXaxis()->SetTitle("Pierre Auger Obs. exposure #geq 32 EeV  [10^{3} km^{2} yr sr]");
		hGraph1->SetStats(0);
		hGraph1->SetDirectory(0);
		hGraph1->SetMinimum(TS_min);
		hGraph1->SetMaximum(TS_max);
		hGraph1->GetYaxis()->SetTitle("Cumulated TS #geq #it{E}_{th}");
		hGraph1->GetYaxis()->SetTitleSize(text_size);
		hGraph1->GetYaxis()->SetLabelSize(text_size);
		hGraph1->GetXaxis()->SetTitleSize(text_size);
		hGraph1->GetXaxis()->SetLabelSize(text_size);
		hGraph1->GetXaxis()->SetTitleOffset(x_offset);
		hGraph1->GetYaxis()->SetTitleOffset(y_offset);
		hGraph1->GetXaxis()->CenterTitle();
		hGraph1->GetYaxis()->CenterTitle();
		
	TLatex l;
		l.SetTextSize(text_size);
		l.SetTextAlign(11);			

	TLegend *leg0 = new TLegend(0.13,0.85-0.06*2,0.33,0.85);
		leg0->SetLineColor(kWhite); 
		leg0->SetFillColor(kWhite);
		leg0->SetMargin(0.3); 
		leg0->SetTextSize(hGraph1->GetXaxis()->GetTitleSize());

		stringstream ssTS;
		ssTS<<"#bf{#color["<<col[cat]<<"]{"<<InputData_vssources[0]<<" - #it{E}_{th} = "<<Eth_TS<<" EeV}}";
		leg0->AddEntry(G_TS,ssTS.str().c_str(),"l");

		stringstream ssCen;
		ssCen<<"#bf{#color["<<col[4]<<"]{Centaurus region - #it{E}_{th} = "<<Eth_Cen<<" EeV}}";
		leg0->AddEntry(Gcen,ssCen.str().c_str(),"l");		

	//Main canvas
	TCanvas *c1 = new TCanvas("c0","c0",1400,600);
	double x0=0, x1=0.67, x2=1;
	double y0=0, y1=0.5, y2=1;
	double fact_side_text = 1.65;
	
	//Plot TS + Cen growth
	c1->cd();
	TPad *cgraph1 = new TPad("cgraph1","cgraph1",x0, y0,x1, y2);
	cgraph1->SetMargin(side_margin,side_margin,bottom_margin,top_margin);
	cgraph1->Draw();
	cgraph1->cd();			
		hGraph1->Draw();
		G->Draw();
		Gup_minor->Draw();
		Gup->Draw();
		leg0->Draw();
		G_TS->Draw("same l");
		Gcen->Draw("same l");
		vector< int > vyears = {2006,2008,2010,2012,2014,2016,2018,2020};
		for(unsigned int i=0; i<vyears.size(); i++) l.DrawLatex(fup->Eval(vyears[i]), TS_max*1.015,to_string(vyears[i]).c_str());
		
	//Plot TS growth
	c1->cd();
		
	TH1D *hGraphTS = new TH1D("hGraphTS","",100,0,expo_max);
		hGraphTS->GetXaxis()->SetTitle("Exposure #geq 32 EeV  [10^{3} km^{2} yr sr]");
		hGraphTS->SetStats(0);
		hGraphTS->SetDirectory(0);
		hGraphTS->SetMinimum(TS_min);
		hGraphTS->SetMaximum(TS_max);
		hGraphTS->GetYaxis()->SetTitle("Cumulated TS #geq #it{E}_{th}");
		hGraphTS->GetYaxis()->SetTitleSize(text_size*fact_side_text);
		hGraphTS->GetYaxis()->SetLabelSize(text_size*fact_side_text);
		hGraphTS->GetXaxis()->SetTitleSize(text_size*fact_side_text);
		hGraphTS->GetXaxis()->SetLabelSize(text_size*fact_side_text);
		hGraphTS->GetXaxis()->SetTitleOffset(x_offset);
		hGraphTS->GetYaxis()->SetTitleOffset(y_offset/TMath::Sqrt(fact_side_text));
		hGraphTS->GetXaxis()->CenterTitle();
		hGraphTS->GetYaxis()->CenterTitle();	

	TLegend *legTS = new TLegend(0.19,0.85-0.06*2*fact_side_text,0.45,0.85);
		legTS->SetLineColor(kWhite); 
		legTS->SetFillColor(kWhite);
		legTS->SetMargin(0.3); 
		legTS->SetTextSize(hGraphTS->GetXaxis()->GetTitleSize());

		for(unsigned int i=0; i<vG_TS_1_2sigma.size(); i++){
			stringstream ssTS;
			ssTS<<"#bf{#color["<<vcolTS[i]<<"]{"<<vCL[i]<<"% C.L.}}";		
			legTS->AddEntry(vG_TS_1_2sigma[i],ssTS.str().c_str(),"f");
		}
		
	TPad *cgraph2 = new TPad("cgraph2","cgraph2",x1, y1,x2, y2);
	cgraph2->SetMargin(side_margin*fact_side_text,0.05,bottom_margin*fact_side_text,top_margin);
	cgraph2->Draw();
	cgraph2->cd();
	cgraph2->SetTickx();
	cgraph2->SetTicky();	
		hGraphTS->Draw();
		legTS->Draw();
		for(int i=vG_TS_1_2sigma.size()-1; i>=0; i--) vG_TS_1_2sigma[i]->Draw("same f");	
		G_TS_side->Draw("same l");		
		
	//Plot Cen growth
	c1->cd();
			
	TH1D *hGraphCen = new TH1D("hGraphCen","",100,0,expo_max);
		hGraphCen->GetXaxis()->SetTitle("Exposure #geq 32 EeV  [10^{3} km^{2} yr sr]");
		hGraphCen->SetStats(0);
		hGraphCen->SetDirectory(0);
		hGraphCen->SetMinimum(TS_min/XS_2_TS);
		hGraphCen->SetMaximum(TS_max/XS_2_TS);
		hGraphCen->GetYaxis()->SetTitle("Centaurus region excess");
		hGraphCen->GetYaxis()->SetTitleColor(col[4]);
		hGraphCen->GetYaxis()->SetLabelColor(col[4]); 
		hGraphCen->GetYaxis()->SetTitleSize(text_size*fact_side_text);
		hGraphCen->GetYaxis()->SetLabelSize(text_size*fact_side_text);
		hGraphCen->GetXaxis()->SetTitleSize(text_size*fact_side_text);
		hGraphCen->GetXaxis()->SetLabelSize(text_size*fact_side_text);
		hGraphCen->GetXaxis()->SetTitleOffset(x_offset);
		hGraphCen->GetYaxis()->SetTitleOffset(y_offset/TMath::Sqrt(fact_side_text));
		hGraphCen->GetXaxis()->CenterTitle();
		hGraphCen->GetYaxis()->CenterTitle();	

	TLegend *legCen = new TLegend(0.19,0.85-0.06*2*fact_side_text,0.45,0.85);
		legCen->SetLineColor(kWhite); 
		legCen->SetFillColor(kWhite);
		legCen->SetMargin(0.3); 
		legCen->SetTextSize(hGraphCen->GetXaxis()->GetTitleSize());

		for(unsigned int i=0; i<vGcen_1_2sigma.size(); i++){
			stringstream ssCen;
			ssCen<<"#bf{#color["<<vcolCen[i]<<"]{"<<vCL[i]<<"% C.L.}}";		
			legCen->AddEntry(vGcen_1_2sigma[i],ssCen.str().c_str(),"f");
		}
		
	TPad *cgraph3 = new TPad("cgraph3","cgraph3",x1, y0,x2, y1);
	cgraph3->SetMargin(side_margin*fact_side_text,0.05,bottom_margin*fact_side_text,top_margin);
	cgraph3->Draw();
	cgraph3->cd();
	cgraph3->SetTickx();
	cgraph3->SetTicky();	
		hGraphCen->Draw();
		legCen->Draw();
		for(int i=vGcen_1_2sigma.size()-1; i>=0; i--) vGcen_1_2sigma[i]->Draw("same f");			
	Gcen_side->Draw("same l");	
	
	c1->SaveAs("fig7.pdf");	
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

	////////////////////////////////////////////////////////////////////////////
	//                                                                        //
	//              Plot the results for a fixed composition                  //
	//                                                                        //
	////////////////////////////////////////////////////////////////////////////

	LoadScan_interpretation();
	int nsimu=1000;//Nsims to obtain the contours -> 1000 for 1 and 2sigma (<5min)
	LoadScan_evolution(nsimu);//Note: hard-coded best-fit entries within that function

	////////////////////////////////////////////////////////////////////////////
	//                                                                        //
	//                         End of the code                                //
	//                                                                        //
	////////////////////////////////////////////////////////////////////////////

	cout << endl <<"------ Program Finished Normally: Good Job! ------" << endl;
	gSystem->Exit(0);
	rint->Run(kTRUE);
}

