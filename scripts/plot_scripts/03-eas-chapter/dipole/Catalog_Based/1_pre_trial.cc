///////////////////////////////////////////////////////
//                                                   //
//   Dev. by Jonathan Biteau (biteau@in2p3.fr)       //
//                       last edit: 2022-04-14       //
///////////////////////////////////////////////////////

// C/C++ classics
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

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
#include "TLegend.h"
#include "TLatex.h"
#include "TF1.h"
#include "TFile.h"
#include "TGaxis.h"

// UTILITIES CATALOG ANALYSIS
#include "utils.h"

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

//Cosmetics//////////////
Color_t col[4] = {kAzure+6, kMagenta+1, kOrange+8, kSpring-8};
int lstyle[4] = {1, 7, 2, 9};
int lwidth[4] = {2, 3, 3, 2};
int mstyle[4] = {21, 25, 20, 27};

//Event axis/////////////
TGraph* Gth_to_nevts;
TGraph* Gnevts_to_th;
double Eth_to_nevts(double Eth){return Gth_to_nevts->Eval(Eth, 0, "S");}
double nevts_to_Eth(double nevts){return Gnevts_to_th->Eval(nevts, 0, "S");}
double nevts_to_th(double *x, double* par){return Gnevts_to_th->Eval(x[0], 0, "S");}
/////////////////////////

void Usage(string sinput)
{
	cout << endl;
	cout << " Synopsis : " << endl;
	cout << sinput << " -h, --help to obtain this message" << endl;
	cout << " Description :" << endl;
	cout << " Scan as a function of threshold energy: output -> fig3.root and fig4.root"<<endl
		<< " -- that's all folks! -- "
		<< endl<< endl;

	exit(0);
}

//Simply returns the best alpha, theta and TS above Eth
void LoadResultsFirstPass(double& alpha_bestfit, double& theta_bestfit, double& TS_bestfit, string fileout, THealpixMap exposureMap, vector<double> vl, vector<double> vb, double Eth = 32, int cat = 0, int compo_model=0, bool fast_computation = false){

	////////////////////////////////////////////////////////////////////////////
	//                                                                        //
	//                      Model loading                                     //
	//                                                                        //
	////////////////////////////////////////////////////////////////////////////
		
	//Select catalog and set cosmetics
	string ssrc_file = InputData_vsfile[cat];

	//Composition model
	stringstream ssfile_compo;
	ssfile_compo<<"Catalogs/ModelsUHECR/"<<ssrc_file<<"/"<<ssrc_file;
	if(InputData_vscenarios[compo_model].length()>1) ssfile_compo<<"_"<<InputData_vscenarios[compo_model];
	ssfile_compo<<"_threshold"<<int(Eth);	

	unsigned int nSide = exposureMap.NSide();	

	////////////////////////////////////////////////////////////////////////////
	//                                                                        //
	//     	                 FIT THE DATA ABOVE Eth                           //
	//                                                                        //
	////////////////////////////////////////////////////////////////////////////

	//Load the model
	vector< string > vSrcName;
	vector< double > vl_model, vb_model, vw_model;
	LoadModel(ssfile_compo.str(),vSrcName, vl_model, vb_model, vw_model);
	vector< THealpixMap > vsmoothedMap = LoadSmoothedModelMaps(nSide, vl_model, vb_model, vw_model);

	//Fit the data
	double theta_bf, alpha_bf, logL_0;
	double TS_simplex = TS_fixed_E(vl,vb, exposureMap, vsmoothedMap, theta_bf, alpha_bf, logL_0);

	////////////////////////////////////////////////////////////////////////////
	//                                                                        //
	//     	                 LOAD TS VS PARAMS                                //
	//                                                                        //
	////////////////////////////////////////////////////////////////////////////

	if(fast_computation){
		alpha_bestfit = alpha_bf;
		theta_bestfit = theta_bf;
		TS_bestfit = TS_simplex;			
	}
	else{
		//Load the best fit likelihood maps
		double alpha_max = 0.5, theta_min = 3, theta_max=50;//range of the parameter space to plot
		int nbins = 250.;
		double dalpha = alpha_max/(nbins-1.);
		double dtheta = (theta_max-theta_min)/(nbins-1.);

		//scans vs alpha, theta
		vector< double > valpha_theta, vtheta, vTS_theta;
		for(int i=0; i<nbins; i++){
			double theta = theta_min+i*dtheta;
			for(int j=0; j<nbins; j++){
				double alpha_tmp = j*dalpha;
				double logLij = LogLikelihood(theta, alpha_tmp);
				if(logLij>-1.E50){
					valpha_theta.push_back(alpha_tmp);
					vtheta.push_back(theta);
					vTS_theta.push_back(2*(logLij-logL_0));
				}
			}
		}

		//Plot TS vs ...
		double dbin = 0.5*dalpha;
		double dbin_theta = 0.5*dtheta;
		//Cosmetics
		double x_offset = 1.05, y_offset = 1.0;
		double top_margin = 0.1, bottom_margin = 0.11, left_margin = 0.11, right_margin = 0.15;
		double text_size = 0.05;	
		
		TH2D *h_theta = LoadHisto("h2D", "TS = 2 ln L(#alpha,#Theta) / L(0,#bullet)", "Signal fraction, #alpha", "Fisher search radius, #Theta [deg]",  nbins,-dbin, alpha_max+dbin, nbins,theta_min-dbin_theta, theta_max+dbin_theta, valpha_theta, vtheta, vTS_theta);
			h_theta->SetContour(100.);
			h_theta->SetTitleSize(text_size);
			h_theta->GetYaxis()->SetTitleSize(text_size);
			h_theta->GetYaxis()->SetLabelSize(text_size);
			h_theta->GetXaxis()->SetTitleSize(text_size);
			h_theta->GetXaxis()->SetLabelSize(text_size);
			h_theta->GetXaxis()->SetTitleOffset(x_offset);
			h_theta->GetYaxis()->SetTitleOffset(y_offset);
			h_theta->GetZaxis()->SetTitleSize(text_size);
			h_theta->GetZaxis()->SetLabelSize(text_size);	
				
		int locmax_x, locmax_y, locmax_z;
		int i_max = h_theta->GetMaximumBin(locmax_x, locmax_y, locmax_z);
		alpha_bestfit = h_theta->GetXaxis()->GetBinCenter(locmax_x);
		theta_bestfit = h_theta->GetYaxis()->GetBinCenter(locmax_y);
		TS_bestfit = h_theta->GetBinContent(i_max);	
		cout<<"Test statistic: "<<TS_bestfit<<endl;
		
		////////////////////////////////////////////////////////////////////////////
		//                                                                        //
		//     	                 PLOT TS VS PARAMS                                //
		//                                                                        //
		////////////////////////////////////////////////////////////////////////////
	
		//Contour levels for 1, 2 sigma assuming 2 dof
		int ndof = 2, ncont_level = 1;
		vector< double > vsig, vTS;
		double TS_max = 25, dTS=0.1;
		vTS.push_back(0.);
		while(vTS.back()<=TS_max) vTS.push_back(vTS.back()+dTS);
		for(unsigned int i=0; i<vTS.size(); i++) vsig.push_back(sqrt(2)*TMath::ErfcInverse(TMath::Prob(vTS[i],ndof)));
		TGraph* G_sig_TS = new TGraph(vsig.size(),&vsig[0],&vTS[0]);
		vector< double > vcont_level;
		for(int i=1; i<=ncont_level; i++) vcont_level.push_back(TS_bestfit-G_sig_TS->Eval(i));

		vector< TGraph* > vG_theta = LoadContoursAndBestFit(h_theta,vcont_level,alpha_bestfit,theta_bestfit);		
		
		TLatex tl;
			tl.SetNDC();
			tl.SetTextAlign(22);
			tl.SetTextSize(h_theta->GetXaxis()->GetTitleSize());
			tl.SetTextFont(h_theta->GetXaxis()->GetTitleFont());
			
		TLegend *leg_TS = new TLegend(0.38,0.13,0.62,0.13+0.044*3);
			leg_TS->SetName("leg");
			leg_TS->SetLineColor(kWhite); 
			leg_TS->SetFillColor(kWhite);
			leg_TS->SetMargin(0.2); 
			leg_TS->SetTextSize(h_theta->GetXaxis()->GetTitleSize()); 
			stringstream ssEth;
			ssEth<<"#it{E}#geq "<<Eth<<" EeV";
			leg_TS->AddEntry(vG_theta[vG_theta.size()-1],("Best-fit at "+ssEth.str()).c_str(),"p");	
			leg_TS->AddEntry(vG_theta[0],"68\% confidence contour","l");			

		//Plot		
		gStyle->SetPalette(53);
		TFile* f_fig4 = TFile::Open(fileout.c_str(),"UPDATE");
		stringstream cname;
		cname<<"c2D_cat"<<cat<<"_Eth"<<Eth;
		TCanvas *cTheta = new TCanvas(cname.str().c_str(),cname.str().c_str());
			cTheta->cd();
			cTheta->SetLeftMargin(left_margin);
			cTheta->SetRightMargin(right_margin);
			cTheta->SetTopMargin(top_margin);
			cTheta->SetBottomMargin(bottom_margin);	
			cTheta->SetTickx();
			cTheta->SetTicky();

			h_theta->Draw("colz");
			string sexp = InputData_vssources[cat]+" - "+ssEth.str(); 
			tl.DrawLatex(0.48,0.94,sexp.c_str());
			for(unsigned int i=0; i<vG_theta.size()-1; i++) vG_theta[i]->Draw("same l");
			vG_theta[vG_theta.size()-1]->Draw("same p");
			leg_TS->Draw();
		cTheta->Write();	
		f_fig4->Close();	
	}
}




//Store TS, alpha, theta scan vs Eth in a root file
void PlotTSvsEth(string fileout, vector< double > vEth, vector< double > vnevts, vector< vector< double > > valpha, vector< vector< double > > vtheta, vector< vector< double > > vTS, int compo_model){

	////////////////////////////////////////////////////////////////////////////
	//                                                                        //
	//            Load the results for a fixed composition                    //
	//                                                                        //
	////////////////////////////////////////////////////////////////////////////

	//Loop on the catalogs
	vector< string > vcat;
	vector< TGraph* > vG_TS;
	vector< TGraph* > vG_alpha, vG_theta;
	for(unsigned int j =0; j<vTS.size(); j++){
		stringstream catname;
		catname<<"#color["<<col[j]+1<<"]{"<<InputData_vssources[j]<<"}";
		vcat.push_back(catname.str());
	
		vector <double> vts = vTS[j], va = valpha[j], vt = vtheta[j];
		TGraph *G_TS = new TGraph(vEth.size(), &vEth[0], &vts[0]);
		TGraph *G_alpha = new TGraph(vEth.size(), &vEth[0], &va[0]);
		TGraph *G_theta = new TGraph(vEth.size(), &vEth[0], &vt[0]);				
		vector< TGraph* > vG_loaded = {G_TS, G_alpha, G_theta}; 
		for(unsigned int i=0; i<vG_loaded.size(); i++){
			vG_loaded[i]->SetLineStyle(lstyle[j]);
			vG_loaded[i]->SetMarkerStyle(mstyle[j]);
			vG_loaded[i]->SetLineColor(col[j]+1);
			vG_loaded[i]->SetLineWidth(lwidth[j]);
			vG_loaded[i]->SetMarkerColor(col[j]+1);
			vG_loaded[i]->SetMarkerSize(0.85);
		}
		stringstream ssTS;
		ssTS<<"G_TS"<<j;
		vG_TS.push_back(vG_loaded[0]);
		vG_TS.back()->SetName(ssTS.str().c_str());
		
		stringstream ssalpha;
		ssalpha<<"G_alpha"<<j;
		vG_alpha.push_back(vG_loaded[1]);
		vG_alpha.back()->SetName(ssalpha.str().c_str());
				
		stringstream sstheta;
		sstheta<<"G_theta"<<j;
		vG_theta.push_back(vG_loaded[2]);		
		vG_theta.back()->SetName(sstheta.str().c_str());		
	}
	
	////////////////////////////////////////////////////////////////////////////
	//                                                                        //
	//            Plot the TS results for a fixed composition                 //
	//                                                                        //
	////////////////////////////////////////////////////////////////////////////
	double x_offset = 1.05, y_offset = 0.75;
	double top_margin = 0.14, bottom_margin = 0.14, side_margin = 0.09;
	double text_size = 0.06;
	
	//Range of the TS plot
	double TS_min = 0, TS_max = 28;
	double Emin = 31, Emax = 81;
		
	TH1D *hGraph_TS = new TH1D("hGraph_TS","",100,Emin, Emax);
		hGraph_TS->SetStats(0);
		hGraph_TS->SetDirectory(0);
		hGraph_TS->SetMinimum(TS_min);
		hGraph_TS->SetMaximum(TS_max);
		hGraph_TS->GetYaxis()->SetTitleSize(text_size);
		hGraph_TS->GetYaxis()->SetLabelSize(text_size);
		hGraph_TS->GetXaxis()->SetTitleSize(text_size);
		hGraph_TS->GetXaxis()->SetTitleSize(text_size);		
		hGraph_TS->GetXaxis()->SetLabelSize(text_size);
		hGraph_TS->GetXaxis()->SetTitleOffset(x_offset);
		hGraph_TS->GetYaxis()->SetTitleOffset(y_offset);		
		hGraph_TS->GetYaxis()->SetTitle("Test statistic, TS");
		hGraph_TS->GetXaxis()->SetTitle("Threshold energy, #it{E}_{th}  [EeV]");
		hGraph_TS->GetXaxis()->CenterTitle();
		hGraph_TS->GetYaxis()->CenterTitle();
		TString s("0 ");
		hGraph_TS->GetYaxis()->ChangeLabel(1,-1,-1,21, -1, -1, "0 ");

	//Event axis
	Gth_to_nevts = new TGraph(vnevts.size(), &vEth[0], &vnevts[0]);	
	Gnevts_to_th = new TGraph(vnevts.size(), &vnevts[0], &vEth[0]);
	TF1 *fnevts_to_th = new TF1("fnevts_to_th",nevts_to_th, Eth_to_nevts(Emin), Eth_to_nevts(Emax), 0);

	TGaxis *Gevts_TS = new TGaxis(Emin,TS_max,Emax,TS_max,"fnevts_to_th",510,"-G");
		Gevts_TS->SetTitle("#it{N}_{events}(#geq #it{E}_{th})");
		Gevts_TS->CenterTitle();
		Gevts_TS->SetTitleSize(text_size);
		Gevts_TS->SetLabelSize(text_size);

	TLegend *leg_TS = new TLegend(0.50,0.83-0.06*vG_TS.size(),0.90,0.83);
		leg_TS->SetLineColor(kWhite); 
		leg_TS->SetFillColor(kWhite);
		leg_TS->SetMargin(0.25); 
		leg_TS->SetTextSize(0.05); 
		for(unsigned int i=0; i<vcat.size(); i++) leg_TS->AddEntry(vG_TS[i],("#bf{"+vcat[i]+"}").c_str(),"lp");
		
	TLatex tl;
		tl.SetNDC();
		tl.SetTextAlign(12);
		tl.SetTextSize(0.05);
		tl.SetTextFont(hGraph_TS->GetXaxis()->GetTitleFont());	

	TLatex tl_axis;
		tl_axis.SetTextAlign(22);
		tl_axis.SetTextSize(hGraph_TS->GetXaxis()->GetTitleSize());
		tl_axis.SetTextFont(hGraph_TS->GetXaxis()->GetTitleFont());	
		
	TLatex tl_ticks;
		tl_ticks.SetTextAlign(23);
		tl_ticks.SetTextSize(0.5*hGraph_TS->GetXaxis()->GetTitleSize());
		tl_ticks.SetTextFont(hGraph_TS->GetXaxis()->GetTitleFont());		
		
	////////////////////////////////////////////////////////////////////////////
	//                                                                        //
	//            Plot the TS results for a fixed composition                 //
	//                                                                        //
	////////////////////////////////////////////////////////////////////////////

	TCanvas *c_res = new TCanvas("c_res","c_res", 700, 1200);
	c_res->Divide(1,3,0.001,0.001);

	c_res->cd(1);
		c_res->cd(1)->SetLeftMargin(side_margin);
		c_res->cd(1)->SetRightMargin(side_margin);
		c_res->cd(1)->SetTopMargin(top_margin);
		c_res->cd(1)->SetBottomMargin(0);		
		c_res->cd(1)->SetTicky();
		hGraph_TS->Draw();
		Gevts_TS->Draw();
		leg_TS->Draw();
		vector< double > vnth = {2000, 1000, 500, 200, 100};
		for(unsigned int i=0; i<vnth.size(); i++){
			stringstream sslab;
			sslab<<int(vnth[i]);
			tl_axis.DrawLatex(nevts_to_Eth(vnth[i]),TS_max*1.04, sslab.str().c_str());
			tl_ticks.DrawLatex(nevts_to_Eth(vnth[i]),TS_max, "#bf{|}");
		}
		tl.DrawLatex(0.12,0.10,"#bf{Attenuation model:}");
		tl.DrawLatex(0.12,0.05,("#bf{"+InputData_vname_scenarios[compo_model]+"}").c_str());
		for(int i=vcat.size()-1; i>-1; i--) vG_TS[i]->Draw("same lp");
		
	////////////////////////////////////////////////////////////////////////////
	//                                                                        //
	//            Plot the alpha results for a fixed composition              //
	//                                                                        //
	////////////////////////////////////////////////////////////////////////////

	//Range of the TS plot
	double alpha_max = 21;
	
	TGaxis *Gevts_alpha = new TGaxis(Emin,alpha_max,Emax,alpha_max,"fnevts_to_th",510,"-G");
		Gevts_alpha->SetTitle("#it{N}_{events}(#geq #it{E}_{th})");
		Gevts_alpha->CenterTitle();
		Gevts_alpha->SetTitleSize(text_size);
		Gevts_alpha->SetLabelSize(text_size);
		
	TH1D *hGraph_alpha = new TH1D("hGraph_alpha","",100,Emin, Emax);
		hGraph_alpha->SetStats(0);
		hGraph_alpha->SetDirectory(0);
		hGraph_alpha->SetMinimum(0);
		hGraph_alpha->SetMaximum(alpha_max);
		hGraph_alpha->GetYaxis()->SetTitle("Signal fraction, #alpha  [%]");
		hGraph_alpha->GetXaxis()->SetTitle("");
		hGraph_alpha->GetYaxis()->SetTitleOffset(y_offset);
		hGraph_alpha->GetXaxis()->CenterTitle();
		hGraph_alpha->GetYaxis()->CenterTitle();
		hGraph_alpha->GetYaxis()->SetTitleSize(text_size);
		hGraph_alpha->GetYaxis()->SetLabelSize(text_size);
		hGraph_alpha->GetXaxis()->SetTitleSize(text_size);
		hGraph_alpha->GetXaxis()->SetLabelSize(0);
		hGraph_alpha->GetXaxis()->SetTitleOffset(x_offset);
		hGraph_alpha->GetYaxis()->SetTitleOffset(y_offset);		

	//Graph of alpha
	c_res->cd(2);
		c_res->cd(2)->SetLeftMargin(side_margin);
		c_res->cd(2)->SetRightMargin(side_margin);
		c_res->cd(2)->SetTopMargin(0.5*top_margin);
		c_res->cd(2)->SetBottomMargin(0.5*top_margin);		
		c_res->cd(2)->SetTicky();		
		c_res->cd(2)->SetTickx();			
		hGraph_alpha->Draw();
		for(int i=vcat.size()-1; i>-1; i--) vG_alpha[i]->Draw("same lpz");		

		
	////////////////////////////////////////////////////////////////////////////
	//                                                                        //
	//            Plot the theta results for a fixed composition              //
	//                                                                        //
	////////////////////////////////////////////////////////////////////////////

	//Range of the TS plot
	double theta_min = 0, theta_max = 26;
	
	TGaxis *Gevts_theta = new TGaxis(Emin,theta_max,Emax,theta_max,"fnevts_to_th",510,"-G");
		Gevts_theta->SetTitle("#it{N}_{events}(#geq #it{E}_{th})");
		Gevts_theta->CenterTitle();
		Gevts_theta->SetTitleSize(text_size);
		Gevts_theta->SetLabelSize(text_size);
		
	TH1D *hGraph_theta = new TH1D("hGraph_theta","",100,Emin, Emax);
		hGraph_theta->SetStats(0);
		hGraph_theta->SetDirectory(0);
		hGraph_theta->SetMinimum(theta_min);
		hGraph_theta->SetMaximum(theta_max);
		hGraph_theta->GetYaxis()->SetTitle("Fisher search radius, #Theta  [deg]");
		hGraph_theta->GetXaxis()->SetTitle("Threshold energy, #it{E}_{th}  [EeV]");
		hGraph_theta->GetYaxis()->SetTitleOffset(y_offset);
		hGraph_theta->GetXaxis()->CenterTitle();
		hGraph_theta->GetYaxis()->CenterTitle();
		hGraph_theta->GetYaxis()->SetTitleSize(text_size);
		hGraph_theta->GetYaxis()->SetLabelSize(text_size);
		hGraph_theta->GetXaxis()->SetTitleSize(text_size);
		hGraph_theta->GetXaxis()->SetLabelSize(text_size);
		hGraph_theta->GetXaxis()->SetTitleOffset(x_offset);
		hGraph_theta->GetYaxis()->SetTitleOffset(y_offset);		

	//Graph of theta
	c_res->cd(3);
		c_res->cd(3)->SetLeftMargin(side_margin);
		c_res->cd(3)->SetRightMargin(side_margin);
		c_res->cd(3)->SetTopMargin(0);
		c_res->cd(3)->SetBottomMargin(bottom_margin);		
		c_res->cd(3)->SetTicky();			
		c_res->cd(3)->SetTickx();			
		hGraph_theta->Draw();
		for(int i=vcat.size()-1; i>-1; i--) vG_theta[i]->Draw("same lpz");

	c_res->SaveAs("fig3.pdf");
	TFile* f_fig3 = TFile::Open(fileout.c_str(),"RECREATE");	
	c_res->Write();
	f_fig3->Close();
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
	vector< double > vEth, vnevts;
	for(unsigned int Eth=32; Eth<81; Eth++){
		//Trim the data below threshold
		Trim(Eth,vE, vZen, vl,vb, vexp);
		
		vEth.push_back((double)Eth);
		vnevts.push_back(vE.size());
	}	
	
	//Composition model
	int compo_model = 0;//by default EPOS-LHC composition model

	//Exposure/////////////////////////////////////////////////////////////////
	THealpixMap exposureMap((char*)InputData_sexpofitsfile.c_str());
	exposureMap.SetCoordSys('G');

	//Runs the fit above each threshold 
	string fileout_fig4 = "fig4.root";
	TFile* f_fig4 = TFile::Open(fileout_fig4.c_str(),"RECREATE");
	f_fig4->Close();
	vector< vector< double > > valpha, vtheta, vTS;
	for(unsigned cat = 0; cat<InputData_vssources.size(); cat++){
		cout<<"-------------- "<<InputData_vssources[cat]<<" --------------"<<endl;

		//Events///////////////////////////////////////////////////////////////////	
		vector< double > vId, vZen, vl, vb, vE, vexp;
		LoadAugerData(InputData_seventfile , vId, vZen, vl, vb, vE, vexp);
		
		vector< double > valpha_singlecat, vtheta_singlecat, vTS_singlecat;
		for(unsigned int i=0; i<vEth.size(); i++){
			Trim(vEth[i],vE, vZen, vl, vb, vexp);//Trim the data below threshold

			double alpha, theta, TS;
			cout<<"-- Threshold: "<<vEth[i]<<" EeV"<<endl;
			LoadResultsFirstPass(alpha, theta, TS, fileout_fig4, exposureMap, vl, vb, vEth[i], cat, compo_model);
			valpha_singlecat.push_back(100.*alpha);
			vtheta_singlecat.push_back(theta);		
			vTS_singlecat.push_back(TS);				
		}
		valpha.push_back(valpha_singlecat);
		vtheta.push_back(vtheta_singlecat);
		vTS.push_back(vTS_singlecat);
	}
	
	//Saves the TS/alpha/theta vs Eth in dedicated root file
	string fileout_fig3 = "fig3.root";
	PlotTSvsEth(fileout_fig3, vEth, vnevts, valpha, vtheta, vTS, compo_model);

	////////////////////////////////////////////////////////////////////////////
	//                                                                        //
	//                         End of the code                                //
	//                                                                        //
	////////////////////////////////////////////////////////////////////////////

	cout << endl <<"------ Scan Step Finished Normally: Good Job! ------" << endl;
	gSystem->Exit(0);
	rint->Run(kTRUE);
}
