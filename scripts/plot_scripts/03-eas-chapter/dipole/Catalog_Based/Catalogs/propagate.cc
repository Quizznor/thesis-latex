///////////////////////////////////////////////////////
//                                                   //
//   Dev. by Jonathan Biteau (jbiteau.pro@gmail.com) //
//   with inputs from M. Unger						 //
//                       last edit: 2021-11-25       //
///////////////////////////////////////////////////////

// C/C++ classics
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <sys/stat.h>
#include <numeric>
#include <time.h>
#include <utility>

// ROOT
#include "TRint.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraph.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TSystem.h"

using namespace std;

#ifdef gcc323
char* operator+( streampos&, char* );
#endif

void Usage(string sinput)
{
	cout << endl;
	cout << " Synopsis : " << endl;
	cout << sinput << " -h, --help to obtain this message" << endl;

	cout << " Description :" << endl;
	cout << " Load the attenuated flux weights, stored in ModelsUHECR/XXX"<<endl
			 << " -- that's all folks! -- "
			 << endl<< endl;

	exit(0);
}

// Scenario structure used for composition
struct Scenario {
	Scenario(const string& name,
			const unsigned int modelIndex,
			const double fP, const double fHe, const double fN,
			const double fSi, const double fFe) :
	fName(name), fModel(modelIndex)
	{
		fFraction[0] = fP;
		fFraction[1] = fHe;
		fFraction[2] = fN;
		fFraction[3] = fSi;
		fFraction[4] = fFe;
	}

	static const unsigned int gnMass = 5;
	string fName;
	unsigned int fModel;
	double fFraction[gnMass];
};

//Load a composition weight
double LoadWeight(TFile* crpFile, double distance, const unsigned int* masses, int threshold, const Scenario& s){
	double sum = 0;
	for(unsigned int iMass = 0; iMass < Scenario::gnMass; ++iMass) {
		ostringstream gname;
		gname << "atten" << masses[iMass] << "_" << threshold - 20 << "_" << s.fModel << "_0";
		TGraph* graph = (TGraph*) crpFile->Get(gname.str().c_str());
		if (!graph) {
			cerr << gname.str() << " does not exist " << endl;
			sum = 0;
			break;
		}
		if (distance < *(graph->GetX() + graph->GetN() - 1)) {
			if (distance <= 0) sum += s.fFraction[iMass];
			else sum += s.fFraction[iMass] * graph->Eval(distance);
		}
	}
	return sum;
}

//Load a composition weight marginalized over distance uncertainty
double LoadWeightMarg(TFile* crpFile, double distance, double eln_distance, const unsigned int* masses, double threshold, const Scenario& s){

	//defines the range of integration
	int npoints = 31;//number of steps between -nsigma and +nsigma for marginalisation
	double nsigma  = 3.;
	vector< double > vdist, vw_dist;
	double ln_dist0 = TMath::Log(distance);
	for(int i=0; i<npoints; i++){
		double ln_dist = ln_dist0 + eln_distance*(-nsigma + 2.*nsigma*i/(npoints-1.));//between -nsig and nsig
		vdist.push_back(TMath::Exp(ln_dist));
		vw_dist.push_back(TMath::Gaus(ln_dist,ln_dist0,eln_distance));
	}

	//computes the marginalized weight
	double weight_marg = 0.;
	double norm_marg = 0.;
	for(unsigned int i=0; i<vdist.size(); i++){
		double att_fac = LoadWeight(crpFile, vdist[i], masses, (int)round(threshold), s);
		weight_marg+=att_fac*vw_dist[i];
		norm_marg+=vw_dist[i];
	}

	return weight_marg/norm_marg;
}

//Load the source model
void LoadModel(string filename, bool mag, vector< string >& vSrcName, vector< double >& vRA, vector< double >& vDec, vector< double >& vw, vector< double >& vdist, vector< double >& velndist, bool verbose=true)
{
	//Following numbers depend on the format of the file!!
	unsigned int nfloat_per_line = 7;
	//SourceName HyperLedaName AssocName Type l(deg) b(deg) DMod eDMod D(Mpc) eD Phi[10-1000GeV](1E-10cm-2.s-1) ePhi InApJL2018?

	//Clear the vectors to file
	vRA.clear();
	vDec.clear();
	vw.clear();
	vdist.clear();
	velndist.clear();

	//Load data from file
	ifstream myfile;	
	myfile.open(filename.c_str());
	if(myfile.is_open()){
		//loops on the lines
		string line0;
		getline (myfile,line0);
		while(myfile.good()){
			string line;
			getline (myfile,line);
			stringstream ss(line);
			string srcName1, srcName2, srcName3, type;
			ss>>srcName1>>srcName2>>srcName3>>type;
			vector< double > vbuf;
			double buf;
			while(ss>>buf) vbuf.push_back(buf);
			if(vbuf.size()>=nfloat_per_line){
				vSrcName.push_back(srcName3);
				vRA.push_back(vbuf[0]);
				vDec.push_back(vbuf[1]);

				vdist.push_back(vbuf[4]);
				velndist.push_back(vbuf[5]);

				double flux = vbuf[6];
				if(mag) flux = pow(10,-0.4*(vbuf[6]-5));//zero point <-> unity flux arbitrary at m=5 for 2MASS

				vw.push_back(flux);
			}
			else if(verbose && line.size()>0) cout<<"Line unread in source file: "<<line<<endl;
		}
		myfile.close();
	}
	else cout<<"Could not open file: "<<filename<<endl;
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
	if(argc != 2) Usage(argv[0]);
	unsigned int threshold = abs(int(atof(argv[1])));

	// Initialize root
	int fargc = 1;
	TRint *rint = NULL;
	rint = new TRint("Test plot", &fargc, argv);
	gROOT->SetStyle("Plain");
	gStyle->SetTitleFont(30,"TITLE");

	// Load the file containing the weights
	const string crpFilename("anaCRP3.root");
	TFile* crpFile = TFile::Open(crpFilename.c_str());
	if (!crpFile) {
		cerr << " error opening " << crpFilename << endl;
		gSystem->Exit(0);
	}

	//Load the hadronic models
	//Names should be consistent with InputData.h
	static const unsigned int masses[Scenario::gnMass] = {1, 4, 14, 28, 56};
	vector<Scenario> vscenarios;
	vscenarios.push_back(Scenario("EPO1st", 0, 0.000, 0.673, 0.281, 0.046, 0));
	vscenarios.push_back(Scenario("EPO2nd", 1, 0.000, 0.000, 0.798, 0.202, 0));
	vscenarios.push_back(Scenario("Sib1st", 2, 0.702, 0.295, 0.003, 0.000, 0));

	////////////////////////////////////////////////////////////////////////////
	//                                                                        //
	//                  Compute the weighted models above threshold           //
	//                                                                        //
	////////////////////////////////////////////////////////////////////////////
	
	vector< string > InputData_vsmodelfile = {"Multiwavelength/Lunardini_SBG.dat","Multiwavelength/Fermi3FHL_AGN.dat","Multiwavelength/SwiftBAT105_AGN.dat","Multiwavelength/2MASS_HyperLEDA.dat"};
	vector< string > InputData_vssources = {"Starburst galaxies (radio)", "Jetted AGN (#gamma-rays)", "All AGN (hard X-rays)", "Galaxies > 1 Mpc (IR)"};//Associated display name
	vector< string > InputData_vsfile = {"sbg", "fermi_lat_agn", "swift_bat_agn", "2mass_hyperleda"} ;//Associated folder
	vector< bool > InputData_vmag = {false, false, false, true};//false if the catalog provides fluxes, true if magnitudes
	
	for(unsigned im=0; im<InputData_vsmodelfile.size(); im++){
		vector< double > vRA_model, vDec_model, vw_model, vdist_model, velndist_model;
		vector< string > vSrcName;
		LoadModel(InputData_vsmodelfile[im], InputData_vmag[im], vSrcName, vRA_model, vDec_model, vw_model, vdist_model, velndist_model);

		cout<<"Computing weights for "<<InputData_vssources[im]<<" above "<<threshold<<" EeV"<<endl; 
		
		//for the three scenarios
		for(unsigned int j=0; j<vscenarios.size(); j++){
			const Scenario& s = vscenarios[j];
			stringstream ssout;
			ssout<<"ModelsUHECR/"<<InputData_vsfile[im]<<"/"<<InputData_vsfile[im]<<"_"<<vscenarios[j].fName<<"_threshold"<<threshold;
			ofstream file_out(ssout.str().c_str());
			for(unsigned int k=0; k<vSrcName.size(); k++){
				double compo_weight;
				if(InputData_vmag[im]) compo_weight = LoadWeight(crpFile,vdist_model[k], masses, threshold, s);
				else compo_weight = LoadWeightMarg(crpFile, vdist_model[k], velndist_model[k], masses, threshold, s);
				double tot_weight = compo_weight*vw_model[k];
				file_out<<vSrcName[k]<<" "<<vRA_model[k]<<" "<<vDec_model[k]<<" "<<tot_weight<<"\n";
			}
			file_out.close();
		}

		//scenario without weight
		stringstream ssout;
		ssout<<"ModelsUHECR/"<<InputData_vsfile[im]<<"/"<<InputData_vsfile[im]<<"_threshold"<<threshold;
		ofstream file_out(ssout.str().c_str());
		for(unsigned int k=0; k<vSrcName.size(); k++) file_out<<vSrcName[k]<<" "<<vRA_model[k]<<" "<<vDec_model[k]<<" "<<vw_model[k]<<"\n";
		file_out.close();
	}

	////////////////////////////////////////////////////////////////////////////
	//                                                                        //
	//                         End of the code                                //
	//                                                                        //
	////////////////////////////////////////////////////////////////////////////

	cout << endl <<"------ Program Finished Normally: Good Job! ------" << endl;
	gSystem->Exit(0);
	rint->Run(kTRUE);
}
