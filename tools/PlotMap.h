// A class to simplify the handling of TH1, TH2 and TProfiles for analysis.
// Basically just an std::map<std::string,TH1*> with wrapping paper.
// Uses map rather than unsorted_map so when you look at the output file
// interactively, ls() gives the histo names in alphabetic order.
// Evidently, TProfile is derived from TH1 by Root v6.
//
// Leo Bellantoni, 2019 copyright FRA.
//
// Use like so:
//
//		TH1::SetDefaultSumw2(true);
//		PlotMap plots;
//		enum {X=0,Y=1};
//		Int_t nBins[2];	Double_t lo[2],hi[2];
//	
//		nBins[X] = nBins[Y] = 150;		lo[X] = lo[Y] = -500.0;		hi[X] = hi[Y] = +500.0;
//		plots.addTH2D("MCvtxAllXY",	"XY, MC primary vertex, \"data\"",
//				nBins[X],lo[X],hi[X],nBins[Y],lo[Y],hi[Y]);
//		.
//		.
//		.
//		Float_t XvertMC, YvertMC, ZvertMC;
//		int nMCvert = MCVertX.size();
//		for (int iMCvert=0; iMCvert<nMCvert; ++iMCvert) {
//			XvertMC = MCVertX.getData(iMCvert);
//			YvertMC = MCVertY.getData(iMCvert);
//			plots.Fill("MCvtxAllXY", XvertMC,YvertMC);
//		}
//		.
//		.
//		.
//		TFile Output(RootOutName.c_str(),"RECREATE");
//		if (Output.IsZombie()) {cout << "Eat you brains!" << endl; return 1;};
//		plots.writePlots();
//		Output.Close();



#include <map>
using std::map;
#include <string>
using std::string;

#include "Rtypes.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"





class PlotMap {
public:
	PlotMap() { TH1::SetDefaultSumw2(true); }



	void addTH1D(string name, string title,
							  Int_t nBins,  Double_t lo,  Double_t hi) {
		table[name] = new TH1D(name.c_str(),title.c_str(),
							  nBins,lo,hi);
	}
	void addTH1D(string name, string title,
							  Int_t nBins,  Double_t* binArray) {
		table[name] = new TH1F(name.c_str(),title.c_str(),
							  nBins,binArray);
	}

	void addTH2D(string name, string title,
							  Int_t nBinsX, Double_t loX, Double_t hiX,
							  Int_t nBinsY, Double_t loY, Double_t hiY) {
		table[name] = new TH2D(name.c_str(),title.c_str(),
							  nBinsX,loX,hiX, nBinsY,loY,hiY);
	}
	void addTH2D(string name, string title,
							  Int_t nBinsX, Double_t* binArrayX,
							  Int_t nBinsY, Double_t loY, Double_t hiY) {
		table[name] = new TH2D(name.c_str(),title.c_str(),
							  nBinsX,binArrayX, nBinsY,loY,hiY);
	}
	void addTH2D(string name, string title,
							  Int_t nBinsX, Double_t loX, Double_t hiX,
							  Int_t nBinsY, Double_t* binArrayY) {
		table[name] = new TH2D(name.c_str(),title.c_str(),
							  nBinsX,loX,hiX, nBinsY,binArrayY);

	}
	void addTH2D(string name, string title,
							  Int_t nBinsX, Double_t* binArrayX,
							  Int_t nBinsY, Double_t* binArrayY) {
		table[name] = new TH2D(name.c_str(),title.c_str(),
							  nBinsX,binArrayX, nBinsY,binArrayY);

	}

	// Default opt is err on mean; "S" is RMS.
	void addProf(string name, string title,
							  Int_t nBins,  Double_t lo,  Double_t hi,
							  string opt="") {
		table[name] = new TProfile(name.c_str(),title.c_str(),
							  nBins,lo,hi, opt.c_str());
	}
	void addProf(string name, string title,
							  Int_t nBinsX, Double_t loX, Double_t hiX,
							  Double_t loY, Double_t hiY,
							  string opt="") {
		table[name] = new TProfile(name.c_str(),title.c_str(),
							  nBinsX,loX,hiX, loY,hiY, opt.c_str());
	}
	void addProf(string name, string title,
							  Int_t nBins,  Double_t* binArray,
							  string opt="") {
		table[name] = new TProfile(name.c_str(),title.c_str(),
							  nBins,binArray, opt.c_str());
	}
	void addProf(string name, string title,
							  Int_t nBins,  Double_t* binArray,
							  Double_t loY, Double_t hiY,
							  string opt="") {
		table[name] = new TProfile(name.c_str(),title.c_str(),
							  nBins,binArray, loY,hiY, opt.c_str());
	}



	TH1* operator[](string name) {return table[name];}

	void writePlots(string outName) {
		string RootOutName = outName;
		RootOutName += ".root";
		TFile output(RootOutName.c_str(),"RECREATE");
		if (output.IsZombie()) {
			cout << "Eat you brains!" << endl;
			output.Close();
			return;
		};
		writePlots();
		output.Close();
	}

	void writePlots() {
		map<string,TH1*>::iterator it;
		for (it = table.begin(); it!=table.end(); ++it) it->second->Write();
	}

private:
map<string,TH1*> table;

};
