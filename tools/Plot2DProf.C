// Plot a 2D histogram and overlay with a profile.  Honestly, not always
// a good plot to show.  Over/Underflow can garble it.  Makes the TProfile 
// from the TH2D.  Default option is to show the RMS rather than the error 
// on the mean.  Resulting TProfile* is public so that it may be cloned or 
// fit or whatever after doItAll is executed.


#include <iostream>
#include <fstream>
#include <sstream>
#include <deque>

#include "Rtypes.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"


//---   class definition here   ----------------------------------------
class Plot2DProf {
public:
	Plot2DProf() {
		ProfCanvas = new TCanvas("Prof2DCanvas","2DProfCanvas",50,50,600,450);
		ProfCanvas->ToggleEventStatus();
		Profile = nullptr;
	};
	Plot2DProf(const char* canvasname, int Xwide=600,int Yhigh=450, int Xpos=50,int Ypos=50) {
		ProfCanvas = new TCanvas(canvasname,canvasname,Xpos,Ypos,Xwide,Yhigh);
		ProfCanvas->ToggleEventStatus();
		Profile = nullptr;
	};

	~Plot2DProf() {
		delete ProfCanvas;
		delete Profile;
	}

	void help() {
		cout
			<< "Plot2DProf(canvasname, Xwide=600,Yhigh=900, Xpos=50,Ypos=50)" << endl;
		cout
			<< "doProf2D(hist, bool RMS=true)" << endl;
		cout
			<< "help()" << endl;
		cout
			<< "\nNeed to implement minYaxis, maxYaxis someday." << endl;
	}

	int doProf2D(char* histname, bool RMS=true) {
		TH2* hist = (TH2*)gDirectory->Get(histname);
		return doItAll(hist, RMS);
	}
	int doProf2D(TH2* hist, bool RMS=true) {
		return doItAll(hist, RMS);
	}

	TCanvas*  ProfCanvas;
	TProfile* Profile;



private:
	int doItAll(TH2* hist, bool RMS) {
		if (Profile!=NULL) delete Profile;
		if (RMS) {
			cout << "RMS option selected" << endl;
			Profile = hist->ProfileX("_pfx",1,-1,"s");
		} else {
			cout << "Err on mean option selected" << endl;
			Profile = hist->ProfileX("_pfx",1,-1);
		}

		gStyle->SetOptStat(0);
		ProfCanvas->Clear();
		ProfCanvas->cd(1);
		hist->Draw("COLZ");
		Profile->Draw("SAME");

		return 0;
	}
};
