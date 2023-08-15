// Data/MC ratio plotter for TH2D.  doRat for errors in quadrature, doEff for binomial
// uncertainties.



#include <iostream>
#include <fstream>
#include <sstream>
#include <deque>

#include "Rtypes.h"
#include "TMath.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"

using namespace std;



//---   class definition here   ----------------------------------------
class PlotRatio2D {
public:
	PlotRatio2D() {
		rawCanvas   = new TCanvas("rawCanvas",  "rawCanvas",50,50,600,900);
		ratioCanvas = new TCanvas("ratioCanvas","ratioCanvas",100,100,600,600);
		rawCanvas->ToggleEventStatus();
		ratioCanvas->ToggleEventStatus();
	};
	PlotRatio2D(const char* canvasname, int Xwide=600,int Yhigh=900, int Xpos=50,int Ypos=50) {
		string rawName   = string("raw_")   + string(canvasname);
		string ratioName = string("ratio_") + string(canvasname);
		rawCanvas   = new TCanvas(rawName.c_str(),  canvasname,Xpos,   Ypos,   Xwide,Yhigh);
		ratioCanvas = new TCanvas(ratioName.c_str(),canvasname,Xpos+50,Ypos+50,Xwide,Xwide);
		rawCanvas->ToggleEventStatus();
		ratioCanvas->ToggleEventStatus();
	};
	~PlotRatio2D() {
		delete rawCanvas;
		delete ratioCanvas;
	}

	void help() {
		cout
			<< "Nondefault constructor:\tcanvasname, Xwide=600,Yhigh=900, Xpos=50,Ypos=50" << endl;
		cout
			<< "doRat:\thistNum, histDen, Log=false, Scale=true, Opts=ZCOL" << endl;
		cout
			<< "doEff:\thistNum, histDen, Log=true,  Scale=false, Opts=ZCOL" << endl;
		cout
			<< "PlotRatio2D* ratioPlotter;" << endl;
		cout
			<< "ratioPlotter->help()" << endl;
		cout
			<< "Need to implement minYaxis, maxYaxis someday." << endl;
	}

	// NULL plotting option is actually ZCOL.  Non-NULL will change how ratio is plotted only.
	int doRat(char* histNum, char* histDen, bool Log=false, bool Scale=true, char* Opts=NULL) {
		TH2D* hNum = (TH2D*)gDirectory->Get(histNum);
		TH2D* hDen = (TH2D*)gDirectory->Get(histDen);
		return doFullRatio(hNum,hDen,false,Log,Scale,Opts);
	}
	int doEff(char* histNum, char* histDen, bool Log=true, bool Scale=false, char* Opts=NULL) {
		TH2D* hNum = (TH2D*)gDirectory->Get(histNum);
		TH2D* hDen = (TH2D*)gDirectory->Get(histDen);
		return doFullRatio(hNum,hDen,true,Log,Scale,Opts);
	}

	int doRat(TH2D* hNum, TH2D* hDen, bool Log=false, bool Scale=true, char* Opts=NULL) {
		return doFullRatio(hNum,hDen,false,Log,Scale,Opts);
	}
	int doEff(TH2D* hNum, TH2D* hDen, bool Log=true, bool Scale=false, char* Opts=NULL) {
		return doFullRatio(hNum,hDen,true,Log,Scale,Opts);
	}


	
private:
	TCanvas* rawCanvas;
	TCanvas* ratioCanvas;



	int doFullRatio(TH2D* hNum, TH2D* hDen,
					bool binomial, bool Log, bool Scale, char* Opts) {

		// Check the binning
		Int_t	  NbinX[2];				Int_t    NbinY[2];
		NbinX[0] = hNum->GetNbinsX();
		NbinX[1] = hDen->GetNbinsX();
		if (NbinX[0] != NbinX[1]) {
			cout << "PlotRatio: Binning mismatch (Nbin; X)" << endl;
			return 1;
		}
		NbinY[0] = hNum->GetNbinsY();
		NbinY[1] = hDen->GetNbinsY();
		if (NbinY[0] != NbinY[1]) {
			cout << "PlotRatio: Binning mismatch (Nbin; Y)" << endl;
			return 1;
		}

		int Ncells = (NbinX[0]+2)*(NbinY[0]+2);

		// I think one needs to use Sumw2() in all scenarios
		if (hNum->GetSumw2N()==0) {
			cout << "Numerator never heard you say Sumw2()\n";
			return 1;
		}
		if (hDen->GetSumw2N()==0) {
			cout << "Denominator never heard you say Sumw2()\n";
			return 1;
		}

		// Create RAT for ratio
		TH2D* RAT = (TH2D*)hNum->Clone();
		RAT->SetName("RAT");
		RAT->Divide(hDen);

		// Create SCL = scaled denominator for plotting, Chi2 only
		TH2D* SCL = (TH2D*)hDen->Clone();		SCL->Sumw2();
		Double_t ScaleFactor = 1.0;
		if (Scale) {
			ScaleFactor = hNum->Integral(1,NbinX[0], 1,NbinY[0])
						 /hDen->Integral(1,NbinX[0], 1,NbinY[0]);
			// cout << "Scaling MC by " << ScaleFactor;
		}
		SCL->Scale(ScaleFactor);


		// Find the correct uncertainties.  Default is sum-in-quad of relative
		// uncertainties; optionally, do a binomial ratio
		if (binomial) {
			// Two and a half weeks I spent working on this and all I can
			// say is that PROBABLY this is the right formula
			Double_t num, den, den2, wrong;
			int badN=0, badD=0;		Int_t binX,binY;
			for (Int_t iGlobal=0; iGlobal<=Ncells; ++iGlobal) {
				hNum->GetBin(iGlobal,binX,binY);

				num   = hNum->GetBinContent(iGlobal);
				den   = hDen->GetBinContent(iGlobal);
				den2  = hDen->GetBinError(iGlobal);
				den2 *= den2;

				if ( num > den ) {
					RAT->SetBinContent(iGlobal, 0.0);
					RAT->SetBinError  (iGlobal, 0.0);
					if ( badN==0 ) {
						cout << "N(numerator) > N(denominator) ";
						cout << "starting at bin (" << binX << ", " << binY << ")" << endl;
						badN = iGlobal;
					}
				}
				if ( den <= 0.0 ) {
					RAT->SetBinContent(iGlobal, 0.0);
					RAT->SetBinError(iGlobal, 0.0);
					if ( !(RAT->IsBinUnderflow(iGlobal)) && !(RAT->IsBinOverflow(iGlobal))
						 && badD==0 ) {
						cout << "N(denominator) <= 0 ";
						cout << "starting at bin (" << binX << ", " << binY << ")" << endl;
						badD = iGlobal;
					}
				}
				if ( num<=den && den>0 ) {
					wrong = den2*(num/den)*(den-num)/den;
					wrong = sqrt(wrong) / den;
					RAT->SetBinError(iGlobal,wrong);
				}
			}
		}


		// Plot the dang thing.
		rawCanvas->Clear();
		rawCanvas->Divide(1,2);
		rawCanvas->cd(1);			// Plot the MC over the data on top
		Double_t macks;

		if ( hNum->GetMaximum() > SCL->GetMaximum() ){
			macks = hNum->GetMaximum();
		} else {
			macks = SCL->GetMaximum();
		}
		if (macks > 0.0 && Log) {
			gPad->SetLogz(1);
		} else {
			gPad->SetLogz(0);
		}
		hNum->SetMaximum(1.2*macks);
		hNum->DrawCopy("ZCOL");
		gStyle->SetOptStat("ou");
		gPad->Update();

		rawCanvas->cd(2);
		if (macks > 0.0 && Log) {
			gPad->SetLogz(1);
		} else {
			gPad->SetLogz(0);
		}
		SCL->SetMaximum(1.2*macks);
		SCL->DrawCopy("ZCOL");
		gStyle->SetOptStat("ou");
		gPad->Update();		

		ratioCanvas->Clear();		// Show the ratio on the other canvas
		ratioCanvas->cd();
		hNum->SetMaximum(macks);
		if (Log) {
			gPad->SetLogz(1);
		} else {
			gPad->SetLogz(0);
		}
		if (Opts) {
			RAT->DrawCopy(Opts);
		} else {
			RAT->DrawCopy("ZCOL");
		}
		gStyle->SetOptStat("");
		gPad->Update();


		Double_t Chi2,pull,err2;	// Print quality of fit message
		Chi2 = 0.0;
		Int_t NDoF = 0;
		if (Scale) ++NDoF;
		Double_t Chi2OU;	Chi2OU = 0.0;
		Int_t NDoFOU = 0;
		Int_t binX,binY;

		// Compute inside the histogram vs under & overflow bins separately
		for (Int_t iGlobal=0; iGlobal<=Ncells; ++iGlobal) {
			if (hNum->GetBinContent(iGlobal)!=0.0 || SCL->GetBinContent(iGlobal)!=0.0) {
				pull =  hNum->GetBinContent(iGlobal) -SCL->GetBinContent(iGlobal);
				err2 =  hNum->GetBinError(iGlobal)*hNum->GetBinError(iGlobal)
						+SCL->GetBinError(iGlobal)*SCL->GetBinError(iGlobal);
				//cout << "Bin " << iGlobal << " Pull " << pull/sqrt(err2) << endl;
				
				if ( !(hNum->IsBinUnderflow(iGlobal)) && !(hNum->IsBinOverflow(iGlobal)) ) {			
					if (err2<=0.0) {
						hNum->GetBin(iGlobal,binX,binY);
						cout << "Zero error, non-zero content at bin (" << binX << ", "
							<< binY << ")" << endl;
						err2 = hNum->GetBinContent(iGlobal) +SCL->GetBinContent(iGlobal);			
					}
					Chi2 += pull*pull/err2;
					++NDoF;
				} else {
					if (err2<=0.0) {
						hNum->GetBin(iGlobal,binX,binY);
						cout << "Zero error, non-zero content at bin (" << binX << ", "
							<< binY << ")" << endl;
						cout << "Which is OVER/UNDERFLOW" << endl;
						err2 = hNum->GetBinContent(iGlobal) +SCL->GetBinContent(iGlobal);			
					}
					Chi2OU += pull*pull/err2;
					++NDoFOU;
				}
			}
		}
		
		cout << "\nPlotRatio: " << hNum->GetName() << " / " << hDen->GetName() << "\n";
		cout << "Chi2/NDoF: " << Chi2 << " / " << NDoF << "   C.L.: " 
			<< TMath::Prob(Chi2,NDoF) << endl;


		cout << "Chi2/NDoF: " << Chi2+Chi2OU << " / " << NDoF+NDoFOU << "   C.L.: " 
			 << TMath::Prob(Chi2+Chi2OU,NDoF+NDoFOU);
		cout << "  with over/underflows\n";
		return 0;
	}
};
