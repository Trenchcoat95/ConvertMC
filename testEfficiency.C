/*
    gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/aliKalman/test/\"")
    gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/MC/\"")
    gSystem->Load("$fastMCKalman/fastMCKalman/aliKalman/test/AliExternalTrackParam.so");
    .L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx++g
    .L $ConvertMC/testEfficiency.C++g
    AliPDG::AddParticlesToPdgDataBase();
    testEfficiency(10000,kTRUE,0,13)
*/

#include "fastSimulation.h"
#include "fastTrackerGAR.h"
#include "TTreeStream.h"
#include "TStopwatch.h"
#include "TRandom.h"
#include "TChain.h"
#include "AliXRDPROOFtoolkit.h"
//#include "AliDrawStyle.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TSystem.h"
#include "TPad.h"
#include "TCanvas.h"
#include "AliPID.h"
#include <algorithm>
#include <iostream>
#include <string>
#include "garutils.h"
#include "TH1F.h"
#include "TEfficiency.h"

const Float_t kDecayFraction=0.5;
const Float_t kRandomPDGFraction=0.5;

TChain * treeFast = 0;
TChain * treeTurn=0;
TChain * treeUnit0=0;
TChain * treeSeed=0;

double GArCenter[3]={0,-150.473,1486}; 

void testEfficiency(Int_t nEv, bool dumpStream=1, size_t FirstEvent=0, int ptype = 13)
{
TStopwatch timer;
  timer.Start();
  





  /////////////////////////////////////////////////////////////////////////////////Read and set the tree source
  TChain* tree_source = new TChain("/anatree/GArAnaTree");    
  tree_source->Add("/home/federico/Documents/Universita/Federico_2020-2021/Aliwork/ConvertMC/MC/Interaction_CustomHit/*");

  std::vector<size_t>  *TPCClusterTrkIDNumber = 0;
  std::vector<int>  *TPCClusterMCindex = 0;
  std::vector<size_t>  *TrackIDNumber = 0;
  std::vector<int>  *TrackMCindex = 0;
  std::vector<int>  *TrajMCPIndex = 0;
  //std::vector<int>  *MCTrkID = 0;
  std::vector<int>     *PDG = 0;
  std::vector<int>     *PDGMother = 0;
  std::vector<float>   *TPCClusterX = 0;
  std::vector<float>   *TPCClusterY = 0;
  std::vector<float>   *TPCClusterZ = 0;
  std::vector<float>   *TPCClusterCovXX = 0;
  std::vector<float>   *TPCClusterCovYY = 0;
  std::vector<float>   *TrajMCPPX = 0;
  std::vector<float>   *TrajMCPPY = 0;
  std::vector<float>   *TrajMCPPZ = 0;
  std::vector<float>   *TrajMCPX = 0;
  std::vector<float>   *TrajMCPY = 0;
  std::vector<float>   *TrajMCPZ = 0;
  std::vector<float> *MCPStartX = 0;
  std::vector<float> *MCPStartY = 0;
  std::vector<float> *MCPStartZ     =0;
  std::vector<float> *MCPStartPX    =0;
  std::vector<float> *MCPStartPY    =0;
  std::vector<float> *MCPStartPZ    =0;
  std::vector<float> *MCPEndX       =0;
  std::vector<float> *MCPEndY       =0;
  std::vector<float> *MCPEndZ       =0;
  std::vector<float> *MCPEndPX      =0;
  std::vector<float> *MCPEndPY      =0;
  std::vector<float> *MCPEndPZ      =0;
  std::vector<float> *TrackStartX = 0;
  std::vector<float> *TrackStartY = 0;
  std::vector<float> *TrackStartZ     =0;
  std::vector<float> *TrackStartPX    =0;
  std::vector<float> *TrackStartPY    =0;
  std::vector<float> *TrackStartPZ    =0;
  std::vector<int>   *TrackStartQ     = 0;
  std::vector<float> *TrackEndX       =0;
  std::vector<float> *TrackEndY       =0;
  std::vector<float> *TrackEndZ       =0;
  std::vector<float> *TrackEndPX      =0;
  std::vector<float> *TrackEndPY      =0;
  std::vector<float> *TrackEndPZ      =0;
  std::vector<int>   *TrackEndQ       =0;
  std::vector<int>   *NTPCClustersOnTrack=0;
  std::vector<float> *VertX      =0;
  std::vector<float> *VertY      =0;
  std::vector<float> *VertZ      =0;
  std::vector<float> *MCVertX      =0;
  std::vector<float> *MCVertY      =0;
  std::vector<float> *MCVertZ      =0;
  std::vector<int>   *GPartPdg   =0;
  std::vector<int>   *CCNC   =0;
  TBranch        *b_TPCClusterTrkIDNumber;   //!
  TBranch        *b_TPCClusterMCindex;
  TBranch        *b_TrackIDNumber;
  TBranch        *b_TrackMCindex;
  TBranch        *b_TrajMCPIndex;
  //TBranch        *b_MCTrkID;
  TBranch        *b_PDG;
  TBranch        *b_PDGMother;
  TBranch        *b_TPCClusterX;   //!
  TBranch        *b_TPCClusterY;   //!
  TBranch        *b_TPCClusterZ;   //!
  TBranch        *b_TPCClusterCovXX;   //!
  TBranch        *b_TPCClusterCovYY;   //!
  TBranch        *b_TrajMCPPX;   //!
  TBranch        *b_TrajMCPPY;   //!
  TBranch        *b_TrajMCPPZ;   //!
  TBranch        *b_TrajMCPX;   //!
  TBranch        *b_TrajMCPY;   //!
  TBranch        *b_TrajMCPZ;   //!
  TBranch *b_TrackStartX ;
  TBranch *b_TrackStartY ;
  TBranch *b_TrackStartZ     ;
  TBranch *b_TrackStartPX    ;
  TBranch *b_TrackStartPY    ;
  TBranch *b_TrackStartPZ    ;
  TBranch  *b_TrackStartQ     ;
  TBranch *b_TrackEndX       ;
  TBranch *b_TrackEndY       ;
  TBranch *b_TrackEndZ       ;
  TBranch *b_TrackEndPX      ;
  TBranch *b_TrackEndPY      ;
  TBranch *b_TrackEndPZ      ;
  TBranch  *b_TrackEndQ       ;
  TBranch *b_NTPCClustersOnTrack;
  TBranch *b_MCPStartX ;
  TBranch *b_MCPStartY ;
  TBranch *b_MCPStartZ     ;
  TBranch *b_MCPStartPX    ;
  TBranch *b_MCPStartPY    ;
  TBranch *b_MCPStartPZ    ;
  TBranch *b_MCPEndX       ;
  TBranch *b_MCPEndY       ;
  TBranch *b_MCPEndZ       ;
  TBranch *b_MCPEndPX      ;
  TBranch *b_MCPEndPY      ;
  TBranch *b_MCPEndPZ      ;
  TBranch *b_VertX      ;
  TBranch *b_VertY      ;
  TBranch *b_VertZ      ;
  TBranch *b_MCVertX      ;
  TBranch *b_MCVertY      ;
  TBranch *b_MCVertZ      ;
  TBranch *b_GPartPdg   ;
  TBranch *b_CCNC       ;
  tree_source->SetBranchAddress("TPCClusterTrkIDNumber", &TPCClusterTrkIDNumber, &b_TPCClusterTrkIDNumber);
  tree_source->SetBranchAddress("TPCClusterMCindex", &TPCClusterMCindex, &b_TPCClusterMCindex);
  tree_source->SetBranchAddress("TrackIDNumber", &TrackIDNumber, &b_TrackIDNumber);
  tree_source->SetBranchAddress("TrackMCindex", &TrackMCindex, &b_TrackMCindex);
  tree_source->SetBranchAddress("TrajMCPIndex", &TrajMCPIndex, &b_TrajMCPIndex);
  //tree_source->SetBranchAddress("MCTrkID", &MCTrkID, &b_MCTrkID);
  tree_source->SetBranchAddress("PDG", &PDG, &b_PDG);
  tree_source->SetBranchAddress("PDGMother", &PDGMother, &b_PDGMother);
  tree_source->SetBranchAddress("TPCClusterX", &TPCClusterX, &b_TPCClusterX);
  tree_source->SetBranchAddress("TPCClusterY", &TPCClusterY, &b_TPCClusterY);
  tree_source->SetBranchAddress("TPCClusterZ", &TPCClusterZ, &b_TPCClusterZ);
  tree_source->SetBranchAddress("TPCClusterCovXX", &TPCClusterCovXX, &b_TPCClusterCovXX);
  tree_source->SetBranchAddress("TPCClusterCovYY", &TPCClusterCovYY, &b_TPCClusterCovYY);
  tree_source->SetBranchAddress("TrajMCPPX", &TrajMCPPX, &b_TrajMCPPX);
  tree_source->SetBranchAddress("TrajMCPPY", &TrajMCPPY, &b_TrajMCPPY);
  tree_source->SetBranchAddress("TrajMCPPZ", &TrajMCPPZ, &b_TrajMCPPZ);
  tree_source->SetBranchAddress("TrajMCPX", &TrajMCPX, &b_TrajMCPX);
  tree_source->SetBranchAddress("TrajMCPY", &TrajMCPY, &b_TrajMCPY);
  tree_source->SetBranchAddress("TrajMCPZ", &TrajMCPZ, &b_TrajMCPZ);
  tree_source->SetBranchAddress("TrackStartX", &TrackStartX, &b_TrackStartX);
  tree_source->SetBranchAddress("TrackStartY", &TrackStartY, &b_TrackStartY);
  tree_source->SetBranchAddress("TrackStartZ", &TrackStartZ, &b_TrackStartZ);
  tree_source->SetBranchAddress("TrackStartPX", &TrackStartPX, &b_TrackStartPX);
  tree_source->SetBranchAddress("TrackStartPY", &TrackStartPY, &b_TrackStartPY);
  tree_source->SetBranchAddress("TrackStartPZ", &TrackStartPZ, &b_TrackStartPZ);
  tree_source->SetBranchAddress("TrackStartQ", &TrackStartQ, &b_TrackStartQ);
  tree_source->SetBranchAddress("TrackEndX", &TrackEndX, &b_TrackEndX);
  tree_source->SetBranchAddress("TrackEndY", &TrackEndY, &b_TrackEndY);
  tree_source->SetBranchAddress("TrackEndZ", &TrackEndZ, &b_TrackEndZ);
  tree_source->SetBranchAddress("TrackEndPX", &TrackEndPX, &b_TrackEndPX);
  tree_source->SetBranchAddress("TrackEndPY", &TrackEndPY, &b_TrackEndPY);
  tree_source->SetBranchAddress("TrackEndPZ", &TrackEndPZ, &b_TrackEndPZ);
  tree_source->SetBranchAddress("TrackEndQ", &TrackEndQ, &b_TrackEndQ);
  tree_source->SetBranchAddress("NTPCClustersOnTrack", &NTPCClustersOnTrack, &b_NTPCClustersOnTrack);
  tree_source->SetBranchAddress("MCPStartX", &MCPStartX, &b_MCPStartX);
  tree_source->SetBranchAddress("MCPStartY", &MCPStartY, &b_MCPStartY);
  tree_source->SetBranchAddress("MCPStartZ", &MCPStartZ, &b_MCPStartZ);
  tree_source->SetBranchAddress("MCPStartPX", &MCPStartPX, &b_MCPStartPX);
  tree_source->SetBranchAddress("MCPStartPY", &MCPStartPY, &b_MCPStartPY);
  tree_source->SetBranchAddress("MCPStartPZ", &MCPStartPZ, &b_MCPStartPZ);
  tree_source->SetBranchAddress("MCPEndX", &MCPEndX, &b_MCPEndX);
  tree_source->SetBranchAddress("MCPEndY", &MCPEndY, &b_MCPEndY);
  tree_source->SetBranchAddress("MCPEndZ", &MCPEndZ, &b_MCPEndZ);
  tree_source->SetBranchAddress("MCPEndPX", &MCPEndPX, &b_MCPEndPX);
  tree_source->SetBranchAddress("MCPEndPY", &MCPEndPY, &b_MCPEndPY);
  tree_source->SetBranchAddress("MCPEndPZ", &MCPEndPZ, &b_MCPEndPZ);
  tree_source->SetBranchAddress("VertX", &VertX, &b_VertX);
  tree_source->SetBranchAddress("VertY", &VertY, &b_VertY);
  tree_source->SetBranchAddress("VertZ", &VertZ, &b_VertZ);
  tree_source->SetBranchAddress("MCVertX", &MCVertX, &b_MCVertX);
  tree_source->SetBranchAddress("MCVertY", &MCVertY, &b_MCVertY);
  tree_source->SetBranchAddress("MCVertZ", &MCVertZ, &b_MCVertZ);
  tree_source->SetBranchAddress("GPartPdg", &GPartPdg, &b_GPartPdg);
  tree_source->SetBranchAddress("CCNC", &CCNC, &b_CCNC);

  int nEvSource = tree_source->GetEntries();
  nEv = std::min(nEv,nEvSource);
  
  



  Int_t nPrimariesMCtot = 0;
  Int_t nPrimariesRecotot = 0;
  TFile* f = new TFile("Eficiency.root","RECREATE");
  TH1F MCptot_MC("MCptot_MC","MCptot_MC",20,0,6);
  TH1F MCptot_Reco("MCptot_Reco","MCptot_Reco",20,0,6);

  for (Int_t i=FirstEvent; i<nEv; i++){
 
      tree_source->GetEntry(i);
      if(TrackIDNumber->size()==0) continue;

      if(CCNC!=0 && GPartPdg->at(0)!=14) continue; //If analysing Interaction file only consider numuCC events   

      if(VertX->size()==0) continue;
      bool vert_in_fid = inFiducial(MCVertZ->at(0)-GArCenter[2],MCVertY->at(0)-GArCenter[1],MCVertX->at(0)-GArCenter[0]);
      if(!vert_in_fid) continue; //If analysng realistic Interaction sample, apply fiducial cut
      
      Int_t nPrimariesMC = 0;
      for(Int_t t=0; t<PDG->size(); t++){
        if(PDG->at(t)==ptype && PDGMother->at(t)==0){
          nPrimariesMC++;
          Float_t pMC = sqrt(pow(MCPStartPX->at(t),2)+pow(MCPStartPY->at(t),2)+pow(MCPStartPY->at(t),2));
          MCptot_MC.Fill(pMC);
        }
      }
 
      Int_t nPrimariesReco = 0;
      for(UInt_t k=0; k<TrackIDNumber->size();k++){
         size_t MCID = TrackMCindex->at(k);
         if(PDG->at(MCID)==ptype && PDGMother->at(MCID)==0){
          nPrimariesReco++;
          Float_t pMC = sqrt(pow(MCPStartPX->at(MCID),2)+pow(MCPStartPY->at(MCID),2)+pow(MCPStartPY->at(MCID),2));
          MCptot_Reco.Fill(pMC);
        }
      }
      std::cout<<"nEv: "<<i<<" "<<tree_source->GetFile()->GetName()<<std::endl;

    
      
  }

  MCptot_MC.Write();
  MCptot_Reco.Write();
  TH1F Ratio_MCptot = MCptot_Reco/MCptot_MC;
  Ratio_MCptot.SetName("Ratio_MCptot");
  Ratio_MCptot.SetTitle("Ratio_MCptot");
  Ratio_MCptot.Write();
  
}