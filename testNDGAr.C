/*

    gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/aliKalman/test/\"")
    gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/MC/\"")
    gSystem->Load("$fastMCKalman/fastMCKalman/aliKalman/test/AliExternalTrackParam.so");
    .L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx++g
    .L $ConvertMC/testNDGAr.C++g

    gSystem->AddIncludePath("-I\"test/\"")
    gSystem->Load("test/AliExternalTrackParam.so");
    .L fastSimulation.cxx++g
    .L testNDGAr.C++g
    
    testNDGAr(50000,"$ConvertMC/MC/Interaction/*",kTRUE,kFALSE,0,kTRUE,kTRUE,0,kTRUE) ////for TKI
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
#include<sstream>
#include "tools/garutils.h"
const Float_t kDecayFraction=0.5;
const Float_t kRandomPDGFraction=0.5;

TChain * treeFast = 0;
TChain * treeTurn=0;
TChain * treeUnit0=0;
TChain * treeSeed=0;


void testNDGAr(Int_t nEv, std::string file, bool dumpStream=1, bool Ideal=kTRUE, size_t FirstEvent=0, bool Interaction=kFALSE, bool RotFrame = kFALSE, double ptype = 13, bool TKISample=false, bool MoveToVertex=false){


  const Int_t   nLayerTPC=278;
  const Int_t   nPoints=nLayerTPC*3;
  const Float_t xx0=8.37758e-04; //1/X0 cm^-1 for ArCH4 at 10 atm
  const Float_t xrho=0.016770000; //rho g/cm^3 for ArCH4 at 10 atm
  const Float_t kMaterialScaling=1;      ////Promote to global variable
  double GArCenter[3]={0,-150.473,1486}; 
  float fSortDistCut = 25.0;
  float fPrintLevel = 0.0;
  const Float_t kMaxResol=0.5;            // point resolution scan max resolution
  const Float_t kMinResol=0.05;           //  point resolution scan min resolution
            


  TStopwatch timer;
  timer.Start();
  

  //////////////////////////////////////////////////////////////////////////////////Create result File
  std::string filename = "FullReco";
  if(Interaction) filename+="Interaction";
  if(Ideal) filename+= "Ideal";
  if(TKISample) filename+="TKI";
  if(MoveToVertex) filename+="Vertex";
  std::string path ="";
  std::string dirname = path+filename;

  std::stringstream stream;
  stream.precision(0);
  stream << std::fixed;
  stream<<ptype;

  if(ptype!=0) filename += stream.str();
  filename+= "_first"+std::to_string(FirstEvent)+"_last"+std::to_string(nEv);
  std::string totalname = dirname +"/"+filename+".root";
  std::string filelist = dirname + "/fastParticle.list";
  const char* cdir = const_cast<char*>(dirname.c_str());
  const char* ctotal = const_cast<char*>(totalname.c_str());
  std::cout<< "Filename: "<<totalname<<std::endl;
  gSystem->mkdir(cdir,1);
  std::ofstream fw(filelist, std::ofstream::out);
  fw << totalname ;
  fw.close();



  /////////////////////////////////////////////////////////////////////////////////Read and set the tree source
  TChain* tree_source = new TChain("/anatree/GArAnaTree");    
  tree_source->Add(file.c_str());
  //else tree_source->Add("/home/federico/Documents/Universita/Federico_2020-2021/Aliwork/ConvertMC/MC/anatree_range*");
  std::vector<float> *TrackLenF = 0;
  std::vector<float> *TrackLenB = 0;
  std::vector<float> *MCNuPx=0;
  std::vector<float> *MCNuPy=0;
  std::vector<float> *MCNuPz=0;
  std::vector<size_t>  *TPCClusterTrkIDNumber = 0;
  std::vector<int>     *TPCClusterMCindex = 0;
  std::vector<size_t>  *TrackIDNumber = 0;
  std::vector<int>     *TrackMCindex = 0;
  std::vector<int>  *TrajMCPIndex = 0;
  //std::vector<int>  *MCTrkID = 0;
  std::vector<int>     *PDG = 0;
  std::vector<int>     *TgtPDG = 0;
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
  TBranch        *b_TrackLenF;
  TBranch        *b_TrackLenB;
  TBranch        *b_MCNuPx;
  TBranch        *b_MCNuPy;
  TBranch        *b_MCNuPz;
  TBranch        *b_TPCClusterTrkIDNumber;   //!
  TBranch        *b_TPCClusterMCindex;
  TBranch        *b_TrackIDNumber;
  TBranch        *b_TrackMCindex;
  TBranch        *b_TrajMCPIndex;
  //TBranch        *b_MCTrkID;
  TBranch        *b_PDG;
  TBranch        *b_TgtPDG;
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
  tree_source->SetBranchAddress("MCNuPx", &MCNuPx, &b_MCNuPx);
  tree_source->SetBranchAddress("MCNuPy", &MCNuPy, &b_MCNuPy);
  tree_source->SetBranchAddress("MCNuPz", &MCNuPz, &b_MCNuPz);
  tree_source->SetBranchAddress("TrackLenF", &TrackLenF, &b_TrackLenF);
  tree_source->SetBranchAddress("TrackLenB", &TrackLenB, &b_TrackLenB);
  tree_source->SetBranchAddress("TPCClusterTrkIDNumber", &TPCClusterTrkIDNumber, &b_TPCClusterTrkIDNumber);
  tree_source->SetBranchAddress("TPCClusterMCindex", &TPCClusterMCindex, &b_TPCClusterMCindex);
  tree_source->SetBranchAddress("TrackIDNumber", &TrackIDNumber, &b_TrackIDNumber);
  tree_source->SetBranchAddress("TrackMCindex", &TrackMCindex, &b_TrackMCindex);
  tree_source->SetBranchAddress("TrajMCPIndex", &TrajMCPIndex, &b_TrajMCPIndex);
  //tree_source->SetBranchAddress("MCTrkID", &MCTrkID, &b_MCTrkID);
  tree_source->SetBranchAddress("PDG", &PDG, &b_PDG);
  tree_source->SetBranchAddress("TgtPDG", &TgtPDG, &b_TgtPDG);
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
  
  




  /////////////////////////////////////////////////////////////////////////////////Build the detector geometry
  fastGeometry geom(nLayerTPC+1);
  geom.fBz=-5;
  float resol[2]={0.0001,0.0001};
  resol[0]=0.4;                   
  resol[1]=0.3;
  geom.setLayerRadiusPower(0,nLayerTPC,1,nLayerTPC,1.0,xx0,xrho,resol);
  for (size_t iLayer=0; iLayer<geom.fLayerX0.size();iLayer++) {
        geom.fLayerX0[iLayer] = xx0 * kMaterialScaling;
        geom.fLayerRho[iLayer] = xrho * kMaterialScaling;
        geom.fLayerResolRPhi[iLayer] =resol[0];
        geom.fLayerResolZ[iLayer] = resol[1];
      }



  ////////////////////////////////////////////////////////////////////////////////Create and reconstruct Particle
  TTreeSRedirector *pcstream = new TTreeSRedirector(ctotal,"recreate");
  TTree * tree = 0;
  for (Int_t i=FirstEvent; i<nEv; i++){
 
      tree_source->GetEntry(i);
      if(TrackIDNumber->size()==0) continue;

      if(Interaction && (CCNC->at(0)!=0 || GPartPdg->at(0)!=14)) continue; //If analysing Interaction file only consider numuCC events   

      //if(Interaction && !Ideal && VertX->size()==0) continue;
      if(Interaction && !Ideal) {
        bool vert_in_fid = inFiducial(MCVertX->at(0)-GArCenter[0],MCVertY->at(0)-GArCenter[1],MCVertZ->at(0)-GArCenter[2]);
        if(!vert_in_fid) continue; //If analysng realistic Interaction sample, apply fiducial cut
      }
      //if(MCNuPx->Size()==0)

      int Target=0;
      float NuPx=0;
      float NuPy=0;
      float NuPz=0;
      if(Interaction){
        Target = TgtPDG->at(0);
        NuPx = MCNuPx->at(0);
        NuPy = MCNuPy->at(0);
        NuPz = MCNuPz->at(0);
      }

      float VertGArx=0;
      float VertGAry=0;
      float VertGArz=0;
      if(Interaction && VertX->size()!=0){
        Target = TgtPDG->at(0);
        VertGArx=VertZ->at(0)-GArCenter[2];
        VertGAry=VertY->at(0)-GArCenter[1];
        VertGArz=VertX->at(0)-GArCenter[0];
      }

      int nPrimaries = 0;
      if(Interaction){
        for(int p=0;p<PDGMother->size();p++){
          if (PDGMother->at(p)==0) nPrimaries++;
          else break;
        }
      }


      size_t ID = 0;
      size_t MCID = 0;
      size_t NTPC = 0;
      Float_t CovXX = 0;
      Float_t CovYY = 0;
      std::vector<size_t> ID_vector;                    ///////For Each event find all the reconstructed tracks and save their IDs
      std::vector<int> MCID_vector;                     ///////For each track find the corresponding MC index to match information with MC truth
      std::vector<int> NTPCClusters_vector;                ///////For each track record the NTPCClusters used by garsoft
      std::vector<std::vector<TVector3>> TrksXYZ;       ///////Vector containing one vector of TPCCluster XYZ points for each track
      std::vector<TVector3> TrkClusterXYZ;              /////Container that will be filled for each track, added to TrksXYZ and then deleted
      std::vector<std::vector<float>> TrksCovXX;       ///////Vector containing one vector of TPCCluster XYZ points for each track
      std::vector<float> TrkClusterCovXX;              /////Container that will be filled for each track, added to TrksXYZ and then deleted
      std::vector<std::vector<float>> TrksCovYY;       ///////Vector containing one vector of TPCCluster XYZ points for each track
      std::vector<float> TrkClusterCovYY;              /////Container that will be filled for each track, added to TrksXYZ and then deleted
                                                        /////Note: for now all is in ND-GAr coordinates
      for(UInt_t k=0; k<TrackIDNumber->size();k++)
      {
          ID = TrackIDNumber->at(k);
          MCID = TrackMCindex->at(k);
          NTPC = NTPCClustersOnTrack->at(k);
          for(UInt_t j=0; j<TPCClusterTrkIDNumber->size();j++)
          {
            TVector3 xyz(TPCClusterX->at(j),TPCClusterY->at(j),TPCClusterZ->at(j));
            if(ID==TPCClusterTrkIDNumber->at(j)) {
              TrkClusterXYZ.push_back(xyz);
              TrkClusterCovXX.push_back(TPCClusterCovXX->at(j));
              TrkClusterCovYY.push_back(TPCClusterCovYY->at(j));
            }
          }
          ID_vector.push_back(ID);
          MCID_vector.push_back(MCID);
          TrksXYZ.push_back(TrkClusterXYZ);
          TrksCovXX.push_back(TrkClusterCovXX);
          TrksCovYY.push_back(TrkClusterCovYY);
          NTPCClusters_vector.push_back(NTPC);
          TrkClusterXYZ.clear();
          TrkClusterCovXX.clear();
          TrkClusterCovYY.clear();
      }

      //if(TKISample && (TrksXYZ.size()!=3 || Target!=1000010010)) continue;
      if(TKISample && (TrksXYZ.size()!=3)) continue;
      
      for(size_t t=0; t<TrksXYZ.size(); t++)   ////Cycle over all the tracks for the event 
                                               ////Note that in the final tree the separation between events is lost and each item is a track)
      {    
        if(Ideal){
          Float_t resolScan= kMinResol + gRandom->Rndm()*(kMaxResol-kMinResol);
          for (size_t iLayer=0; iLayer<geom.fLayerX0.size();iLayer++) {
            geom.fLayerResolRPhi[iLayer] = resolScan;
            geom.fLayerResolZ[iLayer] = resolScan;
          }
        }
        if(!Interaction && t>0) continue; //If analyzing particle gun file, only consider first track
        if(Interaction && Ideal && t>0) continue; //If analyzing particle gun file, only consider first track
        //if(t!=7) continue;
        if(MCID_vector.at(t)>=PDG->size()) continue;
        long PDGcode=PDG->at(MCID_vector.at(t));    ////Save PDGCode for the MC trajectory associated with the ND-GAr reconstructed track
        long PDGMothercode=PDGMother->at(MCID_vector.at(t));    ////Save PDGCode for the MC trajectory associated with the ND-GAr reconstructed track
        // bool contained = inTPC(MCPStartZ->at(MCID_vector.at(t))-GArCenter[2],MCPStartY->at(MCID_vector.at(t))-GArCenter[1],MCPStartX->at(MCID_vector.at(t))-GArCenter[0]);
        // contained *= inTPC(MCPEndZ->at(MCID_vector.at(t))-GArCenter[2],MCPEndY->at(MCID_vector.at(t))-GArCenter[1],MCPEndX->at(MCID_vector.at(t))-GArCenter[0]);
        if(!TKISample && (PDGcode!=ptype || PDGMothercode!=0) ) continue;


          ///////Match MC trajectory that produced the track and save  the xyz's and pxyz's
          std::vector<TVector3> trajxyz; 
          std::vector<float> trajxyzcovxx;
          std::vector<float> trajxyzcovyy;                
          std::vector<TVector3> trajpxyz;
          for(size_t u=0; u<TrajMCPX->size(); u++)          
          {
            if(MCID_vector.at(t)==TrajMCPIndex->at(u))
            {
              TVector3 jxyz(TrajMCPX->at(u),TrajMCPY->at(u),TrajMCPZ->at(u));
              TVector3 pjxyz(TrajMCPPX->at(u),TrajMCPPY->at(u),TrajMCPPZ->at(u));
              double Cluster[] = {TrksXYZ.at(t).at(0).X(),TrksXYZ.at(t).at(0).Y(),TrksXYZ.at(t).at(0).Z()};
              trajxyzcovxx.push_back(0.4);
              trajxyzcovxx.push_back(0.4);
              trajxyz.push_back(jxyz);
              trajpxyz.push_back(pjxyz);
            }
          }

          if(trajxyz.size()==0) continue;       /////Skip tracks that were not successfully associated to a MC trajectory by the back-tracker

          

          if(!Interaction && abs(PDGcode)!=13 ) continue;   ////If analyzing partigle gun file, only consider muons 
          TParticlePDG *particle = TDatabasePDG::Instance()->GetParticle(PDGcode); ///Check for invalid pdgCode
          if (particle == nullptr) continue;
          std::cout<<"ID: "<<ID_vector.at(t)<<" PDG: "<<PDGcode<<" ev: "<< i <<std::endl;




          ///////////////////////////////////////////// Sort TPC Clusters using the ND-GAr method
          std::vector<int> hlf;
          std::vector<int> hlb;
          float lftmp=0;  
          float lbtmp=0;
          sort_TPCClusters_along_track2(TrksXYZ.at(t),hlf,hlb,fPrintLevel,lftmp,lbtmp,fSortDistCut);

          if(hlf.size()!=TrksXYZ.at(t).size()){
            size_t hlfsize = hlf.size();
            size_t trksize = TrksXYZ.at(t).size();
            size_t NTPClus = NTPCClusters_vector[t]; 
            std::cout<<"Ordering failed: "<<hlf.size()<<" != "<<TrksXYZ.at(t).size()<<std::endl;
          }

          std::vector<TVector3> TrkClusterXYZb_NDGAr;
          std::vector<TVector3> TrkClusterXYZf_NDGAr;
          std::vector<float> TrkClusterCovXXb_NDGAr;
          std::vector<float> TrkClusterCovXXf_NDGAr;
          std::vector<float> TrkClusterCovYYb_NDGAr;
          std::vector<float> TrkClusterCovYYf_NDGAr;
          for(size_t k=0;k<hlb.size();k++) {
            TrkClusterXYZb_NDGAr.push_back(TrksXYZ.at(t).at(hlb[k]));
            TrkClusterCovXXb_NDGAr.push_back(TrksCovXX.at(t).at(hlb[k]));
            TrkClusterCovYYb_NDGAr.push_back(TrksCovYY.at(t).at(hlb[k]));
          }
          for(size_t k=0;k<hlf.size();k++) {
            TrkClusterXYZf_NDGAr.push_back(TrksXYZ.at(t).at(hlf[k]));
            TrkClusterCovXXf_NDGAr.push_back(TrksCovXX.at(t).at(hlf[k]));
            TrkClusterCovYYf_NDGAr.push_back(TrksCovYY.at(t).at(hlf[k]));
          }

          double stX = TrkClusterXYZf_NDGAr[0].X();
          double stY = TrkClusterXYZf_NDGAr[0].Y();
          double stZ = TrkClusterXYZf_NDGAr[0].Z();
          double endX = TrkClusterXYZb_NDGAr[0].X();
          double endY = TrkClusterXYZb_NDGAr[0].Y();
          double endZ = TrkClusterXYZb_NDGAr[0].Z();

          ///Now that we have vectors of sorted TPCClusters we can find the closest trajectory point for each cluster 
          ///and save the MC pxyz information for each point
          ///Note that if IsForwardf=kTRUE it means that the forward ordering produced by sort_TPCClusters_along_track2 
          ///matches the actual order of the trajectory points

          Bool_t IsForwardf, isForwardb;
          std::vector<TVector3> trajpxyzb_NDGAr;
          std::vector<TVector3> trajxyzb_NDGAr;
          ClosestTrajectory(trajxyz,trajpxyz,TrkClusterXYZb_NDGAr,trajxyzb_NDGAr,trajpxyzb_NDGAr,isForwardb);
          std::vector<TVector3> trajpxyzf_NDGAr;
          std::vector<TVector3> trajxyzf_NDGAr;
          ClosestTrajectory(trajxyz,trajpxyz,TrkClusterXYZf_NDGAr,trajxyzf_NDGAr,trajpxyzf_NDGAr,IsForwardf);

          ////Build end and start NDGAr reco parameters to be saved 
          ///Start of the track
          AliExternalTrackParam4D paramSt0;
          TVector3 xyzSt(TrackStartX->at(t),TrackStartY->at(t),TrackStartZ->at(t));
          TVector3 pxyzSt(TrackStartPX->at(t),TrackStartPY->at(t),TrackStartPZ->at(t));
          int qSt = TrackStartQ->at(t);
          TVector3 xyzStMC(MCPStartX->at(MCID_vector.at(t)),MCPStartY->at(MCID_vector.at(t)),MCPStartZ->at(MCID_vector.at(t)));
          TVector3 pxyzStMC(MCPStartPX->at(MCID_vector.at(t)),MCPStartPY->at(MCID_vector.at(t)),MCPStartPZ->at(MCID_vector.at(t)));
          int qStMC = 1;
          if(abs(PDGcode)==13) qStMC = (PDGcode>0)?-1:1;
          if(abs(PDGcode)==211) qStMC = (PDGcode<0)?-1:1;
          if(abs(PDGcode)==2212) qStMC = (PDGcode<0)?-1:1;
          BuildParamNDGArReco(paramSt0,GArCenter,xyzSt,pxyzSt,PDGcode,qSt);
          AliExternalTrackParam4D paramStMC;
          BuildParamNDGArReco(paramStMC,GArCenter,xyzStMC,pxyzStMC,PDGcode,qStMC);


          ///End of the track
          AliExternalTrackParam4D paramEnd0;
          TVector3 xyzEnd(TrackEndX->at(t),TrackEndY->at(t),TrackEndZ->at(t));
          TVector3 pxyzEnd(TrackEndPX->at(t),TrackEndPY->at(t),TrackEndPZ->at(t));
          int qEnd = TrackEndQ->at(t);
          TVector3 xyzEndMC(MCPEndX->at(MCID_vector.at(t)),MCPEndY->at(MCID_vector.at(t)),MCPEndZ->at(MCID_vector.at(t)));
          TVector3 pxyzEndMC(MCPEndPX->at(MCID_vector.at(t)),MCPEndPY->at(MCID_vector.at(t)),MCPEndPZ->at(MCID_vector.at(t)));
          int qEndMC=1;
          if(abs(PDGcode)==13) qEndMC = (PDGcode>0)?-1:1;
          if(abs(PDGcode)==211) qEndMC = (PDGcode<0)?-1:1;
          if(abs(PDGcode)==2212) qEndMC = (PDGcode<0)?-1:1;
          BuildParamNDGArReco(paramEnd0,GArCenter,xyzEnd,pxyzEnd,PDGcode,qEnd);
          AliExternalTrackParam4D paramEndMC;
          BuildParamNDGArReco(paramEndMC,GArCenter,xyzEndMC,pxyzEndMC,PDGcode,qStMC);

          AliExternalTrackParam4D paramSt;
          AliExternalTrackParam4D paramEnd;
          double TrackLenSt;
          double TrackLenEnd;

          //NB: In GArSoft what is labeled as TrackStart is positioned at the last TPCCluster,
          //following the sort_TPCClusters_along_track2 forward ordering.
          //If the forward ordering matches the trajectory ordering IsForwardf=kTRUE,
          //hence, what we would consider to be the actual track start (i.e. the extremity closest to the trajectory start) 
          //is labelled by GArSoft as TrackEnd and viceversa

          if(!IsForwardf) 
          {
            TrackLenSt = TrackLenB->at(t);
            TrackLenEnd = TrackLenF->at(t);
            paramSt = paramSt0;
            paramEnd = paramEnd0;
          }
          else
          {
            TrackLenSt = TrackLenF->at(t);
            TrackLenEnd = TrackLenB->at(t);
            paramSt = paramEnd0;
            paramEnd = paramSt0;
          }

          ////Build end and start non rotated NDGAr reco parameters to be saved 
          ///Start of the track
          AliExternalTrackParam4D paramStNoRot0;
          BuildParamNDGArRecoNoRotation(paramStNoRot0,GArCenter,xyzSt,pxyzSt,PDGcode,qSt);
          AliExternalTrackParam4D paramStNoRotMC;
          BuildParamNDGArRecoNoRotation(paramStNoRotMC,GArCenter,xyzStMC,pxyzStMC,PDGcode,qStMC);


          ///End of the track
          AliExternalTrackParam4D paramEndNoRot0;
          BuildParamNDGArRecoNoRotation(paramEndNoRot0,GArCenter,xyzEnd,pxyzEnd,PDGcode,qEnd);
          AliExternalTrackParam4D paramEndNoRotMC;
          BuildParamNDGArRecoNoRotation(paramEndNoRotMC,GArCenter,xyzEndMC,pxyzEndMC,PDGcode,qEndMC);

          AliExternalTrackParam4D paramStNoRot;
          AliExternalTrackParam4D paramEndNoRot;

          if(!IsForwardf) //If Trk points are "forward" ordered (i.e. same order as trajMC) by sort_TPCClusters_along_track2, GArSoft reconstructs the end; it's the opposite in our code
          {
            paramStNoRot = paramStNoRot0;
            paramEndNoRot = paramEndNoRot0;
          }
          else
          {
            paramStNoRot = paramEndNoRot0;
            paramEndNoRot = paramStNoRot0;
          }

          double checkloop =0;
          //if(TrkClusterXYZf_NDGAr.size()>50) continue;
          if(Ideal) {
            ////// Create particle object with ALICE-like global coordinates
            fastParticle particle(trajxyz.size()+1);
            particle.fAddMSsmearing=true;
            particle.fAddPadsmearing=true;
            particle.fUseMCInfo=false;
            particle.fUseGArSeeding=true;
            particle.fgStreamer=pcstream;
            particle.gid=ID_vector.at(t);
            particle.fDecayLength=0;
            double xstart=0;
            double ystart=0;
            if(RotFrame){
              size_t size_dis = 30;
              if (int(trajxyz.size()-1)<30) size_dis = int(trajxyz.size()-1);
              double displacex = (trajxyz.at(size_dis).Z()-trajxyz.at(0).Z());
              double displacey = (trajxyz.at(size_dis).Y()-trajxyz.at(0).Y());
              double displace_mod = sqrt(displacex*displacex+displacey*displacey);
              xstart = -(trajxyz.at(0).Z()-GArCenter[2])+20*displacex/displace_mod;
              ystart = -(trajxyz.at(0).Y()-GArCenter[1])+20*displacey/displace_mod; 
              double Length = 0;
              int Size = trajxyz.size();
              for(int i=1;i<trajxyz.size();i++){
                  Length+=sqrt(pow(trajxyz[i].X()-trajxyz[i-1].X(),2)+pow(trajxyz[i].Y()-trajxyz[i-1].Y(),2)+pow(trajxyz[i].Z()-trajxyz[i-1].Z(),2));
              }
              int Sample = 0.85 * Size/Length;
              if (Sample==0) Sample++;
              BuildParticle(particle,GArCenter,geom,trajxyz,trajpxyz,PDGcode,xstart,ystart,Sample);  //////Built the ALICE particle with the MC trajectory xyz and pxyz points 
              if(particle.fParamMC.size()==0) continue;
              particle.reconstructParticleFull(geom,PDGcode,10000);
              BuildInRotfromIn(particle);
              //particle.reconstructParticleFullOut(geom,PDGcode,10000);
              //particle.refitParticle();
            }else{
              BuildParticleNoRotation(particle,GArCenter,geom,trajxyz,trajpxyz,PDGcode,0);
              if(particle.fParamMC.size()==0) continue;

              particle.reconstructParticleRotate0(geom,PDGcode,10000);
              // if (particle.fParamInRot[0].GetCovariance()[0]==0){
              //   for(int a=1;a<=50;a++)
              //   {
              //     BuildParticleNoRotation(particle,GArCenter,geom,trajxyz,trajpxyz,PDGcode,a*Pi()/50);
              //     particle.reconstructParticleRotate0(geom,PDGcode,10000);
              //     if(particle.fParamInRot[0].GetCovariance()[0]!=0) {
              //       break;
              //     }
              //   }
              // }
            }


            if (dumpStream==kFALSE) continue;
            if (tree) tree->Fill();
            else {
              (*pcstream) << "fastPart" <<
                          "i=" << i <<
                          "t=" <<t<<
                          "xstart="<<xstart<<
                          "ystart="<<ystart<< 
                          //"paramGAr.=" <<&paramGAr<<
                          ///"paramGArbkw.=" <<&paramGArbkw<<
                          "paramSt.="<<&paramSt<<
                          "paramEnd.="<<&paramEnd<<
                          "paramStNoRot.="<<&paramStNoRot<<
                          "paramEndNoRot.="<<&paramEndNoRot<<
                          "paramStNoRotMC.="<<&paramStNoRotMC<<
                          "paramStMC.="<<&paramStMC<<
                          "paramEndMC.="<<&paramEndMC<<
                          "geom.="<<&geom<<
                          "part.=" << &particle <<
                          "\n";
              tree=  ((*pcstream) << "fastPart").GetTree();
            }
          }
          else
          {

            /////// Repeat reconstruction twice because of the energy loss reconstruction
            ////// Create particle object with ALICE-like global coordinates

            //std::cout<<"Forward ordering reconstruction"<<std::endl;
            fastParticle particle_f(hlf.size()+1);
            particle_f.fAddMSsmearing=true;
            particle_f.fAddPadsmearing=false;
            particle_f.fUseGArSeeding=false;
            particle_f.fUseMCInfo=false;
            particle_f.fgStreamer=pcstream;
            particle_f.gid=ID_vector.at(t);
            particle_f.fDecayLength=0;
            fastParticle particle_f_nosmear(hlf.size()+1);
            particle_f_nosmear.fAddMSsmearing=false;
            particle_f_nosmear.fAddPadsmearing=false;
            particle_f_nosmear.fUseGArSeeding=false;
            particle_f_nosmear.fUseMCInfo=false;
            particle_f_nosmear.fgStreamer=pcstream;
            particle_f_nosmear.gid=ID_vector.at(t);
            particle_f_nosmear.fDecayLength=0;
            double xstart=0;
            double ystart=0;
            int Flip = 1;
            double MCVertXloc=0;
            double MCVertYloc=0;
            double MCVertZloc=0;
            
            if(RotFrame){
              if(IsForwardf){
                std::cout<<"Forward"<<std::endl;
                //double sign = (TrkClusterXYZf_NDGAr.at(TrkClusterXYZf_NDGAr.size()-1).Z()-TrkClusterXYZf_NDGAr.at(0).Z())>0?1:-1;
                // for(int a=0;a<=500;a++)
                // {
                  size_t size_dis = 30;
                  if (int(TrkClusterXYZf_NDGAr.size()-1)<30) size_dis = int(TrkClusterXYZf_NDGAr.size()-1);
                  double displacex = (TrkClusterXYZf_NDGAr.at(size_dis).Z()-TrkClusterXYZf_NDGAr.at(0).Z());
                  double displacey = (TrkClusterXYZf_NDGAr.at(size_dis).Y()-TrkClusterXYZf_NDGAr.at(0).Y());
                  double displace_mod = sqrt(displacex*displacex+displacey*displacey);
                  xstart = -(TrkClusterXYZf_NDGAr.at(0).Z()-GArCenter[2])+20*displacex/displace_mod;
                  ystart = -(TrkClusterXYZf_NDGAr.at(0).Y()-GArCenter[1])+20*displacey/displace_mod;           
                  BuildParticle(particle_f,GArCenter,geom,TrkClusterXYZf_NDGAr,trajpxyzf_NDGAr,PDGcode,xstart,ystart);  /////Build the ALICE particle with Track XYZClusters ordered forward and closest MC pxyz
                  //SetParticleLoop(particle_f);
                  if(particle_f.fParamMC.size()==0) continue;
                  BuiltParticleUnsmeared(particle_f_nosmear,GArCenter,TrkClusterXYZf_NDGAr,trajxyzf_NDGAr,xstart,ystart);
                  particle_f.reconstructParticleFull(geom,PDGcode,10000);
                  particle_f.reconstructParticleFullOut(geom,PDGcode,10000);
                //   if(particle_f.fParamIn[0].GetCovariance()[0]!=0) {
                  BuildInRotfromIn(particle_f);
                  BuildParamNDGArReco(paramSt,GArCenter,xyzEnd,pxyzEnd,PDGcode,qEnd,xstart,ystart);
                  Flip = BuildParamNDGArReco(paramStMC,GArCenter,xyzStMC,pxyzStMC,PDGcode,qEnd,xstart,ystart);
                  double xyzVertex[3];
                  paramStMC.GetXYZ(xyzVertex);
                  double checkpxyz[3];
                  paramStMC.GetPxPyPz(checkpxyz);
                  if(MoveToVertex) {MoveParamToPoint(particle_f.fParamIn[0],xyzVertex, xx0, xrho, geom.fLayerResolRPhi[0], geom.fLayerResolZ[0], 
                                                    geom.fBz, 10000/geom.fLayerResolRPhi[0], particle_f.fMassRec);}
                  // double MCVertXloc=MCVertX->at(0)-(GArCenter[2]-xstart);
                  // double MCVertYloc=MCVertY->at(0)-(GArCenter[1]-ystart);
                  // double MCVertZloc=MCVertX->at(0)-GArCenter[0];

                //     break;
                //   }
                // }
                // for(int q=0;q<particle_f.fParamMC.size();q++){
                //   particle_f.fResolRPhi[q] = pow(particle_f.fParamMC[q].GetParameter()[0]-particle_f_nosmear.fParamMC[q].GetParameter()[0],2);
                //   particle_f.fResolZ[q] = pow(particle_f.fParamMC[q].GetParameter()[1]-particle_f_nosmear.fParamMC[q].GetParameter()[1],2);
                // }
                
              }
              else{
                // for(int a=0;a<=280;a++)
                // {
                std::cout<<"Backwards"<<std::endl;
                size_t size_dis = 30;
                if (int(TrkClusterXYZf_NDGAr.size()-1)<30) size_dis = int(TrkClusterXYZb_NDGAr.size()-1);
                double displacex = (TrkClusterXYZb_NDGAr.at(size_dis).Z()-TrkClusterXYZb_NDGAr.at(0).Z());
                double displacey = (TrkClusterXYZb_NDGAr.at(size_dis).Y()-TrkClusterXYZb_NDGAr.at(0).Y());
                double displace_mod = sqrt(displacex*displacex+displacey*displacey);
                xstart = -(TrkClusterXYZb_NDGAr.at(0).Z()-GArCenter[2])+20*displacex/displace_mod;
                ystart = -(TrkClusterXYZb_NDGAr.at(0).Y()-GArCenter[1])+20*displacey/displace_mod; 
                BuildParticle(particle_f,GArCenter,geom,TrkClusterXYZb_NDGAr,trajpxyzb_NDGAr,PDGcode,xstart,ystart);  /////Build the ALICE particle with Track XYZClusters ordered forward and closest MC pxyz
                //SetParticleLoop(particle_f);
                if(particle_f.fParamMC.size()==0) continue;
                BuiltParticleUnsmeared(particle_f_nosmear,GArCenter,TrkClusterXYZb_NDGAr,trajxyzb_NDGAr,xstart,ystart);
                particle_f.reconstructParticleFull(geom,PDGcode,10000);
                particle_f.reconstructParticleFullOut(geom,PDGcode,10000);
                //   if(particle_f.fParamIn[0].GetCovariance()[0]!=0) {
                BuildInRotfromIn(particle_f);
                BuildParamNDGArReco(paramSt,GArCenter,xyzSt,pxyzSt,PDGcode,qSt,xstart,ystart);
                Flip = BuildParamNDGArReco(paramStMC,GArCenter,xyzStMC,pxyzStMC,PDGcode,qSt,xstart,ystart);
                double xyzVertex[3];
                paramStMC.GetXYZ(xyzVertex);
                double checkpxyz[3];
                paramStMC.GetPxPyPz(checkpxyz);
                if(MoveToVertex) {MoveParamToPoint(particle_f.fParamIn[0],xyzVertex, xx0, xrho, geom.fLayerResolRPhi[0], geom.fLayerResolZ[0], 
                                                    geom.fBz, 10000/geom.fLayerResolRPhi[0], particle_f.fMassRec);}
                // double MCVertXloc=MCVertX->at(0)-(GArCenter[2]-xstart);
                // double MCVertYloc=MCVertY->at(0)-(GArCenter[1]-ystart);
                // double MCVertZloc=MCVertX->at(0)-GArCenter[0];
                //     break;
                //   }
                // }
                // for(int q=0;q<particle_f.fParamMC.size();q++){
                //   particle_f.fResolRPhi[q] = pow(particle_f.fParamMC[q].GetParameter()[0]-particle_f_nosmear.fParamMC[q].GetParameter()[0],2);
                //   particle_f.fResolZ[q] = pow(particle_f.fParamMC[q].GetParameter()[1]-particle_f_nosmear.fParamMC[q].GetParameter()[1],2);
                // }
                if(particle_f.fParamMC.size()==0) continue;
              }
              //particle_f.reconstructParticleFull(geom,PDGcode,10000);
              //BuildInRotfromIn(particle_f);
            }
            else{

            if(IsForwardf){
              std::cout<<"Forward"<<std::endl;
              BuildParticleNoRotation(particle_f,GArCenter,geom,TrkClusterXYZf_NDGAr,trajpxyzf_NDGAr,PDGcode);  /////Build the ALICE particle with Track XYZClusters ordered forward and closest MC pxyz
              //SetParticleLoop(particle_f);
              BuiltParticleUnsmearedNoRot(particle_f_nosmear,GArCenter,TrkClusterXYZf_NDGAr,trajxyzf_NDGAr);
              //if(particle_f.fParamMC.size()>50)
              //bool checklarm = CheckLArm(particle_f,75);
              //particle_f.reconstructParticleRotate0(geom,PDGcode,10000);
              //checkloop = CheckLooper(particle_f);
              //if (!checkloop){
                for(int a=0;a<=180;a++)
                {
                  BuildParticleNoRotation(particle_f,GArCenter,geom,TrkClusterXYZf_NDGAr,trajpxyzf_NDGAr,PDGcode,a*Pi()/180.);
                  bool statuscheck = particle_f.reconstructParticleRotate0(geom,PDGcode,10000);
                  BuiltParticleUnsmearedNoRot(particle_f_nosmear,GArCenter,TrkClusterXYZf_NDGAr,trajxyzf_NDGAr,a*Pi()/180.);
                  if(particle_f.fParamInRot[0].GetCovariance()[0]!=0) {
                    break;
                  }
                }
              //}
              checkloop = CheckLooperMC(particle_f);
              if(particle_f.fParamInRot[0].GetCovariance()[0]==0 || particle_f.fParamMC.size()>500){
                BuildParticle(particle_f,GArCenter,geom,TrkClusterXYZf_NDGAr,trajpxyzf_NDGAr,PDGcode,278,0);
                bool statuscheck = particle_f.reconstructParticleFull(geom,PDGcode,10000);
                BuildInRotfromIn(particle_f);
              }
            }
            else{
              std::cout<<"Backwards"<<std::endl;
              BuildParticleNoRotation(particle_f,GArCenter,geom,TrkClusterXYZb_NDGAr,trajpxyzb_NDGAr,PDGcode);  /////Build the ALICE particle with Track XYZClusters ordered forward and closest MC pxyz
              //SetParticleLoop(particle_f);
              BuiltParticleUnsmearedNoRot(particle_f_nosmear,GArCenter,TrkClusterXYZb_NDGAr,trajxyzb_NDGAr);
              //if(particle_f.fParamMC.size()>50) continue;
              // particle_f.reconstructParticleRotate0(geom,PDGcode,10000);
              //bool checklarm = CheckLArm(particle_f,75);
              //checkloop = CheckLooper(particle_f);
              //if (!checkloop){
                  for(int a=0;a<=180;a++)
                  {
                    BuildParticleNoRotation(particle_f,GArCenter,geom,TrkClusterXYZb_NDGAr,trajpxyzb_NDGAr,PDGcode,a*Pi()/180.);
                    bool statuscheck = particle_f.reconstructParticleRotate0(geom,PDGcode,10000);
                    BuiltParticleUnsmearedNoRot(particle_f_nosmear,GArCenter,TrkClusterXYZf_NDGAr,trajxyzf_NDGAr,a*Pi()/180.);
                    if(particle_f.fParamInRot[0].GetCovariance()[0]!=0) {
                      break;
                    }
                  }
              //}
              checkloop = CheckLooperMC(particle_f);
              if(particle_f.fParamInRot[0].GetCovariance()[0]==0 || particle_f.fParamMC.size()>500){
                BuildParticle(particle_f,GArCenter,geom,TrkClusterXYZb_NDGAr,trajpxyzb_NDGAr,PDGcode,278,0);
                bool statuscheck = particle_f.reconstructParticleFull(geom,PDGcode,10000);
                BuildInRotfromIn(particle_f);
              }
            }
            }

            ////// Create particle object with ALICE-like global coordinates
            

            if (dumpStream==kFALSE) continue;
            if (tree) tree->Fill();
            else {
              (*pcstream) << "fastPart" <<
                          "i=" << i <<
                          "t=" <<t<<
                          "Flip="<<Flip<<
                          "nPrimaries="<<nPrimaries<<
                          "PDGMotherCode="<<PDGMothercode<<
                          "xstart="<<xstart<<
                          "ystart="<<ystart<<
                          "VertGArx="<<VertGArx<<
                          "VertGAry="<<VertGAry<<
                          "VertGArz="<<VertGArz<<
                          "NuPx="<<NuPx<<
                          "NuPy="<<NuPy<<
                          "NuPz="<<NuPz<<
                          "Target="<<Target<<
                          "checkloop=" <<checkloop <<
                          "TrackLenEnd=" <<TrackLenEnd<<
                          "TrackLenSt=" <<TrackLenSt<<
                          "paramSt.="<<&paramSt<<
                          "paramEnd.="<<&paramEnd<<
                          "paramStMC.="<<&paramStMC<<
                          "paramEndMC.="<<&paramEndMC<<
                          "paramStNoRot.="<<&paramStNoRot<<
                          "paramEndNoRot.="<<&paramEndNoRot<<
                          "paramStNoRotMC.="<<&paramStNoRotMC<<
                          "paramEndNoRotMC.="<<&paramEndNoRotMC<<
                          "geom.="<<&geom<<
                          "part.=" << &particle_f <<
                          "partnosmear.=" << &particle_f_nosmear <<
                          "\n";
              tree=  ((*pcstream) << "fastPart").GetTree();
            }
          }          
      }
      
  }
  delete pcstream;
  timer.Print();
}


void initTreeFast(const char * inputList="fastParticle.list"){
  const char* inputListPath=gSystem->ExpandPathName(inputList);
  treeFast  = AliXRDPROOFtoolkit::MakeChainRandom(inputListPath,"fastPart",0,10000);
  treeTurn  = AliXRDPROOFtoolkit::MakeChainRandom(inputListPath,"turn",0,10000);
  treeUnit0  = AliXRDPROOFtoolkit::MakeChainRandom(inputListPath,"UnitTestDumpCorrectForMaterial",0,10000);
  treeSeed  = AliXRDPROOFtoolkit::MakeChainRandom(inputListPath,"seedDump",0,10000);
  treeFast->SetMarkerStyle(21);
  treeFast->SetMarkerSize(0.5);
   treeUnit0->SetMarkerStyle(21);
  treeUnit0->SetMarkerSize(0.5);
  treeSeed->BuildIndex("gid");
  treeFast->BuildIndex("gid");
  treeSeed->AddFriend(treeFast,"F");

  //AliDrawStyle::SetDefaults();
  //AliDrawStyle::ApplyStyle("figTemplate");
  gStyle->SetOptTitle(1);
  fastParticle::setAliases(*treeFast);
  //
  treeUnit0->SetAlias("dEdxOutIn","AliExternalTrackParam::BetheBlochSolid(0+paramStepRK.P()/mass)/AliExternalTrackParam::BetheBlochSolid(paramIn.P()/mass)");
  treeUnit0->SetAlias("dEdxIn","AliExternalTrackParam::BetheBlochSolid(paramIn.P()/mass+0)");
  treeUnit0->SetAlias("dEdxOut","AliExternalTrackParam::BetheBlochSolid(paramRK.P()/mass+0)");

  treeFast->SetAlias("ptMCSt","1./TMath::Abs(paramStMC.fP[4])");
  treeFast->SetAlias("rMCSt","TMath::Sqrt((1. - paramStMC.fP[2])*(1. + paramStMC.fP[2]))");
  treeFast->SetAlias("pxMCSt","ptMCSt*(rMCSt*TMath::Cos(paramStMC.fAlpha) - paramStMC.fP[2]*TMath::Sin(paramStMC.fAlpha))");
  treeFast->SetAlias("pyMCSt","ptMCSt*(paramStMC.fP[2]*TMath::Cos(paramStMC.fAlpha) + rMCSt*TMath::Sin(paramStMC.fAlpha))");
  treeFast->SetAlias("pzMCSt","ptMCSt*paramStMC.fP[3]");



  treeFast->SetAlias("ptInSt","1./TMath::Abs(part.fParamIn[0].fP[4])");
  treeFast->SetAlias("rInSt","TMath::Sqrt((1. - part.fParamIn[0].fP[2])*(1. + part.fParamIn[0].fP[2]))");
  treeFast->SetAlias("pxInSt","ptInSt*(rInSt*TMath::Cos(part.fParamIn[0].fAlpha) - part.fParamIn[0].fP[2]*TMath::Sin(part.fParamIn[0].fAlpha))");
  treeFast->SetAlias("pyInSt","ptInSt*(part.fParamIn[0].fP[2]*TMath::Cos(part.fParamIn[0].fAlpha) + rInSt*TMath::Sin(part.fParamIn[0].fAlpha))");
  treeFast->SetAlias("pzInSt","ptInSt*part.fParamIn[0].fP[3]");

  treeFast->SetAlias("ptGArSt","1./TMath::Abs(paramSt.fP[4])");
  treeFast->SetAlias("rGArSt","TMath::Sqrt((1. - paramSt.fP[2])*(1. + paramSt.fP[2]))");
  treeFast->SetAlias("pxGArSt","ptGArSt*(rGArSt*TMath::Cos(paramSt.fAlpha) - paramSt.fP[2]*TMath::Sin(paramSt.fAlpha))");
  treeFast->SetAlias("pyGArSt","ptGArSt*(paramSt.fP[2]*TMath::Cos(paramSt.fAlpha) + rGArSt*TMath::Sin(paramSt.fAlpha))");
  treeFast->SetAlias("pzGArSt","ptGArSt*paramSt.fP[3]");

  // treeFast->SetAlias("yend","part.fParamMC[fParamMC@.size()-1].fX*sin(part.fParamMC[fParamMC@.size()-1].fAlpha)+part.fParamMC[fParamMC@.size()-1].fP[0]*cos(part.fParamMC[fParamMC@.size()-1].fAlpha)");
  // treeFast->SetAlias("ystart","part.fParamMC[0].fX*sin(part.fParamMC[0].fAlpha)+part.fParamMC[0].fP[0]*cos(part.fParamMC[0].fAlpha)");
  // treeFast->SetAlias("xend","part.fParamMC[fParamMC@.size()-1].fX*cos(part.fParamMC[fParamMC@.size()-1].fAlpha)-part.fParamMC[fParamMC@.size()-1].fP[0]*sin(part.fParamMC[fParamMC@.size()-1].fAlpha)");
  // treeFast->SetAlias("xstart","part.fParamMC[0].fX*cos(part.fParamMC[0].fAlpha)-part.fParamMC[0].fP[0]*sin(part.fParamMC[0].fAlpha)");
  // treeFast->SetAlias("lArmMC","sqrt((xend-xstart)*(xend-xstart)+(yend-ystart)*(yend-ystart))");

  treeFast->SetAlias("gxMCdis","part.fParamMC[].fX*cos(part.fParamMC[].fAlpha)-part.fParamMC[].fP[0]*sin(part.fParamMC[].fAlpha)");
  treeFast->SetAlias("gyMCdis","part.fParamMC[].fX*sin(part.fParamMC[].fAlpha)+part.fParamMC[].fP[0]*cos(part.fParamMC[].fAlpha)");
  treeFast->SetAlias("gxMC","gxMCdis-xstart");
  treeFast->SetAlias("gyMC","gyMCdis-ystart");

  
  treeFast->SetAlias("gxIndis","part.fParamIn[].fX*cos(part.fParamIn[].fAlpha)-part.fParamIn[].fP[0]*sin(part.fParamIn[].fAlpha)");
  treeFast->SetAlias("gyIndis","part.fParamIn[].fX*sin(part.fParamIn[].fAlpha)+part.fParamIn[].fP[0]*cos(part.fParamIn[].fAlpha)");
  treeFast->SetAlias("gxIn","gxIndis-xstart");
  treeFast->SetAlias("gyIn","gyIndis-ystart");
  treeFast->SetAlias("rMC","sqrt(gyMC*gyMC+(gxMC-278)*(gxMC-278))");

  treeFast->SetAlias("gxMCSt","(part.fParamMC[0].fX*cos(part.fParamMC[0].fAlpha)-part.fParamMC[0].fP[0]*sin(part.fParamMC[0].fAlpha))-xstart");
  treeFast->SetAlias("gyMCSt","(part.fParamMC[0].fX*sin(part.fParamMC[0].fAlpha)+part.fParamMC[0].fP[0]*cos(part.fParamMC[0].fAlpha))-ystart");
  treeFast->SetAlias("gxMCEnd","(part.fParamMC[part.fParamMC@.size()-1].fX*cos(part.fParamMC[part.fParamMC@.size()-1].fAlpha)-part.fParamMC[part.fParamMC@.size()-1].fP[0]*sin(part.fParamMC[part.fParamMC@.size()-1].fAlpha))-xstart");
  treeFast->SetAlias("gyMCEnd","(part.fParamMC[part.fParamMC@.size()-1].fX*sin(part.fParamMC[part.fParamMC@.size()-1].fAlpha)+part.fParamMC[part.fParamMC@.size()-1].fP[0]*cos(part.fParamMC[part.fParamMC@.size()-1].fAlpha))-ystart");

  treeFast->SetAlias("gxMCgdis","paramStMC.fX*cos(paramStMC.fAlpha)-paramStMC.fP[0]*sin(paramStMC.fAlpha)");
  treeFast->SetAlias("gyMCgdis","paramStMC.fX*sin(paramStMC.fAlpha)+paramStMC.fP[0]*cos(paramStMC.fAlpha)");

  treeFast->SetAlias("gxMCg","gxMCgdis-xstart");
  treeFast->SetAlias("gyMCg","gyMCgdis-ystart");
  treeFast->SetAlias("gzMCg","paramStMC.fP[1]");
  treeFast->SetAlias("rMCg","sqrt(gygMC*gygMC+(gxgMC-278)*(gxgMC-278))");

  treeFast->SetAlias("gxMCno","partnosmear.fParamMC[].fX*cos(partnosmear.fParamMC[].fAlpha)-partnosmear.fParamMC[].fP[0]*sin(partnosmear.fParamMC[].fAlpha)");
  treeFast->SetAlias("gyMCno","partnosmear.fParamMC[].fX*sin(partnosmear.fParamMC[].fAlpha)+partnosmear.fParamMC[].fP[0]*cos(partnosmear.fParamMC[].fAlpha)");
  treeFast->SetAlias("rMCno","sqrt(gyMCno*gyMCno+(gxMCno-278)*(gxMCno-278))");

  treeFast->SetAlias("ResALICE","(fParamMC[0].GetP()-fParamInRot[0].GetP())/fParamMC[0].GetP()");
  treeFast->SetAlias("ResGAr","(fParamMC[0].GetP()-paramSt.GetP())/fParamMC[0].GetP()");

  treeFast->SetAlias("ResALICE","(fParamMC[0].GetP()-fParamInRot[0].GetP())/fParamMC[0].GetP()");
  treeFast->SetAlias("ResGAr","(fParamMC[0].GetP()-paramSt.GetP())/fParamMC[0].GetP()");

}
