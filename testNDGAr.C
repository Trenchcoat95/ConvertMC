/*
    gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/aliKalman/test/\"")
    gSystem->AddIncludePath("-I\"$ConvertMC/Local/\"")
    gSystem->Load("$fastMCKalman/fastMCKalman/aliKalman/test/AliExternalTrackParam.so");

   .L  Local/fastSimulation.cxx++g
    .L testNDGAr.C++g
    AliPDG::AddParticlesToPdgDataBase();
    testNDGAr(300,kTRUE)
    //initTreeFast()
    .> a.log
     testPCStream(5000,kTRUE);
     .>

 */
#include "fastSimulation.h"
#include "TTreeStream.h"
#include "TStopwatch.h"
#include "TRandom.h"
#include "TChain.h"
#include "AliXRDPROOFtoolkit.h"
//#include "AliDrawStyle.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TPad.h"
#include "TCanvas.h"
#include "AliPID.h"
#include <iostream>
#include <string>
#include "garutils.h"
const Float_t kDecayFraction=0.5;
const Float_t kRandomPDGFraction=0.5;

TChain * treeFast = 0;
TChain * treeTurn=0;
TChain * treeUnit0=0;
TChain * treeSeed=0;

void BuildParticle(fastParticle &particle, std::vector<TVector3> ClusterXYZ, double Center[3], fastGeometry geom)
{
          uint fMaxLayer = 0;
          for(size_t k=0;k<ClusterXYZ.size();k++) 
          {
              particle.fDirection.resize(k+1);
              particle.fParamMC.resize(k+1);
              particle.fLayerIndex.resize(k+1);

              TVector3 xyz_conv(ClusterXYZ.at(k).Z()-Center[2],
                                ClusterXYZ.at(k).Y()-Center[1],
                                ClusterXYZ.at(k).X()-Center[0]);

              Double_t alpha=TMath::ATan2(xyz_conv.Y(),xyz_conv.X());
              Double_t radius=sqrt(xyz_conv.Y()*xyz_conv.Y()+xyz_conv.X()*xyz_conv.X());
              Double_t X_loc =  xyz_conv.X()*cos(alpha) + xyz_conv.Y()*sin(alpha);
              Double_t Y_loc = -xyz_conv.X()*sin(alpha) + xyz_conv.Y()*cos(alpha);
              Double_t param[5]={Y_loc,xyz_conv.Z(),0,0,0};
              particle.fParamMC[k].SetParamOnly(X_loc,alpha,param);
              uint indexR = uint(std::upper_bound (geom.fLayerRadius.begin(),geom.fLayerRadius.end(), radius)-geom.fLayerRadius.begin());
              particle.fLayerIndex[k] = indexR;
              if(k!=0)
              {
                TVector3 xyz_prev(ClusterXYZ.at(k-1).Z()-Center[2],
                                  ClusterXYZ.at(k-1).Y()-Center[1],
                                  ClusterXYZ.at(k-1).X()-Center[0]);
                Double_t r_prev = sqrt(xyz_prev.Y()*xyz_prev.Y()+xyz_prev.X()*xyz_prev.X());
                if((radius/r_prev)>1) particle.fDirection[k] = +1;
                else particle.fDirection[k] = -1;
              }
              else particle.fDirection[k] = +1;

              if (indexR>fMaxLayer) fMaxLayer=indexR;
          }
          if(particle.fDirection.size()>1) particle.fDirection[0]=particle.fDirection[1];
}

void testNDGAr(Int_t nEv, bool dumpStream=1){

  const Int_t   nLayerTPC=277;
  const Int_t   nPoints=nLayerTPC*3;
  const Float_t xx0=7.8350968e-05;
  const Float_t xrho=0.0016265266;
  const Float_t kMaterialScaling=10;      ////Promote to global variable
  double GArCenter[3]={0,-150.473,1486}; 
  double GAr_r = 349.9;
  double GAr_L = 669.6;
  float fSortDistCut = 10.0;
  float fPrintLevel = 0.0;
            


  TStopwatch timer;
  timer.Start();
  

  //////////////////////////////////////////////////////////////////////////////////Create result File
  std::string filename = "FullReco";
  std::string path ="/home/federico/Documents/Universita/Federico_2020-2021/Aliwork/ConvertMC/Reco";
  std::string dirname = path+filename;
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
  tree_source->Add("/home/federico/Documents/Universita/Federico_2020-2021/Aliwork/ConvertMC/MC/*");
  std::vector<size_t>  *TPCClusterTrkIDNumber = 0;
  std::vector<float>   *TPCClusterX = 0;
  std::vector<float>   *TPCClusterY = 0;
  std::vector<float>   *TPCClusterZ = 0;
  TBranch        *b_TPCClusterTrkIDNumber;   //!
  TBranch        *b_TPCClusterX;   //!
  TBranch        *b_TPCClusterY;   //!
  TBranch        *b_TPCClusterZ;   //!
  tree_source->SetBranchAddress("TPCClusterTrkIDNumber", &TPCClusterTrkIDNumber, &b_TPCClusterTrkIDNumber);
  tree_source->SetBranchAddress("TPCClusterX", &TPCClusterX, &b_TPCClusterX);
  tree_source->SetBranchAddress("TPCClusterY", &TPCClusterY, &b_TPCClusterY);
  tree_source->SetBranchAddress("TPCClusterZ", &TPCClusterZ, &b_TPCClusterZ);
  int nEvSource = tree_source->GetEntries();
  nEv = std::min(nEv,nEvSource);
  
  




  /////////////////////////////////////////////////////////////////////////////////Build the detector geometry
  fastGeometry geom(nLayerTPC+1);
  geom.fBz=5;
  float resol[2]={0.0001,0.0001};
  resol[0]=0.1;                   
  resol[1]=0.1;
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
  for (Int_t i=0; i<nEv; i++){
 
      tree_source->GetEntry(i);
      if(TPCClusterTrkIDNumber->size()==0) continue;


      size_t ID = TPCClusterTrkIDNumber->at(0);
      std::vector<std::vector<TVector3>> TrksXYZ;
      std::vector<TVector3> TrkClusterXYZ;
      std::vector<size_t> ID_vector;
      ID_vector.push_back(ID);
      for(UInt_t j=0; j<TPCClusterTrkIDNumber->size();j++)
      {
          TVector3 xyz(TPCClusterX->at(j),TPCClusterY->at(j),TPCClusterZ->at(j));

          if(TPCClusterTrkIDNumber->at(j)==ID)
          {
            TrkClusterXYZ.push_back(xyz);
          }
          else
          {
            ID=TPCClusterTrkIDNumber->at(j);
            ID_vector.push_back(ID);
            TrksXYZ.push_back(TrkClusterXYZ);
            TrkClusterXYZ.clear();
            TrkClusterXYZ.push_back(xyz);
          }
      }
      ID_vector.push_back(ID);
      TrksXYZ.push_back(TrkClusterXYZ);

      for(size_t t=0; t<TrksXYZ.size(); t++)
      {
          std::cout<<"ID: "<<ID_vector.at(t)<<std::endl;

          ///////////////////////////////////////////// Sort TPC Clusters
          std::vector<int> hlf;
          std::vector<int> hlb;
          float lftmp=0;  
          float lbtmp=0;
          sort_TPCClusters_along_track2(TrksXYZ.at(t),hlf,hlb,fPrintLevel,lftmp,lbtmp,fSortDistCut);

          std::vector<TVector3> TrkClusterXYZb_NDGAr;
          std::vector<TVector3> TrkClusterXYZf_NDGAr;
          for(size_t k=0;k<hlb.size();k++) TrkClusterXYZb_NDGAr.push_back(TrksXYZ.at(t).at(hlb[k]));
          for(size_t k=0;k<hlf.size();k++) TrkClusterXYZf_NDGAr.push_back(TrksXYZ.at(t).at(hlf[k]));

          ////// Create particle object with ALICE-like global coordinates
          fastParticle particle(hlb.size()+1);
          particle.fAddMSsmearing=true;
          particle.fgStreamer=pcstream;
          particle.gid=ID_vector.at(t);
          particle.fDecayLength=0;
          BuildParticle(particle,TrkClusterXYZb_NDGAr,GArCenter,geom);

          particle.reconstructParticleFull(geom,211,particle.fLayerIndex[0]);

          if (dumpStream==kFALSE) continue;
          if (tree) tree->Fill();
          else {
            (*pcstream) << "fastPart" <<
                        "i=" << i <<
                        "geom.="<<&geom<<
                        "part.=" << &particle <<
                        "\n";
            tree=  ((*pcstream) << "fastPart").GetTree();
          }       
      }
      
  }
  delete pcstream;
  timer.Print();
}

