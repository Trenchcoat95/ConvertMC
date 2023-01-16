/*

    gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/aliKalman/test/\"")
    gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/MC/\"")
    gSystem->Load("$fastMCKalman/fastMCKalman/aliKalman/test/AliExternalTrackParam.so");
    .L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx++g
    .L $ConvertMC/testNDGAr.C++g
    AliPDG::AddParticlesToPdgDataBase();
    testNDGAr(300,kTRUE,kTRUE)

 */
#include "fastSimulation.h"
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
#include <iostream>
#include <string>
#include "garutils.h"
const Float_t kDecayFraction=0.5;
const Float_t kRandomPDGFraction=0.5;

TChain * treeFast = 0;
TChain * treeTurn=0;
TChain * treeUnit0=0;
TChain * treeSeed=0;

size_t ClosestPoint(std::vector<TVector3> trajxyz, TVector3 ClusterXYZ)
{
  Double_t dist = 10000;
  size_t index=0;
  for(size_t i=0; i<trajxyz.size(); i++)
  {
    double cl[] = {ClusterXYZ.X(),ClusterXYZ.Y(),ClusterXYZ.Z()};
    double traj[] = {trajxyz.at(i).X(),trajxyz.at(i).Y(),trajxyz.at(i).Z()};
    TVector3 diff = trajxyz.at(i)-ClusterXYZ;
    Double_t checkdist = diff.Mag2();
    if(checkdist<dist)
    {
      dist=checkdist;
      index=i;
    }

  }
  return index;
}

void ClosestTrajectory(std::vector<TVector3> trajxyz, std::vector<TVector3> trajpxyz, std::vector<TVector3> ClusterXYZ, 
                        std::vector<TVector3> &trajpxyz_closest)
{
  for(size_t t=0;t<ClusterXYZ.size();t++)
  {
    trajpxyz_closest.push_back(trajpxyz[ClosestPoint(trajxyz,ClusterXYZ[t])]);
  }
}

Int_t BuildParticle(fastParticle &particle, double Center[3], fastGeometry geom, 
                   std::vector<TVector3> trajxyz, std::vector<TVector3> trajpxyz, long PDGcode)
{
          uint fMaxLayer = 0;
          TParticlePDG *p = TDatabasePDG::Instance()->GetParticle(PDGcode);
          if (p == nullptr) {
            ::Error("fastParticle::simulateParticle", "Invalid pdgCode %ld", PDGcode);
            return -1;
          }
          Short_t sign = 1 * p->Charge() / 3.;
          Float_t mass = p->Mass();
          particle.fMassMC=mass;

          for(size_t k=0;k<trajxyz.size();k++) 
          {
              Bool_t invert = kFALSE;
              Double_t xyz_conv[3]= {(trajxyz.at(k).Z()-Center[2]),
                                trajxyz.at(k).Y()-Center[1],
                                trajxyz.at(k).X()-Center[0]};

              Double_t pxyz_conv[3]= {(trajpxyz.at(k).Z()),
                                      trajpxyz.at(k).Y(),
                                      trajpxyz.at(k).X()};

              Double_t alpha=TMath::ATan2(xyz_conv[1],xyz_conv[0]);
              Double_t radius=sqrt(xyz_conv[1]*xyz_conv[1]+xyz_conv[0]*xyz_conv[0]);
              Double_t X_loc =  xyz_conv[0]*cos(alpha) + xyz_conv[1]*sin(alpha);
              Double_t Y_loc = -xyz_conv[0]*sin(alpha) + xyz_conv[1]*cos(alpha);
              //Double_t param[5]={Y_loc,xyz_conv.Z(),0,0,0};              
              //particle.fParamMC[k].SetParamOnly(X_loc,alpha,param);
              double covar[21]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
              AliExternalTrackParam param(xyz_conv,pxyz_conv,covar,sign);
              AliExternalTrackParam4D param4D(param,mass,1);
              Double_t pxyz_conv_test[3];
              param4D.GetPxPyPz(pxyz_conv_test);
              Bool_t status = param4D.Rotate(alpha);
              Double_t pxyz_inv[3] = {-pxyz_conv[0],-pxyz_conv[1],-pxyz_conv[2]};

              if(abs(pxyz_conv_test[0]-pxyz_conv[0])>0.00001)
              {
                AliExternalTrackParam param_inv(xyz_conv,pxyz_inv,covar,sign);
                AliExternalTrackParam4D param4D_inv(param_inv,mass,1);
                status = param4D_inv.Rotate(alpha);
                // Double_t pxyz_inv_test[3];
                // param4D_inv.GetPxPyPz(pxyz_inv_test);
                // param4D = param4D_inv;
                // Double_t pxyz_inv_test2[3];
                // param4D.GetPxPyPz(pxyz_inv_test2);
                invert=kTRUE;
              }

              particle.fDirection.resize(k+1);
              particle.fParamMC.resize(k+1);
              particle.fLayerIndex.resize(k+1);
              if(status)
              {
                if(!invert)
                {
                  particle.fParamMC[k].Set(xyz_conv,pxyz_conv,covar,sign);
                  particle.fParamMC[k].Rotate(alpha);
                }
                else
                {
                  particle.fParamMC[k].Set(xyz_conv,pxyz_inv,covar,-sign);
                  particle.fParamMC[k].Rotate(alpha);
                }
              }
              else
              {
                double ptemp[] = {Y_loc,xyz_conv[2],0,0,0};
                particle.fParamMC[k].SetParamOnly(X_loc,alpha,ptemp);
              }
              //particle.fParamMC[k].Rotate(alpha);
              uint indexR = uint(std::upper_bound (geom.fLayerRadius.begin(),geom.fLayerRadius.end(), radius)-geom.fLayerRadius.begin());
              particle.fLayerIndex[k] = indexR;
              if(k!=0)
              {
                TVector3 xyz_prev(trajxyz.at(k-1).Z()-Center[2],
                                  trajxyz.at(k-1).Y()-Center[1],
                                  trajxyz.at(k-1).X()-Center[0]);
                Double_t r_prev = sqrt(xyz_prev.Y()*xyz_prev.Y()+xyz_prev.X()*xyz_prev.X());
                if((radius/r_prev)>1) particle.fDirection[k] = +1;
                else particle.fDirection[k] = -1;
              }
              else particle.fDirection[k] = +1;

              if (indexR>fMaxLayer) fMaxLayer=indexR;
          }
          if(particle.fDirection.size()>1) particle.fDirection[0]=particle.fDirection[1];

          return 1;
}


void CombineParticle(fastParticle &particle, fastParticle part1, fastParticle part2)
{

          particle.fParamMC.resize(part1.fParamMC.size());
          particle.fParamMC = part1.fParamMC;

          particle.fParamIn.resize(part1.fParamRefit.size());
          particle.fParamIn = part1.fParamRefit;
          particle.fStatusMaskIn.resize(part1.fStatusMaskRefit.size());
          particle.fStatusMaskIn = part1.fStatusMaskRefit;

          particle.fParamOut.resize(part2.fParamRefit.size());
          for(size_t i=0; i<part2.fParamRefit.size(); i++) particle.fParamIn[i]=part2.fParamRefit[part2.fParamRefit.size()-1-i];
          particle.fStatusMaskOut.resize(part2.fStatusMaskRefit.size());
          for(size_t i=0; i<part2.fStatusMaskRefit.size(); i++) particle.fStatusMaskOut[i]=part2.fStatusMaskRefit[part2.fStatusMaskRefit.size()-1-i];
}

void testNDGAr(Int_t nEv, bool dumpStream=1, bool Ideal=kTRUE){

  const Int_t   nLayerTPC=278;
  const Int_t   nPoints=nLayerTPC*3;
  const Float_t xx0=8.37758e-04; //1/X0 cm^-1 for ArCH4 at 10 atm
  const Float_t xrho=0.016770000; //rho g/cm^3 for ArCH4 at 10 atm
  const Float_t kMaterialScaling=1;      ////Promote to global variable
  double GArCenter[3]={0,-150.473,1486}; 
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
  std::vector<int>  *TPCClusterMCindex = 0;
  std::vector<size_t>  *TrackIDNumber = 0;
  std::vector<int>  *TrackMCindex = 0;
  std::vector<int>  *TrajMCPIndex = 0;
  std::vector<int>  *MCTrkID = 0;
  std::vector<int>     *PDG = 0;
  std::vector<float>   *TPCClusterX = 0;
  std::vector<float>   *TPCClusterY = 0;
  std::vector<float>   *TPCClusterZ = 0;
  std::vector<float>   *TrajMCPPX = 0;
  std::vector<float>   *TrajMCPPY = 0;
  std::vector<float>   *TrajMCPPZ = 0;
  std::vector<float>   *TrajMCPX = 0;
  std::vector<float>   *TrajMCPY = 0;
  std::vector<float>   *TrajMCPZ = 0;
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
  TBranch        *b_TPCClusterTrkIDNumber;   //!
  TBranch        *b_TPCClusterMCindex;
  TBranch        *b_TrackIDNumber;
  TBranch        *b_TrackMCindex;
  TBranch        *b_TrajMCPIndex;
  TBranch        *b_MCTrkID;
  TBranch        *b_PDG;
  TBranch        *b_TPCClusterX;   //!
  TBranch        *b_TPCClusterY;   //!
  TBranch        *b_TPCClusterZ;   //!
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
  TBranch   *b_TrackStartQ     ;
  TBranch *b_TrackEndX       ;
  TBranch *b_TrackEndY       ;
  TBranch *b_TrackEndZ       ;
  TBranch *b_TrackEndPX      ;
  TBranch *b_TrackEndPY      ;
  TBranch *b_TrackEndPZ      ;
  TBranch   *b_TrackEndQ       ;
  tree_source->SetBranchAddress("TPCClusterTrkIDNumber", &TPCClusterTrkIDNumber, &b_TPCClusterTrkIDNumber);
  tree_source->SetBranchAddress("TPCClusterMCindex", &TPCClusterMCindex, &b_TPCClusterMCindex);
  tree_source->SetBranchAddress("TrackIDNumber", &TrackIDNumber, &b_TrackIDNumber);
  tree_source->SetBranchAddress("TrackMCindex", &TrackMCindex, &b_TrackMCindex);
  tree_source->SetBranchAddress("TrajMCPIndex", &TrajMCPIndex, &b_TrajMCPIndex);
  tree_source->SetBranchAddress("MCTrkID", &MCTrkID, &b_MCTrkID);
  tree_source->SetBranchAddress("PDG", &PDG, &b_PDG);
  tree_source->SetBranchAddress("TPCClusterX", &TPCClusterX, &b_TPCClusterX);
  tree_source->SetBranchAddress("TPCClusterY", &TPCClusterY, &b_TPCClusterY);
  tree_source->SetBranchAddress("TPCClusterZ", &TPCClusterZ, &b_TPCClusterZ);
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

  int nEvSource = tree_source->GetEntries();
  nEv = std::min(nEv,nEvSource);
  
  




  /////////////////////////////////////////////////////////////////////////////////Build the detector geometry
  fastGeometry geom(nLayerTPC+1);
  geom.fBz=-5;
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
      if(TrackIDNumber->size()==0) continue;




      size_t ID = 0;
      size_t MCID = 0;
      std::vector<size_t> ID_vector;                    ///////For Each event find all the reconstructed tracks and save their IDs
      std::vector<int> MCID_vector;                     ///////For each track find the corresponding MC index to match information with MC truth
      std::vector<std::vector<TVector3>> TrksXYZ;       ///////Vector containing one vector of TPCCluster XYZ points for each track
      std::vector<TVector3> TrkClusterXYZ;              /////Container that will be filled for each track, added to TrksXYZ and then deleted
                                                        /////Note: for now all is in ND-GAr coordinates
      for(UInt_t k=0; k<TrackIDNumber->size();k++)
      {
          ID = TrackIDNumber->at(k);
          MCID = TrackMCindex->at(k);
          for(UInt_t j=0; j<TPCClusterTrkIDNumber->size();j++)
          {
            TVector3 xyz(TPCClusterX->at(j),TPCClusterY->at(j),TPCClusterZ->at(j));
            if(ID==TPCClusterTrkIDNumber->at(j)) TrkClusterXYZ.push_back(xyz);
          }
          ID_vector.push_back(ID);
          MCID_vector.push_back(MCID);
          TrksXYZ.push_back(TrkClusterXYZ);
          TrkClusterXYZ.clear();
      }


      
      for(size_t t=0; t<TrksXYZ.size(); t++)   ////Cycle over all the tracks for the event 
                                               ////Note that in the final tree the separation between events is lost and each item is a track)
      {    
          Double_t TStX = TrackStartZ->at(t)- GArCenter[2];         //////Save track's ND-GAr reconstructed quantities in the ALICE coordinate frame
          Double_t TStY = TrackStartY->at(t)- GArCenter[1];
          Double_t TStZ =  TrackStartX->at(t)- GArCenter[0];
          Double_t TStPX = TrackStartPZ->at(t);
          Double_t TStPY = TrackStartPY->at(t);
          Double_t TStPZ =  TrackStartPX->at(t);
          int TStQ = TrackStartQ->at(t);
          Double_t TEndX = TrackEndZ->at(t)- GArCenter[2];
          Double_t TEndY = TrackEndY->at(t)- GArCenter[1];
          Double_t TEndZ =  TrackEndX->at(t)- GArCenter[0];
          Double_t TEndPX = TrackEndPZ->at(t);
          Double_t TEndPY = TrackEndPY->at(t);
          Double_t TEndPZ =  TrackEndPX->at(t);
          int TEndQ = TrackEndQ->at(t);

          std::vector<TVector3> trajxyz;                ///////Find MC trajectory that produced the track and save  the xyz's and pxyz's
          std::vector<TVector3> trajpxyz;
          for(size_t u=0; u<TrajMCPX->size(); u++)          
          {
            if(MCID_vector.at(t)==TrajMCPIndex->at(u))
            {
              TVector3 jxyz(TrajMCPX->at(u),TrajMCPY->at(u),TrajMCPZ->at(u));
              TVector3 pjxyz(TrajMCPPX->at(u),TrajMCPPY->at(u),TrajMCPPZ->at(u));
              double Cluster[] = {TrksXYZ.at(t).at(0).X(),TrksXYZ.at(t).at(0).Y(),TrksXYZ.at(t).at(0).Z()};
              trajxyz.push_back(jxyz);
              trajpxyz.push_back(pjxyz);
            }
          }

          if(trajxyz.size()==0) continue;       /////Skip tracks that were not successfully associated to a MC trajectory by the back-tracker

          
          long PDGcode=PDG->at(MCID_vector.at(t));    ////Save PDGCode for the MC trajectory associated with the ND-GAr reconstructed track

          if(PDGcode==0 || abs(PDGcode)==14 || abs(PDGcode)!=13) continue;   ////For now only checking muons
          std::cout<<"ID: "<<ID_vector.at(t)<<" PDG: "<<PDGcode<<std::endl;


          ///////////////////////////////////////////// Sort TPC Clusters using the ND-GAr method
          std::vector<int> hlf;
          std::vector<int> hlb;
          float lftmp=0;  
          float lbtmp=0;
          sort_TPCClusters_along_track2(TrksXYZ.at(t),hlf,hlb,fPrintLevel,lftmp,lbtmp,fSortDistCut);

          std::vector<TVector3> TrkClusterXYZb_NDGAr;
          std::vector<TVector3> TrkClusterXYZf_NDGAr;
          for(size_t k=0;k<hlb.size();k++) TrkClusterXYZb_NDGAr.push_back(TrksXYZ.at(t).at(hlb[k]));
          for(size_t k=0;k<hlf.size();k++) TrkClusterXYZf_NDGAr.push_back(TrksXYZ.at(t).at(hlf[k]));

          ///Now that we have vectors of sorted TPCClusters we can find the closest trajectory point for each cluster 
          ///and save the MC pxyz information for each point
          std::vector<TVector3> trajpxyzb_NDGAr;
          ClosestTrajectory(trajxyz,trajpxyz,TrkClusterXYZb_NDGAr,trajpxyzb_NDGAr);
          std::vector<TVector3> trajpxyzf_NDGAr;
          ClosestTrajectory(trajxyz,trajpxyz,TrkClusterXYZf_NDGAr,trajpxyzf_NDGAr);

          ////// Create particle object with ALICE-like global coordinates
          fastParticle particle(hlb.size()+1);
          particle.fAddMSsmearing=true;
          particle.fAddPadsmearing=false;
          particle.fUseMCInfo=false;
          particle.fgStreamer=pcstream;
          particle.gid=ID_vector.at(t);
          particle.fDecayLength=0;
        
          if(Ideal) {
            BuildParticle(particle,GArCenter,geom,trajxyz,trajpxyz,PDGcode);  //////Built the ALICE particle with the MC trajectory xyz and pxyz points 
            if(particle.fParamMC.size()==0) continue;
            particle.reconstructParticleFull(geom,PDGcode,10000);
            particle.reconstructParticleFullOut(geom,PDGcode,10000);
            particle.refitParticle();
          }
          else
          {
            BuildParticle(particle,GArCenter,geom,trajxyz,trajpxyz,PDGcode);  /////Build the ALICE particle with Track XYZClusters and closest MC pxyz
            if(particle.fParamMC.size()==0) continue;
            particle.reconstructParticleFull(geom,PDGcode,10000);
            particle.reconstructParticleFullOut(geom,PDGcode,10000);
            particle.refitParticle();
          }


          /*
          fastParticle particle2(hlf.size()+1);
          particle2.fAddMSsmearing=true;
          particle2.fAddPadsmearing=false;
          particle2.fUseMCInfo=false;
          particle2.fgStreamer=pcstream;
          particle2.gid=ID_vector.at(t);
          particle2.fDecayLength=0;
          BuildParticle(particle2,TrkClusterXYZf_NDGAr,GArCenter,geom);

          particle2.reconstructParticleFull(geom,211,particle.fLayerIndex[0]);
          particle2.reconstructParticleFullOut(geom,211,particle.fLayerIndex[0]);
          particle2.refitParticle();

          fastParticle particletot(hlf.size()+1);
          particletot.fAddMSsmearing=true;
          particletot.fAddPadsmearing=false;
          particletot.fUseMCInfo=false;
          particletot.fgStreamer=pcstream;
          particletot.gid=ID_vector.at(t);
          particletot.fDecayLength=0;
          CombineParticle(particletot,particle,particle2);
          particletot.refitParticle();
          */

          if (dumpStream==kFALSE) continue;
          if (tree) tree->Fill();
          else {
            (*pcstream) << "fastPart" <<
                        "i=" << i <<
                        "TStQ=" << TStQ <<
                        "TStX=" << TStX <<
                        "TStY=" << TStY <<
                        "TStZ=" << TStZ <<
                        "TStPX=" << TStPX <<
                        "TStPY=" << TStPY <<
                        "TStPZ=" << TStPZ <<
                        "TEndQ=" << TEndQ <<
                        "TEndX=" << TEndX <<
                        "TEndY=" << TEndY <<
                        "TEndZ=" << TEndZ <<
                        "TEndPX=" << TEndPX <<
                        "TEndPY=" << TEndPY <<
                        "TEndPZ=" << TEndPZ <<
                        "geom.="<<&geom<<
                        "part.=" << &particle <<
                        //"part2.="<< &particle2 <<
                        //"parttot.="<< &particle2 <<
                        "\n";
            tree=  ((*pcstream) << "fastPart").GetTree();
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

  treeFast->SetAlias("ptMCSt","1./TMath::Abs(part.fParamMC[0].fP[4])");
  treeFast->SetAlias("rMCSt","TMath::Sqrt((1. - part.fParamMC[0].fP[2])*(1. + part.fParamMC[0].fP[2]))");
  treeFast->SetAlias("pxMCSt","ptMCSt*(rMCSt*TMath::Cos(part.fParamMC[0].fAlpha) - part.fParamMC[0].fP[2]*TMath::Sin(part.fParamMC[0].fAlpha))");
  treeFast->SetAlias("pyMCSt","ptMCSt*(part.fParamMC[0].fP[2]*TMath::Cos(part.fParamMC[0].fAlpha) + rMCSt*TMath::Sin(part.fParamMC[0].fAlpha))");

}
