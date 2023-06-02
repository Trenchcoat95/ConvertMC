/*

    gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/aliKalman/test/\"")
    gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/MC/\"")
    gSystem->Load("$fastMCKalman/fastMCKalman/aliKalman/test/AliExternalTrackParam.so");
    .L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx++g
    .L $ConvertMC/testNDGAr.C++g
    AliPDG::AddParticlesToPdgDataBase();
    testNDGAr(10000,kTRUE,kTRUE,0,kTRUE,kFALSE)
    testNDGAr(300,kTRUE,kFALSE,0,kTRUE,kTRUE)

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
#include <iostream>
#include <string>
#include "garutils.h"
const Float_t kDecayFraction=0.5;
const Float_t kRandomPDGFraction=0.5;

TChain * treeFast = 0;
TChain * treeTurn=0;
TChain * treeUnit0=0;
TChain * treeSeed=0;

void BuiltParticleUnsmeared(fastParticle &particle, double Center[3], std::vector<TVector3> trajxyz, std::vector<TVector3> trajxyz_uns)
{
  for(size_t k=0;k<trajxyz.size();k++) 
  {
    Double_t xyz_conv[3]= {(trajxyz.at(k).Z()-Center[2]),
                              trajxyz.at(k).Y()-Center[1],
                              trajxyz.at(k).X()-Center[0]};
    Double_t xyz_conv_uns[3]= {(trajxyz_uns.at(k).Z()-Center[2]),
                              trajxyz_uns.at(k).Y()-Center[1],
                              trajxyz_uns.at(k).X()-Center[0]};
    Double_t alpha=TMath::ATan2(xyz_conv[1],xyz_conv[0]);
    Double_t X_loc =  xyz_conv_uns[0]*cos(alpha) + xyz_conv_uns[1]*sin(alpha);
    Double_t Y_loc = -xyz_conv_uns[0]*sin(alpha) + xyz_conv_uns[1]*cos(alpha);
    Double_t Z_loc =  xyz_conv_uns[2];
    double ptemp[] = {Y_loc,Z_loc,0,0,0};
    particle.fParamMC.resize(k+1);
    particle.fParamMC[k].SetParamOnly(X_loc,alpha,ptemp);  
  }
}

void BuiltParticleUnsmearedNoRot(fastParticle &particle, double Center[3], std::vector<TVector3> trajxyz, std::vector<TVector3> trajxyz_uns, double alpha=0)
{
  for(size_t k=0;k<trajxyz.size();k++) 
  {
    Double_t xyz_conv[3]= {(trajxyz.at(k).Z()-Center[2]),
                              trajxyz.at(k).Y()-Center[1],
                              trajxyz.at(k).X()-Center[0]};
    Double_t xyz_conv_uns[3]= {(trajxyz_uns.at(k).Z()-Center[2]),
                              trajxyz_uns.at(k).Y()-Center[1],
                              trajxyz_uns.at(k).X()-Center[0]};
    //Double_t alpha=TMath::ATan2(xyz_conv[1],xyz_conv[0]);
    Double_t X_loc =  xyz_conv_uns[0]*cos(alpha) + xyz_conv_uns[1]*sin(alpha);
    Double_t Y_loc = -xyz_conv_uns[0]*sin(alpha) + xyz_conv_uns[1]*cos(alpha);
    Double_t Z_loc =  xyz_conv_uns[2];
    double ptemp[] = {Y_loc,Z_loc,0,0,0};
    particle.fParamMC.resize(k+1);
    particle.fParamMC[k].SetParamOnly(X_loc,alpha,ptemp);  
  }
}

size_t ClosestPoint(std::vector<TVector3> trajxyz, TVector3 ClusterXYZ)
{
  Double_t dist = 10000;
  size_t index=0;
  for(size_t i=0; i<trajxyz.size(); i++)
  {
    double cl[] = {ClusterXYZ.X(),ClusterXYZ.Y(),ClusterXYZ.Z()};
    double traj[] = {trajxyz.at(i).X(),trajxyz.at(i).Y(),trajxyz.at(i).Z()};
    TVector3 diff = trajxyz.at(i)-ClusterXYZ;
    Double_t checkdist = sqrt(diff.Mag2());
    //checkdist = sqrt(pow(ClusterXYZ.X()-trajxyz.at(i).X(),2)+pow(ClusterXYZ.Y()-trajxyz.at(i).Y(),2));
    if(checkdist<dist)
    {
      dist=checkdist;
      index=i;
    }

  }
  return index;
}

void ClosestTrajectory(std::vector<TVector3> trajxyz, std::vector<TVector3> trajpxyz, std::vector<TVector3> ClusterXYZ, 
                        std::vector<TVector3> &trajxyz_closest, std::vector<TVector3> &trajpxyz_closest,Bool_t &IsForward)
{
  size_t FirstTrajIndex,LastTrajIndex;
  for(size_t t=0;t<ClusterXYZ.size();t++)
  {
    size_t cl_index= ClosestPoint(trajxyz,ClusterXYZ[t]);
    trajpxyz_closest.push_back(trajpxyz[cl_index]);
    trajxyz_closest.push_back(trajxyz[cl_index]);
    if(t==0) FirstTrajIndex=ClosestPoint(trajxyz,ClusterXYZ[t]);
    if(t==(ClusterXYZ.size()-1)) LastTrajIndex=ClosestPoint(trajxyz,ClusterXYZ[t]);
  }
  if(LastTrajIndex>FirstTrajIndex) IsForward=kTRUE;
  else IsForward=kFALSE;
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


Int_t BuildParticleNoRotation(fastParticle &particle, double Center[3], fastGeometry geom, 
                   std::vector<TVector3> trajxyz, std::vector<TVector3> trajpxyz, long PDGcode, float angle = 0.)
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
              Double_t xyz_conv[3]= {(trajxyz.at(k).Z()-(Center[2]-278)),
                                trajxyz.at(k).Y()-Center[1],
                                trajxyz.at(k).X()-Center[0]};

              Double_t pxyz_conv[3]= {(trajpxyz.at(k).Z()),
                                      trajpxyz.at(k).Y(),
                                      trajpxyz.at(k).X()};
              //Double_t totalp = sqrt(pxyz_conv[0]*pxyz_conv[0]+pxyz_conv[0]*pxyz_conv[0]+pxyz_conv[2]*pxyz_conv[2]);
              //if(totalp<0.01) break;

              Double_t alpha=TMath::ATan2(xyz_conv[1],xyz_conv[0]);
              Double_t radius=sqrt((xyz_conv[1])*(xyz_conv[1])+(xyz_conv[0]-278)*(xyz_conv[0]-278));
              Double_t X_loc =  xyz_conv[0]*cos(angle) + xyz_conv[1]*sin(angle);
              Double_t Y_loc = -xyz_conv[0]*sin(angle) + xyz_conv[1]*cos(angle);
              double covar[21]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

              particle.fDirection.resize(k+1);
              particle.fParamMC.resize(k+1);
              particle.fLayerIndex.resize(k+1);
              particle.fLoop.resize(k+1);

              particle.fLoop[k]=0;

              AliExternalTrackParam param(xyz_conv,pxyz_conv,covar,sign);
              bool statuscheck = param.Rotate(angle);
              if(statuscheck){ 
                particle.fParamMC[k].Set(xyz_conv,pxyz_conv,covar,sign);
                particle.fParamMC[k].Rotate(angle);
                bool x = 0;
              }
              else {
                double ptemp[] = {Y_loc,xyz_conv[2],0,0,0};
                particle.fParamMC[k].SetParamOnly(X_loc,angle,ptemp);
                double x = 0;
              }
              //bool statuscheck = particle.fParamMC[k].Rotate(0);
              
              //particle.fParamMC[k].Rotate(alphca);
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

Int_t SetParticleLoop(fastParticle &particle){

 Int_t sign = particle.fParamMC[0].GetParameter()[4]>0? 1:-1;
 Int_t loopcounter = 0;
  for(size_t k=0;k<particle.fParamMC.size();k++) {
     Int_t signtest = particle.fParamMC[k].GetParameter()[4];
     if(signtest>0) signtest=1;
     else if(signtest<0) signtest=-1;
     else if(signtest==0) signtest=0;

     signtest*=sign;
     if(signtest>0){
      particle.fLoop.resize(k+1);
      particle.fLoop[k]=loopcounter;
     }
     if(signtest<0){
      sign*=-1;
      loopcounter++;
      particle.fLoop.resize(k+1);
      particle.fLoop[k]=loopcounter;
     }
     if(signtest==0){
      int zerocounter=0;
      float checkparam=0;
      while(checkparam==0){
        zerocounter++;
        if((k+zerocounter)==particle.fParamMC.size()) return 0;
        checkparam = particle.fParamMC[k+zerocounter].GetParameter()[4];
      }
      for(int i=0;i<=zerocounter;i++){
        if(i<=int(zerocounter/2)){
          particle.fLoop.resize(k+i+1);
          particle.fLoop[k+i]=loopcounter+1;
        }else{
          particle.fLoop.resize(k+i+1);
          particle.fLoop[k+i]=loopcounter+1;
        }
      }
      loopcounter++;
      sign*=-1;
      k+=(zerocounter);
     }

  }
  return 1;
}


Int_t BuildParamNDGArReco(AliExternalTrackParam4D &paramNDGAr4D, double Center[3], 
                          TVector3 trajxyz, TVector3 trajpxyz, long PDGcode,
                          Short_t signr=0)
{
        TParticlePDG *p = TDatabasePDG::Instance()->GetParticle(PDGcode);
        if (p == nullptr) {
          ::Error("fastParticle::simulateParticle", "Invalid pdgCode %ld", PDGcode);
          return -1;
        }
        Float_t mass = p->Mass();
        Short_t sign;
        if(signr==0) sign = 1 * p->Charge() / 3.;
        else sign=signr;

        Bool_t invert = kFALSE;
        Double_t xyz_conv[3]= {(trajxyz.Z()-Center[2]),
                          trajxyz.Y()-Center[1],
                          trajxyz.X()-Center[0]};

        Double_t pxyz_conv[3]= {(trajpxyz.Z()),
                                  trajpxyz.Y(),
                                  trajpxyz.X()};

        Double_t alpha=TMath::ATan2(xyz_conv[1],xyz_conv[0]);
        Double_t radius=sqrt(xyz_conv[1]*xyz_conv[1]+xyz_conv[0]*xyz_conv[0]);
        Double_t X_loc =  xyz_conv[0]*cos(alpha) + xyz_conv[1]*sin(alpha);
        Double_t Y_loc = -xyz_conv[0]*sin(alpha) + xyz_conv[1]*cos(alpha);

        double covar[21]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        AliExternalTrackParam param(xyz_conv,pxyz_conv,covar,sign);
        AliExternalTrackParam4D param4D(param,mass,1);
        Double_t pxyz_conv_test[3];
        param4D.GetPxPyPz(pxyz_conv_test);
        Bool_t status = param4D.Rotate(alpha);
        Double_t pxyz_inv[3] = {-pxyz_conv[0],-pxyz_conv[1],-pxyz_conv[2]};

        if((abs(pxyz_conv_test[0]-pxyz_conv[0])>0.00001) || (!status))
        {
          AliExternalTrackParam param_inv(xyz_conv,pxyz_inv,covar,sign);
          AliExternalTrackParam4D param4D_inv(param_inv,mass,1);
          status = param4D_inv.Rotate(alpha);
          invert=kTRUE;
        }

        if(status)
        {
          if(!invert)
          {
            paramNDGAr4D.Set(xyz_conv,pxyz_conv,covar,sign);
            paramNDGAr4D.Rotate(alpha);
          }
          else
          {
            paramNDGAr4D.Set(xyz_conv,pxyz_inv,covar,-sign);
            paramNDGAr4D.Rotate(alpha);
          }
        }
        else
        {
          double ptemp[] = {Y_loc,xyz_conv[2],0,0,0};
          paramNDGAr4D.SetParamOnly(X_loc,alpha,ptemp);
        }


    return 1;
}


Int_t BuildParamNDGArRecoNoRotation(AliExternalTrackParam4D &paramNDGAr4D, double Center[3], 
                          TVector3 trajxyz, TVector3 trajpxyz, long PDGcode,
                          Short_t signr=0)
{
    TParticlePDG *p = TDatabasePDG::Instance()->GetParticle(PDGcode);
    if (p == nullptr) {
      ::Error("fastParticle::simulateParticle", "Invalid pdgCode %ld", PDGcode);
      return -1;
    }
    Float_t mass = p->Mass();
    Short_t sign;
    if(signr==0) sign = 1 * p->Charge() / 3.;
    else sign=signr;

    Bool_t invert = kFALSE;
    Double_t xyz_conv[3]= {(trajxyz.Z()-(Center[2]-278)),
                      trajxyz.Y()-Center[1],
                      trajxyz.X()-Center[0]};

    Double_t pxyz_conv[3]= {(trajpxyz.Z()),
                              trajpxyz.Y(),
                              trajpxyz.X()};
    Double_t pxyz_inv[3] = {-pxyz_conv[0],-pxyz_conv[1],-pxyz_conv[2]};

    //Double_t alpha=TMath::ATan2(xyz_conv[1],xyz_conv[0]);
    Double_t radius=sqrt(xyz_conv[1]*xyz_conv[1]+(xyz_conv[0]-278)*(xyz_conv[0]-278));
    //Double_t X_loc =  xyz_conv[0]*cos(alpha) + xyz_conv[1]*sin(alpha);
    //Double_t Y_loc = -xyz_conv[0]*sin(alpha) + xyz_conv[1]*cos(alpha);

    double covar[21]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    AliExternalTrackParam param(xyz_conv,pxyz_conv,covar,sign);
    AliExternalTrackParam4D param4D(param,mass,1);
    bool status = param4D.Rotate(0);

    Double_t pxyz_conv_test[3];
    param4D.GetPxPyPz(pxyz_conv_test);
    AliExternalTrackParam4D param4D_inv_perm;

    if((abs(pxyz_conv_test[0]-pxyz_conv[0])>0.00001) || (!status))
    {
          AliExternalTrackParam param_inv(xyz_conv,pxyz_inv,covar,sign);
          AliExternalTrackParam4D param4D_inv(param_inv,mass,1);
          status = param4D_inv.Rotate(0);
          param4D_inv_perm=param4D_inv;
          Double_t pxyz_conv_test_inv[3];
          param4D_inv.GetPxPyPz(pxyz_conv_test_inv);
          invert=kTRUE;
    }

        if(status)
        {
          if(!invert)
          {
            //paramNDGAr4D.Set(xyz_conv,pxyz_conv,covar,sign);
            //paramNDGAr4D.Rotate(0);
            paramNDGAr4D=param4D;
          }
          else
          {
            //paramNDGAr4D.Set(xyz_conv,pxyz_inv,covar,-sign);
            //paramNDGAr4D.Rotate(0);
            paramNDGAr4D=param4D_inv_perm;
            Double_t pxyz_conv_test_inv_again[3];
            paramNDGAr4D.GetPxPyPz(pxyz_conv_test_inv_again);
            int x=0;
          }
        }

    //paramNDGAr4D.Set(xyz_conv,pxyz_conv,covar,-sign);
    //bool statuscheck = paramNDGAr4D.Rotate(0);



    return 1;
}

bool inFiducial(double x, double y, double z) {
	double const Rcut = 400.0 /2.0;		double const Zcut = 439.4 /2.0;
	double r = sqrt( z*z + y*y );
	if ( r > Rcut ) return false;
	if (fabs(x) > Zcut ) return false;
	return true;
};

bool inTPC(double x, double y, double z) {
	double const Rcut = 500 /2.0;		double const Zcut = 470 /2.0;
	double r = sqrt( z*z + y*y );
	if ( r > Rcut ) return false;
	if (fabs(x) > Zcut ) return false;
	return true;
};

void testNDGAr(Int_t nEv, bool dumpStream=1, bool Ideal=kTRUE, size_t FirstEvent=0, bool Interaction=kFALSE, bool RotFrame = kFALSE){

  /// part.fStatusMaskIn[0]&8192>0&&part.fParamMC[0].GetP()<1000 conditions for drawing

  /*
  remake plots
  gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/aliKalman/test/\"")
  gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/MC/\"")
  gSystem->Load("$fastMCKalman/fastMCKalman/aliKalman/test/AliExternalTrackParam.so");
  .L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx++g
  .L $ConvertMC/testNDGAr.C++g
  AliPDG::AddParticlesToPdgDataBase();
  initTreeFast()
  treeFast->SetAlias("deltaP","(part.fParamMC[0].GetP()-part.fParamIn[0].GetP())/part.fParamMC[0].GetP()")
  treeFast->SetAlias("deltaPorig","(part.fParamMC[0].GetP()-paramSt.GetP())/part.fParamMC[0].GetP()")
  treeFast->Draw("deltaP:sqrt(part.fNPointsIn[0]/part.@fParamMC.size())>>hv(10,0,1,100,-0.4,0.4)","part.fStatusMaskIn[0]&8192>0&&part.fParamMC[0].GetP()<1000&&part.fNPointsIn[0]>50","colz")
  hv->FitSlicesY()
  treeFast->Draw("deltaPorig:sqrt(part.fNPointsIn[0]/part.@fParamMC.size())>>hk(10,0,1,100,-0.4,0.4)","part.fStatusMaskIn[0]&8192>0&&part.fParamMC[0].GetP()<1000&&part.fNPointsIn[0]>50","colz")
  hk->FitSlicesY()
  hv_2->Divide(hk_2)
  hv_2->SetTitle("Impact of skipping points on momentum resolution;#sqrt{N_{New}/N_{Old}};res_{New}/res_{Old}")
  hv_2->Draw()

  */
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
  std::string filename = "FullRecoInteraction";
  if(Ideal) filename+= "Ideal";
  std::string path ="/home/federico/Documents/Universita/Federico_2020-2021/Aliwork/ConvertMC/";
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
  if(Interaction) tree_source->Add("/home/federico/Documents/Universita/Federico_2020-2021/Aliwork/ConvertMC/MC/Interaction/*");
  else tree_source->Add("/home/federico/Documents/Universita/Federico_2020-2021/Aliwork/ConvertMC/MC/anatree_range*");
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

      if(Interaction && CCNC!=0 && GPartPdg->at(0)!=14) continue; //If analysing Interaction file only consider numuCC events   

      if(Interaction && !Ideal && VertX->size()==0) continue;
      if(Interaction && !Ideal) {
        bool vert_in_fid = inFiducial(MCVertZ->at(0)-GArCenter[2],MCVertY->at(0)-GArCenter[1],MCVertX->at(0)-GArCenter[0]);
        if(!vert_in_fid) continue; //If analysng realistic Interaction sample, apply fiducial cut
      }


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
        if(!Interaction && t>0) continue; //If analyzing particle gun file, only consider first track
        //if(t!=7) continue;
        long PDGcode=PDG->at(MCID_vector.at(t));    ////Save PDGCode for the MC trajectory associated with the ND-GAr reconstructed track
        long PDGMothercode=PDGMother->at(MCID_vector.at(t));    ////Save PDGCode for the MC trajectory associated with the ND-GAr reconstructed track
        // bool contained = inTPC(MCPStartZ->at(MCID_vector.at(t))-GArCenter[2],MCPStartY->at(MCID_vector.at(t))-GArCenter[1],MCPStartX->at(MCID_vector.at(t))-GArCenter[0]);
        // contained *= inTPC(MCPEndZ->at(MCID_vector.at(t))-GArCenter[2],MCPEndY->at(MCID_vector.at(t))-GArCenter[1],MCPEndX->at(MCID_vector.at(t))-GArCenter[0]);
        if(abs(PDGcode)!=13 || PDGMothercode!=0 ) continue;


          ///////Match MC trajectory that produced the track and save  the xyz's and pxyz's
          std::vector<TVector3> trajxyz;                
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

          

          if(!Interaction && abs(PDGcode)!=13 ) continue;   ////If analyzing partigle gun file, only consider muons 
          std::cout<<"ID: "<<ID_vector.at(t)<<" PDG: "<<PDGcode<<" ev: "<< i <<std::endl;




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
          int qStMC = (PDGcode == 13)?-1:1;
          BuildParamNDGArReco(paramSt0,GArCenter,xyzSt,pxyzSt,PDGcode,qSt);


          ///End of the track
          AliExternalTrackParam4D paramEnd0;
          TVector3 xyzEnd(TrackEndX->at(t),TrackEndY->at(t),TrackEndZ->at(t));
          TVector3 pxyzEnd(TrackEndPX->at(t),TrackEndPY->at(t),TrackEndPZ->at(t));
          int qEnd = TrackEndQ->at(t);
          TVector3 xyzEndMC(MCPEndX->at(MCID_vector.at(t)),MCPEndY->at(MCID_vector.at(t)),MCPEndZ->at(MCID_vector.at(t)));
          TVector3 pxyzEndMC(MCPEndPX->at(MCID_vector.at(t)),MCPEndPY->at(MCID_vector.at(t)),MCPEndPZ->at(MCID_vector.at(t)));
          int qEndMC = (PDGcode == 13)?-1:1;
          BuildParamNDGArReco(paramEnd0,GArCenter,xyzEnd,pxyzEnd,PDGcode,qEnd);

          AliExternalTrackParam4D paramSt;
          AliExternalTrackParam4D paramEnd;

          //NB: In GArSoft what is labeled as TrackStart is positioned at the last TPCCluster,
          //following the sort_TPCClusters_along_track2 forward ordering.
          //If the forward ordering matches the trajectory ordering IsForwardf=kTRUE,
          //hence, what we would consider to be the actual track start (i.e. the extremity closest to the trajectory start) 
          //is labelled by GArSoft as TrackEnd and viceversa

          if(!IsForwardf) 
          {
            paramSt = paramSt0;
            paramEnd = paramEnd0;
          }
          else
          {
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
            if(RotFrame){
              BuildParticle(particle,GArCenter,geom,trajxyz,trajpxyz,PDGcode);  //////Built the ALICE particle with the MC trajectory xyz and pxyz points 
              if(particle.fParamMC.size()==0) continue;
              particle.reconstructParticleFull(geom,PDGcode,10000);
              //particle.reconstructParticleFullOut(geom,PDGcode,10000);
              //particle.refitParticle();
            }else{
              BuildParticleNoRotation(particle,GArCenter,geom,trajxyz,trajpxyz,PDGcode,0);
              if(particle.fParamMC.size()==0) continue;

              particle.reconstructParticleRotate0(geom,PDGcode,10000);
              if (particle.fParamInRot[0].GetCovariance()[0]==0){
                for(int a=1;a<=50;a++)
                {
                  BuildParticleNoRotation(particle,GArCenter,geom,trajxyz,trajpxyz,PDGcode,a*Pi()/50);
                  particle.reconstructParticleRotate0(geom,PDGcode,10000);
                  if(particle.fParamInRot[0].GetCovariance()[0]!=0) {
                    break;
                  }
                }
              }
            }


            if (dumpStream==kFALSE) continue;
            if (tree) tree->Fill();
            else {
              (*pcstream) << "fastPart" <<
                          "i=" << i <<
                          "t=" <<t<<
                          //"paramGAr.=" <<&paramGAr<<
                          ///"paramGArbkw.=" <<&paramGArbkw<<
                          "paramSt.="<<&paramSt<<
                          "paramEnd.="<<&paramEnd<<
                          "paramStNoRot.="<<&paramStNoRot<<
                          "paramEndNoRot.="<<&paramEndNoRot<<
                          "paramStNoRotMC.="<<&paramStNoRotMC<<
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
            particle_f.fAddPadsmearing=true;
            particle_f.fUseGArSeeding=false;
            particle_f.fUseMCInfo=false;
            particle_f.fgStreamer=pcstream;
            particle_f.gid=ID_vector.at(t);
            particle_f.fDecayLength=0;
            fastParticle particle_f_nosmear(hlf.size()+1);
            particle_f_nosmear.fAddMSsmearing=true;
            particle_f_nosmear.fAddPadsmearing=false;
            particle_f_nosmear.fUseGArSeeding=false;
            particle_f_nosmear.fUseMCInfo=false;
            particle_f_nosmear.fgStreamer=pcstream;
            particle_f_nosmear.gid=ID_vector.at(t);
            particle_f_nosmear.fDecayLength=0;
            
            if(RotFrame){
              if(IsForwardf){
                std::cout<<"Forward"<<std::endl;
                BuildParticle(particle_f,GArCenter,geom,TrkClusterXYZf_NDGAr,trajpxyzf_NDGAr,PDGcode);  /////Build the ALICE particle with Track XYZClusters ordered forward and closest MC pxyz
                //SetParticleLoop(particle_f);
                //BuiltParticleUnsmeared(particle_f_nosmear,GArCenter,TrkClusterXYZf_NDGAr,trajxyzf_NDGAr);
                if(particle_f.fParamMC.size()==0) continue;
              }
              else{
                std::cout<<"Backwards"<<std::endl;
                BuildParticle(particle_f,GArCenter,geom,TrkClusterXYZb_NDGAr,trajpxyzb_NDGAr,PDGcode);  /////Build the ALICE particle with Track XYZClusters ordered forward and closest MC pxyz
                //SetParticleLoop(particle_f);
                //BuiltParticleUnsmeared(particle_f_nosmear,GArCenter,TrkClusterXYZb_NDGAr,trajxyzb_NDGAr);
                if(particle_f.fParamMC.size()==0) continue;
              }
              particle_f.reconstructParticleFull(geom,PDGcode,10000);
            }
            else{

            if(IsForwardf){
              std::cout<<"Forward"<<std::endl;
              BuildParticleNoRotation(particle_f,GArCenter,geom,TrkClusterXYZf_NDGAr,trajpxyzf_NDGAr,PDGcode);  /////Build the ALICE particle with Track XYZClusters ordered forward and closest MC pxyz
              //SetParticleLoop(particle_f);
              BuiltParticleUnsmearedNoRot(particle_f_nosmear,GArCenter,TrkClusterXYZf_NDGAr,trajxyzf_NDGAr);
              if(particle_f.fParamMC.size()==0) continue;
              particle_f.reconstructParticleRotate0(geom,PDGcode,10000);
              if (particle_f.fParamInRot[0].GetCovariance()[0]==0){
                for(int a=1;a<=50;a++)
                {
                  BuildParticleNoRotation(particle_f,GArCenter,geom,TrkClusterXYZf_NDGAr,trajpxyzf_NDGAr,PDGcode,a*Pi()/50);
                  particle_f.reconstructParticleRotate0(geom,PDGcode,10000);
                  BuiltParticleUnsmearedNoRot(particle_f_nosmear,GArCenter,TrkClusterXYZf_NDGAr,trajxyzf_NDGAr,a*Pi()/50);
                  if(particle_f.fParamInRot[0].GetCovariance()[0]!=0) {
                    break;
                  }
                }
              }
            }
            else{
              std::cout<<"Backwards"<<std::endl;
              BuildParticleNoRotation(particle_f,GArCenter,geom,TrkClusterXYZb_NDGAr,trajpxyzb_NDGAr,PDGcode);  /////Build the ALICE particle with Track XYZClusters ordered forward and closest MC pxyz
              //SetParticleLoop(particle_f);
              BuiltParticleUnsmearedNoRot(particle_f_nosmear,GArCenter,TrkClusterXYZb_NDGAr,trajxyzb_NDGAr);
              if(particle_f.fParamMC.size()==0) continue;
              particle_f.reconstructParticleRotate0(geom,PDGcode,10000);
              if (particle_f.fParamInRot[0].GetCovariance()[0]==0){
                for(int a=1;a<=50;a++)
                {
                  BuildParticleNoRotation(particle_f,GArCenter,geom,TrkClusterXYZb_NDGAr,trajpxyzb_NDGAr,PDGcode,a*Pi()/50);
                  particle_f.reconstructParticleRotate0(geom,PDGcode,10000);
                  BuiltParticleUnsmearedNoRot(particle_f_nosmear,GArCenter,TrkClusterXYZb_NDGAr,trajxyzb_NDGAr,a*Pi()/50);
                  if(particle_f.fParamInRot[0].GetCovariance()[0]!=0) {
                    break;
                  }
                }
              }
            }
            //particle_f.reconstructParticleRotate0(geom,PDGcode,10000);            
            
            //particle_f.reconstructParticleFullOut(geom,PDGcode,10000);
            //particle_f.refitParticle();
            }

            ////// Create particle object with ALICE-like global coordinates


            if (dumpStream==kFALSE) continue;
            if (tree) tree->Fill();
            else {
              (*pcstream) << "fastPart" <<
                          "i=" << i <<
                          "t=" <<t<<
                          "paramSt.="<<&paramSt<<
                          "paramEnd.="<<&paramEnd<<
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

  treeFast->SetAlias("ptMCSt","1./TMath::Abs(part.fParamMC[0].fP[4])");
  treeFast->SetAlias("rMCSt","TMath::Sqrt((1. - part.fParamMC[0].fP[2])*(1. + part.fParamMC[0].fP[2]))");
  treeFast->SetAlias("pxMCSt","ptMCSt*(rMCSt*TMath::Cos(part.fParamMC[0].fAlpha) - part.fParamMC[0].fP[2]*TMath::Sin(part.fParamMC[0].fAlpha))");
  treeFast->SetAlias("pyMCSt","ptMCSt*(part.fParamMC[0].fP[2]*TMath::Cos(part.fParamMC[0].fAlpha) + rMCSt*TMath::Sin(part.fParamMC[0].fAlpha))");

}
