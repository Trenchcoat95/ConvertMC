/*
  Here we would like to have tests - which checks internal consistence of the fastSimulation
  Dedicated debug streamaers to be used

    gSystem->Load("$fastMCKalman/fastMCKalman/aliKalman/test/AliExternalTrackParam.so");
    .L $ConvertMC/Local/fastSimulation.cxx+g
    .L $ConvertMC/testNDGAr.C+g
    .L $ConvertMC/test_fastSimulation.C+g
    initTreeFast();

    or

    gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/aliKalman/test/\"")
    gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/MC/\"")
    gSystem->Load("$fastMCKalman/fastMCKalman/aliKalman/test/AliExternalTrackParam.so");
    .L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx+g
    .L $ConvertMC/testNDGAr.C+g
    .L $ConvertMC/test_fastSimulation.C+g
    initTreeFast();
 */

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TChain.h"
#include "TTree.h"
#include "TRandom.h"
#include "TH1.h"
#include "TF1.h"
#include "TTreeStream.h"
#include "fastTracker.h"
#include "fastSimulation.h"
#include "TStyle.h"




extern TChain * treeSeed;
extern TChain * treeFast;


/// Test pull of seeds
/// the seeding pulls at the big snp are underestimated - because of rotation - we call test only for smaller |fP[2]|<0.6 further from 1
void testPullsSeed() {
  // Here we should make unit tests - it will work only up to some relative energy loss
  // treeSeed->Draw("log(paramSeed.P()/input.P()):log(input2.P()/input.P())","sign0>0");

  //
  TF1 *mygauss = new TF1("mygauss", "gaus");
  for (int version=0; version<=1; version++) {
    for (int iPar = 0; iPar <= 4; iPar++) {
      treeSeed->Draw(Form("(seed.fP[%d]-input.fP[%d])/sqrt(seed.fC[%d])>>his(100,-6,6)", iPar, iPar, AliExternalTrackParam::GetIndex(iPar, iPar)),
                     Form("version==%d&&abs(seed.fP[2])<0.6",version), "");
      treeSeed->GetHistogram()->Fit("mygauss", "q");
      bool isOK = abs(1 - mygauss->GetParameter(2)) < 5 * mygauss->GetParError(2);
      float rms=treeSeed->GetHistogram()->GetRMS();
      if (isOK) {
        ::Info(Form("testFastTracker seed pull test P%d - version %d",iPar,version), "pullAnalytical - OK - %2.2f\t%2.2f", mygauss->GetParameter(2),rms);
      } else {
        ::Error(Form("testFastTracker seed pull test P%d - version %d",iPar,version), "pullAnalytical- FAILED- %2.2f\t%2.2f", mygauss->GetParameter(2),rms);
      }
    }
  }
}

void drawTrackStatus(int counter, std::string Id = "In", std::string Error = "0x1"){
  treeFast->SetMarkerColor(1);   /// all MC
  treeFast->Draw("gyMC:gxMC","","",1,counter);
  treeFast->SetMarkerColor(4);   /// all reco points
  treeFast->Draw("gyInF:gxInF",Form("part.fStatusMask%s>0",Id.c_str()),"same",1,counter);
  treeFast->SetMarkerColor(2);   /// first MC point
  treeFast->Draw("gyMC:gxMC","Iteration$==0","same",1,counter);
  treeFast->SetMarkerColor(3);   /// trigger problem
  treeFast->Draw("gyInF:gxInF",Form("part.fStatusMask%s==%s",Id.c_str(),Error.c_str()),"same",1,counter);
}

void drawTrackStatus3D(int counter, std::string Id = "In", std::string Error = "0x1"){
  treeFast->SetMarkerColor(1);   /// all MC
  treeFast->Draw("gyMC:gxMC:gzMC","","",1,counter);
  treeFast->SetMarkerColor(4);   /// all reco points
  treeFast->Draw("gyInF:gxInF:gzInF",Form("part.fStatusMask%s>0",Id.c_str()),"same",1,counter);
  treeFast->SetMarkerColor(2);   /// first MC point
  treeFast->Draw("gyMC:gxMC:gzMC","Iteration$==0","same",1,counter);
  treeFast->SetMarkerColor(3);   /// trigger problem
  treeFast->Draw("gyInF:gxInF:gzInF",Form("part.fStatusMask%s==%s",Id.c_str(),Error.c_str()),"same",1,counter);
}

void SetList(std::string Id = "In", std::string Error = "0x1"){

  treeFast->SetAlias("gyInF",Form("sin(part.fParam%s[].fAlpha)*part.fParam%s[].fX",Id.c_str(),Id.c_str()));
  treeFast->SetAlias("gxInF",Form("cos(part.fParam%s[].fAlpha)*part.fParam%s[].fX",Id.c_str(),Id.c_str()));
  treeFast->SetAlias("gzInF",Form("part.fParam%s[].fP[1]",Id.c_str()));

  treeFast->Draw(">>ProblemRot",Form("Sum$(part.fStatusMask%s==%s)",Id.c_str(),Error.c_str()),"entrylist");
  TEntryList* problemList0x1 =(TEntryList*)gDirectory->Get("ProblemRot");
  treeFast->SetEntryList(problemList0x1);
  int counter=0;
  treeFast->SetMarkerSize(1.5);
  gStyle->SetPalette(55);
}

