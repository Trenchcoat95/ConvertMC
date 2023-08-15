/*
    gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/aliKalman/test/\"")
    gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/MC/\"")
    gSystem->Load("$fastMCKalman/fastMCKalman/aliKalman/test/AliExternalTrackParam.so");
    .L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx++g
    .L $ConvertMC/CalculateTKI.C++g
    CalculateTKI("ALICE")
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
#include "tools/PDGutils.h"
#include "tools/AnaFunctions.h"


void CalculateTKI(string Type = "ALICE"){
    std::string folder="/home/federico/Documents/Universita/Federico_2020-2021/Aliwork/ConvertMC/";
    std::string file = "FullRecoInteractionTKI";
    std::string inputList = folder+file+"/fastParticle.list";
    const char* inputListPath=gSystem->ExpandPathName(inputList.c_str());
    TChain * treeFast  = AliXRDPROOFtoolkit::MakeChainRandom(inputListPath,"fastPart",0,10000);
    fastParticle *part = 0;
    AliExternalTrackParam4D *paramSt = 0;
    Int_t i;
    ULong64_t t;
    int Target;
    float NuPx;
    float NuPy;
    float NuPz;
    treeFast->SetBranchAddress("i",&i);
    treeFast->SetBranchAddress("t",&t);
    treeFast->SetBranchAddress("Target",&Target);
    treeFast->SetBranchAddress("NuPx", &NuPx);
    treeFast->SetBranchAddress("NuPy", &NuPy);
    treeFast->SetBranchAddress("NuPz", &NuPz);
    treeFast->SetBranchAddress("part.",&part);
    treeFast->SetBranchAddress("paramSt.",&paramSt);

    // treeFast->GetEntry(0);
    // double NuE = sqrt(NuPx*NuPx+NuPy*NuPy+NuPz*NuPz + massPDG(numuPDG));
    // TLorentzVector* p4_Nu = new TLorentzVector(NuPx,NuPy,NuPz,NuE);
    //int ID_start = i;
    // int A_tgt = GetA_PDG(Target);
    // int Z_tgt = GetZ_PDG(Target);

    std::vector<TLorentzVector> xyz_ev;
    std::vector<TLorentzVector> p4_ev;
    std::vector<int> length_ev;
    std::vector<int> PDG_ev;
    std::vector<size_t> size_ev;
    std::vector<int> mass_ev;
    TH1D dptt_H("dpTT_numuCC_H","dpTT_numuCC_H",80,-400,400);
    TH1D dptt_C("dpTT_numuCC_C","dpTT_numuCC_C",80,-400,400);
    TH1D dptt_Ar("dpTT_numuCC_Ar","dpTT_numuCC_Ar",80,-400,400);

    treeFast->GetEntry(0);
    int check3 = 0;

    int ID_start = i;
    double NuE = sqrt(NuPx*NuPx+NuPy*NuPy+NuPz*NuPz + massPDG(numuPDG));
    TLorentzVector* p4_Nu = new TLorentzVector(NuPx,NuPy,NuPz,NuE);
    int A_tgt = GetA_PDG(Target);
    int Z_tgt = GetZ_PDG(Target);

    if(Type!="ALICE" && Type!="GAr" && Type!="MC"){
        std::cout<<"Didn't choose a valid reconstruction type"<<std::endl;
        return;
    }

    for(size_t e = 0; e<treeFast->GetEntries(); e++){
       treeFast->GetEntry(e);
       if(i==ID_start){

        check3++;
        double pxyz[3];
        if(Type=="ALICE") part->fParamIn[0].GetPxPyPz(pxyz);
        else if (Type=="MC") part->fParamMC[0].GetPxPyPz(pxyz);
        else if (Type=="GAr") paramSt->GetPxPyPz(pxyz);
        double Energy = sqrt(pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]+pxyz[1]*pxyz[1] + part->fMassMC);
        TLorentzVector p4(pxyz[2],pxyz[1],pxyz[0],Energy);
        p4_ev.push_back(p4);
        length_ev.push_back(part->fLengthInRot);
        PDG_ev.push_back(part->fPdgCodeRec);
        mass_ev.push_back(part->fMassMC);
        size_ev.push_back(part->fParamMC.size());
        std::cout<<"ID: "<<i<<" ntracks: "<<check3<<" array size: "<<p4_ev.size()<<" PDG: "<<part->fPdgCodeRec<<" Target: "<<Target<<std::endl;
       }

       if(i!=ID_start || e==(treeFast->GetEntries()-1)){
        if(check3==3) {
            TLorentzVector Scatter(0,0,0,0);
            TLorentzVector Recoil(0,0,0,0);
            int npart = p4_ev.size();
            bool Check = true;
            for(int c=0; c<npart;c++){
                Float_t KE_part = p4_ev[c].Energy()-mass_ev[c];
                Float_t px = p4_ev[c].Px();
                Float_t py = p4_ev[c].Py();
                Float_t pz = p4_ev[c].Pz();
                if(!(abs(px)>0.001 && abs(py)>0.001 && abs(pz)>0.001)) Check = kFALSE;
                if(KE_part<0.003) Check=kFALSE;
                //if(size_ev[c]<50) Check=kFALSE;
            }

            ///rough PID

            int longest = 0;
            float lengthmax = 0;
            for(int r=0; r<length_ev.size(); r++){
            float length = length_ev[r];
                if(length>lengthmax){
                    lengthmax=length;
                    longest=r;
                }
            }

            for(int c=0; c<p4_ev.size();c++){
                if(c==longest) Scatter = p4_ev[c];
                else Recoil += p4_ev.at(c);
            }
            TLorentzVector * p4_Scatter = new TLorentzVector(Scatter);
            TLorentzVector * p4_Recoil = new TLorentzVector(Recoil);

            float pxsc = p4_Scatter->Px();
            float pysc = p4_Scatter->Py();
            float pzsc = p4_Scatter->Pz();
            float pxrec = p4_Recoil->Px();
            float pyrec = p4_Recoil->Py();
            float pzrec = p4_Recoil->Pz();
            float pxnu = p4_Nu->Px();
            float pynu = p4_Nu->Py();
            float pznu = p4_Nu->Pz();


            double dalphat, dphit, dpt, dpTT, beamCalcP, IApN, recoilM, recoilP;
            dalphat = dphit = dpt = dpTT = beamCalcP = IApN = recoilM = recoilP = 0;

            AnaFunctions::getCommonTKI(A_tgt,Z_tgt,p4_Nu,p4_Scatter,p4_Recoil,dalphat,dphit,dpt,dpTT,beamCalcP,IApN,recoilM,recoilP);
            std::cout<<"dpTT value: "<<dpTT<<std::endl;
            if(Check){
                if(Z_tgt==GetZ_PDG(HydrogenPDG)) dptt_H.Fill(dpTT*1000);
                if(Z_tgt==GetZ_PDG(CarbonPDG)) dptt_C.Fill(dpTT*1000);
                if(Z_tgt==GetZ_PDG(ArgonPDG)) dptt_Ar.Fill(dpTT*1000);
            }
        }

        p4_ev.clear();
        mass_ev.clear();
        size_ev.clear();
        length_ev.clear();
        PDG_ev.clear();

        if(e!=(treeFast->GetEntries()-1)){
            ///Set event quantities for new event
            ID_start=i;
            NuE = sqrt(NuPx*NuPx+NuPy*NuPy+NuPz*NuPz + massPDG(numuPDG));
            p4_Nu->SetPxPyPzE(NuPx,NuPy,NuPz,NuE);
            A_tgt = GetA_PDG(Target);
            Z_tgt = GetZ_PDG(Target);
            check3=1;

            ///Set first particle of new event
            double pxyz[3];
            part->fParamMC[0].GetPxPyPz(pxyz);
            double Energy = sqrt(pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]+pxyz[1]*pxyz[1] + part->fMassMC);
            TLorentzVector p4(pxyz[2],pxyz[1],pxyz[0],Energy);
            p4_ev.push_back(p4);
            length_ev.push_back(part->fLengthInRot);
            PDG_ev.push_back(part->fPdgCodeRec);
            mass_ev.push_back(part->fMassMC);
            size_ev.push_back(part->fParamMC.size());
            std::cout<<"ID: "<<i<<" ntracks: "<<check3<<" array size: "<<p4_ev.size()<<" PDG: "<<part->fPdgCodeRec<<" Target: "<<Target<<std::endl;
        }

       }

    }

    std::string out = folder+file+"/TKIFileOut"+Type+".root";
    TFile * FileOut = new TFile(out.c_str(),"RECREATE");
    FileOut->cd();
    FileOut->mkdir("anaplots");
    FileOut->cd("anaplots");

    int Ar = dptt_Ar.GetEntries();
    int H = dptt_H.GetEntries();
    int C = dptt_C.GetEntries();
    std::cout<<"Ar: "<<Ar<<" C: "<<C<<"H: "<<H<<std::endl;

    dptt_Ar.Write();
    dptt_H.Write();
    dptt_C.Write();

    FileOut->Close();
}

 