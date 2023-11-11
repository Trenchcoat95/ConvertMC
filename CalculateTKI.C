/*
    gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/aliKalman/test/\"")
    gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/MC/\"")
    gSystem->Load("$fastMCKalman/fastMCKalman/aliKalman/test/AliExternalTrackParam.so");
    .L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx++g
    .L $ConvertMC/CalculateTKI.C++g

    gSystem->AddIncludePath("-I\"test/\"")
    gSystem->Load("test/AliExternalTrackParam.so");
    .L fastSimulation.cxx++g
    .L CalculateTKI.C++g 

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
#include "TVector.h"
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
#include "TEntryList.h"
#include "TLorentzVector.h"


void CalculateTKI(string Type = "ALICE", bool Track_Cuts = true,  bool PIDcut = false, bool Reco_nuDir = false){
    std::string folder="";
    std::string file = "FullRecoInteractionTKI";
    std::string inputList = folder+file+"/fastParticle.list";
    const char* inputListPath=gSystem->ExpandPathName(inputList.c_str());
    TChain * treeFast  = AliXRDPROOFtoolkit::MakeChainRandom(inputListPath,"fastPart",0,10000);
    fastParticle *part = 0;
    AliExternalTrackParam4D *paramSt = 0;
    AliExternalTrackParam4D *paramStMC = 0;
    Int_t i;
    ULong64_t t;
    int Target;
    int nPrimaries;
    long PDGMotherCode;
    double NuPReco[] = {0,-0.100988589,0.994887584};
    float NuPx;
    float NuPy;
    float NuPz;
    int Flip;
    double TrackLenSt;
    Double_t xstart;
    Double_t ystart;
    float VertGArx;
    float VertGAry;
    float VertGArz;
    treeFast->SetBranchAddress("i",&i);
    treeFast->SetBranchAddress("t",&t);
    treeFast->SetBranchAddress("Target",&Target);
    treeFast->SetBranchAddress("NuPx", &NuPx);
    treeFast->SetBranchAddress("NuPy", &NuPy);
    treeFast->SetBranchAddress("NuPz", &NuPz);
    treeFast->SetBranchAddress("part.",&part);
    treeFast->SetBranchAddress("paramSt.",&paramSt);
    treeFast->SetBranchAddress("paramStMC.",&paramStMC);
    treeFast->SetBranchAddress("TrackLenSt",&TrackLenSt);
    treeFast->SetBranchAddress("xstart",&xstart);
    treeFast->SetBranchAddress("ystart",&ystart);
    treeFast->SetBranchAddress("nPrimaries",&nPrimaries);
    treeFast->SetBranchAddress("Flip",&Flip);
    treeFast->SetBranchAddress("PDGMotherCode",&PDGMotherCode);
    treeFast->SetBranchAddress("VertGArx",&VertGArx);
    treeFast->SetBranchAddress("VertGAry",&VertGAry);
    treeFast->SetBranchAddress("VertGArz",&VertGArz);

    // treeFast->GetEntry(0);
    // double NuE = sqrt(NuPx*NuPx+NuPy*NuPy+NuPz*NuPz + massPDG(numuPDG));
    // TLorentzVector* p4_Nu = new TLorentzVector(NuPx,NuPy,NuPz,NuE);
    //int ID_start = i;
    // int A_tgt = GetA_PDG(Target);
    // int Z_tgt = GetZ_PDG(Target);

    std::vector<TVector3> xyz_ev;
    std::vector<TVector3> xyzstart_ev;
    std::vector<TLorentzVector> p4_ev;
    std::vector<double> length_ev;
    std::vector<double> lArmMC_ev;
    std::vector<int> PDG_ev;
    std::vector<long> PDGMother_ev;
    std::vector<size_t> size_ev;
    std::vector<size_t> nPrimaries_ev;
    std::vector<double> mass_ev;
    std::vector<AliExternalTrackParam4D> param_ev;

    TH1D dptt_H("dpTT_numuCC_H","dpTT_numuCC_H",87,-100,100);
    if(Type=="MCSt") dptt_H = TH1D("dpTT_numuCC_H","dpTT_numuCC_H",101,-0.0005,0.0005);
    TH1D dptt_C("dpTT_numuCC_C","dpTT_numuCC_C",100,-100,100);
    TH1D dptt_Ar("dpTT_numuCC_Ar","dpTT_numuCC_Ar",100,-100,100);

    TH2D dpttVSAvgN_H("dpTTVSAvgN_numuCC_H","dpTTVSAvgN_numuCC_H",8,50,450,25,-60,60);
    TH2D dpttVSTotN_H("dpTTVSTotN_numuCC_H","dpTTVSTotN_numuCC_H",6,0,1200,25,-60,60);
    TH2D dpttVSN2212_H("dpTTVSN2212_numuCC_H","dpTTVSN2212_numuCC_H",8,50,450,25,-60,60);
    TH2D dpttVSN211_H("dpTTVSN211_numuCC_H","dpTTVSN211_numuCC_H",8,50,450,25,-60,60);
    TH2D dpttVSN13_H("dpTTVSN13_numuCC_H","dpTTVSN13_numuCC_H",8,50,450,25,-60,60);

    TH2D dpttVSAvgp_H("dpTTVSAvgp_numuCC_H","dpTTVSAvgp_numuCC_H",9,0.25,2.5,25,-60,60);
    TH2D dpttVSTotp_H("dpTTVSTotp_numuCC_H","dpTTVSTotp_numuCC_H",9,1.6,3.4,25,-60,60);
    TH2D dpttVSp2212_H("dpTTVSp2212_numuCC_H","dpTTVSp2212_numuCC_H",10,0.2,2.2,25,-60,60);
    TH2D dpttVSp211_H("dpTTVSp211_numuCC_H","dpTTVSp211_numuCC_H",6,0,1.2,25,-60,60);
    TH2D dpttVSp13_H("dpTTVSp13_numuCC_H","dpTTVSp13_numuCC_H",7,0,3.5,25,-60,60);

    TH2D dpttVSAvglArmMC_H("dpTTVSAvglArmMC_numuCC_H","dpTTVSAvglArmMC_numuCC_H",5,50,300,25,-60,60);
    TH2D dpttVSTotlArmMC_H("dpTTVSTotlArmMC_numuCC_H","dpTTVSTotlArmMC_numuCC_H",9,100,1000,25,-60,60);
    TH2D dpttVSlArmMC2212_H("dpTTVSlArmMC2212_numuCC_H","dpTTVSlArmMC2212_numuCC_H",8,50,450,25,-60,60);
    TH2D dpttVSlArmMC211_H("dpTTVSlArmMC211_numuCC_H","dpTTVSlArmMC211_numuCC_H",8,50,450,25,-60,60);
    TH2D dpttVSlArmMC13_H("dpTTVSlArmMC13_numuCC_H","dpTTVSlArmMC13_numuCC_H",8,50,450,25,-60,60);


    treeFast->GetEntry(0);
    int check3 = 0;

    int ID_start = i;
    double NuE = sqrt(NuPx*NuPx+NuPy*NuPy+NuPz*NuPz + massPDG(numuPDG));
    TLorentzVector* p4_Nu = new TLorentzVector(NuPx,NuPy,NuPz,NuE);
    if(Reco_nuDir) p4_Nu->SetPxPyPzE(NuPReco[0],NuPReco[1],NuPReco[2],NuE);
    int A_tgt = GetA_PDG(Target);
    int Z_tgt = GetZ_PDG(Target);

    if(Type!="ALICE" && Type!="GAr" && Type!="MC" && Type!="MCSt" && Type!= "ALICEStMC" && Type!= "ALICEStReco" && Type!= "ALICEStRecoGAr"){
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
        else if (Type=="MCSt") {
            paramStMC->GetPxPyPz(pxyz);
            pxyz[0]*=Flip;
            pxyz[1]*=Flip;
            pxyz[2]*=Flip;
        }
        else if (Type=="ALICEStMC") {
            AliExternalTrackParam4D p = part->fParamIn[0];
            double xyz[3];
            paramStMC->GetXYZ(xyz);
            //std::cout<<"Track Start MC: "<<xyz[0]-xstart<<" "<<xyz[1]-ystart<<" "<<xyz[2]<<std::endl;
            float xx0 = 8.37758e-04;
            float xrho = 0.016770000;
            float ResolRPhi = 0.4;
            float ResolZ = 0.3;
            float Bz = -5;
            double chi2Cut= 10000/ResolZ;
            float mass = part->fMassRec;
            MoveParamToPoint(p, xyz, xx0, xrho, ResolRPhi, ResolZ, Bz, chi2Cut, mass);
            p.GetPxPyPz(pxyz);
        }
        else if (Type=="ALICEStRecoGAr") {
            AliExternalTrackParam4D p = part->fParamIn[0];
            double xyz[3] = {VertGArx+xstart,VertGAry+ystart,VertGArz};
            //std::cout<<"Track Start MC: "<<xyz[0]-xstart<<" "<<xyz[1]-ystart<<" "<<xyz[2]<<std::endl;
            float xx0 = 8.37758e-04;
            float xrho = 0.016770000;
            float ResolRPhi = 0.4;
            float ResolZ = 0.3;
            float Bz = -5;
            double chi2Cut= 10000/ResolZ;
            float mass = part->fMassRec;
            MoveParamToPoint(p, xyz, xx0, xrho, ResolRPhi, ResolZ, Bz, chi2Cut, mass);
            p.GetPxPyPz(pxyz);
        }
        else if (Type=="ALICEStReco"){
            param_ev.push_back(part->fParamIn[0]);
            part->fParamIn[0].GetPxPyPz(pxyz);
            TVector3 xyzstartvec(xstart,ystart,0);
            xyzstart_ev.push_back(xyzstartvec);
        }
        double Energy = sqrt(pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]+pxyz[1]*pxyz[1] + part->fMassRec);
        TLorentzVector p4(pxyz[2],pxyz[1],pxyz[0],Energy);
        p4_ev.push_back(p4);
        if(Type!="GAr") length_ev.push_back(part->fLengthInRot);
        else length_ev.push_back(TrackLenSt);
        PDG_ev.push_back(part->fPdgCodeRec);
        PDGMother_ev.push_back(PDGMotherCode);
        mass_ev.push_back(part->fMassRec);
        size_ev.push_back(part->fParamMC.size());
        nPrimaries_ev.push_back(nPrimaries);

        float yend = part->fParamMC[part->fParamMC.size()-1].GetX()*sin(part->fParamMC[part->fParamMC.size()-1].GetAlpha())+part->fParamMC[part->fParamMC.size()-1].GetParameter()[0]*cos(part->fParamMC[part->fParamMC.size()-1].GetAlpha());
        float xend = part->fParamMC[part->fParamMC.size()-1].GetX()*cos(part->fParamMC[part->fParamMC.size()-1].GetAlpha())-part->fParamMC[part->fParamMC.size()-1].GetParameter()[0]*sin(part->fParamMC[part->fParamMC.size()-1].GetAlpha());
        float y0 = part->fParamMC[0].GetX()*sin(part->fParamMC[0].GetAlpha())+part->fParamMC[0].GetParameter()[0]*cos(part->fParamMC[0].GetAlpha());
        float x0 = part->fParamMC[0].GetX()*cos(part->fParamMC[0].GetAlpha())-part->fParamMC[0].GetParameter()[0]*sin(part->fParamMC[0].GetAlpha());
        float lArmMC = sqrt((xend-x0)*(xend-x0)+(yend-y0)*(yend-y0));
        lArmMC_ev.push_back(lArmMC);

        double xyzst[3];
        part->fParamIn[0].GetXYZ(xyzst);
        TVector3 xyz_vec(xyzst[0]-xstart,xyzst[1]-ystart,xyzst[2]);
        xyz_ev.push_back(xyz_vec);
        //std::cout<<"ID: "<<i<<" ntracks: "<<check3<<" array size: "<<p4_ev.size()<<" PDG: "<<part->fPdgCodeRec<<" Target: "<<Target<<std::endl;
       }

       if(i!=ID_start || e==(treeFast->GetEntries()-1)){
        if(check3==3) {
            TLorentzVector Scatter(0,0,0,0);
            TLorentzVector Recoil(0,0,0,0);
            int npart = p4_ev.size();
            bool Check = true;
            bool Check_Cut = true;
            for(int c=0; c<npart;c++){
                Float_t KE_part = p4_ev[c].Energy()-mass_ev[c];
                Float_t px = p4_ev[c].Px();
                Float_t py = p4_ev[c].Py();
                Float_t pz = p4_ev[c].Pz();
                //std::cout<<"Momenta: "<<px<<" "<<py<<" "<<pz<<" "<<std::endl;
                if(Track_Cuts){
                    if(!(abs(px)>0.001 && abs(py)>0.001 && abs(pz)>0.001 && KE_part>0.003)) Check = kFALSE;
                    if(!(lArmMC_ev[c]>50 && size_ev[c]>50)) Check_Cut=kFALSE;
                }
                //if(KE_part<0.003) Check=kFALSE;
                //if(size_ev[c]<50) Check=kFALSE;
            }

            float totsize = 0;
            for(auto s : size_ev) totsize+=s;
            float avgsize=totsize/size_ev.size();

            float totlArmMC = 0;
            for(auto s : lArmMC_ev) totlArmMC+=s;
            float avglArmMC=totlArmMC/size_ev.size();
            //std::cout<<"Avg lArm"<<avglArmMC<<std::endl;

            TLorentzVector totmomentum(0,0,0,0);
            float avgmomentum = 0;
            float totmomentum_mod = 0;
            for(auto s : p4_ev) {
                totmomentum=totmomentum+s;
                avgmomentum+=sqrt(s.Px()*s.Px()+s.Py()*s.Py()+s.Pz()*s.Pz());
            }
            avgmomentum /= p4_ev.size();
            totmomentum_mod = sqrt(totmomentum.Px()*totmomentum.Px()+totmomentum.Py()*totmomentum.Py()+totmomentum.Pz()*totmomentum.Pz());
            //if(!(totmomentum_mod>2)) Check_Cut=kFALSE;
            //std::cout<<"Tot momentum"<<totmomentum_mod<<std::endl;


            ///rough PID

            int longest = 0;
            float lengthmax = 0;
            for(int r=0; r<length_ev.size(); r++){
            float length = length_ev[r];
                if(length>lengthmax){
                    lengthmax=length;
                    longest=r;
                }
                //std::cout<<"Length: "<<length<<" r: "<<r<<std::endl;
            }
            //std::cout<<"longest: "<<longest<<std::endl;

            ////correct PID's cut
            bool ismuon = false;
            bool isproton = false;
            bool ispiplus = false;
            bool allprimaries = true;
            for(auto p:PDG_ev){
                if(p==muonPDG) ismuon = true;
                if(p==pichPDG) ispiplus = true;
                if(p== protonPDG) isproton = true;
            }

            for(auto p:PDGMother_ev){
               if (p!=0) allprimaries=false;
            }

            if(PIDcut && !(ismuon&&isproton&&ispiplus&&nPrimaries_ev[0]==3&&allprimaries)) Check_Cut=kFALSE;


            /////Calculate vertex
            Double_t vert_reco[3];
            //Vertexing3D(p4_ev,xyz_ev,vert_reco);
            //findIntersection(xyz_ev[0],p4_ev[0],xyz_ev[1],p4_ev[1],xyz_ev[2],p4_ev[2],vert_reco);
            double smallestx = xyz_ev[0].X();
            double smallest_index=0;
            for (int q=0;q<xyz_ev.size();q++){
                if(xyz_ev[q].X()<smallestx){
                  smallestx = xyz_ev[q].X();
                  smallest_index = q;
                }
            }
            vert_reco[0] = xyz_ev[smallest_index].X(); 
            vert_reco[1] = xyz_ev[smallest_index].Y(); 
            vert_reco[2] = xyz_ev[smallest_index].Z(); 
            double vertcheck[3];
            paramStMC->GetXYZ(vertcheck);

            //for(auto a:xyz_ev) std::cout<<"Track Start: "<<a.X()<<" "<<a.Y()<<" "<<a.Z()<<std::endl;
            //std::cout<<"Reco Vert: "<<(vert_reco[0])<<" "<<(vert_reco[1])<<" "<<vert_reco[2]<<std::endl<<std::endl;

            if(Type=="ALICEStReco"){
                p4_ev.clear();
                for(int k=0;k<param_ev.size();k++){
                    AliExternalTrackParam4D p = param_ev[k];
                    double xyz[3] = {vert_reco[0]+xyzstart_ev[k].X(),vert_reco[1]+xyzstart_ev[k].Y(),vert_reco[2]+xyzstart_ev[k].Z()};
                    //std::cout<<"Track Start MC: "<<xyz[0]-xstart<<" "<<xyz[1]-ystart<<" "<<xyz[2]<<std::endl;
                    float xx0 = 8.37758e-04;
                    float xrho = 0.016770000;
                    float ResolRPhi = 0.4;
                    float ResolZ = 0.3;
                    float Bz = -5;
                    double chi2Cut= 10000/ResolZ;
                    float mass = mass_ev[k];
                    MoveParamToPoint(p, xyz, xx0, xrho, ResolRPhi, ResolZ, Bz, chi2Cut, mass);
                    double pxyznew[3];
                    p.GetPxPyPz(pxyznew);
                    double Energynew = sqrt(pxyznew[0]*pxyznew[0]+pxyznew[1]*pxyznew[1]+pxyznew[1]*pxyznew[1] + mass_ev[k]);
                    TLorentzVector p4(pxyznew[2],pxyznew[1],pxyznew[0],Energynew);
                    p4_ev.push_back(p4);
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
            if(Check){
                if(Check_Cut){
                std::cout<<"p4_beam: "<<pxnu<<" "<<pynu<<" "<<pznu<<std::endl;
                std::cout<<"p4_scatter: "<<pxsc<<" "<<pysc<<" "<<pzsc<<std::endl;
                std::cout<<"p4_recoil: "<<pxrec<<" "<<pyrec<<" "<<pzrec<<std::endl;
                std::cout<<"dpTT value: "<<dpTT<<" IDStart: "<< ID_start<<std::endl<<std::endl;
                }
                if(Z_tgt==GetZ_PDG(HydrogenPDG)) {
                    if(abs(dpTT)>0.001 && Check_Cut){
                        ofstream myfile;
                        myfile.open ("weird.txt", std::ios_base::app);
                        myfile << ID_start << "\n";
                        myfile.close();
                    }

                    if(Check_Cut) dptt_H.Fill(dpTT*1000);
                    if(Check_Cut){
                        dpttVSAvgN_H.Fill(avgsize,dpTT*1000);
                        dpttVSTotN_H.Fill(totsize,dpTT*1000);
                        dpttVSAvgp_H.Fill(avgmomentum,dpTT*1000);
                        dpttVSTotp_H.Fill(totmomentum_mod,dpTT*1000);
                        dpttVSAvglArmMC_H.Fill(avglArmMC,dpTT*1000);
                        dpttVSTotlArmMC_H.Fill(totlArmMC,dpTT*1000);
                    }
                    if(Check_Cut){
                    for(int p =0; p<PDG_ev.size(); p++){
                            if(abs(PDG_ev[p])==muonPDG) {
                                dpttVSN13_H.Fill(size_ev[p],dpTT*1000); 
                                dpttVSp13_H.Fill(sqrt(p4_ev[p].Px()*p4_ev[p].Px()+p4_ev[p].Py()*p4_ev[p].Py()+p4_ev[p].Pz()*p4_ev[p].Pz()),dpTT*1000);
                                dpttVSlArmMC13_H.Fill(lArmMC_ev[p],dpTT*1000);
                            }
                            if(abs(PDG_ev[p])==pichPDG) {
                                dpttVSN211_H.Fill(size_ev[p],dpTT*1000);
                                dpttVSp211_H.Fill(sqrt(p4_ev[p].Px()*p4_ev[p].Px()+p4_ev[p].Py()*p4_ev[p].Py()+p4_ev[p].Pz()*p4_ev[p].Pz()),dpTT*1000);
                                dpttVSlArmMC211_H.Fill(lArmMC_ev[p],dpTT*1000);
                            }
                            if(abs(PDG_ev[p])==protonPDG) {
                                dpttVSN2212_H.Fill(size_ev[p],dpTT*1000);
                                dpttVSp2212_H.Fill(sqrt(p4_ev[p].Px()*p4_ev[p].Px()+p4_ev[p].Py()*p4_ev[p].Py()+p4_ev[p].Pz()*p4_ev[p].Pz()),dpTT*1000);
                                dpttVSlArmMC2212_H.Fill(lArmMC_ev[p],dpTT*1000);
                            } 
                        }
                    }
                }
                
                if(Z_tgt==GetZ_PDG(CarbonPDG)) dptt_C.Fill(dpTT*1000);
                if(Z_tgt==GetZ_PDG(ArgonPDG)) dptt_Ar.Fill(dpTT*1000);
            }
        }

        p4_ev.clear();
        mass_ev.clear();
        size_ev.clear();
        length_ev.clear();
        PDG_ev.clear();
        PDGMother_ev.clear();
        lArmMC_ev.clear();
        xyz_ev.clear();
        param_ev.clear();
        xyzstart_ev.clear();
        nPrimaries_ev.clear();

        if(e!=(treeFast->GetEntries()-1)){
            ///Set event quantities for new event
            ID_start=i;
            NuE = sqrt(NuPx*NuPx+NuPy*NuPy+NuPz*NuPz + massPDG(numuPDG));
            p4_Nu->SetPxPyPzE(NuPx,NuPy,NuPz,NuE);
            if(Reco_nuDir) p4_Nu->SetPxPyPzE(NuPReco[0],NuPReco[1],NuPReco[2],NuE);
            //p4_Nu->SetPxPyPzE(0,-0.1,0.995,NuE);
            A_tgt = GetA_PDG(Target);
            Z_tgt = GetZ_PDG(Target);
            check3=1;

            ///Set first particle of new event
            double pxyz[3];
            if(Type=="ALICE") part->fParamIn[0].GetPxPyPz(pxyz);
            else if (Type=="MC") part->fParamMC[0].GetPxPyPz(pxyz);
            else if (Type=="GAr") paramSt->GetPxPyPz(pxyz);
            else if (Type=="MCSt") {
                paramStMC->GetPxPyPz(pxyz);
                pxyz[0]*=Flip;
                pxyz[1]*=Flip;
                pxyz[2]*=Flip;
            }
            else if (Type=="ALICEStMC") {
                AliExternalTrackParam4D p = part->fParamIn[0];
                double xyz[3];
                paramStMC->GetXYZ(xyz);
                //std::cout<<"Track Start MC: "<<xyz[0]-xstart<<" "<<xyz[1]-ystart<<" "<<xyz[2]<<std::endl;
                float xx0 = 8.37758e-04;
                float xrho = 0.016770000;
                float ResolRPhi = 0.4;
                float ResolZ = 0.3;
                float Bz = -5;
                double chi2Cut= 10000/ResolZ;
                float mass = part->fMassRec;
                MoveParamToPoint(p, xyz, xx0, xrho, ResolRPhi, ResolZ, Bz, chi2Cut, mass);
                p.GetPxPyPz(pxyz);
            }
            else if (Type=="ALICEStRecoGAr") {
                AliExternalTrackParam4D p = part->fParamIn[0];
                double xyz[3] = {VertGArx+xstart,VertGAry+ystart,VertGArz};
                //std::cout<<"Track Start MC: "<<xyz[0]-xstart<<" "<<xyz[1]-ystart<<" "<<xyz[2]<<std::endl;
                float xx0 = 8.37758e-04;
                float xrho = 0.016770000;
                float ResolRPhi = 0.4;
                float ResolZ = 0.3;
                float Bz = -5;
                double chi2Cut= 10000/ResolZ;
                float mass = part->fMassRec;
                MoveParamToPoint(p, xyz, xx0, xrho, ResolRPhi, ResolZ, Bz, chi2Cut, mass);
                p.GetPxPyPz(pxyz);
            }
            else if (Type=="ALICEStReco"){
                param_ev.push_back(part->fParamIn[0]);
                part->fParamIn[0].GetPxPyPz(pxyz);
                TVector3 xyzstartvec(xstart,ystart,0);
                xyzstart_ev.push_back(xyzstartvec);
            }
            double Energy = sqrt(pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]+pxyz[1]*pxyz[1] + part->fMassRec);
            TLorentzVector p4(pxyz[2],pxyz[1],pxyz[0],Energy);
            p4_ev.push_back(p4);
            if(Type!="GAr") length_ev.push_back(part->fLengthInRot);
            else length_ev.push_back(TrackLenSt);
            PDG_ev.push_back(part->fPdgCodeRec);
            PDGMother_ev.push_back(PDGMotherCode);
            mass_ev.push_back(part->fMassRec);
            size_ev.push_back(part->fParamMC.size());
            nPrimaries_ev.push_back(nPrimaries);

            float yend = part->fParamMC[part->fParamMC.size()-1].GetX()*sin(part->fParamMC[part->fParamMC.size()-1].GetAlpha())+part->fParamMC[part->fParamMC.size()-1].GetParameter()[0]*cos(part->fParamMC[part->fParamMC.size()-1].GetAlpha());
            float xend = part->fParamMC[part->fParamMC.size()-1].GetX()*cos(part->fParamMC[part->fParamMC.size()-1].GetAlpha())-part->fParamMC[part->fParamMC.size()-1].GetParameter()[0]*sin(part->fParamMC[part->fParamMC.size()-1].GetAlpha());
            float y0 = part->fParamMC[0].GetX()*sin(part->fParamMC[0].GetAlpha())+part->fParamMC[0].GetParameter()[0]*cos(part->fParamMC[0].GetAlpha());
            float x0 = part->fParamMC[0].GetX()*cos(part->fParamMC[0].GetAlpha())-part->fParamMC[0].GetParameter()[0]*sin(part->fParamMC[0].GetAlpha());
            float lArmMC = sqrt((xend-x0)*(xend-x0)+(yend-y0)*(yend-y0));
            lArmMC_ev.push_back(lArmMC);

            double xyzst[3];
            part->fParamIn[0].GetXYZ(xyzst);
            TVector3 xyz_vec(xyzst[0]-xstart,xyzst[1]-ystart,xyzst[2]);
            xyz_ev.push_back(xyz_vec);
            //std::cout<<"ID: "<<i<<" ntracks: "<<check3<<" array size: "<<p4_ev.size()<<" PDG: "<<part->fPdgCodeRec<<" Target: "<<Target<<std::endl;
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

    dpttVSAvgN_H.Write();
    dpttVSTotN_H.Write();
    dpttVSN13_H.Write();
    dpttVSN211_H.Write();
    dpttVSN2212_H.Write();

    dpttVSAvgp_H.Write();
    dpttVSTotp_H.Write();
    dpttVSp13_H.Write();
    dpttVSp211_H.Write();
    dpttVSp2212_H.Write();

    dpttVSAvglArmMC_H.Write();
    dpttVSTotlArmMC_H.Write();
    dpttVSlArmMC13_H.Write();
    dpttVSlArmMC211_H.Write();
    dpttVSlArmMC2212_H.Write();

    FileOut->Close();
}


void AnalyzefromTxt(std::string txt, int maxEntries){
    TChain* t = new TChain("/anatree/GArAnaTree");
    t->Add("/pnfs/dune/persistent/users/battisti/TKI/traj/*");
    ifstream myfile(txt,ios::in);
    vector<int> entries;
    int num;

    while(myfile){
       myfile>>num;
       entries.push_back(num);
    }

    std::vector<int>     *PDG = 0;
    std::vector<int>     *TgtPDG = 0;
    std::vector<int>     *PDGMother = 0;
    std::vector<int>   *GPartPdg   =0;
    std::vector<int>   *CCNC   =0;
    std::vector<int>     *TrackMCindex = 0;

    t->SetBranchAddress("PDG",&PDG);
    t->SetBranchAddress("TgtPDG",&TgtPDG);
    t->SetBranchAddress("PDGMother",&PDGMother);
    t->SetBranchAddress("GPartPdg",&GPartPdg);
    t->SetBranchAddress("CCNC",&CCNC);
    t->SetBranchAddress("TrackMCindex",&TrackMCindex);
    
    int count = 0;
    int NCcount = 0;
    int antinumucount = 0;
    int nuecount = 0;
    int mu_proton_piplus_pi0 =0;
    int mu_proton_piplus_2pi0 =0;
    int mu_proton_piplus_3pi0 =0;
    for (auto e:entries){
         count++;
         if(count>maxEntries) break;
         bool status = t->GetEntry(e);
         if(!status) continue;

         if(CCNC->at(0)==1){
            NCcount++;
            continue;
         }

        if(GPartPdg->at(0)==-numuPDG){
            antinumucount++;
            continue;
         }
        
        if(GPartPdg->at(0)==nuePDG){
            nuecount++;
            continue;
         }

        int checkmuon=0;
        int checkproton=0;
        int checkpiplus=0;
        int checkpi0=0;
        int nprimaries=0;

        for(int a = 0; a < PDGMother->size(); a++){
            if (PDGMother->at(a) != 0) continue;
            if(PDG->at(a)==muonPDG) checkmuon++;
            if(PDG->at(a)==pichPDG) checkpiplus++;
            if(PDG->at(a)==pi0PDG)  checkpi0++;
            if(PDG->at(a)==protonPDG) checkproton++;
            nprimaries++;
            //std::cout<<" "<<PDG->at(a);
         }

         if(checkmuon==1 && checkpiplus==1 && checkpi0==1 && checkproton==1 && nprimaries==(checkproton+checkmuon+checkpiplus+checkpi0)){
            mu_proton_piplus_pi0++;
            continue;
         }

        if(checkmuon==1 && checkpiplus==1 && checkpi0==2 && checkproton==1 && nprimaries==(checkproton+checkmuon+checkpiplus+checkpi0)){
            mu_proton_piplus_3pi0++;
            continue;
        }

        if(checkmuon==1 && checkpiplus==1 && checkpi0==3 && checkproton==1 && nprimaries==(checkproton+checkmuon+checkpiplus+checkpi0)){
            mu_proton_piplus_3pi0++;
            continue;
        }

         std::cout<<"EventID: "<<e<<" TgtPDG: "<<TgtPDG->at(0)<<" CCNC: "<<CCNC->at(0)<<" First part PDG: "<<GPartPdg->at(0)<<"\n"; 

         std::cout<<" Primaries PDG's:";
         for(int a = 0; a < PDGMother->size(); a++){
            if (PDGMother->at(a) != 0) continue;
            std::cout<<" "<<PDG->at(a);
         }

         std::cout<<" Tracks PDG's:";
         for(int a = 0; a < TrackMCindex->size(); a++){
            std::cout<<" "<<PDG->at(TrackMCindex->at(a));
         }

        std::cout<<" Tracks PDGMothers:";
        for(int a = 0; a < TrackMCindex->size(); a++){
            std::cout<<" "<<PDGMother->at(TrackMCindex->at(a));
         }

         std::cout<<"\n\n";
    }
    std::cout<<"count: "<<count<<" NCcount: "<<NCcount<<" antinumucount: "<<antinumucount<<" nuecount: "<<nuecount<<std::endl;
    std::cout<<"mu_proton_piplus_pi0: "<<mu_proton_piplus_pi0<<" mu_proton_piplus_2pi0: "<<mu_proton_piplus_2pi0<<" mu_proton_piplus_3pi0: "<<mu_proton_piplus_3pi0<<std::endl;
}

 