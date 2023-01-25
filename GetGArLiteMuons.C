#include "TROOT.h"
#include "TRandom.h"
#include "iostream"
#include "TMath.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TChain.h"

void GetGArLiteMuons()
{
  //freopen("output.txt","w",stdout);
  gRandom->SetSeed(0);
  gSystem->RedirectOutput("muons.txt","w");

  TChain* tree_source = new TChain("/anatree/GArAnaTree");    
  tree_source->Add("/home/federico/Documents/Universita/Federico_2020-2021/OxfordCode/Kalman_Garlite/MCgarlite/6planes/muon_test_ana10k*");
  std::vector<int>     *PDG = 0;
  std::vector<float>   *MCPStartX = 0;
  std::vector<float>   *MCPStartY = 0;
  std::vector<float>   *MCPStartZ = 0;
  std::vector<float>   *MCPStartPX = 0;
  std::vector<float>   *MCPStartPY = 0;
  std::vector<float>   *MCPStartPZ = 0;
  TBranch        *b_PDG;   //!
  TBranch        *b_MCPStartX;   //!
  TBranch        *b_MCPStartY;   //!
  TBranch        *b_MCPStartZ;   //!
  TBranch        *b_MCPStartPX;   //!
  TBranch        *b_MCPStartPY;   //!
  TBranch        *b_MCPStartPZ;   //!
  tree_source->SetBranchAddress("PDG", &PDG, &b_PDG);
  tree_source->SetBranchAddress("MCPStartX", &MCPStartX, &b_MCPStartX);
  tree_source->SetBranchAddress("MCPStartY", &MCPStartY, &b_MCPStartY);
  tree_source->SetBranchAddress("MCPStartZ", &MCPStartZ, &b_MCPStartZ);
  tree_source->SetBranchAddress("MCPStartPX", &MCPStartPX, &b_MCPStartPX);
  tree_source->SetBranchAddress("MCPStartPY", &MCPStartPY, &b_MCPStartPY);
  tree_source->SetBranchAddress("MCPStartPZ", &MCPStartPZ, &b_MCPStartPZ);

  for (int i=0; i<tree_source->GetEntries(); ++i)
    {
      tree_source->GetEntry(i);
      double px=MCPStartPX->at(0);
      double py=MCPStartPY->at(0);
      double pz=MCPStartPZ->at(0);
      double p = TMath::Sqrt(px*px+py*py+pz*pz);

      double pdg = PDG->at(0);
      float x = MCPStartX->at(0);
      float y = MCPStartY->at(0);
      float z = 1486-276;
        
      float m = 0.10566;
      float e = TMath::Sqrt(p*p + m*m);
      std::cout << i << " " <<  1 << std::endl;
      std::cout << "1 "<<pdg<<" 0 0 0 0 " << px << " " << py << " " << pz << " " << e << " "  << m << " " << x << " " << y << " " << z << " " << 0 << std::endl; 
    }
  gApplication->Terminate(0);
  
}