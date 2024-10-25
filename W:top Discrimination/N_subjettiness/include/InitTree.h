//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul 17 21:43:29 2020 by ROOT version 6.12/07
// from TTree trackTree/v1
// found on file: outputMiniAOD.root
//////////////////////////////////////////////////////////




#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

#include <string>

#include <iostream>
#include <iomanip>

#include <vector>
#include "math.h"
#include <numeric>
#include <string>
#include <fstream>
#include <algorithm>
#include <iterator>

#include <TStyle.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TCutG.h"
#include "TRandom3.h"

using namespace std;

#include "Pythia8/Pythia.h"
Pythia8::Pythia pythia;
Pythia8::ParticleData &particleData = pythia.particleData;


//fastjet
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/Nsubjettiness.hh" // In external code, this should be fastjet/contrib/Nsubjettiness.hh
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#include "fastjet/contrib/XConePlugin.hh"

using TMath::ATan;
using TMath::Exp;

#include "1d2d_constants.h"

class Gosling {
public :
   TFile          * fFile;
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   std::vector<std::string> fileList;

   

  
   vector<float>   *genJetEta;
   vector<float>   *genJetPt;
   vector<float>   *genJetPhi;
   vector<int>     *genJetChargedMultiplicity;
   vector<vector<int> > *genDau_chg;
   vector<vector<int> > *genDau_pid;
   vector<vector<float> > *genDau_pt;
   vector<vector<float> > *genDau_eta;
   vector<vector<float> > *genDau_phi;
   // List of branches
   TBranch        *b_genJetEta;   //!
   TBranch        *b_genJetPt;   //!
   TBranch        *b_genJetPhi;   //!
   TBranch        *b_genJetChargedMultiplicity;   //!
   TBranch        *b_genDau_chg;   //!
   TBranch        *b_genDau_pid;   //!
   TBranch        *b_genDau_pt;   //!
   TBranch        *b_genDau_eta;   //!
   TBranch        *b_genDau_phi;   


  

   Gosling(std::vector<std::string> _fileList);
   virtual ~Gosling();
   
   virtual Int_t    GetEntry(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
  
};

Gosling::Gosling(std::vector<std::string> _fileList) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
//   if (tree == 0) {
//      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("total_output_ak8.root");
//      if (!f || !f->IsOpen()) {
//         f = new TFile("total_output_ak8.root");
//      }
//      TDirectory * dir = (TDirectory*)f->Get("total_output_ak8.root:/analyzer");
//      dir->GetObject("trackTree",tree);
//
//   }
    fileList = _fileList;
}

Gosling::~Gosling()
{
   //if (!fChain) return;
   //delete fChain->GetCurrentFile();
}

Int_t Gosling::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

void Gosling::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   genJetEta=0;
   genJetPt=0;
   genJetPhi=0;
   genJetChargedMultiplicity=0;
   genDau_chg=0;
   genDau_pid=0;
   genDau_pt=0;
   genDau_eta=0;
   genDau_phi=0;
   // Set branch addresses and branch pointers
  if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);
   fChain->SetBranchAddress("genJetEta"                ,&genJetEta                ,&b_genJetEta                );
   fChain->SetBranchAddress("genJetPt"                 ,&genJetPt                 ,&b_genJetPt                 );
   fChain->SetBranchAddress("genJetPhi"                ,&genJetPhi                ,&b_genJetPhi                );
   fChain->SetBranchAddress("genJetChargedMultiplicity",&genJetChargedMultiplicity,&b_genJetChargedMultiplicity);
   fChain->SetBranchAddress("genDau_chg"               ,&genDau_chg               ,&b_genDau_chg               );
   fChain->SetBranchAddress("genDau_pid"               ,&genDau_pid               ,&b_genDau_pid               );
   fChain->SetBranchAddress("genDau_pt"                ,&genDau_pt                ,&b_genDau_pt                );
   fChain->SetBranchAddress("genDau_eta"               ,&genDau_eta               ,&b_genDau_eta               );
   fChain->SetBranchAddress("genDau_phi"               ,&genDau_phi               ,&b_genDau_phi               );
   
}















