//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Aug 26 03:28:48 2020 by ROOT version 6.16/00
// from TTree hTree/tpc hit tree for ionization
// found on file: ../Files/slim_G4Hits_sHijing_0-12fm_000000_001000.root
//////////////////////////////////////////////////////////

#ifndef analyzeHits_h
#define analyzeHits_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.


class analyzeHits {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           isOnPlane;
   Float_t         hit_z;
   Float_t         hit_r;
   Float_t         hit_phi;
   Float_t         hit_eion;
   Float_t         ibf_vol;
   Float_t         amp_ele_vol;
   Float_t         event_timestamp;
   Float_t         event_bunchXing;

   //TTree *_tree = 0;
   // List of branches
   TBranch        *b_isOnPlane;   //!
   TBranch        *b_hit_z;   //!
   TBranch        *b_hit_r;   //!
   TBranch        *b_hit_phi;   //!
   TBranch        *b_hit_eion;   //!
   TBranch        *b_ibf_vol;   //!
   TBranch        *b_amp_ele_vol;   //!
   TBranch        *b_event_timestamp;   //!
   TBranch        *b_event_bunchXing;   //!

   TString outputFileName = "outputFile.root";
   std::vector<int> beamXings;
   int beamXing = 734587;
   int beamXingBias = 0;
   int f15kHz = 0;
   int fPrim = 0;
   
   analyzeHits(TString FileName="/sphenix/user/shulga/Work/IBF/DistortionMap/Files/slim_G4Hits_sHijing_0-12fm_000000_001000.root",  TTree *_tree=0);
   virtual ~analyzeHits();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     SetOutputFileName(TString name);
   virtual void     SetBeamXings(std::vector<int> beamXs);
   virtual void     SetBeamXing(int beamX);
   virtual void     SetBeamXingBias(int beamXbias);
   //virtual void     SetTree(TTree *tree);
   virtual void     RunLaser15kHz(int fkHz = 1);
   virtual void     RunPrim(int fprim = 1);

};

#endif

#ifdef analyzeHits_cxx
analyzeHits::analyzeHits(TString FileName, TTree *_tree) : fChain(0) 
{

   // if parameter _tree is not specified (or zero), connect the file
   // used to generate this class and read the Tree.
   //TTree *_tree;

   if (_tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(FileName);
      if (!f || !f->IsOpen()) {
       f = new TFile(FileName);
      }
      f->GetObject("hTree",_tree);
   }
   Init(_tree);
}

analyzeHits::~analyzeHits()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t analyzeHits::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t analyzeHits::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void analyzeHits::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("isOnPlane", &isOnPlane, &b_isOnPlane);
   fChain->SetBranchAddress("hit_z", &hit_z, &b_hit_z);
   fChain->SetBranchAddress("hit_r", &hit_r, &b_hit_r);
   fChain->SetBranchAddress("hit_phi", &hit_phi, &b_hit_phi);
   fChain->SetBranchAddress("hit_eion", &hit_eion, &b_hit_eion);
   fChain->SetBranchAddress("ibf_vol", &ibf_vol, &b_ibf_vol);
   fChain->SetBranchAddress("amp_ele_vol", &amp_ele_vol, &b_amp_ele_vol);
   fChain->SetBranchAddress("event_timestamp", &event_timestamp, &b_event_timestamp);
   fChain->SetBranchAddress("event_bunchXing", &event_bunchXing, &b_event_bunchXing);
   Notify();
}

Bool_t analyzeHits::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void analyzeHits::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
void analyzeHits::SetOutputFileName(TString name){
   //Set up Output File Name
   outputFileName = name;
}
void analyzeHits::SetBeamXings(std::vector<int> beamXs){
   beamXings = beamXs;
}
void analyzeHits::SetBeamXing(int beamX){
   beamXing = beamX;
}
void analyzeHits::SetBeamXingBias(int beamXbias){
   beamXingBias = beamXbias;
}
//void analyzeHits::SetTree(TTree *tree){
//   _tree = tree;
//}
void analyzeHits::RunLaser15kHz(int fkHz){
   f15kHz = fkHz;
}
void analyzeHits::RunPrim(int fprim = 1){
   fPrim = fprim;
}

Int_t analyzeHits::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef analyzeHits_cxx
