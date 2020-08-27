#define analyzeHits_cxx
#include "analyzeHits.h"
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>

#include <TStyle.h>
#include <TCanvas.h>


void analyzeHits::Loop()
{
//   In a ROOT session, you can do:
//      root> .L analyzeHits.C
//      root> analyzeHits t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch





  double f=0.5;//for now, just pick the middle of the hit.  Do better later.
  double ns=1e-9,us=1e-6,ms=1e-3,s=1;
  double um=1e-6, mm=1e-3, cm=1e-2,m=1; //changed to make 'cm' 1.0, for convenience.
  double Hz=1, kHz=1e3, MHz=1e6;
  double V=1;
  //used two ways:  1) to apply units to variables when defined
  //                2) to divide by certain units so that those variables are expressed in those units.

  double ionMobility=3.37*cm*cm/V/s;
  double vIon=ionMobility*400*V/cm;
  //float vIon=16.0*um/us;
  //float mbRate=_freqKhz*kHz;
  //float xingRate = 9.383*MHz;
  //float mean = mbRate/xingRate;
  double z_rdo=105.5*cm;
  double rmin=20*cm;
  double rmax=78*cm;

  double Ne_dEdx = 1.56/cm;   // keV/cm                                                                                                            
  double CF4_dEdx = 7.00/cm;  // keV/cm                                                                                                            
  double Ne_NTotal = 43/cm;    // Number/cm                                                                                                        
  double CF4_NTotal = 100/cm;  // Number/cm                                                                                                        
  double Tpc_NTot = 0.90 * Ne_NTotal + 0.10 * CF4_NTotal;
  double Tpc_dEdx = 0.90 * Ne_dEdx + 0.10 * CF4_dEdx;
  double Tpc_ElectronsPerKeV = Tpc_NTot / Tpc_dEdx;
  double Tpc_ElectronsPerGeV = Tpc_NTot / Tpc_dEdx*1e6; //electrons per gev.                                                                       


  int nr=159;
  int nphi=360;
  int nz=62*2;
 

  char fname[100];
  const int nHist = 10;
  TH3D *hCharge[nHist];
  for(int hi=0;hi<nHist;hi++){
     hCharge[hi]=new TH3D(Form("h_Charge_%d",hi),"SC (ions) per m^3;phi (rad);r (m);z (m)",nphi,0,6.28319,nr,rmin,rmax,nz,0,z_rdo);
  }
   //734352 bunches to fill the TPC
   //int bX=734587;//first one after which TPC has been filled according to timestamps_50kHz.txt 
   int bX = beamXing;
   std::vector<int>::iterator it = std::find(beamXings.begin(), beamXings.end(), bX);
   int index = std::distance(beamXings.begin(), it) + beamXingBias;

   //int evtN=0;
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //cout<<hit_z<<endl;
      // if (Cut(ientry) < 0) continue;
      //cout<<"index="<<index<<endl;
      //beamXings
      for(int hi=0;hi<nHist;hi++){
         bX = beamXings.at(index+hi);
         double phi = hit_phi;
         double r = hit_r;
         double z_prim = hit_z+(bX-event_bunchXing)*106*vIon*ns;
         double z_ibf = 1.05-(bX-event_bunchXing)*106*vIon*ns;
         //evtN++;
         //if(evtN==1){
         //   cout<<z_prim <<"="<< hit_z<<"+("<<bX<<"-"<<event_bunchXing<<")*"<<106*vIon*ns <<endl;
         //   cout<<z_prim <<"="<< hit_z<<"+("<<bX-event_bunchXing<<")*"<<106*vIon*ns <<endl;
         //}
         double w_prim = hit_eion*Tpc_ElectronsPerGeV;
         double w_ibf = ibf_vol;
         if(event_bunchXing<bX){
            hCharge[hi]->Fill(phi,r,z_prim,w_prim);
            if(isOnPlane){
               hCharge[hi]->Fill(phi,r,z_ibf,w_ibf);
            }
         }
      }
   }
   
   TFile outputFile (outputFileName,"RECREATE");
   for(int hi=0;hi<nHist;hi++){
      hCharge[hi]->Write();
   }
   outputFile.Close();


}
