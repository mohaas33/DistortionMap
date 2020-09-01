#define analyzeHits_cxx
#include "analyzeHits.h"
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>

#include <TStyle.h>
#include <TCanvas.h>


void GetPadCenter(double r, double phi, int &rsec, int &nsec, int &padNr, int &padNphi, double &rin, double &cntrr, double &cntrphi);

void GetPadCenter(double r, double phi, int &rsec, int &nsec, int &padNr, int &padNphi, double &rin, double &cntrr, double &cntrphi){
  //these parameters are taken from Feb 12 drawings of frames.
  double tpc_frame_side_gap=0.8;//mm //space between radial line and start of frame
  double tpc_frame_side_width=2.6;//mm //thickness of frame
  double tpc_margin=0.0;//mm // extra gap between edge of frame and start of GEM holes
  
  double tpc_frame_r3_outer=758.4;//mm inner edge of larger-r frame of r3
  double tpc_frame_r3_inner=583.5;//mm outer edge of smaller-r frame of r3
 
  double tpc_frame_r2_outer=574.9;//mm inner edge of larger-r frame of r3
  double tpc_frame_r2_inner=411.4;//mm outer edge of smaller-r frame of r3
 
  double tpc_frame_r1_outer=402.6;//mm inner edge of larger-r frame of r3
  double tpc_frame_r1_inner=221.0;//mm outer edge of smaller-r frame of r3
  
  cntrr=0;
  cntrphi=0;
  padNphi=0;
  padNr=0;
  rsec=0;
  rin=0;
    double pad_r=0;
    if (r>=tpc_frame_r1_inner+tpc_margin && r<=tpc_frame_r1_outer+tpc_margin){
        rsec = 1;
        pad_r = (tpc_frame_r1_outer - tpc_frame_r1_inner)/16;
        rin = tpc_frame_r1_inner;
    }
    if (r>=tpc_frame_r2_inner+tpc_margin && r<=tpc_frame_r2_outer+tpc_margin){
        rsec = 2;
        pad_r = (tpc_frame_r2_outer - tpc_frame_r2_inner)/16;
        rin = tpc_frame_r2_inner;
    }
    if (r>=tpc_frame_r3_inner+tpc_margin and r<=tpc_frame_r3_outer+tpc_margin){
        rsec = 3;
        pad_r = (tpc_frame_r3_outer - tpc_frame_r3_inner)/16;
        rin = tpc_frame_r3_inner;
    }
    padNr= floor((r - rin)/pad_r);
    cntrr = double(padNr)*pad_r+rin+pad_r/2;

    //if the coordinate is within gap+width of a sector boundary, return True:
    //note that this is not a line of constant radius, but a linear distance from a radius.

    //find the two spokes we're between:
    double sectorangle=(M_PI/6);
    int nsectors=phi/sectorangle;
    nsec=floor(nsectors)+1;
    //if(r>tpc_frame_r1_inner)cout<<sectorangle*r/5<<endl;
    double arc_length = (phi-(nsec)*sectorangle)*r;
    padNphi = floor(arc_length/5);
    cntrphi = (nsec)*sectorangle + (double(padNphi)+0.5)*5/cntrr;
 
    
    
    
    nsec++;


}

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
  TH3D *heCharge[nHist];
  TH2D *heChargePad[nHist];
  TH2D *heChargePadXY[nHist];
  for(int hi=0;hi<nHist;hi++){
     hCharge[hi]=new TH3D(Form("h_Charge_%d",hi),"SC (ions) per m^3;phi (rad);r (m);z (m)",nphi,0,6.28319,nr,rmin,rmax,nz,0,z_rdo);
     heCharge[hi]=new TH3D(Form("h_eCharge_%d",hi),"Charge (electrons) per m^3;phi (rad);r (m);z (m)",nphi,0,6.28319,nr,rmin,rmax,nz,0,z_rdo);
     heChargePad[hi]=new TH2D(Form("h_eChargePad_%d",hi),"Charge (electrons) per m^3;phi (rad);r (m);z (m)",nphi,0,6.28319,nr,rmin,rmax);
     heChargePadXY[hi]=new TH2D(Form("h_eChargePadXY_%d",hi),"Charge (electrons) per m^3;phi (rad);r (m);z (m)",1600*4,-0.8,0.8,1600*4,-0.8,0.8);
  }
   //734352 bunches to fill the TPC
   //int bX=734587;//first one after which TPC has been filled according to timestamps_50kHz.txt 
   int bX = beamXing;
   std::vector<int>::iterator it = std::find(beamXings.begin(), beamXings.end(), bX);
   int index; 
   index = std::distance(beamXings.begin(), it) + beamXingBias;
   if(f15kHz==1)index = std::distance(beamXings.begin(), it);
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
         if(f15kHz==1)bX = beamXings.at(index)+(hi+beamXingBias)*38860;
         double phi = hit_phi;
         double r = hit_r;
         double z_prim = hit_z-(bX-event_bunchXing)*106*vIon*ns;
         double z_ibf = 1.05-(bX-event_bunchXing)*106*vIon*ns;
         //cout<<"z_prim = "<<z_prim<<endl;
         //cout<<"bX-event_bunchXing = "<<bX<<"-"<<event_bunchXing<<endl;
         //cout<<"|----------------------------------------------|"<<endl;
         //evtN++;
         //if(evtN==1){
         //   cout<<z_prim <<"="<< hit_z<<"+("<<bX<<"-"<<event_bunchXing<<")*"<<106*vIon*ns <<endl;
         //   cout<<z_prim <<"="<< hit_z<<"+("<<bX-event_bunchXing<<")*"<<106*vIon*ns <<endl;
         //}
         double w_prim = hit_eion*Tpc_ElectronsPerGeV;
         double w_ibf = ibf_vol;
         double w_charge = amp_ele_vol*hit_eion*Tpc_ElectronsPerGeV;
         if(event_bunchXing<bX){
            hCharge[hi]->Fill(phi,r,z_prim,w_prim);
            if(isOnPlane){
               int rsec=0; // number of sector in r
               int nsec=0; // number of sector in phi
               int padNr = 0; // number of pad in r              
               int padNphi = 0; // number of pad in phi    
               double cntrr = 0; // pad center in r
               double cntrphi = 0; // pad center in ph
               double rin=0;
               GetPadCenter(r*1e3, phi, rsec, nsec, padNr, padNphi, rin, cntrr, cntrphi);              
               heCharge[hi]->Fill(phi,r,hit_z,w_charge);
               heChargePad[hi]->Fill(cntrphi,cntrr/1e3,w_charge);
               double x = cntrr/1e3*sin(cntrphi);
               double y = cntrr/1e3*cos(cntrphi);
               heChargePadXY[hi]->Fill(x,y,w_charge);

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
   cout<<outputFileName<<endl;
   TString s_search = "/Files/";
   int s_start = outputFileName.Index(s_search);
   int s_end = outputFileName.Length();
   int s_len = s_search.Length();
   TString outputFileName_e = outputFileName(0,s_start) + "/Files/e_"+outputFileName(s_start+s_len,s_end);
   cout<<outputFileName_e<<endl;

   //std::string subject = outputFileName;
   //std::string 
   //std::string replace = "/Files/e_";
   //size_t pos = 0;
   //while((pos = subject.find(search, pos)) != std::string::npos) {
   //      subject.replace(pos, search.length(), replace);
   //      pos += replace.length();
   //}

   TFile e_outputFile (outputFileName_e ,"RECREATE");
   for(int hi=0;hi<nHist;hi++){
      heCharge[hi]->Write();
      heChargePad[hi]->Write();
      heChargePadXY[hi]->Write();
   }
   e_outputFile.Close();


}
