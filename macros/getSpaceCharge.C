#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TH3.h>

#include <stdio.h>

#include <string>

//TChain *chanT = NULL;

void getSpaceCharge(){
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
  const int nFiles = 1;
  TH3D *hCharge[nFiles];
  TH3D *hCharge_tot=new TH3D("hCharge_tot","SC (ions) per m^3;phi (rad);r (m);z (m)",nphi,0,6.28319,nr,rmin,rmax,nz,0,z_rdo);
  
  for (int i=0;i<nFiles;i++){
  
    TChain *chanT = new TChain("hTree");
    TH3D *hCharge_tmp=new TH3D("hCharge_tmp","SC (ions) per m^3;phi (rad);r (m);z (m)",nphi,0,6.28319,nr,rmin,rmax,nz,0,z_rdo);

    int eventsInFileStart = i*1000;
    int eventsInFileEnd = (i+1)*1000;

    sprintf(fname, "/sphenix/user/shulga/Work/IBF/DistortionMap/Files/slim_G4Hits_sHijing_0-12fm_%06d_%06d.root",eventsInFileStart,eventsInFileEnd);
    //sprintf(fname, "./slim_G4Hits_sHijing_0-12fm_000000_001000.root");
    chanT->Add(fname);
    //734352 bunches to fill the TPC
    int bX=734587;//first one after which TPC has bee filled according to timestamps_50kHz.txt 
    chanT->Draw(Form("hit_z+(%d-event_bunchXing)*106*%.10f:hit_r:hit_phi>> hCharge_tmp",bX,vIon*ns),Form("(hit_eion*%f)*(event_bunchXing<%d)",Tpc_ElectronsPerGeV,bX)); //A map for event right after bunch xing bX has occurred
    chanT->Draw(Form("1.05-(%d-event_bunchXing)*106*%.10f:hit_r:hit_phi>>+hCharge_tmp",bX,vIon*ns),Form("(ibf_vol)*(event_bunchXing<%d && isOnPlane)",bX)); //A map of IBF for event right after bunch xing B has occurred
    hCharge[i]=(TH3D*)hCharge_tmp->Clone("h_Charge_tot");
    if(i==0){
      hCharge_tot=(TH3D*)hCharge_tmp->Clone(Form("h_Charge_%06d_%06d",eventsInFileStart,eventsInFileEnd));
    }else{
      hCharge_tot->Add(hCharge_tmp);
    }
  }
  TFile outputFile ("outputFile.root","RECREATE");
  hCharge_tot->Write();
  outputFile.Close();

  TCanvas *c = new TCanvas();
  hCharge[0]->Draw();
  c->Print("./charges.png");
  
  TCanvas *c_xy = new TCanvas();
  c_xy->GetPad(0)->SetLogz();
  c_xy->SetRightMargin(0.15);
  hCharge[0]->Project3D("xy")->Draw("colz");
  c_xy->Print("./charges_xy.png");
  
  TCanvas *c_z = new TCanvas();
  hCharge[0]->Project3D("z")->Draw("HIST");
  c_z->Print("./charges_z.png");
}
