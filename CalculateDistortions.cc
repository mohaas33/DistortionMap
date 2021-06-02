#include "CalculateDistortions.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/SubsysReco.h>           // for SubsysReco

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>


#include <phool/PHCompositeNode.h>

#include <phool/getClass.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TNtuple.h>
#include <TTree.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <typeinfo>
#include <cstdlib>
#include <map>

using namespace std;


bool IsOverFrame(double r, double phi);

bool IsOverFrame(double r, double phi){
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
 
  //double tpc_sec0_phi=0.0;//get_double_param("tpc_sec0_phi");

  //if the coordinate is in the radial spaces of the frames, return true:
  if (r<tpc_frame_r1_inner+tpc_margin)
    return true;
  if (r>tpc_frame_r1_outer-tpc_margin  && r<tpc_frame_r2_inner+tpc_margin)
    return true;
  if (r>tpc_frame_r2_outer-tpc_margin  && r<tpc_frame_r3_inner+tpc_margin)
    return true;
  if (r>tpc_frame_r3_outer-tpc_margin)
    return true;

  //if the coordinate is within gap+width of a sector boundary, return true:
  //note that this is not a line of constant radius, but a linear distance from a radius.

  //find the two spokes we're between:
  double pi = 2 * acos(0.0);

  float sectorangle=(pi/6);
  float nsectors=phi/sectorangle;
  int nsec=floor(nsectors);
  float reduced_phi=phi-nsec*sectorangle; //between zero and sixty degrees.
  float dist_to_previous=r*sin(reduced_phi);
  float dist_to_next=r*sin(sectorangle-reduced_phi);
  if (dist_to_previous<tpc_frame_side_gap+tpc_frame_side_width+tpc_margin)
    return true;
  if (dist_to_next<tpc_frame_side_gap+tpc_frame_side_width+tpc_margin)
    return true;
  
  return false;
}

//____________________________________________________________________________..
CalculateDistortions::CalculateDistortions(const std::string &name, const std::string &filename):
  SubsysReco(name)
 , hm(nullptr)
 , _filename(filename)
 , outfile(nullptr)
 ,_ampGain(2e3)
 ,_ampIBFfrac(0.02)
 ,_collSyst(0)

{
  cout << "CalculateDistortions::CalculateDistortions(const std::string &name) Calling ctor" << endl;
}

//____________________________________________________________________________..
CalculateDistortions::~CalculateDistortions()
{
  //cout << "CalculateDistortions::~CalculateDistortions() Calling dtor" << endl;
  // delete whatever is created 
  delete hm;
}

//____________________________________________________________________________..
int CalculateDistortions::Init(PHCompositeNode *topNode)
{
  double cm=1e-2; //changed to make 'm' 1.0, for convenience.

  int nr=159;
  int nphi=360;
  int nz=62*2;
  double z_rdo=105.5*cm;
  double rmin=20*cm;
  double rmax=78*cm;
  //cout << "CalculateDistortions::Init(PHCompositeNode *topNode) Initializing" << endl;
  hm = new Fun4AllHistoManager("HITHIST");

  _h_SC_prim = new TH3F("_h_SC_prim","_h_SC_prim;#phi, [rad];R, [m];Z, [m]"  ,nphi,0,6.28319,nr,rmin,rmax,2*nz,-z_rdo,z_rdo);
  _h_SC_ibf  = new TH3F("_h_SC_ibf" ,"_h_SC_ibf;#phi, [rad];R, [m];Z, [m]"   ,nphi,0,6.28319,nr,rmin,rmax,2*nz,-z_rdo,z_rdo);
  _h_hits  = new TH1F("_h_hits" ,"_h_hits;N, [hit]"   ,4000,0,1e6);
  hm->registerHisto(_h_SC_prim);
  hm->registerHisto(_h_SC_ibf );
  hm->registerHisto(_h_hits );

  outfile = new TFile(_filename.c_str(), "RECREATE");
  _event_timestamp = 0;
  _hit_eion  = 0;
  _hit_r   = 0;
  _hit_phi = 0;
  _hit_z = 0;
  _ibf_vol     = 0;
  _amp_ele_vol = 0;
  if(_fSliming==1){
    _rawHits=new TTree("hTree","tpc hit tree for ionization");
    _rawHits->Branch("isOnPlane",&_isOnPlane);
    _rawHits->Branch("hit_z",&_hit_z);
    _rawHits->Branch("hit_r",&_hit_r);
    _rawHits->Branch("hit_phi",&_hit_phi);
    _rawHits->Branch("hit_eion",&_hit_eion);
    _rawHits->Branch("ibf_vol"    ,&_ibf_vol    );
    _rawHits->Branch("amp_ele_vol",&_amp_ele_vol);

    _rawHits->Branch("event_timestamp",&_event_timestamp);
    _rawHits->Branch("event_bunchXing",&_event_bunchXing);
  }
  return 0;
}

//____________________________________________________________________________..
int CalculateDistortions::InitRun(PHCompositeNode *topNode)
{
  //cout << "CalculateDistortions::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << endl;
  std::string line;
  //AA collisions timestamps
  std::string txt_file = "/sphenix/user/shulga/Work/IBF/DistortionMap/timestamps_50kHz.txt";
  int start_line = 3;
  if(_collSyst==1){
    //pp collisions timestamps
    txt_file = "/phenix/u/hpereira/sphenix/work/g4simulations/timestamps_3MHz.txt";
    //txt_file = "/sphenix/user/shulga/Work/IBF/DistortionMap/timestamps_50kHz.txt";
    start_line = 2;
  }
  ifstream InputFile (txt_file);
  if (InputFile.is_open()){
    int n_line=0;
    while ( getline (InputFile,line) )
    {
      n_line++;
      if(n_line>start_line){
        std::istringstream is( line );
        double n[2] = {0,0};
        int i = 0;
        while( is >> n[i] ) {    
            i++;    
        }
        _timestamps[n[0]]=n[1];
        if(n_line<10){
          cout<<n[1]<<endl;
        }
        _keys.push_back(int(n[0]));
      }
    }
    InputFile.close();
  }

  else cout << "Unable to open file:"<<txt_file<<endl; 

  TFile *MapsFile; 
  if(_fUseIBFMap){
    MapsFile = new TFile("/sphenix/user/shulga/Work/IBF/DistortionMap/IBF_Map.root","READ");
    if ( MapsFile->IsOpen() ) printf("File opened successfully\n");
    _h_modules_anode       = (TH2F*)MapsFile ->Get("h_modules_anode")      ->Clone("_h_modules_anode");
    _h_modules_measuredibf = (TH2F*)MapsFile ->Get("h_modules_measuredibf")->Clone("_h_modules_measuredibf");
  }

  _mbRate   = _freqKhz*kHz;
  _xingRate = 9.383*MHz;
  _mean     = mbRate/xingRate;

  return 0;
}

//____________________________________________________________________________..
int CalculateDistortions::process_event(PHCompositeNode *topNode)
{
    
  double bX = _beamxing;
  double z_bias_avg = 0;
  if (_fAvg==1){ 
    z_bias_avg=1.05*(float) rand()/RAND_MAX;
  }
  int bemxingsInFile = _keys.size();
  if (_evtstart>= bemxingsInFile) _evtstart=0;
  int key = _keys.at(_evtstart);
  _event_timestamp = (float)_timestamps[key]*ns;//units in seconds
  _event_bunchXing = key;
  if(_evtstart%100==0) cout<<"_evtstart = "<<_evtstart<<endl;
  _evtstart++;

  ostringstream nodename;
  set<std::string>::const_iterator iter;
  nodename << "G4HIT_TPC";

  PHG4HitContainer *hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());
  int n_hits = 0;
  if (hits){
    PHG4HitContainer::ConstRange hit_range = hits->getHits();
      for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++){
        n_hits++;
        int f_fill_prim=1;
        int f_fill_ibf=1;

        float hit_x0 = hit_iter->second->get_x(0);
        float hit_y0 = hit_iter->second->get_y(0);
        float hit_z0 = hit_iter->second->get_z(0);
        float hit_x1 = hit_iter->second->get_x(1);
        float hit_y1 = hit_iter->second->get_y(1);
        float hit_z1 = hit_iter->second->get_z(1);

        float hit_eion = hit_iter->second->get_eion();
        float N_electrons=hit_eion*Tpc_ElectronsPerGeV;
        float x = (hit_x0 + f * (hit_x1 - hit_x0))*cm;
        float y = (hit_y0 + f * (hit_y1 - hit_y0))*cm;
        float z = (hit_z0 + f * (hit_z1 - hit_z0))*cm;
        
        float r=sqrt(x*x+y*y);
        float phi=atan2(x,y);
        if (phi<0) phi+=2*pi;

        //Reading IBF and Gain weights according to X-Y position
        float w_ibf = 1.;
        float w_gain = 1.;
        if(_fUseIBFMap){
          int bin_x = _h_modules_anode ->GetXaxis()->FindBin(x/mm);
          int bin_y = _h_modules_anode ->GetYaxis()->FindBin(y/mm);
          w_ibf = _h_modules_measuredibf->GetBinContent(bin_x,bin_y);
          w_gain = _h_modules_anode->GetBinContent(bin_x,bin_y);
          }
        float ionsPerEle=w_gain*_ampGain*w_ibf*_ampIBFfrac;

        _isOnPlane = 0;
        if(!IsOverFrame(r/mm,phi)){
          _isOnPlane = 1;
        }
        _hit_z = z;
        _hit_r = r;
        _hit_phi = phi;
        _hit_eion = hit_eion;
        _ibf_vol = N_electrons*ionsPerEle;
        _amp_ele_vol = w_gain*_ampGain;
  	    if(_fSliming==1)_rawHits->Fill();
        double z_prim = -1*1e10;
        double z_ibf =  -1*1e10;

        if(_hit_z>=0){
          if(_fAvg==1){
             z_prim = _hit_z - z_bias_avg;
             z_ibf  = 1.05  - z_bias_avg;
          }else{
            z_prim = _hit_z-(bX-_event_bunchXing)*106*vIon*ns;
            z_ibf = 1.05-(bX-_event_bunchXing)*106*vIon*ns;
          }
          if(z_prim<=0 ){
            f_fill_prim=0;
          }
          if( z_ibf<=0){
            f_fill_ibf=0;
          }
        }
        if(_hit_z<0){
           if(_fAvg==1){
              z_prim = _hit_z + z_bias_avg;
              z_ibf  = -1.05  + z_bias_avg;
           }else{
             z_prim = _hit_z+(bX-_event_bunchXing)*106*vIon*ns;
             z_ibf = -1.05+(bX-_event_bunchXing)*106*vIon*ns;
           }
           if(z_prim>=0 ){
             f_fill_prim=0;
           }
           if( z_ibf>=0){
             f_fill_ibf=0;
           }
        }
        //if(n_hits<5)cout<<"z_bias_avg="<<z_bias_avg<<" IBF ="<<z_ibf<<" prim ="<<z_prim<<endl;

        double w_prim = _hit_eion*Tpc_ElectronsPerGeV;
        if(f_fill_prim==1)_h_SC_prim ->Fill(_hit_phi,_hit_r,z_prim,w_prim);
        if(_isOnPlane && f_fill_ibf==1)_h_SC_ibf  ->Fill(_hit_phi,_hit_r,z_ibf,_ibf_vol);
      }

  }else{
    if(_fSliming==1)_rawHits->Fill();
  }
  _h_hits->Fill(n_hits);
  return 0;
}

//____________________________________________________________________________..
//int CalculateDistortions::ResetEvent(PHCompositeNode *topNode)
//{
//  cout << "CalculateDistortions::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << endl;
//  return Fun4AllReturnCodes::EVENT_OK;
//}

//____________________________________________________________________________..
int CalculateDistortions::EndRun(const int runnumber)
{
  cout << "CalculateDistortions::EndRun(const int runnumber) Ending Run for Run " << runnumber << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CalculateDistortions::End(PHCompositeNode *topNode)
{
  if(_fSliming==1){
    outfile->cd();
    outfile->Write();
    outfile->Close();
    delete outfile;
    _h_SC_prim ->Sumw2( false );
    _h_SC_ibf  ->Sumw2( false );
    _h_hits    ->Sumw2( false );
    hm->dumpHistos(_filename, "UPDATE");
  }else{
    hm->dumpHistos(_filename, "RECREATE");
  }

  return 0;
}

//____________________________________________________________________________..
int CalculateDistortions::Reset(PHCompositeNode *topNode)
{
 cout << "CalculateDistortions::Reset(PHCompositeNode *topNode) being Reset" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void CalculateDistortions::Print(const std::string &what) const
{
  cout << "CalculateDistortions::Print(const std::string &what) const Printing info for " << what << endl;
}

void CalculateDistortions::SetFrequency(int freq){
  _freqKhz = freq;
  cout<<"Frequency is set to: "<<_freqKhz<<" kHz"<<endl;
}
void CalculateDistortions::SetBeamXing(int newBeamXing){
  _beamxing = newBeamXing;
  cout<<"Initial BeamXing is set to: "<<newBeamXing<<endl;

}
void CalculateDistortions::SetEvtStart(int newEvtStart){
  _evtstart = newEvtStart;
  cout<<"Starting event is set to: "<<newEvtStart<<endl;

}

void CalculateDistortions::SetUseIBFMap(bool useIBFMap){
  _fUseIBFMap = useIBFMap;
  cout<<"Using IBF and Gain Maps"<<endl;
}
void CalculateDistortions::SetGain(float ampGain){
  _ampGain = ampGain;
  cout<<"Gain is set to: "<<_ampGain<<endl;
}
void CalculateDistortions::SetIBF(float ampIBFfrac){
  _ampIBFfrac = ampIBFfrac;
  cout<<"IBF is set to: "<<_ampIBFfrac<<endl;
}

void CalculateDistortions::SetCollSyst(int coll_syst){
  _collSyst = coll_syst;
  std::string s_syst[2] = {"AA","pp"};
  cout<<"Collision system is set to: "<<s_syst[_collSyst]<<endl;

}

void CalculateDistortions::SetAvg(int fAvg){
  _fAvg = fAvg;
  std::string s_avg[2] = {"OFF","ON"};
  cout<<"Averaging is set to: "<<s_avg[_fAvg]<<endl;

}
void CalculateDistortions::UseSliming(int fSliming){
  _fSliming = fSliming;
  std::string s_sliming[2] = {"OFF","ON"};
  cout<<"Sliming is set to: "<<s_sliming[_fSliming]<<endl;

}
