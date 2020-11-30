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
 //,hCharge(nullptr)
 //,hChargePlane(nullptr)
 //,hCharge_v1(nullptr)
 //,hCharge_RPhi(nullptr)
 ,_ampGain(2e3)
 ,_ampIBFfrac(0.02)

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
  //cout << "CalculateDistortions::Init(PHCompositeNode *topNode) Initializing" << endl;
  hm = new Fun4AllHistoManager("HITHIST");
  outfile = new TFile(_filename.c_str(), "RECREATE");

  //hCharge=new TH3D("sphenix_minbias_charge","SC (ions) ;phi (rad);r (m);z (m)",nphi,0,6.28319,nr,rmin,rmax,nz,0,z_rdo);
  //hm->registerHisto(hCharge);
  //hChargePlane=new TH3D("sphenix_minbias_charge_plane","SC (ions) ;phi (rad);r (m);z (m)",nphi,0,6.28319,nr,rmin,rmax,nz,0,z_rdo);
  //hm->registerHisto(hChargePlane);
  //hCharge_v1=new TH3D("sphenix_minbias_charge_v1","SC (ions) ;phi (rad);r (m);z (m)",nphi,0,6.28319,nr,rmin,rmax,nz,0,z_rdo);
  //hm->registerHisto(hCharge_v1);
  //hCharge_RPhi=new TH2D("sphenix_minbias_charge_RPhi","SC (ions) ;phi (rad);r (m)",nphi,0,6.28319,nr,rmin,rmax);
  //hm->registerHisto(hCharge_RPhi);
  //TH1 *h1 = new TH1F("edep1GeV", "edep 0-1GeV", 1000, 0, 1);
  //eloss.push_back(h1);
  //h1 = new TH1F("edep100GeV", "edep 0-100GeV", 1000, 0, 100);
  //eloss.push_back(h1);
  //return Fun4AllReturnCodes::EVENT_OK;
  //_beamxing = 0;
  _event_timestamp = 0;
  _hit_eion  = 0;
  _hit_r   = 0;
  _hit_phi = 0;
  _hit_z = 0;
  _ibf_vol     = 0;
  _amp_ele_vol = 0;
  _rawHits=new TTree("hTree","tpc hit tree for ionization");
  //_rawHits->Branch("x",&x);
  //_rawHits->Branch("y",&y);
  //_rawHits->Branch("zorig",&z);
  //_rawHits->Branch("zibf",&zibf);
  _rawHits->Branch("isOnPlane",&_isOnPlane);
  _rawHits->Branch("hit_z",&_hit_z);
  _rawHits->Branch("hit_r",&_hit_r);
  _rawHits->Branch("hit_phi",&_hit_phi);
  _rawHits->Branch("hit_eion",&_hit_eion);
  _rawHits->Branch("ibf_vol"    ,&_ibf_vol    );
  _rawHits->Branch("amp_ele_vol",&_amp_ele_vol);

  _rawHits->Branch("event_timestamp",&_event_timestamp);
  _rawHits->Branch("event_bunchXing",&_event_bunchXing);
  //_rawHits->Branch("ne",&ne);

  return 0;
}

//____________________________________________________________________________..
int CalculateDistortions::InitRun(PHCompositeNode *topNode)
{
  //cout << "CalculateDistortions::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << endl;
  string line;
  string txt_file = "/sphenix/user/shulga/Work/IBF/DistortionMap/timestamps_50kHz.txt";
  ifstream InputFile (txt_file);
  //std::map<int,int> timestamps;
  if (InputFile.is_open()){
    int n_line=0;
    while ( getline (InputFile,line) )
    {
        n_line++;
      //cout << line << '\n';
      if(n_line>3){
        std::istringstream is( line );
        double n[2] = {0,0};
        int i = 0;
        while( is >> n[i] ) {    
            i++;    
        }
        _timestamps[n[0]]=n[1];
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
  int bemxingsInFile = _keys.size();
  if (_beamxing>= bemxingsInFile) _beamxing=0;
  //int n_key = nBeams-_beamxing-1;
  int key = _keys.at(_beamxing);
  _event_timestamp = (float)_timestamps[key]*ns;//units in seconds
  _event_bunchXing = key;
  if(_beamxing%100==0) cout<<"_beamxing = "<<_beamxing<<endl;
  _beamxing++;

  ostringstream nodename;
  set<string>::const_iterator iter;
  //vector<TH1 *>::const_iterator eiter;
  nodename << "G4HIT_TPC";

  PHG4HitContainer *hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());
  if (hits){
    PHG4HitContainer::ConstRange hit_range = hits->getHits();
   //for ( const auto &event_timestamp : _timestamps ) {
      
      //float t0 = _event_timestamp;//(beamxing/xingRate); 
      //float driftedZ=t0*vIon;//drift position in local units

      //double esum = 0;
      for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++){
        //esum += hit_iter->second->get_edep();
        //int hit_layer = hit_iter->second->get_layer();
        //int hit_scint_id = hit_iter->second->get_scint_id();
        float hit_x0 = hit_iter->second->get_x(0);
        float hit_y0 = hit_iter->second->get_y(0);
        float hit_z0 = hit_iter->second->get_z(0);
        float hit_x1 = hit_iter->second->get_x(1);
        float hit_y1 = hit_iter->second->get_y(1);
        float hit_z1 = hit_iter->second->get_z(1);
        //double Edep = hit_iter->second->get_edep();
        float hit_eion = hit_iter->second->get_eion();
        float N_electrons=hit_eion*Tpc_ElectronsPerGeV;
        float x = (hit_x0 + f * (hit_x1 - hit_x0))*cm;
        float y = (hit_y0 + f * (hit_y1 - hit_y0))*cm;
        float z = (hit_z0 + f * (hit_z1 - hit_z0))*cm;
        
        //if (z<0) continue;
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


        //float zprim=z-driftedZ;
        //float zibf=z_rdo-driftedZ;
        //zprim=1;
        //cout<<_freqKhz<<" "<<r<<phi<<zprim<<zibf;
        //int bin=hCharge->GetYaxis()->FindBin(r);
        //float hr=hCharge->GetYaxis()->GetBinLowEdge(bin);
        //float vol=(hzstep*hphistep*(hr+hrstep*0.5)*hrstep);

        //hCharge->Fill(phi,r,zprim,N_electrons/vol); //primary ion, drifted by t0, in cm
        //hCharge_RPhi->Fill(phi,r,N_electrons/vol); //primary ion, drifted by t0, in cm
        _isOnPlane = 0;
        if(!IsOverFrame(r/mm,phi)){
          //hCharge->Fill(phi,r,zibf,N_electrons*ionsPerEle/vol); //amp ion, drifted by t0, in cm
          //hCharge_RPhi->Fill(phi,r,N_electrons*ionsPerEle/vol); ///amp ion, drifted by t0, in cm
          //hChargePlane->Fill(phi,r,zprim,w_gain*_ampGain/vol);
          //if (saveTree){
          _isOnPlane = 1;
          //}
        }
        _hit_z = z;
        _hit_r = r;
        _hit_phi = phi;
        _hit_eion = hit_eion;
        _ibf_vol = N_electrons*ionsPerEle;
        _amp_ele_vol = w_gain*_ampGain;
        //cout<<"hit_eion="<<_hit_eion<<"vs"<<hit_eion<<endl;
  	    _rawHits->Fill();
      }

      //for (eiter = eloss.begin(); eiter != eloss.end(); ++eiter){
      //  (*eiter)->Fill(esum);
      //}
    //}
  }
  //cout << "CalculateDistortions::process_event(PHCompositeNode *topNode) Processing Event" << endl;
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
  //_rawHits->Draw("hit_r:hit_phi:hit_z >> hCharge_v1","(hit_eion)*(event_timestamp>0)");
  //int bX=7343773;//Ne CF4 90:10
  //_rawHits->Draw(Form("hit_r:hit_phi:hit_z-(%d-event_bunchXing)*106*%.10f>> hCharge_v1",bX,vIon/(m/ns)),Form("(hit_eion*%f)*(event_bunchXing<%d)",Tpc_ElectronsPerGeV,bX)); //A map for event right after bunch xing bX has occurred	
  //_rawHits->Draw(Form("hit_r:hit_phi:1.05-(%d-event_bunchXing)*106*%.10f>>+hCharge_v1",bX,vIon/(m/ns)),Form("(ibf_vol)*(event_bunchXing<%d && isOnPlane)",bX)); //A map of IBF for event right after bunch xing B has occurred

  outfile->cd();
  //ntup->Write();
  outfile->Write();
  outfile->Close();
  delete outfile;
  //hm->dumpHistos(_filename, "UPDATE");
  //cout << "CalculateDistortions::End(PHCompositeNode *topNode) This is the End..." << endl;
  //return Fun4AllReturnCodes::EVENT_OK;
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
  //if ( newBeamXing>=nBeams){
  //  _beamxing = newBeamXing - int(newBeamXing/nBeams)*nBeams;
  //  }
  //else{
  //  _beamxing =  newBeamXing;
  //  } 
  _beamxing = newBeamXing;
  cout<<"Initial BeamXing is set to: "<<newBeamXing<<endl;

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
