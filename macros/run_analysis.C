#include "analyzeHits.C"

#include <TH3.h>
#include <TFile.h>

#include <TStyle.h>
#include <TCanvas.h>

std::vector<int> readBeamXings();

std::vector<int> readBeamXings(){
  //cout << "CalculateDistortions::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << endl;
  std::vector<int> bXs;
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
        //_timestamps[n[0]]=n[1];
        bXs.push_back(int(n[0]));
      }
    }
    InputFile.close();
  }
 return bXs;
}


void run_analysis(int bXbias = 10){
    char fname[250];
    const int nFiles = 5;
    const int nHist = 10;
    int bX = 734587;
    //int bXbias = 10;//0;
    //analyzeHits aH[nFiles];
    std::vector<int> bXs = readBeamXings();
    //TH3D *hCharge[nFiles];
    //TFile *f[nFiles];
    TH3D *hCharge_tot[nHist];//=new TH3D("hCharge_tot","SC (ions) per m^3;phi (rad);r (m);z (m)",nphi,0,6.28319,nr,rmin,rmax,nz,0,z_rdo);
    for (int i=0;i<nFiles;i++){       
        int eventsInFileStart = i*1000;
        int eventsInFileEnd = (i+1)*1000;
        //cout<<"bX="<<bXs.at(i)<<endl;
        //sprintf(fname, "/sphenix/user/shulga/Work/IBF/DistortionMap/Files/slim_G4Hits_sHijing_0-12fm_%06d_%06d.root",eventsInFileStart,eventsInFileEnd);
        cout<<Form("/sphenix/user/shulga/Work/IBF/DistortionMap/Files/slim_G4Hits_sHijing_0-12fm_%06d_%06d.root",eventsInFileStart,eventsInFileEnd)<<endl;
        analyzeHits aH(Form("/sphenix/user/shulga/Work/IBF/DistortionMap/Files/slim_G4Hits_sHijing_0-12fm_%06d_%06d.root",eventsInFileStart,eventsInFileEnd));
        //sprintf(fname, "/sphenix/user/shulga/Work/IBF/DistortionMap/Files/outputFile_G4Hits_sHijing_0-12fm_%06d_%06d.root",eventsInFileStart,eventsInFileEnd);
        //sprintf(fname, "/sphenix/user/shulga/Work/IBF/DistortionMap/Files/outputFile_G4Hits_sHijing_0-12fm_%06d_%06d_bX%d_bias%d.root",eventsInFileStart,eventsInFileEnd,bX,bXbias);
        //cout<<Form("/sphenix/user/shulga/Work/IBF/DistortionMap/Files/outputFile_G4Hits_sHijing_0-12fm_%06d_%06d_bX%d_bias%d.root",eventsInFileStart,eventsInFileEnd,bX,bXbias)<<endl;
        aH.SetOutputFileName(Form("/sphenix/user/shulga/Work/IBF/DistortionMap/Files/outputFile_G4Hits_sHijing_0-12fm_%06d_%06d_bX%d_bias%d.root",eventsInFileStart,eventsInFileEnd,bX,bXbias));
        aH.SetBeamXings(bXs);
        aH.SetBeamXing(bX);
        aH.SetBeamXingBias(bXbias);    
        //aH[i] = analyzeHits(fname);
        //cout<<"Loop"<<endl;
        aH.Loop();
    }
    for (int i=0;i<nFiles;i++){       
        int eventsInFileStart = i*1000;
        int eventsInFileEnd = (i+1)*1000;
        //sprintf(fname, "/sphenix/user/shulga/Work/IBF/DistortionMap/Files/outputFile_G4Hits_sHijing_0-12fm_%06d_%06d_bX%d_bias%d.root",eventsInFileStart,eventsInFileEnd,bX,bXbias);
        //std::cout<<fname<<std::endl;
        //TFile *f = new TFile(fname);
        //TFile *f = new TFile("../Files/outputFile_G4Hits_sHijing_0-12fm_000000_001000.root");
        TFile *f = TFile::Open(Form("/sphenix/user/shulga/Work/IBF/DistortionMap/Files/outputFile_G4Hits_sHijing_0-12fm_%06d_%06d_bX%d_bias%d.root",eventsInFileStart,eventsInFileEnd,bX,bXbias));
        for (int hi=0;hi<nHist;hi++){       
            TH3D *hCharge_tmp = (TH3D*)f->Get(Form("h_Charge_%d",hi))->Clone(Form("h_Charge_%d_%d_tmp",i,hi));
            hCharge_tmp->SetDirectory(0);
            if (i==0){
                hCharge_tot[hi] = (TH3D*)hCharge_tmp->Clone(Form("h_Charge_%d_tot",hi));
            }else{
                hCharge_tot[hi]->Add(hCharge_tmp); 
            }
        }

    }

//TCanvas *c = new TCanvas();   
//hCharge_tot[0]->Draw();
//c->Print("./charges.png");
//   
//TCanvas *c_xy = new TCanvas();
//c_xy->GetPad(0)->SetLogz();
//c_xy->SetRightMargin(0.15);
//hCharge_tot[0]->Project3D("xy")->Draw("colz");
//c_xy->Print("./charges_xy.png");
//   
//TCanvas *c_z = new TCanvas();
//hCharge_tot[0]->SetLineColor(1);
//hCharge_tot[1]->SetLineColor(2);
//hCharge_tot[2]->SetLineColor(4);
//hCharge_tot[3]->SetLineColor(6);
//hCharge_tot[0]->SetLineWidth(2);
//hCharge_tot[1]->SetLineWidth(2);
//hCharge_tot[2]->SetLineWidth(2);
//hCharge_tot[3]->SetLineWidth(2);
//
//hCharge_tot[0]->Project3D("z")->Draw("HIST");
//hCharge_tot[1]->Project3D("z")->Draw("HISTsame");
//hCharge_tot[2]->Project3D("z")->Draw("HISTsame");
//hCharge_tot[3]->Project3D("z")->Draw("HISTsame");
//c_z->Print("./charges_z.png");

TFile outputFile (Form("/sphenix/user/shulga/Work/IBF/DistortionMap/Files/mapFile_bX%d_bias%d.root",bX,bXbias),"RECREATE");
for (int hi=0;hi<nHist;hi++){       
    hCharge_tot[hi]->Write();
}
outputFile.Close();



}