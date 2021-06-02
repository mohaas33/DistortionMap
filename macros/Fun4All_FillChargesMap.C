#pragma once
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>

#include <CalculateDistortions.h>

#include <stdio.h>
//#include <sstream>

#include <string>


R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libCalculateDistortions.so)
R__LOAD_LIBRARY(libg4dst.so)

void Fun4All_FillChargesMap(  const int nEvents = 10, const int eventsInFileStart = 0, const string &fname = "/sphenix/sim/sim01/sphnxpro/Micromegas/2/G4Hits_sHijing_0-12fm_000000_001000.root", const string &foutputname = "/sphenix/user/shulga/Work/IBF/DistortionMap/Files/slim_G4Hits_sHijing_0-12fm_000000_001000.root" )
{
  ///////////////////////////////////////////
  // Make the Server
  //////////////////////////////////////////


  //char fname[100];
  //char foutputname[100];
  //for (int i=0;i<nFiles;i++){
  //int eventsInFileStart = i*1000;
  //int eventsInFileEnd = (i+1)*1000;
  //sprintf(fname, "/sphenix/sim/sim01/sphnxpro/Micromegas/2/G4Hits_sHijing_0-12fm_%06d_%06d.root",eventsInFileStart,eventsInFileEnd);
  //sprintf(foutputname, "/sphenix/user/shulga/Work/IBF/DistortionMap/Files/slim_G4Hits_sHijing_0-12fm_%06d_%06d.root",eventsInFileStart,eventsInFileEnd);

  Fun4AllServer *se = Fun4AllServer::instance();
  string cd_name = "CalculateDistortions"+std::to_string(eventsInFileStart);
  //cout<<fname_tmp<<endl;
  CalculateDistortions *dist_calc = new CalculateDistortions(cd_name, foutputname);
  dist_calc->SetFrequency(50);
  dist_calc->SetEvtStart(eventsInFileStart);
  dist_calc->SetBeamXing(1508071); // Set beam crosssing bias
  dist_calc->SetAvg(1); //Set average calculation
  dist_calc->SetUseIBFMap();
  //dist_calc->SetGain(2e3*48.7/71.5);
  dist_calc->SetGain(1400);
  dist_calc->SetIBF(0.004);
  dist_calc->UseSliming(0);//Turn off TTree filling and recording
  //Set pp colliding system
  dist_calc->SetCollSyst(0); //setting pp with = 1

  se->registerSubsystem(dist_calc);
  
  // this (DST) input manager just drives the event loop
  Fun4AllInputManager *in = new Fun4AllDstInputManager("DSTin");
  in->fileopen(fname);
  se->registerInputManager(in);
  // events = 0 => run till end of input file
  if (nEvents <= 0)
  {
    return;
  }
  cout << endl << "Running over " << nEvents << " Events" << endl;
  se->run(nEvents);
  //}
  cout << endl << "Calling End in Fun4All_CalculateDistortions.C" << endl;
  se->End();

  cout << endl << "All done, calling delete Fun4AllServer" << endl;
  delete se;

  cout << endl << "gSystem->Exit(0)" << endl;
  gSystem->Exit(0);
}