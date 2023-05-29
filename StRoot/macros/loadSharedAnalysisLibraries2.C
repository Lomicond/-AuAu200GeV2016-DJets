#include "TSystem.h"

#ifndef __CINT__
#include "StMuDSTMaker/COMMON/macros/loadSharedLibraries.C"
#endif

extern TSystem* gSystem;

void loadSharedAnalysisLibraries2() 
{

  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();


  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StRefMultCorr");
  gSystem->Load("StPicoPrescales");

  gSystem->Load("St_db_Maker");
  gSystem->Load("StDaqLib");

  gSystem->Load("StEmcRawMaker");
  gSystem->Load("StEmcADCtoEMaker");
  gSystem->Load("StEpcMaker");
  gSystem->Load("StTriggerUtilities");

  

  gSystem->Load("StBTofUtil");
  gSystem->Load("StTpcDb");
  gSystem->Load("StDbUtilities");
  gSystem->Load("Sti");
  gSystem->Load("StiUtilities");
  gSystem->Load("StSsdDbMaker");
  gSystem->Load("StSvtDbMaker");
  gSystem->Load("StiMaker");
    gSystem->Load("StDbBroker");
  gSystem->Load("libgeometry_Tables"); //rember, order of loading makers matters
  gSystem->Load("StPicoTowerTest");


	//PYTHIA6 libraries
	//gSystem->Load("libEG.so");
	//gSystem->Load("libEGPythia6.so");
	//gSystem->Load("libPythia6.so");	

  cout << " loading of shared  libraries are done" << endl;
}
