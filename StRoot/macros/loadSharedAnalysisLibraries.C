#include "TSystem.h"

#ifndef __CINT__
#include "StMuDSTMaker/COMMON/macros/loadSharedLibraries.C"
#endif

extern TSystem* gSystem;

void loadSharedAnalysisLibraries() 
{

  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  gSystem->Load("/direct/star+u/robotmon/fastjet1/fastjet_install/lib/libfastjet");
  gSystem->Load("/direct/star+u/robotmon/fastjet1/fastjet_install/lib/libsiscone");
  gSystem->Load("/direct/star+u/robotmon/fastjet1/fastjet_install/lib/libsiscone_spherical"); 
  gSystem->Load("/direct/star+u/robotmon/fastjet1/fastjet_install/lib/libfastjetplugins");
  gSystem->Load("/direct/star+u/robotmon/fastjet1/fastjet_install/lib/libfastjettools");
  gSystem->Load("/direct/star+u/robotmon/fastjet1/fastjet_install/lib/libfastjetcontribfragile"); 

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

  
 gSystem->Load("StPicoCutsBase");
  gSystem->Load("StBTofUtil");
  gSystem->Load("StPicoD0EventMaker");
  gSystem->Load("StPicoHFMaker");
  gSystem->Load("StPicoCuts");
  gSystem->Load("StTpcDb");
  gSystem->Load("StDbUtilities");
  gSystem->Load("Sti");
  gSystem->Load("StiUtilities");
  gSystem->Load("StSsdDbMaker");
  gSystem->Load("StSvtDbMaker");
  gSystem->Load("StiMaker");
  gSystem->Load("StDbBroker");
  gSystem->Load("libgeometry_Tables"); //rember, order of loading makers matters
  gSystem->Load("StPicoD0AnaMaker");


	//PYTHIA6 libraries
	gSystem->Load("libEG.so");
	gSystem->Load("libEGPythia6.so");
	gSystem->Load("libPythia6.so");	

  cout << " loading of shared  libraries are done" << endl;
}
