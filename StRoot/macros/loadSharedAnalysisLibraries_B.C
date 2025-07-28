#include "TSystem.h"

#ifndef __CINT__
#include "StMuDSTMaker/COMMON/macros/loadSharedLibraries.C"
#endif

extern TSystem* gSystem;

void loadSharedAnalysisLibraries() 
{

  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  /*
  gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/fastjet-install/lib/libfastjet");
  gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/fastjet-install/lib/libsiscone");
  gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/fastjet-install/lib/libsiscone_spherical"); 
  gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/fastjet-install/lib/libfastjetplugins");
  gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/fastjet-install/lib/libfastjettools");
  */
    gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/Alex_install/install/fastjet-install/lib/libfastjet");
    gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/Alex_install/install/fastjet-install/lib/libsiscone");
    gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/Alex_install/install/fastjet-install/lib/libsiscone_spherical");
    gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/Alex_install/install/fastjet-install/lib/libfastjetplugins");
    gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/Alex_install/install/fastjet-install/lib/libfastjettools");
    gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/Alex_install/install/fastjet-install/lib/libfastjetcontribfragile");
  //gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/fastjet-install/lib/libfastjetcontribfragile"); 

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
   // gSystem->Load("StPicoBackground");
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
  //gSystem->Load("StPicoD0AnaMaker");
   gSystem->Load("StPicoBackground");

	//PYTHIA6 libraries
	gSystem->Load("libEG.so");
	gSystem->Load("libEGPythia6.so");
	gSystem->Load("libPythia6.so");	

  cout << " loading of shared  libraries are done" << endl;
}
