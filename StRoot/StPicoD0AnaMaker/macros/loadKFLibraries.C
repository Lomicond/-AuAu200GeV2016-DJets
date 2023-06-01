#include "TSystem.h"

extern TSystem* gSystem;


void loadSharedLibraries() {
  // Dynamically link needed shared lib
 	//Dynamically link needed shared libs
	gSystem->Load("St_base");
	gSystem->Load("StChain");
	gSystem->Load("St_Tables");
	gSystem->Load("StUtilities");
	gSystem->Load("StIOMaker");
	gSystem->Load("StarMagField");
	gSystem->Load("StarClassLibrary");
	//gSystem->Load("StTpcDb");
	gSystem->Load("StEvent");
	//gSystem->Load("StEventMaker"); //not needed if event.root branch present
  	gSystem->Load("StEventUtilities");
  	gSystem->Load("StDbLib");
	gSystem->Load("StEmcUtil"); 
	gSystem->Load("StEEmcUtil");
	gSystem->Load("StMcEvent");
	gSystem->Load("StMcEventMaker");
	gSystem->Load("StAssociationMaker");


	gSystem->Load("StTpcDb");
	gSystem->Load("StDbBroker");
	gSystem->Load("StDetectorDbMaker");
	gSystem->Load("StDbUtilities");
	gSystem->Load("St_db_Maker");
	gSystem->Load("StMagF");
        gSystem->Load("StPxlUtil");
        gSystem->Load("StPxlDbMaker");
        gSystem->Load("StIstUtil"); //YF
        gSystem->Load("StIstDbMaker"); //YF

        gSystem->Load("Sti");
	gSystem->Load("StiUtilities");
	gSystem->Load("StSsdDbMaker");
	gSystem->Load("StSvtDbMaker");
	gSystem->Load("StiMaker.so");
	gSystem->Load("libMathMore");
	gSystem->Load("libSpectrum");

 
  
 
  
  cout << " loading of KFVertex libraries done" << endl;
 }
