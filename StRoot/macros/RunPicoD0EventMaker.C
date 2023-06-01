#ifndef __CINT__
#include "TString.h"
#include <iostream>
#include <fstream>
#include <memory>
#include <ctime>
#include <cstdio>
#include <string>
#include <TROOT.h>
#include <TSystem.h>
#include "TChain.h"
#include "StChain/StMaker.h"
#include "StChain/StChain.h"
#include "StRefMultCorr/StRefMultCorr.h"
#include "StPicoD0EventMaker/StPicoD0EventMaker.h"
#include "St_db_Maker/St_db_Maker.h"
#include "StEmcADCtoEMaker/StEmcADCtoEMaker.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoHFMaker/StHFCuts.h"
using namespace std;

#else
class StChain;
#endif
class StMaker;
class StChain;
class StPicoDstMaker;


StChain *chain;

void loadSharedLibraries();
void loadSharedAnalysisLibraries_EM();
void progres(double N, double D);

void RunPicoD0EventMaker(string pico="TestLists/testPico_2016.list",
 string outFileName="Test_EventMaker.root", int pYear = 2014, bool Testing = true)

{ 
 
	 //Check STAR Library. Please set SL_version to the original star library used in the production from http://www.star.bnl.gov/devcgi/dbProdOptionRetrv.pl
   string SL_version = "Unknown";
   string RefMult = "grefmult";
   string Runcode;
   string prodID;
   if (pYear==2016){
      SL_version = "SL20c";
      Runcode = "Run16_AuAu200_VpdMB5";
      prodID = "P16ij";

   } else if (pYear==2014){
      SL_version = "SL22c";
      Runcode = "Run14_AuAu200_VpdMB5";
      prodID = "P16id";
   }else {
      cout << "\033[0;31m Not valid year.\033[0m"<<endl;
      exit(0);
   }

   Int_t nEvents = 1000;

   if(Testing){
   pico=Form("TestLists/testPico_%.d.list",pYear);
   outFileName=Form("Test_EventMaker_%.d",pYear);
   }

   string env_SL = getenv ("STAR");
   if(env_SL.find(SL_version)==string::npos)
   {
      cout << "\033[0;31mEnvironment Star Library does not match the requested library: \033[0m" << "\033[0;32m" << SL_version << "\033[0m";
      cout << "\033[0;31m for run: \033[0m" << "\033[0;32m" << pYear << "\033[0m";
      cout << "\033[0;31m in RunPicoTowerTest_short.C. Exiting...\033[0m" << endl;
      exit(1);
   }
	
   gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();
   gROOT->LoadMacro("StRoot/macros/loadSharedAnalysisLibraries_EM.C");
   loadSharedAnalysisLibraries_EM();

	chain = new StChain();

	StPicoDstMaker* picoDstMaker = new StPicoDstMaker(StPicoDstMaker::IoRead,pico.c_str(),"picoDstMaker");
   StPicoD0EventMaker* picoD0Maker = new StPicoD0EventMaker("picoD0Maker",picoDstMaker,outFileName.c_str(),pYear);

	chain->Init();
	cout<<"chain->Init();"<<endl;

	if(!Testing){
      nEvents = picoDstMaker->chain()->GetEntries();
   }

   cout << " Total entries = " << nEvents << endl;
   //-------------------------------------------
   for (int iEvent = 0; iEvent <= nEvents; ++iEvent)
   {
      chain->Clear();
      int iret = chain->Make();
      if(iEvent%10==0) progres(iEvent,nEvents);
      if (iret)
      {
         cout << "Bad return code!" << iret << endl;
         break;
      }
   }
   if(nEvents%10!=0) progres(nEvents,nEvents);
   //-------------------------------------------
	
	cout << "****************************************** " << endl;
	cout << "Work done... now its time to close up shop!"<< endl;
	cout << "****************************************** " << endl;
	chain->Finish();
	cout << "****************************************** " << endl;
	cout << "total number of events  " << nEvents << endl;
	cout << "****************************************** " << endl;
	
	delete chain;
}
   
   void progres(double citatel, double jmenovatel){
      int Ndilky=50;
      int proc=floor(citatel/jmenovatel*Ndilky);
      
      cout << "\r"<< flush;
      cout << "  │" << flush;
      for (int i=1; i<=proc; i++){
      cout <<"█"<<flush; 
      }
      for (int j=proc+1; j<=Ndilky;j++){
      cout <<  "░"<<flush;
      }
      cout << "│ " << flush;
      if(citatel!=jmenovatel){
         cout << Form("Completed: %.2f ",citatel/jmenovatel*100.)<<"%"<<flush;
      }
               
      else{
         cout << Form("\033[1;32m Completed: %.2f \033[0m",citatel/jmenovatel*100.)<<"\033[1;32m%\033[0m"<<endl;
      }
      return;
    }