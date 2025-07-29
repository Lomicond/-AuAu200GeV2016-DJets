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
#include "StPicoEvent/StPicoMcVertex.h"
#include "StRefMultCorr/CentralityMaker.h"
#include "StRefMultCorr/StRefMultCorr.h"
#include "StPicoD0JetAnaMaker/StPicoD0JetAnaMaker.h"
#include "St_db_Maker/St_db_Maker.h"
#include "StEmcADCtoEMaker/StEmcADCtoEMaker.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
using namespace std;

#else
class StChain;
#endif

StChain *chain;


void loadSharedLibraries();
void loadSharedAnalysisLibraries();
void progres(double N, double D);

void RunPicoD0JetAnaMaker(string pico="testPico.list",
 string outFileName="Test_AnaMaker.root"/*, string badRunListFileName = "picoList_bad_MB.list"*/,int pYear = 2014, bool Testing = true)
{
   //Check STAR Library. Please set SL_version to the original star library used in the production from http://www.star.bnl.gov/devcgi/dbProdOptionRetrv.pl
   string SL_version;
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

   int nEntries = 1000000;

   if(Testing){
   pico=Form("TestLists/testPico_%.d.list",pYear);
   outFileName=Form("Test_AnaMaker_%.d.root",pYear);
   }

   string env_SL = getenv("STAR");
   if (env_SL.find(SL_version) == string::npos)
   {
      cout << "\033[0;31mEnvironment Star Library does not match the requested library: \033[0m" << "\033[0;32m" << SL_version << "\033[0m";
      cout << "\033[0;31m for run: \033[0m" << "\033[0;32m" << pYear << "\033[0m";
      cout << "\033[0;31m in RunPicoTowerTest_short.C. Exiting...\033[0m" << endl;
      exit(1);
   }

   gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
   loadSharedLibraries();
   gSystem->Load("StBTofUtil");
   gSystem->Load("StPicoPrescales");
   gROOT->LoadMacro("StRoot/macros/loadSharedAnalysisLibrariesJet.C");
   loadSharedAnalysisLibraries();
   chain = new StChain();

   StRefMultCorr* grefmultCorrUtil = new StRefMultCorr(RefMult, Runcode, prodID);
   StPicoDstMaker* picoDstMaker = new StPicoDstMaker(StPicoDstMaker::IoRead, pico.c_str(), "picoDstMaker");

   StMessMgr *msg = StMessMgr::Instance();
   msg->SwitchOff("Could not make BEMC detector");
   St_db_Maker *dbMaker = new St_db_Maker("db","MySQL:StarDb","$STAR/StarDb","StarDb");

   StEmcADCtoEMaker *adc = new StEmcADCtoEMaker();
   StPicoD0JetAnaMaker*  picoD0JetAnaMaker = new StPicoD0JetAnaMaker("picoD0JetAnaMaker", outFileName.c_str(), picoDstMaker,grefmultCorrUtil,pYear);

   // -------------- USER variables -------------------------
   picoD0JetAnaMaker->setCutETmin(0.2);
   picoD0JetAnaMaker->setHadronCorr(1.);
   picoD0JetAnaMaker->setOnlyTrackBasedJets(false);
   picoD0JetAnaMaker->setMaxDcaZHadronCorr(3.0); //cm, max DCA_z for global tracks used for hadronic correction 
   picoD0JetAnaMaker->setNJetsRemove(1);
   picoD0JetAnaMaker->setGhostMaxrap(1.0);
   picoD0JetAnaMaker->setMaxNeutralFraction(0.95);
   picoD0JetAnaMaker->setBgSubtraction(2); //0 - Area based, 1 - ICS, 2 - jet shape

   chain->Init();
  
   if(!Testing) nEntries = picoDstMaker->chain()->GetEntries();  
   else if (nEntries>picoDstMaker->chain()->GetEntries()) nEntries = picoDstMaker->chain()->GetEntries();

   cout << "nEntries: " << nEntries << endl;
   //--------------------------------------------------------   
   for (int iEvent = 0; iEvent < nEntries; ++iEvent)
   {
      chain->Clear();
      int iret = chain->Make();
      if(Testing && iEvent%200==0) progres(iEvent,nEntries);
      if (iret)
      {
         cout << "Bad return code!" << iret << endl;
         break;
      }
   }
   if(Testing &&nEntries%200!=0) progres(nEntries,nEntries);
   //--------------------------------------------------------

   chain->Finish();
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
         cout << Form("Completed: %.2f",citatel/jmenovatel*100.)<<"% "<<flush;
      }
               
      else{
         cout << Form("\033[1;32m Completed: %.2f \033[0m",citatel/jmenovatel*100.)<<"\033[1;32m%\033[0m"<<endl;
      }
      return;
    }
