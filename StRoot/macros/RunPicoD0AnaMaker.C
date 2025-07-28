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
#include "StRefMultCorr/CentralityMaker.h"
#include "StRefMultCorr/StRefMultCorr.h"
#include "StPicoD0AnaMaker/StPicoD0AnaMaker.h"
#include "St_db_Maker/St_db_Maker.h"
#include "StEmcADCtoEMaker/StEmcADCtoEMaker.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoHFMaker/StHFCuts.h"
using namespace std;

#else
class StChain;
#endif

StChain *chain;


void loadSharedLibraries();
void loadSharedAnalysisLibraries();
void progres(double N, double D);

void RunPicoD0AnaMaker(string d0list="testD0.list", string pico="testPico.list",
 string outFileName="Test_AnaMaker.root"/*, string badRunListFileName = "picoList_bad_MB.list"*/,int pYear = 2014, bool Testing = true)
{
   //Check STAR Library. Please set SL_version to the original star library used in the production from http://www.star.bnl.gov/devcgi/dbProdOptionRetrv.pl
   string SL_version;
   string RefMult = "grefmult";
   string Runcode;
   string prodID;
   string D0list;
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

   int nEntries = 100000;

   if(Testing){
   pico=Form("TestLists/testPico_%.d.list",pYear);
   d0list=Form("TestLists/testD0_%.d.list",pYear);
   outFileName=Form("Test_AnaMaker_%.d.root",pYear);
   }

   string env_SL = getenv("STAR");
   if (env_SL.find(SL_version) == string::npos)
   {
      cout << "\033[0;31mEnvironment Star Library does not match the requested library: \033[0m" << "\033[0;32m" << SL_version << "\033[0m";
      cout << "\033[0;31m for run: \033[0m" << "\033[0;32m" << pYear << "\033[0m";
      cout << "\033[0;31m in RunPicoTowerTest_short.C. Exiting...\033[0m" << endl;
//      exit(1);
   }

   gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
   loadSharedLibraries();
   gSystem->Load("StBTofUtil");
   gSystem->Load("StPicoPrescales");
   gROOT->LoadMacro("StRoot/macros/loadSharedAnalysisLibraries.C");
   loadSharedAnalysisLibraries();
   chain = new StChain();

   StRefMultCorr* grefmultCorrUtil = new StRefMultCorr(RefMult, Runcode, prodID);
   StPicoDstMaker* picoDstMaker = new StPicoDstMaker(StPicoDstMaker::IoRead, pico.c_str(), "picoDstMaker");

   StMessMgr *msg = StMessMgr::Instance();
   msg->SwitchOff("Could not make BEMC detector");
   St_db_Maker *dbMaker = new St_db_Maker("db","MySQL:StarDb","$STAR/StarDb","StarDb");

   StEmcADCtoEMaker *adc = new StEmcADCtoEMaker();
   StPicoD0AnaMaker*  picoD0AnaMaker = new StPicoD0AnaMaker("picoD0AnaMaker", d0list.c_str(), outFileName.c_str(), picoDstMaker,grefmultCorrUtil,pYear, false, "");


  
   StHFCuts* d0Cuts = new StHFCuts("d0Cuts");

   // -------------- USER variables -------------------------
   picoD0AnaMaker->setHFCuts(d0Cuts);
   picoD0AnaMaker->setCutETmin(0.2);
   picoD0AnaMaker->setHadronCorr(1.);
   picoD0AnaMaker->setOnlyTrackBasedJets(false);
   picoD0AnaMaker->setMaxDcaZHadronCorr(3.0); //cm, max DCA_z for global tracks used for hadronic correction 
   picoD0AnaMaker->setNJetsRemove(1);
   picoD0AnaMaker->setGhostMaxrap(1.0);
   picoD0AnaMaker->setMaxNeutralFraction(0.95);
   // -- File name of bad run list
   //d0Cuts->setBadRunListFileName(badRunListFileName.c_str());

   // add your cuts here.

   // tracking
   d0Cuts->setCutNHitsFitnHitsMax(20);

   // pions
   d0Cuts->setCutTPCNSigmaPion(3.0);

   // kaons
   d0Cuts->setCutTPCNSigmaKaon(2.0);

   // kaonPion pair cuts
   float dcaDaughtersMax = 0.008;  // maximum
   float decayLengthMin  = 0.0030; // minimum
   float decayLengthMax  = 999999; //std::numeric_limits<float>::max();
   float cosThetaMin     = 0.90;   // minimum
   float minMass         = 1.6;
   float maxMass         = 2.1;
   d0Cuts->setCutSecondaryPair(dcaDaughtersMax, decayLengthMin, decayLengthMax, cosThetaMin, minMass, maxMass);

   chain->Init();
  
   if(!Testing){
      nEntries = picoD0AnaMaker->getEntries();
   }else
   {
   if (nEntries>picoD0AnaMaker->getEntries()) nEntries = picoD0AnaMaker->getEntries();
   }

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
