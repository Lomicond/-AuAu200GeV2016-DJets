/* **************************************************
 *  A macro to run StPicoD0AnaMaker
 *
 *  Authors:  **Mustafa Mustafa (mmustafa@lbl.gov)
 *
 *  **Code Maintainer
 *
 * **************************************************
 */

void RunPicoTowerTest_short(TString pico="TestLists/testPico_2016.list",
 TString outFileName="Test_Tower.root", TString badRunListFileName = "picoList_bad_MB.list",  int pYear = 2016)
{
   //Check STAR Library. Please set SL_version to the original star library used in the production from http://www.star.bnl.gov/devcgi/dbProdOptionRetrv.pl
   string SL_version = "Unknown";
   if (pYear==2016){
      SL_version = "SL20c";
   } else if (pYear==2014){
      SL_version = "SL22c";
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
   gROOT->LoadMacro("StRoot/macros/loadSharedAnalysisLibraries2.C");
    loadSharedAnalysisLibraries2();
   chain = new StChain();

    StRefMultCorr* grefmultCorrUtil = new StRefMultCorr("grefmult", "Run16_AuAu200_VpdMB5", "P16ij");

   StPicoDstMaker* picoDstMaker = new StPicoDstMaker(2, pico, "picoDstMaker");

   StMessMgr *msg = StMessMgr::Instance();
   msg->SwitchOff("Could not make BEMC detector");
   St_db_Maker *dbMaker = new St_db_Maker("db","MySQL:StarDb","$STAR/StarDb","StarDb");
   StEmcADCtoEMaker *adc = new StEmcADCtoEMaker();
   StPicoTowerTest*  picoTowerTest = new StPicoTowerTest("picoTowerTest", outFileName.Data(), picoDstMaker,grefmultCorrUtil);


   //StHFCuts* d0Cuts = new StHFCuts("d0Cuts");
   //picoTowerTest->setHFCuts(d0Cuts);

   picoTowerTest->setCutETmin(0.2);
   picoTowerTest->setHadronCorr(1.);
   // -------------- USER variables -------------------------
   picoTowerTest->setMaxDcaZHadronCorr(3.0); //cm, max DCA_z for global tracks used for hadronic correction 
   // -- File name of bad run list
   //d0Cuts->setBadRunListFileName(badRunListFileName);

   // add your cuts here.

   // tracking
   //d0Cuts->setCutNHitsFitMax(20);
   //d0Cuts->setCutNHitsFitnHitsMax(20);

   // pions
   //d0Cuts->setCutTPCNSigmaPion(3.0);

   // kaons
   //d0Cuts->setCutTPCNSigmaKaon(2.0);

   // kaonPion pair cuts
   float dcaDaughtersMax = 0.008;  // maximum
   float decayLengthMin  = 0.0030; // minimum
   float decayLengthMax  = 999999; //std::numeric_limits<float>::max();
   float cosThetaMin     = 0.90;   // minimum
   float minMass         = 1.6;
   float maxMass         = 2.1;
   //d0Cuts->setCutSecondaryPair(dcaDaughtersMax, decayLengthMin, decayLengthMax, cosThetaMin, minMass, maxMass);

   chain->Init();
  
   //int nEntries = picoTowerTest->getEntries();
   int nEntries = 100000000000;
   for (int iEvent = 0; iEvent < nEntries; ++iEvent)
   {
      //if(iEvent%1000==0)
      cout << "Working on eventNumber " << iEvent << endl;
      chain->Clear();
      int iret = chain->Make();
      if (iret)
      {
         cout << "Bad return code!" << iret << endl;
         break;
      }
   }

   chain->Finish();
   delete chain;

   // delete list of picos
  // command = "rm -f correspondingPico.list";
   //gSystem->Exec(command.Data());

}
