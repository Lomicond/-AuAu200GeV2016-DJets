void RunPicoD0AnaMaker_short(TString d0list="testD0.list", TString pico="testPico.list",
 TString outFileName="Test_AnaMaker.root", TString badRunListFileName = "picoList_bad_MB.list",int pYear = 2016, bool Testing = true)
{
   //Check STAR Library. Please set SL_version to the original star library used in the production from http://www.star.bnl.gov/devcgi/dbProdOptionRetrv.pl
   string SL_version = "Unknown";
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
   } 

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
      exit(1);
   }

   int nEntries = 200;

   gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
   loadSharedLibraries();

   gSystem->Load("StBTofUtil");
   gSystem->Load("StPicoPrescales");
   gROOT->LoadMacro("StRoot/macros/loadSharedAnalysisLibraries.C");
   loadSharedAnalysisLibraries();
   chain = new StChain();

   StRefMultCorr* grefmultCorrUtil = new StRefMultCorr(RefMult, Runcode, prodID);
   StPicoDstMaker* picoDstMaker = new StPicoDstMaker(2, pico, "picoDstMaker");

   StMessMgr *msg = StMessMgr::Instance();
   msg->SwitchOff("Could not make BEMC detector");
   St_db_Maker *dbMaker = new St_db_Maker("db","MySQL:StarDb","$STAR/StarDb","StarDb");

   StEmcADCtoEMaker *adc = new StEmcADCtoEMaker();
   StPicoD0AnaMaker*  picoD0AnaMaker = new StPicoD0AnaMaker("picoD0AnaMaker", d0list, outFileName.Data(), picoDstMaker,grefmultCorrUtil);

  
   StHFCuts* d0Cuts = new StHFCuts("d0Cuts");
   //Zkontroluji, jestli se mi to nacte
   picoD0AnaMaker->setHFCuts(d0Cuts);
   picoD0AnaMaker->setCutETmin(0.2);
   picoD0AnaMaker->setHadronCorr(1.);
   // -------------- USER variables -------------------------
   picoD0AnaMaker->setMaxDcaZHadronCorr(3.0); //cm, max DCA_z for global tracks used for hadronic correction 
   // -- File name of bad run list
   d0Cuts->setBadRunListFileName(badRunListFileName);

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
  
   //int nEntries = picoD0AnaMaker->getEntries();
    cout << "nEntries: " << nEntries << endl;
   for (int iEvent = 0; iEvent <= nEntries; ++iEvent)
   {
      chain->Clear();
      int iret = chain->Make();
      if(iEvent%10==0) progres(iEvent,nEntries);
      if (iret)
      {
         cout << "Bad return code!" << iret << endl;
         break;
      }
   }
   if(nEntries%10!=0) progres(iEvent,nEntries);


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
         cout << Form("Completed: %.2f ",citatel/jmenovatel*100.)<<"%"<<flush;
      }
               
      else{
         cout << Form("\033[1;32m Completed: %.2f \033[0m",citatel/jmenovatel*100.)<<"\033[1;32m%\033[0m"<<endl;
      }
      return;
    }