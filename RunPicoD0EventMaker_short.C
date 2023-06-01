#include <TSystem>

class StMaker;
class StChain;
class StPicoDstMaker;


StChain *chain;
//void runPicoD0EventMaker(const Char_t *inputFile="root://xrdstar.rcf.bnl.gov:1095//home/starlib/home/starreco/reco/AuAu_200_production_low_2014/ReversedFullField/P16id.SL22c/2014/095/15095020/st_physics_15095020_raw_1000015.picoDst.root", const Char_t *outputFile="test.root")
void RunPicoD0EventMaker_short(TString pico="TestLists/testPico_2016.list",
 TString outFileName="Test_EventMaker.root", int pYear = 2016)

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
      pico = "TestLists/testPico_2016.list";

   } else if (pYear==2014){
      SL_version = "SL22c";
      Runcode = "Run14_AuAu200_VpdMB5";
      prodID = "P16id";
      pico = "TestLists/testPico_2014.list";
   }else {
      cout << "\033[0;31m Not valid year.\033[0m"
      exit(0);
   }
   outFileName=Form("Test_EventMaker_%.d",pYear);

  string env_SL = getenv ("STAR");
  if(env_SL.find(SL_version)==string::npos)
  {
      cout << "\033[0;31mEnvironment Star Library does not match the requested library: \033[0m" << "\033[0;32m" << SL_version << "\033[0m";
      cout << "\033[0;31m for run: \033[0m" << "\033[0;32m" << pYear << "\033[0m";
      cout << "\033[0;31m in RunPicoTowerTest_short.C. Exiting...\033[0m" << endl;
      exit(1);
  }

  Int_t nEvents = 5000;
	
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();
  gROOT->LoadMacro("StRoot/macros/loadSharedAnalysisLibraries_EM.C");
  loadSharedAnalysisLibraries_EM();

	chain = new StChain();

	StPicoDstMaker* picoDstMaker = new StPicoDstMaker(2,pico,"picoDstMaker");
  StPicoD0EventMaker* picoD0Maker = new StPicoD0EventMaker("picoD0Maker",picoDstMaker,outFileName,pYear);

	chain->Init();
	cout<<"chain->Init();"<<endl;
	int total = picoDstMaker->chain()->GetEntries();
  cout << " Total entries = " << total << endl;
  

	nEvents = 200;
	for (Int_t i=0; i<= nEvents; i++)
  {
		
	  chain->Clear();
	  int iret = chain->Make(i);
		if(i%10==0) progres(i,nEvents);
	  if (iret) { cout << "Bad return code!" << iret << endl; break;}

	  total++;
	}
	 if(nEvents%10!=0) progres(i,nEvents);
	
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