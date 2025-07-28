void runPicoD0AnaMaker(TString d0list="file.list",TString outFileName="test.root")
{
  //Check STAR Library. Please set SL_version to the original star library used in the production from http://www.star.bnl.gov/devcgi/dbProdOptionRetrv.pl
  string SL_version = "SL15c";
  string env_SL = getenv ("STAR");
  if(env_SL.find(SL_version)==string::npos)
  {
      cout<<"Environment Star Library does not match the requested library in runPicoD0EventMaker.C. Exiting..."<<endl;
      exit(1);
  }

  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StPicoPrescales");
  gSystem->Load("StPicoD0EventMaker");
  gROOT->LoadMacro("./loadKFLibraries.C");
  loadSharedLibraries();
	// gSystem->Load("libStPicoKFVertexFitter");
  gSystem->Load("StPicoD0AnaMaker");
  gSystem->Load("StRefMultCorr");
  // gSystem->Load("StPicoHFMaker");
	gSystem->Load("libMinuit");

  chain = new StChain();

  // create list of picoDst files
  //TString command = "sed 's/hft\\\/d0tree/picodsts/g' "+d0list+" >correspondingPico.list";
  TString command = "sed 's/hft\\\/d0tree/picodsts/g' "+d0list+" >correspondingPico.list";
  gSystem->Exec(command.Data());
  command = "sed -i 's/picoD0/picoDst/g;s/kfProd2/physics2/g' correspondingPico.list";
  //command = "sed -i 's/kfTest/physics2/g' correspondingPico.list";
  gSystem->Exec(command.Data());
  StPicoDstMaker* picoDstMaker = new StPicoDstMaker(0,"correspondingPico.list","picoDstMaker");
  StRefMultCorr* grefmultCorrUtil  = CentralityMaker::instance()->getgRefMultCorr();
  StPicoD0AnaMaker*  picoD0AnaMaker = new StPicoD0AnaMaker("picoD0AnaMaker",d0list,outFileName.Data(),picoDstMaker,grefmultCorrUtil);
//  StPicoD0AnaMaker*  picoD0AnaMaker = new StPicoD0AnaMaker("picoD0AnaMaker",d0list,outFileName.Data(), picoDstMaker, grefmultCorrUtil);
  grefmultCorrUtil->setVzForWeight(6, -6.0, 6.0);
  grefmultCorrUtil->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14.txt");

  StHFCuts* d0Cuts = new StHFCuts("d0Cuts");
  picoD0AnaMaker->setHFCuts(d0Cuts);

	// -------------- USER variables -------------------------
	// add your cuts here. 

	// tracking
  // d0Cuts->setCutNHitsFitMax(20);

	// pions
  // d0Cuts->setCutTPCNSigmaPion(3.0);

	// kaons
  // d0Cuts->setCutTPCNSigmaKaon(2.0);
   
	// kaonPion pair cuts
  float dcaDaughtersMax = 0.008;  // maximum
  float decayLengthMin  = 0.0030; // minimum
  float decayLengthMax  = 999999; //std::numeric_limits<float>::max();
  float cosThetaMin     = 0.90;   // minimum
  float minMass         = 1.6;
  float maxMass         = 2.1;
  // d0Cuts->setCutSecondaryPair(dcaDaughtersMax, decayLengthMin, decayLengthMax, cosThetaMin, minMass, maxMass);
	  
  chain->Init();
  int nEntries = picoD0AnaMaker->getEntries();
  cout << " Total entries = " << nEntries << endl;
  time_t t_1,t_2;
  t_1 = clock();
  for(int iEvent = 0; iEvent<nEntries; ++iEvent)
  {
    if(iEvent%100==0)
                cout << "Working on eventNumber " << iEvent << endl;
    chain->Clear();
    int iret = chain->Make();
	  if (iret) { cout << "Bad return code!" << iret << endl; break;}
  }
  t_2 = clock();
  float dtime = 0.000001*difftime(t_2,t_1);
  cout<<"Job is done, using time is "<<dtime<<"s for "<<nEntries<<"events"<<endl;
  chain->Finish();
  delete chain;

  // delete list of picos
//  command = "rm -f correspondingPico.list";
  gSystem->Exec(command.Data());

}
