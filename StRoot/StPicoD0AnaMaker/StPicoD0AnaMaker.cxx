#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TLorentzVector.h"

#include "phys_constants.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoD0EventMaker/StPicoD0Event.h"
#include "StPicoD0EventMaker/StKaonPion.h"
#include "StPicoD0AnaMaker.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StCuts.h"
#include "JetInfo.h"
#include "../StPicoPrescales/StPicoPrescales.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StEmcADCtoEMaker/StBemcData.h"
#include "StEmcADCtoEMaker/StEmcADCtoEMaker.h"
#include "StEmcRawMaker/StBemcRaw.h"
#include "StEmcRawMaker/defines.h"
#include "StEmcRawMaker/StBemcTables.h"
#include "BemcNewCalib.h"
#include "Calibration2016.h"
//////Refit include lib
#include "PhysicalConstants.h"
#include "StThreeVectorF.hh"
#include "StLorentzVectorD.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TFile.h"
#include "StEvent/StDcaGeometry.h"
//
#include <vector>
#include <stdio.h>
#include <time.h>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cmath>
//---------------------------------
#include <fastjet/config.h>
#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
//#ifdef FASTJET_VERSION
#include <fastjet/Selector.hh>
#include <fastjet/tools/Subtractor.hh>
#include <fastjet/tools/Recluster.hh>
#include <fastjet/tools/JetMedianBackgroundEstimator.hh>
//#include <fastjet/contrib/ConstituentSubtractor.hh>
//#include <fastjet/contrib/SoftDrop.hh>  //robotmon 
//#include <fastjet/contrib/RecursiveSoftDrop.hh>
//#include <fastjet/contrib/SoftKiller.hh>
using namespace std;  //robotmon
using namespace fastjet;
//-------------------------------------

ClassImp(StPicoD0AnaMaker)

  StPicoD0AnaMaker::StPicoD0AnaMaker(char const * name,char const * inputFilesList, char const * outName,StPicoDstMaker* picoDstMaker,StRefMultCorr* grefmultCorrUtil, int pYear):
    StMaker(name),mPicoDstMaker(picoDstMaker),mPicoD0Event(NULL), mGRefMultCorrUtil(grefmultCorrUtil),
    mOutFileName(outName), mInputFileList(inputFilesList),mOutputFile(NULL), mChain(NULL), mEventCounter(0),mYear(pYear){}

Int_t StPicoD0AnaMaker::Init()
{
  mPicoD0Event = new StPicoD0Event();

  mChain = new TChain("T");
  std::ifstream listOfFiles(mInputFileList.Data());
  if (listOfFiles.is_open())
  {
    std::string file;
    while (getline(listOfFiles, file))
    {
      LOG_INFO << "StPicoD0AnaMaker - Adding :" << file <<endm;
      mChain->Add(file.c_str());
      LOG_INFO<<" Entries = "<<mChain->GetEntries()<< endm; 
    }
  }
  else
  {
    LOG_ERROR << "StPicoD0AnaMaker - Could not open list of files. ABORT!" << endm;
    return kStErr;
  }

  mPrescales = new StPicoPrescales(mycuts::prescalesFilesDirectoryName);

  mChain->GetBranch("dEvent")->SetAutoDelete(kFALSE);
  mChain->SetBranchAddress("dEvent", &mPicoD0Event);

  mOutputFile = new TFile(mOutFileName.Data(), "RECREATE");
  mOutputFile->cd();

  //Histograms
  vtxz = new TH1F("vtxz",";PVtx.z() [cm]; Count",100,-10,10);
  D0etalike = new TH1D("D0etalike",";#eta; Count",100,-5,5);
  D0etaunlike = new TH1D("D0etaunlike",";#eta; Count",100,-5,5);
  pioneta = new TH1D("pioneta",";#eta; Count",100,-5,5);
  kaoneta =  new TH1D("kaoneta",";#eta; Count",100,-5,5);
  vtxr = new TH2D("vtxr",";PVtx.x() [cm]; PVtx.y() [cm]",100,-3,3,100,-3,3);
  hcentr = new TH1D("hcentr",";C_{ID}",9,-0.5,8.5);
  NEvent = new TH1D("NEvent","",1,-1,1);
  //angledistr = new TH1D("AngleDistr","",100,-7,-7);



  int nRuns = mPrescales->numberOfRuns();
  mh1TotalEventsInRun = new TH1F("mh1TotalEventsInRun","totalEventsInRun;runIndex;totalEventsInRun",nRuns+1,0,nRuns+1);
  mh1TotalGRefMultInRun = new TH1F("mh1TotalGRefMultInRun","totalGRefMultInRun;runIndex;totalGRefMultInRun",nRuns+1,0,nRuns+1);
  mh1TotalKaonsInRun = new TH1F("mh1TotalKaonsInRun","totalKaonsInRun;runIndex;totalKaonsInRun",nRuns+1,0,nRuns+1);
  mh1TotalPionsInRun = new TH1F("mh1TotalPionsInRun","totalPionsInRun;runIndex;totalPionsInRun",nRuns+1,0,nRuns+1);
  mh1TotalD0CandidatesInRun = new TH1F("mh1TotalD0CandidatesInRun","totalD0CandidatesInRun;runIndex;totalD0CandidatesInRun",nRuns+1,0,nRuns+1);
  mh2NKaonsVsNPions = new TH2F("mh2NKaonsVsNPions","nKaonsVsNPions;nPions;nKaons",1000,0,1000,300,0,300);
  mh2KaonDcaVsPt = new TH2F("mh2KaonDcaVsPt","kaonDcaVsPt;p_{T}(K#pi)(GeV/c);K DCA(cm)",120,0,12,50,0,0.05);
  mh2PionDcaVsPt = new TH2F("mh2PionDcaVsPt","pionDcaVsPt;p_{T}(K#pi)(GeV/c);#pi DCA(cm)",120,0,12,50,0,0.05);
  mh2CosThetaVsPt = new TH2F("mh2CosThetaVsPt","cosThetaVsPt;p_{T}(K#pi)(GeV/c);cos(#theta)",120,0,12,500,0,1.0);
  mh2DcaDaughtersVsPt = new TH2F("mh2DcaDaughtersVsPt","dcaDaughtersVsPt;p_{T}(K#pi)(GeV/c);dcaDaughters(cm)",120,0,12,200,0,0.02);
  mh2InvariantMassVsPtUnlike = new TH2F("mh2InvariantMassVsPtUnlike","invariantMassVsPtUnlike;p_{T}(K#pi)(GeV/c);m_{K#pi}(GeV/c^{2})",120,0,12,220,0,2.2);
  mh2InvariantMassVsPtLike = new TH2F("mh2InvariantMassVsPtLike","invariantMassVsPtLike;p_{T}(K#pi)(GeV/c);m_{K#pi}(GeV/c^{2})",120,0,12,220,0,2.2);

  hpT_tr = new TH1F("hpT_tr","hpT_tr;p_{T}(GeV/c);",200,0,100);
  heta_phi_tr = new TH2F("heta_phi_tr","heta_phi_tr;#phi;#eta",120,0,12,50,0,0.05);
  heta_tr = new TH1F("heta_tr","heta_tr;#eta",50,-6,6);
  hphi_tr = new TH1F("hphi_tr","hphi_tr;#phi;",40,-6.2830,6.2830);
  hdca_z_tr = new TH2F("hdca_z_tr","hdca_z_tr;DCA(cm); v_{z}(GeV/c)",120,0,12,50,-10,10);
  hdca_pT = new TH2F("hdca_pT","hdca_pT; DCA(cm); p_{T}(GeV/c) ",120,0,12,100,0,200);
  hdca_tr = new TH1F("hdca_tr","hdca_tr; DCA(cm)",120,0,12);
  hcharged_tr = new TH1F("hcharged_tr","hcharged_tr;Charge;",5,-1,1);

//jet #        rapidity             phi              pt           index
  Jets = new TNtuple("Jets", "Jets", "RunId:EventID:NJet:rapidity:phi:pt:index:D0mass:D0_r:D0_pT:lambda_1_1:z:E");
  //D0_Daughter       = 0 not D0, antiD0 nor daughter
  //                  = 1 D0
  //                  = -1 anti D0
  //                  = 2 daughter pion
  //                  = -2 daughter kaon

  const int xbinSize=100;

  float xbin[101];
    for(int i=0;i<101;i++)
    xbin[i] = 0.1*i;


  float binMass[2001];

  candPt = new TProfile("candPt",";D0Pt [GeV/c];D0Pt [GeV/c]",xbinSize,xbin);

  for(int i=0;i<2001;i++)
    binMass[i] = 0.01*i;
  massPt = new TH2D("massPt",";M_{K#pi} [GeV/c^{2}]; p_{T} [GeV/c]",2000,binMass,xbinSize,xbin);
  massPtLike = new TH2D("massPtLike",";M_{K#pi} [GeV/c^{2}]; p_{T} [GeV/c]",2000,binMass,xbinSize,xbin);
  //angledistrLike = new TH1D("AngleDistrLike","",200,3.14125,3.14165);
  //angledistrUnlike = new TH1D("AngleDistrUnlike","",200,3.14125,3.14165);
  //massAngleLike = new TH2D("massAngleLike",";M_{K#pi} [GeV/c^{2}]; angle",2000,0,2000*0.01,200,3.14125,3.14165);
  //massAngleUnlike = new TH2D("massAngleUnlike",";M_{K#pi} [GeV/c^{2}]; angle",2000,0,2000*0.01,200,3.14125,3.14165);

 
  float ptbin1[12] = {0.225,0.375,0.525,0.675,0.825,0.975,1.12,1.27,1.42,1.58,1.73,1.88};
  
  mOutputFile->cd();
StMaker* maker = GetMaker("Eread");
if (maker) {
  cout << "Maker is valid!" << endl;
} else {
  cout << "Maker is invalid!" << endl;
}
  mADCtoEMaker = dynamic_cast<StEmcADCtoEMaker*>(GetMaker("Eread"));
  if (mADCtoEMaker) {
  cout << "mADCtoEMaker is valid!" << endl;
} else {
  cout << "mADCtoEMaker is invalid!" << endl;
}
  assert(mADCtoEMaker);
  mTables = mADCtoEMaker->getBemcData()->getTables();

  // -------------- USER VARIABLES -------------------------
  mGRefMultCorrUtil = new StRefMultCorr("grefmult");

  return kStOK;
}
//-----------------------------------------------------------------------------
StPicoD0AnaMaker::~StPicoD0AnaMaker()
{
  /*  */
  delete mGRefMultCorrUtil;
}

//-------------------------------------------
struct FourMomentum {
    double E, px, py, pz, D0_antiD0;
};

double delta_R(double eta1, double phi1, double eta2, double phi2) {
    double deta = eta1 - eta2;
    double dphi = TVector2::Phi_mpi_pi(phi1 - phi2);
    return sqrt(deta*deta + dphi*dphi);
}
//-----------------------------------------------------------------------------
Int_t StPicoD0AnaMaker::Finish()
{
  LOG_INFO << " StPicoD0AnaMaker - writing data and closing output file " <<endm;
  fout.close();
  fout1.close();
  mOutputFile->cd();
  // save user variables here
  massPt->Write();
  massPtLike->Write();
    //massPtacut->Write();
  //massPtLikeacut->Write();
  D0etalike->Write();
  D0etaunlike->Write();
  pioneta->Write();
  kaoneta->Write();
  candPt->Write();
  vtxz->Write();
  vtxr->Write();
  NEvent->Write();
  hcentr->Write();
  //angledistr->Write();

    //angledistrLike->Write();
    //angledistrUnlike->Write();
    //massAngleLike->Write();
    //massAngleUnlike->Write();

    mh1TotalEventsInRun->Write();
    mh1TotalGRefMultInRun->Write();
    mh1TotalKaonsInRun->Write();
    mh1TotalPionsInRun->Write();
    mh1TotalD0CandidatesInRun->Write();
    mh2NKaonsVsNPions->Write();

    mh2KaonDcaVsPt->Write();
    mh2PionDcaVsPt->Write();
    mh2CosThetaVsPt->Write();
    mh2DcaDaughtersVsPt->Write();

         hpT_tr->Write();
      heta_phi_tr->Write();
      heta_tr->Write();
      hphi_tr->Write();
      hdca_z_tr->Write();
      hdca_pT->Write();
      hdca_tr->Write();
      hcharged_tr->Write();


 Jets->Write();   

  mOutputFile->Close();
  delete mPrescales;

  return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoD0AnaMaker::Make()
{
  readNextEvent();
  if (!mPicoDstMaker)
  {
    LOG_WARN << " StPicoD0AnaMaker - No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  picoDst = mPicoDstMaker->picoDst();

  if (!picoDst)
  {
    LOG_WARN << "StPicoD0AnaMaker - No PicoDst! Skip! " << endm;
    return kStWarn;
  }
  //cout << mPicoD0Event->runId() << " " << picoDst->event()->runId() << endl;
  //cout << mPicoD0Event->eventId() << " " << picoDst->event()->eventId();
  if(mPicoD0Event->runId() != picoDst->event()->runId() ||
      mPicoD0Event->eventId() != picoDst->event()->eventId())
  {
    LOG_ERROR <<" StPicoD0AnaMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!"<<endm;
    LOG_ERROR <<" StPicoD0AnaMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoEvent and StPicoD0Event are not in sync."<<endm;
    exit(1);
  }

  // I will check if runID is included in array EnergyBadRunList, if yes, I will skip the run. I am not going to use IsBadRun


  // -------------- USER ANALYSIS -------------------------
  TClonesArray const * aKaonPion = mPicoD0Event->kaonPionArray();

  StThreeVectorF pVtx(-999.,-999.,-999.);


  StPicoEvent *event = (StPicoEvent *)picoDst->event();

  
  if(!(isGoodEvent(mYear)))//minBias trigger requires
  {
     return kStOK;
  }

  if(!mGRefMultCorrUtil) {
    LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
    return kStWarn;
  }


  if (mYear==2016 && IsBadEnergyRun(mPicoD0Event->runId())) return kStOK;

  mGRefMultCorrUtil->init(picoDst->event()->runId());
  mGRefMultCorrUtil->initEvent(picoDst->event()->grefMult(),pVtx.z(),picoDst->event()->ZDCx()) ;

  if (mGRefMultCorrUtil->isBadRun(picoDst->event()->runId()))
  {
  //cout<<"This is a bad run from mGRefMultCorrUtil! Skip! " << endl;
  return kStOK;
  }
    
  int centrality  = mGRefMultCorrUtil->getCentralityBin9();
  if(centrality<0) {
    LOG_WARN << "not minBias sample!" << endl;
    return kStOK;
  }
     pVtx = StThreeVectorF(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z());
     vtxz->Fill(pVtx.z());
     vtxr->Fill(pVtx.x(),pVtx.y());
     NEvent->Fill(0);
     hcentr->Fill(centrality);
// << "runID:" << picoDst->event()->runId() << endl;
//---------------------------------
int eventID = mPicoD0Event->eventId();
int RunId = mPicoD0Event->runId();
UInt_t nTracks = picoDst->numberOfTracks();
bool IsThereD0 = false;
std::vector<double> DaughterPionTrackVector;
std::vector<double> DaughterKaonTrackVector;

  // Vypsani hlavicky
  fastjet::ClusterSequence::print_banner();
std::vector<FourMomentum> D0_fourmomentum;
//----------------------------------  


  double reweight = mGRefMultCorrUtil->getWeight();
  //cout << "Number: " << aKaonPion->GetEntries() << endl;
  //cout << "EventId: " << mPicoD0Event->eventId() << endl;

  for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx)
  {
    StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idx);
    StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
    StPicoTrack const* pion = picoDst->track(kp->pionIdx());

    if (!isGoodTrack(kaon) || !isGoodTrack(pion)) continue;
    //cout << "Is good track" << endl;
    if (!isTpcPion(pion)) continue;
    //cout << "Is TPC pion" << endl;
    bool tpcKaon = isTpcKaon(kaon,&pVtx);
    float kBeta = getTofBeta(kaon,&pVtx,picoDst);
    bool tofAvailable = kBeta>0;
    bool tofKaon = tofAvailable && isTofKaon(kaon,kBeta);
    bool goodKaon = (tofAvailable && tofKaon) || (!tofAvailable && tpcKaon);
    if(!goodKaon) continue;
    //cout << "Is good kaon" << endl;

    int charge=0;
    float d0Pt = kp->pt();  
    //cout << "p_{T}: " << d0Pt << endl;
    double dMass = kp->m();
    //cout << "mass: " << dMass << endl;
    if(d0Pt>10) continue;


    double reweight_eff = 1.;
    //if((charge=isD0Pair150(kp))!=0 )
    double pimass = 0.13957018;
    double kaonmass = 0.493677;
      
    if((charge=isD0PairCentrality_pt(kp,centrality))!=0 )
    {
        /*
      TLorentzVector kaonVec(kaon->gMom().Px(), kaon->gMom().Py(), kaon->gMom().Pz(), sqrt(kaonmass*kaonmass+kaon->gMom().Mag()*kaon->gMom().Mag()));
      TLorentzVector pionVec(pion->gMom().Px(), pion->gMom().Py(), pion->gMom().Pz(), sqrt(pimass*pimass+pion->gMom().Mag()*pion->gMom().Mag()));
      kaonVec.Boost(-kp->Px()/kp->Energy(), -kp->Py()/kp->Energy(), -kp->Pz()/kp->Energy());
      pionVec.Boost(-kp->Px()/kp->Energy(), -kp->Py()/kp->Energy(), -kp->Pz()/kp->Energy());
      double agle=kaonVec.Angle(pionVec.Vect());
      angledistr->Fill(kaonVec.Angle(pionVec.Vect()));
      */
      //Angle_like->
      //Angle_unlike->
      //getCorV2(idx,reweight*reweight_eff);//Fill D-hadron 2PC v2 plots
      if(charge==-1)
      {
        if(dMass>1.81&&dMass<1.91)
        candPt->Fill(d0Pt,d0Pt,reweight*reweight_eff);

        //angledistrUnlike->Fill(agle);

        IsThereD0 =true;

        massPt->Fill(dMass,d0Pt,reweight*reweight_eff);
        //massAngleUnlike->Fill(dMass,agle,reweight*reweight_eff);
        //if(agle>3.141582&&agle<3.141602)
        //massPtacut->Fill(dMass,d0Pt,reweight*reweight_eff);
        D0etaunlike->Fill(kp->eta());
      
        //Tracks4Jets->Fill(RunID, eventID, energy, px, py, pz, d0_daughter, charge);
        //Tracks4Jets->Fill(RunId, eventID, kp->Energy(), kp->Px(), kp->Py(), kp->Pz(), pion->charge(), 0);
        //cout << " " << RunId << " " << eventID << " " << kp->Energy() << " " << kp->Px() << " " << kp->Py() << " " << kp->Pz() << " " << pion->charge() << " " << 0 <<endl;
        DaughterPionTrackVector.push_back(pion->id());
        DaughterKaonTrackVector.push_back(kaon->id());
                              // E,             px,       py,       pz,     D0_antiD0;
        FourMomentum D0_actual = {kp->Energy(), kp->Px(), kp->Py(),  kp->Pz(), pion->charge()};
        D0_fourmomentum.push_back(D0_actual);


        pioneta->Fill(pion->pMom().PseudoRapidity());
        kaoneta->Fill(kaon->pMom().PseudoRapidity());
        //----Histograms---------------------------------------
        int runIndex = mPrescales->runIndex(mPicoD0Event->runId());
        mh1TotalEventsInRun->Fill(runIndex);
        mh1TotalGRefMultInRun->Fill(runIndex,picoDst->event()->grefMult());
        mh1TotalKaonsInRun->Fill(runIndex,mPicoD0Event->nKaons());
        mh1TotalPionsInRun->Fill(runIndex,mPicoD0Event->nPions());
        mh1TotalD0CandidatesInRun->Fill(runIndex,mPicoD0Event->nKaonPion());
        mh2NKaonsVsNPions->Fill(mPicoD0Event->nPions(),mPicoD0Event->nKaons());


        mh2KaonDcaVsPt->Fill(kp->pt(),kp->kaonDca());
        mh2PionDcaVsPt->Fill(kp->pt(),kp->pionDca());
        mh2CosThetaVsPt->Fill(kp->pt(),cos(kp->pointingAngle()));
        mh2DcaDaughtersVsPt->Fill(kp->pt(),kp->dcaDaughters());



    //-------------------------------------------------------

    //Tracks4Jets->Fill(eventID, energy, px, py, pz, d0, daughter);

    



      }
      if(charge>0){
         massPtLike->Fill(dMass,d0Pt,reweight*reweight_eff);
        D0etalike->Fill(kp->eta());
        //angledistrLike->Fill(agle);
        //massAngleLike->Fill(dMass,agle,reweight*reweight_eff);
        //if(agle>3.141582&&agle<3.141602)
        //massPtLikeacut->Fill(dMass,d0Pt,reweight*reweight_eff);
      }
       

    }//D loop
  }

//-------------New-one---------------------


  
  
  double pimass = 0.13957018;
  if(IsThereD0){

      //Fast jet
      vector<fastjet::PseudoJet> input_particles;
      vector<fastjet::PseudoJet> chargedjetTracks;
      vector<fastjet::PseudoJet> neutraljetTracks;
      //----------------------------------------------------------
      double R = 0.4;




    //D_0 and antiD_0 candidates
    for (int nD0 = 0; nD0 < D0_fourmomentum.size(); nD0++) {

      fastjet::PseudoJet pj(D0_fourmomentum[nD0].px,D0_fourmomentum[nD0].py,D0_fourmomentum[nD0].pz,D0_fourmomentum[nD0].E);  
      pj.set_user_index(D0_fourmomentum[nD0].D0_antiD0*2);
      input_particles.push_back(pj);
      //cout << D0_fourmomentum[nD0].px << " " << D0_fourmomentum[nD0].py << " " << D0_fourmomentum[nD0].pz << " " << D0_fourmomentum[nD0].E << endl;
    


      //Adding neutral tracks
//-------------------------------------------------------------------------------------------------------------------------------      
      //fill array Sump with momenta of tracks which are matched to BEMC
      GetCaloTrackMomentum(picoDst,TVector3(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z())); 
   
      for (int iTow = 0; iTow < 4800; iTow++){ //get btow info
        
          StPicoBTowHit *towHit = picoDst->btowHit(iTow);
          vector<int> ids = {0,0,0,0,0,0,0,0,0}; 
          if (!towHit || towHit->isBad()) continue; //if the tower is bad or missing info
          int realtowID = towHit->numericIndex2SoftId(iTow);
          double towE;

          if(mYear==2014) {
              if (BadTowerMap[realtowID]) continue; //exclude bad towers
              towE = GetTowerCalibEnergy(iTow+1); //2014 there was some problem with energy calibration
          }
          if(mYear==2016) {
              if (EnergyBadTowerMap[realtowID]) continue; //exclude bad towers
              towE = towHit->energy();
          }

          towE-= fHadronCorr*Sump[iTow];  //fHadronCorr equals 0 or 1, if we want to correct for hadrons or not
          if (towE < 0) towE = 0;
          StEmcGeom* mEmcGeom;
          mEmcGeom = StEmcGeom::getEmcGeom("bemc");
              
          float Toweta_tmp = 0, Towphi = 0;
          mEmcGeom->getEtaPhi(realtowID,Toweta_tmp,Towphi);
          float Toweta = vertexCorrectedEta(Toweta_tmp, event->primaryVertex().z()); //max eta 1.05258 max difference: ET = 0.124452 for E = 0.2, if we cut on |Vz| < 30 cm
          double ET = towE/cosh(Toweta);

          if (ET > 30) {
              TowArr.clear();
              TowEta.clear();
              TowPhi.clear();
              Clusters.clear();
              return kStOK; //discard events with E > 30 GeV towers
              }
              //no clustering
              double px,py,pz;
              //px = towE*cos(Towphi)/cosh(Toweta);
              //py = towE*sin(Towphi)/cosh(Toweta);
              px = ET*cos(Towphi);
              py = ET*sin(Towphi);
              pz = towE*tanh(Toweta);

              PseudoJet inputTower(px, py, pz, towE);
              //cout << "px: " << px << " py: " << py << " pz: " << pz << " towE: " << towE << endl; 
              if (inputTower.perp() > fETmincut){
                //inputTower.set_user_index(0); //default index is -1, 0 means neutral particle
                neutraljetTracks.push_back(inputTower);
                input_particles.push_back(inputTower);
              }
          }

        

//-----------------------------------------------------------------------------------------------------------------------
      //Adding charged tracks
      for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack){
    
      //The i-th track is loaded
      StPicoTrack* trk = picoDst->track(iTrack);


      if(!trk) continue;
      if (!isGoodTrack(trk)) continue;

      double pT = trk->pMom().Perp();
      if(pT != pT) continue; // NaN test.
      float eta = trk->pMom().PseudoRapidity();
      float phi = trk->pMom().Phi();
      float dca = (TVector3(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z()) - trk->origin()).Mag();
      float charged = trk->charge();

        //if the track is not daughter pion nor kion then...
      if (DaughterPionTrackVector[nD0] != trk->id() && DaughterKaonTrackVector[nD0] != trk->id()){
                              //      px,       py,               pz,                                 E, 
      fastjet::PseudoJet pj(trk->pMom().x(),trk->pMom().y(),trk->pMom().z(), sqrt(trk->pMom().Mag()*trk->pMom().Mag()+pimass*pimass));
      //pj.set_user_index(0);
      input_particles.push_back(pj);
      chargedjetTracks.push_back(pj);
      }//End of if (DaughterPionTrackVector[nD0]...


      hpT_tr->Fill(pT, reweight);
      heta_phi_tr->Fill(phi + TMath::Pi(), eta,  reweight);
      heta_tr->Fill(eta, reweight);
      hphi_tr->Fill(phi + TMath::Pi(), reweight); //to shift by pi
      hdca_z_tr->Fill(dca, event->primaryVertex().z(), reweight);
      hdca_pT->Fill(dca, pT, reweight);
      hdca_tr->Fill(dca, reweight);
      hcharged_tr->Fill(charged, reweight);


      
    } //End of track loop
    // background estimation
    JetDefinition jet_def_bkgd(kt_algorithm, R);
    AreaDefinition area_def_bkgd(active_area_explicit_ghosts,GhostedAreaSpec(fGhostMaxrap, 1, 0.01));
    if (centrality == 0 || centrality == 1) nJetsRemove = 2;//remove two hardest jets in central collisions, one in others
    Selector selector = SelectorAbsEtaMax(1.0) * (!SelectorNHardest(nJetsRemove)) * SelectorPtMin(0.01);
    JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);
    bkgd_estimator.set_particles(chargedjetTracks);

    float rho = bkgd_estimator.rho();
    float rho_sigma = bkgd_estimator.sigma();


    //----------------------------------------------------------
    //full jet
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
    AreaDefinition area_def(active_area_explicit_ghosts, GhostedAreaSpec(fGhostMaxrap, 1, 0.01));
    ClusterSequenceArea clust_seq_hard(input_particles, jet_def, area_def);



       /* // run the jet clustering with the above jet definition
        //----------------------------------------------------------
        fastjet::ClusterSequence clust_seq(input_particles, jet_def);
    */
        // get the resulting jets ordered in pt
    //----------------------------------------------------------
    //double ptmin = 5.0;
    double ptmin = 0.0;

    vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq_hard.inclusive_jets(ptmin));
    //To print results
    cout << "Ran " << jet_def.description() << endl;
    printf("%5s %15s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt", "index");
     
    for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
      int user_index = 0;
      double Delta_R_D0 = 0;
      double lambda_alpha = 1;
      double lambda_kappa = 1;
      double zet = 0;
      double lambda = 0;
      double pT_jet = inclusive_jets[i].perp();
      double pT_jet_corr = pT_jet - rho * inclusive_jets[i].area();
      double px_jet = inclusive_jets[i].px();
      double px_jet_corr = px_jet - rho * inclusive_jets[i].area_4vector().px();
      double py_jet = inclusive_jets[i].py();
      double py_jet_corr = py_jet - rho * inclusive_jets[i].area_4vector().py();
      const vector<fastjet::PseudoJet>& constituents = inclusive_jets[i].constituents();

        for (vector<fastjet::PseudoJet>::const_iterator particle = constituents.begin(); particle != constituents.end(); ++particle) {
  
          int index = particle->user_index(); // zjisteni indexu castice v rekonstruovanem PseudoJet
          //cout << "index: " << index << endl;
          double Delta_R =delta_R(inclusive_jets[i].eta(),inclusive_jets[i].phi(),particle->eta(),particle->phi());
          lambda+=pow(particle->pt()/pT_jet_corr,lambda_kappa)*pow( Delta_R /jet_def.R() ,lambda_alpha);
           //Is there D0 in this Jet?
          if (abs(index) == 2 ) {
             user_index=index;
             Delta_R_D0 = Delta_R;
             // z = pT(D0)*^pT(jet)/|pT(jet)|     
             zet = (D0_fourmomentum[nD0].px*px_jet_corr+D0_fourmomentum[nD0].py*py_jet_corr)/(pT_jet_corr*pT_jet_corr);
           }
        } // end loop over jet constituents


      if (abs(user_index) ==2){


      double D0mass = sqrt(D0_fourmomentum[nD0].E*D0_fourmomentum[nD0].E-D0_fourmomentum[nD0].px*D0_fourmomentum[nD0].px-D0_fourmomentum[nD0].py*D0_fourmomentum[nD0].py-D0_fourmomentum[nD0].pz*D0_fourmomentum[nD0].pz);
      double D0_pT = sqrt(D0_fourmomentum[nD0].px*D0_fourmomentum[nD0].px+D0_fourmomentum[nD0].py*D0_fourmomentum[nD0].py);
      //eventID, RunId, inclusive_jets.size(),inclusive_jets[i].rap(), inclusive_jets[i].phi(),
      // inclusive_jets[i].perp(), user_index,D0mass,D0_r,D0_pT,lambda_1_1,z, E
      Jets->Fill( eventID,                                                 
                  RunId, 
                  D0_fourmomentum.size(),
                  inclusive_jets[i].rap(), 
                  inclusive_jets[i].phi(),
                  pT_jet_corr,
                  user_index, 
                  D0mass,
                  Delta_R_D0,
                  D0_pT,
                  lambda,
                  zet,
                  inclusive_jets[i].E()
                );  


        printf("\033[32m%5u %15.8f %15.8f %15.8f %15d\033[0m\n", i, inclusive_jets[i].rap(), inclusive_jets[i].phi(), inclusive_jets[i].perp(), user_index);
      } else{
        printf("%5u %15.8f %15.8f %15.8f %15d\n", i, inclusive_jets[i].rap(), inclusive_jets[i].phi(), inclusive_jets[i].perp(), user_index);
      }
 
  }// end of inclusive jet loop
    //Delete vectors to avoid a pileup
    inclusive_jets.clear();
 cout << "----------------------------------------------" << endl;

    input_particles.clear();
    chargedjetTracks.clear();
    neutraljetTracks.clear();

  }//end of D0 loop


  }//End of IsthereD0 condition

  


//----------------------------------------
DaughterPionTrackVector.clear();
DaughterKaonTrackVector.clear();
//--------------------------------
  return kStOK;
}
//-----------------------------------------------------------------------------

bool StPicoD0AnaMaker::isGoodPair(StKaonPion const* const kp) const
{
  if(!kp) return false;

  StPicoTrack const* kaon = mPicoDstMaker->picoDst()->track(kp->kaonIdx());
  StPicoTrack const* pion = mPicoDstMaker->picoDst()->track(kp->pionIdx());

  StPicoEvent *event = (StPicoEvent *)picoDst->event();
  StThreeVectorF pVtx(-999.,-999.,-999.);
  //pVtx = event->primaryVertex();
  pVtx = StThreeVectorF(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z());

  int charge = kaon->charge() * pion->charge();
  bool pairCuts = kp->m()>1.6 && kp->m()<2.1 &&
    charge==-1;

  return (isTpcKaon(kaon,&pVtx) && isTpcPion(pion) && 
      pairCuts);
}
//----------------------------------------------------------------------------- 
//Correct tower eta for Vz position
//----------------------------------------------------------------------------- 
Double_t StPicoD0AnaMaker::vertexCorrectedEta(double eta, double vz) {
    double tower_theta = 2.0 * atan(exp(-eta));
    double z = 0.0;
    if (eta != 0.0) z = mBarrelRadius / tan(tower_theta);
    double z_diff = z - vz;
    double theta_corr = atan2(mBarrelRadius, z_diff);
    double eta_corr = -log(tan(theta_corr / 2.0));
    return eta_corr;
}
//--------------------------------------------------------------
Bool_t StPicoD0AnaMaker::GetCaloTrackMomentum(StPicoDst *mPicoDst, TVector3 mPrimVtx) {
  //loop over global tracks  - towers

  UInt_t nTracks = mPicoDst->numberOfTracks();

  for (unsigned int itrack = 0; itrack < nTracks; itrack++) {
        StPicoTrack *trk = mPicoDst->track(itrack);
        TVector3 gMom = trk->gMom();
        
        //using global tracks
        double pT = gMom.Perp();
        if(pT != pT || pT < 0.2) continue;
        float eta = gMom.PseudoRapidity();
        if (fabs(eta) > 1) continue;
        float phi = gMom.Phi();
        float nHitsFit = trk->nHitsFit();
        float nHitsMax = trk->nHitsMax();
        if (nHitsFit < 15 || nHitsFit/nHitsMax < 0.52) continue; //some basic QA cuts
        double Bfield = mPicoDst->event()->bField();
        StPicoPhysicalHelix trkhelix = trk->helix(Bfield);
        float vtx_x = mPrimVtx.x();
        float vtx_y = mPrimVtx.y();
        float vtx_z = mPrimVtx.z();

        TVector3 dcaPoint = trkhelix.at(trkhelix.pathLength(vtx_x, vtx_y));
        float dca_z = dcaPoint.z() - vtx_z; //check
        if (fabs(dca_z) > maxdcazhadroncorr) continue; 
        int TowIndex = -99999;
        TowIndex = trk->bemcTowerIndex();
  
        float p = 0;
        
        if (TowIndex > 0) {
          p = gMom.Mag();
          Sump[TowIndex-1] += p;
          //cout << p << endl;
        }
  
  }// END global track loop
  
  return true;
}
//---------------------------------------------------------------------------
Double_t StPicoD0AnaMaker::GetTowerCalibEnergy(Int_t TowerId)
{

  StPicoBTowHit *tower = static_cast<StPicoBTowHit*>(picoDst->btowHit(TowerId-1));
  Float_t pedestal, rms;
  Int_t status;
  mTables->getPedestal(BTOW, TowerId, 0, pedestal, rms);
  mTables->getStatus(BTOW, TowerId, status);  
  Double_t *TowerCoeff;

  //Only for 2014

  if(picoDst->event()->runId() <= 15094020) TowerCoeff = CPre;

  else TowerCoeff = CLowMidHigh;


  TowerCoeff = CLowMidHigh;

  Double_t calibEnergy = TowerCoeff[TowerId-1]*(tower->adc() - pedestal);
  return calibEnergy;
}
//---------------------------------------------------------------------------
int StPicoD0AnaMaker::isD0Pair(StKaonPion const* const kp) const
{

  StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
  StPicoTrack const* pion = picoDst->track(kp->pionIdx());
  TLorentzVector d0Lorentz;
  d0Lorentz.SetPtEtaPhiM(kp->pt(),kp->eta(),kp->phi(),kp->m());

  if(fabs(d0Lorentz.Rapidity())>1.) return 0;

  bool pairCuts = false;
  if(kp->pt()<1)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0061 &&
      kp->pionDca() > 0.0110 && kp->kaonDca() > 0.0103 &&
      kp->dcaDaughters() < 0.0084 && kp->decayLength()>0.0145;  
  }
  else if(kp->pt()<2)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0049 &&
      kp->pionDca() > 0.0111 && kp->kaonDca() > 0.0091 &&
      kp->dcaDaughters() < 0.0066 && kp->decayLength()>0.0181;  
  }
  else if(kp->pt()<3)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0038 &&
      kp->pionDca() > 0.0086 && kp->kaonDca() > 0.0095 &&
      kp->dcaDaughters() < 0.0057 && kp->decayLength()>0.0212;  
  }
  else if(kp->pt()<5)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0038 &&
      kp->pionDca() > 0.0081 && kp->kaonDca() > 0.0079 &&
      kp->dcaDaughters() < 0.0050 && kp->decayLength()>0.0247;  
  }
  else 
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0040 &&
      kp->pionDca() > 0.0062 && kp->kaonDca() > 0.0058 &&
      kp->dcaDaughters() < 0.0060 && kp->decayLength()>0.0259;  
  }

  int charge = kaon->charge() * pion->charge();
  if(charge>0)
    charge = kaon->charge()>0 ? 1:2;


  if(pairCuts)
    return charge;
  else
    return 0;
}
/*
*/
//---------------------------------------------------------------------------
bool StPicoD0AnaMaker::IsBadEnergyRun(int runID) {
    for (int i = 0; i < sizeof(EnergyBadRunList)/sizeof(EnergyBadRunList[0]); i++) {
        if (EnergyBadRunList[i] == runID) {
            return true; // Hodnota se rovná prvku v poli
        }
    }
    return false; // Hodnota se nerovná žádnému prvku v poli
}

//---------------------------------------------------------------------------
int StPicoD0AnaMaker::isD0Pair50(StKaonPion const* const kp) const
{

  StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
  StPicoTrack const* pion = picoDst->track(kp->pionIdx());
  bool pairCuts = false;
  if(kp->pt()<1)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0044 &&
      kp->pionDca() > 0.0120 && kp->kaonDca() > 0.0119 &&
      kp->dcaDaughters() < 0.0069 && kp->decayLength()>0.0144;  
  }
  else if(kp->pt()<2)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0036 &&
      kp->pionDca() > 0.0102 && kp->kaonDca() > 0.0110 &&
      kp->dcaDaughters() < 0.0048 && kp->decayLength()>0.0204;  
  }
  else if(kp->pt()<3)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0031 &&
      kp->pionDca() > 0.0118 && kp->kaonDca() > 0.0109 &&
      kp->dcaDaughters() < 0.0044 && kp->decayLength()>0.0242;  
  }
  else if(kp->pt()<5)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0026 &&
      kp->pionDca() > 0.0109 && kp->kaonDca() > 0.0106 &&
      kp->dcaDaughters() < 0.0049 && kp->decayLength()>0.0245;  
  }
  else 
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0032 &&
      kp->pionDca() > 0.0096 && kp->kaonDca() > 0.0080 &&
      kp->dcaDaughters() < 0.0047 && kp->decayLength()>0.0300;  
  }

  int charge = kaon->charge() * pion->charge();
  if(charge>0)
    charge = kaon->charge()>0 ? 1:2;


  if(pairCuts)
    return charge;
  else
    return 0;
}

int StPicoD0AnaMaker::isD0Pair150(StKaonPion const* const kp) const
{

  StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
  StPicoTrack const* pion = picoDst->track(kp->pionIdx());
  bool pairCuts = false;
  if(kp->pt()<1)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0072 &&
      kp->pionDca() > 0.0092 && kp->kaonDca() > 0.0105 &&
      kp->dcaDaughters() < 0.0077 && kp->decayLength()>0.0110;  
  }
  else if(kp->pt()<2)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0053 &&
      kp->pionDca() > 0.0078 && kp->kaonDca() > 0.0068 &&
      kp->dcaDaughters() < 0.0078 && kp->decayLength()>0.0168;  
  }
  else if(kp->pt()<3)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0047 &&
      kp->pionDca() > 0.0086 && kp->kaonDca() > 0.0080 &&
      kp->dcaDaughters() < 0.0074 && kp->decayLength()>0.0187;  
  }
  else if(kp->pt()<5)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0042 &&
      kp->pionDca() > 0.0065 && kp->kaonDca() > 0.0066 &&
      kp->dcaDaughters() < 0.0068 && kp->decayLength()>0.0199;  
  }
  else 
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0062 &&
      kp->pionDca() > 0.0047 && kp->kaonDca() > 0.0041 &&
      kp->dcaDaughters() < 0.0066 && kp->decayLength()>0.0180;  
  }

  int charge = kaon->charge() * pion->charge();
  if(charge>0)
    charge = kaon->charge()>0 ? 1:2;


  if(pairCuts)
    return charge;
  else
    return 0;
}

int StPicoD0AnaMaker::isD0PairCentrality_pt(StKaonPion const* const kp, int Centrality) const
{

  StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
  StPicoTrack const* pion = picoDst->track(kp->pionIdx());
  bool pairCuts = false;
  
/*
  if(kp->pt()<0.5)
  {
  KPMom = 0;
  }
  else if(kp->pt()<1)
  {
  KPMom = 1;
  }
  else if(kp->pt()<2)
  {
  KPMom = 2;
  }
  else if(kp->pt()<3)
  {
  KPMom = 3;
  }
  else if(kp->pt()<5)
  {
  KPMom = 4;
  }
  else
  {
  KPMom = 5;
  }*/

  // Centr.   0-10%       10-20%         20-40%         40-60%        60-80%
  // C_ID     8,7        6               5,4             3,2         1,0
  // bin       0         1                 2               3           4

 int Centrality2 =  (Centrality == 8 || Centrality == 7) ? 0 :
                    (Centrality == 6) ? 1 :
                    (Centrality == 5 || Centrality == 4) ? 2 :
                    (Centrality == 3 || Centrality == 2) ? 3 :
                    (Centrality == 1 || Centrality == 0) ? 4 : -1;


  int KPMom = (kp->pt() < 0.5) ? 0 : (kp->pt() < 1) ? 1 : (kp->pt() < 2) ? 2 : (kp->pt() < 3) ? 3 : (kp->pt() < 5) ? 4 : 5;

  pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < mycuts::DCA_D0_cut[KPMom][Centrality2] &&
  kp->pionDca() > mycuts::pionDCA_cut[KPMom][Centrality2] && kp->kaonDca() > mycuts::kaonDCA_cut[KPMom][Centrality2] &&
  kp->dcaDaughters() < mycuts::pionkaonDCA_cut[KPMom][Centrality2] && kp->decayLength()> mycuts::D0_decayLength_cut[KPMom][Centrality2];  

  int charge = kaon->charge() * pion->charge();
  if(charge>0)
    charge = kaon->charge()>0 ? 1:2;


  if(pairCuts)
    return charge;
  else
    return 0;
}

int StPicoD0AnaMaker::isD0PairOld(StKaonPion const* const kp) const
{

  StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
  StPicoTrack const* pion = picoDst->track(kp->pionIdx());
  bool pairCuts = false;
  if(kp->pt()<1)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0062 &&
      kp->pionDca() > 0.0109 && kp->kaonDca() > 0.0123 &&
      kp->dcaDaughters() < 0.0082 && kp->decayLength()>0.0149;  
  }
  else if(kp->pt()<2)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0047 &&
      kp->pionDca() > 0.0108 && kp->kaonDca() > 0.0097 &&
      kp->dcaDaughters() < 0.0070 && kp->decayLength()>0.0205;  
  }
  else if(kp->pt()<3)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0040 &&
      kp->pionDca() > 0.0100 && kp->kaonDca() > 0.0091 &&
      kp->dcaDaughters() < 0.0056 && kp->decayLength()>0.0216;  
  }
  else if(kp->pt()<5)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0041 &&
      kp->pionDca() > 0.0074 && kp->kaonDca() > 0.0075 &&
      kp->dcaDaughters() < 0.0065 && kp->decayLength()>0.0233;  
  }
  else 
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0042 &&
      kp->pionDca() > 0.067 && kp->kaonDca() > 0.0053 &&
      kp->dcaDaughters() < 0.0065 && kp->decayLength()>0.0282;  
  }

  int charge = kaon->charge() * pion->charge();


  if(pairCuts)
    return charge;
  else
    return 0;
}



bool StPicoD0AnaMaker::isGoodEvent(int mYear)
{
  StPicoEvent *event = (StPicoEvent *)picoDst->event();
  // return (event->triggerWord() & mycuts::triggerWord) &&
  return (isMBTrigger(mYear) &&
      sqrt(event->primaryVertex().x()*event->primaryVertex().x()+event->primaryVertex().y()*event->primaryVertex().y()) < mycuts::vz &&
      fabs(event->primaryVertex().z()) < mycuts::vz &&
      fabs(event->primaryVertex().z() - event->vzVpd()) < mycuts::vzVpdVz);
  //  return event->triggerWord() & mycuts::triggerWord;
}
bool StPicoD0AnaMaker::isMBTrigger(int mYear)
{
 const std::set<int>* mbTriggers = nullptr;

      if(mYear ==2016) mbTriggers = &mycuts::mbTriggers2016;
      if(mYear ==2014) mbTriggers = &mycuts::mbTriggers2014;

      StPicoEvent* event = static_cast<StPicoEvent*>(mPicoDstMaker->picoDst()->event());
      return std::any_of(mbTriggers->begin(), mbTriggers->end(), [&](int trigger) { return event->isTrigger(trigger); });

}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodTrack(StPicoTrack const * const trk) const
{
  // Require at least one hit on every layer of PXL and IST.
  // It is done here for tests on the preview II data.
  // The new StPicoTrack which is used in official production has a method to check this
  //return trk->gPt() > mycuts::minPt && trk->nHitsFit() >= mycuts::nHitsFit && trk->isHFTTrack();

  bool HFTCondition = (trk->hasPxl1Hit() && trk->hasPxl2Hit()) && (trk->hasSstHit() || trk->hasIstHit());
  return trk->gPt() > mycuts::minPt && trk->nHitsFit() >= mycuts::nHitsFit && HFTCondition;
  //return  trk->nHitsFit() >= mycuts::nHitsFit;
}
bool StPicoD0AnaMaker::isGoodTrack2(StPicoTrack const * const trk) const
{

  return trk->gPt() > mycuts::minPt && trk->nHitsFit() >= mycuts::nHitsFit;

}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodHadron(StPicoTrack const * const trk) const
{
  //return trk->pMom().perp() > mycuts::hadronPtMin &&trk->pMom().perp() < mycuts::hadronPtMax && trk->nHitsFit() >= mycuts::nHitsFit &&fabs(trk->pMom().pseudoRapidity())<1.&&fabs(trk->nSigmaElectron())>3 && (1.0*trk->nHitsFit()/trk->nHitsMax())>0.52;
  return trk->pMom().Perp() > mycuts::hadronPtMin &&trk->pMom().Perp() < mycuts::hadronPtMax && trk->nHitsFit() >= 15 &&fabs(trk->pMom().PseudoRapidity())<1. && (1.0*trk->nHitsFit()/trk->nHitsMax())>0.52;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTpcPion(StPicoTrack const * const trk) const
{
  return fabs(trk->nSigmaPion()) < mycuts::nSigmaPion;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTpcKaon(StPicoTrack const * const trk, StThreeVectorF const* const pVtx) const
{
  return fabs(trk->nSigmaKaon()) < mycuts::nSigmaKaon;
  //      || tofKaon;
}
//-----------------------------------------------------------------------------
float StPicoD0AnaMaker::getTofBeta(StPicoTrack const * const trk, StThreeVectorF const* const pVtx,StPicoDst const* const picoDst) const
{

  int index2tof = trk->bTofPidTraitsIndex();

  float beta = std::numeric_limits<float>::quiet_NaN();

  if(index2tof >= 0)
  {
    StPicoBTofPidTraits *tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);

    if(tofPid)
    {
      beta = tofPid->btofBeta();

      if (beta < 1e-4)
      {
        StThreeVectorF const btofHitPos = StThreeVectorF(tofPid->btofHitPos().x(),tofPid->btofHitPos().y(),tofPid->btofHitPos().z());
        //StPhysicalHelixD helix = trk->helix();
        StPicoPhysicalHelix helix = trk->helix(picoDst->event()->bField());

        float L = tofPathLength(pVtx, &btofHitPos, helix.curvature());
        float tof = tofPid->btof();
        if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
        else beta = std::numeric_limits<float>::quiet_NaN();
      }
    }
  }

  return beta;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTofKaon(StPicoTrack const * const trk, float beta) const
{
  bool tofKaon = false;

  if(beta>0)
  {
    //double ptot = trk->dcaGeometry().momentum().mag();
    double ptot = trk->gMom().Mag();
    float beta_k = ptot/sqrt(ptot*ptot+M_KAON_PLUS*M_KAON_PLUS);
    tofKaon = fabs(1/beta - 1/beta_k) < mycuts::kTofBetaDiff ? true : false;
  }

  return tofKaon;
}


bool StPicoD0AnaMaker::getCorHadron(float eta,vector<float> &hadronsPhi, vector<unsigned int> index1, float phi, float etaCut) 
{
  for(unsigned int i=0;i<picoDst->numberOfTracks();++i)
  {
    StPicoTrack const* hadron = picoDst->track(i);
    if(!hadron)  continue; 
    if(hadron->pMom().Perp()<0.2) continue;
    etaPhi_Hadron_all->Fill(hadron->pMom().Phi(),hadron->pMom().PseudoRapidity());
    vector<unsigned int>::iterator it_index;
    it_index = find(index1.begin(),index1.end(),i);
    if(it_index!=index1.end())  continue;
    if(!isGoodHadron(hadron)) continue;
    float dEta = fabs(hadron->pMom().PseudoRapidity() - eta);
    float dPhi = (hadron->pMom().Phi() - phi);
    if(etaCut<0.001)
    {
      dEtaDHadron->Fill(dEta);
      hEtaD->Fill(eta);
      hEtaHadron->Fill(hadron->pMom().PseudoRapidity());
    }
    //if(dPhi>3.1416) dPhi = 2*3.1416-dPhi;
    if(dEta< etaCut|| dEta > mycuts::corDetaMax)  continue;
    etaPhi->Fill(dPhi,dEta);
    etaPhi_Hadron->Fill(hadron->pMom().Phi(),hadron->pMom().PseudoRapidity());
    hadronsPhi.push_back(hadron->pMom().Phi());
  }
  //  fixPhi(hadronsPhi);
  return true;

}

float StPicoD0AnaMaker::sumCos(float phi,vector<float> &hadronsPhi) 
{
  float sumOfCos = 0;
  for(unsigned int i=0;i<hadronsPhi.size();++i)
  {
    sumOfCos += cos(2*(phi-hadronsPhi[i]));
  }
  return sumOfCos;
}

bool StPicoD0AnaMaker::fixPhi(vector<float> &phi) 
{
  if(phi.size() == 0) return false;
  float sumPhi = 0;
  for(unsigned int i=0;i<phi.size();i++)
    sumPhi+=phi[i];
  float meanPhi = sumPhi/phi.size();
  for(unsigned int i=0;i<phi.size();i++)
    phi[i] = phi[i]-meanPhi;  
  return true;
}


bool StPicoD0AnaMaker::isEtaGap(double dEta,double mGap,double hEta)
{
  if(mGap == 0) return true;
  //double range =  2. - mGap*2;
  // if(dEta> (1.-2*mGap))
  //   return hEta<(dEta-mGap) && hEta>(dEta-mGap-range);
  // else if(dEta<(-1.+2*mGap))
  //   return hEta>(dEta+mGap) && hEta<(dEta+mGap+range);
  // else 
  //   return (hEta>(dEta+mGap) || hEta<(dEta-mGap));
  if(dEta>0)
    return hEta<-mGap;
  else
    return hEta>mGap;
}







