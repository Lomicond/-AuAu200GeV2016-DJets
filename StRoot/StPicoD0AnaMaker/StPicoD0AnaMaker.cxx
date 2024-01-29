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
using namespace std;
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

  //The prescales are used only for rescaling the RunID to the lower numbers
  mPrescales = new StPicoPrescales(mycuts::prescalesFilesDirectoryName);

  //Loading the information from the D0EventMaker
  mChain->GetBranch("dEvent")->SetAutoDelete(kFALSE);
  mChain->SetBranchAddress("dEvent", &mPicoD0Event);

  //Output file
  mOutputFile = new TFile(mOutFileName.Data(), "RECREATE");
  mOutputFile->cd();

  //Number of Runs
  int nRuns = mPrescales->numberOfRuns();

  //Events histograms:
  vtxz = new TH1F("vtxz",";PVtx.z() [cm]; Count",100,-10,10);
  vtxr = new TH2D("vtxr",";PVtx.x() [cm]; PVtx.y() [cm]",100,-3,3,100,-3,3);
  hcentr = new TH1D("hcentr",";C_{ID}",9,-0.5,8.5);
  NEventsCuts = new TH1F("NEventsCuts", "NEventsCuts;Cuts;Count", 9, 0, 9);
      NEventsCuts->GetXaxis()->SetBinLabel(1, "All events");
      NEventsCuts->GetXaxis()->SetBinLabel(2, "Triggers");
      NEventsCuts->GetXaxis()->SetBinLabel(3, "V_{r}");
      NEventsCuts->GetXaxis()->SetBinLabel(4, "V_{z}");
      NEventsCuts->GetXaxis()->SetBinLabel(5, "|V_{z} - V_{z}^{VPD}|");
      NEventsCuts->GetXaxis()->SetBinLabel(6, "Good E run");
      NEventsCuts->GetXaxis()->SetBinLabel(7, "Centrality");
      NEventsCuts->GetXaxis()->SetBinLabel(8, "Good D^{0}");
      NEventsCuts->GetXaxis()->SetBinLabel(9, "Cal.: E_{T} < 30 GeV");
  mh1TotalEventsInRun = new TH1F("mh1TotalEventsInRun","totalEventsInRun;runIndex;totalEventsInRun",nRuns+1,0,nRuns+1);
  mh1TotalGRefMultInRun = new TH1F("mh1TotalGRefMultInRun","totalGRefMultInRun;runIndex;totalGRefMultInRun",nRuns+1,0,nRuns+1);

  //D0 histograms:
  D0etalike = new TH1D("D0etalike",";#eta; Count",100,-5,5);
  D0etaunlike = new TH1D("D0etaunlike",";#eta; Count",100,-5,5);
  pioneta = new TH1D("pioneta",";#eta; Count",100,-5,5);
  kaoneta =  new TH1D("kaoneta",";#eta; Count",100,-5,5);
  mh1TotalKaonsInRun = new TH1F("mh1TotalKaonsInRun","totalKaonsInRun;runIndex;totalKaonsInRun",nRuns+1,0,nRuns+1);
  mh1TotalPionsInRun = new TH1F("mh1TotalPionsInRun","totalPionsInRun;runIndex;totalPionsInRun",nRuns+1,0,nRuns+1);
  mh1TotalD0CandidatesInRun = new TH1F("mh1TotalD0CandidatesInRun","totalD0CandidatesInRun;runIndex;totalD0CandidatesInRun",nRuns+1,0,nRuns+1);
  mh2NKaonsVsNPions = new TH2F("mh2NKaonsVsNPions","nKaonsVsNPions;nPions;nKaons",1000,0,1000,300,0,300);
  mh2KaonDcaVsPt = new TH2F("mh2KaonDcaVsPt","kaonDcaVsPt;p_{T}(K#pi)(GeV/c);K DCA(cm)",120,0,12,50,0,0.05);
  mh2PionDcaVsPt = new TH2F("mh2PionDcaVsPt","pionDcaVsPt;p_{T}(K#pi)(GeV/c);#pi DCA(cm)",120,0,12,50,0,0.05);
  mh2CosThetaVsPt = new TH2F("mh2CosThetaVsPt","cosThetaVsPt;p_{T}(K#pi)(GeV/c);cos(#theta)",120,0,12,500,0,1.0);
  mh2DcaDaughtersVsPt = new TH2F("mh2DcaDaughtersVsPt","dcaDaughtersVsPt;p_{T}(K#pi)(GeV/c);dcaDaughters(cm)",120,0,12,200,0,0.02);

  //Jet tracks
  JetTracksdEdx = new TH2D("JetTracksdEdx","JetTracksdEdx;charge #times |p| [GeV/c]; dE/dx [KeV/cm]",200,-4,4,200,0,8);
  JetTracksdEdxCut = new TH2D("JetTracksdEdxCut","JetTracksdEdxCut;charge #times |p| [GeV/c]; dE/dx [KeV/cm]",200,-4,4,200,0,8);
  hpT_tr = new TH1F("hpT_tr","hpT_tr;p_{T}(GeV/c);",200,0,100);
  heta_phi_tr = new TH2F("heta_phi_tr","heta_phi_tr;#phi;#eta",120,0,12,50,0,0.05);
  heta_tr = new TH1F("heta_tr","heta_tr;#eta",50,-6,6);
  hphi_tr = new TH1F("hphi_tr","hphi_tr;#phi;",40,-6.2830,6.2830);
  hdca_z_tr = new TH2F("hdca_z_tr","hdca_z_tr;DCA(cm); v_{z}(GeV/c)",120,0,12,50,-10,10);
  hdca_pT = new TH2F("hdca_pT","hdca_pT; DCA(cm); p_{T}(GeV/c) ",120,0,12,100,0,200);
  hdca_tr = new TH1F("hdca_tr","hdca_tr; DCA(cm)",120,0,12);
  hcharged_tr = new TH1F("hcharged_tr","hcharged_tr;Charge;",5,-1,1);

  //Jet background
  Jet_grefmult_pt_background = new TH2D("Jet_rho_vs_grefmult","Jet_rho_vs_grefmult;grefmult; p_{T} (GeV/c)",1000,0,1000,100,-1,32);
  Jet_D0pT_vs_D0rapidity = new TH2D("Jet_D0pT_vs_D0rapidity","Jet_D0pT_vs_D0rapidity;rapidity; p_{T} (GeV/c)",100,-1.5,1.5,100,0,10);
  Jet_D0pT_vs_Jetrapidity = new TH2D("Jet_D0pT_vs_JetRapidity","Jet_D0pT_vs_JetRapidity;rapidity; p_{T} (GeV/c)",100,-1.5,1.5,100,0,10);
  Jet_phi = new TH1D("Jet_phi","Jet_phi;#phi;",40,-6.2830,6.2830);

  //TNtuple D0-jets
  Jets = new TNtuple("Jets", "Jets", "RunId:centrality:centr_weight:NJet:pseudorapidity:jet_pt:jet_pt_corr:D0mass:D0_r:D0_pT:lambda_1_1:lambda_1_1half:lambda_1_2:lambda_1_3:z");
      //RunId: sgn(RunID) = 1 -> D0, sgn(RunID) = -1 -> antiD0
      //Centrality: 0 -> 70-80%, 1 -> 60-70%, 2 -> 50-60%, 3 -> 40-50%, 4 -> 30-40%, 5 -> 20-30%, 6 -> 10-20%, 7 -> 5-10%, 8 -> 0-5%
      //Centr_weight: centrality weight
      //NJet: Number of D0 in one events.
      //pseudorapidity: pseudorapidity of D0-jet
      //jet_pt: pt of D0-jet
      //jet_pt_corr: pt of D0-jet after the background subtraction
      //D0mass: mass of D0
      //D0_r: r of D0
      //D0_pT: pt of D0
      //lambda_1_1: angularity kappa=1 and alpha=1 after the background subtraction
      //lambda_1_1half: angularity kappa=1 and alpha=1.5 after the background subtraction
      //lambda_1_2: angularity kappa=1 and alpha=2 after the background subtraction
      //lambda_1_3: angularity kappa=1 and alpha=3 after the background subtraction
      //z: z of D0-jet after the background subtraction

  //TNtuple additional info
  Jets2 = new TNtuple("Jets2", "Jets2", "centrality:centr_weight:jet_phi:grefmult:bg_dens:jet_area:jet_rap");
        // Centrality: centrality
        // centr_weight: centrality weight
        // jet_phi: phi of D0-jet
        // grefmult: grefmult
        // bg_dens: background density
        // Jet_area: area of D0-jet
        // jet_rap: rapidity of D0-jet


    //Bin size of 2D D0 mass-pt like-sign and unlike-sign histograms
  const int xbinSize=100;
  float binMass[2001];

  //Calculating of bin edges
  float xbin[101];
  for(int i=0;i<101;i++) xbin[i] = 0.1*i;
  for(int i=0;i<2001;i++) binMass[i] = 0.01*i;

  //2D D0 mass-pt like-sign and unlike-sign histograms
  massPt = new TH2D("massPt",";M_{K#pi} [GeV/c^{2}]; p_{T} [GeV/c]",2000,binMass,xbinSize,xbin);
  massPtLike = new TH2D("massPtLike",";M_{K#pi} [GeV/c^{2}]; p_{T} [GeV/c]",2000,binMass,xbinSize,xbin);

  //Loading of BEMC tables
  StMaker* maker = GetMaker("Eread");
  mADCtoEMaker = dynamic_cast<StEmcADCtoEMaker*>(GetMaker("Eread"));
  assert(mADCtoEMaker);
  mTables = mADCtoEMaker->getBemcData()->getTables();

  return kStOK;
}

//-----------------------------------------------------------------------------
StPicoD0AnaMaker::~StPicoD0AnaMaker(){
  //Destructor
  delete mGRefMultCorrUtil;
}

//-----------------------------------------------------------------------------
struct FourMomentum {
  //Four-momentum of the reconstructed particle
  double E, px, py, pz, D0_antiD0, D0Mass;
  // D0_antiD0: 1 -> D0, -1 -> antiD0
};

//-----------------------------------------------------------------------------
double delta_R(double eta1, double phi1, double eta2, double phi2) {
  //Calculating of delta R = sqrt(delta eta^2 + delta phi^2)
  double deta = eta1 - eta2;
  double dphi = TVector2::Phi_mpi_pi(phi1 - phi2);

  return sqrt(deta*deta + dphi*dphi);

  //Function returns delta R
}

//-----------------------------------------------------------------------------
Int_t StPicoD0AnaMaker::Finish(){
  //Write histograms and close the output file, if you create a new histogram it has to be added here
  LOG_INFO << " StPicoD0AnaMaker - writing data and closing output file " <<endm;
  fout.close();
  mOutputFile->cd();

  //Events histograms:
  vtxz->Write();
  vtxr->Write();
  hcentr->Write();
  NEventsCuts->Write();
  mh1TotalEventsInRun->Write();
  mh1TotalGRefMultInRun->Write();

  //D0 histograms:
  D0etalike->Write();
  D0etaunlike->Write();
  pioneta->Write();
  kaoneta->Write();
  mh1TotalKaonsInRun->Write();
  mh1TotalPionsInRun->Write();
  mh1TotalD0CandidatesInRun->Write();
  mh2NKaonsVsNPions->Write();
  mh2KaonDcaVsPt->Write();
  mh2PionDcaVsPt->Write();
  mh2CosThetaVsPt->Write();
  mh2DcaDaughtersVsPt->Write();

  //Jet Tracks
  JetTracksdEdx->Write();
  JetTracksdEdxCut->Write();
  hpT_tr->Write();
  heta_phi_tr->Write();
  heta_tr->Write();
  hphi_tr->Write();
  hdca_z_tr->Write();
  hdca_pT->Write();
  hdca_tr->Write();
  hcharged_tr->Write();

  //Jet background
  Jet_grefmult_pt_background->Write();
  Jet_D0pT_vs_D0rapidity->Write();
  Jet_D0pT_vs_Jetrapidity->Write();
  Jet_phi->Write();

  //2D D0 mass-pt like-sign and unlike-sign histograms
  massPt->Write();
  massPtLike->Write();

  //TNtuple D0-jets
  Jets->Write();
  Jets2->Write();

  //Closing of the output file
  mOutputFile->Close();

  //Closing of the input file including prescales
  delete mPrescales;

  return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoD0AnaMaker::Make()
{
  //Main function where the analysis is done, each "event" is analyzed here
  readNextEvent();

  //Check if everything is loaded properly
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

  //Check if the loaded picodsts are consistent with raw D0 reconstructed data from StPicoD0Event
  if(mPicoD0Event->runId() != picoDst->event()->runId() ||
      mPicoD0Event->eventId() != picoDst->event()->eventId())
  {
    LOG_ERROR <<" StPicoD0AnaMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!"<<endm;
    LOG_ERROR <<" StPicoD0AnaMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoEvent and StPicoD0Event are not in sync."<<endm;
    exit(1);
  }

  // -------------- USER ANALYSIS -------------------------

  //Loading of the raw D0 daughters
  TClonesArray const * aKaonPion = mPicoD0Event->kaonPionArray();

  //Vertex of the event
  StThreeVectorF pVtx(-999.,-999.,-999.);

  //Loading of the event
  StPicoEvent *event = (StPicoEvent *)picoDst->event();

  //Loading of the primary vertex
  pVtx = StThreeVectorF(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z());

  //NEventsCuts: All events
  NEventsCuts->Fill(0);

  //Check if the event is good (vertex, pile-up, MB trigger)
  if(!(isGoodEvent(mYear,NEventsCuts))) return kStOK;

  //Check if the GRefMultCorr information is available
  if(!mGRefMultCorrUtil) {
    LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
    return kStWarn;
  }

  //Check if the event is in a bad run (due to BEMC towers)
  if (mYear==2016 && IsBadEnergyRun(mPicoD0Event->runId())) return kStOK;

  //NEventsCuts: Good E
  NEventsCuts->Fill(5);

  //Loadig of the GRefMultCorr information
  mGRefMultCorrUtil->init(picoDst->event()->runId());
  mGRefMultCorrUtil->initEvent(picoDst->event()->grefMult(),pVtx.z(),picoDst->event()->ZDCx()) ;

  //Check if the event is in a bad run (due to centrality estimation)
  if (mGRefMultCorrUtil->isBadRun(picoDst->event()->runId())) return kStOK;

  //Loading of the centrality
  int centrality  = mGRefMultCorrUtil->getCentralityBin9();

  //Check if the centrality is in the range
  if(centrality<0) return kStOK;

  //NEventsCuts: Centrality
  NEventsCuts->Fill(6);

  //Filling events histograms
  vtxz->Fill(pVtx.z());
  vtxr->Fill(pVtx.x(),pVtx.y());
  hcentr->Fill(centrality);

  //Loading event information
  //int eventID = mPicoD0Event->eventId();
  int RunId = mPicoD0Event->runId();

  //Loading of the number of tracks in the event
  UInt_t nTracks = picoDst->numberOfTracks();

  //Variable checking if there is a good D0 candidate in the event, changed to true, if the candidate is found
  bool IsThereD0 = false;

  //Preparation of the vector of the daughter candidates
  std::vector<double> DaughterPionTrackVector;
  std::vector<double> DaughterKaonTrackVector;

  //Print of the Fastjet banner
  fastjet::ClusterSequence::print_banner();

  //Preparation of the four-vector of the D0 candidates
  std::vector<FourMomentum> D0_fourmomentum;

  //Loading of the weight for the centrality correction
  double reweight = mGRefMultCorrUtil->getWeight();

  //For cyklus of all raw D0 candidates
  for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx){

    //Loading of the raw D0 candidates and there daughters
    StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idx);
    StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
    StPicoTrack const* pion = picoDst->track(kp->pionIdx());

//----------TPC-and-TOF-identification-of-the-daughter-tracks----------
    //Check if the daughter tracks are good.
    if (!isGoodTrack(kaon) || !isGoodTrack(pion)) continue;

    //Check if the daughter pions tracks are good.
    if (!isTpcPion(pion)) continue;

    //Check if the daughter kaon tracks are good.
    bool tpcKaon = isTpcKaon(kaon,&pVtx);

    //Calculating of the beta for the kaon (TOF)
    float kBeta = getTofBeta(kaon,&pVtx,picoDst);

    //Check if the there is a TOF information for the kaon
    bool tofAvailable = kBeta>0;

    //Check if the kaon is good w.r.t the expected beta value (TOF)
    bool tofKaon = tofAvailable && isTofKaon(kaon,kBeta);

    //Final check if the kaon is good. If the TOF information is not available the more strict TPC check is used.
    bool goodKaon = (tofAvailable && tofKaon) || (!tofAvailable && tpcKaon);

    //Check if the kaon is good
    if(!goodKaon) continue;

//----------Pair-cuts--------------------------------------------------
    //Initialisation of the charge
    int charge=0;

    //Loading of the D0 charge and mass
    float d0Pt = kp->pt();  
    double dMass = kp->m();

    //Upper cut on the D0 pT
    if(d0Pt>10) continue;

    //Initialisation of the centrality weight
    double reweight_eff = 1.;

    //Check if all pair cuts conditions are met
    if((charge=isD0PairCentrality_pt(kp,centrality, mYear))!=0 ){
        //Charge = -1 -> Unlike-sign (pi+K- or pi-K+), Charge = 1 -> Like-sign (pi+K+), Charge = 2 -> Like-sign (pi-K-)

        //If pair is Unlike-sign
        if(charge==-1){

            //Filling of the D0 histograms
            massPt->Fill(dMass,d0Pt,reweight*reweight_eff);
            D0etaunlike->Fill(kp->eta());
      
            //Saving of the daughter tracks
            DaughterPionTrackVector.push_back(pion->id());
            DaughterKaonTrackVector.push_back(kaon->id());

            //Check ff the mass is in the D0 mass window
            //if(dMass>1.81&&dMass<1.91){

                //The event is noted
                IsThereD0 = true;

                //Saving of D0 four-momenta // E,       px,       py,        pz,         D0_antiD0,            D0 mass;
                FourMomentum D0_actual = {kp->Energy(), kp->Px(), kp->Py(),  kp->Pz(), (double)pion->charge(), dMass};
                D0_fourmomentum.push_back(D0_actual);
            //}

            //Loading of the rescaled RunID
            int runIndex = mPrescales->runIndex(mPicoD0Event->runId());

            //Filling of the histograms
            pioneta->Fill(pion->pMom().PseudoRapidity());
            kaoneta->Fill(kaon->pMom().PseudoRapidity());
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

        } //end of the Unlike-sign if

        //If pair is Like-sign
        if(charge>0){

            //Filling of the D0 histograms
            massPtLike->Fill(dMass,d0Pt,reweight*reweight_eff);
            D0etalike->Fill(kp->eta());

        } //end of the Like-sign if
       
    } //end of the pair cuts if

  } //end of the D0 candidate loop

//---------------------------------------------------------------------
//-------------------JET-RECONSTRUCTION-PART---------------------------
//---------------------------------------------------------------------

  //Table constants
  double pimass = 0.13957018;

  //Check if there is a D0 candidate in the event
  if(IsThereD0){

      //NEventsCuts: Good D0 candidate
      NEventsCuts->Fill(7);

      //Initialisation of the input particle vectors for FastJet
      vector<fastjet::PseudoJet> input_particles;
      vector<fastjet::PseudoJet> chargedjetTracks;
      vector<fastjet::PseudoJet> neutraljetTracks;

      //Radius of the jet
      double R = 0.4;

      //Loop over all D0 candidates in the event.
      //If there are more than one, the jet reconstruction is done for each D0 candidate separately
      //ignoring the other not reconstructed D0 candidates in the event.
      for (unsigned int nD0 = 0; nD0 < D0_fourmomentum.size(); nD0++) {

//-----------D0-track--------------------------------------------------------

        //Defining the four-momentum of the D0 candidate
        fastjet::PseudoJet pj(D0_fourmomentum[nD0].px,D0_fourmomentum[nD0].py,D0_fourmomentum[nD0].pz,D0_fourmomentum[nD0].E);

        //Set flag to 2 if the D0 candidate is a D0 and to -2 if it is a anti-D0
        //It cannot be -1 or 1 because the default flag in FastJet is -1.
        pj.set_user_index(D0_fourmomentum[nD0].D0_antiD0*2);

        //Add the D0 candidate to the inclusive particle vector
        input_particles.push_back(pj);

//-----------Neutral-tracks--------------------------------------------------

        //Fill array Sump with momenta of tracks which are matched to BEMC towers
        GetCaloTrackMomentum(picoDst,TVector3(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z()));

        //Loop over all tracks in the event
        for (int iTow = 0; iTow < 4800; iTow++){

           //Get the tower hit
           StPicoBTowHit *towHit = picoDst->btowHit(iTow);

           //Check if the tower is bad or missing information
           if (!towHit || towHit->isBad()) continue;

           //Get the alternative counting method for the tower ID
           //NumericIndex - Counting from 0 to 4799
           //SoftId - Counting from 1 to 4800
           int realtowID = towHit->numericIndex2SoftId(iTow);

           //Initialize the tower energy
           double towE;

           //Calculation of the tower energy depending on the year
           //In 2014, there was a problem with the energy calibration, so the energy has to be corrected
           if(mYear==2014) {
              //Exclude bad towers, saved in JetInfo.h
              if (BadTowerMap[realtowID]) continue;
                  //Calculate the tower energy
                  towE = GetTowerCalibEnergy(iTow+1);
              }
              if(mYear==2016) {
                  //Exclude bad towers, saved in Calibration2016.h
                  if (EnergyBadTowerMap[realtowID]) continue;
                  //Get the tower energy
                  towE = towHit->energy();
              }

              //Subtraction of the hadronic background, fHadronCorr is set in RunPicoD0AnaMaker.C (0 - no or 1 - yes)
              towE-= fHadronCorr*Sump[iTow];

           //If the tower energy is negative, set it to 0
           if (towE < 0) towE = 0;

           //Initialize the tower geometry
           StEmcGeom* mEmcGeom;
           mEmcGeom = StEmcGeom::getEmcGeom("bemc");
           float Toweta_tmp = 0, Towphi = 0;

           //Get the eta and phi of the tower
           mEmcGeom->getEtaPhi(realtowID,Toweta_tmp,Towphi);

           //Correct the eta of the tower for the vertex position
           //Because the loaded eta is w.r.t. the center of the TPC, but the vertex do not have to be in the center
           float Toweta = vertexCorrectedEta(Toweta_tmp, event->primaryVertex().z());

           //Calculate the transverse energy
           //max eta 1.05258 max difference: ET = 0.124452 for E = 0.2, if we cut on |Vz| < 30 cm
           double ET = towE/cosh(Toweta);

           //If the transverse energy is greater than 30 GeV, discard the event
           if (ET > 30) {
                TowArr.clear();
                TowEta.clear();
                TowPhi.clear();
                Clusters.clear();
                return kStOK;
           }

           //Initialize and calculate the momentum components
           double px,py,pz;
           px = ET*cos(Towphi);
           py = ET*sin(Towphi);
           pz = towE*tanh(Toweta);

           //Create a jet with the calculated momentum components
           PseudoJet inputTower(px, py, pz, towE);

           //Discarding of the towers with low transverse energy (Defined in RunPicoD0AnaMaker.C as setCutETmin)
           //TrackBasedJets = charged tracks + D0 (Defined in RunPicoD0AnaMaker.C as setOnlyTrackBasedJets)
           if (inputTower.perp() > fETmincut && OnlyTrackBasedJets == 0){
               //Set the flag to 10 if the particle is neutral
               inputTower.set_user_index(10);
               //Add the neutral particle to the neutral particle vector
               neutraljetTracks.push_back(inputTower);
               //Add the neutral particle to the inclusive particle vector
               input_particles.push_back(inputTower);
           } //End of minimum ET cut

        } //End of loop over all towers

        //NEventsCuts: Cal.: E_{T} < 30 GeV
        //Only for the first D0 otherwise it is counted multiple times
        if (nD0==0) NEventsCuts->Fill(8);

//-----------Charged-tracks--------------------------------------------------

        //Loop over all tracks in the event
        for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack){

            //The i-th track is loaded
            StPicoTrack* trk = picoDst->track(iTrack);

            //Check if the track exists
            if(!trk) continue;

            //Filling jet track histogram before any cuts
            JetTracksdEdx->Fill(trk->gPtot()*trk->charge(),trk->dEdx());

            //Check if the track is a good track
            if (!isGoodJetTrack(trk,event)) continue;

            //Loading of the pT
            double pT = trk->gMom().Perp();
            //Check if the pT is above 0.2 GeV/c or if it is NaN, because NaN!=NaN
            if(pT != pT) continue; // NaN test.
            //Loading of the eta, phi, dca and charge
            float eta = trk->gMom().PseudoRapidity();
            float phi = trk->gMom().Phi();
            float dca = (TVector3(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z()) - trk->origin()).Mag();
            float charged = trk->charge();

            //Filling jet track histogram before after cuts
            JetTracksdEdxCut->Fill(trk->gPtot()*trk->charge(),trk->dEdx());

            //If the track is not daughter pion nor kion then...
            if (DaughterPionTrackVector[nD0] != trk->id() && DaughterKaonTrackVector[nD0] != trk->id()){

                //Defining the four-momentum of the charged particle, assumed pi+ mass
                //                       px,       py,               pz,                                 E = sqrt(p^2 + m^2),
                fastjet::PseudoJet pj(trk->gMom().x(),trk->gMom().y(),trk->gMom().z(), sqrt(trk->gMom().Mag()*trk->gMom().Mag()+pimass*pimass));

                //Add the charged particle to the charged particle vector
                chargedjetTracks.push_back(pj);
                //Add the charged particle to the inclusive particle vector
                input_particles.push_back(pj);

            }//End of if the track is not daughter pion nor kion

            //Filling the track histograms
            hpT_tr->Fill(pT, reweight);
            heta_phi_tr->Fill(phi + TMath::Pi(), eta,  reweight);
            heta_tr->Fill(eta, reweight);
            hphi_tr->Fill(phi + TMath::Pi(), reweight); //to shift by pi
            hdca_z_tr->Fill(dca, event->primaryVertex().z(), reweight);
            hdca_pT->Fill(dca, pT, reweight);
            hdca_tr->Fill(dca, reweight);
            hcharged_tr->Fill(charged, reweight);

        } //End of loop over all tracks

//-----------Background-estimation--------------------------------------------------

        //Definition of jets for background estimation
        //Contrary to the inclusive jets, the background jets are reconstructed with the kt algorithm (recommended choice)
        JetDefinition jet_def_bkgd(kt_algorithm, R);

        //Definition of the area for background estimation626
        AreaDefinition area_def_bkgd(active_area_explicit_ghosts,GhostedAreaSpec(fGhostMaxrap, 1, 0.01));

        //Remove two hardest jets in central collisions, one in others
        if (centrality == 7 || centrality == 8) nJetsRemove = 2;

        //Definition of the selector for background estimation (eta and pt cut + remove the n hardest jets)
        Selector selector = SelectorAbsEtaMax(1.0) * (!SelectorNHardest(nJetsRemove)) * SelectorPtMin(0.01);

        //Definition of the background estimator
        JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);

        //Estimation of the background using only charged tracks
        bkgd_estimator.set_particles(chargedjetTracks);

        //Calculation of the rho (median) and sigma (fluctuations of the median) for the background
        float rho = bkgd_estimator.rho();
        //float rho_sigma = bkgd_estimator.sigma();

//-----------Jet-reconstruction-and-variable-calculations----------------------------

        //Inclusive (or track based) jet definition
        fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);

        //Definition of the area for jet reconstruction
        AreaDefinition area_def(active_area_explicit_ghosts, GhostedAreaSpec(fGhostMaxrap, 1, 0.01));

        //Definition of the clustering
        ClusterSequenceArea clust_seq_hard(input_particles, jet_def, area_def);

        //Jet minimum pT cut
        double ptmin = 0.0;

        //Sorting of the jets by pT
        vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq_hard.inclusive_jets(ptmin));

        //Print the results
        /*
        cout << "Ran " << jet_def.description() << endl;
        printf("%5s %15s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt", "index");
        */

        //Loop over all jets
        for (unsigned int i = 0; i < inclusive_jets.size(); i++) {

            //Exclude jets with |eta| > 1 - R
            if (abs(inclusive_jets[i].pseudorapidity()) > (1.0 - R)) continue;

            //Initialize the variables (angularities, z-value)
            int user_index = 0;
            double Delta_R_D0 = 0;
            double lambda_alpha_1 = 1.;
            double lambda_alpha_1half = 1.5;
            double lambda_alpha_2 = 2.;
            double lambda_alpha_3 = 3.;
            double lambda_kappa = 1.;
            double zet = 0;
            double lambda_1_1 = 0;
            double lambda_1_1half = 0;
            double lambda_1_2 = 0;
            double lambda_1_3 = 0;
            double neutralpT = 0;

            //Calculate the jet pT + background subtraction (= _corr)
            //pT(sub) = pT - rho * A_jet
            //Since z-value requires px and py, it is calculated as well
            double pT_jet = inclusive_jets[i].perp();
            double pT_jet_corr = pT_jet - rho * inclusive_jets[i].area();
            double px_jet = inclusive_jets[i].px();
            double px_jet_corr = px_jet - rho * inclusive_jets[i].area_4vector().px();
            double py_jet = inclusive_jets[i].py();
            double py_jet_corr = py_jet - rho * inclusive_jets[i].area_4vector().py();

            //Loading the constituents of the jet
            const vector<fastjet::PseudoJet>& constituents = inclusive_jets[i].constituents();

            //Loop over all constituents of the i-th jet
            for (vector<fastjet::PseudoJet>::const_iterator particle = constituents.begin(); particle != constituents.end(); ++particle) {

                //Fraction of neutral particles
                if (particle->user_index() == 10) neutralpT += particle->pt();

                //Loading of the particle index
                int index = particle->user_index();

                //Calculating the delta R = sqrt(delta eta^2 + delta phi^2)
                double Delta_R =delta_R(inclusive_jets[i].eta(),inclusive_jets[i].phi(),particle->eta(),particle->phi());

                //Angularities are calculated only for track based particles (charged + D0)
                if (particle->user_index() != 10) {
                    lambda_1_1+=pow(particle->pt()/pT_jet_corr,lambda_kappa)*pow( Delta_R /jet_def.R() ,lambda_alpha_1);
                    lambda_1_1half+=pow(particle->pt()/pT_jet_corr,lambda_kappa)*pow( Delta_R /jet_def.R() ,lambda_alpha_1half);
                    lambda_1_2+=pow(particle->pt()/pT_jet_corr,lambda_kappa)*pow( Delta_R /jet_def.R() ,lambda_alpha_2);
                    lambda_1_3+=pow(particle->pt()/pT_jet_corr,lambda_kappa)*pow( Delta_R /jet_def.R() ,lambda_alpha_3);
                }

                //Check if the constituent is D0 (D0 = 2, antiD0 = -2)
                if (abs(index) == 2 ) {

                    //user_index is used as a flag to check if the jet contains D0
                    user_index = index;

                    //Delta R for D0
                    Delta_R_D0 = Delta_R;

                    //z-value, z = pT(D0)*^pT(jet)/|pT(jet)|
                    zet = (D0_fourmomentum[nD0].px*px_jet_corr+D0_fourmomentum[nD0].py*py_jet_corr)/(pT_jet_corr*pT_jet_corr);
                }

            } //End of loop over all constituents of the i-th jet

            //Calculation of the fraction of neutral particles
            double nfraction = neutralpT/pT_jet;

            //If the fraction is too high, jet is rejected (Parameter is saved in RunPicoD0AnaMaker.C as setMaxNeutralFraction)
            if (nfraction > maxneutralfrac) continue;

            //If the jet contains D0
            if (abs(user_index) ==2){

                //Calculation of the D0 mass
                double D0mass = D0_fourmomentum[nD0].D0Mass;
                //Calculation of the D0 pT (pT=sqrt(px^2+py^2))
                double D0_pT = sqrt(D0_fourmomentum[nD0].px*D0_fourmomentum[nD0].px+D0_fourmomentum[nD0].py*D0_fourmomentum[nD0].py);

                //Fill the TNtuple
                Jets->Fill( RunId*user_index/2.,                      //EventID, positive for D0, negative for antiD0
                            centrality,                               //Centrality
                            reweight,                                 //Centrality reweighting factor
                            D0_fourmomentum.size(),                   //Number of D0 in the event
                            inclusive_jets[i].pseudorapidity(),       //Jet eta
                            pT_jet,                                   //Jet pT
                            pT_jet_corr,                              //Jet pT after background subtraction
                            D0mass,                                   //D0 mass
                            Delta_R_D0,                               //Delta R between D0 and jet axis
                            D0_pT,                                    //D0 pT
                            lambda_1_1,                               //Angularity lambda_1_1
                            lambda_1_1half,                           //Angularity lambda_1_1.5
                            lambda_1_2,                               //Angularity lambda_1_2
                            lambda_1_3,                               //Angularity lambda_1_3
                            zet                                       //z-value
                            );

                Jets2->Fill(    centrality,                           //Centrality
                                reweight,                             //Centrality reweighting factor
                                inclusive_jets[i].phi(),              //Jet phi
                                picoDst->event()->grefMult(),         //grefMult
                                rho,                                  //density of the background
                                inclusive_jets[i].area(),             //area of the jet
                                inclusive_jets[i].rap()               //rapidity of the jet
                                );


                //Fill the histogram (pt vs background density)
                Jet_grefmult_pt_background->Fill(picoDst->event()->grefMult(),rho);
                //Rapidity calculations and filling histogram
                double D0_rapidity = 1./2.*log((D0_fourmomentum[nD0].E+D0_fourmomentum[nD0].pz)/(D0_fourmomentum[nD0].E-D0_fourmomentum[nD0].pz));
                Jet_D0pT_vs_D0rapidity->Fill(D0_rapidity,D0_pT);
                Jet_D0pT_vs_Jetrapidity->Fill(inclusive_jets[i].rapidity(),D0_pT);
                Jet_phi->Fill(inclusive_jets[i].phi());


                //Print colorfully the D0-jet information
                /*printf("\033[32m%5u %15.8f %15.8f %15.8f %15d\033[0m\n", i, inclusive_jets[i].rap(), inclusive_jets[i].phi(), inclusive_jets[i].perp(), user_index);*/

            } else{

                //Print the jet information
                /*printf("%5u %15.8f %15.8f %15.8f %15d\n", i, inclusive_jets[i].rap(), inclusive_jets[i].phi(), inclusive_jets[i].perp(), user_index);*/

            } //End of if the jet contains D0

        } //End of loop over all jets

        //Delete vectors to avoid a pileup
        inclusive_jets.clear();
        input_particles.clear();
        chargedjetTracks.clear();
        neutraljetTracks.clear();

      } //End of D0 loop

  } //End of IsthereD0 condition

  //Delete daughter particles vectors
  DaughterPionTrackVector.clear();
  DaughterKaonTrackVector.clear();

  //End of the event
  return kStOK;
}

//---------------------------------------------------------------------
//-----------------------------FUNCTIONS-------------------------------
//---------------------------------------------------------------------

Double_t StPicoD0AnaMaker::vertexCorrectedEta(double eta, double vz) {
    //Function to correct the eta value of a track for the z-position of the primary vertex

    //eta = -log(tan(theta/2)) => theta = 2*atan(exp(-eta))
    double tower_theta = 2.0 * atan(exp(-eta));

    //If eta = 0 then z = 0
    //Else calculate z position
    double z = 0.0;
    if (eta != 0.0) z = mBarrelRadius / tan(tower_theta);

    //Difference between the z position of the track and the z position of the primary vertex
    double z_diff = z - vz;

    //Calculate the corrected theta value
    double theta_corr = atan2(mBarrelRadius, z_diff);

    //Calculate the corrected eta value
    double eta_corr = -log(tan(theta_corr / 2.0));

    return eta_corr;

    //Function returns the corrected eta value
}
//---------------------------------------------------------------------------
Bool_t StPicoD0AnaMaker::GetCaloTrackMomentum(StPicoDst *mPicoDst, TVector3 mPrimVtx) {
  //Function to calculate the momentum of the track matched to the calorimeter tower

  //Loading of the number of tracks in the event
  UInt_t nTracks = mPicoDst->numberOfTracks();

  //Loop over all tracks in the event
  for (unsigned int itrack = 0; itrack < nTracks; itrack++) {

      //Loading the track
      StPicoTrack *trk = mPicoDst->track(itrack);
      //Loading of the global momentum
      TVector3 gMom = trk->gMom();

      //Loading of the pT
      double pT = gMom.Perp();
      //Check if the pT is above 0.2 GeV/c or if it is NaN, because NaN!=NaN
      if(pT != pT || pT < 0.2) continue;
      //Loading of eta
      float eta = gMom.PseudoRapidity();
      //Exclude tracks outside of the TPC acceptance
      if (fabs(eta) > 1) continue;
      //Loading of phi
      //float phi = gMom.Phi();

      //Loading of the number of hits
      float nHitsFit = trk->nHitsFit();
      //Loading of the number of hits possible
      float nHitsMax = trk->nHitsMax();
      //Exclude tracks with less than 15 hits or with a ratio of less than 0.52
      if (nHitsFit < 15 || nHitsFit/nHitsMax < 0.52) continue;

      //Loading of the value of the magnetic field
      double Bfield = mPicoDst->event()->bField();
      //Loading of the helix
      StPicoPhysicalHelix trkhelix = trk->helix(Bfield);

      //Loading of the primary vertex
      float vtx_x = mPrimVtx.x();
      float vtx_y = mPrimVtx.y();
      float vtx_z = mPrimVtx.z();

      //Calculation of the DCA to the primary vertex
      TVector3 dcaPoint = trkhelix.at(trkhelix.pathLength(vtx_x, vtx_y));
      //Calculation of the DCA in the x-y plane
      float dca_z = dcaPoint.z() - vtx_z;
      //Exclude tracks with a DCA to the primary vertex in z of more than maxdcazhadroncorr (in RunPicoD0AnaMaker.C)
      if (fabs(dca_z) > maxdcazhadroncorr) continue;

      //Initialization and loading of the tower index
      int TowIndex = -99999;
      TowIndex = trk->bemcTowerIndex();

      //Initialization of the momentum
      float p = 0;

      //Check if the track is matched to a tower
      if (TowIndex > 0) {

        //Loading of the momentum
        p = gMom.Mag();
        //Summing up the momentum of all tracks matched to the same tower
        Sump[TowIndex-1] += p;

      }
  
  } //End of track loop

  return true;

  //Function returns true if it was successful and Sump filled with the momentum of all tracks matched to a tower
}
//---------------------------------------------------------------------------
Double_t StPicoD0AnaMaker::GetTowerCalibEnergy(Int_t TowerId){
  //Function calculates the calibrated energy of a tower

  //Loading of the tower
  StPicoBTowHit *tower = static_cast<StPicoBTowHit*>(picoDst->btowHit(TowerId-1));

  //Initialization of the pedestal, rms and status
  Float_t pedestal, rms;
  Int_t status;

  //Loading of the pedestal, rms and status (it does not work, if you use root instead of root4star)
  mTables->getPedestal(BTOW, TowerId, 0, pedestal, rms);
  mTables->getStatus(BTOW, TowerId, status);

  //Initialization of the tower coefficients
  Double_t *TowerCoeff;

  //Tower coefficients for the different runs, parameters are saved in BemcNewCalib.h
  if(picoDst->event()->runId() <= 15094020) TowerCoeff = CPre;
  else TowerCoeff = CLowMidHigh;

  //Calculation of the calibrated energy E=C*(ADC-Pedestal)
  Double_t calibEnergy = TowerCoeff[TowerId-1]*(tower->adc() - pedestal);

  return calibEnergy;

  //Function returns the calibrated energy of the tower
}
//---------------------------------------------------------------------------
bool StPicoD0AnaMaker::IsBadEnergyRun(int runID) {
    // Check if the run is in the list of BEMC bad runs

    for (unsigned int i = 0; i < sizeof(EnergyBadRunList)/sizeof(EnergyBadRunList[0]); i++) {
        if (EnergyBadRunList[i] == runID) {
            return true;
        }
    }
    return false;

    //Function returns true if the run is in the list of bad runs
}
//---------------------------------------------------------------------------
int StPicoD0AnaMaker::isD0PairCentrality_pt(StKaonPion const* const kp, int Centrality, int mYear) const{
    //Check if the pair passes the cuts for D0

    //Loading the daughter particles tracks
    StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
    StPicoTrack const* pion = picoDst->track(kp->pionIdx());

    //If pseudorapidity of either daughter is greater than 1, it returns 0
    if(fabs(kaon->gMom().PseudoRapidity())>1 || fabs(pion->gMom().PseudoRapidity())>1)  return 0;

    //Initialisation of the pairCuts boolean
    bool pairCuts = false;

    //Recalculation of the centrality binning
    int Centrality2 =  (Centrality == 8 || Centrality == 7) ? 0 :
                       (Centrality == 6) ? 1 :
                       (Centrality == 5 || Centrality == 4) ? 2 :
                       (Centrality == 3 || Centrality == 2) ? 3 :
                       (Centrality == 1 || Centrality == 0) ? 4 : -1;
    // Centr.   0-10%       10-20%         20-40%         40-60%        60-80%
    // C_ID     8,7        6               5,4             3,2         1,0
    // bin       0         1                 2               3           4


    //Recalculation of the momentum binning
    int KPMom = (kp->pt() < 0.5) ? 0 : (kp->pt() < 1) ? 1 : (kp->pt() < 2) ? 2 : (kp->pt() < 3) ? 3 : (kp->pt() < 5) ? 4 : 5;

    //Check if the pair passes the particular cuts
    //Parameters are saved in StCuts.cxx
    if (mYear == 2014){
        pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < mycuts::DCA_D0_cut_2014[KPMom][Centrality2] &&
                    kp->pionDca() > mycuts::pionDCA_cut_2014[KPMom][Centrality2] && kp->kaonDca() > mycuts::kaonDCA_cut_2014[KPMom][Centrality2] &&
                    kp->dcaDaughters() < mycuts::pionkaonDCA_cut_2014[KPMom][Centrality2] && kp->decayLength()> mycuts::D0_decayLength_cut_2014[KPMom][Centrality2] &&
                    cos(kp->pointingAngle()) > mycuts::cosTheta_2014;
    } else if (mYear == 2016){
        pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < mycuts::DCA_D0_cut_2016[KPMom][Centrality2] &&
                    kp->pionDca() > mycuts::pionDCA_cut_2016[KPMom][Centrality2] && kp->kaonDca() > mycuts::kaonDCA_cut_2016[KPMom][Centrality2] &&
                    kp->dcaDaughters() < mycuts::pionkaonDCA_cut_2016[KPMom][Centrality2] && kp->decayLength()> mycuts::D0_decayLength_cut_2016[KPMom][Centrality2] &&
                    cos(kp->pointingAngle()) > mycuts::cosTheta_2016;
    }

    //Calculation of the product of the daughter charges
    int charge = kaon->charge() * pion->charge();

    //If the daughter charges are the same, it changes charge to 1 (K+) or 2 (K-)
    if(charge>0) charge = kaon->charge()>0 ? 1:2;

    //If the pair passes the cuts, it returns the charge, else 0.
    if(pairCuts) return charge;
    else return 0;

    //Function returns:
    //0  - the pair does not pass the cuts
    //-1 - the unlike-sign pair passes the cuts
    //1  - the like-sign pair passes the cuts (pi+K+)
    //2  - the like-sign pair passes the cuts (pi-K-)
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodEvent(int mYear,  TH1F* NEventsCuts){
  //Check if the event passes the cuts (MB triggers, vr, vz, vzVpdVz)

  //Loading event
  StPicoEvent *event = (StPicoEvent *)picoDst->event();

  //Checking triggers
  if (!isMBTrigger(mYear)) return false;
  //NEventsCuts: Triggers
  NEventsCuts->Fill(1);
  //Checking vr = sqrt(vx^2+vy^2)
  if (!(sqrt(event->primaryVertex().x()*event->primaryVertex().x()+event->primaryVertex().y()*event->primaryVertex().y()) < mycuts::vr)) return false;
  //NEventsCuts: v_r
  NEventsCuts->Fill(2);
  //Checking vz
  if (!(fabs(event->primaryVertex().z()) < mycuts::vz)) return false;
  //NEventsCuts: v_z
  NEventsCuts->Fill(3);
  //Checking vzVpdVz
  if (!(fabs(event->primaryVertex().z() - event->vzVpd()) < mycuts::vzVpdVz)) return false;
  //NEventsCuts: |v_z - v_z_vpd|
  NEventsCuts->Fill(4);

  return true;

  //Function returns true if the event passes the cuts
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isMBTrigger(int mYear){
 //Function checks if the event is triggered by the MB trigger

 //Initialization of the set of triggers
 const std::set<int>* mbTriggers = nullptr;

 //Different triggers for different years, saved in StCuts.cxx
 if(mYear ==2016) mbTriggers = &mycuts::mbTriggers2016;
 if(mYear ==2014) mbTriggers = &mycuts::mbTriggers2014;

 //Loading event and checking if it is triggered by the MB trigger
 StPicoEvent* event = static_cast<StPicoEvent*>(mPicoDstMaker->picoDst()->event());
 return std::any_of(mbTriggers->begin(), mbTriggers->end(), [&](int trigger) { return event->isTrigger(trigger); });

 //Function returns true if the event is triggered by the MB trigger
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodTrack(StPicoTrack const * const trk) const{
    // Require at least one hit on every layer of PXL and IST.
    // It is done here for tests on the preview II data.
    // The new StPicoTrack which is used in official production has a method to check this

    //Check if the track meets the HFT requirement
    //2014 - Require at least one hit on every layer of PXL and IST
    //2016 - Require at least one hit on every layer of PXL and (IST or SST)
    //Both can be written in te same way
    bool HFTCondition = (trk->hasPxl1Hit() && trk->hasPxl2Hit()) && (trk->hasSstHit() || trk->hasIstHit());

    //In StCuts.cxx is defined if the HFT is required and the nHitsFit and minPt values.
    return trk->gPt() > mycuts::minPt && trk->nHitsFit() >= mycuts::nHitsFit && (HFTCondition || !mycuts::requireHFT);

    //Function returns true if track is good
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodJetTrack(StPicoTrack const * const trk,StPicoEvent const *const myEvent) const{
    //Check if the track meets the jet track cuts
    //Parameters saved in StCuts.cxx

    //pT range cut
    bool pTTrackJetCut = trk->gPt() > mycuts::jetTrackPtMin && trk->gPt() < mycuts::jetTrackPtMax;
    //eta range cut
    bool etaTrackJetCut = fabs(trk->gMom().PseudoRapidity()) < mycuts::jetTrackEta;
    //nHitsFit cut
    bool nHitsTrackJetCut = trk->nHitsFit() >= mycuts::jetTracknHitsFit;
    //nHitsRatio cut
    bool nHitsRatioTrackJetCut = (1.0*trk->nHitsFit()/trk->nHitsMax())>=mycuts::jetTracknHitsRatio;
    //DCA cut
    bool dcaTrackJetCut = fabs(trk->gDCA(myEvent->primaryVertex().x(),myEvent->primaryVertex().y(),myEvent->primaryVertex().z())) < mycuts::jetTrackDCA;

    return pTTrackJetCut && etaTrackJetCut && nHitsTrackJetCut && nHitsRatioTrackJetCut && dcaTrackJetCut;

    //Return true if all the cuts are passed
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTpcPion(StPicoTrack const * const trk) const{
    //Check if the track meets nsigma (TPC) requirement for pion
    //Parameters are saved in StCuts.cxx

    return fabs(trk->nSigmaPion()) < mycuts::nSigmaPion;

    //Function returns true if track is good
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTpcKaon(StPicoTrack const * const trk, StThreeVectorF const* const pVtx) const{
    //Check if the track meets nsigma (TPC) requirement for kaon
    //Parameters are saved in StCuts.cxx

    return fabs(trk->nSigmaKaon()) < mycuts::nSigmaKaon;

    //Function returns true if track is good
}
//-----------------------------------------------------------------------------
float StPicoD0AnaMaker::getTofBeta(StPicoTrack const * const trk, StThreeVectorF const* const pVtx,StPicoDst const* const picoDst) const{
    //Calculation of beta for the track

    //index2tof is index of the track in the StPicoBTofPidTraits array
    int index2tof = trk->bTofPidTraitsIndex();

    //Initialization of beta
    float beta = std::numeric_limits<float>::quiet_NaN();

    //If index2tof is positive, than the track has a match in the TOF
    if(index2tof >= 0){

        //Getting the pointer to the StPicoBTofPidTraits object
        StPicoBTofPidTraits *tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);

        //If pointer is not null, than we can calculate beta
        if(tofPid){

            //Calculation of beta
            beta = tofPid->btofBeta();

            //If for some reason beta is negative, than we can try to calculate beta using the pathlength and tof
            if (beta < 1e-4){

                //Getting the hit position in the TOF
                StThreeVectorF const btofHitPos = StThreeVectorF(tofPid->btofHitPos().x(),tofPid->btofHitPos().y(),tofPid->btofHitPos().z());

                //Getting the helix of the track
                StPicoPhysicalHelix helix = trk->helix(picoDst->event()->bField());

                //Calculation of the pathlength
                float L = tofPathLength(pVtx, &btofHitPos, helix.curvature());

                //Calculation of the time of flight
                float tof = tofPid->btof();

                //Calculation of beta for positive values of tof
                if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
                //...else beta is not defined
                else beta = std::numeric_limits<float>::quiet_NaN();

            } //End of beta < 1e-4

        } //End of tofPid != NULL

    } //End of index2tof >= 0

    return beta;

    //Function returns beta of the track
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTofKaon(StPicoTrack const * const trk, float beta) const{
    //Check if the track meets |1/beta-1/beta_K| (TOF) requirement for kaon

    //Initialization of tofKaon
    bool tofKaon = false;

    //If beta is positive, than we can calculate |1/beta-1/beta_K|
    if(beta>0){

        //Calculation of the global total momentum
        double ptot = trk->gMom().Mag();

        //Calculation of the expected beta for kaons
        float beta_k = ptot/sqrt(ptot*ptot+M_KAON_PLUS*M_KAON_PLUS);

        //Check if the track meets the TOF requirement
        //Parameters are saved in StCuts.cxx
        tofKaon = fabs(1/beta - 1/beta_k) < mycuts::kTofBetaDiff ? true : false;

    } //End of beta > 0

    return tofKaon;

    //Function returns true if track is good based on TOF information
}








