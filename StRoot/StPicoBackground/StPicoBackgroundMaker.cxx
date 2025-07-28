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
#include "StPicoBackgroundMaker.h"
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
#include "TF1.h"
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
#include "TRandom.h"
#include "TRandom3.h"
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
#include <fastjet/PseudoJet.hh>  // Pro třídu PseudoJet
#include <fastjet/ClusterSequence.hh>  // Pro klasické metody clusteringu (např. k nejbližší sousedství)
#include <fastjet/contrib/ConstituentSubtractor.hh>  // Pro ConstituentSubtractor
#include <fastjet/Selector.hh>  // Pro práci s selektory
#include <fastjet/JetDefinition.hh>  // Pro definici jetů (třeba k-means, anti-kt atd.)
#include <fastjet/AreaDefinition.hh>  // Pro definici plochy (area) v analýzách
#include <fastjet/tools/JetMedianBackgroundEstimator.hh> 

#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/contrib/IterativeConstituentSubtractor.hh" 

//#include "/gpfs01/star/pwg/lomicond/Ondrej/Jets/Alex_install/install/fastjet-contrib-install/include/fastjet/contrib/SignalFreeBackgroundEstimator.hh"
using namespace std;
using namespace fastjet;
//-------------------------------------

ClassImp(StPicoBackgroundMaker)

StPicoBackgroundMaker::StPicoBackgroundMaker(char const * name, /*char const * inputFilesList,*/ char const * outName, StPicoDstMaker* picoDstMaker, StRefMultCorr* grefmultCorrUtil, int pYear, bool Sim, char const *filename)
        : StMaker(name),
          mOutFileName(outName), // Initialize in the same order as declared
        //  mChain(NULL),
          mOutputFile(NULL),
          mGRefMultCorrUtil(grefmultCorrUtil),
          mPicoDstMaker(picoDstMaker),
          mPicoD0Event(NULL),
         /* mInputFileList(inputFilesList),*/
          mEventCounter(0),
          mYear(pYear),
          picoDst(NULL),
          mPrescales(NULL),
          mHFCuts(NULL),
          mSimulation(Sim) {
 
           // fMCFileListName = filename;
          }

Int_t StPicoBackgroundMaker::Init()
{

 

  //mPicoD0Event = new StPicoD0Event();

  //mChain = new TChain("T");
 
 /*
  std::ifstream listOfFiles(mInputFileList.Data());
  if (listOfFiles.is_open())
  {
    std::string file;
    while (getline(listOfFiles, file))
    {
      LOG_INFO << "StPicoBackgroundMaker - Adding :" << file <<endm;
      mChain->Add(file.c_str());
      LOG_INFO<<" Entries = "<<mChain->GetEntries()<< endm; 
    }
  }
  else
  {
    LOG_ERROR << "StPicoBackgroundMaker - Could not open list of files. ABORT!" << endm;
    return kStErr;
  }
  
*/
  //The prescales are used only for rescaling the RunID to the lower numbers
  mPrescales = new StPicoPrescales(mycuts::prescalesFilesDirectoryName);

  //Loading the information from the D0EventMaker
 // mChain->GetBranch("dEvent")->SetAutoDelete(kFALSE);
 // mChain->SetBranchAddress("dEvent", &mPicoD0Event);

  //Output file
  mOutputFile = new TFile(mOutFileName.Data(), "RECREATE");
  mOutputFile->cd();

  //Number of Runs
  int nRuns = mPrescales->numberOfRuns();

  //Events histograms:
  vtxz = new TH1F("vtxz",";PVtx.z() [cm]; Count",100,-10,10);
  vtxr = new TH2D("vtxr",";PVtx.x() [cm]; PVtx.y() [cm]",100,-3,3,100,-3,3);
  hcentr = new TH1D("hcentr",";C_{ID}",9,-0.5,8.5);
  hcentrW = new TH1D("hcentrW",";C_{ID}",9,-0.5,8.5);
  NEventsCuts = new TH1F("NEventsCuts", "NEventsCuts;Cuts;Count", 10, 0, 10);
      NEventsCuts->GetXaxis()->SetBinLabel(1, "All events");
      NEventsCuts->GetXaxis()->SetBinLabel(2, "Triggers");
      NEventsCuts->GetXaxis()->SetBinLabel(3, "V_{r}");
      NEventsCuts->GetXaxis()->SetBinLabel(4, "V_{z}");
      NEventsCuts->GetXaxis()->SetBinLabel(5, "|V_{z} - V_{z}^{VPD}|");
      NEventsCuts->GetXaxis()->SetBinLabel(6, "|V_{x,y,z}|!=0");
      NEventsCuts->GetXaxis()->SetBinLabel(7, "Good E run");
      NEventsCuts->GetXaxis()->SetBinLabel(8, "Centrality");
      NEventsCuts->GetXaxis()->SetBinLabel(9, "Good D^{0}");
      NEventsCuts->GetXaxis()->SetBinLabel(10, "Cal.: E_{T} < 30 GeV");
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
  Jets = new TNtuple("Jets", "Jets", "RunId:EventID:centrality:centr_weight:cone_eta:cone_phi:conePTSumCharged:conePxSumCharged:conePySumCharged:conePzSumCharged:conePTSumInclusive:conePxSumInclusive:conePySumInclusive:conePzSumInclusive:conePTSumNeutral:conePxSumNeutral:conePySumNeutral:conePzSumNeutral:rho:Q1_vec:Q2_vec:Psi_2:Q1_vec_rec:Q2_vec_rec:Psi_2_shifted");




VariableJets = {
    {"RunId", 0},        
    {"EventID", 1},    
    {"centrality", 2},   // 0 -> 70-80%, 1 -> 60-70%, 2 -> 50-60%, 3 -> 40-50%, 4 -> 30-40%, 5 -> 20-30%, 6 -> 10-20%, 7 -> 5-10%, 8 -> 0-5%
    {"centr_weight", 3},   // centrality weight       
    {"cone_eta", 4},       
    {"cone_phi", 5},      
    {"conePTSumCharged", 6},       
    {"conePxSumCharged", 7},     
    {"conePySumCharged", 8},      
    {"conePzSumCharged", 9},       
    {"conePTSumInclusive", 10},  
    {"conePxSumInclusive", 11},      
    {"conePySumInclusive", 12},       
    {"conePzSumInclusive", 13},   
    {"conePTSumNeutral", 14},  
    {"conePxSumNeutral", 15},      
    {"conePySumNeutral", 16},       
    {"conePzSumNeutral", 17},         
    {"bg_dens", 18},
    {"Q1_vec", 19},
    {"Q2_vec", 20},
    {"Psi_2", 21},
    {"Q1_vec_rec", 22},
    {"Q2_vec_rec", 23},
    {"Psi_2_shifted", 24}
};

////EventStats = new TNtuple("EventStats", "EventStats","RunId:centrality:centr_weight");

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
  ////StMaker* maker = GetMaker("Eread");
  mADCtoEMaker = dynamic_cast<StEmcADCtoEMaker*>(GetMaker("Eread"));
  assert(mADCtoEMaker);
  mTables = mADCtoEMaker->getBemcData()->getTables();
  
  ifstream filelistforMCEvents(fMCFileListName.Data());
  
  if(mSimulation){
  
	  string line;
	  while (getline(filelistforMCEvents, line))
	  {
	    TString s(line);
	    filenamesforHIOverlay.push_back(s);
	  }
  
  } else {
  
  
  }



  return kStOK;
}

//-----------------------------------------------------------------------------
StPicoBackgroundMaker::~StPicoBackgroundMaker(){
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
//--------Event-plane---------------
Int_t StPicoBackgroundMaker::EP_IsGoodTrack(StPicoTrack *trk, TVector3 pVertex){

	enum Position {West = -1, East = 1};
	Position Side = West;
	
	TVector3 mTrkMom;
	
	//Primary tracks
   	mTrkMom = trk->pMom(); //If not exists, it equals 0, 0, 0

	//0.2 < pt < 2 GeV
	float pt = mTrkMom.Perp();
	if (pt!=pt||pt <= 0.2 || pt >= 2.0) return 0;

	//0.05 < |eta| < 1.00 (West/East)
	float eta = mTrkMom.PseudoRapidity();
	if (abs(eta) <= 0.05 || abs(eta) >= 1) return 0;
	if (eta > 0) Side = East;

	//nHitsFit > 15
	float nHitsFit = trk->nHitsFit();
	float nHitsMax = trk->nHitsMax();
	
	//nHitsFit/nHitsMax => 0.52
	float nHitsRatio = 1.0*nHitsFit/nHitsMax;
	if (nHitsFit < 15 || nHitsRatio < 0.52) return 0;
	
	
	//DCA < 1 cm
	float dca = trk->gDCA(pVertex).Mag();
	if (dca >= 1) return 0;

return Side;
}


void StPicoBackgroundMaker::CalculateEventPlane(){


   StPicoEvent *eventQ = (StPicoEvent *)picoDst->event();

    //Primary vertex
    TVector3 prVertex = eventQ->primaryVertex();
    //Q-vectors
    double Q_1 = 0;
    double Q_2 = 0;
    
     for (UInt_t iTrack = 0; iTrack < picoDst->numberOfTracks(); iTrack++){
      
        StPicoTrack *trk = static_cast<StPicoTrack *>(picoDst->track(iTrack));
        if (!trk) continue;
        
        //{West = -1, East = 1}
        double Goodtrack = EP_IsGoodTrack(trk,prVertex);
        if (!abs(Goodtrack)) continue;
        
        TVector3 trackMom = trk->pMom();
        
        double pPt = trackMom.Perp();
        double phi = trackMom.Phi();
        
        //Q-vectors calculating
        Q_1 += 1.*pPt*cos(2*phi);
        Q_2 += 1.*pPt*sin(2*phi);
        
     }
    
   // cout << "Q_1: " << Q_1 << " Q_2: " << Q_2 << endl;
    
    fQ_1 = Q_1;
    fQ_2 = Q_2;
    
    //cout << Q_1 << " " << Q_2 << " Q1 a Q2" << endl;
    
    //Recentering
    //-------------------
    //You have to run the code for "all" events to get the mean value of Q vector
    double Q_1rc = -0.567;
    double Q_2rc = 1.574;
    
    
    double Q_1corr = Q_1 - Q_1rc;
    double Q_2corr = Q_2 - Q_2rc;
    
    fQ_1_rec = Q_1corr;
    fQ_2_rec = Q_2corr;    
    //-------------------
    
    double Psi_2 = 1./2 * TMath::ATan2(Q_2corr, Q_1corr);
    
    
    
    //-------------------
    //Psi shift	
    //You have to run the code for "all" events to get the mean value of Q vector (again)
    double DeltaPsi2 = 0;
    
std::vector<double> A_2 = {-0.0610273, -0.0120242, -0.00426917, 0.00316563, 0.00544594, 0.00601391, 0.0024604, 0.00164042, 0.0011545, -0.00280563, -0.0016235, 0.00195694, 0.00333251, 0.00114701, 9.20402e-05, 0.00145072, -0.00515215, -0.00140633, -0.00178873, -0.00189488, -0.00376466};
std::vector<double> B_2 = {0.0225684, -0.0170903, -0.00714623, -0.00231669, 0.00213819, -0.0082172, -0.00416794, 0.0054867, -0.0029049, -0.00133012, -0.00428737, -0.00368623, 0.00300743, 0.00132147, 0.000450483, -0.00332185, -0.00244103, 0.000236871, 0.00101393, 0.000561244, -0.00627737};



    double CorrectedPsi2 = Psi_2;
    for (int i = 1; i <= 21; i++){
         CorrectedPsi2 +=(1.0 / 2) * (2.0 / i) *(-A_2[i - 1] * cos(2 * i * Psi_2) + B_2[i - 1] * sin(2 * i * Psi_2));
    }
    
    if (Psi_2!=Psi_2 || CorrectedPsi2!=CorrectedPsi2) return;
    
	//Force the range (-pi/2,pi/2)
    CorrectedPsi2 = TMath::ATan2(TMath::Sin(2 * CorrectedPsi2), TMath::Cos(2 * CorrectedPsi2)) / 2.;
    Psi_2 = TMath::ATan2(TMath::Sin(2 * Psi_2), TMath::Cos(2 * Psi_2)) / 2.;
    
    //-------------------
    fPsi_2 = Psi_2;
    fPsi_2_shifted = CorrectedPsi2;
    
    //v2 calc
    
    
    
    

return;
}
//-----------------------------------------------------------------------------
Int_t StPicoBackgroundMaker::Finish(){
  //Write histograms and close the output file, if you create a new histogram it has to be added here
  LOG_INFO << " StPicoBackgroundMaker - writing data and closing output file " <<endm;
  fout.close();
  mOutputFile->cd();
/*
  //Events histograms:
  vtxz->Write();
  vtxr->Write();
  hcentr->Write();
  hcentrW->Write();
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
*/
  //2D D0 mass-pt like-sign and unlike-sign histograms
 // massPt->Write();
 // massPtLike->Write();

  //TNtuple D0-jets
  Jets->Write();
  ////EventStats->Write();
  
  //Closing of the output file
  mOutputFile->Close();

  //Closing of the input file including prescales
  delete mPrescales;

  return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoBackgroundMaker::Make()
{
  //Main function where the analysis is done, each "event" is analyzed here
 // readNextEvent();



  //Check if everything is loaded properly
  if (!mPicoDstMaker)
  {
    LOG_WARN << " StPicoBackgroundMaker - No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  picoDst = mPicoDstMaker->picoDst();
  if (!picoDst)
  {
    LOG_WARN << "StPicoBackgroundMaker - No PicoDst! Skip! " << endm;
    return kStWarn;
  }

  //if (picoDst->event()->eventId()!=2698540) return kStOK;

  //if (picoDst->event()->eventId()!=2804657)return kStOK;

  // -------------- USER ANALYSIS -------------------------

  //Loading of the raw D0 daughters
  TClonesArray const * aKaonPion;
 // if (!mSimulation) aKaonPion = mPicoD0Event->kaonPionArray();

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
  if (mYear==2016 && IsBadEnergyRun(picoDst->event()->runId())) return kStOK;

  //2014 Good run list
  bool isGoodRun = true;
  const std::set<int>* goodRunList = nullptr;
  if (mYear == 2014){
	int runIdLoad = picoDst->event()->runId();
	goodRunList = &mycuts::goodRun2014; 
	isGoodRun = goodRunList->find(runIdLoad) != goodRunList->end();
  }
	
  if (!isGoodRun) return kStOK;
  ///////	

    // Inicializace generátoru náhodných čísel
    TRandom3 randGen(0);

    // Generování náhodné polohy v eta-phi prostoru
    double eta_rnd = randGen.Uniform(-0.6, 0.6);
    double phi_rnd = randGen.Uniform(-TMath::Pi(), TMath::Pi());
//cout << phi_rnd << "=?" << endl;


  //NEventsCuts: Good E
  NEventsCuts->Fill(6);

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
  NEventsCuts->Fill(7);

  //Loading of the weight for the centrality correction
  double reweight = mGRefMultCorrUtil->getWeight();

  //Filling events histograms
  vtxz->Fill(pVtx.z());
  vtxr->Fill(pVtx.x(),pVtx.y());
  hcentr->Fill(centrality);
  hcentrW->Fill(centrality,reweight);
  
  ////EventStats->Fill(picoDst->event()->runId(),centrality,reweight); 

  //Loading event information
  int EventID = picoDst->event()->eventId();
  int RunId = picoDst->event()->runId();

  //Loading of the number of tracks in the event
  UInt_t nTracks = picoDst->numberOfTracks();


  
  
      double conePTSumCharged = 0;
      double conePxSumCharged = 0;
      double conePySumCharged = 0;
      double conePzSumCharged = 0;
      
      double conePTSumInclusive = 0;
      double conePxSumInclusive = 0;
      double conePySumInclusive = 0;
      double conePzSumInclusive = 0;    
      
      double conePTSumNeutral = 0;
      double conePxSumNeutral = 0;
      double conePySumNeutral = 0;
      double conePzSumNeutral = 0;   

CalculateEventPlane();
//----------------------------------------------------------------

//---------------------------------------------------------------------
//-------------------JET-RECONSTRUCTION-PART---------------------------
//---------------------------------------------------------------------

      //NEventsCuts: Good D0 candidate
      NEventsCuts->Fill(8);

      //Initialisation of the input particle vectors for FastJet
      vector<fastjet::PseudoJet> input_particles;
      vector<fastjet::PseudoJet> chargedjetTracks;
      vector<fastjet::PseudoJet> neutraljetTracks;

      //Radius of the jet
      double R = 0.4;



	//Delete all energies calculated for hadr. corr. from previous event
	SumE.fill(0);
	
//-----------Charged-tracks--------------------------------------------------

/*

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

           
       
           //Defining the four-momentum of the charged particle, assumed pi+ mass
           //                       px,       py,               pz,                                 E = sqrt(p^2 + m^2),
           fastjet::PseudoJet pj(trk->gMom().x(),trk->gMom().y(),trk->gMom().z(), sqrt(trk->gMom().Mag()*trk->gMom().Mag()+M_PION_PLUS*M_PION_PLUS));
           
           cout << "input_particles.push_back("<<trk->gMom().x()<<","<<trk->gMom().y()<<","<<trk->gMom().z()<<","<<sqrt(trk->gMom().Mag()*trk->gMom().Mag()+M_PION_PLUS*M_PION_PLUS)<<"); //charged" << endl;

	   //Set the flag to 3 if the particle is charged
	   pj.set_user_index(3);
           //Add the charged particle to the charged particle vector
           chargedjetTracks.push_back(pj);
           //Add the charged particle to the inclusive particle vector
           input_particles.push_back(pj);
	//    bool isInRndCone(double eta, double phi, double eta_rnd, double phi_rnd) const;
	if(isInRndCone(eta, phi, eta_rnd, phi_rnd)){
	
	conePTSumCharged+=pT;
	conePxSumCharged+=trk->gMom().X();
  	conePySumCharged+=trk->gMom().Y();
  	conePzSumCharged+=trk->gMom().Z();
  	conePTSumInclusive+=pT;
	conePxSumInclusive+=trk->gMom().X();
  	conePySumInclusive+=trk->gMom().Y();
  	conePzSumInclusive+=trk->gMom().Z();
  	
           }
        } //End of loop over all tracks
        
        
        */
//-----------Neutral-tracks--------------------------------------------------

/* //TEST
        //Fill array SumE with momenta of tracks which are matched to BEMC towers
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
                  if (BadTowerMap[realtowID-1]) continue;

                  //Calculate the tower energy
                  towE = GetTowerCalibEnergy(realtowID);           
                  ////towE = towHit->energy(); //Only test
              }
              if(mYear==2016) {
                  //Exclude bad towers, saved in Calibration2016.h
                  if (EnergyBadTowerMap[realtowID-1]) continue; //!!! Check -1
                  //Get the tower energy
                  towE = towHit->energy();
              }
              towE-= fHadronCorr*SumE[iTow];
              
           //If the tower energy is negative, set it to 0
           if (towE < 0) towE = 0;

           //Initialize the tower geometry
           //StEmcGeom* mEmcGeom;
           //mEmcGeom = StEmcGeom::getEmcGeom("bemc");
           StEmcPosition* mEmcPosition;
           mEmcPosition = new StEmcPosition();

           //Correct the eta of the tower for the vertex position
           //Because the loaded eta is w.r.t. the center of the TPC, but the vertex do not have to be in the center
           StThreeVectorF towerPosition = mEmcPosition->getPosFromVertex(StThreeVectorF(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z()), realtowID);
           float Toweta = towerPosition.pseudoRapidity();
           float Towphi = towerPosition.phi();

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

	  // Double_t p = 1.0*TMath::Sqrt(towE*towE - M_PION_PLUS*M_PION_PLUS);
	   Double_t p = 1.0*TMath::Sqrt(towE*towE - 0*0);
	   double posX = towerPosition.x();
  	   double posY = towerPosition.y();
  	   double posZ = towerPosition.z();
  	     Double_t r = TMath::Sqrt(posX*posX + posY*posY + posZ*posZ) ;
  	                double px,py,pz;
  	   px = p*posX/r*1.;
    	   py = p*posY/r*1.;
    	   pz = p*posZ/r*1.;

           //Create a jet with the calculated momentum components
           PseudoJet inputTower(px, py, pz, towE);
		//if (realtowID==4706) cout << "px: " << px << " py: " << py << " pz: " << pz << " ET: " << ET << " sumE: " << SumE[iTow] << endl;
           //Discarding of the towers with low transverse energy (Defined in RunPicoD0AnaMaker.C as setCutETmin)
           //TrackBasedJets = charged tracks + D0 (Defined in RunPicoD0AnaMaker.C as setOnlyTrackBasedJets)
           if (ET > fETmincut && OnlyTrackBasedJets == 0){
                   
               //Set the flag to 10 if the particle is neutral
               inputTower.set_user_index(10);
               //Add the neutral particle to the neutral particle vector
               neutraljetTracks.push_back(inputTower);
               //Add the neutral particle to the inclusive particle vector
               input_particles.push_back(inputTower);
               
         cout << "input_particles.push_back("<<px<<","<<py<<","<<pz<<","<<towE<<"); //neutral" << endl;
           if(isInRndCone(Toweta, Towphi, eta_rnd, phi_rnd)){
               	conePTSumInclusive+=TMath::Sqrt(px*px+py*py);
		conePxSumInclusive+=px;
  		conePySumInclusive+=py;
  		conePzSumInclusive+=pz;
  		conePTSumNeutral+=TMath::Sqrt(px*px+py*py);
		conePxSumNeutral+=px;
  		conePySumNeutral+=py;
  		conePzSumNeutral+=pz;
           }
               
           } //End of minimum ET cut

        } //End of loop over all towers

        //NEventsCuts: Cal.: E_{T} < 30 GeV
        //Only for the first D0 otherwise it is counted multiple times

        */




//-----------Background-estimation--------------------------------------------------


	/* //TEST
        //Definition of jets for background estimation
        //Contrary to the inclusive jets, the background jets are reconstructed with the kt algorithm (recommended choice)
        JetDefinition jet_def_bkgd(kt_algorithm, R, E_scheme, Best);

        //Definition of the area for background estimation
        AreaDefinition area_def_bkgd(active_area_explicit_ghosts,GhostedAreaSpec(1., 1, 0.01)); //Comp. test
       // AreaDefinition area_def_bkgd(active_area_explicit_ghosts,GhostedAreaSpec(1.2, 1, 0.005)); //Comp. test
        
        

        /////
        //unsigned int seed1 = 12345;
    	//unsigned int seed2 = 56789;
	//std::vector<int> seeds = { static_cast<int>(seed1), static_cast<int>(seed2) };
    	//fastjet::AreaDefinition area_def_bkgd = initial_area_def.with_fixed_seed(seeds);
	//////


        //Remove two hardest jets in central collisions, one in others
        if (centrality == 7 || centrality == 8) nJetsRemove = 2; //Comp. test
	//nJetsRemove = 2; //Comp. test

        //Definition of the selector for background estimation (eta and pt cut + remove the n hardest jets)
        //Selector selector = (!SelectorNHardest(nJetsRemove)) * SelectorAbsEtaMax(1.) * SelectorPtMin(0.01); //Comp. test
        Selector selector = (!SelectorNHardest(nJetsRemove)) * SelectorAbsEtaMax(0.6); //Comp. test

        //Definition of the background estimator
        JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);

        //Estimation of the background using only charged tracks
        bkgd_estimator.set_particles(input_particles); //Bug previously: chargedjetTracks instead of input_particles

        //Calculation of the rho (median) and sigma (fluctuations of the median) for the background
        float rho = bkgd_estimator.rho();
        //float emptyjets = bkgd_estimator.n_empty_jets(); //not used
       
       
       */ //TEST
   
//-----------Jet-reconstruction-and-variable-calculations----------------------------
/*
        //Inclusive (or track based) jet definition
        fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R, E_scheme, Best);

        //Definition of the area for jet reconstruction
	fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, GhostedAreaSpec(1., 1, 0.01)); //Comp. test
	//fastjet::AreaDefinition initial_area_def2(fastjet::active_area_explicit_ghosts, GhostedAreaSpec(1.2, 1, 0.005)); //Comp. test
	
	/////
        //unsigned int seed1 = 12345;
    	//unsigned int seed2 = 56789;
	//std::vector<int> seeds = { static_cast<int>(seed1), static_cast<int>(seed2) };
    	fastjet::AreaDefinition area_def = initial_area_def2.with_fixed_seed(seeds);
	//////
	
        //Definition of the clustering
        ClusterSequenceArea clust_seq_hard(input_particles, jet_def, area_def);

        //Jet minimum pT cut
        double ptmin = 0.0;

        //Sorting of the jets by pT
        vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq_hard.inclusive_jets(ptmin));
*/


                //Fill the TNtuple
		TupleVariables[VariableJets["RunId"]] = RunId;           		// RunID, positive for D0, negative for antiD0
		TupleVariables[VariableJets["EventID"]] = EventID; 
		TupleVariables[VariableJets["centrality"]] = centrality;                    		// Centrality
		TupleVariables[VariableJets["centr_weight"]] = reweight;  
		TupleVariables[VariableJets["cone_eta"]] = eta_rnd;
		TupleVariables[VariableJets["cone_phi"]] = phi_rnd; 
		TupleVariables[VariableJets["conePTSumCharged"]] =  conePTSumCharged;		
		TupleVariables[VariableJets["conePxSumCharged"]] =  conePxSumCharged;			
		TupleVariables[VariableJets["conePySumCharged"]] =  conePySumCharged;
		TupleVariables[VariableJets["conePzSumCharged"]] =  conePzSumCharged;
		TupleVariables[VariableJets["conePTSumInclusive"]] = conePTSumInclusive; 		
		TupleVariables[VariableJets["conePxSumInclusive"]] = conePxSumInclusive;  				
		TupleVariables[VariableJets["conePySumInclusive"]] = conePySumInclusive; 		
		TupleVariables[VariableJets["conePzSumInclusive"]] = conePzSumInclusive;
		TupleVariables[VariableJets["conePTSumNeutral"]] = conePTSumNeutral; 		
		TupleVariables[VariableJets["conePxSumNeutral"]] = conePxSumNeutral;  				
		TupleVariables[VariableJets["conePySumNeutral"]] = conePySumNeutral; 		
		TupleVariables[VariableJets["conePzSumNeutral"]] = conePzSumNeutral;
		TupleVariables[VariableJets["Q1_vec"]] = fQ_1;                               		
		TupleVariables[VariableJets["Q2_vec"]] = fQ_2;                               		
		TupleVariables[VariableJets["Psi_2"]] = fPsi_2;        
		TupleVariables[VariableJets["Q1_vec_rec"]] = fQ_1_rec;                               		
		TupleVariables[VariableJets["Q2_vec_rec"]] = fQ_2_rec;                               		
		TupleVariables[VariableJets["Psi_2_shifted"]] = fPsi_2_shifted;  
		                       		
                
                Jets->Fill(TupleVariables);
	

        //Delete vectors to avoid a pileup
        //inclusive_jets.clear();
        input_particles.clear();
        chargedjetTracks.clear();
        neutraljetTracks.clear();

      //} //End of D0 loop

  //} //End of IsthereD0 condition

  //Delete daughter particles vectors
 // DaughterPionTrackVector.clear();
 // DaughterKaonTrackVector.clear();

  //End of the event
  return kStOK;
}

//---------------------------------------------------------------------
//-----------------------------FUNCTIONS-------------------------------
//---------------------------------------------------------------------
 //Not used anymore
Double_t StPicoBackgroundMaker::vertexCorrectedEta(double eta, double vz) {
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
Bool_t StPicoBackgroundMaker::GetCaloTrackMomentum(StPicoDst *mPicoDst, TVector3 mPrimVtx) {
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
      ////float dca_z = dcaPoint.z() - vtx_z; //Test
      
      float dca_z = trk->gDCAz(vtx_z); //? TO DO
      
      //Exclude tracks with a DCA to the primary vertex in z of more than maxdcazhadroncorr (in RunPicoD0AnaMaker.C)
      if (fabs(dca_z) > maxdcazhadroncorr) continue;

      //Initialization and loading of the tower index
      int TowIndex = -99999;
      TowIndex = trk->bemcTowerIndex(); //ID
      float p = 0;
      
      //Check if the track is matched to a tower
      if (TowIndex >= 0) {

        //Loading of the momentum
        p = gMom.Mag();          
	double TrackEnergy = 1.0*TMath::Sqrt(p*p + M_PION_PLUS*M_PION_PLUS);
        
        //Summing up the energy of all tracks matched to the same tower //Previously neglected pion mass 
        SumE[TowIndex] += TrackEnergy;
      }
  
  } //End of track loop

  return true;

  //Function returns true if it was successful and SumE filled with the momentum of all tracks matched to a tower
}
//---------------------------------------------------------------------------
Double_t StPicoBackgroundMaker::GetTowerCalibEnergy(Int_t TowerId){
  //Function calculates the calibrated energy of a tower

  //Loading of the tower
  StPicoBTowHit *tower = static_cast<StPicoBTowHit*>(picoDst->btowHit(TowerId-1));

  //Initialization of the pedestal, rms and status
  Float_t pedestal, rms;
  Int_t status;

  //Loading of the pedestal, rms and status (it does not work, if you use root instead of root4star)
  mTables->getPedestal(BTOW, TowerId, 0, pedestal, rms); //softID
  mTables->getStatus(BTOW, TowerId, status); //softID

  //Initialization of the tower coefficients
  Double_t *TowerCoeff;

  //Tower coefficients for the different runs, parameters are saved in BemcNewCalib.h
  if(picoDst->event()->runId() <= 15094020) TowerCoeff = CPre;
  else TowerCoeff = CLowMidHigh;

  //Calculation of the calibrated energy E=C*(ADC-Pedestal)
  Double_t calibEnergy = TowerCoeff[TowerId-1]*(tower->adc() - pedestal); //softID

  return calibEnergy;

  //Function returns the calibrated energy of the tower
}
//---------------------------------------------------------------------------
bool StPicoBackgroundMaker::IsBadEnergyRun(int runID) {
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
bool StPicoBackgroundMaker::isInRndCone(double eta, double phi, double eta_rnd, double phi_rnd) const{


    // Výpočet rozdílů v eta a phi
    double dEta = eta - eta_rnd;
    double dPhi = phi - phi_rnd;

    // Zajištění správné periody phi v rozsahu [-pi, pi]
    if (dPhi > TMath::Pi()) {
        dPhi -= 2 * TMath::Pi();
    } else if (dPhi < -TMath::Pi()) {
        dPhi += 2 * TMath::Pi();
    }

    // Výpočet metrické vzdálenosti v eta-phi prostoru
    double deltaR = TMath::Sqrt(dEta * dEta + dPhi * dPhi);

    // Kontrola vzdálenosti
    return (deltaR < 0.4);
}

//---------------------------------------------------------------------------
int StPicoBackgroundMaker::isD0PairCentrality_pt(StKaonPion const* const kp, int Centrality, int mYear) const{
    //Check if the pair passes the cuts for D0

    //Loading the daughter particles tracks
    StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
    StPicoTrack const* pion = picoDst->track(kp->pionIdx());

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
bool StPicoBackgroundMaker::isGoodEvent(int mYear,  TH1F* NEventsCuts){
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
  
  //Check on suspicious all-0 position
  bool nonezeroVertex = (event->primaryVertex().x()!=0 && event->primaryVertex().y()!=0 && event->primaryVertex().z()!=0);
  if (!nonezeroVertex) return false;
  NEventsCuts->Fill(5);

  return true;

  //Function returns true if the event passes the cuts
}
//-----------------------------------------------------------------------------
bool StPicoBackgroundMaker::isMBTrigger(int mYear){
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
bool StPicoBackgroundMaker::isGoodTrack(StPicoTrack const * const trk) const{
    // Require at least one hit on every layer of PXL and IST.
    // It is done here for tests on the preview II data.
    // The new StPicoTrack which is used in official production has a method to check this

    //Check if the track meets the HFT requirement
    //2014 - Require at least one hit on every layer of PXL and IST
    //2016 - Require at least one hit on every layer of PXL and (IST or SST)
    //Both can be written in te same way
    bool HFTCondition = (trk->hasPxl1Hit() && trk->hasPxl2Hit()) && (trk->hasSstHit() || trk->hasIstHit());

    //Check if |eta| < 1
    bool EtaCondition = abs(trk->gMom().PseudoRapidity()) < 1;

    //In StCuts.cxx is defined if the HFT is required and the nHitsFit and minPt values.
    return (trk->gPt() > mycuts::minPt && trk->nHitsFit() >= mycuts::nHitsFit && (HFTCondition || !mycuts::requireHFT) && EtaCondition);



    //Function returns true if track is good
}
//-----------------------------------------------------------------------------
bool StPicoBackgroundMaker::isGoodJetTrack(StPicoTrack const * const trk,StPicoEvent const *const myEvent) const{
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
bool StPicoBackgroundMaker::isGoodJetTrackSim(TVector3 simTrack, Int_t reco, Double_t simDca) const{
    //Check if the track meets the jet track cuts
    //Parameters saved in StCuts.cxx

    //pT range cut
    bool pTTrackJetCut = simTrack.Perp() > mycuts::jetTrackPtMin && simTrack.Perp() < mycuts::jetTrackPtMax;
    //eta range cut
    bool etaTrackJetCut = fabs(simTrack.PseudoRapidity()) < mycuts::jetTrackEta;
    //nHitsFit cut
    bool nHitsTrackJetCut = abs(Track_mNHitsFit[reco]) >= mycuts::jetTracknHitsFit;
    //nHitsRatio cut
    bool nHitsRatioTrackJetCut = (1.0 * abs(Double_t(Track_mNHitsFit[reco]) / Double_t(Track_mNHitsMax[reco]))) >= mycuts::jetTracknHitsRatio;
    //DCA cut
    bool dcaTrackJetCut = simDca < mycuts::jetTrackDCA;

    return pTTrackJetCut && etaTrackJetCut && nHitsTrackJetCut && nHitsRatioTrackJetCut && dcaTrackJetCut;

    //Return true if all the cuts are passed
}
//-----------------------------------------------------------------------------
bool StPicoBackgroundMaker::isTpcPion(StPicoTrack const * const trk) const{
    //Check if the track meets nsigma (TPC) requirement for pion
    //Parameters are saved in StCuts.cxx

    return fabs(trk->nSigmaPion()) < mycuts::nSigmaPion;

    //Function returns true if track is good
}
//-----------------------------------------------------------------------------
bool StPicoBackgroundMaker::isTpcKaon(StPicoTrack const * const trk) const{
    //Check if the track meets nsigma (TPC) requirement for kaon
    //Parameters are saved in StCuts.cxx

    return fabs(trk->nSigmaKaon()) < mycuts::nSigmaKaon; //In D0 event maker it is set to 2

    //Function returns true if track is good
}
//-----------------------------------------------------------------------------
float StPicoBackgroundMaker::getTofBeta(StPicoTrack const * const trk, StThreeVectorF const* const pVtx,StPicoDst const* const picoDst) const{
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
bool StPicoBackgroundMaker::isTofKaon(StPicoTrack const * const trk, float beta) const{
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
//-----------------------------------------------------------------------------
bool StPicoBackgroundMaker::isTofPion(StPicoTrack const * const trk, float beta) const{
    //Check if the track meets |1/beta-1/beta_pi| (TOF) requirement for pion

    //Initialization of tofPion
    bool tofPion = false;
    
    //If beta is positive, than we can calculate |1/beta-1/beta_K|
    if(beta>0){

        //Calculation of the global total momentum
        double ptot = trk->gMom().Mag();
        //Calculation of the expected beta for kaons
        float beta_pi = ptot/sqrt(ptot*ptot+M_PION_PLUS*M_PION_PLUS);
        //Check if the track meets the TOF requirement
        //Parameters are saved in StCuts.cxx
        tofPion = fabs(1/beta - 1/beta_pi) < mycuts::pTofBetaDiff ? true : false;
    } //End of beta > 0

    return tofPion;

    //Function returns true if track is good based on TOF information
}





