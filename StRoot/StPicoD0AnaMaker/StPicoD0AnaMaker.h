#ifndef StPicoD0AnaMaker_h
#define StPicoD0AnaMaker_h

/***********************************************************************************
 **
 ** Taken from D0CorrelationV2Analyser Author: Leon He
 ** Modified by: Ondrej Lomicky
 **
 ************************************************************************************
 **
 ** Description: 
 **
 ************************************************************************************
 **
 ** Log:
 **
 ********************************************
 *  A Maker to read a StPicoEvent and StPicoD0Event
 *  simultaneously and do analysis. 
 *
 *  Please write your analysis in the ::Make() function.
 *
 *  Authors:  Xin Dong        (xdong@lbl.gov)
 *            Michael Lomnitz (mrlomnitz@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *            Jochen Thaeder  (jmthader@lbl.gov)   
 *
 * **************************************************
 */
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/projection/StEmcPosition.h"
#include "TChain.h"
#include "StMaker.h"
//
#include "StThreeVectorF.hh"
#include "TSpectrum.h"
#include "Math/Functor.h"
#include "Math/GSLMinimizer1D.h"
#include "TCanvas.h"
#include "TH1K.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TProfile.h"
#include "StiMaker/StKFVerticesCollection.h"
#include "StPhysicalHelixD.hh"
#include "TLorentzVector.h"


////#include <fastjet/PseudoJet.hh>  // Pro práci s PseudoJet
///#include <fastjet/ClusterSequence.hh>  // Pro clustering
///#include <fastjet/Selector.hh>  // Pro selektory
///#include <fastjet/tools/JetMedianBackgroundEstimator.hh>  // Pro Background Estimation

///using namespace fastjet;  // Alias pro jmenný prostor fastjet



class StPrimaryVertex; 
class StEvent;
class StDcaGeometry; 
//// 

class TString;
class TFile;
class TNtuple;
class StPicoEvent;
class StPicoD0Event;
class StKaonPion;
class StPicoDstMaker;
class StPicoDst;
class StPicoTrack;
class StHFCuts;
class StPicoPrescales;
class StRefMultCorr;
class StEmcADCtoEMaker;
class StBemcTables;



class StPicoD0AnaMaker : public StMaker
{
  private:

    TString mOutFileName;
    TChain* mChain;
    TFile* mOutputFile;
    StRefMultCorr* mGRefMultCorrUtil;
    StPicoDstMaker* mPicoDstMaker;
    StPicoD0Event* mPicoD0Event;
    TString mInputFileList;
    int mEventCounter;
    int mYear;
    StPicoDst *picoDst;
    StPicoPrescales* mPrescales;
    StHFCuts* mHFCuts;
    bool mSimulation;

  public:
    StPicoD0AnaMaker(
            char const * name,
            char const * inputFilesList,
            char const * outName,
            StPicoDstMaker* picoDstMaker,
            StRefMultCorr* grefmultCorrUtil,
            int pYear,
            bool Sim,
            char const *filename
    );
    virtual ~StPicoD0AnaMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();

    void setMaxDcaZHadronCorr(float max);
    int getEntries() const;
    void setHadronCorr(float corr);
    void setOnlyTrackBasedJets(bool onlyTrackBasedJets);
    void setHFCuts(StHFCuts* cuts);   
    void setCutETmin(float min);
    void setGhostMaxrap(float fGhostMaxrap);
    void setNJetsRemove(int nJetsRemove);
    void setMaxNeutralFraction(float max);
  void CalculateEventPlane();
  Int_t EP_IsGoodTrack(StPicoTrack *trk, TVector3 pVertex);
    virtual Bool_t GetCaloTrackMomentum(StPicoDst *mPicoDst, TVector3 mPrimVtx);
    virtual Double_t GetTowerCalibEnergy(Int_t TowerId);
    virtual Double_t vertexCorrectedEta(double eta, double vz);

    StEmcADCtoEMaker *mADCtoEMaker;
    StBemcTables     *mTables;

  private:

    StPicoD0AnaMaker() {}
    void readNextEvent();
    ofstream fout;
    int isD0PairCentrality_pt(StKaonPion const* const kp, int Centrality, int mYear) const;
    bool isGoodEvent(int mYear, TH1F* NEventsCuts);
    bool isMBTrigger(int mYear);
    bool isGoodTrack(StPicoTrack const*) const;
    bool isGoodJetTrack(StPicoTrack const*,StPicoEvent const*) const; //my
    bool isGoodJetTrackSim(TVector3 simTrack, Int_t reco, Double_t simDca) const; //my
    bool LoadMcChargedTracks(const Int_t &iD0) const; //my
    bool isTpcPion(StPicoTrack const*) const;
    bool isTpcKaon(StPicoTrack const*) const;
    bool isTofKaon(StPicoTrack const* const, float beta) const;
    bool isTofPion(StPicoTrack const* const, float beta) const;
    float getTofBeta(StPicoTrack const*,StThreeVectorF const * pVtx,StPicoDst const* const picoDst) const;
    bool IsBadEnergyRun(int);
    void ReadTreeMc();
    TVector3 FastSimMom(TVector3 p, Int_t pid) const;
    Bool_t KeepTrack(const Int_t & particleid, const Int_t & centralitybin, const Double_t &  pt) const;
    void GetAllTracksFromVertex(const Int_t &vertexid, vector<Int_t> &trackvec) const;
    ///void JetReconstruction(vector<fastjet::PseudoJet> &fInputVectors, int fCentrality) const;

    float fGhostMaxrap;
    int nJetsRemove;
    float maxneutralfrac;
    float fETmincut;
    float fHadronCorr;
    float OnlyTrackBasedJets;
    float maxdcazhadroncorr;
    const double mBarrelRadius = 225.405;

    //Events histograms:
    TH1F *vtxz;
    TH2D *vtxr;
    TH1D *hcentr;
    TH1D *hcentrW;
    TH1F* NEventsCuts;
    TH1F* mh1TotalEventsInRun;
    TH1F* mh1TotalHftTracksInRun;
    TH1F* mh1TotalGRefMultInRun;

    //D0 histograms:
    TH1D *D0etalike;
    TH1D *D0etaunlike;
    TH1D *pioneta;
    TH1D *kaoneta;
    TH1F* mh1TotalKaonsInRun;
    TH1F* mh1TotalPionsInRun;
    TH1F* mh1TotalD0CandidatesInRun;
    TH2F* mh2NKaonsVsNPions;
    TH2F* mh2KaonDcaVsPt;
    TH2F* mh2PionDcaVsPt;
    TH2F* mh2CosThetaVsPt;
    TH2F* mh2DcaDaughtersVsPt;

    //Jet Tracks
    TH2D* JetTracksdEdx;
    TH2D* JetTracksdEdxCut;
    TH1F* hpT_tr;
    TH2F* heta_phi_tr;
    TH1F* heta_tr;
    TH1F* hphi_tr;
    TH2F* hdca_z_tr;
    TH2F* hdca_pT;
    TH1F* hdca_tr;
    TH1F* hcharged_tr;
    TH1D* hDcaSign[9][6][2];

    //Jet backround
    TH2D* Jet_grefmult_pt_background;
    TH2D* Jet_D0pT_vs_D0rapidity;
    TH2D* Jet_D0pT_vs_Jetrapidity;
    TH1D* Jet_phi;

    //TNtuple D0-jets
    Float_t TupleVariables[27];
    TNtuple* Jets[4];
    TNtuple* EventStats;
    std::map<std::string, int> VariableJets;

    //2D D0 mass-pt like-sign and unlike-sign histograms
    TH2D *massPt;
    TH2D *massPtLike;
    
    //Simulation
    vector<TString> filenamesforHIOverlay;
    TFile *f;
    TTree *fMCPico;
        
      static const Int_t kMaxNumberOfD0Events = 100;
  vector<TLorentzVector> fRecoMcEventTracks[kMaxNumberOfD0Events];     // For charged tracks
    vector<TLorentzVector> fRecoMcEventTowers[kMaxNumberOfD0Events];     // For neutral towers
      vector<TLorentzVector> fMcEventTracks[kMaxNumberOfD0Events];         // For charged particles
   vector<TLorentzVector> fMcEventTowers[kMaxNumberOfD0Events];         // For neutral particles   
     pair<TVector3, TVector3> fMcD0Information[kMaxNumberOfD0Events];     // Pion momenta, Kaon momenta (3 components each)
  pair<TVector3, TVector3> fMcRecoD0Information[kMaxNumberOfD0Events]; // Pion momenta, Kaon momenta (3 components each)
    TVector3 fOrigin[kMaxNumberOfD0Events];                              // MC Event Origin Information
     pair<int, int> fMcEventInfo[kMaxNumberOfD0Events];                   // RunID, EventID
     
     
       // variables
  Int_t fRunNumber;
  Int_t centralitybinforefficiency;
  Double_t fRhoVal;
  vector<Int_t> vertexids;
  vector<Int_t> pionids;
  vector<Int_t> kaonids;
  vector<Int_t> matchedpionids;
  vector<Int_t> matchedkaonids;
  map<int, vector<Int_t>> fVertexToTracks;

  Double_t fPsi_2;
  Double_t fPsi_2_shifted;
  Double_t fQ_1;
  Double_t fQ_2;
  Double_t fQ_1_rec;
  Double_t fQ_2_rec;

    TF1 *fKaonMomResolution;
    TF1 *fPionMomResolution;
    TF1 *fProtonMomResolution;

    TGraph *fPionWeight[3];
    TGraph *fKaonWeight[3];
    TGraph *fProtonWeight[3];
    TGraph *fAProtonWeight[3];

    vector<Int_t> fDroppedMCTracks; // Dropped MC Tracks (This will be tracks which came from the KPi from D0 decaying. We don't want them.)
      // Fixed size dimensions of array or collections stored in the TTree if any.
	  static const Int_t kMaxEvent = 1;
	  static const Int_t kMaxTrack = 1000;
	  static const Int_t kMaxEmcTrigger = 89;
	  static const Int_t kMaxMtdTrigger = 1;
	  static const Int_t kMaxBTowHit = 4800;

	  static const Int_t kMaxMcVertex = 6000;
	  static const Int_t kMaxMcTrack = 6000;
    
	  Int_t Event_;
	  Int_t Event_mRunId[kMaxEvent];   //[Event_]
	  Int_t Event_mEventId[kMaxEvent]; //[Event_]
	  Float_t Event_mPrimaryVertexX[kMaxEvent];
	  Float_t Event_mPrimaryVertexY[kMaxEvent];
	  Float_t Event_mPrimaryVertexZ[kMaxEvent];

	  Int_t Track_;
	  Float_t Track_mGMomentumX[kMaxTrack]; //[Track_]
	  Float_t Track_mGMomentumY[kMaxTrack]; //[Track_]
	  Float_t Track_mGMomentumZ[kMaxTrack]; //[Track_]
	  Float_t Track_mOriginX[kMaxTrack];    //[Track_]
	  Float_t Track_mOriginY[kMaxTrack];    //[Track_]
	  Float_t Track_mOriginZ[kMaxTrack];    //[Track_]

	  Char_t Track_mNHitsFit[kMaxTrack];  //[Track_]
	  UChar_t Track_mNHitsMax[kMaxTrack]; //[Track_]

	  Short_t Track_mNSigmaPion[kMaxTrack];            //[Track_]
	  Short_t Track_mNSigmaKaon[kMaxTrack];            //[Track_]
	  Short_t Track_mNSigmaProton[kMaxTrack];          //[Track_]
	  Short_t Track_mNSigmaElectron[kMaxTrack];        //[Track_]
	  Short_t Track_mBEmcMatchedTowerIndex[kMaxTrack]; //[Track_]
	  UShort_t Track_mIdTruth[kMaxTrack];              //[Track_]
	  UShort_t Track_mQATruth[kMaxTrack];              //[Track_]

	  Int_t BTowHit_;
	  Short_t BTowHit_mE[kMaxBTowHit]; //[BTowHit_]

	  Int_t McVertex_;
	  Int_t McVertex_mId[kMaxMcVertex];             //[McVertex_]
	  UShort_t McVertex_mNoDaughters[kMaxMcVertex]; //[McVertex_]
	  Int_t McVertex_mIdParTrk[kMaxMcVertex];       //[McVertex_]
	  Int_t McVertex_mIsInterm[kMaxMcVertex];       //[McVertex_]
	  Float_t McVertex_mTime[kMaxMcVertex];         //[McVertex_]
	  Float_t McVertex_mVx[kMaxMcVertex];           //[McVertex_]
	  Float_t McVertex_mVy[kMaxMcVertex];           //[McVertex_]
	  Float_t McVertex_mVz[kMaxMcVertex];           //[McVertex_]
	  Int_t McTrack_;
	  UShort_t McTrack_mId[kMaxMcTrack];         //[McTrack_]
	  Int_t McTrack_mGePid[kMaxMcTrack];         //[McTrack_]
	  Char_t McTrack_mCharge[kMaxMcTrack];       //[McTrack_]
	  UChar_t McTrack_mHits[kMaxMcTrack][22];    //[McTrack_]
	  Float_t McTrack_mPx[kMaxMcTrack];          //[McTrack_]
	  Float_t McTrack_mPy[kMaxMcTrack];          //[McTrack_]
	  Float_t McTrack_mPz[kMaxMcTrack];          //[McTrack_]
	  Float_t McTrack_mE[kMaxMcTrack];           //[McTrack_]
	  Bool_t McTrack_mIsFromShower[kMaxMcTrack]; //[McTrack_]
	  Short_t McTrack_mIdVtxStart[kMaxMcTrack];  //[McTrack_]
	  Short_t McTrack_mIdVtxStop[kMaxMcTrack];   //[McTrack_]
	  Short_t McTrack_mIdVtxItrmd[kMaxMcTrack];  //[McTrack_]

	  // List of branches
	  TBranch *b_Event_;                //!
	  TBranch *b_Event_mRunId;          //!
	  TBranch *b_Event_mEventId;        //!
	  TBranch *b_Event_mPrimaryVertexX; //!
	  TBranch *b_Event_mPrimaryVertexY; //!
	  TBranch *b_Event_mPrimaryVertexZ; //!

	  TBranch *b_Track_;            //!
	  TBranch *b_Track_mGMomentumX; //!
	  TBranch *b_Track_mGMomentumY; //!
	  TBranch *b_Track_mGMomentumZ; //!
	  TBranch *b_Track_mOriginX;    //!
	  TBranch *b_Track_mOriginY;    //!
	  TBranch *b_Track_mOriginZ;    //!

	  TBranch *b_Track_mNHitsFit;              //!
	  TBranch *b_Track_mNHitsMax;              //!
	  TBranch *b_Track_mNSigmaPion;            //!
	  TBranch *b_Track_mNSigmaKaon;            //!
	  TBranch *b_Track_mNSigmaProton;          //!
	  TBranch *b_Track_mNSigmaElectron;        //!
	  TBranch *b_Track_mTopologyMap;           //!
	  TBranch *b_Track_mBEmcMatchedTowerIndex; //!
	  TBranch *b_Track_mIdTruth;               //!
	  TBranch *b_Track_mQATruth;               //!

	  TBranch *b_BTowHit_;   //!
	  TBranch *b_BTowHit_mE; //!

	  TBranch *b_McVertex_;             //!
	  TBranch *b_McVertex_mId;          //!
	  TBranch *b_McVertex_mNoDaughters; //!
	  TBranch *b_McVertex_mIdParTrk;    //!
	  TBranch *b_McVertex_mIsInterm;    //!
	  TBranch *b_McVertex_mTime;        //!
	  TBranch *b_McVertex_mVx;          //!
	  TBranch *b_McVertex_mVy;          //!
	  TBranch *b_McVertex_mVz;          //!
	  TBranch *b_McTrack_;              //!
	  TBranch *b_McTrack_mId;           //!
	  TBranch *b_McTrack_mGePid;        //!
	  TBranch *b_McTrack_mCharge;       //!
	  TBranch *b_McTrack_mHits;         //!
	  TBranch *b_McTrack_mPx;           //!
	  TBranch *b_McTrack_mPy;           //!
	  TBranch *b_McTrack_mPz;           //!
	  TBranch *b_McTrack_mE;            //!
	  TBranch *b_McTrack_mIsFromShower; //!
	  TBranch *b_McTrack_mIdVtxStart;   //!
	  TBranch *b_McTrack_mIdVtxStop;    //!
	  TBranch *b_McTrack_mIdVtxItrmd;   //!

  TString fMCFileListName;

Int_t fNumberOfEventsToOverLay = 1; ///THIS


    ClassDef(StPicoD0AnaMaker, 1)

};

inline int StPicoD0AnaMaker::getEntries() const {
  return mChain? mChain->GetEntries() : 0;
}

inline void StPicoD0AnaMaker::setCutETmin(float min){
    fETmincut = min;
}

inline void StPicoD0AnaMaker::setMaxDcaZHadronCorr(float max) {
    maxdcazhadroncorr = max;
}

inline void StPicoD0AnaMaker::setGhostMaxrap(float fGhostMaxrap) {
    StPicoD0AnaMaker::fGhostMaxrap = fGhostMaxrap;
}

inline void StPicoD0AnaMaker::setNJetsRemove(int nJetsRemove) {
    StPicoD0AnaMaker::nJetsRemove = nJetsRemove;
}

inline void StPicoD0AnaMaker::readNextEvent(){
    mChain->GetEntry(mEventCounter++);
}

inline void StPicoD0AnaMaker::setMaxNeutralFraction(float max) {
    maxneutralfrac = max;
}

inline void StPicoD0AnaMaker::setHadronCorr(float corr) {
    fHadronCorr = corr;
}

inline void StPicoD0AnaMaker::setHFCuts(StHFCuts* cuts) {
    mHFCuts = cuts;
}

inline void StPicoD0AnaMaker::setOnlyTrackBasedJets(bool OTBJets) {
    OnlyTrackBasedJets = OTBJets;
}

#endif
