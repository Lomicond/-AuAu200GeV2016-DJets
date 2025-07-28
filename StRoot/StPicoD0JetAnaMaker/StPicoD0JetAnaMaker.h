#ifndef StPicoD0JetAnaMaker_h
#define StPicoD0JetAnaMaker_h

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
 **s
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

//#ifdef FASTJET_VERSION

#ifndef __CINT__
namespace fastjet { 
	class PseudoJet; 
	class ClusterSequenceArea;
}
#endif


class StPrimaryVertex; 
class StEvent;
class StDcaGeometry; 

class TString;
class TFile;
class StPicoEvent;
class StKaonPion;
class StPicoDstMaker;
class StPicoDst;
class StPicoTrack;
class StRefMultCorr;
class StEmcADCtoEMaker;
class StBemcTables;



class StPicoD0JetAnaMaker : public StMaker
{
  private:

    TString mOutFileName;
    TFile* mOutputFile;
    StRefMultCorr* mGRefMultCorrUtil;
    StPicoDstMaker* mPicoDstMaker;
    Int_t mYear;
    StPicoDst *picoDst;

  public:
    StPicoD0JetAnaMaker(
            char const * name,
            char const * outName,
            StPicoDstMaker* picoDstMaker,
            StRefMultCorr* grefmultCorrUtil,
            Int_t pYear
    );
    virtual ~StPicoD0JetAnaMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();

    void setMaxDcaZHadronCorr(Float_t max);
    Int_t getEntries() const;
    void setHadronCorr(Float_t corr);
    void setOnlyTrackBasedJets(Bool_t onlyTrackBasedJets);
    void setCutETmin(Float_t min);
    void setGhostMaxrap(Float_t fGhostMaxrap);
    void setNJetsRemove(Int_t nJetsRemove);
    void setMaxNeutralFraction(Float_t max);
  void CalculateEventPlane();
    Int_t EP_IsGoodTrack(StPicoTrack *trk, TVector3 pVertex);
    virtual Bool_t GetCaloTrackMomentum(StPicoDst *mPicoDst, TVector3 mPrimVtx);
    virtual Double_t GetTowerCalibEnergy(Int_t TowerId);
    virtual Double_t vertexCorrectedEta(Double_t eta, Double_t vz);

    StEmcADCtoEMaker *mADCtoEMaker;
    StBemcTables     *mTables;

  private:

    StPicoD0JetAnaMaker() {}
    void readNextEvent();
    Int_t isD0PairCentrality_pt(StKaonPion const & kp, Int_t Centrality, Int_t mYear) const;
    Bool_t isGoodEvent(Int_t mYear, TH1D* hEventsCuts);
    Bool_t isMBTrigger(Int_t mYear);
    Bool_t isGoodTrack(StPicoTrack const*) const;
    Bool_t isGoodJetTrack(StPicoTrack const*,StPicoEvent const*) const; //my
    Bool_t isTpcPion(StPicoTrack const*) const;
    Bool_t isTpcKaon(StPicoTrack const*) const;
    Bool_t isTofKaon(StPicoTrack const* const, Float_t beta) const;
    Bool_t isTofPion(StPicoTrack const* const, Float_t beta) const;
    Float_t getTofBeta(StPicoTrack const*,StThreeVectorF const * pVtx,StPicoDst const* const picoDst) const;
    Bool_t IsBadEnergyRun(Int_t);
    #ifndef __CINT__
    std::vector<fastjet::PseudoJet> JetReconstructionICS(std::vector<fastjet::PseudoJet> fInputVectors, Int_t fCentrality, Double_t fR, Bool_t fBackSub, Double_t fEP_psi2, Int_t difiter);
    fastjet::ClusterSequenceArea* fClustSeq = nullptr;
    #endif
    Float_t fGhostMaxrap;
    Int_t nJetsRemove;
    Float_t maxneutralfrac;
    Float_t fETmincut;
    Float_t fHadronCorr;
    Float_t OnlyTrackBasedJets;
    Float_t maxdcazhadroncorr;
    const Double_t mBarrelRadius = 225.405;

TTree* Jets;


//Event
Int_t runId;
Int_t eventId;
Int_t centrality;
Float_t weightCentrality;
Int_t gRefMult;
Float_t backgroundDensity;
Float_t psi2;

//D0 meson
Int_t d0PdgSign;
Float_t d0Mass;
Float_t d0Pt;
Float_t d0Rapidity;
Float_t d0Eta;

//Jet obsevables
Float_t jetEta;
Float_t jetPhi;
Float_t jetRapidity;
Float_t jetArea;
Float_t jetPt; 
Float_t lambda1_0_5;
Float_t lambda1_1;
Float_t lambda1_1_5;
Float_t lambda1_2;
Float_t lambda1_3;
Float_t momDisp;
Float_t z;
Int_t nJetConst;
Int_t nJetsInEvent;
Float_t jetD0DeltaR;
Float_t jetNeutralPtFrac;
Float_t jetTrackPtSum;

// Event histograms:
TH1D* hVtxZ;
TH2D* hVtxR;
TH1D* hVzDiff;
TH1D* hCentrality;
TH1D* hCentralityW;
TH1D* hEventsCuts;

// D0 histograms:
TH2D* hD0MassPtUnlike;
TH2D* hD0MassPtLike;
TH1D* hD0EtaUnlike;
TH1D* hD0EtaLike;
TH2D* hPionEtaVsPt;
TH2D* hKaonEtaVsPt;
TH2D* hNKaonsVsNPions;

// Topological cuts:
TH2D* hKaonDcaVsPtD0;
TH2D* hPionDcaVsPtD0;
TH2D* hDcaDaughtersVsPt;
TH2D* hD0DcaVsPt;
TH2D* hCosThetaVsPt;
TH2D* hD0DecayLengthVsPt;

// Daughter PID:
TH2D* hTofBetaDiffKaonVsPt;
TH2D* hTofBetaDiffPionVsPt;
TH2D* hBetaVsSPKaon;
TH2D* hBetaVsSPPion;
TH2D* hTpcNsigmaKaonVsPt;
TH2D* hTpcNsigmaPionVsPt;
TH2D* hDedxVsSPKaon;
TH2D* hDedxVsSPPion;
TH2D* hNHitsFitKaonVsD0Pt;
TH2D* hNHitsFitPionVsD0Pt;

//Jet constituents:
TH2D* hJetTracksDedx;
TH2D* hJetTracksDedxAfterCuts;
TH1D* hJetTracksPt;
TH2D* hJetTracksEtaPhi;
TH1D* hJetTracksNHitsFit;
TH1D* hJetTracksNHitsRatio;
TH1D* hJetTracksDCA;
TH1D* hJetNeutralPt;
TH2D* hJetNeutralEtaPhi;
TH2D* hJetNeutralEtBefAftHC;
TH2D* hJetNeutralECalibBefAft;
TH1D* hJetConstCharge;
TH2D* hJetConstRapPhi;
TH2D* hJetConstRapPhiICS;
TH2D* hJetConstEtaPhi;
TH2D* hJetConstEtaPhiICS;

//Neutral particles hadronic correction:
TH1D* hJetHadrCorrNHitsFit;
TH1D* hJetHadrCorrNHitsRatio;
TH1D* hJetHadrCorrDcaZ;
TH2D* hJetHadrCorrEtaVsPt;
TH1D* hJetHadrCorrE;


//Event plane calculation
Double_t fPsi_2;
Double_t fPsi_2_shifted;
Double_t fQ_1;
Double_t fQ_2;
Double_t fQ_1_rec;
Double_t fQ_2_rec;

ClassDef(StPicoD0JetAnaMaker, 1)

};

inline void StPicoD0JetAnaMaker::setCutETmin(Float_t min){
    fETmincut = min;
}

inline void StPicoD0JetAnaMaker::setMaxDcaZHadronCorr(Float_t max) {
    maxdcazhadroncorr = max;
}

inline void StPicoD0JetAnaMaker::setGhostMaxrap(Float_t fGhostMaxrap) {
    StPicoD0JetAnaMaker::fGhostMaxrap = fGhostMaxrap;
}

inline void StPicoD0JetAnaMaker::setNJetsRemove(Int_t nJetsRemove) {
    StPicoD0JetAnaMaker::nJetsRemove = nJetsRemove;
}

inline void StPicoD0JetAnaMaker::setMaxNeutralFraction(Float_t max) {
    maxneutralfrac = max;
}

inline void StPicoD0JetAnaMaker::setHadronCorr(Float_t corr) {
    fHadronCorr = corr;
}

inline void StPicoD0JetAnaMaker::setOnlyTrackBasedJets(Bool_t OTBJets) {
    OnlyTrackBasedJets = OTBJets;
}

#endif
