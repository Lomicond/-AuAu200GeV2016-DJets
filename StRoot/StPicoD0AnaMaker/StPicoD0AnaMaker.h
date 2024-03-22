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
#include "TProfile.h"
#include "StiMaker/StKFVerticesCollection.h"
#include "StPhysicalHelixD.hh"
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
    int mEventCounter;
    StPicoDstMaker* mPicoDstMaker;
    int mYear;
    StPicoDst *picoDst;
    TString mInputFileList;
    TFile* mOutputFile;
    StPicoD0Event* mPicoD0Event;
    StPicoPrescales* mPrescales;
    StRefMultCorr* mGRefMultCorrUtil;
    StHFCuts* mHFCuts;

  public:
    StPicoD0AnaMaker(
            char const * name,
            char const * inputFilesList,
            char const * outName,
            StPicoDstMaker* picoDstMaker,
            StRefMultCorr* grefmultCorrUtil,
            int pYear
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
    bool isTpcPion(StPicoTrack const*) const;
    bool isTpcKaon(StPicoTrack const*,StThreeVectorF const * pVtx) const;
    bool isTofKaon(StPicoTrack const* const, float beta) const;
    float getTofBeta(StPicoTrack const*,StThreeVectorF const * pVtx,StPicoDst const* const picoDst) const;
    bool IsBadEnergyRun(int);





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

    //Jet backround
    TH2D* Jet_grefmult_pt_background;
    TH2D* Jet_D0pT_vs_D0rapidity;
    TH2D* Jet_D0pT_vs_Jetrapidity;
    TH1D* Jet_phi;

    //TNtuple D0-jets
    Float_t TupleVariables[24];
    TNtuple* Jets;
    std::map<std::string, int> VariableJets;

    //2D D0 mass-pt like-sign and unlike-sign histograms
    TH2D *massPt;
    TH2D *massPtLike;

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
