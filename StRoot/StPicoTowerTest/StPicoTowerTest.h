#ifndef StPicoTowerTest_h
#define StPicoTowerTest_h

/***********************************************************************************
 **
 ** D0CorrelationV2Analyser
 **
 ** Author: Leon He
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
class StPicoDstMaker;
class StPicoDst;
class StPicoTrack;
class StHFCuts;
class StPicoPrescales;
class StRefMultCorr;
class StEmcADCtoEMaker;
class StBemcTables;

class StPicoTowerTest : public StMaker
{
  public:
    StPicoTowerTest(char const * name, 
        char const * outName,StPicoDstMaker* picoDstMaker,StRefMultCorr* grefmultCorrUtil);
    virtual ~StPicoTowerTest();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();
  void setMaxDcaZHadronCorr(float max);
    int getEntries() const;
    void setHadronCorr(float corr);
    void setHFCuts(StHFCuts* cuts);   
    void setCutETmin(float min); 
    ofstream fout1;
        virtual Double_t GetTowerCalibEnergy(Int_t TowerId);
          virtual Double_t vertexCorrectedEta(double eta, double vz);
  StEmcADCtoEMaker *mADCtoEMaker;
  StBemcTables     *mTables;
  private:
    StPicoTowerTest() {}
    void readNextEvent();
    ofstream fout;

    //-------------------------------------------
    int D0Reco(StThreeVectorF *);
    bool isGoodEvent();
    bool isMBTrigger();
    bool  isGoodTrack(StPicoTrack const*) const;
    bool  isGoodTrack2(StPicoTrack const*) const; //my
    bool  isGoodHadron(StPicoTrack const*) const;
    bool  isTpcPion(StPicoTrack const*) const;
    bool  isTpcKaon(StPicoTrack const*,StThreeVectorF const * pVtx) const;
    bool isTofKaon(StPicoTrack const* const, float beta) const;
    float getTofBeta(StPicoTrack const*,StThreeVectorF const * pVtx,StPicoDst const* const picoDst) const;
    bool getCorHadron(float eta, vector<float> &hadronsPhi, const vector<unsigned int>, float, float);
    float sumCos(float phi, vector<float> &hadronsPhi);
    bool fixPhi(vector<float> &phi);
    bool getHadronCorV2(int );
    //bool getCorV2(int , double, int &);
    bool getCorV2(int , double);
    bool isEtaGap(double, double ,double);

    StPicoDstMaker* mPicoDstMaker;
    StPicoD0Event* mPicoD0Event;
    StPicoPrescales* mPrescales;
    StRefMultCorr* mGRefMultCorrUtil;
    //StPicoDstMaker *
    StPicoDst *picoDst;

    TString mOutFileName;
    TString mInputFileList;
    TFile* mOutputFile;
    //TFile* mPhi;
    TChain* mChain;
    int mEventCounter;
      const double mBarrelRadius = 225.405;
    StHFCuts* mHFCuts;
    float fETmincut;
  // hadronic correction fraction
  float fHadronCorr;
    // -------------- USER variables -------------------------
    // add your member variables here. 
    // Remember that ntuples size can be really big, use histograms where appropriate

    TNtuple *mEventTuple;
    TNtuple *mDTuple;
    TNtuple *mHadronTuple;
    TH1F *dEtaDHadron;
    TH1F *hEtaD;
    TH1F *hEtaHadron;
    TH2F *hPhiHadron[8][3];
    TH2F *hPhiD[8][3];
    TH1F *vtxz;
    TH2D *vtxr;
    TH1D *NEvent;
    TH1D *hcentr;
    TH1D *D0etalike;
    TH1D *D0etaunlike;
    TNtuple* Jets;
   TProfile* TProfileTower_ADC;
TH2D *MapTower;
TH2D *MapTower_ADC;
TH2D* MapTower_RunID_NHits_TowerID;
TH2D* MapTower_RunID_Energy_TowerID;
TH2D* MapTower_RunID_Energy_TowerID2;

TNtuple* MapTower_RunID_ADC_TowerID;
TH1D *MapTower_RunID_NEvents;
   TH2D *MapTower_pedestal; 
   TH2D *MapTower_ADCmPedestal; 
   TH2D *MapTower_corrected; 
   TH3D *MapTower_vs_RunID;

   TH1F* mh1TotalEventsInRun;
   TH1F* mh1TotalHftTracksInRun;
   TH1F* mh1TotalGRefMultInRun;
   TH1F* mh1TotalKaonsInRun;
   TH1F* mh1TotalPionsInRun;
   TH1F* mh1TotalD0CandidatesInRun;
   TH2F* mh2NKaonsVsNPions;
   TH2F* mh2KaonDcaVsPt;
   TH2F* mh2PionDcaVsPt;
   TH2F* mh2CosThetaVsPt;
   TH2F* mh2DcaDaughtersVsPt;
   TH2F* mh2InvariantMassVsPtUnlike;
   TH2F* mh2InvariantMassVsPtLike;

 TH1F* hpT_tr;
 TH2F* heta_phi_tr;
 TH1F* heta_tr;
 TH1F* hphi_tr;
 TH2F* hdca_z_tr;
 TH2F* hdca_pT;
 TH1F* hdca_tr;
 TH1F* hcharged_tr;

    TH1D *pioneta;
    TH1D *kaoneta;
    TH2F *etaPhi;
    TH2F *etaPhi_D;
    TH2F *etaPhi_Hadron;
    TH2F *etaPhi_Hadron_all;
    TProfile *profV2[8][5][3];//i.S or B; j.flatten; k. differetn etaGap
    TH1D *hadronV2[5][3];
    TH1D *hadronV2_sum[5][3];
    TH1D *hadronV2_excl[5][9][3];
    TH1D *hadronV2_excl_sum[5][9][3];
    TH2D *fitPhi[6];
    TH2D *massPt;
    TH2D *massPtLike;
    TH2D *massLike;
    TH2D *massLike2;
    TH2D *massUnlike;
    TH2D *v2Weight[8][3];
    TH2D *likeV2Mass[6][5];
    TH2D *likeV2Mass2[6][5];
    TH2D *unlikeV2Mass[6][5];
    TProfile *V2Mass[2][6][5];
    TProfile *candPt;
    TNtuple *checkNew;
    TNtuple *checkOld;

    TH2D *checkPeak;

    double efficiency[4][6];
    ClassDef(StPicoTowerTest, 1)

    float maxdcazhadroncorr;
};

inline int StPicoTowerTest::getEntries() const 
{
  return mChain? mChain->GetEntries() : 0;
}
inline void StPicoTowerTest::setCutETmin(float min)            { fETmincut = min;}

inline void StPicoTowerTest::setMaxDcaZHadronCorr(float max)           { maxdcazhadroncorr = max;}

inline void StPicoTowerTest::readNextEvent()
{
  mChain->GetEntry(mEventCounter++);
}
inline void StPicoTowerTest::setHadronCorr(float corr)   { fHadronCorr = corr;} 

inline void StPicoTowerTest::setHFCuts(StHFCuts* cuts)   
{ 
  mHFCuts = cuts; 
}

#endif
