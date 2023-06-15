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
#include "StPicoPrescales/StPicoPrescales.h"
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
        char const * outName,StPicoDstMaker* picoDstMaker,StRefMultCorr* grefmultCorrUtil, int pYear);
    virtual ~StPicoTowerTest();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();
    void setMaxDcaZHadronCorr(float max);
    int getEntries() const;
    void setHadronCorr(float corr);
    void setCutETmin(float min); 
    ofstream fout1;
    virtual Double_t vertexCorrectedEta(double eta, double vz);
    StEmcADCtoEMaker *mADCtoEMaker;
    StBemcTables     *mTables;
    private:
    int mYear;
    StPicoTowerTest() {}
    void readNextEvent();
    ofstream fout;

    //-------------------------------------------

    bool isGoodEvent();
    bool isMBTrigger(int mYear);
    bool isGoodTrack(StPicoTrack const*) const;

    StPicoDstMaker* mPicoDstMaker;
    StPicoPrescales* mPrescales;
    StRefMultCorr* mGRefMultCorrUtil;
    StPicoDst *picoDst;

    TString mOutFileName;
    TString mInputFileList;
    TFile* mOutputFile;
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

    TH2D *vtxr;
    TH2D* MapTower_RunID_NHits_TowerID;
    TH2D* MapTower_RunID_Energy_TowerID;
    TH1D *MapTower_RunID_NEvents;


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



#endif
