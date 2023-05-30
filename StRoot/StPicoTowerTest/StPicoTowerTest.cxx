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

#include "StPicoTowerTest.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoBTowHit.h"
#include "StCuts.h"

#include "../StPicoPrescales/StPicoPrescales.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StEmcADCtoEMaker/StBemcData.h"
#include "StEmcADCtoEMaker/StEmcADCtoEMaker.h"
#include "StEmcRawMaker/StBemcRaw.h"
#include "StEmcRawMaker/defines.h"
#include "StEmcRawMaker/StBemcTables.h"
#include "BemcNewCalib.h"
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
#include <set>
#include <vector>
#include <stdio.h>
#include <time.h>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cmath>
using namespace std; 

//-------------------------------------

ClassImp(StPicoTowerTest)

  StPicoTowerTest::StPicoTowerTest(char const * name, char const * outName,StPicoDstMaker* picoDstMaker,StRefMultCorr* grefmultCorrUtil,int pYear): 
    StMaker(name),mPicoDstMaker(picoDstMaker), mGRefMultCorrUtil(grefmultCorrUtil),
    mOutFileName(outName),mOutputFile(NULL), mChain(NULL), mEventCounter(0),mYear(pYear){}

Int_t StPicoTowerTest::Init()
{

  mPrescales = new StPicoPrescales(mycuts::prescalesFilesDirectoryName);


  mOutputFile = new TFile(mOutFileName.Data(), "RECREATE");
  mOutputFile->cd();

  int nRuns = mPrescales->numberOfRuns();
  MapTower_RunID_NHits_TowerID = new TH2D("MapTower_RunID_NHits_TowerID","MapTower_RunID_NHits_TowerID;runIndex;towerID",nRuns+1,0,nRuns+1,4800,0,4800);
  MapTower_RunID_Energy_TowerID = new TH2D("MapTower_RunID_Energy_TowerID","MapTower_RunID_Energy_TowerID;runIndex;towerID",nRuns+1,0,nRuns+1,4800,0,4800);
  MapTower_RunID_NEvents = new TH1D("MapTower_RunID_NEvents", "MapTower_RunID_NEvents", nRuns+1,0,nRuns+1); // 80 x 60 bins
  vtxr = new TH2D("vtxr",";PVtx.x() [cm]; PVtx.y() [cm]",100,-3,3,100,-3,3);

 
 
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
  //mGRefMultCorrUtil = new StRefMultCorr("grefmult");

  return kStOK;
}
//-----------------------------------------------------------------------------
StPicoTowerTest::~StPicoTowerTest()
{
  /*  */
  delete mGRefMultCorrUtil;
}

//-----------------------------------------------------------------------------
Int_t StPicoTowerTest::Finish()
{
  LOG_INFO << " StPicoTowerTest - writing data and closing output file " <<endm;
  fout.close();
  fout1.close();
  mOutputFile->cd();

  vtxr->Write();
  MapTower_RunID_NHits_TowerID->Write();
  MapTower_RunID_Energy_TowerID->Write();
  MapTower_RunID_NEvents->Write();

  mOutputFile->Close();
  delete mPrescales;

  return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoTowerTest::Make()
{

  //readNextEvent();
  if (!mPicoDstMaker)
  {
    LOG_WARN << " StPicoTowerTest - No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  picoDst = mPicoDstMaker->picoDst();

  if (!picoDst)
  {
    LOG_WARN << "StPicoTowerTest - No PicoDst! Skip! " << endm;
    return kStWarn;
  }

  //cout << "rok: " << mYear << endl;

  StThreeVectorF pVtx(-999.,-999.,-999.);
  StPicoEvent *event = (StPicoEvent *)picoDst->event();

  if(!(isGoodEvent()))//minBias trigger requires
  {
    return kStOK;
  }
  
  if(!mGRefMultCorrUtil) {
    LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
    return kStWarn;
  }

  mGRefMultCorrUtil->init(picoDst->event()->runId());
  mGRefMultCorrUtil->initEvent(picoDst->event()->grefMult(),pVtx.z(),picoDst->event()->ZDCx()) ;

  int centrality  = mGRefMultCorrUtil->getCentralityBin9();

  if(centrality<0) {
    //LOG_WARN << "not minBias sample!" << endl;
    return kStOK;
  }

  pVtx = StThreeVectorF(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z());
  vtxr->Fill(pVtx.x(),pVtx.y());

  int runIndex = mPrescales->runIndex(picoDst->event()->runId());
  MapTower_RunID_NEvents->Fill(runIndex);

      for (int iTow = 0; iTow < 4800; iTow++){ //get btow info
        
          StPicoBTowHit *towHit = picoDst->btowHit(iTow);
          vector<int> ids = {0,0,0,0,0,0,0,0,0}; 
          if (!towHit) continue; 
          //if (!towHit || towHit->isBad()) continue; //if the tower is bad or missing info
          int realtowID = towHit->numericIndex2SoftId(iTow);


          StPicoBTowHit *tower = static_cast<StPicoBTowHit*>(picoDst->btowHit(iTow));
          Float_t pedestal, rms;
          Int_t status;
          mTables->getPedestal(BTOW, iTow+1, 0, pedestal, rms);
          mTables->getStatus(BTOW, iTow+1, status);  
          Double_t *TowerCoeff;

          if(picoDst->event()->runId() <= 15094020) TowerCoeff = CPre;
          else TowerCoeff = CLowMidHigh;
  
          StEmcGeom* mEmcGeom;
          mEmcGeom = StEmcGeom::getEmcGeom("bemc");
              
          float Toweta_tmp = 0, Towphi = 0;
          mEmcGeom->getEtaPhi(realtowID,Toweta_tmp,Towphi);

          float Toweta = vertexCorrectedEta(Toweta_tmp, event->primaryVertex().z());
          double E_T = tower->energy()/cosh(Toweta);
          // cout << "E_T: " << E_T << endl;
          if(E_T>0.2){
             MapTower_RunID_NHits_TowerID->Fill(runIndex,iTow);
             MapTower_RunID_Energy_TowerID->Fill(runIndex, iTow, tower->energy());
          }
  
      }

  return kStOK;
}
//-----------------------------------------------------------------------------

//----------------------------------------------------------------------------- 
//Correct tower eta for Vz position
//----------------------------------------------------------------------------- 
Double_t StPicoTowerTest::vertexCorrectedEta(double eta, double vz) {
    double tower_theta = 2.0 * atan(exp(-eta));
    double z = 0.0;
    if (eta != 0.0) z = mBarrelRadius / tan(tower_theta);
    double z_diff = z - vz;
    double theta_corr = atan2(mBarrelRadius, z_diff);
    double eta_corr = -log(tan(theta_corr / 2.0));
    return eta_corr;
}



bool StPicoTowerTest::isGoodEvent()
{
  StPicoEvent *event = (StPicoEvent *)picoDst->event();
  // return (event->triggerWord() & mycuts::triggerWord) &&
  //cout << "Test: " << isMBTrigger(mYear) << endl;
  return (isMBTrigger(mYear) &&
      sqrt(event->primaryVertex().x()*event->primaryVertex().x()+event->primaryVertex().y()*event->primaryVertex().y()) < mycuts::vz &&
      fabs(event->primaryVertex().z()) < mycuts::vz &&
      fabs(event->primaryVertex().z() - event->vzVpd()) < mycuts::vzVpdVz);
  //  return event->triggerWord() & mycuts::triggerWord;
}
bool StPicoTowerTest::isMBTrigger(int mYear)
{

  const std::set<int>* mbTriggers = nullptr;

  if(mYear ==2016){
      static const std::set<int> mbTriggers2016 = {
        520001, 520011, 520021, 520031, 520041, 520051,           // VPDMB-5-p-sst
        570002, 570001,                                           // VPDMB-5-nosst  (production 2, nosst stream), VPDMB-5-sst (production 2, sst stream )
        520201, 520211, 520221, 520231, 520241, 520251, 520261,   // BHT1*VPDMB-10
        520203,                                                   // BHT
        520101, 520111, 520121, 520131, 520141,                   // central-5  
        520007, 520017, 520027, 520037,                           // vpdmb-10
        520003, 520013, 520023, 520033, 520043,                   // VPDMB-5
        520802, 520812, 520822, 520832, 520842,                   // VPDMB-5-p-hlt
        520002, 520012, 520022, 520032, 520042                    // VPDMB-5-p-nosst
      };
      mbTriggers = &mbTriggers2016;
  }

  if(mYear ==2014){  
      static const std::set<int> mbTriggers2014 = {
      450050, 450060, 450005, 450015, 450025
      };
      mbTriggers = &mbTriggers2014;
  }

    StPicoEvent* event = static_cast<StPicoEvent*>(picoDst->event());
  return std::any_of(mbTriggers->begin(), mbTriggers->end(), [&](int trigger) { return event->isTrigger(trigger); });
}
//-----------------------------------------------------------------------------
bool StPicoTowerTest::isGoodTrack(StPicoTrack const * const trk) const
{
  // Require at least one hit on every layer of PXL and IST.
  // It is done here for tests on the preview II data.
  // The new StPicoTrack which is used in official production has a method to check this
  //return trk->gPt() > mycuts::minPt && trk->nHitsFit() >= mycuts::nHitsFit && trk->isHFTTrack();

  bool HFTCondition = (trk->hasPxl1Hit() && trk->hasPxl2Hit()) && (trk->hasSstHit() || trk->hasIstHit());
  return trk->gPt() > mycuts::minPt && trk->nHitsFit() >= mycuts::nHitsFit && HFTCondition;
  //return  trk->nHitsFit() >= mycuts::nHitsFit;
}
















