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
#include <vector>
#include <stdio.h>
#include <time.h>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cmath>
using namespace std;  //robotmon

//-------------------------------------

ClassImp(StPicoTowerTest)

  StPicoTowerTest::StPicoTowerTest(char const * name, char const * outName,StPicoDstMaker* picoDstMaker,StRefMultCorr* grefmultCorrUtil): 
    StMaker(name),mPicoDstMaker(picoDstMaker), mGRefMultCorrUtil(grefmultCorrUtil),
    mOutFileName(outName),mOutputFile(NULL), mChain(NULL), mEventCounter(0){}

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
    LOG_WARN << "not minBias sample!" << endl;
    return kStWarn;
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
  return (isMBTrigger() &&
      sqrt(event->primaryVertex().x()*event->primaryVertex().x()+event->primaryVertex().y()*event->primaryVertex().y()) < mycuts::vz &&
      fabs(event->primaryVertex().z()) < mycuts::vz &&
      fabs(event->primaryVertex().z() - event->vzVpd()) < mycuts::vzVpdVz);
  //  return event->triggerWord() & mycuts::triggerWord;
}
bool StPicoTowerTest::isMBTrigger()
{
  StPicoEvent *event = (StPicoEvent *)picoDst->event();
  return (event->isTrigger(520001) ||// VPDMB-5-p-sst (production 1, physics stream)||
      event->isTrigger(520011)|| // VPDMB-5-p-sst
      event->isTrigger(520021)|| // VPDMB-5-p-sst
      event->isTrigger(520031)|| // VPDMB-5-p-sst
      event->isTrigger(520041)|| // VPDMB-5-p-sst
      event->isTrigger(520051)|| // VPDMB-5-p-sst
      event->isTrigger(570002)|| // VPDMB-5-nosst (production 2, nosst stream)
      event->isTrigger(570001)||  // VPDMB-5-sst (production 2, sst stream )
      event->isTrigger(520201)|| //BHT1*VPDMB-10
      event->isTrigger(520211)|| //BHT1*VPDMB-10
      event->isTrigger(520221)|| //BHT1*VPDMB-10
      event->isTrigger(520231)|| //BHT1*VPDMB-10
      event->isTrigger(520241)|| //BHT1*VPDMB-10
      event->isTrigger(520251)|| //BHT1*VPDMB-10
      event->isTrigger(520261)|| //BHT1*VPDMB-10
      event->isTrigger(520203)|| //BHT
      event->isTrigger(520101)|| //central-5
      event->isTrigger(520111)|| //central-5
      event->isTrigger(520121)|| //central-5
      event->isTrigger(520131)|| //central-5
      event->isTrigger(520141)|| //central-5
      event->isTrigger(520007)|| //vpdmb-10
      event->isTrigger(520017)|| //vpdmb-10
      event->isTrigger(520027)|| //vpdmb-10
      event->isTrigger(520037)|| //vpdmb-10
      event->isTrigger(520003)|| //VPDMB-5
      event->isTrigger(520013)|| //VPDMB-5
      event->isTrigger(520023)|| //VPDMB-5
      event->isTrigger(520033)|| //VPDMB-5
      event->isTrigger(520043)|| //VPDMB-5
      event->isTrigger(520802)|| //VPDMB-5-p-hlt
      event->isTrigger(520812)|| //VPDMB-5-p-hlt
      event->isTrigger(520822)|| //VPDMB-5-p-hlt
      event->isTrigger(520832)|| //VPDMB-5-p-hlt
      event->isTrigger(520842)|| //VPDMB-5-p-hlt
      event->isTrigger(520002)|| //VPDMB-5-p-nosst
      event->isTrigger(520012)|| //VPDMB-5-p-nosst
      event->isTrigger(520022)|| //VPDMB-5-p-nosst
      event->isTrigger(520032)|| //VPDMB-5-p-nosst
      event->isTrigger(520042)); //VPDMB-5-p-nosst            
      
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
















