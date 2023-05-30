#include <vector>
#include <cmath>
#include <algorithm>
#include <string>

#include <set>
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "StThreeVectorF.hh"
#include "StLorentzVectorF.hh"
#include "StPhysicalHelixD.hh"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoD0Event.h"
#include "StPicoD0Hists.h"
#include "StKaonPion.h"
#include "StCuts.h"

#include "StPicoD0EventMaker.h"

ClassImp(StPicoD0EventMaker)



StPicoD0EventMaker::StPicoD0EventMaker(char const* makerName, StPicoDstMaker* picoMaker, char const* fileBaseName,int pYear)
   : StMaker(makerName), mPicoDstMaker(picoMaker), mPicoEvent(NULL), mPicoD0Hists(NULL), mOutputFile(NULL), mTree(NULL), mPicoD0Event(NULL),mYear(pYear)

{
   mPicoD0Event = new StPicoD0Event();

   TString baseName(fileBaseName);
   mOutputFile = new TFile(Form("%s.picoD0.root",fileBaseName), "RECREATE");
   mOutputFile->SetCompressionLevel(1);
   int BufSize = (int)pow(2., 16.);
   int Split = 1;
   mTree = new TTree("T", "T", BufSize);
   mTree->SetAutoSave(1000000); // autosave every 1 Mbytes
   mTree->Branch("dEvent", "StPicoD0Event", &mPicoD0Event, BufSize, Split);

   mPicoD0Hists = new StPicoD0Hists(fileBaseName);
}

StPicoD0EventMaker::~StPicoD0EventMaker()
{
   /* mTree is owned by mOutputFile directory, it will be destructed once
    * the file is closed in ::Finish() */
   delete mPicoD0Hists;
}

Int_t StPicoD0EventMaker::Init()
{
   return kStOK;
}

Int_t StPicoD0EventMaker::Finish()
{
   mOutputFile->cd();
   mOutputFile->Write();
   mOutputFile->Close();
   mPicoD0Hists->closeFile();
   //mKfVertexEvent.closeFile();
   return kStOK;
}

void StPicoD0EventMaker::Clear(Option_t *opt)
{
   mPicoD0Event->clear("C");
}

Int_t StPicoD0EventMaker::Make()
{
  // cout << "huh " <<mYear<< endl;
   if (!mPicoDstMaker)
   {
      LOG_WARN << " No PicoDstMaker! Skip! " << endm;
      return kStWarn;
   }

   StPicoDst const * picoDst = mPicoDstMaker->picoDst();
   if (!picoDst)
   {
      LOG_WARN << " No PicoDst! Skip! " << endm;
      return kStWarn;
   }
   //cout << "huh2" << endl;
   mPicoEvent = picoDst->event();


   unsigned int nHftTracks = 0;
   
   //cout << "hmmm" << endl;
   if (isGoodEvent(mYear))
   {
      //cout << "hmmm2" << endl;
      UInt_t nTracks = picoDst->numberOfTracks();

      std::vector<unsigned short> idxPicoKaons;
      std::vector<unsigned short> idxPicoPions;

      StThreeVectorF const pVtx = StThreeVectorF(mPicoEvent->primaryVertex().x(),mPicoEvent->primaryVertex().y(),mPicoEvent->primaryVertex().z());

      for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack)
      {
         StPicoTrack* trk = picoDst->track(iTrack);
         //cout << "aaa" << endl;
         if(!trk) continue;

      
         if (!isGoodTrack(trk)) continue;
         ++nHftTracks; 

         if (isPion(trk)) idxPicoPions.push_back(iTrack);
         if (isKaon(trk)) idxPicoKaons.push_back(iTrack);

      } // .. end tracks loop

        mPicoD0Event->nKaons(idxPicoKaons.size());
        mPicoD0Event->nPions(idxPicoPions.size());

        float const bField = mPicoEvent->bField();
        for (unsigned short ik = 0; ik < idxPicoKaons.size(); ++ik)
        {

          StPicoTrack const * kaon = picoDst->track(idxPicoKaons[ik]);

          // make Kπ pairs
          for (unsigned short ip = 0; ip < idxPicoPions.size(); ++ip)
          {
            
            if (idxPicoKaons[ik] == idxPicoPions[ip]) continue;
            StPicoTrack const * pion = picoDst->track(idxPicoPions[ip]);

            StKaonPion kaonPion(kaon, pion, idxPicoKaons[ik], idxPicoPions[ip], pVtx, bField);

            if (!isGoodPair(kaonPion)) continue;
            if(isGoodMass(kaonPion)) mPicoD0Event->addKaonPion(&kaonPion);
            bool fillMass = isGoodQaPair(&kaonPion,*kaon,*pion);

            /*
            if (fillMass){
            cout<<"⣿⣿⣿⣿⣿⣿⠿⢋⣥⣴⣶⣶⣶⣬⣙⠻⠟⣋⣭⣭⣭⣭⡙⠻⣿⣿⣿⣿⣿"<<endl;
            cout<<"⣿⣿⣿⣿⡿⢋⣴⣿⣿⠿⢟⣛⣛⣛⠿⢷⡹⣿⣿⣿⣿⣿⣿⣆⠹⣿⣿⣿⣿"<<endl;
            cout<<"⣿⣿⣿⡿⢁⣾⣿⣿⣴⣿⣿⣿⣿⠿⠿⠷⠥⠱⣶⣶⣶⣶⡶⠮⠤⣌⡙⢿⣿"<<endl;
            cout<<"⣿⡿⢛⡁⣾⣿⣿⣿⡿⢟⡫⢕⣪⡭⠥⢭⣭⣉⡂⣉⡒⣤⡭⡉⠩⣥⣰⠂⠹"<<endl;
            cout<<"⡟⢠⣿⣱⣿⣿⣿⣏⣛⢲⣾⣿⠃⠄⠐⠈⣿⣿⣿⣿⣿⣿⠄⠁⠃⢸⣿⣿⡧"<<endl;
            cout<<"⢠⣿⣿⣿⣿⣿⣿⣿⣿⣇⣊⠙⠳⠤⠤⠾⣟⠛⠍⣹⣛⣛⣢⣀⣠⣛⡯⢉⣰"<<endl;
            cout<<"⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⡶⠶⢒⣠⣼⣿⣿⣛⠻⠛⢛⣛⠉⣴⣿⣿"<<endl;
            cout<<"⣿⣿⣿⣿⣿⣿⣿⡿⢛⡛⢿⣿⣿⣶⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⡈⢿⣿"<<endl;
            cout<<"⣿⣿⣿⣿⣿⣿⣿⠸⣿⡻⢷⣍⣛⠻⠿⠿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠿⢇⡘⣿"<<endl;
            cout<<"⣿⣿⣿⣿⣿⣿⣿⣷⣝⠻⠶⣬⣍⣛⣛⠓⠶⠶⠶⠤⠬⠭⠤⠶⠶⠞⠛⣡⣿"<<endl;
            cout<<"⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⣶⣬⣭⣍⣙⣛⣛⣛⠛⠛⠛⠿⠿⠿⠛⣠⣿⣿"<<endl;
            cout<<"⣦⣈⠉⢛⠻⠿⠿⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⠿⠛⣁⣴⣾⣿⣿⣿⣿"<<endl;
            cout<<"⣿⣿⣿⣶⣮⣭⣁⣒⣒⣒⠂⠠⠬⠭⠭⠭⢀⣀⣠⣄⡘⠿⣿⣿⣿⣿⣿⣿⣿"<<endl;
            cout<<"⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣦⡈⢿⣿⣿⣿⣿⣿"<<endl;
            }
            */

            bool unlike = kaon->charge() * pion->charge() < 0 ? true : false;
            //cout << "fillMass: " << fillMass << " unlike: " << unlike << endl;
            if(fillMass || unlike) mPicoD0Hists->addKaonPion(&kaonPion,fillMass, unlike);

          } // .. end make Kπ pairs
        } // .. end of kaons loop
   } //.. end of good event fill


   mPicoD0Event->addPicoEvent(*mPicoEvent);
   mPicoD0Hists->addEvent(*mPicoEvent,*mPicoD0Event,nHftTracks);

   // This should never be inside the good event block
   // because we want to save header information about all events, good or bad
   mTree->Fill();
   mPicoD0Event->clear("C");

   return kStOK;
}

bool StPicoD0EventMaker::isGoodEvent(int mYear)
{
   /*
   return (mPicoEvent->triggerWord() & cuts::triggerWord) &&
          fabs(mPicoEvent->primaryVertex().z()) < cuts::vz &&
          fabs(mPicoEvent->primaryVertex().z() - mPicoEvent->vzVpd()) < cuts::vzVpdVz;
          */


   return   isMinBiasTrigger(mYear) &&
            sqrt(mPicoEvent->primaryVertex().x()*mPicoEvent->primaryVertex().x()+mPicoEvent->primaryVertex().y()*mPicoEvent->primaryVertex().y()) < cuts::vz &&
            fabs(mPicoEvent->primaryVertex().z()) < cuts::vz &&
            fabs(mPicoEvent->primaryVertex().z() - mPicoEvent->vzVpd()) < cuts::vzVpdVz;

/*
      return  sqrt(mPicoEvent->primaryVertex().x()*mPicoEvent->primaryVertex().x()+mPicoEvent->primaryVertex().y()*mPicoEvent->primaryVertex().y()) < cuts::vz &&
            fabs(mPicoEvent->primaryVertex().z()) < cuts::vz &&
            fabs(mPicoEvent->primaryVertex().z() - mPicoEvent->vzVpd()) < cuts::vzVpdVz;
*/
}

bool StPicoD0EventMaker::isMinBiasTrigger(int mYear)
{

      const std::set<int>* mbTriggers = nullptr;

      if(mYear ==2016) mbTriggers = &cuts::mbTriggers2016;
      if(mYear ==2014) mbTriggers = &cuts::mbTriggers2014;

      StPicoEvent* event = static_cast<StPicoEvent*>(mPicoDstMaker->picoDst()->event());
      return std::any_of(mbTriggers->begin(), mbTriggers->end(), [&](int trigger) { return event->isTrigger(trigger); });
   

}

bool StPicoD0EventMaker::isGoodTrack(StPicoTrack const * const trk) const
{
   // Require at least one hit on every layer of PXL and IST.
   // It is done here for tests on the preview II data.
   // The new StPicoTrack which is used in official production has a method to check this
   //cout << " trk->isHft() " << trk->isHft() << endl;
   //cout << " trk->nHitsFit() " << trk->nHitsFit() << ">=" << cuts::nHitsFit <<endl;
/*
   cout << " trk->hasPxl1Hit() " << trk->hasPxl1Hit() << endl;
   cout << " trk->hasPxl2Hit() " << trk->hasPxl2Hit() << endl;
   cout << " trk->hasIstHit() " << trk->hasIstHit() << endl;
   cout << " trk->hasSstHit() " << trk->hasSstHit() << endl;
   */
   bool HFTCondition = false;
   if(mYear ==2016) HFTCondition = (trk->hasPxl1Hit() && trk->hasPxl2Hit()) && (trk->hasSstHit() || trk->hasIstHit());
   if(mYear ==2014) HFTCondition = trk->isHFTTrack();                             //Bool_t isHFTTrack()const { return hasPxl1Hit() && hasPxl2Hit() && (hasIstHit() || hasSstHit()); }
   return (!cuts::requireHFT || HFTCondition) &&  trk->nHitsFit() >= cuts::nHitsFit;

}
bool StPicoD0EventMaker::isPion(StPicoTrack const * const trk) const
{
  // cout <<"fabs(trk->nSigmaPion()): " <<fabs(trk->nSigmaPion())<< endl;
   return fabs(trk->nSigmaPion()) < cuts::nSigmaPion;
}
bool StPicoD0EventMaker::isKaon(StPicoTrack const * const trk) const
{
   return fabs(trk->nSigmaKaon()) < cuts::nSigmaKaon;
}
bool StPicoD0EventMaker::isGoodPair(StKaonPion const & kp) const
{
  return std::cos(kp.pointingAngle()) > cuts::cosTheta &&
         kp.decayLength() > cuts::decayLength &&
         kp.dcaDaughters() < cuts::dcaDaughters;
/*
   cout << "std::cos(kp.pointingAngle()): " << std::cos(kp.pointingAngle()) << endl;
   cout << "> cuts::cosTheta: " << cuts::cosTheta << endl;

   cout << "kp.decayLength(): " << kp.decayLength() << endl;
 cout << "cuts::decayLength: " << cuts::decayLength << endl;

cout << "kp.dcaDaughters(): " << kp.dcaDaughters() << endl;
cout << "cuts::dcaDaughters: " << cuts::dcaDaughters << endl;

cout << "--------------" << endl;*/
   //return 1;
}

bool StPicoD0EventMaker::isGoodMass(StKaonPion const & kp) const
{
/*
   cout << "kp.m(): " << kp.m() << endl;
   cout << "cuts::minMass: " << cuts::cosTheta << endl;

   cout << "kp.m(): " << kp.m() << endl;
 cout << "cuts::maxMass: " << cuts::maxMass << endl;
cout << "--------------" << endl;
*/
   return kp.m() > cuts::minMass && kp.m() < cuts::maxMass;
}

bool  StPicoD0EventMaker::isGoodQaPair(StKaonPion const& kp, StPicoTrack const& kaon,StPicoTrack const& pion)
{
   /*
   cout << "pion.gPt(): " << pion.gPt() << endl;
   cout << ">= " << cuts::qaPt << endl;
   cout << "kaon.gPt(): " << kaon.gPt() << endl;
   cout << ">= " << cuts::qaPt << endl;
   cout << "pion.nHitsFit(): " << pion.nHitsFit() << endl;
   cout <<  ">= " << cuts::qaNHitsFit << endl;
   cout << "kaon.nHitsFit(): " << kaon.nHitsFit() << endl;
   cout << ">= "<< cuts::qaNHitsFit << endl;
   cout << "fabs(kaon.nSigmaKaon()): " << fabs(kaon.nSigmaKaon()) << endl;
   cout << "< "<< cuts::qaNSigmaKaon << endl;
   cout << "cos(kp.pointingAngle()): " << cos(kp.pointingAngle()) << endl;
   cout << "> "<< cuts::qaCosTheta << endl;
   cout << "kp.pionDca(): " << kp.pionDca() << endl;
   cout << "> "<< cuts::qaPDca << endl;
   cout << "kp.kaonDca()(): " << kp.kaonDca() << endl;
   cout << "> "<< cuts::qaKDca << endl;
   cout << "kp.dcaDaughters(): " << kp.dcaDaughters() << endl;
   cout << "< "<< cuts::qaDcaDaughters << endl;
   cout << "-----------------" << endl;*/
  return pion.gPt() >= cuts::qaPt && kaon.gPt() >= cuts::qaPt && 
         pion.nHitsFit() >= cuts::qaNHitsFit && kaon.nHitsFit() >= cuts::qaNHitsFit &&
         fabs(kaon.nSigmaKaon()) < cuts::qaNSigmaKaon && 
         cos(kp.pointingAngle()) > cuts::qaCosTheta &&
         kp.pionDca() > cuts::qaPDca && kp.kaonDca() > cuts::qaKDca &&
         kp.dcaDaughters() < cuts::qaDcaDaughters;
}
