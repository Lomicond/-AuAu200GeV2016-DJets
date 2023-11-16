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
   : StMaker(makerName), mPicoDstMaker(picoMaker), mPicoEvent(NULL), mPicoD0Hists(NULL), mOutputFile(NULL), mTree(NULL), mPicoD0Event(NULL),mYear(pYear){

   mPicoD0Event = new StPicoD0Event();

   //Setting up the name of the output file
   TString baseName(fileBaseName);

   //Creating the output file
   mOutputFile = new TFile(Form("%s.picoD0.root",fileBaseName), "RECREATE");

   //Creating the output tree
   mOutputFile->SetCompressionLevel(1);
   int BufSize = (int)pow(2., 16.);
   int Split = 1;
   mTree = new TTree("T", "T", BufSize);

   //Setting the branch address
   mTree->SetAutoSave(1000000); // autosave every 1 Mbytes
   mTree->Branch("dEvent", "StPicoD0Event", &mPicoD0Event, BufSize, Split);

   //Creating the histograms
   mPicoD0Hists = new StPicoD0Hists(fileBaseName);
}

StPicoD0EventMaker::~StPicoD0EventMaker(){
   /* mTree is owned by mOutputFile directory, it will be destructed once
    * the file is closed in ::Finish() */
   delete mPicoD0Hists;
}

Int_t StPicoD0EventMaker::Init(){
   return kStOK;
}

Int_t StPicoD0EventMaker::Finish(){
   //Saving all root files.
   mOutputFile->cd();
   mOutputFile->Write();
   mOutputFile->Close();
   mPicoD0Hists->closeFile();

   return kStOK;
}

void StPicoD0EventMaker::Clear(Option_t *opt){
   mPicoD0Event->clear("C");
}

Int_t StPicoD0EventMaker::Make(){

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
   mPicoEvent = picoDst->event();

   //Set initial number of HFT tracks to 0
   unsigned int nHftTracks = 0;

   //Check if the event is a good event
   if (isGoodEvent(mYear)){

      //Get number of tracks in event
      UInt_t nTracks = picoDst->numberOfTracks();

      //Vectors to hold indices of kaons and pions
      std::vector<unsigned short> idxPicoKaons;
      std::vector<unsigned short> idxPicoPions;

      //Calculation of primary vertex variables
      StThreeVectorF const pVtx = StThreeVectorF(mPicoEvent->primaryVertex().x(),mPicoEvent->primaryVertex().y(),mPicoEvent->primaryVertex().z());

      //Loop over all tracks in PicoDst to find good kaon and pion candidates
      for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack)
      {
         //Loading track information
         StPicoTrack* trk = picoDst->track(iTrack);

         //If the track is null, skip it
         if(!trk) continue;

         //Check if the track is good for use in analysis
         if (!isGoodTrack(trk)) continue;

         //HFT track counter
         ++nHftTracks; 

         //Check if the track is a kaon or pion
         //If you want to find different particles, you can change the isPion() and isKaon() functions there and cuts in StCuts.cxx
         if (isPion(trk)) idxPicoPions.push_back(iTrack);
         if (isKaon(trk)) idxPicoKaons.push_back(iTrack);

      } // ... end tracks loop

      //Set number of kaon and pion candidates in event
      mPicoD0Event->nKaons(idxPicoKaons.size());
      mPicoD0Event->nPions(idxPicoPions.size());

      //Get magnetic field
      float const bField = mPicoEvent->bField();

      //Loops over all kaon and pion candidates to make D0 candidates
      //Loop over all kaon candidates
      for (unsigned short ik = 0; ik < idxPicoKaons.size(); ++ik){

        //Get kaon track information
        StPicoTrack const * kaon = picoDst->track(idxPicoKaons[ik]);

        //Loop over all pion candidates
        for (unsigned short ip = 0; ip < idxPicoPions.size(); ++ip)
        {

            //If pion and kaon are the same track, skip it
            //Sometimes it can happen that the same track is identified as both kaon and pion
            //We do not want to pair particle with itself
            if (idxPicoKaons[ik] == idxPicoPions[ip]) continue;

            //Get pion track information
            StPicoTrack const * pion = picoDst->track(idxPicoPions[ip]);

            //Make a D0 candidate
            StKaonPion kaonPion(kaon, pion, idxPicoKaons[ik], idxPicoPions[ip], pVtx, bField);

            //General check if the D0 candidate is good
            if (!isGoodPair(kaonPion)) continue;

            //Check if the D0 candidate is good for use in analysis
            if (isGoodMass(kaonPion)){

                //Save D0 candidate in EvenetMaker.root.picoD0.root
                mPicoD0Event->addKaonPion(&kaonPion);
            }

//-----------------------------------------------------------------
//The rest of the code is only for checking if there are any D0 candidates
//The other analysis of D0 candidates is done in StPicoD0AnaMaker.cxx
//Therefore any cuts and histograms AFTER THIS LINE are NOT used in the analysis
//-----------------------------------------------------------------

            //Check if the D0 candidate is good for use in analysis
            bool fillMass = isGoodQaPair(&kaonPion,*kaon,*pion);

            /*
            //Pepe the Frog appears if there is a good D0 candidate
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

            //Check if kaon and pion have opposite charge (=unlike)
            bool unlike = kaon->charge() * pion->charge() < 0 ? true : false;

            //Good candidates are saved in EvetMaker.root.picoD0.hists.root
            if(fillMass || unlike) mPicoD0Hists->addKaonPion(&kaonPion,fillMass, unlike);

        } // ... end of loop over all pion candidates

      } // ... end of loop over all pion candidates

   } //.. end of good event fill

   //Fill event information
   mPicoD0Event->addPicoEvent(*mPicoEvent);
   mPicoD0Hists->addEvent(*mPicoEvent,*mPicoD0Event,nHftTracks);

   // This should never be inside the good event block
   // because we want to save header information about all events, good or bad
   mTree->Fill();
   mPicoD0Event->clear("C");

   return kStOK;
}

bool StPicoD0EventMaker::isGoodEvent(int mYear){
   //Check if good run and good event

   return   isMinBiasTrigger(mYear) && //Check on trigger
            sqrt(mPicoEvent->primaryVertex().x()*mPicoEvent->primaryVertex().x()+mPicoEvent->primaryVertex().y()*mPicoEvent->primaryVertex().y()) < cuts::vr && //Check on vertex radius
            fabs(mPicoEvent->primaryVertex().z()) < cuts::vz && //Check on vertex z
            fabs(mPicoEvent->primaryVertex().z() - mPicoEvent->vzVpd()) < cuts::vzVpdVz; //Check on vertex z - vzVpd

   //It returns true if event is good
}

bool StPicoD0EventMaker::isMinBiasTrigger(int mYear){
    //Check if event is minBias trigger

    //List of minBias triggers is taken from StCuts.h
    const std::set<int>* mbTriggers = nullptr;

    //Different triggers for different years
    if(mYear ==2016) mbTriggers = &cuts::mbTriggers2016;
    if(mYear ==2014) mbTriggers = &cuts::mbTriggers2014;

    //Check if event has one of the minBias triggers
    StPicoEvent* event = static_cast<StPicoEvent*>(mPicoDstMaker->picoDst()->event());
    return std::any_of(mbTriggers->begin(), mbTriggers->end(), [&](int trigger) { return event->isTrigger(trigger); });

    //It returns true if event is minBias trigger
}

bool StPicoD0EventMaker::isGoodTrack(StPicoTrack const * const trk) const{
   // Require at least one hit on every layer of PXL and IST.
   // It is done here for tests on the preview II data.
   // The new StPicoTrack which is used in official production has a method to check this

   //Check if the track meets the HFT requirement
   bool HFTCondition = false;

   //Different HFT conditions for different years because of malfunction
   //2014 - Require at least one hit on every layer of PXL and IST
   if(mYear ==2016) HFTCondition = (trk->hasPxl1Hit() && trk->hasPxl2Hit()) && (trk->hasSstHit() || trk->hasIstHit());

   //2016 - Require at least one hit on every layer of PXL and (IST or SST)
   if(mYear ==2014) HFTCondition = trk->isHFTTrack();
   //Bool_t isHFTTrack()const { return hasPxl1Hit() && hasPxl2Hit() && (hasIstHit() || hasSstHit()); }

   //In StCuts.h is defined if the HFT is required and the nHitsFit value
   return (!cuts::requireHFT || HFTCondition) &&  trk->nHitsFit() >= cuts::nHitsFit;

   //Function returns true if track is good
}

bool StPicoD0EventMaker::isPion(StPicoTrack const * const trk) const{
    //Check if track is pion
    //Parameters are saved in StCuts.h
    return fabs(trk->nSigmaPion()) < cuts::nSigmaPion;

    //Function returns true if track is pion
}

bool StPicoD0EventMaker::isKaon(StPicoTrack const * const trk) const{
    //Check if track is kaon
    //Parameters are saved in StCuts.h
    return fabs(trk->nSigmaKaon()) < cuts::nSigmaKaon;

    //Function returns true if track is kaon
}

bool StPicoD0EventMaker::isGoodPair(StKaonPion const & kp) const{
    //Basic check for kaon-pion pair
    //Parameters are saved in StCuts.h
    //Those cuts are very general. The more strict cuts are applied in StPicoD0AnaMaker.cxx

    return  std::cos(kp.pointingAngle()) > cuts::cosTheta &&   //cosTheta cut
            kp.decayLength() > cuts::decayLength &&            //decayLength cut
            kp.dcaDaughters() < cuts::dcaDaughters;            //dcaDaughters cut

    //Function returns true if pair is good
}

bool StPicoD0EventMaker::isGoodMass(StKaonPion const & kp) const{
    //Check if mass is in the range
    //Parameters are saved in StCuts.h

    return kp.m() > cuts::minMass && kp.m() < cuts::maxMass;

    //Function returns true if mass is good
}

bool  StPicoD0EventMaker::isGoodQaPair(StKaonPion const& kp, StPicoTrack const& kaon,StPicoTrack const& pion){
    //Check if pair is good for QA
    //Parameters are saved in StCuts.h
    //Those cuts are ONLY for QA

    return pion.gPt() >= cuts::qaPt && kaon.gPt() >= cuts::qaPt &&
           pion.nHitsFit() >= cuts::qaNHitsFit && kaon.nHitsFit() >= cuts::qaNHitsFit &&
           fabs(kaon.nSigmaKaon()) < cuts::qaNSigmaKaon &&
           cos(kp.pointingAngle()) > cuts::qaCosTheta &&
           kp.pionDca() > cuts::qaPDca && kp.kaonDca() > cuts::qaKDca &&
           kp.dcaDaughters() < cuts::qaDcaDaughters;

    //Function returns true if pair is good for QA
}
