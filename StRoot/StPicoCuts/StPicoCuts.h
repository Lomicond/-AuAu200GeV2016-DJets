#ifndef STPICOCUTSBASE_H
#define STPICOCUTSBASE_H

#include <limits>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

#include "TVector3.h"
#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "StarClassLibrary/SystemOfUnits.h"
#include "StarClassLibrary/StLorentzVectorF.hh"
#include "StBTofUtil/tofPathLength.hh"
#include "StBTofUtil/StV0TofCorrection.h"
#include "phys_constants.h"

#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoBEmcPidTraits.h"
#include "StPicoEvent/StPicoBTowHit.h"

#include "TNamed.h"
#include "TString.h"

class StPicoCuts : public TNamed
{
 public:
  
  StPicoCuts();
  StPicoCuts(const Char_t *name);
  ~StPicoCuts();
  
  void initBase();
  virtual void init() { initBase(); }


  bool isGoodEvent(StPicoDst const * const picoDst, int *aEventCuts = NULL);
  bool isGoodRun(StPicoEvent const * const picoEvent) const;
  //bool isGoodRunHFT(StPicoEvent const * const picoEvent) const;
	
	
  bool isGoodTrigger(StPicoEvent const * const picoEvent) const;
  bool isGoodTrack(StPicoTrack const * const trk) const;
  bool isGoodPrimaryTrack(StPicoTrack const * const trk) const;
	bool isGoodTowHit(StPicoBTowHit const * const towHit) const;

	//bool isGoodBTow(StPicoTrack const * const picoTrack) const;

	//bool isGoodStripEta(StPicoTrack const * const picoTrack) const;
	//bool isGoodStripPhi(StPicoTrack const * const picoTrack) const;
  
  const unsigned int&  eventStatMax()  const { return mEventStatMax; }

  void setBadRunListFileName(const char* fileName);
  void addTriggerId(unsigned int triggerId);

	//void setHotTowerListFileName(const char* fileName);

	//void setHotStripEtaListFileName(const char* fileName);
	//void setHotStripPhiListFileName(const char* fileName);

  void setCutVzMax(float f);
  void setCutVzVpdVzMax(float f);

  void setCutNHitsFitMin(int i);
  void setCutRequireHFT(bool b);
  void setCutNHitsFitnHitsMax(float f);

  void setCutPtRange(float min, float max);
  void setCutERange(float min, float max);

  void setCutRefMult(float min, float max);
  void setCutDcaMin(float min);
  void setCutTPCNSigma(float f);
  void setCutEta(float f);

  float getCutEta();

private:
  
  StPicoCuts(StPicoCuts const &);       
  StPicoCuts& operator=(StPicoCuts const &);


  TVector3 mPrimVtx;   // primary vertex of current event
  const StPicoDst*  mPicoDst;   //! ptr to picoDst

  unsigned int mEventStatMax;   // number of event cuts


  // -- bad run list
  TString mBadRunListFileName;
  std::vector<int> mVecBadRunList;

  // -- trigger id list
  std::vector<unsigned int> mVecTriggerIdList;

	// -- hot tower list

/*	TString mHotTowerListFileName;
  std::vector<int> mVecHotTowerList;*/


  // -- event cuts
  float mVzMax;
  float mVzVpdVzMax;
  float mRefMultMin;
  float mRefMultMax;

	// -- tracking
  int   mNHitsFitMin;
  bool  mRequireHFT; 
  float mNHitsFitnHitsMax;
  float mEta;

  // -- acceptance
  float mPtRange[2];
  float mERangeMin, mERangeMax;

  // -- dca to primary vertex
  float mDcaMin;

  float mTPCNSigmaMax;

  ClassDef(StPicoCuts,1)
};

inline void StPicoCuts::setBadRunListFileName(const char* fileName) { mBadRunListFileName = fileName; }
inline void StPicoCuts::addTriggerId(unsigned int triggerId) {mVecTriggerIdList.push_back(triggerId);}

//inline void StPicoCuts::setHotTowerListFileName(const char* fileName) { mHotTowerListFileName = fileName; }

//inline void StPicoCuts::setHotStripEtaListFileName(const char* fileName) { mHotStripEtaListFileName = fileName; }
//inline void StPicoCuts::setHotStripPhiListFileName(const char* fileName) { mHotStripPhiListFileName = fileName; }

inline void StPicoCuts::setCutVzMax(float f)              { mVzMax            = f; }
inline void StPicoCuts::setCutVzVpdVzMax(float f)         { mVzVpdVzMax       = f; }

inline void StPicoCuts::setCutNHitsFitMin(int i)          { mNHitsFitMin      = i; }
inline void StPicoCuts::setCutRequireHFT(bool b)          { mRequireHFT       = b; }
inline void StPicoCuts::setCutNHitsFitnHitsMax(float f)   { mNHitsFitnHitsMax = f; }


inline void StPicoCuts::setCutPtRange(float min, float max)            { mPtRange[0] = min; mPtRange[1] = max; }
inline void StPicoCuts::setCutERange(float min, float max)            { mERangeMin = min; mERangeMax = max; }

inline void StPicoCuts::setCutDcaMin(float min)                        { mDcaMin = min; }
inline void StPicoCuts::setCutTPCNSigma(float f)                       { mTPCNSigmaMax = f; }
inline void StPicoCuts::setCutEta(float f)                           { mEta = f; }
inline void StPicoCuts::setCutRefMult(float min, float max)          { mRefMultMin = min; mRefMultMax = max; }

inline float StPicoCuts::getCutEta() {return mEta; }

#endif
