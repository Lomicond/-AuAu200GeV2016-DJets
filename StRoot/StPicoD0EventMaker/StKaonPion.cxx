#include <limits>
#include <cmath>

#ifdef __ROOT__
#include "StKaonPion.h"

#include "StLorentzVectorF.hh"
#include "StThreeVectorF.hh"
#include "StPhysicalHelixD.hh"
#include "phys_constants.h"
#include "SystemOfUnits.h"
#include "StPicoEvent/StPicoTrack.h"


ClassImp(StKaonPion)


StKaonPion::StKaonPion(): mLorentzVector(),
   mPointingAngle(std::numeric_limits<float>::quiet_NaN()), mDecayLength(std::numeric_limits<float>::quiet_NaN()),
   mKaonDca(std::numeric_limits<float>::quiet_NaN()), mPionDca(std::numeric_limits<float>::quiet_NaN()),
   mKaonIdx(std::numeric_limits<unsigned short>::quiet_NaN()), mPionIdx(std::numeric_limits<unsigned short>::quiet_NaN()),
   mDcaDaughters(std::numeric_limits<float>::quiet_NaN()), mCosThetaStar(std::numeric_limits<float>::quiet_NaN())

{
}
//------------------------------------
StKaonPion::StKaonPion(StKaonPion const * t) : mLorentzVector(t->mLorentzVector),
   mPointingAngle(t->mPointingAngle), mDecayLength(t->mDecayLength),
   mKaonDca(t->mKaonDca), mPionDca(t->mPionDca),
   mKaonIdx(t->mKaonIdx), mPionIdx(t->mPionIdx),
   mDcaDaughters(t->mDcaDaughters), mCosThetaStar(t->mCosThetaStar)
{
}
//------------------------------------
StKaonPion::StKaonPion(StPicoTrack const * const kaon, StPicoTrack const * const pion,
                       unsigned short const kIdx, unsigned short const pIdx,
                       StThreeVectorF const & vtx, float const bField) : mLorentzVector(),
   mPointingAngle(std::numeric_limits<float>::quiet_NaN()), mDecayLength(std::numeric_limits<float>::quiet_NaN()),
   mKaonDca(std::numeric_limits<float>::quiet_NaN()), mPionDca(std::numeric_limits<float>::quiet_NaN()),
   mKaonIdx(kIdx), mPionIdx(pIdx),
   mDcaDaughters(std::numeric_limits<float>::quiet_NaN()), mCosThetaStar(std::numeric_limits<float>::quiet_NaN())
{
   if ((!kaon || !pion) || (kaon->id() == pion->id()))
   {
      mKaonIdx = std::numeric_limits<unsigned short>::quiet_NaN();
      mPionIdx = std::numeric_limits<unsigned short>::quiet_NaN();
      return;
   }

   /// prefixes code:
   ///   k means kaon
   ///   p means pion
   ///   kp means kaon-pion pair


   // to be used for testing with preview II pico production
   //StPhysicalHelixD kHelix = kaon->dcaGeometry().helix();
   StPicoPhysicalHelix kHelix = kaon->helix(bField);
   //StPhysicalHelixD pHelix = pion->dcaGeometry().helix();
   StPicoPhysicalHelix pHelix = pion->helix(bField);

   // move origins of helices to the primary vertex origin
   //kHelix.moveOrigin(kHelix.pathLength(vtx));
   kHelix.moveOrigin(kHelix.pathLength(TVector3(vtx.x(),vtx.y(),vtx.z())));
   //pHelix.moveOrigin(pHelix.pathLength(vtx));
   pHelix.moveOrigin(pHelix.pathLength(TVector3(vtx.x(),vtx.y(),vtx.z())));

   // use straight lines approximation to get point of DCA of kaon-pion pair
   //StThreeVectorF const kMom = kHelix.momentum(bField * kilogauss);
   StThreeVectorF const kMom = StThreeVectorF(kHelix.momentum(bField * kilogauss).x(),kHelix.momentum(bField * kilogauss).y(),kHelix.momentum(bField * kilogauss).z());
   //StThreeVectorF const pMom = pHelix.momentum(bField * kilogauss);
   StThreeVectorF const pMom = StThreeVectorF(pHelix.momentum(bField * kilogauss).x(),pHelix.momentum(bField * kilogauss).y(),pHelix.momentum(bField * kilogauss).z());
   //StPhysicalHelixD const kStraightLine(kMom, kHelix.origin(), 0, kaon->charge());
   StPhysicalHelixD const kStraightLine(kMom, StThreeVectorF(kHelix.origin().x(),kHelix.origin().y(),kHelix.origin().z()), 0, kaon->charge());

   //StPhysicalHelixD const pStraightLine(pMom, pHelix.origin(), 0, pion->charge());
   StPhysicalHelixD const pStraightLine(pMom, StThreeVectorF(pHelix.origin().x(),pHelix.origin().y(),pHelix.origin().z()), 0, pion->charge());

   pair<double, double> const ss = kStraightLine.pathLengths(pStraightLine);
   StThreeVectorF const kAtDcaToPion = kStraightLine.at(ss.first);
   StThreeVectorF const pAtDcaToKaon = pStraightLine.at(ss.second);

   // calculate DCA of pion to kaon at their DCA
   mDcaDaughters = (kAtDcaToPion - pAtDcaToKaon).mag();

   // calculate Lorentz vector of kaon-pion pair
   //StThreeVectorF const kMomAtDca = kHelix.momentumAt(ss.first, bField * kilogauss);
   StThreeVectorF const kMomAtDca = StThreeVectorF((kHelix.momentumAt(ss.first, bField * kilogauss)).x(),(kHelix.momentumAt(ss.first, bField * kilogauss)).y(),(kHelix.momentumAt(ss.first, bField * kilogauss)).z());
   //StThreeVectorF const pMomAtDca = pHelix.momentumAt(ss.second, bField * kilogauss);
   StThreeVectorF const pMomAtDca = StThreeVectorF((pHelix.momentumAt(ss.second, bField * kilogauss)).x(),(pHelix.momentumAt(ss.second, bField * kilogauss)).y(),(pHelix.momentumAt(ss.second, bField * kilogauss)).z());

   StLorentzVectorF const kFourMom(kMomAtDca, kMomAtDca.massHypothesis(M_KAON_PLUS));
   StLorentzVectorF const pFourMom(pMomAtDca, pMomAtDca.massHypothesis(M_PION_PLUS));

   mLorentzVector = kFourMom + pFourMom;

   // calculate cosThetaStar
   StLorentzVectorF const kpFourMomReverse(-mLorentzVector.px(), -mLorentzVector.py(), -mLorentzVector.pz(), mLorentzVector.e());
   StLorentzVectorF const kFourMomStar = kFourMom.boost(kpFourMomReverse);
   mCosThetaStar = std::cos(kFourMomStar.vect().angle(mLorentzVector.vect()));

   // calculate pointing angle and decay length
   StThreeVectorF const vtxToV0 = (kAtDcaToPion + pAtDcaToKaon) * 0.5 - vtx;
   mPointingAngle = vtxToV0.angle(mLorentzVector.vect());
   mDecayLength = vtxToV0.mag();
  
     // calculate DCA of tracks to primary vertex
   // mKaonDca = (kHelix.origin() - vtx).mag();
   mKaonDca = (StThreeVectorF(kHelix.origin().x(),kHelix.origin().y(),kHelix.origin().z()) - vtx).mag();
   // mPionDca = (pHelix.origin() - vtx).mag();
   mPionDca = (StThreeVectorF(pHelix.origin().x(),pHelix.origin().y(),pHelix.origin().z()) - vtx).mag();
}
#endif // __ROOT__
