#include <limits>
#include <cmath>

#include "StHFQuadruplet.h"

#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "StarClassLibrary/SystemOfUnits.h"
#include "StPicoEvent/StPicoTrack.h"

ClassImp(StHFQuadruplet)

// _________________________________________________________
StHFQuadruplet::StHFQuadruplet(): mLorentzVector(StLorentzVectorF()), mDecayVertex(StThreeVectorF()),
  mPointingAngle(std::numeric_limits<float>::quiet_NaN()), mDecayLength(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Dca(std::numeric_limits<float>::quiet_NaN()), mParticle2Dca(std::numeric_limits<float>::quiet_NaN()), 
  mParticle3Dca(std::numeric_limits<float>::quiet_NaN()),  mParticle4Dca(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Idx(std::numeric_limits<unsigned short>::max()), mParticle2Idx(std::numeric_limits<unsigned short>::max()),
  mParticle3Idx(std::numeric_limits<unsigned short>::max()),mParticle4Idx(std::numeric_limits<unsigned short>::max()),
  mDcaDaughters12(std::numeric_limits<float>::max()),  mDcaDaughters13(std::numeric_limits<float>::max()),
  mDcaDaughters14(std::numeric_limits<float>::max()),  mDcaDaughters23(std::numeric_limits<float>::max()),
  mDcaDaughters24(std::numeric_limits<float>::max()),  mDcaDaughters34(std::numeric_limits<float>::max()),
  mCosThetaStar(std::numeric_limits<float>::quiet_NaN())
{
}

// _________________________________________________________
StHFQuadruplet::StHFQuadruplet(StHFQuadruplet const * t) : 
  mLorentzVector(t->mLorentzVector), mDecayVertex(t->mDecayVertex),
  mPointingAngle(t->mPointingAngle), mDecayLength(t->mDecayLength), 
  mParticle1Dca(t->mParticle1Dca), mParticle2Dca(t->mParticle2Dca), mParticle3Dca(t->mParticle3Dca),mParticle4Dca(t->mParticle4Dca),
  mParticle1Idx(t->mParticle1Idx), mParticle2Idx(t->mParticle2Idx), mParticle3Idx(t->mParticle3Idx),mParticle4Idx(t->mParticle4Idx),
  mDcaDaughters12(t->mDcaDaughters12), mDcaDaughters13(t->mDcaDaughters13), 
  mDcaDaughters14(t->mDcaDaughters14), mDcaDaughters23(t->mDcaDaughters23), 
  mDcaDaughters24(t->mDcaDaughters24), mDcaDaughters34(t->mDcaDaughters34), 
  mCosThetaStar(t->mCosThetaStar)
{
}
// _________________________________________________________
StHFQuadruplet::StHFQuadruplet(StPicoTrack const * const particle1, StPicoTrack const * const particle2, StPicoTrack const * const particle3, StPicoTrack const * particle4,
			       float p1MassHypo, float p2MassHypo, float p3MassHypo,float p4MassHypo,
			       unsigned short const p1Idx, unsigned short const p2Idx, unsigned short const p3Idx, unsigned short const p4Idx,
			 StThreeVectorF const & vtx, float const bField)  : 
  mLorentzVector(StLorentzVectorF()), mDecayVertex(StThreeVectorF()),
  mPointingAngle(std::numeric_limits<float>::quiet_NaN()), mDecayLength(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Dca(std::numeric_limits<float>::quiet_NaN()), mParticle2Dca(std::numeric_limits<float>::quiet_NaN()), 
  mParticle3Dca(std::numeric_limits<float>::quiet_NaN()),  mParticle4Dca(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Idx(p1Idx), mParticle2Idx(p2Idx),  mParticle3Idx(p3Idx),mParticle4Idx(p4Idx),
  mDcaDaughters12(std::numeric_limits<float>::max()), mDcaDaughters13(std::numeric_limits<float>::max()),  
  mDcaDaughters14(std::numeric_limits<float>::max()),  mDcaDaughters23(std::numeric_limits<float>::max()),
  mDcaDaughters24(std::numeric_limits<float>::max()),  mDcaDaughters34(std::numeric_limits<float>::max()),
  mCosThetaStar(std::numeric_limits<float>::min())
{
  // -- Create quadruplet out of 4 tracks
  //     prefixes code:
  //      p1 means particle 1
  //      p2 means particle 2
  //      p3 means particle 3
  //      p4 means particle 4
  //      pair means particle1-particle2 pair||  particle2-particle3 pair ||  particle1-particle3 pair
  //      quadruplet particle1-particle2-particle3-particle4

  if ((!particle1 || !particle2 || !particle3 || !particle4) || 
      (particle1->id() == particle2->id() || particle1->id() == particle3->id() || particle1->id() == particle4->id() || particle2->id() == particle3->id() || particle2->id() == particle4->id() || particle3->id() == particle4->id())) {
    mParticle1Idx = std::numeric_limits<unsigned short>::max();
    mParticle2Idx = std::numeric_limits<unsigned short>::max();
    mParticle3Idx = std::numeric_limits<unsigned short>::max();
    mParticle4Idx = std::numeric_limits<unsigned short>::max();
    return;
  }
  
  //StPhysicalHelixD p1Helix = particle1->dcaGeometry().helix();
  StPicoPhysicalHelix p1Helix = particle1->helix(bField);
  //StPhysicalHelixD p2Helix = particle2->dcaGeometry().helix();
  StPicoPhysicalHelix p2Helix = particle2->helix(bField);
  //StPhysicalHelixD p3Helix = particle3->dcaGeometry().helix();
  StPicoPhysicalHelix p3Helix = particle3->helix(bField);
  //StPhysicalHelixD p4Helix = particle4->dcaGeometry().helix();
  StPicoPhysicalHelix p4Helix = particle4->helix(bField);
   // -- move origins of helices to the primary vertex origin
  //p1Helix.moveOrigin(p1Helix.pathLength(vtx));
  p1Helix.moveOrigin(p1Helix.pathLength(vtx.x(),vtx.y()));
  //p2Helix.moveOrigin(p2Helix.pathLength(vtx));
  p2Helix.moveOrigin(p2Helix.pathLength(vtx.x(),vtx.y()));
  //p3Helix.moveOrigin(p3Helix.pathLength(vtx));
  p3Helix.moveOrigin(p3Helix.pathLength(vtx.x(),vtx.y()));
  //p4Helix.moveOrigin(p4Helix.pathLength(vtx));
  p4Helix.moveOrigin(p4Helix.pathLength(vtx.x(),vtx.y()));
  
  // -- use straight lines approximation to get point of DCA of particle1-particle2 pair
  //StThreeVectorF const p1Mom = p1Helix.momentum(bField * kilogauss);
  //StThreeVectorF const p2Mom = p2Helix.momentum(bField * kilogauss);
  //StThreeVectorF const p3Mom = p3Helix.momentum(bField * kilogauss);
  //StThreeVectorF const p4Mom = p4Helix.momentum(bField * kilogauss);
  
  StThreeVectorF const p1Mom = StThreeVectorF(  (p1Helix.momentum(bField * kilogauss)).x(),
                                                (p1Helix.momentum(bField * kilogauss)).y(),
                                                (p1Helix.momentum(bField * kilogauss)).z()
                                                                                  );
  StThreeVectorF const p2Mom = StThreeVectorF(  (p2Helix.momentum(bField * kilogauss)).x(),
                                                (p2Helix.momentum(bField * kilogauss)).y(),
                                                (p2Helix.momentum(bField * kilogauss)).z()
                                                                                  );
  StThreeVectorF const p3Mom = StThreeVectorF(  (p3Helix.momentum(bField * kilogauss)).x(),
                                                (p3Helix.momentum(bField * kilogauss)).y(),
                                                (p3Helix.momentum(bField * kilogauss)).z()
                                                                                  );
  StThreeVectorF const p4Mom = StThreeVectorF(  (p4Helix.momentum(bField * kilogauss)).x(),
                                                (p4Helix.momentum(bField * kilogauss)).y(),
                                                (p4Helix.momentum(bField * kilogauss)).z()
                                                                                  );
  
  //StPhysicalHelixD const p1StraightLine(p1Mom, p1Helix.origin(), 0, particle1->charge());
  //StPhysicalHelixD const p2StraightLine(p2Mom, p2Helix.origin(), 0, particle2->charge());
  //StPhysicalHelixD const p3StraightLine(p3Mom, p3Helix.origin(), 0, particle3->charge());
  //StPhysicalHelixD const p4StraightLine(p4Mom, p4Helix.origin(), 0, particle4->charge());
  StPhysicalHelixD const p1StraightLine(p1Mom, StThreeVectorF(p1Helix.origin().x(),p1Helix.origin().y(),p1Helix.origin().z()), 0, particle1->charge());
  StPhysicalHelixD const p2StraightLine(p2Mom, StThreeVectorF(p2Helix.origin().x(),p2Helix.origin().y(),p2Helix.origin().z()), 0, particle2->charge());
  StPhysicalHelixD const p3StraightLine(p3Mom, StThreeVectorF(p3Helix.origin().x(),p3Helix.origin().y(),p3Helix.origin().z()), 0, particle3->charge());
  StPhysicalHelixD const p4StraightLine(p4Mom, StThreeVectorF(p4Helix.origin().x(),p4Helix.origin().y(),p4Helix.origin().z()), 0, particle4->charge());


  pair<double, double> const ss12 = p1StraightLine.pathLengths(p2StraightLine);
  //pair<double, double> const ss12 = p1StraightLine.pathLengths(p2StraightLine.x(),p2StraightLine.y());
  StThreeVectorF const p1AtDcaToP2 = p1StraightLine.at(ss12.first);
  StThreeVectorF const p2AtDcaToP1 = p2StraightLine.at(ss12.second);

  pair<double, double> const ss13 = p1StraightLine.pathLengths(p3StraightLine);
  //pair<double, double> const ss13 = p1StraightLine.pathLengths(p3StraightLine.x(),p3StraightLine.y());
  StThreeVectorF const p1AtDcaToP3 = p1StraightLine.at(ss13.first);
  StThreeVectorF const p3AtDcaToP1 = p3StraightLine.at(ss13.second);

  pair<double, double> const ss14 = p1StraightLine.pathLengths(p4StraightLine);
  //pair<double, double> const ss14 = p1StraightLine.pathLengths(p4StraightLine.x(),p4StraightLine.y());
  StThreeVectorF const p1AtDcaToP4 = p1StraightLine.at(ss14.first);
  StThreeVectorF const p4AtDcaToP1 = p4StraightLine.at(ss14.second);

  pair<double, double> const ss23 = p2StraightLine.pathLengths(p3StraightLine);
  //pair<double, double> const ss23 = p2StraightLine.pathLengths(p3StraightLine.x(),p3StraightLine.y());
  StThreeVectorF const p2AtDcaToP3 = p2StraightLine.at(ss23.first);
  StThreeVectorF const p3AtDcaToP2 = p3StraightLine.at(ss23.second);

  pair<double, double> const ss24 = p2StraightLine.pathLengths(p4StraightLine);
  //pair<double, double> const ss24 = p2StraightLine.pathLengths(p4StraightLine.x(),p4StraightLine.y());
  StThreeVectorF const p2AtDcaToP4 = p2StraightLine.at(ss24.first);
  StThreeVectorF const p4AtDcaToP2 = p4StraightLine.at(ss24.second);
  
  pair<double, double> const ss34 = p3StraightLine.pathLengths(p4StraightLine);
  //pair<double, double> const ss34 = p3StraightLine.pathLengths(p4StraightLine.x(),p4StraightLine.y());
  StThreeVectorF const p3AtDcaToP4 = p3StraightLine.at(ss34.first);
  StThreeVectorF const p4AtDcaToP3 = p4StraightLine.at(ss34.second);
  
  // -- calculate DCA of particle1 to particl2 at their DCA
  mDcaDaughters12 = (p1AtDcaToP2 - p2AtDcaToP1).mag();

  // -- calculate DCA of particle1 to particl3 at their DCA
  mDcaDaughters13 = (p1AtDcaToP3 - p3AtDcaToP1).mag();
  
  // -- calculate DCA of particle1 to particl4 at their DCA
  mDcaDaughters14 = (p1AtDcaToP4 - p4AtDcaToP1).mag();

  // -- calculate DCA of particle2 to particl3 at their DCA
  mDcaDaughters23 = (p2AtDcaToP3 - p3AtDcaToP2).mag();

  // -- calculate DCA of particle2 to particl3 at their DCA
  mDcaDaughters24 = (p2AtDcaToP4 - p4AtDcaToP2).mag();
  
  // -- calculate DCA of particle2 to particl3 at their DCA
  mDcaDaughters34 = (p3AtDcaToP4 - p4AtDcaToP3).mag();

  
  // -- calculate decay vertex (secondary)
  StThreeVectorF mDecayVertex = ( p1AtDcaToP2 + p2AtDcaToP1 + p1AtDcaToP3 + p3AtDcaToP1 + p1AtDcaToP4 + p4AtDcaToP1 + p2AtDcaToP3 + p3AtDcaToP2 + p2AtDcaToP4 + p4AtDcaToP2 +p3AtDcaToP4 + p4AtDcaToP3)/12.0;
  
  // -- constructing mother daughter four momentum. Need helix (not straight line) for each daughter
 // double const p1AtV0 = p1Helix.pathLength( mDecayVertex );
  double const p1AtV0 = p1Helix.pathLength( mDecayVertex.x(),mDecayVertex.y() );
  //StThreeVectorF const p1MomAtDca = p1Helix.momentumAt(p1AtV0 ,  bField * kilogauss);
  StThreeVectorF const p1MomAtDca = StThreeVectorF( p1Helix.momentumAt(p1AtV0 ,  bField * kilogauss).x(),
                                                    p1Helix.momentumAt(p1AtV0 ,  bField * kilogauss).y(),
                                                    p1Helix.momentumAt(p1AtV0 ,  bField * kilogauss).z()
                                                                                                  );

  //double const p2AtV0 = p2Helix.pathLength( mDecayVertex );
  double const p2AtV0 = p2Helix.pathLength( mDecayVertex.x(),mDecayVertex.y()  );

//  StThreeVectorF const p2MomAtDca = p2Helix.momentumAt(p2AtV0 ,  bField * kilogauss);
    StThreeVectorF const p2MomAtDca = StThreeVectorF( p2Helix.momentumAt(p2AtV0 ,  bField * kilogauss).x(),
                                                    p2Helix.momentumAt(p2AtV0 ,  bField * kilogauss).y(),
                                                    p2Helix.momentumAt(p2AtV0 ,  bField * kilogauss).z()
                                                                                                  );


 // double const p3AtV0 = p3Helix.pathLength( mDecayVertex );
  double const p3AtV0 = p3Helix.pathLength( mDecayVertex.x(),mDecayVertex.y()  );

  //StThreeVectorF const p3MomAtDca = p3Helix.momentumAt(p3AtV0 ,  bField * kilogauss);
  StThreeVectorF const p3MomAtDca = StThreeVectorF( p3Helix.momentumAt(p3AtV0 ,  bField * kilogauss).x(),
                                                    p3Helix.momentumAt(p3AtV0 ,  bField * kilogauss).y(),
                                                    p3Helix.momentumAt(p3AtV0 ,  bField * kilogauss).z()
                                                                                                  );

  //double const p4AtV0 = p4Helix.pathLength( mDecayVertex );
  double const p4AtV0 = p4Helix.pathLength( mDecayVertex.x(),mDecayVertex.y()  );

  //StThreeVectorF const p4MomAtDca = p4Helix.momentumAt(p4AtV0 ,  bField * kilogauss);
    StThreeVectorF const p4MomAtDca = StThreeVectorF( p4Helix.momentumAt(p4AtV0 ,  bField * kilogauss).x(),
                                                    p4Helix.momentumAt(p4AtV0 ,  bField * kilogauss).y(),
                                                    p4Helix.momentumAt(p4AtV0 ,  bField * kilogauss).z()
                                                                                                  );


  StLorentzVectorF const p1FourMom(p1MomAtDca, p1MomAtDca.massHypothesis(p1MassHypo));
  StLorentzVectorF const p2FourMom(p2MomAtDca, p2MomAtDca.massHypothesis(p2MassHypo));
  StLorentzVectorF const p3FourMom(p3MomAtDca, p3MomAtDca.massHypothesis(p3MassHypo));
  StLorentzVectorF const p4FourMom(p4MomAtDca, p4MomAtDca.massHypothesis(p4MassHypo));
  
  mLorentzVector = p1FourMom + p2FourMom + p3FourMom + p4FourMom;
   
  // -- calculate cosThetaStar
  //    ->> Lomnitz: Need to rethink theta star
  StLorentzVectorF const quadrupletFourMomReverse(-mLorentzVector.px(), -mLorentzVector.py(), -mLorentzVector.pz(), mLorentzVector.e());
  StLorentzVectorF const p1FourMomStar = p1FourMom.boost(quadrupletFourMomReverse);
  mCosThetaStar = std::cos(p1FourMomStar.vect().angle(mLorentzVector.vect()));
  
  // -- calculate pointing angle and decay length
  StThreeVectorF const vtxToV0 = mDecayVertex - vtx;
  mPointingAngle = vtxToV0.angle(mLorentzVector.vect());
  mDecayLength = vtxToV0.mag();
  
  // --- calculate DCA of tracks to primary vertex
  //mParticle1Dca = (p1Helix.origin() - vtx).mag();
  //mParticle2Dca = (p2Helix.origin() - vtx).mag();
  //mParticle3Dca = (p3Helix.origin() - vtx).mag();
  //mParticle4Dca = (p4Helix.origin() - vtx).mag();
  mParticle1Dca = (p1Helix.origin() - TVector3(vtx.x(),vtx.y(),vtx.z())).Mag();
  mParticle2Dca = (p2Helix.origin() - TVector3(vtx.x(),vtx.y(),vtx.z())).Mag();
  mParticle3Dca = (p3Helix.origin() - TVector3(vtx.x(),vtx.y(),vtx.z())).Mag();
  mParticle4Dca = (p4Helix.origin() - TVector3(vtx.x(),vtx.y(),vtx.z())).Mag();
}
// _________________________________________________________
StHFQuadruplet::StHFQuadruplet(StPicoTrack const * const particle1, StPicoTrack const * const particle2, StPicoTrack const * const particle3, StHFPair const * particle4,
			       float p1MassHypo, float p2MassHypo, float p3MassHypo,float p4MassHypo,
			       unsigned short const p1Idx, unsigned short const p2Idx, unsigned short const p3Idx, unsigned short const p4Idx,
			 StThreeVectorF const & vtx, float const bField)  : 
  mLorentzVector(StLorentzVectorF()), mDecayVertex(StThreeVectorF()),
  mPointingAngle(std::numeric_limits<float>::quiet_NaN()), mDecayLength(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Dca(std::numeric_limits<float>::quiet_NaN()), mParticle2Dca(std::numeric_limits<float>::quiet_NaN()), 
  mParticle3Dca(std::numeric_limits<float>::quiet_NaN()),  mParticle4Dca(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Idx(p1Idx), mParticle2Idx(p2Idx),  mParticle3Idx(p3Idx),mParticle4Idx(p4Idx),
  mDcaDaughters12(std::numeric_limits<float>::max()), mDcaDaughters13(std::numeric_limits<float>::max()),  
  mDcaDaughters14(std::numeric_limits<float>::max()),  mDcaDaughters23(std::numeric_limits<float>::max()),
  mDcaDaughters24(std::numeric_limits<float>::max()),  mDcaDaughters34(std::numeric_limits<float>::max()),
  mCosThetaStar(std::numeric_limits<float>::min())
{
  // -- Create quadruplet out of 4 tracks
  //     prefixes code:
  //      p1 means particle 1
  //      p2 means particle 2
  //      p3 means particle 3
  //      p4 means particle 4
  //      pair means particle1-particle2 pair||  particle2-particle3 pair ||  particle1-particle3 pair
  //      triplet particle1-particle2-particle3

  if ((!particle1 || !particle2 || !particle3 || !particle4) || 
      (particle1->id() == particle2->id() || particle1->id() == particle3->id() || particle1->id() == (particle4->particle1Idx() || particle4->particle2Idx()) || particle2->id() == particle3->id() || particle2->id() == (particle4->particle1Idx() || particle4->particle2Idx()) || particle3->id() == (particle4->particle1Idx() || particle4->particle2Idx()))) {
    mParticle1Idx = std::numeric_limits<unsigned short>::max();
    mParticle2Idx = std::numeric_limits<unsigned short>::max();
    mParticle3Idx = std::numeric_limits<unsigned short>::max();
    mParticle4Idx = std::numeric_limits<unsigned short>::max();
    return;
  }
  
  //StPhysicalHelixD p1Helix = particle1->dcaGeometry().helix();
  //StPhysicalHelixD p2Helix = particle2->dcaGeometry().helix();
  //StPhysicalHelixD p3Helix = particle3->dcaGeometry().helix();
  StPicoPhysicalHelix p1Helix = particle1->helix(bField);
  StPicoPhysicalHelix p2Helix = particle2->helix(bField);
  StPicoPhysicalHelix p3Helix = particle3->helix(bField);
  
   // -- move origins of helices to the primary vertex origin
  //p1Helix.moveOrigin(p1Helix.pathLength(vtx));
  //p2Helix.moveOrigin(p2Helix.pathLength(vtx));
  //p3Helix.moveOrigin(p3Helix.pathLength(vtx));
  p1Helix.moveOrigin(p1Helix.pathLength(vtx.x(),vtx.y()));
  p2Helix.moveOrigin(p2Helix.pathLength(vtx.x(),vtx.y()));
  p3Helix.moveOrigin(p3Helix.pathLength(vtx.x(),vtx.y()));
  // No need to move p4 to origin already defined as such

  // -- use straight lines approximation to get point of DCA of particle1-particle2 pair
  //StThreeVectorF const p1Mom = p1Helix.momentum(bField * kilogauss);
  //StThreeVectorF const p2Mom = p2Helix.momentum(bField * kilogauss);
  //StThreeVectorF const p3Mom = p3Helix.momentum(bField * kilogauss);
  StThreeVectorF const p1Mom = StThreeVectorF(p1Helix.momentum(bField * kilogauss).x(),
                                              p1Helix.momentum(bField * kilogauss).y(),
                                              p1Helix.momentum(bField * kilogauss).z());
  StThreeVectorF const p2Mom = StThreeVectorF(p2Helix.momentum(bField * kilogauss).x(),
                                              p2Helix.momentum(bField * kilogauss).y(),
                                              p2Helix.momentum(bField * kilogauss).z());
  StThreeVectorF const p3Mom = StThreeVectorF(p3Helix.momentum(bField * kilogauss).x(),
                                              p3Helix.momentum(bField * kilogauss).y(),
                                              p3Helix.momentum(bField * kilogauss).z());
  StThreeVectorF const p4Mom(particle4->px(), particle4->py(), particle4->pz());
  
  // Build p4 helix from pair: assuming pair will always be neutral charge and origin is set to the given vtx 
  StPhysicalHelixD p4Helix(p4Mom, vtx, bField * kilogauss, 0);

  //StPhysicalHelixD const p1StraightLine(p1Mom, p1Helix.origin(), 0, particle1->charge());
  //StPhysicalHelixD const p2StraightLine(p2Mom, p2Helix.origin(), 0, particle2->charge());
  //StPhysicalHelixD const p3StraightLine(p3Mom, p3Helix.origin(), 0, particle3->charge());

  StPhysicalHelixD const p1StraightLine(p1Mom, StThreeVectorF(p1Helix.origin().x(),p1Helix.origin().y(),p1Helix.origin().z()), 0, particle1->charge());
  StPhysicalHelixD const p2StraightLine(p2Mom, StThreeVectorF(p2Helix.origin().x(),p2Helix.origin().y(),p2Helix.origin().z()), 0, particle2->charge());
  StPhysicalHelixD const p3StraightLine(p3Mom, StThreeVectorF(p3Helix.origin().x(),p3Helix.origin().y(),p3Helix.origin().z()), 0, particle3->charge());
  
  // Four paticle, i.e. pair, is already straight line (charge neutral)
  
  pair<double, double> const ss12 = p1StraightLine.pathLengths(p2StraightLine);
  StThreeVectorF const p1AtDcaToP2 = p1StraightLine.at(ss12.first);
  StThreeVectorF const p2AtDcaToP1 = p2StraightLine.at(ss12.second);

  pair<double, double> const ss13 = p1StraightLine.pathLengths(p3StraightLine);
  StThreeVectorF const p1AtDcaToP3 = p1StraightLine.at(ss13.first);
  StThreeVectorF const p3AtDcaToP1 = p3StraightLine.at(ss13.second);

  pair<double, double> const ss14 = p1StraightLine.pathLengths(p4Helix);
  StThreeVectorF const p1AtDcaToP4 = p1StraightLine.at(ss14.first);
  StThreeVectorF const p4AtDcaToP1 = p4Helix.at(ss14.second);

  pair<double, double> const ss23 = p2StraightLine.pathLengths(p3StraightLine);
  StThreeVectorF const p2AtDcaToP3 = p2StraightLine.at(ss23.first);
  StThreeVectorF const p3AtDcaToP2 = p3StraightLine.at(ss23.second);

  pair<double, double> const ss24 = p2StraightLine.pathLengths(p4Helix);
  StThreeVectorF const p2AtDcaToP4 = p2StraightLine.at(ss24.first);
  StThreeVectorF const p4AtDcaToP2 = p4Helix.at(ss24.second);
  
  pair<double, double> const ss34 = p3StraightLine.pathLengths(p4Helix);
  StThreeVectorF const p3AtDcaToP4 = p3StraightLine.at(ss34.first);
  StThreeVectorF const p4AtDcaToP3 = p4Helix.at(ss34.second);
  
  // -- calculate DCA of particle1 to particl2 at their DCA
  mDcaDaughters12 = (p1AtDcaToP2 - p2AtDcaToP1).mag();

  // -- calculate DCA of particle1 to particl3 at their DCA
  mDcaDaughters13 = (p1AtDcaToP3 - p3AtDcaToP1).mag();
  
  // -- calculate DCA of particle1 to particl4 at their DCA
  mDcaDaughters14 = (p1AtDcaToP4 - p4AtDcaToP1).mag();

  // -- calculate DCA of particle2 to particl3 at their DCA
  mDcaDaughters23 = (p2AtDcaToP3 - p3AtDcaToP2).mag();

  // -- calculate DCA of particle2 to particl3 at their DCA
  mDcaDaughters24 = (p2AtDcaToP4 - p4AtDcaToP2).mag();
  
  // -- calculate DCA of particle2 to particl3 at their DCA
  mDcaDaughters34 = (p3AtDcaToP4 - p4AtDcaToP3).mag();

  
  // -- calculate decay vertex (secondary)
  StThreeVectorF mDecayVertex = ( p1AtDcaToP2 + p2AtDcaToP1 + p1AtDcaToP3 + p3AtDcaToP1 + p1AtDcaToP4 + p4AtDcaToP1 + p2AtDcaToP3 + p3AtDcaToP2 + p2AtDcaToP4 + p4AtDcaToP2 +p3AtDcaToP4 + p4AtDcaToP3)/12.0;
  
  // -- constructing mother daughter four momentum. Need helix (not straight line) for each daughter
  //double const p1AtV0 = p1Helix.pathLength( mDecayVertex );
  double const p1AtV0 = p1Helix.pathLength( mDecayVertex.x(),mDecayVertex.y() );
  //StThreeVectorF const p1MomAtDca = p1Helix.momentumAt(p1AtV0 ,  bField * kilogauss);
  StThreeVectorF const p1MomAtDca = StThreeVectorF( p1Helix.momentumAt(p1AtV0 ,  bField * kilogauss).x(),
                                                    p1Helix.momentumAt(p1AtV0 ,  bField * kilogauss).y(),
                                                    p1Helix.momentumAt(p1AtV0 ,  bField * kilogauss).z()
                                                                                                  );

  //double const p2AtV0 = p2Helix.pathLength( mDecayVertex );
  double const p2AtV0 = p2Helix.pathLength( mDecayVertex.x(),mDecayVertex.y() );
  //StThreeVectorF const p2MomAtDca = p2Helix.momentumAt(p2AtV0 ,  bField * kilogauss);
  StThreeVectorF const p2MomAtDca = StThreeVectorF( p2Helix.momentumAt(p2AtV0 ,  bField * kilogauss).x(),
                                                    p2Helix.momentumAt(p2AtV0 ,  bField * kilogauss).y(),
                                                    p2Helix.momentumAt(p2AtV0 ,  bField * kilogauss).z()
                                                                                                  );
  
  //double const p3AtV0 = p3Helix.pathLength( mDecayVertex );
  double const p3AtV0 = p3Helix.pathLength( mDecayVertex.x(),mDecayVertex.y() );
  //StThreeVectorF const p3MomAtDca = p3Helix.momentumAt(p3AtV0 ,  bField * kilogauss);
    StThreeVectorF const p3MomAtDca = StThreeVectorF( p3Helix.momentumAt(p3AtV0 ,  bField * kilogauss).x(),
                                                    p3Helix.momentumAt(p3AtV0 ,  bField * kilogauss).y(),
                                                    p3Helix.momentumAt(p3AtV0 ,  bField * kilogauss).z()
                                                                                                  );

  //double const p4AtV0 = p4Helix.pathLength( mDecayVertex );
  double const p4AtV0 = p4Helix.pathLength( mDecayVertex.x(),mDecayVertex.y() );

  //StThreeVectorF const p4MomAtDca = p4Helix.momentumAt(p4AtV0 ,  bField * kilogauss);
  StThreeVectorF const p4MomAtDca = StThreeVectorF( p4Helix.momentumAt(p4AtV0 ,  bField * kilogauss).x(),
                                                    p4Helix.momentumAt(p4AtV0 ,  bField * kilogauss).y(),
                                                    p4Helix.momentumAt(p4AtV0 ,  bField * kilogauss).z()
                                                                                                  );

  StLorentzVectorF const p1FourMom(p1MomAtDca, p1MomAtDca.massHypothesis(p1MassHypo));
  StLorentzVectorF const p2FourMom(p2MomAtDca, p2MomAtDca.massHypothesis(p2MassHypo));
  StLorentzVectorF const p3FourMom(p3MomAtDca, p3MomAtDca.massHypothesis(p3MassHypo));
  StLorentzVectorF const p4FourMom = particle4->lorentzVector();
  
  mLorentzVector = p1FourMom + p2FourMom + p3FourMom + p4FourMom;
   
  // -- calculate cosThetaStar
  //    ->> Lomnitz: Need to rethink theta star
  StLorentzVectorF const quadrupletFourMomReverse(-mLorentzVector.px(), -mLorentzVector.py(), -mLorentzVector.pz(), mLorentzVector.e());
  StLorentzVectorF const p1FourMomStar = p1FourMom.boost(quadrupletFourMomReverse);
  mCosThetaStar = std::cos(p1FourMomStar.vect().angle(mLorentzVector.vect()));
  
  // -- calculate pointing angle and decay length
  StThreeVectorF const vtxToV0 = mDecayVertex - vtx;
  mPointingAngle = vtxToV0.angle(mLorentzVector.vect());
  mDecayLength = vtxToV0.mag();
  
  // --- calculate DCA of tracks to primary vertex
  //mParticle1Dca = (p1Helix.origin() - vtx).mag();
  //mParticle2Dca = (p2Helix.origin() - vtx).mag();
  //mParticle3Dca = (p3Helix.origin() - vtx).mag();
  //mParticle4Dca = (p4Helix.origin() - vtx).mag();

  mParticle1Dca = (StThreeVectorF(p1Helix.origin().x(),p1Helix.origin().y(),p1Helix.origin().z()) - vtx).mag();
  mParticle2Dca = (StThreeVectorF(p2Helix.origin().x(),p2Helix.origin().y(),p2Helix.origin().z()) - vtx).mag();
  mParticle3Dca = (StThreeVectorF(p3Helix.origin().x(),p3Helix.origin().y(),p3Helix.origin().z()) - vtx).mag();
  mParticle4Dca = (StThreeVectorF(p4Helix.origin().x(),p4Helix.origin().y(),p4Helix.origin().z()) - vtx).mag();

  //mParticle1Dca = (p1Helix.origin() - TVector3(vtx.x(),vtx.y(),vtx.z())).Mag();
  //mParticle2Dca = (p2Helix.origin() - TVector3(vtx.x(),vtx.y(),vtx.z())).Mag();
  //mParticle3Dca = (p3Helix.origin() - TVector3(vtx.x(),vtx.y(),vtx.z())).Mag();
  //mParticle4Dca = (p4Helix.origin() - TVector3(vtx.x(),vtx.y(),vtx.z())).Mag();
}

