#ifndef CUTS_H
#define CUTS_H

/* **************************************************
 *  Cuts namespace.
 *
 *  Authors:  Xin Dong        (xdong@lbl.gov)
 *            Michael Lomnitz (mrlomnitz@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *            Jochen Thaeder  (jmthader@lbl.gov)   
 *
 * **************************************************
 */

#include "Rtypes.h"
#include <string>

namespace mycuts
{
   extern std::string const prescalesFilesDirectoryName;
   //event
   extern UShort_t const triggerWord;
   extern float const vz;
   extern float const vzVpdVz;

   //tracking
   extern int const nHitsFit;
   extern bool const requireHFT;
   extern float const minPt;

   //pions
   extern float const nSigmaPion;

   //kaons
   extern float const nSigmaKaon;
   extern float const kTofBetaDiff;

   // tree kaonPion pair cuts
   extern float const cosTheta;
   extern float const dcaDaughters;
   extern float const decayLength;
   extern float const minMass;
   extern float const maxMass;
   extern float const kDca;
   extern float const pDca;

   // histograms kaonPion pair cuts
   extern float const qaNHitsFit;
   extern float const qaNSigmaKaon;
   extern float const qaCosTheta;
   extern float const qaDcaDaughters;
   extern float const qaKDca;
   extern float const qaPDca;
   //Hadron cuts
   extern float const hadronPtMin;
   extern float const hadronPtMax;
   extern float const corDetaMin;
   extern float const corDetaMax;

   extern float const pionDCA_cut[6][5];
   extern float const kaonDCA_cut[6][5];
   extern float const DCA_D0_cut[6][5];
   extern float const D0_decayLength_cut[6][5];
   extern float const pionkaonDCA_cut[6][5];
}
#endif
