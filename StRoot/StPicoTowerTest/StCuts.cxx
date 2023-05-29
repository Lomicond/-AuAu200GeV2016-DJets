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

#include "StCuts.h"

namespace mycuts
{
   // path to lists of triggers prescales
   // lists are obtained from http://www.star.bnl.gov/protected/common/common2014/trigger2014/plots_au200gev/
   std::string const prescalesFilesDirectoryName = "./run14AuAu200GeVPrescales";
   //event
   UShort_t const triggerWord = 0x1F; //first five bits see http://rnc.lbl.gov/~xdong/SoftHadron/picoDst.html
   float const vz = 6.0;// cm.
   float const vzVpdVz = 3.0; // 3 cm.

   float const vr = 2.0; // cm.

   //tracking
   int const nHitsFit = 20;
   bool const requireHFT = true;
   float const minPt = 0.6;

   //pions
   float const nSigmaPion = 3.0;

   //kaons
   float const nSigmaKaon = 2.0;
   float const kTofBetaDiff = 0.03;

   // tree kaonPion pair cuts
   float const cosTheta = 0.995; // minimum
   float const dcaDaughters = 0.0050; // maximum
   float const decayLength = 0.0030; // minimum
   float const minMass = 1.6;
   float const maxMass = 2.1;

   // histograms kaonPion pair cuts
   float const qaNHitsFit = 20;
   float const qaNSigmaKaon = 2.0;
   float const kDca = 0.008; // minimum
   float const pDca = 0.008;
   // Hadron cuts
   float const hadronPtMin= 0.2;
   float const hadronPtMax = 2.0;
   float const corDetaMin = 0.75;
   float const corDetaMax = 10;

   //Pt-Centrality-Cuts
   //Centality                          0-10%       10-20%         20-40%         40-60%        60-80%
float const pionDCA_cut[6][5] = {       {0.0131,     0.0141,        0.0131,      0.0145,       0.0098},      //pT = 0-0.5
                                        {0.0105,     0.0100,        0.0113,      0.0128,       0.0098},      //pT = 0.5-1
                                        {0.0093,     0.0074,        0.0099,      0.0072,       0.0083},      //pT = 1-2
                                        {0.0097,     0.0077,        0.0106,      0.0079,       0.0073},      //pT = 2-3
                                        {0.0067,     0.0066,        0.0065,      0.0060,       0.0056},      //pT = 3-5
                                        {0.0055,     0.0052,        0.0052,      0.0051,       0.0050}  };   //pT = 5-10                                                              }   //pT = 0.5-1

   //Centality                          0-10%       10-20%         20-40%         40-60%        60-80%
float const kaonDCA_cut[6][5] = {       {0.0138,     0.0145,        0.0151,      0.0140,       0.0106},      //pT = 0-0.5
                                        {0.0109,     0.0113,        0.0102,      0.0100,       0.0106},      //pT = 0.5-1
                                        {0.0082,     0.0094,        0.0104,      0.0075,       0.0069},      //pT = 1-2
                                        {0.0094,     0.0089,        0.0099,      0.0072,       0.0068},      //pT = 2-3
                                        {0.0076,     0.0069,        0.0063,      0.0060,       0.0050},      //pT = 3-5
                                        {0.0054,     0.0050,        0.0050,      0.0050,       0.0050}  };   //pT = 5-10       

   //Centality                          0-10%       10-20%         20-40%         40-60%        60-80%
float const DCA_D0_cut[6][5] =     {    {0.0062,     0.0063,        0.0066,      0.0072,       0.0076},      //pT = 0-0.5
                                        {0.0055,     0.0047,        0.0055,      0.0057,       0.0076},      //pT = 0.5-1
                                        {0.0040,     0.0045,        0.0053,      0.0058,       0.0053},      //pT = 1-2
                                        {0.0040,     0.0046,        0.0046,      0.0049,       0.0054},      //pT = 2-3
                                        {0.0040,     0.0042,        0.0041,      0.0049,       0.0054},      //pT = 3-5
                                        {0.0044,     0.0044,        0.0050,      0.0047,       0.0042}  };   //pT = 5-10   

   //Centality                          0-10%       10-20%         20-40%         40-60%        60-80%
float const D0_decayLength_cut[6][5] ={ {0.0100,     0.0172,        0.0178,      0.0171,       0.0175},      //pT = 0-0.5
                                        {0.0199,     0.0215,        0.0206,      0.0196,       0.0175},      //pT = 0.5-1
                                        {0.0227,     0.0252,        0.0221,      0.0210,       0.0187},      //pT = 1-2
                                        {0.0232,     0.0232,        0.0209,      0.0187,       0.0178},      //pT = 2-3
                                        {0.0236,     0.0236,        0.0219,      0.0190,       0.0184},      //pT = 3-5
                                        {0.0255,     0.0237,        0.0240,      0.0214,       0.0187}  };   //pT = 5-10     

   //Centality                          0-10%       10-20%         20-40%         40-60%        60-80%
float const pionkaonDCA_cut[6][5] ={    {0.0071,     0.0076,        0.0078,      0.0080,       0.0077},      //pT = 0-0.5
                                        {0.0064,     0.0078,        0.0073,      0.0083,       0.0077},      //pT = 0.5-1
                                        {0.0070,     0.0092,        0.0080,      0.0092,       0.0094},      //pT = 1-2
                                        {0.0063,     0.0072,        0.0093,      0.0081,       0.0078},      //pT = 2-3
                                        {0.0082,     0.0086,        0.0096,      0.0094,       0.0081},      //pT = 3-5
                                        {0.0080,     0.0085,        0.0103,      0.0106,       0.0120}  };   //pT = 5-10    








}
