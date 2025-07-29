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
#include "StPicoD0EventMaker/StPicoD0Event.h"
#include "StPicoD0EventMaker/StKaonPion.h"
#include "StPicoD0JetAnaMaker.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StCuts.h"
#include "JetInfo.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StEmcADCtoEMaker/StBemcData.h"
#include "StEmcADCtoEMaker/StEmcADCtoEMaker.h"
#include "StEmcRawMaker/StBemcRaw.h"
#include "StEmcRawMaker/defines.h"
#include "StEmcRawMaker/StBemcTables.h"
#include "BemcNewCalib.h"
#include "Calibration2016.h"
//////Refit include lib
#include "PhysicalConstants.h"
#include "StThreeVectorF.hh"
#include "StLorentzVectorD.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TF1.h"
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
//---------------------------------
#include "TRandom.h"
#include "TRandom3.h"
//---------------------------------
#include <fastjet/config.h>
#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
//#ifdef FASTJET_VERSION
#include <fastjet/Selector.hh>
#include <fastjet/tools/Subtractor.hh>
#include <fastjet/tools/Recluster.hh>
#include <fastjet/tools/JetMedianBackgroundEstimator.hh>
#include <fastjet/PseudoJet.hh>  // Pro třídu PseudoJet
#include <fastjet/ClusterSequence.hh>  // Pro klasické metody clusteringu (např. k nejbližší sousedství)
#include <fastjet/contrib/ConstituentSubtractor.hh>  // Pro ConstituentSubtractor
#include <fastjet/Selector.hh>  // Pro práci s selektory
#include <fastjet/JetDefinition.hh>  // Pro definici jetů (třeba k-means, anti-kt atd.)
#include <fastjet/AreaDefinition.hh>  // Pro definici plochy (area) v analýzách
#include <fastjet/tools/JetMedianBackgroundEstimator.hh> 

#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/contrib/IterativeConstituentSubtractor.hh" 
#include "fastjet/contrib/RescalingClasses.hh"
#include "fastjet/contrib/GenericSubtractor.hh"
#include "fastjet/contrib/ShapeWithPartition.hh"


// jet includes
#include "FJ_includes.h"
#define Error(...) ::Error(__VA_ARGS__)
#include "StJetShape.h"
#undef Error
#include "StJetPicoDefinitions.h"
#include "AngularDef.h"
//#include "/gpfs01/star/pwg/lomicond/Ondrej/Jets/Alex_install/install/fastjet-contrib-install/include/fastjet/contrib/SignalFreeBackgroundEstimator.hh"
using namespace std;
//using namespace fastjet;
namespace fj = fastjet;


//-------------------------------------

ClassImp(StPicoD0JetAnaMaker)

StPicoD0JetAnaMaker::StPicoD0JetAnaMaker(char const * name, char const * outName, StPicoDstMaker* picoDstMaker, StRefMultCorr* grefmultCorrUtil, Int_t pYear)
        : StMaker(name),
          mOutFileName(outName),
          mOutputFile(NULL),
          mGRefMultCorrUtil(grefmultCorrUtil),
          mPicoDstMaker(picoDstMaker),
          mYear(pYear),
          picoDst(NULL){}

///////////////////////////////////////////////////////////////////////////////////////////////////////
//--------Event-plane---------------
Int_t StPicoD0JetAnaMaker::EP_IsGoodTrack(StPicoTrack *trk, TVector3 pVertex){

	enum Position {West = -1, East = 1};
	Position Side = West;
	
	TVector3 mTrkMom;
	
	//Primary tracks
   	mTrkMom = trk->pMom(); //If not exists, it equals 0, 0, 0

	//0.2 < pt < 2 GeV
	Float_t pt = mTrkMom.Perp();
	if (pt!=pt||pt <= 0.2 || pt >= 2.0) return 0;

	//0.05 < |eta| < 1.00 (West/East)
	Float_t eta = mTrkMom.PseudoRapidity();
	if (abs(eta) <= 0.05 || abs(eta) >= 1) return 0;
	if (eta > 0) Side = East;

	//nHitsFit > 15
	Float_t nHitsFit = trk->nHitsFit();
	Float_t nHitsMax = trk->nHitsMax();
	
	//nHitsFit/nHitsMax => 0.52
	Float_t nHitsRatio = 1.0*nHitsFit/nHitsMax;
	if (nHitsFit < 15 || nHitsRatio < 0.52) return 0;
	
	
	//DCA < 1 cm
	Float_t dca = trk->gDCA(pVertex).Mag();
	if (dca >= 1) return 0;

return Side;
}

void StPicoD0JetAnaMaker::CalculateEventPlane(){



    //Primary vertex
    TVector3 prVertex = picoDst->event()->primaryVertex();
    //Q-vectors
    Double_t Q_1 = 0;
    Double_t Q_2 = 0;
    
 
    

    for (UInt_t iTrack = 0; iTrack < picoDst->numberOfTracks(); iTrack++){
      
        StPicoTrack *trk = static_cast<StPicoTrack *>(picoDst->track(iTrack));
        if (!trk) continue;
        
        Double_t Goodtrack = EP_IsGoodTrack(trk,prVertex);
        if (!abs(Goodtrack)) continue;
        
        TVector3 trackMom = trk->pMom();
        
        Double_t pPt = trackMom.Perp();
        Double_t phi = trackMom.Phi();
        
        //Q-vectors calculating
        Q_1 += 1.*pPt*cos(2*phi);
        Q_2 += 1.*pPt*sin(2*phi);
    }
    
   // cout << "Q_1: " << Q_1 << " Q_2: " << Q_2 << endl;
    
    fQ_1 = Q_1;
    fQ_2 = Q_2;
    
    //Recentering
    //-------------------
    //You have to run the code for "all" events to get the mean value of Q vector
        //7.2.2025
    Double_t Q_1rc = -0.5976;
    Double_t Q_2rc = 1.842;
    
    Double_t Q_1corr = Q_1 - Q_1rc;
    Double_t Q_2corr = Q_2 - Q_2rc;
    
    fQ_1_rec = Q_1corr;
    fQ_2_rec = Q_2corr;
    //-------------------
    //7.2.2025

    Double_t Psi_2 = 1./2 * TMath::ATan2(Q_2corr, Q_1corr);
   
    //-------------------
    //Psi shift	
    //You have to run the code for "all" events to get the mean value of Q vector (again)
    std::vector<Double_t> A_2 = {-0.0471131, -0.0114311, -0.000157011, 0.00128572, -0.000253291, 0.000446538, 0.00223142, 0.000574021, 0.00133264, 0.00157648, 0.00169576, -0.00023281, 0.00174849, 0.00194597, 0.000343813, 0.00104523, 0.00194369, -0.000649944, 0.00073013, -0.00138147, -0.000686723};
    std::vector<Double_t> B_2 = {0.0204715, -0.0151115, -0.00306008, -0.000966112, 0.000103265, -0.0014271, -0.00244014, 0.000801171, -0.00235772, 0.00173085, -8.08274e-05, -0.000880916, 0.000967381, 0.00182509, 0.00167076, -6.81227e-05, 8.50044e-05, -0.000525437, -0.000912995, 0.000685826, -0.000467597};

    Double_t CorrectedPsi2 = Psi_2;
    for (Int_t i = 1; i <= 21; i++){
         CorrectedPsi2 +=(1.0 / 2) * (2.0 / i) *(-A_2[i - 1] * cos(2 * i * Psi_2) + B_2[i - 1] * sin(2 * i * Psi_2));
    }
    
    if (Psi_2!=Psi_2 || CorrectedPsi2!=CorrectedPsi2) return;
    
    fPsi_2 = Psi_2;
    
    //Force the range (-pi/2,pi/2)
    CorrectedPsi2 = TMath::ATan2(TMath::Sin(2 * CorrectedPsi2), TMath::Cos(2 * CorrectedPsi2)) / 2.;
    Psi_2 = TMath::ATan2(TMath::Sin(2 * Psi_2), TMath::Cos(2 * Psi_2)) / 2.;

    
    //-------------------
    psi2 = CorrectedPsi2;
   



return;
}
//////
std::vector<fj::PseudoJet> StPicoD0JetAnaMaker::JetReconstructionShape(vector<fj::PseudoJet> fInputVectors, Int_t fCentrality, Double_t fR, Bool_t fBackSub, Double_t fEP_psi2, Int_t difiter){


	bool fPhiModulation = true;
	    //Scaling    //v2 settings
    double v2 = 0.0;
//https://journals.aps.org/prc/pdf/10.1103/PhysRevC.77.054901
//https://www.hepdata.net/record/ins777954
    switch (fCentrality) {

    case 0: v2 = 0.0696; break; //70-80%
    case 1: v2 = 0.0723; break; //60-70%
    case 2: v2 = 0.0744; break; //50-60%
    case 3: v2 = 0.074;  break; //40-50%
    case 4: v2 = 0.0703; break; //30-40%
    case 5: v2 = 0.0618; break; //20-30%
    case 6: v2 = 0.0476; break; //10-20%
    case 7: v2 = 0.0339; break; //5-10%
    case 8: v2 = 0.0232; break; //0-5%
 
    default: v2 = 0.0; break;
}
    
    double v3 = 0;
    double v4 = 0;
    double psi = fEP_psi2;

    //BackgroundRescalingYPhiUsingVectorForY(double v2, double v3, double v4, double psi, std::vector<double> values, std::vector<double> rap_binning);
    //BackgroundRescalingYPhi(): _v2(0), _v3(0), _v4(0), _psi(0), _a1(1), _sigma1(1000), _a2(0), _sigma2(1000), _use_rap(false), _use_phi(false) {}
    fj::contrib::BackgroundRescalingYPhi rescaling(v2,v3,v4,psi,0.,1.,0.,1.);
    rescaling.use_rap_term(false);    // this is useful to check if the vectors with rapidity dependence have good sizes, if one wants to use also the rapidity rescaling.
    rescaling.use_phi_term(true);
    
	

	//-----------------------------------------------------------
	
	// Define jet algorithm for actual jet reconstruction (Anti-kt with radius fR)
	fj::JetDefinition jet_def(fj::kt_algorithm, fR, fj::E_scheme, fj::Best);

	// Define area for jet reconstruction
	fj::AreaDefinition area_def_jet(fj::active_area_explicit_ghosts, fj::GhostedAreaSpec(1.2, 1, 0.01));
	////fastjet::AreaDefinition area_def_jet = area_def_jet2.with_fixed_seed(seeds);

	// Perform jet clustering with the background-subtracted (or original) particles
	fClustSeq = new fj::ClusterSequenceArea(fInputVectors, jet_def, area_def_jet);
	
	// Sort the jets by transverse momentum (pt)
	double ptmin = -1;
	vector<fastjet::PseudoJet> fInclusiveJets;
	fInclusiveJets.clear();
	fInclusiveJets = sorted_by_pt(fClustSeq->inclusive_jets(ptmin));
	
	//----------------------------------------------------------
	//--- Background Estimation ---
	fj::JetDefinition jet_def_bkgd(fj::kt_algorithm, fR, fj::E_scheme, fj::Best); // Define jet algorithm for background estimation
	fj::AreaDefinition area_def_bkgd(fj::active_area_explicit_ghosts, fj::GhostedAreaSpec(1.2, 1, 0.01)); // Define the area for background estimation

	// Decide how many jets to remove based on centrality
	int nJetsRemove = 1;
	if (fCentrality == 7 || fCentrality == 8) nJetsRemove = 2;

	// Selector to choose the hardest jets based on centrality and eta
	fj::Selector selector = (!fj::SelectorNHardest(nJetsRemove)) * fj::SelectorAbsEtaMax(0.6) * fj::SelectorPtMin(0.01);

	// Create background estimator using the previously defined selector and jet algorithm
	fj::JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);
	
	//Estimation of the background using only charged tracks
	if (fPhiModulation) bkgd_estimator.set_rescaling_class(&rescaling);
	bkgd_estimator.set_particles(fInputVectors);
	
	if (fBackSub){
	   	backgroundDensity = bkgd_estimator.rho();
	   	backgroundDensityM = bkgd_estimator.rho_m(); 
	}
	//--------------------------------------------------------------
	
	fastjet::contrib::GenericSubtractor gensub(&bkgd_estimator);
	gensub.set_common_bge_for_rho_and_rhom(true);
	fastjet::contrib::GenericSubtractorInfo info;

	Angularity my_angularity_10half(1,0.5,fR);
	Angularity my_angularity_11(1,1,fR);
	Angularity my_angularity_11half(1,1.5,fR);
	Angularity my_angularity_12(1,2,fR);
	Angularity my_angularity_13(1,3,fR);
	Angularity my_angularity_Disp(2,0,fR);
	
	fastjet::FunctionOfPseudoJet<double>* shape10half = &my_angularity_10half;
	fastjet::FunctionOfPseudoJet<double>* shape11 = &my_angularity_11;
	fastjet::FunctionOfPseudoJet<double>* shape11half = &my_angularity_11half;
	fastjet::FunctionOfPseudoJet<double>* shape12 = &my_angularity_12;
	fastjet::FunctionOfPseudoJet<double>* shape13 = &my_angularity_13;
	fastjet::FunctionOfPseudoJet<double>* shapeDisp = &my_angularity_Disp;
	
	for (unsigned int i = 0; i < fInclusiveJets.size(); i++) {
	    const fastjet::PseudoJet &jet = fInclusiveJets[i];

	    // Nejprve zkontrolujeme constituenty jetu, zda obsahuje D⁰ meson
	    std::vector<fastjet::PseudoJet> constituents = jet.constituents();
	    int d0_count = 0;
	    for (size_t j = 0; j < constituents.size(); j++) {
		// Předpokládáme, že D⁰ meson má hmotnost větší než 1 GeV
		if (constituents[j].m() > 1.0) {
		    d0_count++;
		}
	    }

	    // Pokud je v jetu právě jeden D⁰ meson, provedeme odečet angularity
	    if (d0_count == 1) {
		fAngul10half = gensub(*shape10half, jet, info);
		fAngul11 = gensub(*shape11, jet, info);
		fAngul11half = gensub(*shape11half, jet, info);
		fAngul12 = gensub(*shape12, jet, info);
		fAngul13 = gensub(*shape13, jet, info);
		fAngulDisp = sqrt(gensub(*shapeDisp, jet, info));
	    } 
	    
	}
	
	return fInclusiveJets;
	
	
}

//////
std::vector<fj::PseudoJet> StPicoD0JetAnaMaker::JetReconstructionICS(vector<fj::PseudoJet> fInputVectors, Int_t fCentrality, Double_t fR, Bool_t fBackSub, Double_t fEP_psi2, Int_t difiter){


Int_t NiterA[4] = {4,3,2,2};
Double_t R_max1A[4] = {0.05,0.125,0.100,0.150};
Double_t R_max2A[4] = {0.005,0.005,0.175,0.100};
Bool_t fMassiveTest = true;
Bool_t fPhiModulation = true;
/*****/
   //Scaling    //v2 settings
    Double_t v2 = 0.0;
//https://journals.aps.org/prc/pdf/10.1103/PhysRevC.77.054901
//https://www.hepdata.net/record/ins777954
    switch (fCentrality) {

    case 0: v2 = 0.0696; break; //70-80%
    case 1: v2 = 0.0723; break; //60-70%
    case 2: v2 = 0.0744; break; //50-60%
    case 3: v2 = 0.074;  break; //40-50%
    case 4: v2 = 0.0703; break; //30-40%
    case 5: v2 = 0.0618; break; //20-30%
    case 6: v2 = 0.0476; break; //10-20%
    case 7: v2 = 0.0339; break; //5-10%
    case 8: v2 = 0.0232; break; //0-5%
 
    default: v2 = 0.0; break;
}
    
    Double_t v3 = 0;
    Double_t v4 = 0;
    Double_t psi = fEP_psi2;

    //BackgroundRescalingYPhiUsingVectorForY(Double_t v2, Double_t v3, Double_t v4, Double_t psi, std::vector<Double_t> values, std::vector<Double_t> rap_binning);
    //BackgroundRescalingYPhi(): _v2(0), _v3(0), _v4(0), _psi(0), _a1(1), _sigma1(1000), _a2(0), _sigma2(1000), _use_rap(false), _use_phi(false) {}
    fj::contrib::BackgroundRescalingYPhi rescaling(v2,v3,v4,psi,0.,1.,0.,1.);
    rescaling.use_rap_term(false);    // this is useful to check if the vectors with rapidity dependence have good sizes, if one wants to use also the rapidity rescaling.
    rescaling.use_phi_term(true);
    
	//--- Background Estimation ---
	fj::JetDefinition jet_def_bkgd(fj::kt_algorithm, fR, fj::E_scheme, fj::Best); // Define jet algorithm for background estimation
	fj::AreaDefinition area_def_bkgd(fj::active_area_explicit_ghosts, fj::GhostedAreaSpec(1.2, 1, 0.01)); // Define the area for background estimation

	// Decide how many jets to remove based on centrality
	Int_t nJetsRemove = 1;
	if (fCentrality == 7 || fCentrality == 8) nJetsRemove = 2;

	// Selector to choose the hardest jets based on centrality and eta
	fj::Selector selector = (!fj::SelectorNHardest(nJetsRemove)) * fj::SelectorAbsEtaMax(0.6) * fj::SelectorPtMin(0.01);

	// Define seeds for random number generator (used for area calculation)
	/*
	unsigned Int_t seed1 = 12345;
	unsigned Int_t seed2 = 56789;
	std::vector<Int_t> seeds = { static_cast<Int_t>(seed1), static_cast<Int_t>(seed2) };
	*/

	// Define the area for background estimation with fixed seeds
	////fastjet::AreaDefinition area_def_bkgd = area_def_bkgd.with_fixed_seed(seeds);

	// Create background estimator using the previously defined selector and jet algorithm
	fj::JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);

	// Set the maximum eta value for background estimation
	Double_t max_eta = 1;
        
	//--- Jet Reconstruction ---
	fj::contrib::IterativeConstituentSubtractor subtractor; // Set up the background subtraction algorithm
	subtractor.set_distance_type(fj::contrib::ConstituentSubtractor::deltaR); // Set distance type for subtraction

  	// Set parameters for the background subtraction (maximum distance and alpha values)
	vector<Double_t> max_distances;
	vector<Double_t> alphas;
	
	//For higher iterations always used the same R_1max and alpha
	if (difiter < 0 || difiter >= 4) {
	    std::cerr << "ERROR: difiter out of range! difiter = " << difiter << std::endl;
	    exit(1);
	}
	if(NiterA[difiter] > 2){
	
		for (Int_t it = 0; it < NiterA[difiter] ; it++){
			max_distances.push_back(R_max1A[difiter]);
			alphas.push_back(0);
		}
	} else //else two iterations
	{
		max_distances.push_back(R_max1A[difiter]);
		max_distances.push_back(R_max2A[difiter]);
		alphas.push_back(0);
		alphas.push_back(0);
	}
	

	
	if (fBackSub&&false){
		cout << NiterA[difiter] << "test" << endl;
		for (auto& a : alphas) {
			cout << "alpha: " << a << endl;
		}
		for (auto& rr : max_distances) {	
			cout << "max_distances: " << rr << endl;
		}
	}
	
	//exclude D0
	fj::Selector notD02 = fastjet::SelectorMassMax(1); // 1 GeV max mass
	std::vector<fastjet::PseudoJet> particles_without_D0;
	std::vector<fastjet::PseudoJet> D0_pseudojet;
	std::vector<fastjet::PseudoJet> all_vectors;

	for (auto& p : fInputVectors) {
	
	if (std::isnan(p.E()) || std::isnan(p.px()) || std::isnan(p.py()) || std::isnan(p.pz())) {
	    std::cerr << "Found invalid input PseudoJet: px=" << p.px() << " E=" << p.E() << std::endl;
	}

	
	    if (notD02(p)) {
	    	fastjet::PseudoJet temp_jet; 

	    	if (fMassiveTest) temp_jet = fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.E());
	    	else temp_jet = fastjet::PseudoJet(p.px(), p.py(), p.pz(), sqrt(p.px()*p.px()+p.py()*p.py()+p.pz()*p.pz()));
	    	
		temp_jet.set_user_index(p.user_index());
	    	particles_without_D0.push_back(temp_jet);
		all_vectors.push_back(temp_jet);
	    } else {
		    fastjet::PseudoJet temp_jet; 
		    
		    if (fMassiveTest) temp_jet = fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.E());
		    else temp_jet = fastjet::PseudoJet(p.px(), p.py(), p.pz(), sqrt(p.px()*p.px()+p.py()*p.py()+p.pz()*p.pz()));
		    
		    temp_jet.set_user_index(p.user_index());
		    all_vectors.push_back(temp_jet);
		    D0_pseudojet.push_back(temp_jet);
	    }
	}
	///////////////////////////////////
	
	// Apply the subtraction parameters
	subtractor.set_parameters(max_distances, alphas);
	subtractor.set_ghost_removal(true); // Enable ghost removal
	subtractor.set_ghost_area(0.005); // Set ghost area value
	subtractor.set_max_eta(max_eta); // Set maximum eta for particles
	
	
	subtractor.set_background_estimator(&bkgd_estimator); // Link the background estimator
	
	if(fMassiveTest) subtractor.set_common_bge_for_rho_and_rhom(true); // Set common background estimation for rho and rhom
	if(fMassiveTest) subtractor.set_keep_original_masses(); // Keep the original masses of particles
	subtractor.set_scale_fourmomentum();
        
       	// Selector for particles with mass below 1 GeV (likely excludes certain particles)
	fj::Selector notD0 = fastjet::SelectorMassMax(1); // 1 GeV max mass
	subtractor.set_particle_selector(&notD0);  

	subtractor.initialize(); // Initialize the background subtraction algorithm

	// Define jet algorithm for actual jet reconstruction (Anti-kt with radius fR)
	fj::JetDefinition jet_def(fj::antikt_algorithm, fR, fj::E_scheme, fj::Best);

	// Define area for jet reconstruction
	fj::AreaDefinition area_def_jet(fj::active_area, fj::GhostedAreaSpec(1.2, 1, 0.01));
	////fastjet::AreaDefinition area_def_jet = area_def_jet2.with_fixed_seed(seeds);

	// Set particles for background estimator (background estimation for the input particles)
	////bkgd_estimator.set_particles(fInputVectors); //all_vectors
	//Rescaling
	if (fPhiModulation) bkgd_estimator.set_rescaling_class(&rescaling);
	bkgd_estimator.set_particles(all_vectors);
	
	// Minimum pt for jets (below this, jets are excluded)
	Double_t ptmin = -0.01;
	vector<fj::PseudoJet> corrected_event;

	// Apply background subtraction if enabled
	if (fBackSub) corrected_event = subtractor.subtract_event(particles_without_D0); 
	else corrected_event = particles_without_D0; // Otherwise, use original input vectors
	
	if (fBackSub){
	/*
		    // Odstranit částice s |eta| > 1 z vektoru corrected_event
		    corrected_event.erase(
			std::remove_if(corrected_event.begin(), corrected_event.end(),
				       [](const fastjet::PseudoJet& p) {
				           return (fabs(p.eta()) > 1.0 || p.perp() < 0.2);
				         //  return (fabs(p.eta()) > 1.0);
				       }),
			corrected_event.end()
		    );
		    
		    */
	
        for (vector<fastjet::PseudoJet>::const_iterator particle = corrected_event.begin(); particle != corrected_event.end(); ++particle) {
		 hJetConstRapPhiICS->Fill(particle->phi_std(), particle->rap(),weightCentrality);
		 hJetConstEtaPhiICS->Fill(particle->phi_std(), particle->eta(),weightCentrality);		
	}
	
	}
	
	// Return back D0
	//corrected_event.push_back(D0_pseudojet.back());
	if (!D0_pseudojet.empty()) corrected_event.push_back(D0_pseudojet.back());
	else cout << "D0 missing!!!" << endl;

	
	// Perform jet clustering with the background-subtracted (or original) particles
	if (fClustSeq) {
	    	delete fClustSeq;
	    	fClustSeq = nullptr;
	}
	fClustSeq = new fastjet::ClusterSequenceArea(corrected_event, jet_def, area_def_jet);
	
	// Sort the jets by transverse momentum (pt)
	vector<fastjet::PseudoJet> fInclusiveJets;
	fInclusiveJets.clear();
	fInclusiveJets = sorted_by_pt(fClustSeq->inclusive_jets(ptmin));

	backgroundDensity = bkgd_estimator.rho();



	// Return the list of inclusive jets
	return fInclusiveJets;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
Int_t StPicoD0JetAnaMaker::Init(){

  	//Output file
  	mOutputFile = new TFile(mOutFileName.Data(), "RECREATE");
  	mOutputFile->cd();

 	Jets = new TTree("jets", "Tree containing event and jet information");

	// Event
	Jets->Branch("runId", &runId, "runId/I");
	Jets->Branch("eventId", &eventId, "eventId/I");
	Jets->Branch("centrality", &centrality, "centrality/I");
	Jets->Branch("weightCentrality", &weightCentrality, "weightCentrality/F");
	Jets->Branch("gRefMult", &gRefMult, "gRefMult/I");
	Jets->Branch("backgroundDensity", &backgroundDensity, "backgroundDensity/F");
	Jets->Branch("backgroundDensityM", &backgroundDensityM, "backgroundDensityM/F");
	Jets->Branch("psi2", &psi2, "psi2/F");

	// D0 meson
	Jets->Branch("d0PdgSign", &d0PdgSign, "d0PdgSign/I");
	Jets->Branch("d0Mass", &d0Mass, "d0Mass/F");
	Jets->Branch("d0Pt", &d0Pt, "d0Pt/F");
	Jets->Branch("d0Rapidity", &d0Rapidity, "d0Rapidity/F");
	Jets->Branch("d0Eta", &d0Eta, "d0Eta/F");

	// Jet observables
	Jets->Branch("jetEta", &jetEta, "jetEta/F");
	Jets->Branch("jetPhi", &jetPhi, "jetPhi/F");
	Jets->Branch("jetRapidity", &jetRapidity, "jetRapidity/F");
	Jets->Branch("jetArea", &jetArea, "jetArea/F");
	Jets->Branch("jetPt", &jetPt, "jetPt/F");
	Jets->Branch("lambda1_0_5", &lambda1_0_5, "lambda1_0_5/F");
	Jets->Branch("lambda1_1", &lambda1_1, "lambda1_1/F");
	Jets->Branch("lambda1_1_5", &lambda1_1_5, "lambda1_1_5/F");
	Jets->Branch("lambda1_2", &lambda1_2, "lambda1_2/F");
	Jets->Branch("lambda1_3", &lambda1_3, "lambda1_3/F");
	Jets->Branch("momDisp", &momDisp, "momDisp/F");
	Jets->Branch("z", &z, "z/F");
	Jets->Branch("nJetConst", &nJetConst, "nJetConst/I");
	Jets->Branch("nJetsInEvent", &nJetsInEvent, "nJetsInEvent/I");
	Jets->Branch("jetD0DeltaR", &jetD0DeltaR, "jetD0DeltaR/F");
	Jets->Branch("jetNeutralPtFrac", &jetNeutralPtFrac, "jetNeutralPtFrac/F");
	Jets->Branch("jetTrackPtSum", &jetTrackPtSum, "jetTrackPtSum/F");

	// Event histograms:
	hVtxZ = new TH1D("hVtxZ", ";PVtx.z() [cm]; Count", 100, -10, 10);
	hVtxR = new TH2D("hVtxR", ";PVtx.x() [cm]; PVtx.y() [cm]", 100, -3, 3, 100, -3, 3);
	hVzDiff = new TH1D("hVzDiff", "V_{z} - V_{z}^{VPD};V_{z} - V_{z}^{VPD} [cm];Count", 80, -4.0, 4.0);
	hCentrality = new TH1D("hCentrality", ";C_{ID}", 9, -0.5, 8.5);
	hCentralityW = new TH1D("hCentralityW", ";C_{ID}", 9, -0.5, 8.5);
	hEventsCuts = new TH1D("hEventsCuts", "hEventsCuts;Cuts;Count", 9, 0, 9);
	hEventsCuts->GetXaxis()->SetBinLabel(1, "All events");
	hEventsCuts->GetXaxis()->SetBinLabel(2, "Triggers & Bad runs");
	hEventsCuts->GetXaxis()->SetBinLabel(3, "V_{r}");
	hEventsCuts->GetXaxis()->SetBinLabel(4, "V_{z}");
	hEventsCuts->GetXaxis()->SetBinLabel(5, "|V_{z} - V_{z}^{VPD}|");
	hEventsCuts->GetXaxis()->SetBinLabel(6, "|V_{x,y,z}|!=0");
	hEventsCuts->GetXaxis()->SetBinLabel(7, "Centrality");
	hEventsCuts->GetXaxis()->SetBinLabel(8, "Good D^{0}");
	hEventsCuts->GetXaxis()->SetBinLabel(9, "Cal.: E_{T} < 30 GeV");

	// D0 histograms:
	hD0MassPtUnlike       = new TH2D("hMassPtUnlike",";M_{K#pi} [GeV/c^{2}]; p_{T} [GeV/c]", 200, 1.6, 2.1, 100, 0, 10);
        hD0MassPtLike         = new TH2D("hMassPtLike",";M_{K#pi} [GeV/c^{2}]; p_{T} [GeV/c]", 200, 1.6, 2.1, 100, 0, 10);
        hD0EtaUnlike          = new TH1D("hD0EtaUnlike",          "Unlike-sign;#eta;Count", 100, -5, 5);
	hD0EtaLike            = new TH1D("hD0EtaLike",            "Like-sign;#eta;Count", 100, -5, 5);
	hPionEtaVsPt 	      = new TH2D("hPionEtaVsPt", "#pi pseudorapidity vs p_{T};p_{T}(#pi) [GeV/c];#eta", 150, 0, 30, 100, -5, 5);
	hKaonEtaVsPt          = new TH2D("hKaonEtaVsPt", "K pseudorapidity vs p_{T};p_{T}(K) [GeV/c];#eta", 150, 0, 30, 100, -5, 5);
	hNKaonsVsNPions       = new TH2D("hNKaonsVsNPions",       "Number of kaons vs pions;n_{#pi};n_{K}", 300, 0, 300, 300, 0, 300);

	// Topological cuts:
	hKaonDcaVsPtD0       = new TH2D("hKaonDcaVsPtD0",       "Kaon DCA vs p_{T}^{D^{0}};p_{T}(D^{0}) [GeV/c];K DCA [cm]",              120, 0, 12, 50, 0, 0.05);
	hPionDcaVsPtD0       = new TH2D("hPionDcaVsPtD0",       "Pion DCA vs p_{T}^{D^{0}};p_{T}(D^{0}) [GeV/c];#pi DCA [cm]",           120, 0, 12, 50, 0, 0.05);
	hDcaDaughtersVsPt    = new TH2D("hDcaDaughtersVsPt",    "Daughter DCA vs p_{T}^{D^{0}};p_{T}(D^{0}) [GeV/c];DCA_{K#pi} [cm]",    120, 0, 12, 300, 0, 0.03);
	hD0DcaVsPt           = new TH2D("hD0DcaVsPt",           "D^{0} DCA vs p_{T}^{D^{0}};p_{T}(D^{0}) [GeV/c];DCA_{D^{0}} [cm]",      120, 0, 12, 300, 0, 0.03);
	hCosThetaVsPt        = new TH2D("hCosThetaVsPt",        "cos(#theta) vs p_{T}^{D^{0}};p_{T}(D^{0}) [GeV/c];cos(#theta)",         120, 0, 12, 100, 0, 1.0);
	hD0DecayLengthVsPt   = new TH2D("hD0DecayLengthVsPt", "D^{0} decay length vs p_{T}^{D^{0}};p_{T}(D^{0}) [GeV/c];Decay length [cm]", 120, 0, 12, 300, 0, 0.6);

	// Daughter PID:
	hTofBetaDiffKaonVsPt = new TH2D("hTofBetaDiffKaonVsPt", "|#Delta#beta^{-1}_{K}| vs p_{T}^{D^{0}};p_{T}(D^{0}) [GeV/c];|#Delta#beta^{-1}_{K}|", 120, 0, 12, 100, 0, 0.1);
	hTofBetaDiffPionVsPt = new TH2D("hTofBetaDiffPionVsPt", "|#Delta#beta^{-1}_{#pi}| vs p_{T}^{D^{0}};p_{T}(D^{0}) [GeV/c];|#Delta#beta^{-1}_{#pi}|", 120, 0, 12, 100, 0, 0.1);
	hBetaVsSPKaon        = new TH2D("hBetaVsSPKaon",   "#beta vs sign(q) #times p (K);sign(q) #times p [GeV/c];#beta", 300, -30, 30, 100, 0, 1.2);
	hBetaVsSPPion        = new TH2D("hBetaVsSPPion",   "#beta vs sign(q) #times p (#pi);sign(q) #times p [GeV/c];#beta", 300, -30, 30, 100, 0, 1.2);
	hTpcNsigmaKaonVsPt   = new TH2D("hTpcNsigmaKaonVsPt",   "TPC n#sigma_{K} vs p_{T}^{D^{0}};p_{T}(D^{0}) [GeV/c];n#sigma_{K}",    120, 0, 12, 100, -5, 5);
	hTpcNsigmaPionVsPt   = new TH2D("hTpcNsigmaPionVsPt",   "TPC n#sigma_{#pi} vs p_{T}^{D^{0}};p_{T}(D^{0}) [GeV/c];n#sigma_{#pi}", 120, 0, 12, 100, -5, 5);
	hDedxVsSPKaon        = new TH2D("hDedxVsSPKaon",   "dE/dx vs sign(q) #times p (K);sign(q) #times p [GeV/c];dE/dx [keV/cm]", 300, -30, 30, 300, 0, 15);
	hDedxVsSPPion        = new TH2D("hDedxVsSPPion",   "dE/dx vs sign(q) #times p (#pi);sign(q) #times p [GeV/c];dE/dx [keV/cm]", 300, -30, 30, 300, 0, 15);
	hNHitsFitKaonVsD0Pt  = new TH2D("hNHitsFitKaonVsD0Pt", "nHitsFit (K) vs p_{T}^{D^{0}};p_{T}(D^{0}) [GeV/c];nHitsFit (K)", 120, 0, 12, 51, -0.5, 50.5);
	hNHitsFitPionVsD0Pt  = new TH2D("hNHitsFitPionVsD0Pt", "nHitsFit (#pi) vs p_{T}^{D^{0}};p_{T}(D^{0}) [GeV/c];nHitsFit (#pi)", 120, 0, 12, 51, -0.5, 50.5);

	//Jet constituents:
	hJetTracksDedx           = new TH2D("hJetTracksDedx",           "Track dE/dx vs signed |p|; sign(q) #times |p| [GeV/c];dE/dx [keV/cm]", 200, -4, 4, 400, 0, 20);
	hJetTracksDedxAfterCuts  = new TH2D("hJetTracksDedxAfterCuts",  "Track dE/dx vs signed |p| (after cuts);sign(q) #times |p| [GeV/c];dE/dx [keV/cm]", 400, -4, 4, 200, 0, 20);
	hJetTracksPt             = new TH1D("hJetTracksPt",         "Jet charged constituent p_{T};p_{T} [GeV/c];Count", 200, 0, 40);
	hJetTracksEtaPhi         = new TH2D("hJetTracksEtaPhi",      "Jet charged constituents #eta vs #phi;#phi;#eta", 120, -TMath::Pi(), TMath::Pi(), 100, -1.1, 1.1);
	hJetTracksNHitsFit       = new TH1D("hJetTracksNHitsFit",     "Jet charged constituents nHitsFit;nHitsFit;Count", 51, -0.5, 50.5);
	hJetTracksNHitsRatio     = new TH1D("hJetTracksNHitsRatio",   "Jet charged constituents nHitsRatio;nHitsRatio;Count", 101, 0.0, 1.01);
	hJetTracksDCA            = new TH1D("hJetTracksDCA",           "Jet charged constituents DCA;DCA [cm];Count", 100, 0, 4);
	 
	hJetNeutralPt         = new TH1D("hJetNeutralPt",        "Jet neutral constituent p_{T};p_{T} [GeV/c];Count", 200, 0, 40);
	hJetNeutralEtaPhi     = new TH2D("hJetNeutralEtaPhi",    "Jet neutral constituent #eta vs #phi;#phi;#eta", 120, -TMath::Pi(), TMath::Pi(), 100, -1.1, 1.1);
	hJetNeutralEtBefAftHC = new TH2D("hJetNeutralEtBefAftHC", "Tower E_{T} before vs after HC;E_{T}^{before} [GeV];E_{T}^{after} [GeV]", 350, 0, 35, 350, 0, 35);
	hJetNeutralECalibBefAft = new TH2D("hJetNeutralECalibBefAft","Tower E before vs after calibration;E^{before calib} [GeV];E^{after calib} [GeV]", 350, 0, 35, 350, 0, 35);
	
	hJetConstCharge = new TH1D("hJetConstCharge", "Jet constituent charge;charge;Count", 3, -1.5, 1.5);
	hJetConstRapPhi = new TH2D("hJetConstRapPhi", "Jet constituent (w/o D^{0}) y vs #phi;#phi;rapidity", 120, -TMath::Pi(), TMath::Pi(), 240, -1.2, 1.2);
	hJetConstRapPhiICS = new TH2D("hJetConstRapPhiICS", "Jet constituent (w/o D^{0}) y vs #phi after ICS;#phi;rapidity", 120, -TMath::Pi(), TMath::Pi(), 240, -1.2, 1.2);
	hJetConstEtaPhi = new TH2D("hJetConstEtaPhi", "Jet constituent (w/o D^{0}) #eta vs #phi;#phi;#eta", 120, -TMath::Pi(), TMath::Pi(), 200, -5, 5);
	hJetConstEtaPhiICS = new TH2D("hJetConstEtaPhiICS", "Jet constituent (w/o D^{0}) #eta vs #phi after ICS;#phi;#eta", 120, -TMath::Pi(), TMath::Pi(), 200, -5, 5);
	
	//Neutral particles hadronic correction:
	hJetHadrCorrNHitsFit   = new TH1D("hJetHadrCorrNHitsFit",   "Hadron. corr. track nHitsFit;nHitsFit;Count", 51, -0.5, 50.5);
	hJetHadrCorrNHitsRatio = new TH1D("hJetHadrCorrNHitsRatio", "Hadron. corr. track nHitsRatio;nHitsRatio;Count", 101, 0.0, 1.01);
	hJetHadrCorrDcaZ       = new TH1D("hJetHadrCorrDcaZ",       "Hadron. corr. track DCA_{z};DCA_{z} [cm];Count", 100, 0, 4);
	hJetHadrCorrEtaVsPt    = new TH2D("hJetHadrCorrEtaVsPt",    "Hadron. corr. track #eta vs p_{T};p_{T} [GeV/c];#eta", 200, 0, 40, 150, -1.5, 1.5);
	hJetHadrCorrE          = new TH1D("hJetHadrCorrE",    "Hadron. corr. track energy;E [GeV];Count", 200, 0, 40);
        //Loading of BEMC tables
        mADCtoEMaker = dynamic_cast<StEmcADCtoEMaker*>(GetMaker("Eread"));
        assert(mADCtoEMaker);
        mTables = mADCtoEMaker->getBemcData()->getTables();
    
        return kStOK;
}

//-----------------------------------------------------------------------------
StPicoD0JetAnaMaker::~StPicoD0JetAnaMaker(){

//Destructor
delete mGRefMultCorrUtil;
delete fClustSeq;
}

//-----------------------------------------------------------------------------
struct FourMomentum {
  //Four-momentum of the reconstructed particle
  Double_t E, px, py, pz, D0_antiD0, D0Mass;
  // D0_antiD0: 2 -> D0, -2 -> antiD0
};

//-----------------------------------------------------------------------------
Double_t delta_R(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2) {
  //Calculating of delta R = sqrt(delta eta^2 + delta phi^2)
  Double_t deta = eta1 - eta2;
  Double_t dphi = TVector2::Phi_mpi_pi(phi1 - phi2);

  return sqrt(deta*deta + dphi*dphi);

  //Function returns delta R
}

//-----------------------------------------------------------------------------
Int_t StPicoD0JetAnaMaker::Finish(){
  	//Write histograms and close the output file, if you create a new histogram it has to be added here
  	LOG_INFO << " StPicoD0JetAnaMaker - writing data and closing output file " <<endm;
  	mOutputFile->cd();


	// Event histograms:
	TDirectory* dirEvent = mOutputFile->mkdir("event");
	dirEvent->cd();
	hVtxZ->Write();
	hVtxR->Write();
	hVzDiff->Write();
	hCentrality->Write();
	hCentralityW->Write();
	hEventsCuts->Write();
	mOutputFile->cd();

	// D0 histograms:
	TDirectory* dirD0Meson = mOutputFile->mkdir("d0Meson");
	dirD0Meson->cd();
	hD0MassPtUnlike->Write();
	hD0MassPtLike->Write();
	hD0EtaUnlike->Write();
	hD0EtaLike->Write();
	hPionEtaVsPt->Write();
	hKaonEtaVsPt->Write();
	hNKaonsVsNPions->Write();
	mOutputFile->cd();

	// Topological cuts:
	TDirectory* dirTopologicalCuts = mOutputFile->mkdir("topologicalCuts");
	dirTopologicalCuts->cd();
	hKaonDcaVsPtD0->Write();
	hPionDcaVsPtD0->Write();
	hDcaDaughtersVsPt->Write();
	hD0DcaVsPt->Write();
	hCosThetaVsPt->Write();
	hD0DecayLengthVsPt->Write();
	mOutputFile->cd();

	// Daughter PID:
	TDirectory* dirDaughterPid = mOutputFile->mkdir("daughterPid");
	dirDaughterPid->cd();
	hTofBetaDiffKaonVsPt->Write();
	hTofBetaDiffPionVsPt->Write();
	hBetaVsSPKaon->Write();
	hBetaVsSPPion->Write();
	hTpcNsigmaKaonVsPt->Write();
	hTpcNsigmaPionVsPt->Write();
	hDedxVsSPKaon->Write();
	hDedxVsSPPion->Write();
	hNHitsFitKaonVsD0Pt->Write();
	hNHitsFitPionVsD0Pt->Write();
	mOutputFile->cd();

	//Jet constituents:
	TDirectory* dirJetConstituents = mOutputFile->mkdir("jetConstituents");
	dirJetConstituents->cd();
	hJetTracksDedx->Write();
	hJetTracksDedxAfterCuts->Write();
	hJetTracksPt->Write();
	hJetTracksEtaPhi->Write();
	hJetTracksNHitsFit->Write();
	hJetTracksNHitsRatio->Write();
	hJetTracksDCA->Write();
	hJetNeutralPt->Write();
	hJetNeutralEtaPhi->Write();
	hJetNeutralEtBefAftHC->Write();
	hJetNeutralECalibBefAft->Write();
	hJetConstCharge->Write();
	hJetConstEtaPhi->Write();
	hJetConstEtaPhiICS->Write();  
	hJetConstRapPhi->Write();
	hJetConstRapPhiICS->Write();     
	mOutputFile->cd();

	//Neutral particles hadronic correction:
	TDirectory* dirHadronCorr = mOutputFile->mkdir("hadronCorr");
	dirHadronCorr->cd();
	hJetHadrCorrEtaVsPt->Write();
	hJetHadrCorrE->Write();
	hJetHadrCorrNHitsFit->Write();
	hJetHadrCorrNHitsRatio->Write();
	hJetHadrCorrDcaZ->Write();
	mOutputFile->cd();

	Jets->Write();

  
	//Closing of the output file
	mOutputFile->Close();

        return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoD0JetAnaMaker::Make()
{
  //Check if everything is loaded properly
  if (!mPicoDstMaker){
    LOG_WARN << " StPicoD0JetAnaMaker - No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  picoDst = mPicoDstMaker->picoDst();
  if (!picoDst){
    LOG_WARN << "StPicoD0JetAnaMaker - No PicoDst! Skip! " << endm;
    return kStWarn;
  }
 
  // -------------- USER ANALYSIS -------------------------
  //Vertex of the event
  StThreeVectorF pVtx(-999.,-999.,-999.);

  //Loading of the event
  StPicoEvent *event = (StPicoEvent *)picoDst->event();

  //Loading of the primary vertex
  pVtx = StThreeVectorF(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z());

  //hEventsCuts: All events
  hEventsCuts->Fill(0);
  
  Bool_t isBadRun = false;
  runId = event->runId();

  //2014 good run list
  if (mYear == 2014 && mycuts::AllBadRunList2014.count(runId)) {
      isBadRun = true;
  }

        
  if(isBadRun) return kStOK;
  
  //Check if the event is good (vertex, pile-up, MB trigger)
  if(!(isGoodEvent(mYear,hEventsCuts))) return kStOK;

  //Check if the GRefMultCorr information is available
  if(!mGRefMultCorrUtil) {
    LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
    return kStWarn;
  }

  //Check if the event is in a bad run (due to BEMC towers)
  if (mYear==2016 && IsBadEnergyRun(runId)) return kStOK;

  //Loadig of the GRefMultCorr information
  mGRefMultCorrUtil->init(runId);
  mGRefMultCorrUtil->initEvent(event->grefMult(),pVtx.z(),event->ZDCx()) ;

  //Check if the event is in a bad run (due to centrality estimation)
  if (mGRefMultCorrUtil->isBadRun(runId)) return kStOK;

  //Loading of the centrality
  centrality = mGRefMultCorrUtil->getCentralityBin9();

  //Check if the centrality is in the range
  if(centrality<0) return kStOK;

  //hEventsCuts: Centrality
  hEventsCuts->Fill(6);

  //Loading of the weight for the centrality correction
  weightCentrality = mGRefMultCorrUtil->getWeight();

  //Filling events histograms
  hVtxZ->Fill(pVtx.z());
  hVtxR->Fill(pVtx.x(),pVtx.y());
  hVzDiff->Fill(event->primaryVertex().z() - event->vzVpd());
  hCentrality->Fill(centrality);
  hCentralityW->Fill(centrality,weightCentrality);
  
  ////EventStats->Fill(picoDst->event()->runId(),centrality,weightCentrality); 

  //Loading event information
  eventId = event->eventId();
  gRefMult = event->grefMult();

  //Loading of the number of tracks in the event
  UInt_t nTracks = picoDst->numberOfTracks();

//----------------------------------------------------------------

  //Variable checking if there is a good D0 candidate in the event, changed to true, if the candidate is found
  Bool_t IsThereD0 = false;

  //Preparation of the vector of the daughter candidates
  std::vector<Double_t> DaughterPionTrackVector;
  std::vector<Double_t> DaughterKaonTrackVector;

  //Print of the Fastjet banner
  fastjet::ClusterSequence::print_banner();

  //Looking for good daughter candidates
  //Vectors to hold indices of kaons and pions
  std::vector<unsigned short> idxPicoKaons;
  std::vector<unsigned short> idxPicoPions;
  
  //Loop over all tracks in PicoDst to find good kaon and pion candidates
  for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack){
  
  	//Loading track information
         StPicoTrack* trk = picoDst->track(iTrack);

         //If the track is null, skip it
         if(!trk) continue;

         //Check if the track is good for use in analysis
         if (!isGoodTrack(trk)) continue;

         //Check if the daughter pions tracks are good.
	 Bool_t tpcPion = isTpcPion(trk);

	 //Check if the daughter kaon tracks are good.
	 Bool_t tpcKaon = isTpcKaon(trk);
	 
	 //TPC information is always required.
	 if (!tpcPion && !tpcKaon) continue;

	 //Calculating of the beta for the kaon (TOF)
	 Float_t trkBeta = getTofBeta(trk,&pVtx,picoDst);

	 //Check if the there is a TOF information for the kaon
	 Bool_t trkTofAvailable = trkBeta>0;
	 
	 Bool_t tofKaon = false;
	 Bool_t tofPion = false;
	 
	 if (trkTofAvailable) {
	 	//Check if the kaon is good w.r.t the expected beta value (TOF)
    		tofKaon = isTofKaon(trk, trkBeta);
    		//Check if the pion is good w.r.t the expected beta value (TOF)
    		tofPion = isTofPion(trk, trkBeta);
	 }

	 //Final check if the kaon is good. If the TOF information is not available the more strict TPC check is used.
	 Bool_t goodKaon = tpcKaon && (!trkTofAvailable || tofKaon);
	    
	 //Final check if the pion is good. If the TOF information is not available the more strict TPC check is used.
	 Bool_t goodPion = tpcPion && (!trkTofAvailable || tofPion);

	 //Check if the kaon is good
	 if(goodKaon) idxPicoKaons.push_back(iTrack);
	    
         //Check if the pion is good
	 if(goodPion) idxPicoPions.push_back(iTrack);

  } // ... end tracks loop



	//Preparation of the four-vector of the D0 candidates
	std::vector<FourMomentum> D0_fourmomentum;

	//Get magnetic field
        Float_t const bField = event->bField();  
	  
	//Loop over all pion candidates 
	for (unsigned short ik = 0; ik < idxPicoKaons.size(); ++ik){

		//Get kaon track information
		StPicoTrack const * kaon = picoDst->track(idxPicoKaons[ik]);

		//Loop over all pion candidates
		for (unsigned short ip = 0; ip < idxPicoPions.size(); ++ip){
		
		    //If pion and kaon are the same track, skip it
		    //Sometimes it can happen that the same track is identified as both kaon and pion
		    //We do not want to pair particle with itself
		    if (idxPicoKaons[ik] == idxPicoPions[ip]) continue;
		    
		    //Get pion track information
                    StPicoTrack const * pion = picoDst->track(idxPicoPions[ip]);
                    
                    //Make a D0 candidate
                    StKaonPion kaonPion(kaon, pion, idxPicoKaons[ik], idxPicoPions[ip], pVtx, bField);

		//----------Pair-cuts--------------------------------------------------
		    //Initialisation of the charge
		    Int_t charge=0;

		    //Loading of the D0 charge and mass
		    Float_t d0Pt = kaonPion.pt();  
		    Double_t d0Mass = kaonPion.m();

		    //Upper cut on the D0 pT
		    if(d0Pt>10) continue;

		    //Check if the mass is in the D0 mass window
		    if(d0Mass<mycuts::minMass||d0Mass>mycuts::maxMass) continue;

		    //Check if all pair cuts conditions are met
		    if((charge=isD0PairCentrality_pt(kaonPion,centrality, mYear))!=0 ){
			//Charge = -1 -> Unlike-sign (pi+K- or pi-K+), Charge = 1 -> Like-sign (pi+K+), Charge = 2 -> Like-sign (pi-K-)

			//If pair is Unlike-sign
			if(charge==-1){

			    //Filling of the D0 histograms
			    hD0MassPtUnlike->Fill(d0Mass,d0Pt,weightCentrality);
			    hD0EtaUnlike->Fill(kaonPion.eta(),weightCentrality);
		      
			    //Saving of the daughter tracks
			    DaughterPionTrackVector.push_back(pion->id());
			    DaughterKaonTrackVector.push_back(kaon->id());

			    //The event is noted
			    IsThereD0 = true;

			    //Saving of D0 four-momenta // E,       px,       py,        pz,         D0_antiD0,            D0 mass;
			    FourMomentum D0_actual = {kaonPion.Energy(), kaonPion.Px(), kaonPion.Py(),  kaonPion.Pz(), (Double_t)pion->charge(), d0Mass};
			    D0_fourmomentum.push_back(D0_actual);
			 

			    //Loading of the rescaled RunID
			    //Int_t runIndex = mPrescales->runIndex(mPicoD0Event->runId());

			    //Filling of the histograms
			    hPionEtaVsPt->Fill(pion->gMom().Perp(), pion->gMom().PseudoRapidity(),weightCentrality);
			    hKaonEtaVsPt->Fill(kaon->gMom().Perp(), kaon->gMom().PseudoRapidity(),weightCentrality);
			    hNKaonsVsNPions->Fill(idxPicoPions.size(),idxPicoKaons.size(),weightCentrality);
			    hKaonDcaVsPtD0->Fill(kaonPion.pt(),kaonPion.kaonDca(),weightCentrality);
			    hPionDcaVsPtD0->Fill(kaonPion.pt(),kaonPion.pionDca(),weightCentrality);
			    hDcaDaughtersVsPt->Fill(kaonPion.pt(),kaonPion.dcaDaughters(),weightCentrality);
			    hD0DcaVsPt->Fill(kaonPion.pt(),sin(kaonPion.pointingAngle())*kaonPion.decayLength(),weightCentrality);
			    hCosThetaVsPt->Fill(kaonPion.pt(),cos(kaonPion.pointingAngle()),weightCentrality);
			    hD0DecayLengthVsPt->Fill(kaonPion.pt(),kaonPion.decayLength(),weightCentrality);
			    
			    // Daughter PID:
			    Double_t kaonMomTot = kaon->gMom().Mag();    
		            Double_t pionMomTot = pion->gMom().Mag();

			    //Calculating of the beta for the kaon (TOF)
			    Float_t kBeta = getTofBeta(kaon,&pVtx,picoDst);
				    
			    //Calculating of the beta for the pion (TOF)
			    Float_t pBeta = getTofBeta(pion,&pVtx,picoDst);
			    
			    //Check if the there is a TOF information for the kaon
			    Bool_t kTofAvailable = kBeta>0;
		    
			    //Check if the there is a TOF information for the pion
			    Bool_t pTofAvailable = pBeta>0;

			    if (kTofAvailable){
			    Double_t beta_k = kaonMomTot/sqrt(kaonMomTot*kaonMomTot+M_KAON_PLUS*M_KAON_PLUS);
			    hTofBetaDiffKaonVsPt->Fill(d0Pt,abs(1/kBeta-1/beta_k),weightCentrality);
			    hBetaVsSPKaon->Fill(kaonMomTot*kaon->charge(), kBeta, weightCentrality);
			    }
			    
			    if (pTofAvailable){
			    Double_t beta_p = pionMomTot/sqrt(pionMomTot*pionMomTot+M_PION_PLUS*M_PION_PLUS);
			    hTofBetaDiffPionVsPt->Fill(d0Pt,abs(1/pBeta-1/beta_p),weightCentrality);
			    hBetaVsSPPion->Fill(pionMomTot*pion->charge(), pBeta, weightCentrality);
			    }
			    
			    hTpcNsigmaKaonVsPt->Fill(d0Pt, kaon->nSigmaKaon(), weightCentrality);
			    hTpcNsigmaPionVsPt->Fill(d0Pt, pion->nSigmaPion(), weightCentrality);
			    hDedxVsSPKaon->Fill(kaonMomTot*kaon->charge(), kaon->dEdx(), weightCentrality);
		 	    hDedxVsSPPion->Fill(pionMomTot*pion->charge(), pion->dEdx(), weightCentrality);
		 	    hNHitsFitKaonVsD0Pt->Fill(d0Pt, kaon->nHitsFit(), weightCentrality);
			    hNHitsFitPionVsD0Pt->Fill(d0Pt, pion->nHitsFit(), weightCentrality);

			} //end of the Unlike-sign if

			//If pair is Like-sign
			if(charge>0){

			    //Filling of the D0 histograms
			    hD0MassPtLike->Fill(d0Mass,d0Pt,weightCentrality);
			    hD0EtaLike->Fill(kaonPion.eta(), weightCentrality);

			} //end of the Like-sign if
		       
		    } //end of the pair cuts if

		  } //end of the D0 candidate loop pions
	} //end of the D0 candidate loop kaons

//---------------------------------------------------------------------
//-------------------JET-RECONSTRUCTION-PART---------------------------
//---------------------------------------------------------------------

  //Check if there is a D0 candidate in the event
  if(IsThereD0){

  CalculateEventPlane();

      //hEventsCuts: Good D0 candidate
      hEventsCuts->Fill(7);

      //Initialisation of the input particle vectors for FastJet
      vector<fastjet::PseudoJet> input_particles;
      vector<fastjet::PseudoJet> chargedjetTracks;
      vector<fastjet::PseudoJet> neutraljetTracks;

      //Radius of the jet
      Double_t R = 0.4;

      //Initialize the tower geometry
      //StEmcPosition* mEmcPosition;
      //mEmcPosition = new StEmcPosition();
      StEmcPosition mEmcPosition;

      //Loop over all D0 candidates in the event.
      //If there are more than one, the jet reconstruction is done for each D0 candidate separately
      //ignoring the other not reconstructed D0 candidates in the event.
      for (UInt_t nD0 = 0; nD0 < D0_fourmomentum.size(); nD0++) {

 	//for (Int_t iTow = 0; iTow < 4800; iTow++) SumE[iTow] = 0;

	//Delete all energies calculated for hadr. corr. from previous event
	SumE.fill(0);
	
	//for (Int_t iTow = 0; iTow < 4800; iTow++) cout << SumE[iTow] << endl;
//-----------D0-track--------------------------------------------------------

        //Defining the four-momentum of the D0 candidate
        fastjet::PseudoJet pj(D0_fourmomentum[nD0].px,D0_fourmomentum[nD0].py,D0_fourmomentum[nD0].pz,D0_fourmomentum[nD0].E); //comp. testing

        //Set flag to 2 if the D0 candidate is a D0 and to -2 if it is a anti-D0
        //It cannot be -1 or 1 because the default flag in FastJet is -1.
        pj.set_user_index(D0_fourmomentum[nD0].D0_antiD0*2);
  	
        //Add the D0 candidate to the inclusive particle vector
        input_particles.push_back(pj);
        
//-----------Neutral-tracks--------------------------------------------------

        //Fill array SumE with momenta of tracks which are matched to BEMC towers
        GetCaloTrackMomentum(picoDst,TVector3(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z()));

        //Loop over all tracks in the event
        for (Int_t iTow = 0; iTow < 4800; iTow++){

           //Get the tower hit
           StPicoBTowHit *towHit = picoDst->btowHit(iTow);

           //Check if the tower is bad or missing information
           if (!towHit || towHit->isBad()) continue;

           //Get the alternative counting method for the tower ID
           //NumericIndex - Counting from 0 to 4799
           //SoftId - Counting from 1 to 4800
           Int_t realtowID = towHit->numericIndex2SoftId(iTow);

           //Initialize the tower energy
           Double_t towE = 0;
           Double_t towEUncal = 0;

           //Calculation of the tower energy depending on the year
           //In 2014, there was a problem with the energy calibration, so the energy has to be corrected
           
           if(mYear==2014) {
                  //Exclude bad towers, saved in JetInfo.h
                  if (BadTowerMap[realtowID-1]) continue;

                  //Calculate the tower energy
                  towE = GetTowerCalibEnergy(realtowID); 
                  towEUncal = towHit->energy();         

              }
           if(mYear==2016) {
                  //Exclude bad towers, saved in Calibration2016.h
                  if (EnergyBadTowerMap2016[realtowID-1]<0) continue;
                 
                  //Get the tower energy
                  towE = towHit->energy();
                  towEUncal = towE;
           }
           
           Double_t towEUncorr = towE;
           
           towE-= fHadronCorr*SumE[iTow];
           
           //If the tower energy is negative, set it to 0
           if (towE < 0) towE = 0;

           //Correct the eta of the tower for the vertex position
           //Because the loaded eta is w.r.t. the center of the TPC, but the vertex do not have to be in the center
           StThreeVectorF towerPosition = mEmcPosition.getPosFromVertex(StThreeVectorF(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z()), realtowID);
           Float_t Toweta = towerPosition.pseudoRapidity();
           //Float_t Towphi = towerPosition.phi();

           //Calculate the transverse energy
           //max eta 1.05258 max difference: ET = 0.124452 for E = 0.2, if we cut on |Vz| < 30 cm
           Double_t ET = towE/cosh(Toweta);
           Double_t ETUncorr = towEUncorr/cosh(Toweta);

           //If the transverse energy is greater than 30 GeV, discard the event
           if (ET > 30) {
                return kStOK;
           }
////////////////////////////////////////////////////////////////////////////////// comp. test
	   Double_t p = 1.0*TMath::Sqrt(towE*towE - 0*0);
	   Double_t posX = towerPosition.x();
  	   Double_t posY = towerPosition.y();
  	   Double_t posZ = towerPosition.z();
  	     Double_t r = TMath::Sqrt(posX*posX + posY*posY + posZ*posZ) ;
  	                Double_t px,py,pz;
  	   px = p*posX/r*1.;
    	   py = p*posY/r*1.;
    	   pz = p*posZ/r*1.;

           //Create a jet with the calculated momentum components
           PseudoJet inputTower(px, py, pz, towE);
           //Discarding of the towers with low transverse energy (Defined in RunPicoD0AnaMaker.C as setCutETmin)
           //TrackBasedJets = charged tracks + D0 (Defined in RunPicoD0AnaMaker.C as setOnlyTrackBasedJets)
           if (ET > fETmincut && OnlyTrackBasedJets == 0){
                   
               //Set the flag to 10 if the particle is neutral
               inputTower.set_user_index(10);
               //Add the neutral particle to the neutral particle vector
               neutraljetTracks.push_back(inputTower);
               //Add the neutral particle to the inclusive particle vector
               input_particles.push_back(inputTower);
               
               ////cout << "input_particles.push_back({"<<px<<","<<py<<","<<pz<<","<<towE<<"}); //neutral" << endl;
               //cout << px <<"," << py << "," << pz << endl;
               
               TVector3 NeutralPart(px, py, pz);
               hJetNeutralPt->Fill(NeutralPart.Perp(), weightCentrality);
	       hJetNeutralEtaPhi->Fill(NeutralPart.Phi(), NeutralPart.Eta(), weightCentrality);
	       hJetNeutralEtBefAftHC->Fill(ETUncorr, ET, weightCentrality);
	       hJetNeutralECalibBefAft->Fill(towEUncal, towEUncorr, weightCentrality);
	       hJetConstCharge->Fill(0., weightCentrality);
	       hJetConstRapPhi->Fill(inputTower.phi_std(),inputTower.rap(), weightCentrality); 
	       hJetConstEtaPhi->Fill(inputTower.phi_std(),inputTower.eta(), weightCentrality); 
           } //End of minimum ET cut

        } //End of loop over all towers

        //hEventsCuts: Cal.: E_{T} < 30 GeV
        //Only for the first D0 otherwise it is counted multiple times
        if (nD0==0) hEventsCuts->Fill(8);
        

//-----------Charged-tracks--------------------------------------------------




        //Loop over all tracks in the event
        for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack){

            //The i-th track is loaded
            StPicoTrack* trk = picoDst->track(iTrack);

            //Check if the track exists
            if(!trk) continue;

            //Filling jet track histogram before any cuts
            hJetTracksDedx->Fill(trk->gPtot()*trk->charge(),trk->dEdx(),weightCentrality);
            

            //Check if the track is a good track
            if (!isGoodJetTrack(trk,event)) continue;

            //Loading of the pT
            Double_t pT = trk->gMom().Perp();
            //Check if the pT is above 0.2 GeV/c or if it is NaN, because NaN!=NaN
            if(pT != pT) continue; // NaN test.
            //Loading of the eta, phi, dca and charge
            Float_t eta = trk->gMom().PseudoRapidity();
            Float_t phi = trk->gMom().Phi();
            Float_t dca = (TVector3(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z()) - trk->origin()).Mag();
            Float_t charged = trk->charge();

            //Filling jet track histogram before after cuts
            hJetTracksDedxAfterCuts->Fill(trk->gPtot()*trk->charge(),trk->dEdx(),weightCentrality);

            //If the track is not daughter pion nor kion then...
            if (DaughterPionTrackVector[nD0] != trk->id() && DaughterKaonTrackVector[nD0] != trk->id()){

                //Defining the four-momentum of the charged particle, assumed pi+ mass
                //                       px,       py,               pz,                                 E = sqrt(p^2 + m^2),
                fastjet::PseudoJet pj(trk->gMom().x(),trk->gMom().y(),trk->gMom().z(), sqrt(trk->gMom().Mag()*trk->gMom().Mag()+M_PION_PLUS*M_PION_PLUS));
       
		//Set the flag to 3 if the particle is charged
		pj.set_user_index(3);
                //Add the charged particle to the charged particle vector
                chargedjetTracks.push_back(pj);
                //Add the charged particle to the inclusive particle vector
                input_particles.push_back(pj);
                
                //QA histograms
                hJetTracksPt->Fill(pT,weightCentrality);
                hJetTracksEtaPhi->Fill(phi, eta, weightCentrality);
                hJetTracksNHitsFit->Fill(trk->nHitsFit(), weightCentrality);
		hJetTracksNHitsRatio->Fill(1.0*trk->nHitsFit()/trk->nHitsMax(), weightCentrality);
		hJetTracksDCA->Fill(dca, weightCentrality);
		hJetConstCharge->Fill(charged, weightCentrality);
		hJetConstRapPhi->Fill(pj.phi_std(),pj.rap(), weightCentrality);
		hJetConstEtaPhi->Fill(pj.phi_std(),pj.eta(), weightCentrality);

            }//End of if the track is not daughter pion nor kion

     } //End of loop over all tracks

//-----------Background-estimation--------------------------------------------------
        //vector<fastjet::PseudoJet> input_particles2 =input_particles;
        //vector<fastjet::PseudoJet> input_particles3 =input_particles;
/*

        //Definition of jets for background estimation
        //Contrary to the inclusive jets, the background jets are reconstructed with the kt algorithm (recommended choice)
        JetDefinition jet_def_bkgd(kt_algorithm, R, E_scheme, Best);
        
        //Definition of the area for background estimation
        AreaDefinition initial_area_def(active_area,GhostedAreaSpec(1., 1, 0.01)); //Comp. test
        //AreaDefinition initial_area_def(active_area_explicit_ghosts,GhostedAreaSpec(1.2, 1, 0.005)); //Comp. test
        /////
        
        unsigned Int_t seed1 = 12345;
    	unsigned Int_t seed2 = 56789;
	std::vector<Int_t> seeds = { static_cast<Int_t>(seed1), static_cast<Int_t>(seed2) };
    	fastjet::AreaDefinition area_def_bkgd = initial_area_def.with_fixed_seed(seeds);
    	
	//////

        //Remove two hardest jets in central collisions, one in others
        if (centrality == 7 || centrality == 8) nJetsRemove = 2; //Comp. test
	//nJetsRemove = 2; //Comp. test

        //Definition of the selector for background estimation (eta and pt cut + remove the n hardest jets)
        //Selector selector = (!SelectorNHardest(nJetsRemove)) * SelectorAbsEtaMax(1.) * SelectorPtMin(0.01); //Comp. test
        Selector selector = (!SelectorNHardest(nJetsRemove)) * SelectorAbsEtaMax(0.6); //Comp. test
        
        //Definition of the background estimator
        JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);
        ////JetMedianBackgroundEstimator bkgd_estimator2(selector, jet_def_bkgd, area_def_bkgd);
        //Estimation of the background using only charged tracks
        bkgd_estimator.set_particles(input_particles); //Bug previously: chargedjetTracks instead of input_particles

        //Calculation of the rho (median) and sigma (fluctuations of the median) for the background
        Float_t rho = bkgd_estimator.rho();
        //Float_t emptyjets = bkgd_estimator.n_empty_jets(); //not used
        //Float_t alljets = bkgd_estimator.n_jets_used(); //not used
        Float_t rhom = bkgd_estimator.rho_m(); 
        
        //iterativ background estimation
        Double_t max_eta=100;   
       	Double_t max_eta_jet=3; // the maximal pseudorapidity for selected jets. Not important for the subtraction.


 	contrib::IterativeConstituentSubtractor subtractor;
  	subtractor.set_distance_type(contrib::ConstituentSubtractor::deltaR); 




        vector<Double_t> max_distances;
	max_distances.push_back(0.15);
	//max_distances.push_back(0.2);
    
	vector<Double_t> alphas;
	alphas.push_back(1);
	//alphas.push_back(1);
	subtractor.set_parameters(max_distances,alphas); // in this example, 2 CS corrections will be performed: 1. correction with max_distance of 0.1 and alpha of 0, 2. correction with max_distance of 0.15 and alpha of 0. After the first correction, the scalar sum of pt from remaining ghosts is evaluated and redistributed uniformly accross the event.
	
	subtractor.set_ghost_removal(true);  
  	subtractor.set_ghost_area(0.005);
	subtractor.set_max_eta(max_eta);
	subtractor.set_background_estimator(&bkgd_estimator);
        subtractor.set_common_bge_for_rho_and_rhom(true);
        subtractor.set_keep_original_masses();    
        
        // Example selector for ConstituentSubtractor:                     
        fastjet::Selector notD0 = fastjet::SelectorMassMax(1); //1 GeV max
        subtractor.set_particle_selector(&notD0);  

        
        subtractor.initialize();
         
       //  cout << subtractor.description() << endl;
        


 
 
//------Alternative-background-subtraction-------------------------------------------

//-----------Jet-reconstruction-and-variable-calculations----------------------------

        //Inclusive (or track based) jet definition
        fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R, E_scheme, Best);

        //Definition of the area for jet reconstruction
	fastjet::AreaDefinition initial_area_def2(fastjet::active_area_explicit_ghosts, GhostedAreaSpec(1., 1, 0.01)); //Comp. test
	//fastjet::AreaDefinition initial_area_def2(fastjet::active_area_explicit_ghosts, GhostedAreaSpec(1.2, 1, 0.005)); //Comp. test
	
	/////
        //unsigned Int_t seed1 = 12345;
    	//unsigned Int_t seed2 = 56789;
	//std::vector<Int_t> seeds = { static_cast<Int_t>(seed1), static_cast<Int_t>(seed2) };
    	fastjet::AreaDefinition area_def = initial_area_def2.with_fixed_seed(seeds);
	//////
	
	
        //Definition of the clustering
        ClusterSequenceArea clust_seq_hard(input_particles, jet_def, area_def);
        
        //Jet minimum pT cut
        Double_t ptmin = -0.1;
        //  fFilteredJets = fClustSeqSA->inclusive_jets(fMinJetPt-0.1); //becasue this is < not <=
        

        //Sorting of the jets by pT
        vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq_hard.inclusive_jets(ptmin-0.01));


*/
	
	Int_t difiter = 2;
	
	vector<PseudoJet> corrected_jets;
	
	if (fBgSubtraction == 1) corrected_jets = JetReconstructionICS(input_particles, centrality, R, true, psi2, difiter);
	if (fBgSubtraction == 0 || fBgSubtraction == 2) corrected_jets = JetReconstructionShape(input_particles, centrality, R, true, psi2, difiter);

        //Loop over all jets
        for (UInt_t i = 0; i < corrected_jets.size(); i++) {
          
      	    //Exclude jets with |eta| > 1 - R
            ////if (abs(inclusive_jets[i].pseudorapidity()) > (1.0 - R)) continue; Postponed to local analysis

	    //Constituent counter
	    Int_t ConstCounter = 0;

            //Initialize the variables (angularities, z-value)
            Int_t user_index = 0;
            Double_t Delta_R_D0 = 0;
            Double_t zet = 0;
            Double_t lambda_1_0half = 0;
            Double_t lambda_1_1 = 0;
            Double_t lambda_1_1half = 0;
            Double_t lambda_1_2 = 0;
            Double_t lambda_1_3 = 0;
            Double_t lambda_2_0 = 0;
            Double_t neutralpT = 0;

            //Calculate the jet pT + background subtraction (= _corr)
            //pT(sub) = pT - rho * A_jet
            //Since z-value requires px and py, it is calculated as well
            Double_t pT_jet = corrected_jets[i].perp();
            Double_t pT_jet_corr = -999;
            if (fBgSubtraction == 1) pT_jet_corr = pT_jet; //ICS
            else if (fBgSubtraction == 0 || fBgSubtraction == 2) pT_jet_corr = pT_jet - backgroundDensity * corrected_jets[i].area(); //Area+Shape
            
            
            Double_t px_jet = corrected_jets[i].px();
            Double_t px_jet_corr = -999;
            
            if (fBgSubtraction == 1) px_jet_corr = px_jet; //ICS
            else if (fBgSubtraction == 0 || fBgSubtraction == 2) px_jet_corr = px_jet - backgroundDensity * corrected_jets[i].area_4vector().px(); //Area+Shape
            

            Double_t py_jet = corrected_jets[i].py();
            Double_t py_jet_corr = -999;
            
            if (fBgSubtraction == 1) py_jet_corr = py_jet; //ICS
            else if (fBgSubtraction == 0 || fBgSubtraction == 2) py_jet_corr = py_jet - backgroundDensity * corrected_jets[i].area_4vector().py(); //Area+Shape
                        
 
	    Double_t pT_jet_scalar = 0;
            const vector<fastjet::PseudoJet>& constituents = corrected_jets[i].constituents();
		
            //Loop over all constituents of the i-th jet
            for (vector<fastjet::PseudoJet>::const_iterator particle = constituents.begin(); particle != constituents.end(); ++particle) {

                //Loading of the particle index
                Int_t index = particle->user_index();
                
                //Number of constituents (+-2 = D0, 3 = charged, 10 = neutral, -1 ghost)
                if (index == -1) continue;
                //if (particle->pt() < 0.01); continue;
                
                //Scalar
		pT_jet_scalar += particle->pt();
		
                //Constituent counter
                if(particle->pt() > 0.001) ConstCounter++;
                
                //Fraction of neutral particles
                if (index == 10) neutralpT += particle->pt();
                
                //Calculating the delta R = sqrt(delta eta^2 + delta phi^2)
                Double_t Delta_R =delta_R(corrected_jets[i].eta(),corrected_jets[i].phi(),particle->eta(),particle->phi());

                //Angularities are calculated only for track based particles (charged + D0)
                if (particle->user_index() != 10) {
                    lambda_1_0half+=	pow(particle->pt()/pT_jet_corr, 1)*		pow( Delta_R /R ,0.5);
                    lambda_1_1+=	pow(particle->pt()/pT_jet_corr, 1)*		pow( Delta_R /R ,1. );
                    lambda_1_1half+=	pow(particle->pt()/pT_jet_corr, 1)*		pow( Delta_R /R ,1.5);
                    lambda_1_2+=	pow(particle->pt()/pT_jet_corr, 1)*		pow( Delta_R /R ,2. );
                    lambda_1_3+=	pow(particle->pt()/pT_jet_corr, 1)*		pow( Delta_R /R ,3. );
 		    lambda_2_0+=	pow(particle->pt()/pT_jet_corr, 2)*		pow( Delta_R /R ,0. );
                }

                //Check if the constituent is D0 (D0 = 2, antiD0 = -2)
                if (abs(index) == 2 ) {
		    
                    //user_index is used as a flag to check if the jet contains D0
                    user_index = index;

                    //Delta R for D0
                    Delta_R_D0 = Delta_R;

                    //z-value, z = pT(D0)*^pT(jet)/|pT(jet)|
                    zet = (D0_fourmomentum[nD0].px*px_jet_corr+D0_fourmomentum[nD0].py*py_jet_corr)/(pT_jet_corr*pT_jet_corr);
                    
                }

            } //End of loop over all constituents of the i-th jet

            //Calculation of the fraction of neutral particles
            Double_t nfraction = 1.*neutralpT/pT_jet;

            //If the fraction is too high, jet is rejected (Parameter is saved in RunPicoD0AnaMaker.C as setMaxNeutralFraction)
            ////if (nfraction > maxneutralfrac) continue; //Postponed to local analysis
            
            //if (area_jet < fAcuts) continue; //Not implemented

            //If the jet contains D0
            if (abs(user_index) ==2){

                //Calculation of the D0 mass
                Double_t D0mass = D0_fourmomentum[nD0].D0Mass;
                //Calculation of the D0 pT (pT=sqrt(px^2+py^2))
                Double_t D0_pT = sqrt(D0_fourmomentum[nD0].px*D0_fourmomentum[nD0].px+D0_fourmomentum[nD0].py*D0_fourmomentum[nD0].py);
                
                //Rapidity calculations and filling histogram
                TLorentzVector v(D0_fourmomentum[nD0].px,D0_fourmomentum[nD0].py,D0_fourmomentum[nD0].pz,D0_fourmomentum[nD0].E);
                //Double_t D0_rapidity = 1./2.*log((D0_fourmomentum[nD0].E+D0_fourmomentum[nD0].pz)/(D0_fourmomentum[nD0].E-D0_fourmomentum[nD0].pz));
  		Double_t D0_rapidity = v.Rapidity();
                Double_t D0_pseudorapidity = v.PseudoRapidity();
                
		// D0 meson
		d0PdgSign     = Int_t(user_index/2.);  //1 pro D0, -1 pro anti-D0
		d0Mass        = D0mass;
		d0Pt          = D0_pT;
		d0Rapidity    = D0_rapidity;
		d0Eta         = D0_pseudorapidity;

		// Jet observables
		jetEta            = corrected_jets[i].pseudorapidity();
		jetPhi            = corrected_jets[i].phi();
		jetRapidity       = corrected_jets[i].rap();
		jetArea           = corrected_jets[i].area();
		jetPt             = pT_jet;
		if (fBgSubtraction == 0 || fBgSubtraction ==1){ //Area+ICS
			lambda1_0_5       = lambda_1_0half;
			lambda1_1         = lambda_1_1;
			lambda1_1_5       = lambda_1_1half;
			lambda1_2         = lambda_1_2;
			lambda1_3         = lambda_1_3;
			momDisp           = sqrt(lambda_2_0);
		} else if (fBgSubtraction ==2){
			lambda1_0_5       = fAngul10half;
			lambda1_1         = fAngul11;
			lambda1_1_5       = fAngul11half;
			lambda1_2         = fAngul12;
			lambda1_3         = fAngul13;
			momDisp           = fAngulDisp;	//sqrt done previously	
		}
		z                 = zet;
		nJetConst         = ConstCounter;
		nJetsInEvent      = D0_fourmomentum.size();
		jetD0DeltaR       = Delta_R_D0;
		jetNeutralPtFrac  = nfraction;
		jetTrackPtSum     = pT_jet_scalar;

		// Nakonec naplň strom
		Jets->Fill();
   
            } else{

                //Print the jet information
               //// printf("%5u %15.8f %15.8f %15.8f %15d\n", i, corrected_jets[i].rap(), corrected_jets[i].phi(), corrected_jets[i].perp(), user_index);

            } //End of if the jet contains D0

        } //End of loop over all jets

        corrected_jets.clear();

        //Delete vectors to avoid a pileup
        ////inclusive_jets.clear();
        input_particles.clear();
        chargedjetTracks.clear();
        neutraljetTracks.clear();


      } //End of D0 loop

  } //End of IsthereD0 condition

  //Delete daughter particles vectors
  DaughterPionTrackVector.clear();
  DaughterKaonTrackVector.clear();

  //End of the event
  return kStOK;
}

//---------------------------------------------------------------------
//-----------------------------FUNCTIONS-------------------------------
//---------------------------------------------------------------------
 //Not used anymore
Double_t StPicoD0JetAnaMaker::vertexCorrectedEta(Double_t eta, Double_t vz) {
    //Function to correct the eta value of a track for the z-position of the primary vertex

    //eta = -log(tan(theta/2)) => theta = 2*atan(exp(-eta))
    Double_t tower_theta = 2.0 * atan(exp(-eta));

    //If eta = 0 then z = 0
    //Else calculate z position
    Double_t z = 0.0;
    if (eta != 0.0) z = mBarrelRadius / tan(tower_theta);

    //Difference between the z position of the track and the z position of the primary vertex
    Double_t z_diff = z - vz;

    //Calculate the corrected theta value
    Double_t theta_corr = atan2(mBarrelRadius, z_diff);

    //Calculate the corrected eta value
    Double_t eta_corr = -log(tan(theta_corr / 2.0));

    return eta_corr;

    //Function returns the corrected eta value
}
//---------------------------------------------------------------------------
Bool_t StPicoD0JetAnaMaker::GetCaloTrackMomentum(StPicoDst *mPicoDst, TVector3 mPrimVtx) {
  //Function to calculate the momentum of the track matched to the calorimeter tower

  //Loading of the number of tracks in the event
  UInt_t nTracks = mPicoDst->numberOfTracks();

  //Loop over all tracks in the event
  for (UInt_t itrack = 0; itrack < nTracks; itrack++) {

      //Loading the track
      StPicoTrack *trk = mPicoDst->track(itrack);
      //Loading of the global momentum
      TVector3 gMom = trk->gMom();

      //Loading of the pT
      Double_t pT = gMom.Perp();
           
      //Check if the pT is above 0.2 GeV/c or if it is NaN, because NaN!=NaN
      if(pT != pT || pT < 0.2) continue;
      //Loading of eta
      Float_t eta = gMom.PseudoRapidity();
      //Exclude tracks outside of the TPC acceptance
      if (fabs(eta) > 1) continue;
      //Loading of phi
      //Float_t phi = gMom.Phi();

      //Loading of the number of hits
      Float_t nHitsFit = trk->nHitsFit();
      //Loading of the number of hits possible
      Float_t nHitsMax = trk->nHitsMax();
      //Exclude tracks with less than 15 hits or with a ratio of less than 0.52
      if (nHitsFit < 15 || 1.*nHitsFit/nHitsMax < 0.52) continue;

      //Loading of the value of the magnetic field
      Double_t Bfield = mPicoDst->event()->bField();
      //Loading of the helix
      StPicoPhysicalHelix trkhelix = trk->helix(Bfield);

      //Loading of the primary vertex
      Float_t vtx_x = mPrimVtx.x();
      Float_t vtx_y = mPrimVtx.y();
      Float_t vtx_z = mPrimVtx.z();

      //Calculation of the DCA to the primary vertex
      TVector3 dcaPoint = trkhelix.at(trkhelix.pathLength(vtx_x, vtx_y));
      //Calculation of the DCA in the x-y plane
      ////Float_t dca_z = dcaPoint.z() - vtx_z; //Test
      
      Float_t dca_z = trk->gDCAz(vtx_z); //? TO DO
      
      //Exclude tracks with a DCA to the primary vertex in z of more than maxdcazhadroncorr (in RunPicoD0AnaMaker.C)
      if (fabs(dca_z) > maxdcazhadroncorr) continue;

      //Initialization and loading of the tower index
      Int_t TowIndex = -99999;
      TowIndex = trk->bemcTowerIndex(); //ID
      Float_t p = 0;
      
      //Check if the track is matched to a tower
      if (TowIndex >= 0) {

        //Loading of the momentum
        p = gMom.Mag();          
	Double_t TrackEnergy = 1.0*TMath::Sqrt(p*p + M_PION_PLUS*M_PION_PLUS);
	
	hJetHadrCorrNHitsFit->Fill(nHitsFit, weightCentrality);
	hJetHadrCorrNHitsRatio->Fill(1.*nHitsFit/nHitsMax, weightCentrality);
	hJetHadrCorrDcaZ->Fill(dca_z, weightCentrality);
	hJetHadrCorrEtaVsPt->Fill(pT, eta, weightCentrality);
	hJetHadrCorrE->Fill(TrackEnergy, weightCentrality);
        
        //Summing up the energy of all tracks matched to the same tower //Previously neglected pion mass 
        SumE[TowIndex] += TrackEnergy;
      }
  
  } //End of track loop

  return true;

  //Function returns true if it was successful and SumE filled with the momentum of all tracks matched to a tower
}
//---------------------------------------------------------------------------
Double_t StPicoD0JetAnaMaker::GetTowerCalibEnergy(Int_t TowerId){
  //Function calculates the calibrated energy of a tower

  //Loading of the tower
  StPicoBTowHit *tower = static_cast<StPicoBTowHit*>(picoDst->btowHit(TowerId-1));

  //Initialization of the pedestal, rms and status
  Float_t pedestal, rms;
  Int_t status;

  //Loading of the pedestal, rms and status (it does not work, if you use root instead of root4star)
  mTables->getPedestal(BTOW, TowerId, 0, pedestal, rms); //softID
  mTables->getStatus(BTOW, TowerId, status); //softID

  //Initialization of the tower coefficients
  Double_t *TowerCoeff;

  //Tower coefficients for the different runs, parameters are saved in BemcNewCalib.h
  if(picoDst->event()->runId() <= 15094020) TowerCoeff = CPre;
  else TowerCoeff = CLowMidHigh;

  //Calculation of the calibrated energy E=C*(ADC-Pedestal)
  Double_t calibEnergy = TowerCoeff[TowerId-1]*(tower->adc() - pedestal); //softID

  return calibEnergy;

  //Function returns the calibrated energy of the tower
}
//---------------------------------------------------------------------------
Bool_t StPicoD0JetAnaMaker::IsBadEnergyRun(Int_t runID) {
    // Check if the run is in the list of BEMC bad runs

    for (UInt_t i = 0; i < sizeof(EnergyBadRunList2016)/sizeof(EnergyBadRunList2016[0]); i++) {
        if (EnergyBadRunList2016[i] == runID) {
            return true;
        }
    }
    return false;

    //Function returns true if the run is in the list of bad runs
}
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
Int_t StPicoD0JetAnaMaker::isD0PairCentrality_pt(StKaonPion const & kp, Int_t Centrality, Int_t mYear) const{
    //Check if the pair passes the cuts for D0

    //Loading the daughter particles tracks
    StPicoTrack const* kaon = picoDst->track(kp.kaonIdx());
    StPicoTrack const* pion = picoDst->track(kp.pionIdx());

    //Initialisation of the pairCuts boolean
    Bool_t pairCuts = false;

    //Recalculation of the centrality binning
    Int_t Centrality2 =  (Centrality == 8 || Centrality == 7) ? 0 :
                       (Centrality == 6) ? 1 :
                       (Centrality == 5 || Centrality == 4) ? 2 :
                       (Centrality == 3 || Centrality == 2) ? 3 :
                       (Centrality == 1 || Centrality == 0) ? 4 : -1;
    // Centr.   0-10%       10-20%         20-40%         40-60%        60-80%
    // C_ID     8,7        6               5,4             3,2         1,0
    // bin       0         1                 2               3           4


    //Recalculation of the momentum binning
    Int_t KPMom = (kp.pt() < 0.5) ? 0 : (kp.pt() < 1) ? 1 : (kp.pt() < 2) ? 2 : (kp.pt() < 3) ? 3 : (kp.pt() < 5) ? 4 : 5;

    //Check if the pair passes the particular cuts
    //Parameters are saved in StCuts.cxx
    if (mYear == 2014){
        pairCuts =  sin(kp.pointingAngle())*kp.decayLength() < mycuts::DCA_D0_cut_2014[KPMom][Centrality2] &&
                    kp.pionDca() > mycuts::pionDCA_cut_2014[KPMom][Centrality2] && kp.kaonDca() > mycuts::kaonDCA_cut_2014[KPMom][Centrality2] &&
                    kp.dcaDaughters() < mycuts::pionkaonDCA_cut_2014[KPMom][Centrality2] && kp.decayLength()> mycuts::D0_decayLength_cut_2014[KPMom][Centrality2] &&
                    cos(kp.pointingAngle()) > mycuts::cosTheta_2014;
    } else if (mYear == 2016){
        pairCuts =  sin(kp.pointingAngle())*kp.decayLength() < mycuts::DCA_D0_cut_2016[KPMom][Centrality2] &&
                    kp.pionDca() > mycuts::pionDCA_cut_2016[KPMom][Centrality2] && kp.kaonDca() > mycuts::kaonDCA_cut_2016[KPMom][Centrality2] &&
                    kp.dcaDaughters() < mycuts::pionkaonDCA_cut_2016[KPMom][Centrality2] && kp.decayLength()> mycuts::D0_decayLength_cut_2016[KPMom][Centrality2] &&
                    cos(kp.pointingAngle()) > mycuts::cosTheta_2016;
    }

    //Calculation of the product of the daughter charges
    Int_t charge = kaon->charge() * pion->charge();

    //If the daughter charges are the same, it changes charge to 1 (K+) or 2 (K-)
    if(charge>0) charge = kaon->charge()>0 ? 1:2;

    //If the pair passes the cuts, it returns the charge, else 0.
    if(pairCuts) return charge;
    else return 0;

    //Function returns:
    //0  - the pair does not pass the cuts
    //-1 - the unlike-sign pair passes the cuts
    //1  - the like-sign pair passes the cuts (pi+K+)
    //2  - the like-sign pair passes the cuts (pi-K-)
}
//-----------------------------------------------------------------------------
Bool_t StPicoD0JetAnaMaker::isGoodEvent(Int_t mYear,  TH1D* hEventsCuts){
  //Check if the event passes the cuts (MB triggers, vr, vz, vzVpdVz)

  //Loading event
  StPicoEvent *event = (StPicoEvent *)picoDst->event();



  //Checking triggers
  if (!isMBTrigger(mYear)) return false;

  Double_t ShiftVz = 0;
  if (mYear == 2014) ShiftVz = 0; //2.1486; //It could be used but due to HFT we will not shift it

  //hEventsCuts: Triggers
  hEventsCuts->Fill(1);
  //Checking vr = sqrt(vx^2+vy^2)
  if (!(sqrt(event->primaryVertex().x()*event->primaryVertex().x()+event->primaryVertex().y()*event->primaryVertex().y()) < mycuts::vr)) return false;
  //hEventsCuts: v_r
  hEventsCuts->Fill(2);
  //Checking vz
  if (!(fabs(event->primaryVertex().z()-ShiftVz) < mycuts::vz)) return false;
  //hEventsCuts: v_z
  hEventsCuts->Fill(3);
  //Checking vzVpdVz
  if (!(fabs(event->primaryVertex().z() - event->vzVpd()) < mycuts::vzVpdVz)) return false;
  //hEventsCuts: |v_z - v_z_vpd|
  hEventsCuts->Fill(4);
  
  //Check on suspicious all-0 position
  Bool_t nonezeroVertex = (event->primaryVertex().x()!=0 && event->primaryVertex().y()!=0 && event->primaryVertex().z()!=0);
  if (!nonezeroVertex) return false;
  hEventsCuts->Fill(5);

  return true;

  //Function returns true if the event passes the cuts
}
//-----------------------------------------------------------------------------
Bool_t StPicoD0JetAnaMaker::isMBTrigger(Int_t mYear){
 //Function checks if the event is triggered by the MB trigger

 //Initialization of the set of triggers
 const std::set<Int_t>* mbTriggers = nullptr;

 //Different triggers for different years, saved in StCuts.cxx
 if(mYear ==2016) mbTriggers = &mycuts::mbTriggers2016;
 if(mYear ==2014) mbTriggers = &mycuts::mbTriggers2014;

 //Loading event and checking if it is triggered by the MB trigger
 StPicoEvent* event = static_cast<StPicoEvent*>(mPicoDstMaker->picoDst()->event());
 return std::any_of(mbTriggers->begin(), mbTriggers->end(), [&](Int_t trigger) { return event->isTrigger(trigger); });

 //Function returns true if the event is triggered by the MB trigger
}
//-----------------------------------------------------------------------------
Bool_t StPicoD0JetAnaMaker::isGoodTrack(StPicoTrack const * const trk) const{
    // Require at least one hit on every layer of PXL and IST.
    // It is done here for tests on the preview II data.
    // The new StPicoTrack which is used in official production has a method to check this

    //Check if the track meets the HFT requirement
    //2014 - Require at least one hit on every layer of PXL and IST
    //2016 - Require at least one hit on every layer of PXL and (IST or SST)
    //Both can be written in te same way
    Bool_t HFTCondition = (trk->hasPxl1Hit() && trk->hasPxl2Hit()) && (trk->hasSstHit() || trk->hasIstHit());

    //Check if |eta| < 1
    Bool_t EtaCondition = abs(trk->gMom().PseudoRapidity()) < 1;

    //In StCuts.cxx is defined if the HFT is required and the nHitsFit and minPt values.
    return (trk->gPt() > mycuts::minPt && trk->nHitsFit() >= mycuts::nHitsFit && (HFTCondition || !mycuts::requireHFT) && EtaCondition);



    //Function returns true if track is good
}
//-----------------------------------------------------------------------------
Bool_t StPicoD0JetAnaMaker::isGoodJetTrack(StPicoTrack const * const trk,StPicoEvent const *const myEvent) const{
    //Check if the track meets the jet track cuts
    //Parameters saved in StCuts.cxx

    //pT range cut
    Bool_t pTTrackJetCut = trk->gPt() > mycuts::jetTrackPtMin && trk->gPt() < mycuts::jetTrackPtMax;
    //eta range cut
    Bool_t etaTrackJetCut = fabs(trk->gMom().PseudoRapidity()) < mycuts::jetTrackEta;
    //nHitsFit cut
    Bool_t nHitsTrackJetCut = trk->nHitsFit() >= mycuts::jetTracknHitsFit;
    //nHitsRatio cut
    Bool_t nHitsRatioTrackJetCut = (1.0*trk->nHitsFit()/trk->nHitsMax())>=mycuts::jetTracknHitsRatio;
    //DCA cut
    Bool_t dcaTrackJetCut = fabs(trk->gDCA(myEvent->primaryVertex().x(),myEvent->primaryVertex().y(),myEvent->primaryVertex().z())) < mycuts::jetTrackDCA;

    return pTTrackJetCut && etaTrackJetCut && nHitsTrackJetCut && nHitsRatioTrackJetCut && dcaTrackJetCut;

    //Return true if all the cuts are passed
}
//-----------------------------------------------------------------------------
Bool_t StPicoD0JetAnaMaker::isTpcPion(StPicoTrack const * const trk) const{
    //Check if the track meets nsigma (TPC) requirement for pion
    //Parameters are saved in StCuts.cxx

    return fabs(trk->nSigmaPion()) < mycuts::nSigmaPion;

    //Function returns true if track is good
}
//-----------------------------------------------------------------------------
Bool_t StPicoD0JetAnaMaker::isTpcKaon(StPicoTrack const * const trk) const{
    //Check if the track meets nsigma (TPC) requirement for kaon
    //Parameters are saved in StCuts.cxx

    return fabs(trk->nSigmaKaon()) < mycuts::nSigmaKaon; //In D0 event maker it is set to 2

    //Function returns true if track is good
}
//-----------------------------------------------------------------------------
Float_t StPicoD0JetAnaMaker::getTofBeta(StPicoTrack const * const trk, StThreeVectorF const* const pVtx,StPicoDst const* const picoDst) const{
    //Calculation of beta for the track

    //index2tof is index of the track in the StPicoBTofPidTraits array
    Int_t index2tof = trk->bTofPidTraitsIndex();

    //Initialization of beta
    Float_t beta = std::numeric_limits<Float_t>::quiet_NaN();

    //If index2tof is positive, than the track has a match in the TOF
    if(index2tof >= 0){

        //Getting the pointer to the StPicoBTofPidTraits object
        StPicoBTofPidTraits *tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);

        //If pointer is not null, than we can calculate beta
        if(tofPid){

            //Calculation of beta
            beta = tofPid->btofBeta();

            //If for some reason beta is negative, than we can try to calculate beta using the pathlength and tof
            if (beta < 1e-4){

                //Getting the hit position in the TOF
                StThreeVectorF const btofHitPos = StThreeVectorF(tofPid->btofHitPos().x(),tofPid->btofHitPos().y(),tofPid->btofHitPos().z());

                //Getting the helix of the track
                StPicoPhysicalHelix helix = trk->helix(picoDst->event()->bField());

                //Calculation of the pathlength
                Float_t L = tofPathLength(pVtx, &btofHitPos, helix.curvature());

                //Calculation of the time of flight
                Float_t tof = tofPid->btof();

                //Calculation of beta for positive values of tof
                if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
                //...else beta is not defined
                else beta = std::numeric_limits<Float_t>::quiet_NaN();

            } //End of beta < 1e-4

        } //End of tofPid != NULL

    } //End of index2tof >= 0

    return beta;

    //Function returns beta of the track
}
//-----------------------------------------------------------------------------
Bool_t StPicoD0JetAnaMaker::isTofKaon(StPicoTrack const * const trk, Float_t beta) const{
    //Check if the track meets |1/beta-1/beta_K| (TOF) requirement for kaon

    //Initialization of tofKaon
    Bool_t tofKaon = false;

    //If beta is positive, than we can calculate |1/beta-1/beta_K|
    if(beta>0){

        //Calculation of the global total momentum
        Double_t ptot = trk->gMom().Mag();

        //Calculation of the expected beta for kaons
        Float_t beta_k = ptot/sqrt(ptot*ptot+M_KAON_PLUS*M_KAON_PLUS);

        //Check if the track meets the TOF requirement
        //Parameters are saved in StCuts.cxx
        tofKaon = fabs(1/beta - 1/beta_k) < mycuts::kTofBetaDiff ? true : false;

    } //End of beta > 0

    return tofKaon;

    //Function returns true if track is good based on TOF information
}
//-----------------------------------------------------------------------------
Bool_t StPicoD0JetAnaMaker::isTofPion(StPicoTrack const * const trk, Float_t beta) const{
    //Check if the track meets |1/beta-1/beta_pi| (TOF) requirement for pion

    //Initialization of tofPion
    Bool_t tofPion = false;
    
    //If beta is positive, than we can calculate |1/beta-1/beta_K|
    if(beta>0){

        //Calculation of the global total momentum
        Double_t ptot = trk->gMom().Mag();
        //Calculation of the expected beta for kaons
        Float_t beta_pi = ptot/sqrt(ptot*ptot+M_PION_PLUS*M_PION_PLUS);
        //Check if the track meets the TOF requirement
        //Parameters are saved in StCuts.cxx
        tofPion = fabs(1/beta - 1/beta_pi) < mycuts::pTofBetaDiff ? true : false;
    } //End of beta > 0

    return tofPion;

    //Function returns true if track is good based on TOF information
}





