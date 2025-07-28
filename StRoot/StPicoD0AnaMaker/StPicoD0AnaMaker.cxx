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
#include "StPicoD0AnaMaker.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StCuts.h"
#include "JetInfo.h"
#include "../StPicoPrescales/StPicoPrescales.h"
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
#include "fastjet/contrib/IterativeConstituentSubtractor.hh" 
#include "fastjet/contrib/RescalingClasses.hh"
#include "fastjet/contrib/GenericSubtractor.hh"
#include "fastjet/contrib/ShapeWithPartition.hh"
//#include "/gpfs01/star/pwg/lomicond/Ondrej/Jets/Alex_install/install/fastjet-contrib-install/include/fastjet/contrib/SignalFreeBackgroundEstimator.hh"
using namespace std;
using namespace fastjet;
namespace fj = fastjet;
//-------------------------------------

ClassImp(StPicoD0AnaMaker)

StPicoD0AnaMaker::StPicoD0AnaMaker(char const * name, char const * inputFilesList, char const * outName, StPicoDstMaker* picoDstMaker, StRefMultCorr* grefmultCorrUtil, int pYear, bool Sim, char const *filename)
        : StMaker(name),
          mOutFileName(outName), // Initialize in the same order as declared
          mChain(NULL),
          mOutputFile(NULL),
          mGRefMultCorrUtil(grefmultCorrUtil),
          mPicoDstMaker(picoDstMaker),
          mPicoD0Event(NULL),
          mInputFileList(inputFilesList),
          mEventCounter(0),
          mYear(pYear),
          picoDst(NULL),
          mPrescales(NULL),
          mHFCuts(NULL),
          mSimulation(Sim) {
          
            fMCFileListName = filename;
          }
///////////////////////////////////////////////////////////////////////////////////////////////////////
//--------Event-plane---------------
Int_t StPicoD0AnaMaker::EP_IsGoodTrack(StPicoTrack *trk, TVector3 pVertex){

	enum Position {West = -1, East = 1};
	Position Side = West;
	
	TVector3 mTrkMom;
	
	//Primary tracks
   	mTrkMom = trk->pMom(); //If not exists, it equals 0, 0, 0

	//0.2 < pt < 2 GeV
	float pt = mTrkMom.Perp();
	if (pt!=pt||pt <= 0.2 || pt >= 2.0) return 0;

	//0.05 < |eta| < 1.00 (West/East)
	float eta = mTrkMom.PseudoRapidity();
	if (abs(eta) <= 0.05 || abs(eta) >= 1) return 0;
	if (eta > 0) Side = East;

	//nHitsFit > 15
	float nHitsFit = trk->nHitsFit();
	float nHitsMax = trk->nHitsMax();
	
	//nHitsFit/nHitsMax => 0.52
	float nHitsRatio = 1.0*nHitsFit/nHitsMax;
	if (nHitsFit < 15 || nHitsRatio < 0.52) return 0;
	
	
	//DCA < 1 cm
	float dca = trk->gDCA(pVertex).Mag();
	if (dca >= 1) return 0;

return Side;
}

void StPicoD0AnaMaker::CalculateEventPlane(){



    //Primary vertex
    TVector3 prVertex = picoDst->event()->primaryVertex();
    //Q-vectors
    double Q_1 = 0;
    double Q_2 = 0;
    
 
    

    for (UInt_t iTrack = 0; iTrack < picoDst->numberOfTracks(); iTrack++){
      
        StPicoTrack *trk = static_cast<StPicoTrack *>(picoDst->track(iTrack));
        if (!trk) continue;
        
        double Goodtrack = EP_IsGoodTrack(trk,prVertex);
        if (!abs(Goodtrack)) continue;
        
        TVector3 trackMom = trk->pMom();
        
        double pPt = trackMom.Perp();
        double phi = trackMom.Phi();
        
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
    double Q_1rc = -0.5976;
    double Q_2rc = 1.842;
    
    double Q_1corr = Q_1 - Q_1rc;
    double Q_2corr = Q_2 - Q_2rc;
    
    fQ_1_rec = Q_1corr;
    fQ_2_rec = Q_2corr;
    //-------------------
    //7.2.2025

    double Psi_2 = 1./2 * TMath::ATan2(Q_2corr, Q_1corr);
   
    //-------------------
    //Psi shift	
    //You have to run the code for "all" events to get the mean value of Q vector (again)
    std::vector<double> A_2 = {-0.0471131, -0.0114311, -0.000157011, 0.00128572, -0.000253291, 0.000446538, 0.00223142, 0.000574021, 0.00133264, 0.00157648, 0.00169576, -0.00023281, 0.00174849, 0.00194597, 0.000343813, 0.00104523, 0.00194369, -0.000649944, 0.00073013, -0.00138147, -0.000686723};
    std::vector<double> B_2 = {0.0204715, -0.0151115, -0.00306008, -0.000966112, 0.000103265, -0.0014271, -0.00244014, 0.000801171, -0.00235772, 0.00173085, -8.08274e-05, -0.000880916, 0.000967381, 0.00182509, 0.00167076, -6.81227e-05, 8.50044e-05, -0.000525437, -0.000912995, 0.000685826, -0.000467597};

    double CorrectedPsi2 = Psi_2;
    for (int i = 1; i <= 21; i++){
         CorrectedPsi2 +=(1.0 / 2) * (2.0 / i) *(-A_2[i - 1] * cos(2 * i * Psi_2) + B_2[i - 1] * sin(2 * i * Psi_2));
    }
    
    if (Psi_2!=Psi_2 || CorrectedPsi2!=CorrectedPsi2) return;
    
    fPsi_2 = Psi_2;
    
    //Force the range (-pi/2,pi/2)
    CorrectedPsi2 = TMath::ATan2(TMath::Sin(2 * CorrectedPsi2), TMath::Cos(2 * CorrectedPsi2)) / 2.;
    Psi_2 = TMath::ATan2(TMath::Sin(2 * Psi_2), TMath::Cos(2 * Psi_2)) / 2.;

    
    //-------------------
    fPsi_2_shifted = CorrectedPsi2;
   



return;
}
//////
std::vector<fj::PseudoJet> JetReconstructionICS(vector<fj::PseudoJet> fInputVectors, int fCentrality, double fR, bool fBackSub, double fEP_psi2, int difiter){


int NiterA[4] = {4,3,2,2};
double R_max1A[4] = {0.05,0.125,0.100,0.150};
double R_max2A[4] = {0.005,0.005,0.175,0.100};
bool fMassiveTest = true;
bool fPhiModulation = true;
/*****/
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
    
	//--- Background Estimation ---
	fj::JetDefinition jet_def_bkgd(fj::kt_algorithm, fR, fj::E_scheme, fj::Best); // Define jet algorithm for background estimation
	fj::AreaDefinition area_def_bkgd(fj::active_area_explicit_ghosts, fj::GhostedAreaSpec(1., 1, 0.01)); // Define the area for background estimation

	// Decide how many jets to remove based on centrality
	int nJetsRemove = 1;
	if (fCentrality == 7 || fCentrality == 8) nJetsRemove = 2;

	// Selector to choose the hardest jets based on centrality and eta
	fj::Selector selector = (!fj::SelectorNHardest(nJetsRemove)) * fj::SelectorAbsEtaMax(0.6);

	// Define seeds for random number generator (used for area calculation)
	/*
	unsigned int seed1 = 12345;
	unsigned int seed2 = 56789;
	std::vector<int> seeds = { static_cast<int>(seed1), static_cast<int>(seed2) };
	*/

	// Define the area for background estimation with fixed seeds
	////fastjet::AreaDefinition area_def_bkgd = area_def_bkgd.with_fixed_seed(seeds);

	// Create background estimator using the previously defined selector and jet algorithm
	fj::JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);

	// Set the maximum eta value for background estimation
	double max_eta = 1;
        
	//--- Jet Reconstruction ---
	fj::contrib::IterativeConstituentSubtractor subtractor; // Set up the background subtraction algorithm
	subtractor.set_distance_type(fj::contrib::ConstituentSubtractor::deltaR); // Set distance type for subtraction

  	// Set parameters for the background subtraction (maximum distance and alpha values)
	vector<double> max_distances;
	vector<double> alphas;
	
	//For higher iterations always used the same R_1max and alpha
	if(NiterA[difiter] > 2){
	
		for (int it = 0; it < NiterA[difiter] ; it++){
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
		cout << NiterA[difiter] << " je tp 11?" << endl;
	for (auto& a : alphas) {
	
	cout << "alpha: " << a << endl;
	}
	for (auto& rr : max_distances) {
	
	cout << "max_distances: " << rr << endl;
	}
}

	//Parama
	///max_distances.push_back(frMax2);





		
	//Parama
	//alphas.push_back(fAlpha2);
	
	
	//exclude D0
	fj::Selector notD02 = fastjet::SelectorMassMax(1); // 1 GeV max mass
	std::vector<fastjet::PseudoJet> particles_without_D0;
	std::vector<fastjet::PseudoJet> D0_pseudojet;
	std::vector<fastjet::PseudoJet> all_vectors;

	
	for (auto& p : fInputVectors) {
	////if ( abs(p.user_index()) > 9999) continue;
	
	    if (notD02(p)) {
	    	////fastjet::PseudoJet temp_jet(p.px(), p.py(), p.pz(), sqrt(p.px()*p.px()+p.py()*p.py()+p.pz()*p.pz()));
	    	
	    	fastjet::PseudoJet temp_jet; 


	    	if (fMassiveTest) temp_jet = fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.E());
	    	else temp_jet = fastjet::PseudoJet(p.px(), p.py(), p.pz(), sqrt(p.px()*p.px()+p.py()*p.py()+p.pz()*p.pz()));
	    	
	    	
		temp_jet.set_user_index(p.user_index());
	    	particles_without_D0.push_back(temp_jet);
		all_vectors.push_back(temp_jet);
	    } else {
	    fastjet::PseudoJet temp_jet; 
	    ////fastjet::PseudoJet temp_jet(p.px(), p.py(), p.pz(), sqrt(p.px()*p.px()+p.py()*p.py()+p.pz()*p.pz()));
	    if (fMassiveTest) temp_jet = fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.E());
	    else temp_jet = fastjet::PseudoJet(p.px(), p.py(), p.pz(), sqrt(p.px()*p.px()+p.py()*p.py()+p.pz()*p.pz()));
	    temp_jet.set_user_index(p.user_index());
	    all_vectors.push_back(temp_jet);
	    D0_pseudojet.push_back(temp_jet);
	    }
	}
	///////////////////////////////////
	///////////////////////////////////
	///////////////////////////////////
	/*
	 for (vector<fastjet::PseudoJet>::const_iterator particle = D0_pseudojet.begin(); particle != D0_pseudojet.end(); ++particle) {
	                ////int index = particle->user_index();
                cerr << "#  rap   phi   perp   m   index" << endl;
		printf("%5u %15.8f %15.8f %15.8f %15.8f %5u\n", 0, particle->rap(), particle->phi(), particle->perp(), particle->m(), particle->user_index());
	 }
	*/
	// Apply the subtraction parameters
	subtractor.set_parameters(max_distances, alphas);
	subtractor.set_ghost_removal(true); // Enable ghost removal
	subtractor.set_ghost_area(0.005); // Set ghost area value
	subtractor.set_max_eta(max_eta); // Set maximum eta for particles
	
	
	
	subtractor.set_background_estimator(&bkgd_estimator); // Link the background estimator
	
	if(fMassiveTest) subtractor.set_common_bge_for_rho_and_rhom(true); // Set common background estimation for rho and rhom
	if(fMassiveTest) subtractor.set_keep_original_masses(); // Keep the original masses of particles
	
        
       	// Selector for particles with mass below 1 GeV (likely excludes certain particles)
	fj::Selector notD0 = fastjet::SelectorMassMax(1); // 1 GeV max mass
	subtractor.set_particle_selector(&notD0);  

	subtractor.initialize(); // Initialize the background subtraction algorithm

	// Define jet algorithm for actual jet reconstruction (Anti-kt with radius fR)
	fj::JetDefinition jet_def(fj::antikt_algorithm, fR, fj::E_scheme, fj::Best);

	// Define area for jet reconstruction
	fj::AreaDefinition area_def_jet(fj::active_area, fj::GhostedAreaSpec(1., 1, 0.01));
	////fastjet::AreaDefinition area_def_jet = area_def_jet2.with_fixed_seed(seeds);

	// Set particles for background estimator (background estimation for the input particles)
	////bkgd_estimator.set_particles(fInputVectors); //all_vectors
	//Rescaling
	if (fPhiModulation) bkgd_estimator.set_rescaling_class(&rescaling);

	bkgd_estimator.set_particles(all_vectors);

	
	// Minimum pt for jets (below this, jets are excluded)
	double ptmin = -0.01;
	vector<fj::PseudoJet> corrected_event;

	// Apply background subtraction if enabled
	////if (fBackSub) corrected_event = subtractor.subtract_event(fInputVectors); 
	////else corrected_event = fInputVectors; // Otherwise, use original input vectors
	if (fBackSub) corrected_event = subtractor.subtract_event(particles_without_D0); 
	else corrected_event = particles_without_D0; // Otherwise, use original input vectors
	
	//Reutrn back D0
	corrected_event.push_back(D0_pseudojet.back());
	
	// Perform jet clustering with the background-subtracted (or original) particles
	ClusterSequenceArea *fClustSeq = new fastjet::ClusterSequenceArea(corrected_event, jet_def, area_def_jet);
	
	// Sort the jets by transverse momentum (pt)
	vector<fastjet::PseudoJet> fInclusiveJets;
	fInclusiveJets.clear();
	fInclusiveJets = sorted_by_pt(fClustSeq->inclusive_jets(ptmin));
/*****/
	// Return the list of inclusive jets
	return fInclusiveJets;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
Int_t StPicoD0AnaMaker::Init()
{

  if (mSimulation){
  cout << "***************************************" << endl;
  cout << "**                                   **" << endl;
  cout << "**            SIMULATION             **" << endl;
  cout << "**                                   **" << endl;
  cout << "***************************************" << endl;


  }


  mPicoD0Event = new StPicoD0Event();

  mChain = new TChain("T");
  std::ifstream listOfFiles(mInputFileList.Data());
  if (listOfFiles.is_open())
  {
    std::string file;
    while (getline(listOfFiles, file))
    {
      LOG_INFO << "StPicoD0AnaMaker - Adding :" << file <<endm;
      mChain->Add(file.c_str());
      LOG_INFO<<" Entries = "<<mChain->GetEntries()<< endm; 
    }
  }
  else
  {
    LOG_ERROR << "StPicoD0AnaMaker - Could not open list of files. ABORT!" << endm;
    return kStErr;
  }
  

  //The prescales are used only for rescaling the RunID to the lower numbers
  mPrescales = new StPicoPrescales(mycuts::prescalesFilesDirectoryName);

  //Loading the information from the D0EventMaker
  mChain->GetBranch("dEvent")->SetAutoDelete(kFALSE);
  mChain->SetBranchAddress("dEvent", &mPicoD0Event);

  //Output file
  mOutputFile = new TFile(mOutFileName.Data(), "RECREATE");
  mOutputFile->cd();

  //Number of Runs
  int nRuns = mPrescales->numberOfRuns();

  //Events histograms:
  vtxz = new TH1F("vtxz",";PVtx.z() [cm]; Count",100,-10,10);
  vtxr = new TH2D("vtxr",";PVtx.x() [cm]; PVtx.y() [cm]",100,-3,3,100,-3,3);
  hcentr = new TH1D("hcentr",";C_{ID}",9,-0.5,8.5);
  hcentrW = new TH1D("hcentrW",";C_{ID}",9,-0.5,8.5);
  NEventsCuts = new TH1F("NEventsCuts", "NEventsCuts;Cuts;Count", 10, 0, 10);
      NEventsCuts->GetXaxis()->SetBinLabel(1, "All events");
      NEventsCuts->GetXaxis()->SetBinLabel(2, "Triggers");
      NEventsCuts->GetXaxis()->SetBinLabel(3, "V_{r}");
      NEventsCuts->GetXaxis()->SetBinLabel(4, "V_{z}");
      NEventsCuts->GetXaxis()->SetBinLabel(5, "|V_{z} - V_{z}^{VPD}|");
      NEventsCuts->GetXaxis()->SetBinLabel(6, "|V_{x,y,z}|!=0");
      NEventsCuts->GetXaxis()->SetBinLabel(7, "Good E run");
      NEventsCuts->GetXaxis()->SetBinLabel(8, "Centrality");
      NEventsCuts->GetXaxis()->SetBinLabel(9, "Good D^{0}");
      NEventsCuts->GetXaxis()->SetBinLabel(10, "Cal.: E_{T} < 30 GeV");
  mh1TotalEventsInRun = new TH1F("mh1TotalEventsInRun","totalEventsInRun;runIndex;totalEventsInRun",nRuns+1,0,nRuns+1);
  mh1TotalGRefMultInRun = new TH1F("mh1TotalGRefMultInRun","totalGRefMultInRun;runIndex;totalGRefMultInRun",nRuns+1,0,nRuns+1);

  //D0 histograms:
  D0etalike = new TH1D("D0etalike",";#eta; Count",100,-5,5);
  D0etaunlike = new TH1D("D0etaunlike",";#eta; Count",100,-5,5);
  pioneta = new TH1D("pioneta",";#eta; Count",100,-5,5);
  kaoneta =  new TH1D("kaoneta",";#eta; Count",100,-5,5);
  mh1TotalKaonsInRun = new TH1F("mh1TotalKaonsInRun","totalKaonsInRun;runIndex;totalKaonsInRun",nRuns+1,0,nRuns+1);
  mh1TotalPionsInRun = new TH1F("mh1TotalPionsInRun","totalPionsInRun;runIndex;totalPionsInRun",nRuns+1,0,nRuns+1);
  mh1TotalD0CandidatesInRun = new TH1F("mh1TotalD0CandidatesInRun","totalD0CandidatesInRun;runIndex;totalD0CandidatesInRun",nRuns+1,0,nRuns+1);
  mh2NKaonsVsNPions = new TH2F("mh2NKaonsVsNPions","nKaonsVsNPions;nPions;nKaons",1000,0,1000,300,0,300);
  mh2KaonDcaVsPt = new TH2F("mh2KaonDcaVsPt","kaonDcaVsPt;p_{T}(K#pi)(GeV/c);K DCA(cm)",120,0,12,50,0,0.05);
  mh2PionDcaVsPt = new TH2F("mh2PionDcaVsPt","pionDcaVsPt;p_{T}(K#pi)(GeV/c);#pi DCA(cm)",120,0,12,50,0,0.05);
  mh2CosThetaVsPt = new TH2F("mh2CosThetaVsPt","cosThetaVsPt;p_{T}(K#pi)(GeV/c);cos(#theta)",120,0,12,500,0,1.0);
  mh2DcaDaughtersVsPt = new TH2F("mh2DcaDaughtersVsPt","dcaDaughtersVsPt;p_{T}(K#pi)(GeV/c);dcaDaughters(cm)",120,0,12,200,0,0.02);

  //Jet tracks
  JetTracksdEdx = new TH2D("JetTracksdEdx","JetTracksdEdx;charge #times |p| [GeV/c]; dE/dx [KeV/cm]",200,-4,4,200,0,8);
  JetTracksdEdxCut = new TH2D("JetTracksdEdxCut","JetTracksdEdxCut;charge #times |p| [GeV/c]; dE/dx [KeV/cm]",200,-4,4,200,0,8);
  hpT_tr = new TH1F("hpT_tr","hpT_tr;p_{T}(GeV/c);",200,0,100);
  heta_phi_tr = new TH2F("heta_phi_tr","heta_phi_tr;#phi;#eta",120,0,12,50,0,0.05);
  heta_tr = new TH1F("heta_tr","heta_tr;#eta",50,-6,6);
  hphi_tr = new TH1F("hphi_tr","hphi_tr;#phi;",40,-6.2830,6.2830);
  hdca_z_tr = new TH2F("hdca_z_tr","hdca_z_tr;DCA(cm); v_{z}(GeV/c)",120,0,12,50,-10,10);
  hdca_pT = new TH2F("hdca_pT","hdca_pT; DCA(cm); p_{T}(GeV/c) ",120,0,12,100,0,200);
  hdca_tr = new TH1F("hdca_tr","hdca_tr; DCA(cm)",120,0,12);
  hcharged_tr = new TH1F("hcharged_tr","hcharged_tr;Charge;",5,-1,1);
  
  /*for (int cent = 0; cent < 9; cent++){
	  for (int ptbin = 0; ptbin < 6; ptbin++){
	  	for (int ch = 0; ch <= 1; ch++){
	  		hDcaSign[cent][ptbin][ch] = new TH1D(Form("hDcaSign_%i_%i_%i", cent, ptbin, ch), Form("hDcaSign_%i_%i_%i", cent, ptbin, ch), 200, -5, 5);
	  	}
	  }
  }
  */

  //Jet background
  Jet_grefmult_pt_background = new TH2D("Jet_rho_vs_grefmult","Jet_rho_vs_grefmult;grefmult; p_{T} (GeV/c)",1000,0,1000,100,-1,32);
  Jet_D0pT_vs_D0rapidity = new TH2D("Jet_D0pT_vs_D0rapidity","Jet_D0pT_vs_D0rapidity;rapidity; p_{T} (GeV/c)",100,-1.5,1.5,100,0,10);
  Jet_D0pT_vs_Jetrapidity = new TH2D("Jet_D0pT_vs_JetRapidity","Jet_D0pT_vs_JetRapidity;rapidity; p_{T} (GeV/c)",100,-1.5,1.5,100,0,10);
  Jet_phi = new TH1D("Jet_phi","Jet_phi;#phi;",40,-6.2830,6.2830);


  //TNtuple D0-jets
  for (int i = 0; i < 4; i++){
  Jets[i] = new TNtuple(Form("Jets_%i",i), Form("Jets_%i",i), "RunId:centrality:centr_weight:NJet:jet_eta:jet_phi:grefmult:bg_dens:jet_area:jet_rap:jet_pt:jet_pt_corr:D0mass:D0_r:D0_pT:lambda_1_0half:lambda_1_1:lambda_1_1half:lambda_1_2:lambda_1_3:lambda_2_0:z:NConst:NpTfraction:D0_rap:D0_eta:Psi2");
}
VariableJets = {
    {"RunId", 0},          // sgn(RunID) = 1 -> D0, sgn(RunID) = -1 -> antiD0
    {"centrality", 1},     // 0 -> 70-80%, 1 -> 60-70%, 2 -> 50-60%, 3 -> 40-50%, 4 -> 30-40%, 5 -> 20-30%, 6 -> 10-20%, 7 -> 5-10%, 8 -> 0-5%
    {"centr_weight", 2},   // centrality weight
    {"NJet", 3},           // Number of D0 in one event
    {"jet_eta", 4},        // pseudorapidity of D0-jet
    {"jet_phi", 5},        // phi of D0-jet
    {"grefmult", 6},       // grefmult
    {"bg_dens", 7},        // background density
    {"jet_area", 8},       // area of D0-jet
    {"jet_rap", 9},        // rapidity of D0-jet
    {"jet_pt", 10},        // pt of D0-jet
    {"jet_pt_corr", 11},   // pt of D0-jet after background subtraction
    {"D0mass", 12},        // mass of D0
    {"D0_r", 13},          // r of D0
    {"D0_pT", 14},         // pt of D0
    {"lambda_1_0half", 15},// angularity kappa=1 and alpha=0.5 after background subtraction
    {"lambda_1_1", 16},    // angularity kappa=1 and alpha=1 after background subtraction
    {"lambda_1_1half", 17},// angularity kappa=1 and alpha=1.5 after background subtraction
    {"lambda_1_2", 18},    // angularity kappa=1 and alpha=2 after background subtraction
    {"lambda_1_3", 19},    // angularity kappa=1 and alpha=3 after background subtraction
    {"lambda_2_0", 20},    // angularity kappa=2 and alpha=0 after background subtraction
    {"z", 21},             // z of D0-jet after background subtraction
    {"NConst", 22},        // Number of constituents
    {"NpTfraction", 23},   // Neutral pT fraction of D0-jet
    {"D0_rap", 24},   	   // D0 meson rapidity
    {"D0_eta", 25},   	   // D0 meson pseudorapidity
    {"Psi2", 26}

};

////EventStats = new TNtuple("EventStats", "EventStats","RunId:centrality:centr_weight");

  //Bin size of 2D D0 mass-pt like-sign and unlike-sign histograms
  const int xbinSize=100;
  float binMass[2001];

  //Calculating of bin edges
  float xbin[101];
  for(int i=0;i<101;i++) xbin[i] = 0.1*i;
  for(int i=0;i<2001;i++) binMass[i] = 0.01*i;

  //2D D0 mass-pt like-sign and unlike-sign histograms
  massPt = new TH2D("massPt",";M_{K#pi} [GeV/c^{2}]; p_{T} [GeV/c]",2000,binMass,xbinSize,xbin);
  massPtLike = new TH2D("massPtLike",";M_{K#pi} [GeV/c^{2}]; p_{T} [GeV/c]",2000,binMass,xbinSize,xbin);

  //Loading of BEMC tables
  ////StMaker* maker = GetMaker("Eread");
  mADCtoEMaker = dynamic_cast<StEmcADCtoEMaker*>(GetMaker("Eread"));
  assert(mADCtoEMaker);
  mTables = mADCtoEMaker->getBemcData()->getTables();
  
  ifstream filelistforMCEvents(fMCFileListName.Data());
  
  if(mSimulation){
  
	  string line;
	  while (getline(filelistforMCEvents, line))
	  {
	    TString s(line);
	    filenamesforHIOverlay.push_back(s);
	  }
  
  } else {
  
  
  }



  return kStOK;
}

//-----------------------------------------------------------------------------
StPicoD0AnaMaker::~StPicoD0AnaMaker(){
  //Destructor
  delete mGRefMultCorrUtil;
}

//-----------------------------------------------------------------------------
struct FourMomentum {
  //Four-momentum of the reconstructed particle
  double E, px, py, pz, D0_antiD0, D0Mass;
  // D0_antiD0: 1 -> D0, -1 -> antiD0
};

//-----------------------------------------------------------------------------
double delta_R(double eta1, double phi1, double eta2, double phi2) {
  //Calculating of delta R = sqrt(delta eta^2 + delta phi^2)
  double deta = eta1 - eta2;
  double dphi = TVector2::Phi_mpi_pi(phi1 - phi2);

  return sqrt(deta*deta + dphi*dphi);

  //Function returns delta R
}

//-----------------------------------------------------------------------------
Int_t StPicoD0AnaMaker::Finish(){
  //Write histograms and close the output file, if you create a new histogram it has to be added here
  LOG_INFO << " StPicoD0AnaMaker - writing data and closing output file " <<endm;
  fout.close();
  mOutputFile->cd();

  //Events histograms:
  vtxz->Write();
  vtxr->Write();
  hcentr->Write();
  hcentrW->Write();
  NEventsCuts->Write();
  mh1TotalEventsInRun->Write();
  mh1TotalGRefMultInRun->Write();

  //D0 histograms:
  D0etalike->Write();
  D0etaunlike->Write();
  pioneta->Write();
  kaoneta->Write();
  mh1TotalKaonsInRun->Write();
  mh1TotalPionsInRun->Write();
  mh1TotalD0CandidatesInRun->Write();
  mh2NKaonsVsNPions->Write();
  mh2KaonDcaVsPt->Write();
  mh2PionDcaVsPt->Write();
  mh2CosThetaVsPt->Write();
  mh2DcaDaughtersVsPt->Write();

  //Jet Tracks
  JetTracksdEdx->Write();
  JetTracksdEdxCut->Write();
  hpT_tr->Write();
  heta_phi_tr->Write();
  heta_tr->Write();
  hphi_tr->Write();
  hdca_z_tr->Write();
  hdca_pT->Write();
  hdca_tr->Write();
  hcharged_tr->Write();
  /*
  for (int cent = 0; cent < 1; cent++){
	  for (int ptbin = 0; ptbin < 1; ptbin++){
	  	for (int ch = 0; ch <= 1; ch=ch++){
	  		hDcaSign[cent][ptbin][ch]->Write();
	  	}
	  }
  }*/

  //Jet background
  Jet_grefmult_pt_background->Write();
  Jet_D0pT_vs_D0rapidity->Write();
  Jet_D0pT_vs_Jetrapidity->Write();
  Jet_phi->Write();

  //2D D0 mass-pt like-sign and unlike-sign histograms
  massPt->Write();
  massPtLike->Write();

  //TNtuple D0-jets
  for (int i = 0; i <4; i++){
  	Jets[i]->Write();
  }
  ////EventStats->Write();
  
  //Closing of the output file
  mOutputFile->Close();

  //Closing of the input file including prescales
  delete mPrescales;

  return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoD0AnaMaker::Make()
{
  //Main function where the analysis is done, each "event" is analyzed here
  readNextEvent();

  //Check if everything is loaded properly
  if (!mPicoDstMaker)
  {
    LOG_WARN << " StPicoD0AnaMaker - No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  picoDst = mPicoDstMaker->picoDst();
  if (!picoDst)
  {
    LOG_WARN << "StPicoD0AnaMaker - No PicoDst! Skip! " << endm;
    return kStWarn;
  }
 

  //Check if the loaded picodsts are consistent with raw D0 reconstructed data from StPicoD0Event
  //Only if it is not simulation
  if(mPicoD0Event->runId() != picoDst->event()->runId() || mPicoD0Event->eventId() != picoDst->event()->eventId())
  {
    LOG_ERROR <<" StPicoD0AnaMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!"<<endm;
    LOG_ERROR <<" StPicoD0AnaMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoEvent and StPicoD0Event are not in sync."<<endm;
    exit(1);
  }
  

 ////if (picoDst->event()->eventId()!=2698540)return kStOK;
 ////cout << "event: " << picoDst->event()->eventId() << endl;

  // -------------- USER ANALYSIS -------------------------

  //Loading of the raw D0 daughters
  TClonesArray const * aKaonPion = mPicoD0Event->kaonPionArray();

  //Vertex of the event
  StThreeVectorF pVtx(-999.,-999.,-999.);

  //Loading of the event
  StPicoEvent *event = (StPicoEvent *)picoDst->event();

  //Loading of the primary vertex
  pVtx = StThreeVectorF(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z());

  //NEventsCuts: All events
  NEventsCuts->Fill(0);
  
  bool isBadRun = false;
  int runID = mPicoD0Event->runId();
  /*
  //Is bad run
  if (mYear == 2016){
	for (size_t i = 0; i < sizeof(EnergyBadRunList2016)/sizeof(int); ++i) {
	    if (EnergyBadRunList2016[i] == runID) {
		isBadRun = true;
		break;
	    }
	}
	for (size_t i = 0; i < sizeof(BadRunListHFT2016)/sizeof(int); ++i) {
	    if (BadRunListHFT2016[i] == runID) {
		isBadRun = true;
		break;
	    }
	}
  }
*/  
if (mYear == 2014 && mycuts::AllBadRunList2014.count(runID)) {
  isBadRun = true;
}

        
  if(isBadRun) return kStOK;
  
  //Check if the event is good (vertex, pile-up, MB trigger)
  if(!(isGoodEvent(mYear,NEventsCuts))) return kStOK;

  //Check if the GRefMultCorr information is available
  if(!mGRefMultCorrUtil) {
    LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
    return kStWarn;
  }

  //Check if the event is in a bad run (due to BEMC towers)
  if (mYear==2016 && IsBadEnergyRun(mPicoD0Event->runId())) return kStOK;

  //2014 Good run list
  
  bool isGoodRun = true;
  const std::set<int>* goodRunList = nullptr;
  /*
  if (mYear == 2014){
	int runIdLoad = mPicoD0Event->runId();
	goodRunList = &mycuts::goodRun2014; 
	isGoodRun = goodRunList->find(runIdLoad) != goodRunList->end();
  }
  
  if (!isGoodRun) return kStOK;
  */
  ///////	

  //NEventsCuts: Good E
  NEventsCuts->Fill(6);

  //Loadig of the GRefMultCorr information
  mGRefMultCorrUtil->init(picoDst->event()->runId());
  mGRefMultCorrUtil->initEvent(picoDst->event()->grefMult(),pVtx.z(),picoDst->event()->ZDCx()) ;

  //Check if the event is in a bad run (due to centrality estimation)
  if (mGRefMultCorrUtil->isBadRun(picoDst->event()->runId())) return kStOK;

  //Loading of the centrality
  int centrality  = mGRefMultCorrUtil->getCentralityBin9();

  //Check if the centrality is in the range
  if(centrality<0) return kStOK;

  //NEventsCuts: Centrality
  NEventsCuts->Fill(7);

  //Loading of the weight for the centrality correction
  double reweight = mGRefMultCorrUtil->getWeight();

  //Filling events histograms
  vtxz->Fill(pVtx.z());
  vtxr->Fill(pVtx.x(),pVtx.y());
  hcentr->Fill(centrality);
  hcentrW->Fill(centrality,reweight);
  
  ////EventStats->Fill(picoDst->event()->runId(),centrality,reweight); 

  //Loading event information
  //int eventID = mPicoD0Event->eventId();
  int RunId = mPicoD0Event->runId();

  //Loading of the number of tracks in the event
  UInt_t nTracks = picoDst->numberOfTracks();



//---------------------Sim----------------------------------------

//----------------------------------------------------------------

  //Variable checking if there is a good D0 candidate in the event, changed to true, if the candidate is found
  bool IsThereD0 = false;

  //Preparation of the vector of the daughter candidates
  std::vector<double> DaughterPionTrackVector;
  std::vector<double> DaughterKaonTrackVector;

  //Print of the Fastjet banner
  fastjet::ClusterSequence::print_banner();

  //Preparation of the four-vector of the D0 candidates
  std::vector<FourMomentum> D0_fourmomentum;



	  for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx){

	    //Loading of the raw D0 candidates and there daughters
	    StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idx);
	    StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
	    StPicoTrack const* pion = picoDst->track(kp->pionIdx());

	//----------TPC-and-TOF-identification-of-the-daughter-tracks----------
	    //Check if the daughter tracks are good.
	    if (!isGoodTrack(kaon) || !isGoodTrack(pion)) continue;

	    //Check if the daughter pions tracks are good.
	    bool tpcPion = isTpcPion(pion);

	    //Check if the daughter kaon tracks are good.
	    bool tpcKaon = isTpcKaon(kaon);

	    //Calculating of the beta for the kaon (TOF)
	    float kBeta = getTofBeta(kaon,&pVtx,picoDst);
	    
	    //Calculating of the beta for the pion (TOF)
	    float pBeta = getTofBeta(pion,&pVtx,picoDst);

	    //Check if the there is a TOF information for the kaon
	    bool kTofAvailable = kBeta>0;
	    
	    //Check if the there is a TOF information for the pion
	    bool pTofAvailable = pBeta>0; //comp. check

	    //Check if the kaon is good w.r.t the expected beta value (TOF)
	    bool tofKaon = kTofAvailable && isTofKaon(kaon,kBeta);
	    
	    //Check if the pion is good w.r.t the expected beta value (TOF)
	    bool tofPion = pTofAvailable && isTofPion(pion,pBeta); //comp. check

	    //Final check if the kaon is good. If the TOF information is not available the more strict TPC check is used.
	    bool goodKaon = (kTofAvailable && tofKaon) || (!kTofAvailable && tpcKaon);
	    
	    //Final check if the pion is good. If the TOF information is not available the more strict TPC check is used.
	    //Hybrid pion (not used)
	    bool goodPion = (pTofAvailable && tofPion) || (!pTofAvailable && tpcPion); //comp. check
	    //TPC pion
	    //bool goodPion = tpcPion;

	    //Check if the kaon is good
	    if(!goodKaon) continue;
	    

	    //Check if the pion is good
	    if(!goodPion) continue;


	//----------Pair-cuts--------------------------------------------------
	    //Initialisation of the charge
	    int charge=0;

	    //Loading of the D0 charge and mass
	    float d0Pt = kp->pt();  
	    double dMass = kp->m();

	    //Upper cut on the D0 pT
	    if(d0Pt>10) continue;

	    //Initialisation of the centrality weight
	    double reweight_eff = 1.;

	    //Check if all pair cuts conditions are met
	    if((charge=isD0PairCentrality_pt(kp,centrality, mYear))!=0 ){
		//Charge = -1 -> Unlike-sign (pi+K- or pi-K+), Charge = 1 -> Like-sign (pi+K+), Charge = 2 -> Like-sign (pi-K-)

		//If pair is Unlike-sign
		if(charge==-1){

		    //Filling of the D0 histograms
		    massPt->Fill(dMass,d0Pt,reweight*reweight_eff);
		    D0etaunlike->Fill(kp->eta());
	      
		    //Saving of the daughter tracks
		    DaughterPionTrackVector.push_back(pion->id());
		    DaughterKaonTrackVector.push_back(kaon->id());

		    //Check ff the mass is in the D0 mass window
		    //if(dMass>1.81&&dMass<1.91){

		        //The event is noted
		        IsThereD0 = true;

		        //Saving of D0 four-momenta // E,       px,       py,        pz,         D0_antiD0,            D0 mass;
		        FourMomentum D0_actual = {kp->Energy(), kp->Px(), kp->Py(),  kp->Pz(), (double)pion->charge(), dMass};
		        D0_fourmomentum.push_back(D0_actual);
		    //}

		    //Loading of the rescaled RunID
		    int runIndex = mPrescales->runIndex(mPicoD0Event->runId());

		    //Filling of the histograms
		    pioneta->Fill(pion->pMom().PseudoRapidity());
		    kaoneta->Fill(kaon->pMom().PseudoRapidity());
		    mh1TotalEventsInRun->Fill(runIndex);
		    mh1TotalGRefMultInRun->Fill(runIndex,picoDst->event()->grefMult());
		    mh1TotalKaonsInRun->Fill(runIndex,mPicoD0Event->nKaons());
		    mh1TotalPionsInRun->Fill(runIndex,mPicoD0Event->nPions());
		    mh1TotalD0CandidatesInRun->Fill(runIndex,mPicoD0Event->nKaonPion());
		    mh2NKaonsVsNPions->Fill(mPicoD0Event->nPions(),mPicoD0Event->nKaons());
		    mh2KaonDcaVsPt->Fill(kp->pt(),kp->kaonDca());
		    mh2PionDcaVsPt->Fill(kp->pt(),kp->pionDca());
		    mh2CosThetaVsPt->Fill(kp->pt(),cos(kp->pointingAngle()));
		    mh2DcaDaughtersVsPt->Fill(kp->pt(),kp->dcaDaughters());

		} //end of the Unlike-sign if

		//If pair is Like-sign
		if(charge>0){

		    //Filling of the D0 histograms
		    massPtLike->Fill(dMass,d0Pt,reweight*reweight_eff);
		    D0etalike->Fill(kp->eta());

		} //end of the Like-sign if
	       
	    } //end of the pair cuts if

	  } //end of the D0 candidate loop

//---------------------------------------------------------------------
//-------------------JET-RECONSTRUCTION-PART---------------------------
//---------------------------------------------------------------------

//fake d0 for tuning
/*
cout << "IsThereD0: " << IsThereD0 << endl;
                FourMomentum D0_actual = {1, 1, 100,  0, 1, 100};
                D0_fourmomentum.push_back(D0_actual);
            DaughterPionTrackVector.push_back(1000000);
            DaughterKaonTrackVector.push_back(1000000);
IsThereD0=1;
*/
  //Check if there is a D0 candidate in the event
  if(IsThereD0){

  CalculateEventPlane();
//cout << "psi_2: " << fPsi_2_shifted << endl;

      //NEventsCuts: Good D0 candidate
      NEventsCuts->Fill(8);

      //Initialisation of the input particle vectors for FastJet
      vector<fastjet::PseudoJet> input_particles;
      vector<fastjet::PseudoJet> chargedjetTracks;
      vector<fastjet::PseudoJet> neutraljetTracks;

      //Radius of the jet
      double R = 0.4;

      //Loop over all D0 candidates in the event.
      //If there are more than one, the jet reconstruction is done for each D0 candidate separately
      //ignoring the other not reconstructed D0 candidates in the event.
      for (unsigned int nD0 = 0; nD0 < D0_fourmomentum.size(); nD0++) {

 	//for (int iTow = 0; iTow < 4800; iTow++) SumE[iTow] = 0;

	//Delete all energies calculated for hadr. corr. from previous event
	SumE.fill(0);
	
	//for (int iTow = 0; iTow < 4800; iTow++) cout << SumE[iTow] << endl;
//-----------D0-track--------------------------------------------------------

	//double d0pionmass = sqrt(D0_fourmomentum[nD0].px*D0_fourmomentum[nD0].px+D0_fourmomentum[nD0].py*D0_fourmomentum[nD0].py+D0_fourmomentum[nD0].pz*D0_fourmomentum[nD0].pz+M_PION_PLUS*M_PION_PLUS);
	//comp. testing

        //Defining the four-momentum of the D0 candidate
        fastjet::PseudoJet pj(D0_fourmomentum[nD0].px,D0_fourmomentum[nD0].py,D0_fourmomentum[nD0].pz,D0_fourmomentum[nD0].E); //comp. testing
        //fastjet::PseudoJet pj(D0_fourmomentum[nD0].px,D0_fourmomentum[nD0].py,D0_fourmomentum[nD0].pz,d0pionmass);


        //Set flag to 2 if the D0 candidate is a D0 and to -2 if it is a anti-D0
        //It cannot be -1 or 1 because the default flag in FastJet is -1.
        pj.set_user_index(D0_fourmomentum[nD0].D0_antiD0*2);
	TLorentzVector v(D0_fourmomentum[nD0].px,D0_fourmomentum[nD0].py,D0_fourmomentum[nD0].pz,D0_fourmomentum[nD0].E);
  	
        //Add the D0 candidate to the inclusive particle vector
        //if(abs(v.PseudoRapidity())>1.0) continue;  //TEST
        input_particles.push_back(pj);
        
        
       // cout << D0_fourmomentum[nD0].px << "," << D0_fourmomentum[nD0].py <<","<< D0_fourmomentum[nD0].pz << endl;
        
//-----------Neutral-tracks--------------------------------------------------

        //Fill array SumE with momenta of tracks which are matched to BEMC towers
        GetCaloTrackMomentum(picoDst,TVector3(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z()));
//cout << "vertex: " << event->primaryVertex().x() << "," << event->primaryVertex().y() << "," << event->primaryVertex().z() << endl;
        //Loop over all tracks in the event
        for (int iTow = 0; iTow < 4800; iTow++){

           //Get the tower hit
           StPicoBTowHit *towHit = picoDst->btowHit(iTow);

           //Check if the tower is bad or missing information
           if (!towHit || towHit->isBad()) continue;

           //Get the alternative counting method for the tower ID
           //NumericIndex - Counting from 0 to 4799
           //SoftId - Counting from 1 to 4800
           int realtowID = towHit->numericIndex2SoftId(iTow);

           //Initialize the tower energy
           double towE = 0;

           //Calculation of the tower energy depending on the year
           //In 2014, there was a problem with the energy calibration, so the energy has to be corrected
           
           if(mYear==2014) {
                  //Exclude bad towers, saved in JetInfo.h
                  if (BadTowerMap[realtowID-1]) continue;

                  //Calculate the tower energy
                  towE = GetTowerCalibEnergy(realtowID);           
                  ////towE = towHit->energy(); //Only test
              }
              if(mYear==2016) {
                  //Exclude bad towers, saved in Calibration2016.h
                  if (EnergyBadTowerMap2016[realtowID-1]<0) continue; //!!! Check -1
                  //Get the tower energy
                  towE = towHit->energy();
              }
              towE-= fHadronCorr*SumE[iTow];
              
           //If the tower energy is negative, set it to 0
           if (towE < 0) towE = 0;

           //Initialize the tower geometry
           //StEmcGeom* mEmcGeom;
           //mEmcGeom = StEmcGeom::getEmcGeom("bemc");
           StEmcPosition* mEmcPosition;
           mEmcPosition = new StEmcPosition();

           //Correct the eta of the tower for the vertex position
           //Because the loaded eta is w.r.t. the center of the TPC, but the vertex do not have to be in the center
           StThreeVectorF towerPosition = mEmcPosition->getPosFromVertex(StThreeVectorF(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z()), realtowID);
           float Toweta = towerPosition.pseudoRapidity();
           float Towphi = towerPosition.phi();

           //Calculate the transverse energy
           //max eta 1.05258 max difference: ET = 0.124452 for E = 0.2, if we cut on |Vz| < 30 cm
           double ET = towE/cosh(Toweta);

           //If the transverse energy is greater than 30 GeV, discard the event
           if (ET > 30) {
                TowArr.clear();
                TowEta.clear();
                TowPhi.clear();
                Clusters.clear();
                return kStOK;
           }
////////////////////////////////////////////////////////////////////////////////// comp. test
	   Double_t p = 1.0*TMath::Sqrt(towE*towE - 0*0);
	   double posX = towerPosition.x();
  	   double posY = towerPosition.y();
  	   double posZ = towerPosition.z();
  	     Double_t r = TMath::Sqrt(posX*posX + posY*posY + posZ*posZ) ;
  	                double px,py,pz;
  	   px = p*posX/r*1.;
    	   py = p*posY/r*1.;
    	   pz = p*posZ/r*1.;
    	  //// cout << "R: " << r << endl;
////////////////////////////////////////////////////////////////////////////////// comp. test
           //Initialize and calculate the momentum components
          /* comp. test
           double px,py,pz;
           px = ET*cos(Towphi);
           py = ET*sin(Towphi);
           pz = towE*tanh(Toweta); 
          */
           //Create a jet with the calculated momentum components
           PseudoJet inputTower(px, py, pz, towE);
		//if (realtowID==4706) cout << "px: " << px << " py: " << py << " pz: " << pz << " ET: " << ET << " sumE: " << SumE[iTow] << endl;
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
           } //End of minimum ET cut

        } //End of loop over all towers

        //NEventsCuts: Cal.: E_{T} < 30 GeV
        //Only for the first D0 otherwise it is counted multiple times
        if (nD0==0) NEventsCuts->Fill(9);
        

//-----------Charged-tracks--------------------------------------------------

        //Loop over all tracks in the event
        for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack){

            //The i-th track is loaded
            StPicoTrack* trk = picoDst->track(iTrack);

            //Check if the track exists
            if(!trk) continue;

            //Filling jet track histogram before any cuts
            JetTracksdEdx->Fill(trk->gPtot()*trk->charge(),trk->dEdx());

            //Check if the track is a good track
            if (!isGoodJetTrack(trk,event)) continue;

            //Loading of the pT
            double pT = trk->gMom().Perp();
            //Check if the pT is above 0.2 GeV/c or if it is NaN, because NaN!=NaN
            if(pT != pT) continue; // NaN test.
            //Loading of the eta, phi, dca and charge
            float eta = trk->gMom().PseudoRapidity();
            float phi = trk->gMom().Phi();
            float dca = (TVector3(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z()) - trk->origin()).Mag();
            float charged = trk->charge();

            //Filling jet track histogram before after cuts
            JetTracksdEdxCut->Fill(trk->gPtot()*trk->charge(),trk->dEdx());

            //If the track is not daughter pion nor kion then...
            if (DaughterPionTrackVector[nD0] != trk->id() && DaughterKaonTrackVector[nD0] != trk->id()){

                //Defining the four-momentum of the charged particle, assumed pi+ mass
                //                       px,       py,               pz,                                 E = sqrt(p^2 + m^2),
                fastjet::PseudoJet pj(trk->gMom().x(),trk->gMom().y(),trk->gMom().z(), sqrt(trk->gMom().Mag()*trk->gMom().Mag()+M_PION_PLUS*M_PION_PLUS));
       // cout << "input_particles.push_back({"<<trk->gMom().x()<<","<<trk->gMom().y()<<","<<trk->gMom().z()<<","<<sqrt(trk->gMom().Mag()*trk->gMom().Mag()+M_PION_PLUS*M_PION_PLUS)<<"}); //charged" << endl;
		//Set the flag to 3 if the particle is charged
		pj.set_user_index(3);
                //Add the charged particle to the charged particle vector
                chargedjetTracks.push_back(pj);
                //Add the charged particle to the inclusive particle vector
                input_particles.push_back(pj);

/*
		float bField = picoDst->event()->bField();
		StPicoPhysicalHelix helix = trk->helix(bField*kilogauss);
*/

            }//End of if the track is not daughter pion nor kion

            //Filling the track histograms
            hpT_tr->Fill(pT, reweight);
            heta_phi_tr->Fill(phi + TMath::Pi(), eta,  reweight);
            heta_tr->Fill(eta, reweight);
            hphi_tr->Fill(phi + TMath::Pi(), reweight); //to shift by pi
            hdca_z_tr->Fill(dca, event->primaryVertex().z(), reweight);
            hdca_pT->Fill(dca, pT, reweight);
            hdca_tr->Fill(dca, reweight);
            hcharged_tr->Fill(charged, reweight);
            /*
	    int ptBin = pT < 1  ? 0 :
		        pT < 4  ? 1 :
			pT < 10 ? 2 :
			pT < 15 ? 3 :
			pT < 20 ? 4 : 5;
*/
	    //TVector3 pVertex = TVector3(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z());
	    
            //double dcaSign = (1.*trk->gMom().x()/trk->gPt())*trk->gDCAy(pVertex.y()) - (1.*trk->gMom().y()/trk->gPt())*trk->gDCAx(pVertex.x());
            //int charge = charged < 0 ? 0 : 1;
            
           // hDcaSign[centrality][ptBin][charge]->Fill(dcaSign,reweight);

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
        
        unsigned int seed1 = 12345;
    	unsigned int seed2 = 56789;
	std::vector<int> seeds = { static_cast<int>(seed1), static_cast<int>(seed2) };
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
        float rho = bkgd_estimator.rho();
        //float emptyjets = bkgd_estimator.n_empty_jets(); //not used
        //float alljets = bkgd_estimator.n_jets_used(); //not used
        float rhom = bkgd_estimator.rho_m(); 
        
        //iterativ background estimation
        double max_eta=100;   
       	double max_eta_jet=3; // the maximal pseudorapidity for selected jets. Not important for the subtraction.


 	contrib::IterativeConstituentSubtractor subtractor;
  	subtractor.set_distance_type(contrib::ConstituentSubtractor::deltaR); 




        vector<double> max_distances;
	max_distances.push_back(0.15);
	//max_distances.push_back(0.2);
    
	vector<double> alphas;
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
        //unsigned int seed1 = 12345;
    	//unsigned int seed2 = 56789;
	//std::vector<int> seeds = { static_cast<int>(seed1), static_cast<int>(seed2) };
    	fastjet::AreaDefinition area_def = initial_area_def2.with_fixed_seed(seeds);
	//////
	
	
        //Definition of the clustering
        ClusterSequenceArea clust_seq_hard(input_particles, jet_def, area_def);
        
        //Jet minimum pT cut
        double ptmin = -0.1;
        //  fFilteredJets = fClustSeqSA->inclusive_jets(fMinJetPt-0.1); //becasue this is < not <=
        

        //Sorting of the jets by pT
        vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq_hard.inclusive_jets(ptmin-0.01));





//***************
        ////bkgd_estimator2.set_particles(input_particles); //Bug previously: chargedjetTracks instead of input_particles

	vector<PseudoJet> corrected_event=subtractor.subtract_event(input_particles); 

	// clustering of the corrected event
	ClusterSequence clust_seq_corr(corrected_event, jet_def);
	Selector sel_jets = SelectorAbsEtaMax(100);
	
	
	*/
	
	for (int difiter = 0; difiter < 4; difiter++){
	
	
	
	vector<PseudoJet> corrected_jets = JetReconstructionICS(input_particles, centrality, R, true, fPsi_2_shifted, difiter);
  		double rho = 0;

//***************


        //Loop over all jets
      //  for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
        //cout << "number: " << corrected_jets.size() << endl;
        for (unsigned int i = 0; i < corrected_jets.size(); i++) {
          
        

      	    //double ptA = inclusive_jets[i].pt();
    	    //double areaA = inclusive_jets[i].area();  // metoda 'area()' vyžaduje použití 'ClusterSequenceArea'
            //double pt_area_ratio = ptA / areaA;

            //Exclude jets with |eta| > 1 - R
            ////if (abs(inclusive_jets[i].pseudorapidity()) > (1.0 - R)) continue; Postponed to local analysis

	    //Constituent counter
	    int ConstCounter = 0;

            //Initialize the variables (angularities, z-value)
            int user_index = 0;
            double Delta_R_D0 = 0;
            double lambda_alpha_0 = 0.;
            double lambda_alpha_0half = 0.5;
            double lambda_alpha_1 = 1.;
            double lambda_alpha_1half = 1.5;
            double lambda_alpha_2 = 2.;
            double lambda_alpha_3 = 3.;
            double lambda_kappa_1 = 1.;
            double lambda_kappa_2 = 2.;
            double zet = 0;
            double lambda_1_0half = 0;
            double lambda_1_1 = 0;
            double lambda_1_1half = 0;
            double lambda_1_2 = 0;
            double lambda_1_3 = 0;
            double lambda_2_0 = 0;
            double neutralpT = 0;

            //Calculate the jet pT + background subtraction (= _corr)
            //pT(sub) = pT - rho * A_jet
            //Since z-value requires px and py, it is calculated as well
            double pT_jet = corrected_jets[i].perp();
            double pT_jet_corr = pT_jet;
            //double pT_jet_corr = pT_jet - rho * inclusive_jets[i].area();
            double px_jet = corrected_jets[i].px();
            double px_jet_corr = px_jet;
            //double px_jet_corr = px_jet - rho * inclusive_jets[i].area_4vector().px();
            double py_jet = corrected_jets[i].py();
            double py_jet_corr = py_jet;
            //double py_jet_corr = py_jet - rho * inclusive_jets[i].area_4vector().py();



            //Loading the constituents of the jet
         //   const vector<fastjet::PseudoJet>& constituents = inclusive_jets[i].constituents();

            const vector<fastjet::PseudoJet>& constituents = corrected_jets[i].constituents();
		
	  //  bool ExitFlag = false;		
		
            //Loop over all constituents of the i-th jet
            for (vector<fastjet::PseudoJet>::const_iterator particle = constituents.begin(); particle != constituents.end(); ++particle) {

                //Loading of the particle index
                int index = particle->user_index();
                
                //Number of constituents (+-2 = D0, 3 = charged, 10 = neutral, -1 ghost)
                if (index == -1) continue;
                //if (particle->pt() < 0.01); continue;
                
                //Constituent counter
                if(particle->pt() > 0.001) ConstCounter++;
                
                //Fraction of neutral particles
                if (index == 10) neutralpT += particle->pt();
                
                //Calculating the delta R = sqrt(delta eta^2 + delta phi^2)
                double Delta_R =delta_R(corrected_jets[i].eta(),corrected_jets[i].phi(),particle->eta(),particle->phi());
                double Delta_R2 =delta_R(corrected_jets[i].rap(),corrected_jets[i].phi(),particle->rap(),particle->phi());
                //Angularities are calculated only for track based particles (charged + D0)
                if (particle->user_index() != 10) {
                    lambda_1_0half+=	pow(particle->pt()/pT_jet_corr,lambda_kappa_1)*		pow( Delta_R /R ,lambda_alpha_0half);
                    lambda_1_1+=	pow(particle->pt()/pT_jet_corr,lambda_kappa_1)*		pow( Delta_R /R ,lambda_alpha_1);
                    lambda_1_1half+=	pow(particle->pt()/pT_jet_corr,lambda_kappa_1)*		pow( Delta_R /R ,lambda_alpha_1half);
                    lambda_1_2+=	pow(particle->pt()/pT_jet_corr,lambda_kappa_1)*		pow( Delta_R /R ,lambda_alpha_2);
                    lambda_1_3+=	pow(particle->pt()/pT_jet_corr,lambda_kappa_1)*		pow( Delta_R /R ,lambda_alpha_3);
 		    lambda_2_0+=	pow(particle->pt()/pT_jet_corr,lambda_kappa_2)*		pow( Delta_R /R ,lambda_alpha_0);
                }

                //Check if the constituent is D0 (D0 = 2, antiD0 = -2)
                if (abs(index) == 2 ) {
		    
                    //user_index is used as a flag to check if the jet contains D0
                    user_index = index;

                    //Delta R for D0
                    Delta_R_D0 = Delta_R;

                    //z-value, z = pT(D0)*^pT(jet)/|pT(jet)|
                    zet = (D0_fourmomentum[nD0].px*px_jet_corr+D0_fourmomentum[nD0].py*py_jet_corr)/(pT_jet_corr*pT_jet_corr);
                    
                    if (abs(index) == 2 && Delta_R > 1.0) {
    std::cout << std::fixed << std::setprecision(3)
              << "⚠️ Suspicious D0:\n"
              << "  ΔR        = " << Delta_R << "\n"
              << "  ΔR2        = " << Delta_R2 << "\n"
              << "  jet η     = " << corrected_jets[i].eta()
              << ", D0 η      = " << particle->eta() << "\n"
              << "  jet φ     = " << TVector2::Phi_0_2pi(corrected_jets[i].phi())
              << ", D0 φ      = " << TVector2::Phi_0_2pi(particle->phi()) << "\n"
              << "  user_index = " << index << "\n"
              << std::endl;
              std::cout << "  jet pt     = " << corrected_jets[i].pt()
          << ", D0 pt       = " << particle->pt() << "\n"
          <<", event ID     = " << picoDst->event()->eventId() << endl;

}

                    
                }
                    


            } //End of loop over all constituents of the i-th jet

            //Calculation of the fraction of neutral particles
            double nfraction = neutralpT/pT_jet;

            //If the fraction is too high, jet is rejected (Parameter is saved in RunPicoD0AnaMaker.C as setMaxNeutralFraction)
            ////if (nfraction > maxneutralfrac) continue; //Postponed to local analysis
            
            //if (area_jet < fAcuts) continue; //Not implemented

            //If the jet contains D0
            if (abs(user_index) ==2){
/*
            for (vector<fastjet::PseudoJet>::const_iterator particle = constituents.begin(); particle != constituents.end(); ++particle) {
	            if (particle->user_index() ==-1) continue;
       	                            cout << "px: " << particle->px()  << ", py: " << particle->py()  << ", pz: " << particle->pz() << endl;
            }
          */  
		//cout << "rho local: " << bge_rho.rho_m(inclusive_jets[i]) << endl;
                //Calculation of the D0 mass
                double D0mass = D0_fourmomentum[nD0].D0Mass;
                //Calculation of the D0 pT (pT=sqrt(px^2+py^2))
                double D0_pT = sqrt(D0_fourmomentum[nD0].px*D0_fourmomentum[nD0].px+D0_fourmomentum[nD0].py*D0_fourmomentum[nD0].py);

		//Fill the histogram (pt vs background density)
                Jet_grefmult_pt_background->Fill(picoDst->event()->grefMult(),rho);
                
                //Rapidity calculations and filling histogram
                TLorentzVector v(D0_fourmomentum[nD0].px,D0_fourmomentum[nD0].py,D0_fourmomentum[nD0].pz,D0_fourmomentum[nD0].E);
                //double D0_rapidity = 1./2.*log((D0_fourmomentum[nD0].E+D0_fourmomentum[nD0].pz)/(D0_fourmomentum[nD0].E-D0_fourmomentum[nD0].pz));
  		double D0_rapidity = v.Rapidity();
                double D0_pseudorapidity = v.PseudoRapidity();

                Jet_D0pT_vs_D0rapidity->Fill(D0_rapidity,D0_pT);
                Jet_D0pT_vs_Jetrapidity->Fill(corrected_jets[i].rapidity(),D0_pT);
                Jet_phi->Fill(corrected_jets[i].phi());
		
                //Fill the TNtuple
		TupleVariables[VariableJets["RunId"]] = RunId * user_index / 2.;             		// RunID, positive for D0, negative for antiD0
		TupleVariables[VariableJets["centrality"]] = centrality;                    		// Centrality
		TupleVariables[VariableJets["centr_weight"]] = reweight;                    		// Centrality reweighting factor
		TupleVariables[VariableJets["NJet"]] = D0_fourmomentum.size();              		// Number of D0 in the event
		TupleVariables[VariableJets["jet_eta"]] = corrected_jets[i].pseudorapidity(); 	    // Jet eta
		TupleVariables[VariableJets["jet_phi"]] = corrected_jets[i].phi();         		    // Jet phi
		TupleVariables[VariableJets["grefmult"]] = picoDst->event()->grefMult();    		// grefMult
		TupleVariables[VariableJets["bg_dens"]] = rho;                               		// density of the background
		TupleVariables[VariableJets["jet_area"]] = -999; /*corrected_jets[i].area();*/         		// area of the jet
		TupleVariables[VariableJets["jet_rap"]] = corrected_jets[i].rap();           		// rapidity of the jet
		TupleVariables[VariableJets["jet_pt"]] = pT_jet;                             		// Jet pT
		TupleVariables[VariableJets["jet_pt_corr"]] = pT_jet_corr;                   		// Jet pT after background subtraction
		TupleVariables[VariableJets["D0mass"]] = D0mass;                             		// D0 mass
		TupleVariables[VariableJets["D0_r"]] = Delta_R_D0;                           		// Delta R between D0 and jet axis
		TupleVariables[VariableJets["D0_pT"]] = D0_pT;                               		// D0 pT
		TupleVariables[VariableJets["lambda_1_0half"]] = lambda_1_0half;             		// Angularity lambda_1_0.5
		TupleVariables[VariableJets["lambda_1_1"]] = lambda_1_1;                     		// Angularity lambda_1_1
		TupleVariables[VariableJets["lambda_1_1half"]] = lambda_1_1half;             		// Angularity lambda_1_1.5
		TupleVariables[VariableJets["lambda_1_2"]] = lambda_1_2;                     		// Angularity lambda_1_2
		TupleVariables[VariableJets["lambda_1_3"]] = lambda_1_3;                     		// Angularity lambda_1_3
		TupleVariables[VariableJets["lambda_2_0"]] = lambda_2_0;                     		// Angularity lambda_2_0
		TupleVariables[VariableJets["z"]] = zet;                                     		// zet
		TupleVariables[VariableJets["NConst"]] = ConstCounter; 					// Number of constituents
		TupleVariables[VariableJets["NpTfraction"]] = nfraction;                     		// Neutral pT fraction
		TupleVariables[VariableJets["D0_rap"]] = D0_rapidity;                     		// D0 meson rapidity
        	TupleVariables[VariableJets["D0_eta"]] = D0_pseudorapidity;                    		// D0 meson pseudorapidity
        	TupleVariables[VariableJets["Psi2"]] = fPsi_2_shifted;                    		// D0 meson pseudorapidity
		//if (abs(inclusive_jets[i].pseudorapidity())>1.0) continue;
		
                //Fill the array to TNtuple                
                Jets[difiter]->Fill(TupleVariables);
                
                //if (abs(inclusive_jets[i].pseudorapidity()) > 1) continue;
              //  if (D0_pseudorapidity > 1) continue;
                     //                    cout << "uncorrected pt: "<< pT_jet <<" subtracted jet: " << subtracted_jet.perp() << " jet->Pt() " << pT_jet_corr << endl;
                     /*
		cout << "runID: " << RunId << endl;
		cout << "eventID: " << picoDst->event()->eventId() << endl;
		cout << "mass: " << D0mass << " ptD0: " << D0_pT << " pt: " << pT_jet << " pt_corr: " << pT_jet_corr << " rho: " << rho << /*" area: " << corrected_jets[i].area() <<*//* endl;
	        cout << "rap: " << corrected_jets[i].rap() << " eta: " << corrected_jets[i].pseudorapidity() << endl;
	        cout << "z: " <<zet << " D0eta: " << D0_pseudorapidity << " D0rap " << D0_rapidity <<endl;
	        */
                //Print colorfully the D0-jet information
                //printf("\033[32m%5u %15.8f %15.8f %15.8f %15d\033[0m\n", i, inclusive_jets[i].rap(), inclusive_jets[i].phi(), inclusive_jets[i].perp(), user_index);
//cout << D0mass << " " << D0_pT << " ptjet: " << pT_jet << " ptjetcorr: " << pT_jet_corr << " " << inclusive_jets[i].area() << " " << rho << endl;		
		//ExitFlag = true;
            } else{

                //Print the jet information
               //// printf("%5u %15.8f %15.8f %15.8f %15d\n", i, corrected_jets[i].rap(), corrected_jets[i].phi(), corrected_jets[i].perp(), user_index);

            } //End of if the jet contains D0
	  //if (ExitFlag) break;
        } //End of loop over all jets

        corrected_jets.clear();
	} //End of 4 different methods

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
Double_t StPicoD0AnaMaker::vertexCorrectedEta(double eta, double vz) {
    //Function to correct the eta value of a track for the z-position of the primary vertex

    //eta = -log(tan(theta/2)) => theta = 2*atan(exp(-eta))
    double tower_theta = 2.0 * atan(exp(-eta));

    //If eta = 0 then z = 0
    //Else calculate z position
    double z = 0.0;
    if (eta != 0.0) z = mBarrelRadius / tan(tower_theta);

    //Difference between the z position of the track and the z position of the primary vertex
    double z_diff = z - vz;

    //Calculate the corrected theta value
    double theta_corr = atan2(mBarrelRadius, z_diff);

    //Calculate the corrected eta value
    double eta_corr = -log(tan(theta_corr / 2.0));

    return eta_corr;

    //Function returns the corrected eta value
}
//---------------------------------------------------------------------------
Bool_t StPicoD0AnaMaker::GetCaloTrackMomentum(StPicoDst *mPicoDst, TVector3 mPrimVtx) {
  //Function to calculate the momentum of the track matched to the calorimeter tower

  //Loading of the number of tracks in the event
  UInt_t nTracks = mPicoDst->numberOfTracks();

  //Loop over all tracks in the event
  for (unsigned int itrack = 0; itrack < nTracks; itrack++) {

      //Loading the track
      StPicoTrack *trk = mPicoDst->track(itrack);
      //Loading of the global momentum
      TVector3 gMom = trk->gMom();

      //Loading of the pT
      double pT = gMom.Perp();
           
      //Check if the pT is above 0.2 GeV/c or if it is NaN, because NaN!=NaN
      if(pT != pT || pT < 0.2) continue;
      //Loading of eta
      float eta = gMom.PseudoRapidity();
      //Exclude tracks outside of the TPC acceptance
      if (fabs(eta) > 1) continue;
      //Loading of phi
      //float phi = gMom.Phi();

      //Loading of the number of hits
      float nHitsFit = trk->nHitsFit();
      //Loading of the number of hits possible
      float nHitsMax = trk->nHitsMax();
      //Exclude tracks with less than 15 hits or with a ratio of less than 0.52
      if (nHitsFit < 15 || 1.*nHitsFit/nHitsMax < 0.52) continue;

      //Loading of the value of the magnetic field
      double Bfield = mPicoDst->event()->bField();
      //Loading of the helix
      StPicoPhysicalHelix trkhelix = trk->helix(Bfield);

      //Loading of the primary vertex
      float vtx_x = mPrimVtx.x();
      float vtx_y = mPrimVtx.y();
      float vtx_z = mPrimVtx.z();

      //Calculation of the DCA to the primary vertex
      TVector3 dcaPoint = trkhelix.at(trkhelix.pathLength(vtx_x, vtx_y));
      //Calculation of the DCA in the x-y plane
      ////float dca_z = dcaPoint.z() - vtx_z; //Test
      
      float dca_z = trk->gDCAz(vtx_z); //? TO DO
      
      //Exclude tracks with a DCA to the primary vertex in z of more than maxdcazhadroncorr (in RunPicoD0AnaMaker.C)
      if (fabs(dca_z) > maxdcazhadroncorr) continue;

      //Initialization and loading of the tower index
      int TowIndex = -99999;
      TowIndex = trk->bemcTowerIndex(); //ID
      float p = 0;
      
      //Check if the track is matched to a tower
      if (TowIndex >= 0) {

        //Loading of the momentum
        p = gMom.Mag();          
	double TrackEnergy = 1.0*TMath::Sqrt(p*p + M_PION_PLUS*M_PION_PLUS);
        
        //Summing up the energy of all tracks matched to the same tower //Previously neglected pion mass 
        SumE[TowIndex] += TrackEnergy;
      }
  
  } //End of track loop

  return true;

  //Function returns true if it was successful and SumE filled with the momentum of all tracks matched to a tower
}
//---------------------------------------------------------------------------
Double_t StPicoD0AnaMaker::GetTowerCalibEnergy(Int_t TowerId){
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
bool StPicoD0AnaMaker::IsBadEnergyRun(int runID) {
    // Check if the run is in the list of BEMC bad runs

    for (unsigned int i = 0; i < sizeof(EnergyBadRunList2016)/sizeof(EnergyBadRunList2016[0]); i++) {
        if (EnergyBadRunList2016[i] == runID) {
            return true;
        }
    }
    return false;

    //Function returns true if the run is in the list of bad runs
}
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
int StPicoD0AnaMaker::isD0PairCentrality_pt(StKaonPion const* const kp, int Centrality, int mYear) const{
    //Check if the pair passes the cuts for D0

    //Loading the daughter particles tracks
    StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
    StPicoTrack const* pion = picoDst->track(kp->pionIdx());

    //Initialisation of the pairCuts boolean
    bool pairCuts = false;

    //Recalculation of the centrality binning
    int Centrality2 =  (Centrality == 8 || Centrality == 7) ? 0 :
                       (Centrality == 6) ? 1 :
                       (Centrality == 5 || Centrality == 4) ? 2 :
                       (Centrality == 3 || Centrality == 2) ? 3 :
                       (Centrality == 1 || Centrality == 0) ? 4 : -1;
    // Centr.   0-10%       10-20%         20-40%         40-60%        60-80%
    // C_ID     8,7        6               5,4             3,2         1,0
    // bin       0         1                 2               3           4


    //Recalculation of the momentum binning
    int KPMom = (kp->pt() < 0.5) ? 0 : (kp->pt() < 1) ? 1 : (kp->pt() < 2) ? 2 : (kp->pt() < 3) ? 3 : (kp->pt() < 5) ? 4 : 5;

    //Check if the pair passes the particular cuts
    //Parameters are saved in StCuts.cxx
    if (mYear == 2014){
        pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < mycuts::DCA_D0_cut_2014[KPMom][Centrality2] &&
                    kp->pionDca() > mycuts::pionDCA_cut_2014[KPMom][Centrality2] && kp->kaonDca() > mycuts::kaonDCA_cut_2014[KPMom][Centrality2] &&
                    kp->dcaDaughters() < mycuts::pionkaonDCA_cut_2014[KPMom][Centrality2] && kp->decayLength()> mycuts::D0_decayLength_cut_2014[KPMom][Centrality2] &&
                    cos(kp->pointingAngle()) > mycuts::cosTheta_2014;
    } else if (mYear == 2016){
        pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < mycuts::DCA_D0_cut_2016[KPMom][Centrality2] &&
                    kp->pionDca() > mycuts::pionDCA_cut_2016[KPMom][Centrality2] && kp->kaonDca() > mycuts::kaonDCA_cut_2016[KPMom][Centrality2] &&
                    kp->dcaDaughters() < mycuts::pionkaonDCA_cut_2016[KPMom][Centrality2] && kp->decayLength()> mycuts::D0_decayLength_cut_2016[KPMom][Centrality2] &&
                    cos(kp->pointingAngle()) > mycuts::cosTheta_2016;
    }

    //Calculation of the product of the daughter charges
    int charge = kaon->charge() * pion->charge();

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
bool StPicoD0AnaMaker::isGoodEvent(int mYear,  TH1F* NEventsCuts){
  //Check if the event passes the cuts (MB triggers, vr, vz, vzVpdVz)

  //Loading event
  StPicoEvent *event = (StPicoEvent *)picoDst->event();



  //Checking triggers
  if (!isMBTrigger(mYear)) return false;

  double ShiftVz = 0;
  if (mYear == 2014) ShiftVz = 0; //2.1486; //It could be used

  //NEventsCuts: Triggers
  NEventsCuts->Fill(1);
  //Checking vr = sqrt(vx^2+vy^2)
  if (!(sqrt(event->primaryVertex().x()*event->primaryVertex().x()+event->primaryVertex().y()*event->primaryVertex().y()) < mycuts::vr)) return false;
  //NEventsCuts: v_r
  NEventsCuts->Fill(2);
  //Checking vz
  if (!(fabs(event->primaryVertex().z()-ShiftVz) < mycuts::vz)) return false;
  //NEventsCuts: v_z
  NEventsCuts->Fill(3);
  //Checking vzVpdVz
  if (!(fabs(event->primaryVertex().z() - event->vzVpd()) < mycuts::vzVpdVz)) return false;
  //NEventsCuts: |v_z - v_z_vpd|
  NEventsCuts->Fill(4);
  
  //Check on suspicious all-0 position
  bool nonezeroVertex = (event->primaryVertex().x()!=0 && event->primaryVertex().y()!=0 && event->primaryVertex().z()!=0);
  if (!nonezeroVertex) return false;
  NEventsCuts->Fill(5);

  return true;

  //Function returns true if the event passes the cuts
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isMBTrigger(int mYear){
 //Function checks if the event is triggered by the MB trigger

 //Initialization of the set of triggers
 const std::set<int>* mbTriggers = nullptr;

 //Different triggers for different years, saved in StCuts.cxx
 if(mYear ==2016) mbTriggers = &mycuts::mbTriggers2016;
 if(mYear ==2014) mbTriggers = &mycuts::mbTriggers2014;

 //Loading event and checking if it is triggered by the MB trigger
 StPicoEvent* event = static_cast<StPicoEvent*>(mPicoDstMaker->picoDst()->event());
 return std::any_of(mbTriggers->begin(), mbTriggers->end(), [&](int trigger) { return event->isTrigger(trigger); });

 //Function returns true if the event is triggered by the MB trigger
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodTrack(StPicoTrack const * const trk) const{
    // Require at least one hit on every layer of PXL and IST.
    // It is done here for tests on the preview II data.
    // The new StPicoTrack which is used in official production has a method to check this

    //Check if the track meets the HFT requirement
    //2014 - Require at least one hit on every layer of PXL and IST
    //2016 - Require at least one hit on every layer of PXL and (IST or SST)
    //Both can be written in te same way
    bool HFTCondition = (trk->hasPxl1Hit() && trk->hasPxl2Hit()) && (trk->hasSstHit() || trk->hasIstHit());

    //Check if |eta| < 1
    bool EtaCondition = abs(trk->gMom().PseudoRapidity()) < 1;

    //In StCuts.cxx is defined if the HFT is required and the nHitsFit and minPt values.
    return (trk->gPt() > mycuts::minPt && trk->nHitsFit() >= mycuts::nHitsFit && (HFTCondition || !mycuts::requireHFT) && EtaCondition);



    //Function returns true if track is good
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodJetTrack(StPicoTrack const * const trk,StPicoEvent const *const myEvent) const{
    //Check if the track meets the jet track cuts
    //Parameters saved in StCuts.cxx

    //pT range cut
    bool pTTrackJetCut = trk->gPt() > mycuts::jetTrackPtMin && trk->gPt() < mycuts::jetTrackPtMax;
    //eta range cut
    bool etaTrackJetCut = fabs(trk->gMom().PseudoRapidity()) < mycuts::jetTrackEta;
    //nHitsFit cut
    bool nHitsTrackJetCut = trk->nHitsFit() >= mycuts::jetTracknHitsFit;
    //nHitsRatio cut
    bool nHitsRatioTrackJetCut = (1.0*trk->nHitsFit()/trk->nHitsMax())>=mycuts::jetTracknHitsRatio;
    //DCA cut
    bool dcaTrackJetCut = fabs(trk->gDCA(myEvent->primaryVertex().x(),myEvent->primaryVertex().y(),myEvent->primaryVertex().z())) < mycuts::jetTrackDCA;

    return pTTrackJetCut && etaTrackJetCut && nHitsTrackJetCut && nHitsRatioTrackJetCut && dcaTrackJetCut;

    //Return true if all the cuts are passed
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodJetTrackSim(TVector3 simTrack, Int_t reco, Double_t simDca) const{
    //Check if the track meets the jet track cuts
    //Parameters saved in StCuts.cxx

    //pT range cut
    bool pTTrackJetCut = simTrack.Perp() > mycuts::jetTrackPtMin && simTrack.Perp() < mycuts::jetTrackPtMax;
    //eta range cut
    bool etaTrackJetCut = fabs(simTrack.PseudoRapidity()) < mycuts::jetTrackEta;
    //nHitsFit cut
    bool nHitsTrackJetCut = abs(Track_mNHitsFit[reco]) >= mycuts::jetTracknHitsFit;
    //nHitsRatio cut
    bool nHitsRatioTrackJetCut = (1.0 * abs(Double_t(Track_mNHitsFit[reco]) / Double_t(Track_mNHitsMax[reco]))) >= mycuts::jetTracknHitsRatio;
    //DCA cut
    bool dcaTrackJetCut = simDca < mycuts::jetTrackDCA;

    return pTTrackJetCut && etaTrackJetCut && nHitsTrackJetCut && nHitsRatioTrackJetCut && dcaTrackJetCut;

    //Return true if all the cuts are passed
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTpcPion(StPicoTrack const * const trk) const{
    //Check if the track meets nsigma (TPC) requirement for pion
    //Parameters are saved in StCuts.cxx

    return fabs(trk->nSigmaPion()) < mycuts::nSigmaPion;

    //Function returns true if track is good
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTpcKaon(StPicoTrack const * const trk) const{
    //Check if the track meets nsigma (TPC) requirement for kaon
    //Parameters are saved in StCuts.cxx

    return fabs(trk->nSigmaKaon()) < mycuts::nSigmaKaon; //In D0 event maker it is set to 2

    //Function returns true if track is good
}
//-----------------------------------------------------------------------------
float StPicoD0AnaMaker::getTofBeta(StPicoTrack const * const trk, StThreeVectorF const* const pVtx,StPicoDst const* const picoDst) const{
    //Calculation of beta for the track

    //index2tof is index of the track in the StPicoBTofPidTraits array
    int index2tof = trk->bTofPidTraitsIndex();

    //Initialization of beta
    float beta = std::numeric_limits<float>::quiet_NaN();

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
                float L = tofPathLength(pVtx, &btofHitPos, helix.curvature());

                //Calculation of the time of flight
                float tof = tofPid->btof();

                //Calculation of beta for positive values of tof
                if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
                //...else beta is not defined
                else beta = std::numeric_limits<float>::quiet_NaN();

            } //End of beta < 1e-4

        } //End of tofPid != NULL

    } //End of index2tof >= 0

    return beta;

    //Function returns beta of the track
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTofKaon(StPicoTrack const * const trk, float beta) const{
    //Check if the track meets |1/beta-1/beta_K| (TOF) requirement for kaon

    //Initialization of tofKaon
    bool tofKaon = false;

    //If beta is positive, than we can calculate |1/beta-1/beta_K|
    if(beta>0){

        //Calculation of the global total momentum
        double ptot = trk->gMom().Mag();

        //Calculation of the expected beta for kaons
        float beta_k = ptot/sqrt(ptot*ptot+M_KAON_PLUS*M_KAON_PLUS);

        //Check if the track meets the TOF requirement
        //Parameters are saved in StCuts.cxx
        tofKaon = fabs(1/beta - 1/beta_k) < mycuts::kTofBetaDiff ? true : false;

    } //End of beta > 0

    return tofKaon;

    //Function returns true if track is good based on TOF information
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTofPion(StPicoTrack const * const trk, float beta) const{
    //Check if the track meets |1/beta-1/beta_pi| (TOF) requirement for pion

    //Initialization of tofPion
    bool tofPion = false;
    
    //If beta is positive, than we can calculate |1/beta-1/beta_K|
    if(beta>0){

        //Calculation of the global total momentum
        double ptot = trk->gMom().Mag();
        //Calculation of the expected beta for kaons
        float beta_pi = ptot/sqrt(ptot*ptot+M_PION_PLUS*M_PION_PLUS);
        //Check if the track meets the TOF requirement
        //Parameters are saved in StCuts.cxx
        tofPion = fabs(1/beta - 1/beta_pi) < mycuts::pTofBetaDiff ? true : false;
    } //End of beta > 0

    return tofPion;

    //Function returns true if track is good based on TOF information
}





