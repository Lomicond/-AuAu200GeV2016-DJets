#ifndef PICOPRESCALESCONSTANTS_H
#define PICOPRESCALESCONSTANTS_H

#include <vector>

namespace PicoPrescalesConstants
{
  std::vector<unsigned int> const triggerId = {
	// st_physics stream
  
	450050,    // vpdmb-5-p-nobsmd-hlt (production_mid_2014, production_low_2014)
	450060,    // vpdmb-5-p-nobsmd-hlt (production_mid_2014, production_low_2014)
	450005,    // vpdmb-5-p-nobsmd (production_2014)
	450015,    // vpdmb-5-p-nobsmd (production_2014, production_mid_2014, production_low_2014)
	450025,    // vpdmb-5-p-nobsmd (production_mid_2014, production_low_2014)
	450014,    // VPDMB-5-nobsmd
	450024,    // VPDMB-5-nobsmd
	450008,    // VPDMB-5 (production_2014, production_mid_2014, production_low_2014)
	450018,    // VPDMB-5 (production_2014, production_mid_2014, production_low_2014)
	450010,    // VPDMB-30 (production_2014, production_mid_2014, production_low_2014)
	450020,    // VPDMB-30 (production_2014, production_mid_2014, production_low_2014)
	450013,    // VPD-ZDC-novtx-mon (production_2014, production_mid_2014, production_low_2014)
	450023,    // VPD-ZDC-novtx-mon (production_2014, production_mid_2014, production_low_2014)
	450009,    // vpdmb-5-p-nobsmd-ssd-hlt (production_mid_2014, production_low_2014)
	450012,    // ZDC-mon (production_2014, production_mid_2014, production_low_2014)
	450022,    // ZDC-mon (production_2014, production_mid_2014, production_low_2014)
	450011,    // MB-mon (production_2014, production_mid_2014, production_low_2014)
	450021,    // MB-mon (production_2014, production_mid_2014, production_low_2014)
	450103,    // Central-5 (production_2014, production_mid_2014, production_low_2014)
	450201,    // BHT1*VPDMB-30 (production_2014, production_mid_2014, production_low_2014)
	450211,    // BHT1*VPDMB-30 (production_2014, production_mid_2014, production_low_2014)
	450202,    // BHT2*VPDMB-30 (production_2014, production_mid_2014, production_low_2014)
	450212,    // BHT2*VPDMB-30 (production_2014, production_mid_2014, production_low_2014)
	450203,    // BHT3 (production_2014, production_mid_2014, production_low_2014)
	450213,    // BHT3 (production_2014, production_mid_2014, production_low_2014)
	
  520001,    // VPDMB-5-p-sst (production 1, physics stream)
  520011,    // VPDMB-5-p-sst 
  520021,    // VPDMB-5-p-sst
  520031,    // VPDMB-5-p-sst
  520041,    // VPDMB-5-p-sst
  520051,    // VPDMB-5-p-sst

  570002,    // VPDMB-5-nosst (production 2, nosst stream)
  570001, // VPDMB-5-sst (production 2, sst stream )

  520201, 520211, 520221, 520231, 520241, 520251, 520261,   // BHT1*VPDMB-10
  520203,                                                   // BHT3
  520101, 520111, 520121, 520131, 520141,                   // central-5  
  520007, 520017, 520027, 520037,                           // vpdmb-10
  520003, 520013, 520023, 520033, 520043,                   // VPDMB-5
  520802, 520812, 520822, 520832, 520842,                   // VPDMB-5-p-hlt
  520002, 520012, 520022, 520032, 520042                    // VPDMB-5-p-nosst
  };

  std::vector<unsigned int> const triggerIdMtd = {
  // st_mtd stream
  450601, // dimuon
  450611, // dimuon
  450621, // dimuon
  450631, // dimuon
  450641, // dimuon
  450604,    // dimuon-30-hft (production_2014)
  450605,    // dimuon-5-hft (production_mid_2014, production_low_2014)
  450606,    // dimuon-5-hft (production_mid_2014)
  450602,    // e-mu (production_2014)
  450612,    // e-mu (production_2014, production_low_2014)
  450622,    // e-mu (production_2014, production_low_2014)
  450632,    // e-mu (production_mid_2014)
  450642,    // e-mu (production_2014, production_low_2014)
  450600,    // single-muon (production_2014)
  450610,    // single-muon (production_2014, production_low_2014)
  450620,    // single-muon (production_2014, production_low_2014)
  450630,    // single-muon (production_mid_2014)
  450640     // single-muon (production_2014, production_low_2014)
  }; 
}
#endif
