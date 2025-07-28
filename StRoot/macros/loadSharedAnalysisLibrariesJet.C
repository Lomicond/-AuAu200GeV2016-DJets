#include "TSystem.h"

#ifndef __CINT__
#include "StMuDSTMaker/COMMON/macros/loadSharedLibraries.C"
#endif

extern TSystem* gSystem;

void loadSharedAnalysisLibraries() 
{

  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  /*
  gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/fastjet-install/lib/libfastjet");
  gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/fastjet-install/lib/libsiscone");
  gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/fastjet-install/lib/libsiscone_spherical"); 
  gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/fastjet-install/lib/libfastjetplugins");
  gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/fastjet-install/lib/libfastjettools");
  */
  
    gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/Alex_install/install/fastjet-install/lib/libfastjet");
    gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/Alex_install/install/fastjet-install/lib/libsiscone");
    gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/Alex_install/install/fastjet-install/lib/libsiscone_spherical");
    gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/Alex_install/install/fastjet-install/lib/libfastjetplugins");
    gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/Alex_install/install/fastjet-install/lib/libfastjettools");
   gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/Alex_install/install/fastjet-install/lib/libfastjetcontribfragile");
  //gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/fastjet-install/lib/libfastjetcontribfragile"); 
  
  //
   // gSystem->Load("/cvmfs/star.sdcc.bnl.gov/star-spack/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/fastjet-3.3.4-j5evuymea6juu4tkqxxim6nj3z6ldbg3/lib/libfastjet");
  //  gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/Alex_install/install/fastjet-install/lib/libsiscone");
   // gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/Alex_install/install/fastjet-install/lib/libsiscone_spherical");
    //gSystem->Load("/cvmfs/star.sdcc.bnl.gov/star-spack/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/fastjet-3.3.4-j5evuymea6juu4tkqxxim6nj3z6ldbg3/lib/libfastjetplugins");
   // gSystem->Load("/cvmfs/star.sdcc.bnl.gov/star-spack/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/fastjet-3.3.4-j5evuymea6juu4tkqxxim6nj3z6ldbg3/lib/libfastjettools");
  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StRefMultCorr");
  gSystem->Load("StPicoPrescales");
  gSystem->Load("StStrangeMuDstMaker");
  gSystem->Load("St_db_Maker");
  gSystem->Load("StDaqLib");

  gSystem->Load("StEmcRawMaker");
  gSystem->Load("StEmcADCtoEMaker");
  gSystem->Load("StEpcMaker");
  gSystem->Load("StTriggerUtilities");

  
 gSystem->Load("StPicoCutsBase");
  gSystem->Load("StBTofUtil");
  gSystem->Load("StPicoD0EventMaker");
   // gSystem->Load("StPicoBackground");
  gSystem->Load("StPicoHFMaker");
  gSystem->Load("StPicoCuts");
  gSystem->Load("StTpcDb");
  gSystem->Load("StDbUtilities");
  gSystem->Load("Sti");
  gSystem->Load("StiUtilities");
  gSystem->Load("StSsdDbMaker");
  gSystem->Load("StSvtDbMaker");
  gSystem->Load("StiMaker");
  gSystem->Load("StDbBroker");
  gSystem->Load("libgeometry_Tables"); //rember, order of loading makers matters
  gSystem->Load("StPicoD0JetAnaMaker");
  // gSystem->Load("StPicoBackground");

	

  cout << " loading of shared  libraries are done" << endl;
}
