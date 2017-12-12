///
/// \file TestEMCALReconstruction.C
/// \ingroup EMCAL_TestSimRec
/// \brief Simple macro to test EMCAL Reconstruction
///
/// Simple macro to test EMCAL reconstruction
///
/// \author Jenn Klay, LLNL
///

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TObjectTable.h>
#include <TStopwatch.h>
#include <TSystem.h>

#include "AliReconstruction.h"
#include "AliLog.h"

#endif

///
/// Main execution method
///
/// \param nev: number of events
///
void TestEMCALReconstruction(Int_t nev =-1) 
{
  // Uncomment to debug with different verbosity levels
  //  AliLog::SetModuleDebugLevel("EMCAL",100);
  //  AliLog::SetGlobalDebugLevel(1000);
  
  AliReconstruction rec;
  
  rec.SetRunTracking("");
  rec.SetRunVertexFinder(kFALSE);
  
  //calls local reconstruction of EMCAL and filling of ESD
  rec.SetRunLocalReconstruction("EMCAL");  //only do emcal
  rec.SetFillESD("EMCAL");
  rec.SetEventRange(0,nev);
  rec.SetRunQA(":");
  
  // Decomment this line in case of real data,
  // add the proper name of the file
  //rec.SetInput("raw.root");
  
  //OCDB settings
  rec.SetDefaultStorage("local://$ALIROOT_OCDB_ROOT/OCDB");
  
  // Decomment this line in case of anchored MC runs or data,
  // with the appropriate year
  //rec.SetDefaultStorage("alien://Folder=/alice/data/2011/OCDB");
  
  rec.SetSpecificStorage("GRP/GRP/Data",
                         Form("local://%s",gSystem->pwd()));
  
  TStopwatch timer;
  timer.Start();
  
  rec.Run();
  
  timer.Stop();
  timer.Print();
  
  gObjectTable->Print();
  
}
