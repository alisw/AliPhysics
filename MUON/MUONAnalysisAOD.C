#if !defined(__CINT__) || defined(__MAKECINT__)
// ROOT includes
#include <TSystem.h>
#include <Riostream.h>

#include <AliESD.h>

#include "AliAnalysis.h"
#include "AliAODRun.h"
#include "AliAOD.h"
#include "AliAODParticle.h"
#include "AliReaderESDTree.h"
#include "AliRunAnalysis.h"
#include "AliMuonAnalysis.h"
#endif


/////////////////////////////////////////////////////////////////////////////
// first version for analysing the invariant mass plot in the AOD framework
// since the framework is in progress, this macro may change as a consequence
// don't forget the include with ".includepath $ALICE_ROOT/ANALYSIS"
// if you compile this macro comment out 
// gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libANALYSIS.so");
// and launch this command in Aliroot outside the macro
//
// Ch. Finck
/////////////////////////////////////////////////////////////////////////////

void analysisAOD(Int_t first = 0, Int_t last = 99999, const Char_t* esdfilename = "AliESDs.root")
{

  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libANALYSIS.so");

  AliVAODParticle::SetDebug(0);
  AliRunAnalysis* analysis = new AliRunAnalysis();
  
  // set reader for ESD
  AliReaderESDTree* reader = new AliReaderESDTree(esdfilename);
  // active muon ESD reader
  reader->SetReadMuon(kTRUE);
  // disable central barrel (default = kTRUE)
  reader->SetReadCentralBarrel(kFALSE);
  // disable simulated data (not implemented yet)
  reader->ReadSimulatedData(kFALSE);
  // number of event to be read
  reader->ReadEventsFromTo(first,last);

  // set MUON analysis
  AliMuonAnalysis* muon = new AliMuonAnalysis();
  analysis->SetReader(reader);
  
  analysis->Add(muon);
  analysis->Run();
  delete analysis;
}
