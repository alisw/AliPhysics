/* $Id$ */

#include "AlidNdEtaAnalysisSelector.h"

#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TH2F.h>
#include <TChain.h>
#include <TFile.h>

#include <AliLog.h>
#include <AliGenEventHeader.h>
#include <AliHeader.h>

#include "dNdEtaAnalysis.h"
#include "dNdEtaCorrection.h"

ClassImp(AlidNdEtaAnalysisSelector)

AlidNdEtaAnalysisSelector::AlidNdEtaAnalysisSelector() :
  AliSelector(),
  fdNdEtaAnalysis(0),
  fdNdEtaCorrection(0)
{
  //
  // Constructor. Initialization of pointers
  //
}

AlidNdEtaAnalysisSelector::~AlidNdEtaAnalysisSelector()
{
  //
  // Destructor
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
}

void AlidNdEtaAnalysisSelector::SlaveBegin(TTree * tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  AliSelector::SlaveBegin(tree);

  fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta", "dndeta");
}

void AlidNdEtaAnalysisSelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  AliSelector::SlaveTerminate();

  // Add the histograms to the output on each slave server
  if (!fOutput)
  {
    AliDebug(AliLog::kError, Form("ERROR: Output list not initialized."));
    return;
  }

  fOutput->Add(fdNdEtaAnalysis);
}

void AlidNdEtaAnalysisSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliSelector::Terminate();

  fdNdEtaAnalysis = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndeta"));

  if (!fdNdEtaAnalysis)
  {
    AliDebug(AliLog::kError, Form("ERROR: Histograms not available %p", (void*) fdNdEtaAnalysis));
    return;
  }

  fdNdEtaAnalysis->Finish(fdNdEtaCorrection);

  TFile* fout = new TFile("out.root","RECREATE");
  WriteObjects();
  fout->Write();
  fout->Close();

  fdNdEtaAnalysis->DrawHistograms();
}

void AlidNdEtaAnalysisSelector::WriteObjects()
{
  // Write objects to output file
  // this is an extra function to be overloaded...
  //

  fdNdEtaAnalysis->SaveHistograms();
}
