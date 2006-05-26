/* $Id$ */

#include "AlidNdEtaAnalysisSelector.h"

#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TH2F.h>

#include <AliLog.h>
#include <AliGenEventHeader.h>

#include "dNdEtaAnalysis.h"

ClassImp(AlidNdEtaAnalysisSelector)

AlidNdEtaAnalysisSelector::AlidNdEtaAnalysisSelector() :
  AliSelector(),
  fdNdEtaAnalysis(0),
  fdNdEtaAnalysisFinal(0)
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

  fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta");
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

  fOutput->Add(fdNdEtaAnalysis->GetEtaVsVtxHistogram());
  fOutput->Add(fdNdEtaAnalysis->GetEtaVsVtxUncorrectedHistogram());
  fOutput->Add(fdNdEtaAnalysis->GetVtxHistogram());

  fdNdEtaAnalysis->GetVtxHistogram()->Print();
  fOutput->Print();
}

void AlidNdEtaAnalysisSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliSelector::Terminate();

  TH2F* etaVsVtxHistogram = dynamic_cast<TH2F*> (fOutput->FindObject("eta_vs_vtx"));
  TH2F* etaVsVtxUncorrectedHistogram = dynamic_cast<TH2F*> (fOutput->FindObject("eta_vs_vtx_uncorrected"));
  TH1D* vtxHistogram = dynamic_cast<TH1D*> (fOutput->FindObject("vtx"));

  if (!etaVsVtxHistogram || !vtxHistogram || !etaVsVtxUncorrectedHistogram)
  {
     AliDebug(AliLog::kError, Form("ERROR: Histograms not available %p %p %p", (void*) etaVsVtxHistogram, (void*) etaVsVtxUncorrectedHistogram, (void*) vtxHistogram));
    return;
  }

  fdNdEtaAnalysisFinal = new dNdEtaAnalysis("dNdEtaResult");

  fdNdEtaAnalysisFinal->SetEtaVsVtxHistogram(etaVsVtxHistogram);
  fdNdEtaAnalysisFinal->SetEtaVsVtxUncorrectedHistogram(etaVsVtxUncorrectedHistogram);
  fdNdEtaAnalysisFinal->SetVtxHistogram(vtxHistogram);

  fdNdEtaAnalysisFinal->Finish();

  TFile* fout = new TFile("out.root","RECREATE");
  WriteObjects();
  fout->Write();
  fout->Close();

  fdNdEtaAnalysisFinal->DrawHistograms();
}

void AlidNdEtaAnalysisSelector::WriteObjects()
{
  // Write objects to output file
  // this is an extra function to be overloaded...
  //

  fdNdEtaAnalysisFinal->SaveHistograms();
}
