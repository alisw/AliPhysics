/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

// The ESD is available as member fESD
//
// The Process function is nearly empty. Implement your analysis there and look at the other listed below functions you
// might need.
//
// The following methods can be overrriden. Please do not forgot to call the base class function.
//
//    Begin():        called everytime a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Init():         called for each new tree. Enable/Disable branches here.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
//  Author: Jan.Fiete.Grosse-Oetringhaus@cern.ch

#include "AliVertexSelector.h"

#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TAxis.h>

#include <AliLog.h>
#include <AliESD.h>
#include <AliHeader.h>
#include <AliGenEventHeader.h>

#include <AliPWG0Helper.h>
#include "esdTrackCuts/AliESDtrackCuts.h"

ClassImp(AliVertexSelector)

AliVertexSelector::AliVertexSelector() :
  AliSelectorRL(),
  fEsdTrackCuts(0),
  fVertexMC(0),
  fVertexESD(0),
  fVertexCorr(0),
  fVertexCorr2(0),
  fFakes(0)
{
  //
  // Constructor. Initialization of pointers
  //
}

AliVertexSelector::~AliVertexSelector()
{
  //
  // Destructor
  //
}

void AliVertexSelector::SlaveBegin(TTree* tree)
{
  AliSelectorRL::SlaveBegin(tree);

  fVertexMC = new TH1F("fVertexMC", "fVertexMC;MC vtx-z;Count", 501, -10, 10);
  fVertexESD = new TH1F("fVertexESD", "fVertexESD;ESD vtx-z;Count", 501, -10, 10);
  fVertexCorr = new TH2F("fVertexCorr", "fVertexCorr;MC vtx-z;ESD vtx-z - MC vtx-z", 40, -10, 10, 200, -1, 1);
  fVertexCorr2 = new TH2F("fVertexCorr2", "fVertexCorr2;MC vtx-z;ESD vtx-z - MC vtx-z", 40, -10, 10, 200, -1, 1);
  fFakes = new TH2F("fFakes", "fFakes;type;label", 5, 0, 5, 1001, -500.5, 500.5);

  if (!fEsdTrackCuts && fInput)
    fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (fInput->FindObject("AliESDtrackCuts"));
}

void AliVertexSelector::Init(TTree* tree)
{
  // read the user objects

  AliSelectorRL::Init(tree);

  // Enable only the needed branches
  if (tree)
  {
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("fTriggerMask", 1);
    tree->SetBranchStatus("fSPDVertex*", 1);
    tree->SetBranchStatus("fTracks.fLabel", 1);

    AliESDtrackCuts::EnableNeededBranches(tree);
  }
}

Bool_t AliVertexSelector::Process(Long64_t entry)
{
  //
  // Implement your analysis here. Do not forget to call the parent class Process by
  // if (AliSelector::Process(entry) == kFALSE)
  //   return kFALSE;
  //

  if (AliSelectorRL::Process(entry) == kFALSE)
    return kFALSE;

  // Check prerequisites
  if (!fESD)
  {
    AliDebug(AliLog::kError, "ESD branch not available");
    return kFALSE;
  }

  if (!fEsdTrackCuts)
  {
    AliDebug(AliLog::kError, "fESDTrackCuts not available");
    return kFALSE;
  }

  if (AliPWG0Helper::IsEventTriggered(fESD) == kFALSE)
  {
    AliDebug(AliLog::kDebug+1, Form("Skipping event %d because it was not triggered", (Int_t) entry));
    return kTRUE;
  }

  if (AliPWG0Helper::IsVertexReconstructed(fESD) == kFALSE)
  {
    AliDebug(AliLog::kDebug+1, Form("Skipping event %d because its vertex was not reconstructed", (Int_t) entry));
    return kTRUE;
  }

  AliHeader* header = GetHeader();
  if (!header)
  {
    AliDebug(AliLog::kError, "Header not available");
    return kFALSE;
  }

  // get the MC vertex
  AliGenEventHeader* genHeader = header->GenEventHeader();

  TArrayF vtxMC(3);
  genHeader->PrimaryVertex(vtxMC);

  // ########################################################
  // get the ESD vertex
  const AliESDVertex* vtxESD = fESD->GetVertex();
  Double_t vtx[3];
  vtxESD->GetXYZ(vtx);

  fVertexMC->Fill(vtxMC[2]);
  fVertexESD->Fill(vtx[2]);
  fVertexCorr->Fill(vtxMC[2], vtx[2] - vtxMC[2]);

  Float_t correctedVertex = vtx[2] + 2.951034e-03 + 6.859620e-04 * vtxMC[2];

  fVertexCorr2->Fill(vtxMC[2], correctedVertex - vtxMC[2]);

  const Int_t max = 10000;

  Bool_t used[max];
  for (Int_t i=0; i<max; ++i)
    used[i] = kFALSE;

  Int_t fakeTracks = 0;

  Int_t nTracks = fESD->GetNumberOfTracks();
  for (Int_t t=0; t<nTracks; t++)
  {
    AliESDtrack* esdTrack = fESD->GetTrack(t);

    if (!esdTrack)
      continue;


    if (!fEsdTrackCuts->AcceptTrack(esdTrack))
      continue;

    UInt_t status = esdTrack->GetStatus();
    if ((!status & AliESDtrack::kTPCrefit))
      continue;

    Int_t label = TMath::Abs(esdTrack->GetLabel());
    if (label < max)
    {
      if (used[label])
      {
        fFakes->Fill(1.5, esdTrack->GetLabel());
        fakeTracks++;
      }

      used[label] = kTRUE;
    }

    fFakes->Fill(0.5, esdTrack->GetLabel());
  }

  if (fakeTracks > 10)
    printf("In event %lld we have %d fake tracks.\n", entry, fakeTracks);

  return kTRUE;
}

void AliVertexSelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  AliSelectorRL::SlaveTerminate();

  fOutput->Add(fVertexMC);
  fOutput->Add(fVertexESD);
  fOutput->Add(fVertexCorr);
  fOutput->Add(fVertexCorr2);
  fOutput->Add(fFakes);
}

void AliVertexSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliSelectorRL::Terminate();

  fVertexMC = dynamic_cast<TH1F*> (fOutput->FindObject("fVertexMC"));
  fVertexESD = dynamic_cast<TH1F*> (fOutput->FindObject("fVertexESD"));
  fVertexCorr = dynamic_cast<TH2F*> (fOutput->FindObject("fVertexCorr"));
  fVertexCorr2 = dynamic_cast<TH2F*> (fOutput->FindObject("fVertexCorr2"));
  fFakes = dynamic_cast<TH2F*> (fOutput->FindObject("fFakes"));

  if (!fVertexMC || !fVertexESD || !fVertexCorr || !fVertexCorr2 || !fFakes)
    return;

  TFile* file = TFile::Open("vertex.root", "RECREATE");
  fVertexMC->Write();
  fVertexESD->Write();
  fVertexCorr->Write();
  fVertexCorr2->Write();
  fFakes->Write();
  file->Close();

  TCanvas* canvas = new TCanvas("canvas", "canvas", 1000, 1000);
  canvas->Divide(2, 2);

  canvas->cd(1);
  fVertexCorr->DrawCopy("COLZ");

  canvas->cd(2);
  fVertexCorr->ProfileX()->DrawCopy();

  canvas->cd(3);
  fVertexCorr2->DrawCopy("COLZ");

  canvas->cd(4);
  fFakes->DrawCopy();

  printf("%f tracks, %f with label > 0, %f with label == 0, %f with label < 0\n", fFakes->Integral(1, 1, 0, fFakes->GetYaxis()->GetNbins()), fFakes->Integral(1, 1, fFakes->GetYaxis()->FindBin(0.0), fFakes->GetYaxis()->GetNbins()), fFakes->GetBinContent(1, fFakes->GetYaxis()->FindBin(0.0)), fFakes->Integral(1, 1, 0, fFakes->GetYaxis()->FindBin(0.0)));
  printf("%f tracks, %f with label > 0, %f with label == 0, %f with label < 0\n", fFakes->Integral(2, 2, 0, fFakes->GetYaxis()->GetNbins()), fFakes->Integral(2, 2, fFakes->GetYaxis()->FindBin(0.0), fFakes->GetYaxis()->GetNbins()), fFakes->GetBinContent(2, fFakes->GetYaxis()->FindBin(0.0)), fFakes->Integral(2, 2, 0, fFakes->GetYaxis()->FindBin(0.0)));
}
