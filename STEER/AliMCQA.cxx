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

/*
$Log$
*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <strings.h>

#include "TObjArray.h"
#include "TH1.h"
#include "TList.h"
#include "TROOT.h"
#include "TBrowser.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TMath.h"

#include "AliMCQA.h"
#include "AliRun.h"
#include "AliModule.h"
#include "AliMC.h"

ClassImp(AliMCQA)


//_____________________________________________________________________________
AliMCQA::AliMCQA(Int_t ndets): fNdets(ndets)
{
  //
  // Constructor, creates the list of lists of histograms
  //
  TList *list;

  fQAList = new TObjArray(ndets);
  TObjArray &hist = *fQAList;

  char title[100];
  //
  TObjArray &mods = *(gAlice->Modules());
  AliModule *mod;
  TList *dir = gDirectory->GetList();
  for (Int_t i=0; i<ndets; i++) {
    hist[i] = list = new TList();
    mod = (AliModule *) mods[i];

    // Energy Spectrum
    sprintf(title,"Spectrum entering: %s ",mod->GetName());
    list->Add(new TH1F("hEnIn",strcat(title,mod->GetTitle()),100,-4,2));
    dir->Remove(dir->FindObject("hEnIn"));

    sprintf(title,"Spectrum exiting: %s ",mod->GetName());
    list->Add(new TH1F("hEnOut",strcat(title,mod->GetTitle()),100,-4,2));
    dir->Remove(dir->FindObject("hEnOut"));

    // Z position
    sprintf(title,"Z coordinate entering: %s ",mod->GetName());
    list->Add(new TH1F("hZIn",strcat(title,mod->GetTitle()),100,
		       mod->ZMin(),mod->ZMax()));
    dir->Remove(dir->FindObject("hZIn"));

    sprintf(title,"Z coordinate exiting: %s ",mod->GetName());
    list->Add(new TH1F("hZOut",strcat(title,mod->GetTitle()),100,
		       mod->ZMin(),mod->ZMax()));
    dir->Remove(dir->FindObject("hZOut"));
  }
  //
  gROOT->GetListOfBrowsables()->Add(this,"AliMCQA");

  fDetDone = new Int_t[fNdets];
}

//_____________________________________________________________________________
void AliMCQA::Browse(TBrowser *b)
{
  //
  // Called when the item "Run" is clicked on the left pane
  // of the Root browser.
  // It displays the Root Trees and all detectors.
  //

  TIter next(fQAList);
  TList *histos;
  TH1 *hist;
  while((histos = (TList*)next())) {
    TIter next1(histos);
    while((hist = (TH1*)next1())) {
      b->Add(hist,hist->GetTitle());
    }
  }
}

//_____________________________________________________________________________
void AliMCQA::PreTrack()
{
  fOldId=-1;
  for(Int_t i=0; i<fNdets; i++) fDetDone[i]=0;
}

//_____________________________________________________________________________
void AliMCQA::StepManager(Int_t id)
{
  if(fOldId != id) {
    TH1F *hist;
    TLorentzVector p, x;
    gMC->TrackMomentum(p);
    gMC->TrackPosition(x);
    Double_t energy = TMath::Max(
      p[3]-gAlice->PDGDB()->GetParticle(gMC->TrackPid())->Mass(),1.e-12);
    if(fOldId > -1) {
      if(!fDetDone[fOldId] && !gMC->IsNewTrack()) {
	TList *histold = (TList*) (*fQAList)[fOldId];
	hist = (TH1F*) histold->FindObject("hEnOut");
	hist->Fill(TMath::Log10(energy));
	hist = (TH1F*) histold->FindObject("hZOut");
	hist->Fill(x[2]);
	fDetDone[fOldId]=1;
      }
    }
    if(!fDetDone[id] && !gMC->IsNewTrack()) {
      TList *histnew = (TList*) (*fQAList)[id];
      hist = (TH1F*) histnew->FindObject("hEnIn");
      hist->Fill(TMath::Log10(energy));
      hist = (TH1F*) histnew->FindObject("hZIn");
      hist->Fill(x[2]);
    }
    fOldId=id;
  }
}
