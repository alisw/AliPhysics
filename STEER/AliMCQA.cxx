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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Quality Assurance class                                                   //
// Provides standard histograms for every element of the detector            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <strings.h>

#include <TBrowser.h>
#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TExec.h>
#include <TH1.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPaveLabel.h>
#include <TROOT.h>
#include <TVirtualMC.h>

#include "AliMC.h"
#include "AliMCQA.h"
#include "AliModule.h"
#include "AliRun.h"

ClassImp(AliMCQA)

//_______________________________________________________________________
AliMCQA::AliMCQA():
  fNdets(0),
  fNvolumes(0),
  fQAList(0),
  fOldId(0),
  fDetDone(0),
  fQAHist(0),
  fVolNames(0),
  fModNames(0),
  fMPaveLabel(0),
  fVPaveLabel(0)
{
  //
  // Default constructor
  //
}

//_______________________________________________________________________
AliMCQA::AliMCQA(const AliMCQA &qa):
  TObject(qa),
  fNdets(0),
  fNvolumes(0),
  fQAList(0),
  fOldId(0),
  fDetDone(0),
  fQAHist(0),
  fVolNames(0),
  fModNames(0),
  fMPaveLabel(0),
  fVPaveLabel(0)
{
  //
  // Copy constructor
  //
  qa.Copy(*this);
}

//_______________________________________________________________________
AliMCQA::AliMCQA(Int_t ndets):
  fNdets(ndets),
  fNvolumes(gMC->NofVolumes()),
  fQAList(new TObjArray(ndets)),
  fOldId(0),
  fDetDone(new Int_t[ndets]),
  fQAHist(new TObjArray(2)),
  fVolNames(new TObjArray(fNvolumes)),
  fModNames(new TObjArray(fNdets)),
  fMPaveLabel(0),
  fVPaveLabel(0)
{
  //
  // Constructor, creates the list of lists of histograms
  //
  TList *list;
  TH1F* h;
  Int_t i;
  
  TObjArray &hist = *fQAList;

  char title[100];
  //
  TObjArray &mods = *(gAlice->Modules());
  TList *dir = gDirectory->GetList();
  for (i=0; i<ndets; i++) {
    hist[i] = list = new TList();
    AliModule *mod = dynamic_cast<AliModule*>(mods[i]);

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

  //
  // Global QA histograms
  //
  
  fQAHist->Add(new TH1F("hMCVcalls","Monte Carlo calls per volume",
			fNvolumes, 0.5, fNvolumes+0.5));
  h = dynamic_cast<TH1F*>(dir->FindObject("hMCVcalls"));
  h->GetListOfFunctions()->Add(new TExec("ex","gAlice->GetMCQA()->AddVolumeName()"));
  dir->Remove(dir->FindObject("hMCVcalls"));
  //
  // Build list of volume names
  //
  AliMC * mc = gAlice->GetMCApp();
  for(i=0;i<fNvolumes;++i) {
    AliModule *mod = dynamic_cast<AliModule*>
      ((*gAlice->Modules())[mc->DetFromMate(gMC->VolId2Mate(i+1))]);
    (*fVolNames)[i]=new TNamed(gMC->VolName(i+1),mod->GetName());
  }

  fQAHist->Add(new TH1F("hMCMcalls","Monte Carlo calls per module",
			fNdets, -0.5, fNdets-0.5));
  h = dynamic_cast<TH1F*>(dir->FindObject("hMCMcalls"));
   h->GetListOfFunctions()->Add(new TExec("ex","gAlice->GetMCQA()->AddModuleName()"));

  dir->Remove(dir->FindObject("hMCMcalls"));
  //
  // Build list of module names
  //
  for(i=0;i<fNdets;++i) 
    (*fModNames)[i]=
      new TNamed((dynamic_cast<AliModule*>((*gAlice->Modules())[i]))->GetName(),"");
}

//_______________________________________________________________________
void AliMCQA::Copy(TObject &) const
{
  Fatal("Copy ctor","Not implemented!\n");
}

//_______________________________________________________________________
AliMCQA::~AliMCQA() 
 {
  //
  // Destructor
  //
  gROOT->GetListOfBrowsables()->Remove(this);
  //if program crashes here - it probobly means that 
  //one of added browsables was deleted and not removed previously from that list
  //skowron
  
  if (fQAList) {
    fQAList->Delete();
    delete fQAList;
    fQAList = 0;
  }
  if (fQAHist) {
    fQAHist->Delete();
    delete fQAHist;
    fQAHist = 0;
  }
  if (fVolNames) {
    fVolNames->Delete();
    delete fVolNames;
    fVolNames = 0;
  }
  if (fModNames) {
    fModNames->Delete();
    delete fModNames;
    fModNames = 0;
  }
  delete [] fDetDone;
  delete fMPaveLabel;
  delete fVPaveLabel;
}

//_______________________________________________________________________
void AliMCQA::Browse(TBrowser *b)
{
  //
  // Called when the item "Run" is clicked on the left pane
  // of the Root browser.
  // It displays the Root Trees and all detectors.
  //
  TH1 *hist;
  //
  // Global histos first
  //
  TIter global(fQAHist);
  while((hist = dynamic_cast<TH1*>(global())))
    b->Add(hist,hist->GetTitle());
  //
  // Module histograms now
  //
  TIter next(fQAList);
  TList *histos;
  while((histos = dynamic_cast<TList*>(next()))) {
    TIter next1(histos);
    while((hist = dynamic_cast<TH1*>(next1())))
      b->Add(hist,hist->GetTitle());
  }
}

//_______________________________________________________________________
void AliMCQA::PreTrack()
{
  //
  // Called before each track
  //
  fOldId=-1;
  for(Int_t i=0; i<fNdets; i++) fDetDone[i]=0;
}

//_______________________________________________________________________
void AliMCQA::StepManager(Int_t id)
{
  //
  // Called at each step
  //
  TH1F *hist;
  Int_t copy;
  //
  // Fill Global histograms first
  //


  static TH1F* mcvcalls = dynamic_cast<TH1F*>(fQAHist->FindObject("hMCVcalls"));
  mcvcalls->Fill(gMC->CurrentVolID(copy));
  static TH1F* mcmcalls = dynamic_cast<TH1F*>(fQAHist->FindObject("hMCMcalls"));
  mcmcalls->Fill(id);

  //
  // Now the step manager histograms
  //
  if(fOldId != id) {
    static  Double_t mpi0=0;
    static  Double_t mpip=0;
    static  Double_t mpim=0;
    static  Double_t mep=0;
    static  Double_t mem=0;
    Double_t mass = 0;
    Int_t num = gMC->TrackPid();

    switch (num) {
    case 111: 
      if (mpi0==0) mpi0=gAlice->PDGDB()->GetParticle(num)->Mass(); 
      mass=mpi0;
      break;
    case 211:
      if (mpip==0) mpip=gAlice->PDGDB()->GetParticle(num)->Mass(); 
      mass=mpip;
      break;
    case -211:
      if (mpim==0) mpim=gAlice->PDGDB()->GetParticle(num)->Mass(); 
      mass=mpim;
      break;
    case 11:
      if (mep==0) mep=gAlice->PDGDB()->GetParticle(num)->Mass(); 
      mass=mep;
      break;
    case -11:
      if (mem==0) mem=gAlice->PDGDB()->GetParticle(num)->Mass(); 
      mass=mem;
      break;
    default:
	if (gAlice->PDGDB()->GetParticle(num))
	    mass = gAlice->PDGDB()->GetParticle(num)->Mass();
      break; 
    }

    static TLorentzVector p, x;
    gMC->TrackMomentum(p);
    gMC->TrackPosition(x);
    Double_t energy = TMath::Max(p[3]-mass,1.e-12);
    if(fOldId > -1) {
      if(!fDetDone[fOldId] && !gMC->IsNewTrack()) {
	TList *histold = dynamic_cast<TList*>((*fQAList)[fOldId]);
	hist = dynamic_cast<TH1F*>(histold->FindObject("hEnOut"));
	hist->Fill(TMath::Log10(energy));
	hist = dynamic_cast<TH1F*>(histold->FindObject("hZOut"));
	hist->Fill(x[2]);
	fDetDone[fOldId]=1;
      }
    }
    if(!fDetDone[id] && !gMC->IsNewTrack()) {
      TList *histnew = dynamic_cast<TList*>((*fQAList)[id]);
      hist = dynamic_cast<TH1F*>(histnew->FindObject("hEnIn"));
      hist->Fill(TMath::Log10(energy));
      hist = dynamic_cast<TH1F*>(histnew->FindObject("hZIn"));
      hist->Fill(x[2]);
    }
    fOldId=id;
  }
}

//_______________________________________________________________________
void AliMCQA::AddModuleName()
{
  //
  // Add function DrawModuleName to the module frequency histogram pad
  //
  static TVirtualPad* oldpad=0;
  if(oldpad!=gPad) {
    gPad->GetCanvas()->FeedbackMode(kTRUE);
    if(gPad) gPad->AddExec("ex","gAlice->GetMCQA()->DrawModuleName()");
    DrawPaveLabel(fMPaveLabel);
    oldpad=gPad;
  }
}

//_______________________________________________________________________
void AliMCQA::AddVolumeName()
{
  //
  // Add function DrawVolumeName to the volume frequency histogram pad
  //
  static TVirtualPad* oldpad=0;
  if(oldpad!=gPad) {
    gPad->GetCanvas()->FeedbackMode(kTRUE);
    if(gPad) gPad->AddExec("ex","gAlice->GetMCQA()->DrawVolumeName()");
    DrawPaveLabel(fVPaveLabel);
    oldpad=gPad;
  }
}

//_______________________________________________________________________
void AliMCQA::DrawPaveLabel(TPaveLabel *&pv)
{
  //
  // Draws the PaveLabel with the meaning of the bin
  //
  float uxmin = gPad->GetUxmin();
  float uxmax = gPad->GetUxmax();
  float uymin = gPad->GetUymin();
  float uymax = gPad->GetUymax();
  float lx = uxmax-uxmin;
  float ly = uymax-uymin;

  if(pv) delete pv;
  pv 
    = new TPaveLabel(uxmin+0.05*lx,uymax-(0.05+0.1)*ly,
		     uxmin+(0.05+0.2)*lx,uymax-0.05*ly,
		     "");
  pv->Draw();
}

//_______________________________________________________________________
Int_t AliMCQA::GetHBin(const char* hname)
{
  //
  // Get the bin where the cursor is
  //
  TList *dir = gDirectory->GetList();
  TH1 *h=dynamic_cast<TH1*>(dir->FindObject(hname));
  

  int px = gPad->GetEventX();
  Float_t upx = gPad->AbsPixeltoX(px);
  Float_t x = gPad->PadtoX(upx);
    
  return h->GetXaxis()->FindBin(x);
}

//_______________________________________________________________________
void AliMCQA::DrawModuleName()
{
  //
  // Writes the name of the module of the bin where the cursor is
  //
  TObject *select = gPad->GetSelected();
  if(!select) return;
  
  Int_t binx = GetHBin("hMCMcalls");
  if(0<binx && binx<=fNdets) {
    char lab[15];
    strcpy(lab,dynamic_cast<TNamed*>((*fModNames)[binx-1])->GetName());
    fMPaveLabel->SetLabel(lab);
  
    gPad->Modified();
    gPad->Update();
  }
}

//_______________________________________________________________________
void AliMCQA::DrawVolumeName()
{
  //
  // Writes the name of the volume:module of the bin where the cursor is
  //
  TObject *select = gPad->GetSelected();
  if(!select) return;

  Int_t binx = GetHBin("hMCVcalls");
  if(0<binx && binx<=fNvolumes) {
    char lab[20];
    sprintf(lab,"%s: %s",dynamic_cast<TNamed*>((*fVolNames)[binx-1])->GetName(),
	    dynamic_cast<TNamed*>((*fVolNames)[binx-1])->GetTitle());
    fVPaveLabel->SetLabel(lab);
    
    gPad->Modified();
    gPad->Update();
  }
}


