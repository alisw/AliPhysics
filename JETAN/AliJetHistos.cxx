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

#include <TList.h>
#include <TClonesArray.h>
#include <TH1I.h>
#include <TH1F.h>
#include <TMath.h>

#include "AliAODJet.h"
#include "AliJetHistos.h"

ClassImp(AliJetHistos)

AliJetHistos::AliJetHistos():
  fNJetsH(0x0),
  fPtH(0x0),
  fEtaH(0x0),
  fEneH(0x0),
  fPhiH(0x0)
{
  // Default constructor

  fNJetsH = new TH1I("NJetsH","Number of Jets",12,0,11);
  SetProperties(fNJetsH,"Number of jets","Entries");

  fPtH = new TH1F("PtH","Pt of Jets",50,0.,200.);
  SetProperties(fPtH,"P_{#perp} [GeV]","Entries");

  fEtaH = new TH1F("EtaH","Pseudorapidity of Jets",30,-1.5,1.5);
  SetProperties(fEtaH,"#eta","Entries");

  fEneH = new TH1F("EneH","Energy of Jets",50,0.,200.);
  SetProperties(fEneH,"Energy [GeV]","Entries");

  fPhiH = new TH1F("PhiH","Azimuthal angle of Jets",
                   60,0.,2.0*TMath::Pi());
  SetProperties(fPhiH,"#phi","Entries");
}

AliJetHistos::~AliJetHistos()
{
// Destructor
  delete fNJetsH;
  delete fPtH;
  delete fEtaH;
  delete fEneH;
  delete fPhiH;
}

void AliJetHistos::SetProperties(TH1* h,const char* x, const char* y) const
{
// Sets the histogram style properties
  h->SetMarkerStyle(20);
  h->SetMarkerSize(.5);
  h->SetMarkerColor(2);
  h->SetXTitle(x);
  h->SetYTitle(y);
  h->Sumw2();
}

void AliJetHistos::AddHistosToList(TList *list)
{
// Add histos to the list
  list->Add(fNJetsH);
  list->Add(fPtH);
  list->Add(fEtaH);
  list->Add(fEneH);
  list->Add(fPhiH);
}


void AliJetHistos::FillHistos(TClonesArray *jets)
{
// Fill histograms
  if(!jets)return;
  Int_t nj = jets->GetEntries();
  fNJetsH->Fill(nj,1);

  if (nj == 0 ) return;

  AliAODJet *j;
  for (Int_t i=0;i<nj;i++) {
    j = (AliAODJet *) jets->At(i);
    fPtH->Fill(j->Pt(),1);
    fEtaH->Fill(j->Eta(),1);
    fEneH->Fill(j->E(),1);
    fPhiH->Fill(j->Phi(),1);
  }
  
}
