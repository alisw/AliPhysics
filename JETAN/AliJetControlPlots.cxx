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
  
//---------------------------------------------------------------------
// Jet Control Plots class 
// manages histograms with control plots of jet searching
// Stores the output of a jet algorithm
// Author: jgcn@mda.cinvestav.mx
//---------------------------------------------------------------------

#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TFile.h> 
#include <TClonesArray.h>
#include <TH1I.h>
#include <TH1D.h>
#include "AliJetReader.h"
#include "AliJet.h"

#include "AliJetControlPlots.h"
ClassImp(AliJetControlPlots)
  
////////////////////////////////////////////////////////////////////////

AliJetControlPlots::AliJetControlPlots()
{
  //
  // Constructor
  //
  fNJetT=0;

  fNJetsH = new TH1I("fNJetsH","Number of Jets",12,0,11);
  SetProperties(fNJetsH,"Number of jets","Entries");

  fMultH = new TH1I("fMultH","Multiplicity of Jets",22,0,21);
  SetProperties(fMultH,"Multiplicity of jets","Entries");

  fPtH = new TH1D("fPtH","Pt of Jets",50,0.,200.);
  SetProperties(fPtH,"P_{t} [GeV]","Entries");

  fEtaH = new TH1D("fEtaH","Pseudorapidity of Jets",30,-1.5,1.5);
  SetProperties(fEtaH,"#eta","Entries");

  fEneH = new TH1D("fEneH","Energy of Jets",50,0.,200.);
  SetProperties(fEneH,"Energy [GeV]","Entries");

  fFragH = new TH1D("fFragH","Jet Fragmentation",20,0.,1.);
  SetProperties(fFragH,"x=E_{i}/E_{JET}","1/N_{JET}dN_{ch}/dx");

  fFragLnH = new TH1D("fFragLnH","Jet Fragmentation",20,0.,10);
  SetProperties(fFragLnH,"#xi=ln(E_{JET}/E_{i}","1/N_{JET}dN_{ch}/d#xi");

  fPhiH = new TH1D("fPhiH","Azimuthal angle of Jets",
		   60,-TMath::Pi(),TMath::Pi());
  SetProperties(fPhiH,"#phi","Entries");
  
  fInJetH = new  TH1D("fInJetH","Percentage of particles in jets",
		      50,0,1);
  SetProperties(fInJetH,"percentage of particles in jets","Entries");
}

////////////////////////////////////////////////////////////////////////

AliJetControlPlots::~AliJetControlPlots()
{
  //
  // Destructor
  //
  delete fNJetsH;
  delete fMultH;
  delete fPtH;
  delete fEtaH;
  delete fEneH;
  delete fFragH;
  delete fFragLnH;
  delete fPhiH;
  delete fInJetH;
}

////////////////////////////////////////////////////////////////////////

void AliJetControlPlots::FillHistos(AliJet *j, AliJetReader *r)
{
// Fills the histograms

  Int_t nj = j->GetNJets();
  fNJetsH->Fill(nj);
  fNJetT+=nj;

  // kinematics, occupancy and multiplicities
  TArrayI mj = j->GetMultiplicities();
  Int_t mjTot=0;
  for (Int_t i=0;i<nj;i++) {
    mjTot+=mj[i];
    fMultH->Fill(mj[i]);    
    fPtH->Fill(j->GetPt(i));
    fEtaH->Fill(j->GetEta(i));
    fEneH->Fill(j->GetE(i));
    fPhiH->Fill(j->GetPhi(i));
  }
  if (nj>0) 
    fInJetH->Fill(((Double_t) mjTot)/((Double_t) j->GetNinput()));

  // fragmentation 
  TClonesArray *lvArray = r->GetMomentumArray();
  TArrayI inJet = j->GetInJet();
  Int_t nIn = j->GetNinput();
  for(Int_t i=0;i<nIn;i++) {
    if (inJet[i] == -1) continue;
    TLorentzVector *lv = (TLorentzVector*) lvArray->At(i);
    Double_t xe = (lv->E())/(j->GetE(inJet[i]));
    fFragH->Fill(xe);
    fFragLnH->Fill(TMath::Log(1.0/xe));
  }
}

////////////////////////////////////////////////////////////////////////

void AliJetControlPlots::Normalize()
{
  if (fNJetT == 0) return;
  fFragH->Sumw2();
  fFragH->Scale(20.0/((Double_t) fNJetT));
  fFragLnH->Sumw2();
  fFragLnH->Scale(2.0/((Double_t) fNJetT));
}

////////////////////////////////////////////////////////////////////////

void AliJetControlPlots::PlotHistos()
{
  // style
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptTitle(kFALSE);
  
  TCanvas *c = new TCanvas("c","AliJetControlPlots",50,200,900,700);
  c->Divide(3,3);
  c->cd(1);  fNJetsH->Draw("e1");
  c->cd(2);  fMultH->Draw("e1");
  c->cd(3);  fInJetH->Draw("e1");
  c->cd(4);  fPtH->Draw("e1");
  c->cd(5);  fEtaH->Draw("e1");
  c->cd(6);  fPhiH->Draw("e1");
  c->cd(7);  fEneH->Draw("e1");
  c->cd(8);  fFragH->Draw("e1");
  c->cd(9);  fFragLnH->Draw("e1");
}

////////////////////////////////////////////////////////////////////////

void AliJetControlPlots::SetProperties(TH1* h,const char* x, const char* y) const
{
//
// Sets the histogram style properties
  h->SetMarkerStyle(20);
  h->SetMarkerSize(.5);
  h->SetMarkerColor(2);
  h->SetXTitle(x);
  h->SetYTitle(y);
}


