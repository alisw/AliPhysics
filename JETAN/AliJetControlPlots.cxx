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

//---------------------------------------------------------------------
// Jet Control Plots class 
// manages histograms with control plots of jet searching
// Stores the output of a jet algorithm
// Author: jgcn@mda.cinvestav.mx
//---------------------------------------------------------------------

#include <TCanvas.h>
#include <TH1D.h>
#include <TH1I.h>
#include <TMath.h>
#include <TStyle.h>

#include "AliJetReader.h"
#include "AliJet.h"

#include "AliJetControlPlots.h"
ClassImp(AliJetControlPlots)
  
////////////////////////////////////////////////////////////////////////

AliJetControlPlots::AliJetControlPlots():
  fNJetsH(0),
  fMultH(0),
  fPtH(0),
  fEtaH(0),
  fEneH(0),
  fFragH(0),
  fFragLnH(0),
  fFragrH(0),
  fFragLnrH(0),
  fShapeH(0),
  fShaperH(0),
  fPhiH(0),
  fInJetH(0),
  fNJetT(0)

{
  // Default constructor

  // general properties
  fNJetsH = new TH1I("fNJetsH","Number of Jets",12,0,11);
  SetProperties(fNJetsH,"Number of jets","Entries");

  fMultH = new TH1I("fMultH","Multiplicity of Jets",22,0,21);
  SetProperties(fMultH,"Multiplicity of jets","Entries");

  fInJetH = new  TH1D("fInJetH","Percentage of particles in jets",50,0,1);
  SetProperties(fInJetH,"percentage of particles in jets","Entries");

  // kinematics
  fPtH = new TH1D("fPtH","Pt of Jets",50,0.,200.);
  SetProperties(fPtH,"P_{t} [GeV]","Entries");

  fEtaH = new TH1D("fEtaH","Pseudorapidity of Jets",30,-1.5,1.5);
  SetProperties(fEtaH,"#eta","Entries");

  fEneH = new TH1D("fEneH","Energy of Jets",50,0.,200.);
  SetProperties(fEneH,"Energy [GeV]","Entries");

  fPhiH = new TH1D("fPhiH","Azimuthal angle of Jets",
		   60,0.,2.0*TMath::Pi());
  SetProperties(fPhiH,"#phi","Entries");

  // fragmentation 
  fFragH = new TH1D("fFragH","Leading Jet Fragmentation (selected tracks)",
		    20,0.,1.);
  SetProperties(fFragH,"x=E_{i}/E_{JET}","1/N_{JET}dN_{ch}/dx");

  fFragrH = new TH1D("fFragrH","Leading Jet Fragmentation (rejected tracks)",
		    20,0.,1.);
  SetProperties(fFragrH,"x=E_{i}/E_{JET}","1/N_{JET}dN_{ch}/dx");

  // fragmentation log
  fFragLnH = new TH1D("fFragLnH","Leading Jet Fragmentation (selected tracks)",
		      20,0.,10);
  SetProperties(fFragLnH,"#xi=ln(E_{JET}/E_{i})","1/N_{JET}dN_{ch}/d#xi");

  fFragLnrH = new TH1D("fFragLnrH",
		       "Leading Jet Fragmentation (rejected tracks)",
		       20,0.,10);
  SetProperties(fFragLnrH,"#xi=ln(E_{JET}/E_{i})","1/N_{JET}dN_{ch}/d#xi");

  // jet shape
  fShapeH = new TH1D("fShapeH","Leading Jet Shape (selected tracks)",
		     20,0.,1.);
  SetProperties(fShapeH,"r","1/N_{JET}#Sigma P_{T}(0,r)/P_{T}(0,R)");
      
  fShaperH = new TH1D("fShaperH","Leading Jet Shape (rejected tracks)",
		     20,0.,1.);
  SetProperties(fShaperH,"r","1/N_{JET}#Sigma P_{T}(0,r)/P_{T}(0,R)");
}

////////////////////////////////////////////////////////////////////////

AliJetControlPlots::~AliJetControlPlots()
{
  // Destructor
  delete fNJetsH;
  delete fMultH;
  delete fPtH;
  delete fEtaH;
  delete fEneH;
  delete fPhiH;
  delete fInJetH;
  delete fFragH;
  delete fFragLnH;
  delete fFragrH;
  delete fFragLnrH;
  delete fShapeH;
  delete fShaperH;
}

////////////////////////////////////////////////////////////////////////

void AliJetControlPlots::FillHistos(AliJet *j)
{
  // Fills the histograms
  Int_t nj = j->GetNJets();
  fNJetsH->Fill(nj,1);
  if (nj == 0) return;
  
  // kinematics, occupancy and multiplicities
  TArrayI mj = j->GetMultiplicities();
  Int_t mjTot=0;
  for (Int_t i=0;i<nj;i++) {
    mjTot+=mj[i];
    fMultH->Fill(mj[i],1);    
    fPtH->Fill(j->GetPt(i),1);
    fEtaH->Fill(j->GetEta(i),1);
    fEneH->Fill(j->GetE(i),1);
    fPhiH->Fill(j->GetPhi(i),1);
  }
  fInJetH->Fill(((Double_t) mjTot)/((Double_t) j->GetNinput()),1);
  
  // fragmentation of leading jet 
  TArrayI inJet = j->GetInJet();
  TArrayF etain = j->GetEtaIn();
  TArrayF ptin = j->GetPtIn();
  for(Int_t i=0;i<(inJet.GetSize());i++) {
    Float_t t1 = TMath::Tan(2.0*TMath::ATan(TMath::Exp(-etain[i])));
    Float_t ene = ptin[i]*TMath::Sqrt(1.+1./(t1*t1));
    if (inJet[i] == 1) {
      fFragH->Fill((Float_t) ene/(j->GetE(0)),1);
      fFragLnH->Fill((Float_t) TMath::Log((j->GetE(0))/ene),1);
    }
    if (inJet[i] == -1) {
      fFragrH->Fill((Float_t) ene/(j->GetE(0)),1);
      fFragLnrH->Fill((Float_t) TMath::Log((j->GetE(0))/ene),1);
    }
  }
  
  // shape of leading jet 
  // * CAREFUL: depends on number of bins and bin sizes
  //   of shape histo. HARDWIRED at the moment *
  //   Plot against dr NOT dr/R_jet!!!
  TArrayF phiin = j->GetPhiIn();
  Float_t rptS[20], rptR[20];
  for(Int_t i=0;i<20;i++) rptS[i]=rptR[i]=0.0;
  
  for(Int_t i=0;i<inJet.GetSize();i++) {
    if (inJet[i] == 1 || inJet[i] == -1) {
      Float_t deta = etain[i] - j->GetEta(0);
      Float_t dphi = phiin[i] - j->GetPhi(0);
      if (dphi < -TMath::Pi()) dphi= -dphi - 2.0 * TMath::Pi();
      if (dphi > TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
      Float_t dr = TMath::Sqrt(deta * deta + dphi * dphi);
      if (dr>1) continue;
      Int_t bin = (Int_t) TMath::Floor(dr/0.05);
      if (inJet[i] == 1) rptS[bin]+=ptin[i]/(j->GetPt(0));
      if (inJet[i] == -1) rptR[bin]+=ptin[i]/(j->GetPt(0));	  
    }
  }
  
  Float_t ptS,ptR,r;
  ptS=ptR=0.0;
  for (Int_t i=0;i<20;i++) {
    ptS+=rptS[i];
    ptR+=rptR[i];
    r=(i+1)*0.05-0.025;
    fShapeH->Fill(r,ptS);
    fShaperH->Fill(r,ptR);      
  }
  fNJetT++;
}

////////////////////////////////////////////////////////////////////////

void AliJetControlPlots::Normalize()
{
// *CAREFUL: depends on histogram number of bins 
  if (fNJetT == 0) return;
  fFragH->Scale(20.0/((Double_t) fNJetT));
  fFragLnH->Scale(2.0/((Double_t) fNJetT));
  fFragrH->Scale(20.0/((Double_t) fNJetT));
  fFragLnrH->Scale(2.0/((Double_t) fNJetT));
  fShapeH->Scale(1.0/((Double_t) fNJetT));
  fShaperH->Scale(1.0/((Double_t) fNJetT));
}

////////////////////////////////////////////////////////////////////////

void AliJetControlPlots::PlotHistos()
{
  // style
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptTitle(kFALSE);
  
  TCanvas *c = new TCanvas("c","AliJetControlPlots",50,200,1200,700);
  c->Divide(4,3);
  c->cd(1);  fNJetsH->Draw("e1");
  c->cd(2);  fMultH->Draw("e1");
  c->cd(3);  fInJetH->Draw("e1");
  c->cd(4);  fPtH->Draw("e1");
  c->cd(5);  fEtaH->Draw("e1");
  c->cd(6);  fPhiH->Draw("e1");
  c->cd(7);  fEneH->Draw("e1");
  c->cd(8);  fFragLnH->Draw("e1");
  c->cd(9);  fFragLnrH->Draw("e1");
  c->cd(10);  fShapeH->Draw("e1");
  c->cd(11);  fShaperH->Draw("e1");
}

////////////////////////////////////////////////////////////////////////

void AliJetControlPlots::SetProperties(TH1* h,const char* x, const char* y) const
{
// Sets the histogram style properties
  h->SetMarkerStyle(20);
  h->SetMarkerSize(.5);
  h->SetMarkerColor(2);
  h->SetXTitle(x);
  h->SetYTitle(y);
  // h->Sumw2();
}


