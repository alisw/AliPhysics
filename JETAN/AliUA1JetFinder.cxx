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
// UA1 Jet finder 
// manages the search for jets 
// Author: jgcn@mda.cinvestav.mx
// (code adapted from EMCAL directory)
//---------------------------------------------------------------------

#include <Riostream.h>
#include <TLorentzVector.h>
#include <TH2F.h>
#include "AliUA1JetFinder.h"
#include "AliUA1JetHeader.h"
#include "UA1Common.h"
#include "AliJetReaderHeader.h"
#include "AliJetReader.h"
#include "AliJet.h"


ClassImp(AliUA1JetFinder)

////////////////////////////////////////////////////////////////////////

AliUA1JetFinder::AliUA1JetFinder()

{
  //
  // Constructor
  //
  fHeader = 0x0;
  fLego   = 0x0;
}

////////////////////////////////////////////////////////////////////////

AliUA1JetFinder::~AliUA1JetFinder()

{
  //
  // destructor
  //

  // delete fLego;            
  
  // reset and delete header
}

////////////////////////////////////////////////////////////////////////

#ifndef WIN32
# define ua1_jet_finder ua1_jet_finder_
# define hf1 hf1_
# define type_of_call

#else
# define ua1_jet_finder UA1_JET_FINDER
# define hf1 HF1
# define type_of_call _stdcall
#endif

extern "C" void type_of_call
ua1_jet_finder(Int_t& ncell, Int_t& ncell_tot,
               Float_t  etc[60000],  Float_t etac[60000],
               Float_t  phic[60000],
               Float_t& min_move, Float_t& max_move, Int_t& mode,
               Float_t& prec_bg,  Int_t& ierror);

extern "C" void type_of_call hf1(Int_t& id, Float_t& x, Float_t& wgt);

////////////////////////////////////////////////////////////////////////

void AliUA1JetFinder::FindJets()

{
  // initialize event, look for jets, download jet info


  // initialization in 2 steps
  // 1) transform input to pt,eta,phi plus lego
  // 2) dump lego
  
  // 1) transform input to pt,eta,phi plus lego
  TClonesArray *lvArray = fReader->GetMomentumArray();
  Int_t nIn =  lvArray->GetEntries();
  if (nIn == 0) return;
    
  // local arrays for input
  Float_t* ptT  = new Float_t[nIn];
  Float_t* etaT = new Float_t[nIn];
  Float_t* phiT = new Float_t[nIn];

  // load input vectors
  for (Int_t i = 0; i < nIn; i++){
      TLorentzVector *lv = (TLorentzVector*) lvArray->At(i);
      ptT[i]  = lv->Pt();
      etaT[i] = lv->Eta();
      phiT[i] = ((lv->Phi() < 0) ? (lv->Phi()) + 2 * TMath::Pi() : lv->Phi());
      fLego->Fill(etaT[i], phiT[i], ptT[i]);
  }
  fJets->SetNinput(nIn);
  
  // 2) dump lego
  // check enough space! *to be done*
  Float_t etCell[60000];   //! Cell Energy
  Float_t etaCell[60000];  //! Cell eta
  Float_t phiCell[60000];  //! Cell phi

  Int_t nCell = 0;
  TAxis* xaxis = fLego->GetXaxis();
  TAxis* yaxis = fLego->GetYaxis();
  Float_t e = 0.0;
  for (Int_t i = 1; i <= fHeader->GetNbinEta(); i++) {
      for (Int_t j = 1; j <= fHeader->GetNbinPhi(); j++) {
	  e = fLego->GetBinContent(i,j);
	  if (e > 0.0) e -= fHeader->GetMinCellEt();
	  if (e < 0.0) e = 0.;
	  Float_t eta  = xaxis->GetBinCenter(i);
	  Float_t phi  = yaxis->GetBinCenter(j);	    
	  etCell[nCell]  = e;
	  etaCell[nCell] = eta;
	  phiCell[nCell] = phi;
	  nCell++;
      }
  }

  // run the algo. Parameters from header
  Int_t nTot      = (fHeader->GetNbinEta())*(fHeader->GetNbinPhi());
  Float_t minmove = fHeader->GetMinMove();
  Float_t maxmove = fHeader->GetMaxMove();
  Int_t mode      = fHeader->GetMode();
  Float_t precbg  = fHeader->GetPrecBg();
  Int_t ierror;
  ua1_jet_finder(nCell, nTot, etCell, etaCell, phiCell, 
		 minmove, maxmove, mode, precbg, ierror);


  // sort jets
  Int_t * idx  = new Int_t[UA1JETS.njet];
  TMath::Sort(UA1JETS.njet, UA1JETS.etj, idx);

  // download jet info.   
  AliJetReaderHeader* rheader = fReader->GetReaderHeader();
  for(Int_t i = 0; i < UA1JETS.njet; i++) {
    // reject events outside acceptable eta range
    if (((UA1JETS.etaj[1][idx[i]])> (rheader->GetFiducialEtaMax()))
	|| ((UA1JETS.etaj[1][idx[i]]) < (rheader->GetFiducialEtaMin())))
      continue;
    Float_t px, py,pz,en; // convert to 4-vector
    px = UA1JETS.etj[idx[i]] * TMath::Cos(UA1JETS.phij[1][idx[i]]);
    py = UA1JETS.etj[idx[i]] * TMath::Sin(UA1JETS.phij[1][idx[i]]);
    pz = UA1JETS.etj[idx[i]] /
      TMath::Tan(2.0 * TMath::ATan(TMath::Exp(-UA1JETS.etaj[1][idx[i]])));
    en = TMath::Sqrt(px * px + py * py + pz * pz);
    fJets->AddJet(px, py, pz, en);
  }
  
  // find multiplicities and relationship jet-particle
  // find also percentage of pt from pythia
  Int_t* injet = new Int_t[nIn];
  for (Int_t i = 0; i < nIn; i++) injet[i]= -1;
  Int_t* mult  = new Int_t[fJets->GetNJets()];
  Float_t* percentage  = new Float_t[fJets->GetNJets()];

  for(Int_t i = 0; i < (fJets->GetNJets()); i++) {
    Float_t pt_sig = 0.0;
    mult[i] = 0;
    for (Int_t j = 0; j < nIn; j++) {
      Float_t deta = etaT[j] - fJets->GetEta(i);
      Float_t dphi = phiT[j] - fJets->GetPhi(i);
      if (dphi < -TMath::Pi()) dphi= -dphi - 2.0 * TMath::Pi();
      if (dphi > TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
      Float_t dr = TMath::Sqrt(deta * deta + dphi * dphi);
      if (dr < fHeader->GetRadius() && injet[j] == -1) {
	if (fReader->GetSignalFlag(j) == 1) pt_sig+=ptT[j];
	injet[j] = i;
	mult[i]++;
      }
    }
    percentage[i] = pt_sig/((Double_t) fJets->GetPt(i));    
  }
  fJets->SetPtFromSignal(percentage);
  fJets->SetMultiplicities(mult);
  fJets->SetInJet(injet);
}

////////////////////////////////////////////////////////////////////////

void AliUA1JetFinder::Reset()

{
  fLego->Reset();
  fJets->ClearJets();
}

////////////////////////////////////////////////////////////////////////

void AliUA1JetFinder::WriteJHeaderToFile()
{
  fOut->cd();
  fHeader->Write();
}


////////////////////////////////////////////////////////////////////////

void AliUA1JetFinder::Init()
{
  // initializes some variables
  Float_t dEta, dPhi;
  dEta = (fHeader->GetEtaMax()-fHeader->GetEtaMin())
    /((Float_t) fHeader->GetNbinEta());
  dPhi = (fHeader->GetPhiMax()-fHeader->GetPhiMin())
    /((Float_t) fHeader->GetNbinPhi());

  UA1CELL.etaCellSize = dEta;
  UA1CELL.phiCellSize = dPhi;

  // jet parameters
  UA1PARA.coneRad = fHeader->GetRadius();
  UA1PARA.etSeed  = fHeader->GetEtSeed();
  UA1PARA.ejMin   = fHeader->GetMinJetEt();
  UA1PARA.etMin   = fHeader->GetMinCellEt();

  // book lego
  fLego = new 
    TH2F("legoH","eta-phi",
	 fHeader->GetNbinEta(), fHeader->GetEtaMin(), fHeader->GetEtaMax(),
	 fHeader->GetNbinPhi(), fHeader->GetPhiMin(), fHeader->GetPhiMax());
}
