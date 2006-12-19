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
#include <TFile.h>
#include <TH2F.h>
#include "AliUA1JetFinder.h"
#include "AliUA1JetHeader.h"
#include "UA1Common.h"
#include "AliJetReaderHeader.h"
#include "AliJetReader.h"
#include "AliJet.h"


ClassImp(AliUA1JetFinder)

////////////////////////////////////////////////////////////////////////

AliUA1JetFinder::AliUA1JetFinder():
  fHeader(0x0),
  fLego(0x0)
{
  // Default constructor
}

////////////////////////////////////////////////////////////////////////

AliUA1JetFinder::~AliUA1JetFinder()

{
  // destructor
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

void AliUA1JetFinder::FindJetsTPC()

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
  Float_t* enT  = new Float_t[nIn];
  Float_t* ptT  = new Float_t[nIn];
  Float_t* etaT = new Float_t[nIn];
  Float_t* phiT = new Float_t[nIn];

  // load input vectors
  for (Int_t i = 0; i < nIn; i++){
    TLorentzVector *lv = (TLorentzVector*) lvArray->At(i);
    enT[i]  = lv->Energy();
    ptT[i]  = lv->Pt();
    etaT[i] = lv->Eta();
    phiT[i] = ((lv->Phi() < 0) ? (lv->Phi()) + 2 * TMath::Pi() : lv->Phi());
    if (fReader->GetCutFlag(i) == 1) 
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
  for (Int_t i = 1; i <= fHeader->GetLegoNbinEta(); i++) {
      for (Int_t j = 1; j <= fHeader->GetLegoNbinPhi(); j++) {
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
  Int_t nTot      = (fHeader->GetLegoNbinEta())*(fHeader->GetLegoNbinPhi());
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
  for(Int_t i = 0; i < UA1JETS.njet; i++) {
    // reject events outside acceptable eta range
    if (((UA1JETS.etaj[1][idx[i]])> (fHeader->GetJetEtaMax()))
	|| ((UA1JETS.etaj[1][idx[i]]) < (fHeader->GetJetEtaMin())))
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
  Int_t* sflag = new Int_t[nIn];
  for (Int_t i = 0; i < nIn; i++) {injet[i]= 0;sflag[i]=0;}
  Int_t* mult  = new Int_t[fJets->GetNJets()];
  Int_t* ncell  = new Int_t[fJets->GetNJets()];
  Float_t* percentage  = new Float_t[fJets->GetNJets()];

  for(Int_t i = 0; i < (fJets->GetNJets()); i++) {
      Float_t pt_sig = 0.0;
      mult[i] = 0;
      ncell[i] = UA1JETS.ncellj[i];
      for (Int_t j = 0; j < nIn; j++) {
	  Float_t deta = etaT[j] - fJets->GetEta(i);
	  Float_t dphi = phiT[j] - fJets->GetPhi(i);
	  if (dphi < -TMath::Pi()) dphi= -dphi - 2.0 * TMath::Pi();
	  if (dphi > TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
	  Float_t dr = TMath::Sqrt(deta * deta + dphi * dphi);
	  if (dr < fHeader->GetRadius() && injet[j] == 0) {
	      injet[j] = -(i+1);
	      if ((fReader->GetCutFlag(j)) == 1 &&
		  (etaT[j] < fHeader->GetLegoEtaMax()) &&
		  (etaT[j] > fHeader->GetLegoEtaMin())) {
		  injet[j] = i+1;
		  mult[i]++;
		  if (fReader->GetSignalFlag(j) == 1) {
		    pt_sig+=ptT[j];
		    sflag[j]=1;
		  }
	      }
	  }
      }
      percentage[i] = (pt_sig-ncell[i]*UA1JETS.etavg)/
	((Double_t) fJets->GetPt(i));    
  }
  fJets->SetNCells(ncell);
  fJets->SetPtFromSignal(percentage);
  fJets->SetMultiplicities(mult);
  fJets->SetInJet(injet);
  fJets->SetEtaIn(etaT);
  fJets->SetPhiIn(phiT);
  fJets->SetPtIn(ptT);
  fJets->SetEtAvg(UA1JETS.etavg);
  delete ncell;
  delete enT;
  delete ptT;
  delete etaT;
  delete phiT;
  delete injet;
  delete idx;
  delete mult;
  delete percentage;
}

void AliUA1JetFinder::FindJets()
{ // Find jets with TPC or EMCAL or TPC + EMCAL information
  // initialize event, look for jets, download jet info

  // 1) transform input to pt,eta,phi plus lego
    AliJetUnitArray *fUnit = fReader->GetUnitArray();
    Int_t nIn = fUnit->GetUnitEntries();
    Int_t fOption = fReader->GetReaderHeader()->GetDetector();
    Int_t fDebug = fReader->GetReaderHeader()->GetDebug();
    
    if(fDebug>1){ 
	printf("Inside FindJets function ! \n");
	printf("fOption : %d \n",fOption);
    }
    
    if(fDebug>10){ 
	cout << "fUnit : " << fUnit << endl;
	printf("nIn = fUnit->GetUnitEntries() : %d \n",nIn);
    }
    
    if(fDebug>10){
	printf("     -----------------------------------------------------------------\n");
	printf("     --- All inputs in fUnitArray ---\n");
	for(Int_t i=0; i<nIn ; i++){
	    if(fUnit[i].GetUnitEnergy()!=0){
		printf("     -----------------------------------------------------------------\n");
		cout << "     |  ID : " << fUnit[i].GetUnitID() << "  |  Eta : " << fUnit[i].GetUnitEta() << "  |  Phi : " << fUnit[i].GetUnitPhi() << "  |  Energy : " << fUnit[i].GetUnitEnergy() << " |" <<endl;
	    }
	}
    }
    
    if (nIn == 0) return;
    
    Int_t nCandidateCut = 0;
    Int_t nCandidateCut2 = 0;
    Int_t nCandidate = 0;
    for (Int_t i = 0; i < nIn; i++){
	if(fUnit[i].GetUnitEnergy()>0. && fUnit[i].GetUnitCutFlag()==1) {
	    // Number of good tracks/digits which have passed the pt cut
	    nCandidateCut += 1;
	}
    if(fUnit[i].GetUnitEnergy()>0.) {
	// Number of good tracks/digits without pt cut
	nCandidate += 1;
    }
    }
    
    if(fDebug>=3){
    cout << "nCandidate : " << nCandidate << endl;
    cout << "nCandidateCut : " << nCandidateCut << endl;
    cout << "nCandidateCut2 : " << nCandidateCut2 << endl;
    }
    
  // local arrays for input No Cuts
  // Both pt < ptMin and pt > ptMin
    Float_t* enT       = new Float_t[nCandidate];
    Float_t* ptT       = new Float_t[nCandidate];
    Float_t* etaT      = new Float_t[nCandidate];
    Float_t* phiT      = new Float_t[nCandidate];
    Float_t* detaT     = new Float_t[nCandidate];
    Float_t* dphiT     = new Float_t[nCandidate];
    Float_t* cFlagT    = new Float_t[nCandidate];
    Float_t* cClusterT = new Float_t[nCandidate];
    Float_t* idT       = new Float_t[nCandidate];
    Int_t loop1 = 0;
    
    fJets->SetNinput(nCandidate);
    
    if(fDebug>3){
	cout << "nCandidate : " << nCandidate << endl;
	//    cout << "nMultCandidate : " << nMultCandidate << endl;
    }
    
    Float_t *etCell = new Float_t[nIn];   //! Cell Energy - Extracted from UnitArray
    Float_t *etaCell = new Float_t[nIn];  //! Cell eta - Extracted from UnitArray
    Float_t *phiCell = new Float_t[nIn];  //! Cell phi - Extracted from UnitArray
    
    // Information extracted from fUnitArray
    for(Int_t i=0; i<nIn; i++) 
    {
	if(fUnit[i].GetUnitCutFlag()==1){ 
	    etCell[i] = fUnit[i].GetUnitEnergy();
	    if (etCell[i] > 0.0) etCell[i] -= fHeader->GetMinCellEt();
	    if (etCell[i] < 0.0) etCell[i] = 0.;
	    etaCell[i] = fUnit[i].GetUnitEta();
	    phiCell[i] = fUnit[i].GetUnitPhi();
	}
	else {
	    etCell[i]  = 0.;
	    etaCell[i] = fUnit[i].GetUnitEta();
	    phiCell[i] = fUnit[i].GetUnitPhi();
	}
	
	if(fUnit[i].GetUnitEnergy()>0.){
	    ptT[loop1]   = fUnit[i].GetUnitEnergy();
	    enT[loop1]   = fUnit[i].GetUnitEnergy();
	    etaT[loop1]  = fUnit[i].GetUnitEta();
	    phiT[loop1]  = fUnit[i].GetUnitPhi();
	    detaT[loop1] = fUnit[i].GetUnitDeta();
	    dphiT[loop1] = fUnit[i].GetUnitDphi();
	    cFlagT[loop1]= fUnit[i].GetUnitCutFlag();
	    idT[loop1]   = fUnit[i].GetUnitID();
	    loop1++;
	}
    }
    
    
    if(fDebug > 40) // For comparison
    {
	for(Int_t j=0; j<nIn; j++) {
	    if(etCell[j]>0){
		cout << "etCell[" << j << "] : " << etCell[j] << endl;
		cout << "etaCell[" << j << "] : " << etaCell[j] << endl;
		cout << "phiCell[" << j << "] : " << phiCell[j] << endl;
	    }
	}
    }
    
    // Run the algo. Parameters from header
    //  Int_t nTot      = (fHeader->GetLegoNbinEta())*(fHeader->GetLegoNbinPhi());
    Int_t nTot      = nIn;
    Float_t minmove = fHeader->GetMinMove();
    Float_t maxmove = fHeader->GetMaxMove();
    Int_t mode      = fHeader->GetMode();
    Float_t precbg  = fHeader->GetPrecBg();
    Int_t ierror;
    
    ua1_jet_finder(nIn, nTot, etCell, etaCell, phiCell, 
		   minmove, maxmove, mode, precbg, ierror);
    
    // sort jets
    Int_t * idx  = new Int_t[UA1JETS.njet];
    TMath::Sort(UA1JETS.njet, UA1JETS.etj, idx);
    
    if(fDebug > 20)
    {
	for(Int_t i = 0; i < UA1JETS.njet; i++) 
	{
	    cout << "Number of jets found, UA1JETS.njet : " << UA1JETS.njet << endl;
	    cout << "UA1JETS.etj : " << UA1JETS.etj << endl;
	    cout << "idx[" << i << "] : " << idx[i] << endl;
	    cout << "UA1JETS.etaj[1][" << idx[i] << "] : " << UA1JETS.etaj[1][idx[i]] << endl;
	    cout << "UA1JETS.phij[1][" << idx[i] << "] : " << UA1JETS.phij[1][idx[i]] << endl;
	    cout << "UA1JETS.etj[" << idx[i] << "] : " << UA1JETS.etj[idx[i]] << endl;
	}
    }
    
    // download jet info.   
    for(Int_t i = 0; i < UA1JETS.njet; i++) {
	// reject events outside acceptable eta range
	if (((UA1JETS.etaj[1][idx[i]])> (fHeader->GetJetEtaMax()))
	    || ((UA1JETS.etaj[1][idx[i]]) < (fHeader->GetJetEtaMin())))
	{
	    continue;
	}
	
	Float_t px, py,pz,en,pT; // convert to 4-vector
	px = UA1JETS.etj[idx[i]] * TMath::Cos(UA1JETS.phij[1][idx[i]]);
	py = UA1JETS.etj[idx[i]] * TMath::Sin(UA1JETS.phij[1][idx[i]]);
	pz = UA1JETS.etj[idx[i]] /
	    TMath::Tan(2.0 * TMath::ATan(TMath::Exp(-UA1JETS.etaj[1][idx[i]])));
	en = TMath::Sqrt(px * px + py * py + pz * pz);
	pT = TMath::Sqrt(px * px + py * py);
	
	fJets->AddJet(px, py, pz, en);
	
    }

    // find multiplicities and relationship jet-particle
    // find also percentage of pt from pythia
    
    Int_t* injet = new Int_t[nCandidate];
    Int_t* sflag = new Int_t[nCandidate];
    for (Int_t i = 0; i < nCandidate; i++) {injet[i]= 0;sflag[i]=0;}
    Int_t* mult  = new Int_t[fJets->GetNJets()];
    Int_t* ncell  = new Int_t[fJets->GetNJets()];
    Float_t* percentage  = new Float_t[fJets->GetNJets()];
    
    // Instead of using etaT below, it would be interesting to use the previous fUnitArray object
    // With the particle ID, it is possible to to have access to its physical properties and one can,
    // for example, set if the corresponding particle is inside or outside the jet with the flag 
    // kOutJet/kInJet, other possibilities...
    
    for(Int_t i = 0; i < (fJets->GetNJets()); i++) {
	Float_t pt_sig = 0.0;
	mult[i] = 0;
	ncell[i] = UA1JETS.ncellj[i];
	for (Int_t j = 0; j < nCandidate; j++) {
	    Float_t deta = etaT[j] - fJets->GetEta(i);
	    Float_t dphi = phiT[j] - fJets->GetPhi(i);
	    if (dphi < -TMath::Pi()) dphi= -dphi - 2.0 * TMath::Pi();
	    if (dphi > TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
	    Float_t dr = TMath::Sqrt(deta * deta + dphi * dphi);
	    if (dr < fHeader->GetRadius() && injet[j] == 0) {
		injet[j] = -(i+1);
		if(cFlagT[j] == 1 &&
		   (etaT[j] < fHeader->GetLegoEtaMax()) &&
		   (etaT[j] > fHeader->GetLegoEtaMin())) {
		    injet[j] = i+1;
		    mult[i]++;
		    pt_sig+=enT[j];
		    sflag[j]=1;
		}
	    }
	    if(fDebug>10){
		cout << "mult[" << i << "] : " << mult[i] << endl;
		cout << "ncell[" << i << "] : " << ncell[i] << endl;
	    }
	}
	percentage[i] = (pt_sig-ncell[i]*UA1JETS.etavg)/
	    ((Double_t) fJets->GetPt(i));    
    }
    
    fJets->SetNCells(ncell);
    fJets->SetPtFromSignal(percentage);
    fJets->SetMultiplicities(mult);
    fJets->SetInJet(injet);
    fJets->SetEtaIn(etaT);
    if(fDebug>10){
	for(Int_t i=0; i<nCandidate ; i++){
	    cout << "phiT[" << i << "] : " << phiT[i] << endl;
	    cout << "etaT[" << i << "] : " << etaT[i] << endl;
	}
    }
    fJets->SetPhiIn(phiT);
    fJets->SetPtIn(enT);
    fJets->SetEtAvg(UA1JETS.etavg);
    delete etCell;
    delete etaCell;
    delete phiCell;
    delete ncell;
    delete cFlagT;
    delete cClusterT;
    delete enT;
    delete ptT;
    delete etaT;
    delete phiT;
    delete detaT;
    delete dphiT;
    delete injet;
    delete idx;
    delete mult;
    delete percentage;
    
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
  dEta = (fHeader->GetLegoEtaMax()-fHeader->GetLegoEtaMin())
    /((Float_t) fHeader->GetLegoNbinEta());
  dPhi = (fHeader->GetLegoPhiMax()-fHeader->GetLegoPhiMin())
    /((Float_t) fHeader->GetLegoNbinPhi());

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
	 fHeader->GetLegoNbinEta(), fHeader->GetLegoEtaMin(), 
	 fHeader->GetLegoEtaMax(),  fHeader->GetLegoNbinPhi(), 
	 fHeader->GetLegoPhiMin(), fHeader->GetLegoPhiMax());
}
