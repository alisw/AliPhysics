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
// UA1 Cone Algorithm Jet finder for charged + neutral jet studies
// manages the search for jets using charged particle momentum and 
// neutral cell energy information
// Based on UA1 V1 (from R. Diaz)
// Author: magali.estienne@subatech.in2p3.fr
//---------------------------------------------------------------------

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRefArray.h>
#include "TFile.h"

#include "AliUA1JetFinderV2.h"
#include "AliUA1JetHeaderV1.h"
#include "AliJetUnitArray.h"
#include "AliJetReaderHeader.h"
#include "AliJetReader.h"
#include "AliJetHeader.h"

class TArrayF;
class TFile;
class AliJetReader;
class AliAODJet;

ClassImp(AliUA1JetFinderV2)


////////////////////////////////////////////////////////////////////////
AliUA1JetFinderV2::AliUA1JetFinderV2() :
  AliJetFinder(),
  fLego(0),  
  fOpt(0)
{
  //
  // Constructor
  //
}

////////////////////////////////////////////////////////////////////////
AliUA1JetFinderV2::~AliUA1JetFinderV2()
{
  //
  // Destructor
  //
}

////////////////////////////////////////////////////////////////////////
void AliUA1JetFinderV2::FindJetsC()
{ 
  // 
  //  Used to find jets using charged particle momentum information
  //
  //  1) Fill cell map array
  //  2) calculate total energy and fluctuation level
  //  3) Run algorithm
  //     3.1) look centroides in cell map
  //     3.2) calculate total energy in cones
  //     3.3) flag as a possible jet
  //     3.4) reorder cones by energy
  //  4) subtract backg in accepted jets
  //  5) fill AliJet list

  //  Transform input to pt,eta,phi plus lego
    
  AliUA1JetHeaderV1* header  = (AliUA1JetHeaderV1*) fHeader;
  TClonesArray*      lvArray = fReader->GetMomentumArray();
  Int_t              nIn     = lvArray->GetEntries();
  fDebug   = fHeader->GetDebug();
  
  if (nIn == 0) return;
  
  // local arrays for input
  Float_t* ptT    = new Float_t[nIn];
  Float_t* etaT   = new Float_t[nIn];
  Float_t* phiT   = new Float_t[nIn];
  Int_t*   cFlagT = new Int_t[nIn]; // Temporarily added
  Int_t*   sFlagT = new Int_t[nIn]; // Temporarily added
  Int_t*   injet  = new Int_t[nIn];

  for (Int_t i = 0; i < nIn; i++) {
    ptT[i]    = 0.;
    etaT[i]   = 0.;
    phiT[i]   = 0.;
    cFlagT[i] = 0; 
    sFlagT[i] = 0; 
    injet[i]  = 0;
  }
  //total energy in array
  Float_t  etbgTotal = 0.0;
  TH1F*    hPtTotal  = new TH1F("hPt","Pt distribution of all particles ",100,0.0,15.0);
  
  // load input vectors and calculate total energy in array
  for (Int_t i = 0; i < nIn; i++){
    TLorentzVector *lv = (TLorentzVector*) lvArray->At(i);
    ptT[i]  = lv->Pt();
    etaT[i] = lv->Eta();
    phiT[i] = ((lv->Phi() < 0) ? (lv->Phi()) + 2 * TMath::Pi() : lv->Phi());
    cFlagT[i] = fReader->GetCutFlag(i); 
    sFlagT[i] = fReader->GetSignalFlag(i); 
    
    if (fReader->GetCutFlag(i) != 1) continue;
    fLego->Fill(etaT[i], phiT[i], ptT[i]);
    hPtTotal->Fill(ptT[i]);
    etbgTotal+= ptT[i];
  }
  
  
  // calculate total energy and fluctuation in map
  Double_t meanpt   = hPtTotal->GetMean();
  Double_t ptRMS    = hPtTotal->GetRMS();
  Double_t npart    = hPtTotal->GetEntries();
  Double_t dEtTotal = (TMath::Sqrt(npart))*TMath::Sqrt(meanpt * meanpt + ptRMS*ptRMS);
  
  // arrays to hold jets
  Float_t etaJet[30];  // eta jet
  Float_t phiJet[30];  // phi jet
  Float_t etJet[30];  // et jet
  Float_t etsigJet[30];  // signal et in jet
  Float_t etallJet[30];  // total et in jet (tmp variable)
  Int_t   ncellsJet[30];
  Int_t   multJet[30];
  //--- Added for jet reordering at the end of the jet finding procedure
  Float_t etaJetOk[30];
  Float_t phiJetOk[30];
  Float_t etJetOk[30];
  Float_t etsigJetOk[30];  // signal et in jet
  Float_t etallJetOk[30];  // total et in jet (tmp variable)
  Int_t   ncellsJetOk[30];
  Int_t   multJetOk[30];
  //--------------------------
  Int_t nJets; // to hold number of jets found by algorithm
  Int_t nj;    // number of jets accepted
  Float_t prec  = header->GetPrecBg();
  Float_t bgprec = 1;
  while(bgprec > prec){
    //reset jet arrays in memory
    memset(etaJet,0,sizeof(Float_t)*30);
    memset(phiJet,0,sizeof(Float_t)*30);
    memset(etJet,0,sizeof(Float_t)*30);
    memset(etallJet,0,sizeof(Float_t)*30);
    memset(etsigJet,0,sizeof(Float_t)*30);
    memset(ncellsJet,0,sizeof(Int_t)*30);
    memset(multJet,0,sizeof(Int_t)*30);
    //--- Added for jet reordering at the end of the jet finding procedure
    memset(etaJetOk,0,sizeof(Float_t)*30);
    memset(phiJetOk,0,sizeof(Float_t)*30);
    memset(etJetOk,0,sizeof(Float_t)*30);
    memset(etallJetOk,0,sizeof(Float_t)*30);
    memset(etsigJetOk,0,sizeof(Float_t)*30);
    memset(ncellsJetOk,0,sizeof(Int_t)*30);
    memset(multJetOk,0,sizeof(Int_t)*30);
    //--------------------------
    nJets = 0;
    nj = 0;
    
    // reset particles-jet array in memory
    memset(injet,-1,sizeof(Int_t)*nIn);
    //run cone algorithm finder
    RunAlgoritmC(etbgTotal,dEtTotal,nJets,etJet,etaJet,phiJet,etallJet,ncellsJet);
    
    //run background subtraction
    if(nJets > header->GetNAcceptJets()) // limited number of accepted jets per event
      nj = header->GetNAcceptJets();
    else
      nj = nJets;
    //subtract background
    Float_t etbgTotalN = 0.0; //new background
    if(header->GetBackgMode() == 1) // standard
      SubtractBackgC(nIn,nj,etbgTotalN,ptT,etaT,phiT,etJet,etaJet,phiJet,etsigJet,multJet,injet);
    if(header->GetBackgMode() == 2) //cone
      SubtractBackgCone(nIn,nj,etbgTotalN,ptT,etaT,phiT,cFlagT,sFlagT,etJet,etaJet,phiJet,etsigJet,multJet,injet);
    if(header->GetBackgMode() == 3) //ratio
      SubtractBackgRatio(nIn,nj,etbgTotalN,ptT,etaT,phiT,cFlagT,sFlagT,etJet,etaJet,phiJet,etsigJet,multJet,injet);
    if(header->GetBackgMode() == 4) //statistic
      SubtractBackgStat(nIn,nj,etbgTotalN,ptT,etaT,phiT,cFlagT,sFlagT,etJet,etaJet,phiJet,etsigJet,multJet,injet);
    //calc precision
    if(TMath::Abs(etbgTotalN) > 0.001)
      bgprec = (etbgTotal - etbgTotalN)/etbgTotalN;
    else
      bgprec = 0;
    etbgTotal = etbgTotalN; // update with new background estimation
  } //end while
  
  // add jets to list
  Int_t* idxjets = new Int_t[nj];
  Int_t nselectj = 0;
  printf("Found %d jets \n", nj);

  // Reorder jets by et in cone
  Int_t * idx  = new Int_t[nJets];
  TMath::Sort(nJets, etJet, idx);
  for(Int_t p = 0; p < nJets; p++){
    etaJetOk[p]    = etaJet[idx[p]];
    phiJetOk[p]    = phiJet[idx[p]];
    etJetOk[p]     = etJet[idx[p]];
    etallJetOk[p]  = etJet[idx[p]];
    etsigJetOk[p]  = etsigJet[idx[p]];
    ncellsJetOk[p] = ncellsJet[idx[p]];
    multJetOk[p]   = multJet[idx[p]];
  }
  
  for(Int_t kj=0; kj<nj; kj++)
    {
      if ((etaJetOk[kj] > (header->GetJetEtaMax())) ||
	  (etaJetOk[kj] < (header->GetJetEtaMin())) ||
	  (etJetOk[kj] < header->GetMinJetEt())) continue; // acceptance eta range and etmin
      Float_t px, py,pz,en; // convert to 4-vector
      px = etJetOk[kj] * TMath::Cos(phiJetOk[kj]);
      py = etJetOk[kj] * TMath::Sin(phiJetOk[kj]);
      pz = etJetOk[kj] / TMath::Tan(2.0 * TMath::ATan(TMath::Exp(-etaJetOk[kj])));
      en = TMath::Sqrt(px * px + py * py + pz * pz);
      
      AliAODJet jet(px, py, pz, en);
      jet.Print("");
      
      AddJet(jet);
      
      idxjets[nselectj] = kj;
      nselectj++;
    }

  //add signal percentage and total signal  in AliJets for analysis tool
  Float_t* percentage  = new Float_t[nselectj];
  Int_t* ncells      = new Int_t[nselectj];
  Int_t* mult        = new Int_t[nselectj];
  for(Int_t i = 0; i< nselectj; i++)
    {
      percentage[i] = etsigJetOk[idxjets[i]]/etJetOk[idxjets[i]];
      ncells[i] = ncellsJetOk[idxjets[i]];
      mult[i] = multJetOk[idxjets[i]];
    }

  //add particle-injet relationship ///
  for(Int_t bj = 0; bj < nIn; bj++)
    {
      if(injet[bj] == -1) continue; //background particle
      Int_t bflag = 0;
      for(Int_t ci = 0; ci< nselectj; ci++){
	if(injet[bj] == idxjets[ci]){
	  injet[bj]= ci;
	  bflag++;
	  break;
	}
      }
      if(bflag == 0) injet[bj] = -1; // set as background particle
    }

 
  //delete
  delete[] ptT;
  delete[] etaT;
  delete[] phiT;
  delete[] cFlagT;
  delete[] sFlagT;
  delete[] injet;
  delete hPtTotal;
  delete[] idxjets;
  delete[] idx;

  delete[] percentage;
  delete[] ncells;
  delete[] mult;
  //--------------------------

}

////////////////////////////////////////////////////////////////////////
void AliUA1JetFinderV2::FindJets()
{
  // 
  //  Used to find jets using charged particle momentum information 
  //  & neutral energy from calo cells
  //
  //  1) Fill cell map array
  //  2) calculate total energy and fluctuation level
  //  3) Run algorithm
  //     3.1) look centroides in cell map
  //     3.2) calculate total energy in cones
  //     3.3) flag as a possible jet
  //     3.4) reorder cones by energy
  //  4) subtract backg in accepted jets
  //  5) fill AliJet list

  // transform input to pt,eta,phi plus lego
    
  AliUA1JetHeaderV1* header   = (AliUA1JetHeaderV1*) fHeader;
  TClonesArray*      fUnit    = fReader->GetUnitArray();
  Int_t              nCand    = fReader->GetNumCandidate();
  Int_t              nCandCut = fReader->GetNumCandidateCut();
  Int_t              nIn      = fUnit->GetEntries();
  Float_t            ptMin   = fReader->GetReaderHeader()->GetPtCut();

  if (nIn == 0) return;

  Int_t nCandidateCut = 0;
  Int_t nCandidate = 0;
  
  nCandidate = nCand;
  nCandidateCut = nCandCut;

  // local arrays for input No Cuts
  // Both pt < ptMin and pt > ptMin
  Float_t*   ptT       = new Float_t[nCandidate];
  Float_t*   en2T      = new Float_t[nCandidate];
  Float_t*   pt2T      = new Float_t[nCandidate];
  Int_t*     detT      = new Int_t[nCandidate]; 
  Float_t*   etaT      = new Float_t[nCandidate];
  Float_t*   phiT      = new Float_t[nCandidate];
  Int_t*     cFlagT    = new Int_t[nCandidate];
  Int_t*     cFlag2T   = new Int_t[nCandidate];
  Int_t*     sFlagT    = new Int_t[nCandidate];
  Float_t*   cClusterT = new Float_t[nCandidate];
  Int_t*     vectT     = new Int_t[nCandidate];
  Int_t      loop1     = 0;
  Int_t*     injet     = new Int_t[nCandidate];
  Int_t*     sflag     = new Int_t[nCandidate];
  
  for(Int_t i = 0; i < nCandidate; i++) {
      ptT[i]       = 0;
      en2T[i]      = 0;
      pt2T[i]      = 0;
      detT[i]      = 0;
      etaT[i]      = 0;
      phiT[i]      = 0;
      cFlagT[i]    = 0;
      cFlag2T[i]   = 0;
      sFlagT[i]    = 0;
      cClusterT[i] = 0;
      vectT[i]     = 0;
      injet[i]     = 0;
      sflag[i]     = 0;
}

  TRefArray* trackRef  = new TRefArray();

  //total energy in array
  Float_t etbgTotal = 0.0;
  TH1F*    hPtTotal  = new TH1F("hPt","Pt distribution of all particles ",100,0.0,15.0);

  // Input cell info
  Float_t *etCell    = new Float_t[nIn];   //! Cell Energy - Extracted from UnitArray
  Float_t *etaCell   = new Float_t[nIn];  //! Cell eta - Extracted from UnitArray
  Float_t *phiCell   = new Float_t[nIn];  //! Cell phi - Extracted from UnitArray
  Int_t   *flagCell  = new Int_t[nIn];   //! Cell phi - Extracted from UnitArray
  Float_t *etCell2   = new Float_t[nIn];  //! Cell Energy - Extracted from UnitArray
  Float_t *etaCell2  = new Float_t[nIn]; //! Cell eta - Extracted from UnitArray
  Float_t *phiCell2  = new Float_t[nIn]; //! Cell phi - Extracted from UnitArray
  Int_t   *flagCell2 = new Int_t[nIn];  //! Cell phi - Extracted from UnitArray
  for(Int_t i = 0; i < nIn; i++) {
    etCell[i]    = 0.;
    etaCell[i]   = 0.;
    phiCell[i]   = 0.;
    flagCell[i]  = 0;
    etCell2[i]   = 0.;
    etaCell2[i]  = 0.;
    phiCell2[i]  = 0.;
    flagCell2[i] = 0;
  }
  // Information extracted from fUnitArray
  // Load input vectors and calculate total energy in array
  for(Int_t i=0; i<nIn; i++) 
    {
      // Recover particle information from UnitArray
      
      AliJetUnitArray *uArray = (AliJetUnitArray*)fUnit->At(i);
      TRefArray* ref = uArray->GetUnitTrackRef();
      Int_t nRef = ref->GetEntries();

      if(uArray->GetUnitEnergy()>0.){

	for(Int_t jpart=0; jpart<nRef;jpart++)
	  trackRef->Add((AliVTrack*)ref->At(jpart));
	ptT[loop1]   = uArray->GetUnitEnergy();
        detT[loop1]  = uArray->GetUnitDetectorFlag();
	etaT[loop1]  = uArray->GetUnitEta();
	phiT[loop1]  = uArray->GetUnitPhi();
	cFlagT[loop1]= uArray->GetUnitCutFlag();   // pt cut tpc
	cFlag2T[loop1]= uArray->GetUnitCutFlag2(); // pt cut emcal
	sFlagT[loop1]= uArray->GetUnitSignalFlag();
	vectT[loop1] = nRef;
	if(cFlagT[loop1] == 1 || cFlag2T[loop1] == 1) {
	  pt2T[loop1] = 0.;
	  en2T[loop1] = 0.;
	  if(detT[loop1]==1){
	    en2T[loop1] = ptT[loop1] - header->GetMinCellEt();
            if(en2T[loop1] < 0) en2T[loop1]=0;
	    hPtTotal->Fill(en2T[loop1]);
	    etbgTotal += en2T[loop1];
	  }
	  if(detT[loop1]==0){ // TPC+ITS
	    Float_t pt = 0.;
	    for(Int_t j=0; j<nRef;j++){
	      Float_t x=0.;  Float_t y=0.;  Float_t z=0.;
	      x = ((AliVTrack*)ref->At(j))->Px();
	      y = ((AliVTrack*)ref->At(j))->Py();
	      z = ((AliVTrack*)ref->At(j))->Pz();
	      pt = TMath::Sqrt(x*x+y*y);
	      if(pt>ptMin) {
		pt2T[loop1] += pt;
		en2T[loop1] += pt;
		hPtTotal->Fill(pt);
		etbgTotal+= pt;
	      }
	    }
	  }
	  if(detT[loop1]==2) { // EMCal
	    Float_t ptCTot = 0.;
	    Float_t pt = 0.;
	    Float_t enC = 0.;
	    for(Int_t j=0; j<nRef;j++){
	      Float_t x=0.;  Float_t y=0.;  Float_t z=0.;
	      x = ((AliVTrack*)ref->At(j))->Px();
	      y = ((AliVTrack*)ref->At(j))->Py();
	      z = ((AliVTrack*)ref->At(j))->Pz();
	      pt = TMath::Sqrt(x*x+y*y);
	      if(pt>ptMin) {
		pt2T[loop1]+=pt;
		en2T[loop1]+=pt;
		hPtTotal->Fill(pt);
		etbgTotal+= pt;
	      }	
	      ptCTot += pt;
	    }
	    enC = ptT[loop1] - ptCTot - header->GetMinCellEt();
            if(enC < 0.) enC=0.;
	    en2T[loop1] += enC;
	    hPtTotal->Fill(enC);
	    etbgTotal+= enC;
	  }
	}	
	loop1++;
      }

      if(uArray->GetUnitCutFlag()==1) {
        if(uArray->GetUnitDetectorFlag()==1){ // EMCal case
          etCell[i] = uArray->GetUnitEnergy() - header->GetMinCellEt();
          if ((uArray->GetUnitEnergy() - header->GetMinCellEt()) < 0.0) etCell[i]=0.;
          etaCell[i] = uArray->GetUnitEta();
          phiCell[i] = uArray->GetUnitPhi();
          flagCell[i] = 0; // default
          etCell2[i] = etCell[i];
          etaCell2[i] = uArray->GetUnitEta();
          phiCell2[i] = uArray->GetUnitPhi();
          flagCell2[i] = 0; // default
        }
        if(uArray->GetUnitDetectorFlag()==0){ // TPC case
          Float_t pt = 0.; Float_t et1 = 0.; Float_t et2 = 0.;
          for(Int_t j=0; j<nRef;j++)
            {
              Float_t x=0.;  Float_t y=0.;  Float_t z=0.;
	      x = ((AliVTrack*)ref->At(j))->Px();
	      y = ((AliVTrack*)ref->At(j))->Py();
	      z = ((AliVTrack*)ref->At(j))->Pz();
              pt = TMath::Sqrt(x*x+y*y);
              if(pt>ptMin) {
                et1 += pt;
                et2 += pt;
              }
            }
          etCell[i] = et1;
          etCell2[i] = et2;
          if(et1 < 0.) etCell[i] = etCell2[i] = 0.;
          etaCell[i] = uArray->GetUnitEta();
          phiCell[i] = uArray->GetUnitPhi();
          flagCell[i] = 0; // default
          etaCell2[i] = uArray->GetUnitEta();
          phiCell2[i] = uArray->GetUnitPhi();
          flagCell2[i] = 0; // default
        }
        if(uArray->GetUnitDetectorFlag()==2){ // TPC + EMCal case
          Float_t ptCTot = 0.;
          Float_t pt = 0.; Float_t et1 = 0.; Float_t et2 = 0.;
          Float_t enC = 0.;
          for(Int_t j=0; j<nRef;j++)
            {
              Float_t x=0.;  Float_t y=0.;  Float_t z=0.;
	      x = ((AliVTrack*)ref->At(j))->Px();
	      y = ((AliVTrack*)ref->At(j))->Py();
	      z = ((AliVTrack*)ref->At(j))->Pz();
              pt = TMath::Sqrt(x*x+y*y);
              if(pt>ptMin) {
                et1 += pt;
                et2 += pt;
              }
              ptCTot += pt;
            }
          enC = uArray->GetUnitEnergy() - ptCTot;
          etCell[i] = et1 + enC - header->GetMinCellEt();
          etCell2[i] = et2 + enC - header->GetMinCellEt();
          if((enC + et1 - header->GetMinCellEt()) < 0.) etCell[i] = etCell2[i] = 0.;
          etaCell[i] = uArray->GetUnitEta();
          phiCell[i] = uArray->GetUnitPhi();
          flagCell[i] = 0; // default
          etaCell2[i] = uArray->GetUnitEta();
          phiCell2[i] = uArray->GetUnitPhi();
          flagCell2[i] = 0; // default
        }
      }
      else {
        etCell[i]  = 0.;
        etaCell[i] = uArray->GetUnitEta();
        phiCell[i] = uArray->GetUnitPhi();
        flagCell[i] = 0;
        etCell2[i]  = 0.;
        etaCell2[i] = uArray->GetUnitEta();
        phiCell2[i] = uArray->GetUnitPhi();
        flagCell2[i] = 0;
      }
    } // end loop on nCandidate


  // calculate total energy and fluctuation in map
  Double_t meanpt = hPtTotal->GetMean();
  Double_t ptRMS = hPtTotal->GetRMS();
  Double_t npart = hPtTotal->GetEntries();
  Double_t dEtTotal = (TMath::Sqrt(npart))*TMath::Sqrt(meanpt * meanpt + ptRMS*ptRMS);

  // arrays to hold jets
  Float_t etaJet[30];
  Float_t phiJet[30];
  Float_t etJet[30];
  Float_t etsigJet[30]; //signal et in jet
  Float_t etallJet[30];  // total et in jet (tmp variable)
  Int_t   ncellsJet[30];
  Int_t   multJet[30];
  //--- Added by me for jet reordering at the end of the jet finding procedure
  Float_t etaJetOk[30];
  Float_t phiJetOk[30];
  Float_t etJetOk[30];
  Float_t etsigJetOk[30]; //signal et in jet
  Float_t etallJetOk[30];  // total et in jet (tmp variable)
  Int_t   ncellsJetOk[30];
  Int_t   multJetOk[30];
  //--------------------------
  Int_t    nJets; // to hold number of jets found by algorithm
  Int_t    nj;    // number of jets accepted
  Float_t  prec  = header->GetPrecBg();
  Float_t  bgprec = 1;

  while(bgprec > prec){

    //reset jet arrays in memory
    memset(etaJet,0,sizeof(Float_t)*30);
    memset(phiJet,0,sizeof(Float_t)*30);
    memset(etJet,0,sizeof(Float_t)*30);
    memset(etallJet,0,sizeof(Float_t)*30);
    memset(etsigJet,0,sizeof(Float_t)*30);
    memset(ncellsJet,0,sizeof(Int_t)*30);
    memset(multJet,0,sizeof(Int_t)*30);
    //--- Added by me for jet reordering at the end of the jet finding procedure
    memset(etaJetOk,0,sizeof(Float_t)*30);
    memset(phiJetOk,0,sizeof(Float_t)*30);
    memset(etJetOk,0,sizeof(Float_t)*30);
    memset(etallJetOk,0,sizeof(Float_t)*30);
    memset(etsigJetOk,0,sizeof(Float_t)*30);
    memset(ncellsJetOk,0,sizeof(Int_t)*30);
    memset(multJetOk,0,sizeof(Int_t)*30);

    nJets = 0;  
    nj = 0;

    // reset particles-jet array in memory
    memset(injet,-1,sizeof(Int_t)*nCandidate);
    //run cone algorithm finder
    RunAlgoritm(nIn,etCell,etaCell,phiCell,flagCell,etCell2,etaCell2,phiCell2,
		flagCell2,etbgTotal,dEtTotal,nJets,etJet,etaJet,phiJet,
		etallJet,ncellsJet);

    //run background subtraction
    if(nJets > header->GetNAcceptJets()) // limited number of accepted jets per event
      nj = header->GetNAcceptJets();
    else
      nj = nJets;

    //subtract background
    Float_t etbgTotalN = 0.0; //new background
    if(header->GetBackgMode() == 1) // standard
      SubtractBackg(nCandidate,nj,etbgTotalN,en2T,vectT,etaT,phiT,cFlagT,cFlag2T,sFlagT,etJet,etaJet,phiJet,etsigJet,multJet,injet);
    // To be modified ------------------------
    if(header->GetBackgMode() == 2) //cone
      SubtractBackgCone(nCandidate,nj,etbgTotalN,ptT,etaT,phiT,cFlagT,sFlagT,etJet,etaJet,phiJet,etsigJet,multJet,injet);
    if(header->GetBackgMode() == 3) //ratio
      SubtractBackgRatio(nCandidate,nj,etbgTotalN,ptT,etaT,phiT,cFlagT,sFlagT,etJet,etaJet,phiJet,etsigJet,multJet,injet);
    if(header->GetBackgMode() == 4) //statistic
      SubtractBackgStat(nCandidate,nj,etbgTotalN,ptT,etaT,phiT,cFlagT,sFlagT,etJet,etaJet,phiJet,etsigJet,multJet,injet);
    //----------------------------------------
    //calc precision
    if(etbgTotalN != 0.0)
      bgprec = (etbgTotal - etbgTotalN)/etbgTotalN;
    else
      bgprec = 0;
    etbgTotal = etbgTotalN; // update with new background estimation
  } //end while
  
  // add jets to list
  Int_t* idxjets = new Int_t[nj];
  Int_t nselectj = 0;
  printf("Found %d jets \n", nj);
  
  // Reorder jets by et in cone
  // Sort jets by energy
  Int_t * idx  = new Int_t[nJets];
  TMath::Sort(nJets, etJet, idx);
  for(Int_t p = 0; p < nJets; p++)
    {
      etaJetOk[p]    = etaJet[idx[p]];
      phiJetOk[p]    = phiJet[idx[p]];
      etJetOk[p]     = etJet[idx[p]];
      etallJetOk[p]  = etJet[idx[p]];
      etsigJetOk[p]  = etsigJet[idx[p]];
      ncellsJetOk[p] = ncellsJet[idx[p]];
      multJetOk[p]   = multJet[idx[p]];
    }

  TRefArray *refs = 0;
  Bool_t fromAod = !strcmp(fReader->ClassName(),"AliJetAODReader");
  if (fromAod) refs = fReader->GetReferences();
  Int_t nTracks = 0;
  if (fromAod) nTracks = ((TRefArray*)refs)->GetEntries();
  Int_t* trackinjet     = new Int_t[nTracks];
  for(Int_t it=0; it<nTracks; it++) trackinjet[it]=-1;

  for(Int_t kj=0; kj<nj; kj++)
    {
      if ((etaJetOk[kj] > (header->GetJetEtaMax())) ||
	  (etaJetOk[kj] < (header->GetJetEtaMin())) ||
	  (etJetOk[kj] < header->GetMinJetEt())) continue; // acceptance eta range and etmin
      Float_t px, py,pz,en; // convert to 4-vector
      px = etJetOk[kj] * TMath::Cos(phiJetOk[kj]);
      py = etJetOk[kj] * TMath::Sin(phiJetOk[kj]);
      pz = etJetOk[kj] / TMath::Tan(2.0 * TMath::ATan(TMath::Exp(-etaJetOk[kj])));
      en = TMath::Sqrt(px * px + py * py + pz * pz);

      AliAODJet jet(px, py, pz, en);
      jet.Print("");
      
      if (fromAod){
	for(Int_t jpart = 0; jpart < nTracks; jpart++) { // loop for all particles in array
	  Float_t deta = ((AliAODTrack*)refs->At(jpart))->Eta() - etaJetOk[kj];
	  Float_t dphi = ((AliAODTrack*)refs->At(jpart))->Phi() - phiJetOk[kj];
	  if (dphi < -TMath::Pi()) dphi= -dphi - 2.0 * TMath::Pi();
	  if (dphi > TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
	  
	  Float_t dr = TMath::Sqrt(deta * deta + dphi * dphi);
	  if(dr <= header->GetRadius() && fReader->GetCutFlag(jpart) == 1) {
	      // particles inside this cone
	      if(trackinjet[jpart]==-1) {
		  trackinjet[jpart] = kj;
	      } else if(fDebug>10) {
		  printf("The track already belongs to jet %d \n",trackinjet[jpart]);
	      }
	  }
	  if(trackinjet[jpart]==kj)
	      jet.AddTrack(refs->At(jpart));  // check if the particle belongs to the jet and add the ref
	}
      }

      AddJet(jet);
      
      idxjets[nselectj] = kj;
      nselectj++;
    }

  //add signal percentage and total signal  in AliJets for analysis tool
  Float_t* percentage  = new Float_t[nselectj];
  Int_t* ncells      = new Int_t[nselectj];
  Int_t* mult        = new Int_t[nselectj];
  for(Int_t i = 0; i< nselectj; i++)
    {
      percentage[i] = etsigJetOk[idxjets[i]]/etJetOk[idxjets[i]];
      ncells[i]     = ncellsJetOk[idxjets[i]];
      mult[i]       = multJetOk[idxjets[i]];
    }

  //add particle-injet relationship ///
  for(Int_t bj = 0; bj < nCandidate; bj++)
    {
      if(injet[bj] == -1) continue; //background particle
      Int_t bflag = 0;
      for(Int_t ci = 0; ci< nselectj; ci++){
	if(injet[bj] == idxjets[ci]){
	  injet[bj]= ci;
	  bflag++;
	  break;
	}
      }
    if(bflag == 0) injet[bj] = -1; // set as background particle
  }


  //delete
  delete [] ptT;
  delete [] en2T;
  delete [] pt2T;
  delete [] etaT;
  delete [] phiT;
  delete [] detT;
  delete [] cFlagT;
  delete [] cFlag2T;
  delete [] sFlagT;
  delete [] cClusterT;
  delete [] vectT;
  delete [] injet;
  delete [] sflag;
  trackRef->Delete();
  delete trackRef;

  delete hPtTotal;
  delete [] etCell;
  delete [] etaCell;
  delete [] phiCell;
  delete [] flagCell;
  delete [] etCell2;
  delete [] etaCell2;
  delete [] phiCell2;
  delete [] flagCell2;
  //--------------------------
  delete [] idxjets;
  delete [] idx;
  delete [] trackinjet;

  delete [] percentage;
  delete [] ncells;
  delete [] mult;

}

////////////////////////////////////////////////////////////////////////
void AliUA1JetFinderV2::RunAlgoritm(Int_t nIn, Float_t* etCell, Float_t* const etaCell, Float_t* phiCell, 
		   Int_t* const flagCell, const Float_t* etCell2, const Float_t* etaCell2, const Float_t* phiCell2, 
		   const Int_t* flagCell2, Float_t etbgTotal, Double_t dEtTotal, 
		   Int_t& nJets, Float_t* const etJet, Float_t* const etaJet, Float_t* const phiJet,
		   Float_t* const etallJet, Int_t* const ncellsJet)
{
  //
  // Main method for jet finding
  // UA1 base cone finder
  //

  Int_t nCell  = nIn;

  // Dump lego
  // Check enough space! *to be done*
  AliUA1JetHeaderV1* header = (AliUA1JetHeaderV1*) fHeader;
  for(Int_t i=0; i<nCell; i++){
    etCell[i]   = etCell2[i];
    etaCell[i]  = etaCell2[i];
    phiCell[i]  = phiCell2[i];
    flagCell[i] = flagCell2[i];
  }

  // Parameters from header
  Float_t minmove = header->GetMinMove();
  Float_t maxmove = header->GetMaxMove();
  Float_t rc      = header->GetRadius();
  Float_t etseed  = header->GetEtSeed();

  // Tmp array of jets form algoritm
  Float_t etaAlgoJet[30]    = {0.};
  Float_t phiAlgoJet[30]    = {0.};
  Float_t etAlgoJet[30]     = {0.};
  Int_t   ncellsAlgoJet[30] = {0};

  // Run algorithm//

  // Sort cells by et
  Int_t * index  = new Int_t[nCell];
  TMath::Sort(nCell, etCell, index);

  // Variable used in centroide loop
  Float_t eta   = 0.0;
  Float_t phi   = 0.0;
  Float_t eta0  = 0.0;
  Float_t phi0  = 0.0;
  Float_t etab  = 0.0;
  Float_t phib  = 0.0;
  Float_t etas  = 0.0;
  Float_t phis  = 0.0;
  Float_t ets   = 0.0;
  Float_t deta  = 0.0;
  Float_t dphi  = 0.0;
  Float_t dr    = 0.0;
  Float_t etsb  = 0.0;
  Float_t etasb = 0.0;
  Float_t phisb = 0.0;
  Float_t dphib = 0.0;

  for(Int_t icell = 0; icell < nCell; icell++)
    {
      Int_t jcell = index[icell];
      if(etCell[jcell] <= etseed) continue; // if cell energy is low et seed
      if(flagCell[jcell] != 0) continue; // if cell was used before

      eta  = etaCell[jcell];
      phi  = phiCell[jcell];
      eta0 = eta;
      phi0 = phi;
      etab = eta;
      phib = phi;
      ets  = etCell[jcell];
      etas = 0.0;
      phis = 0.0;
      etsb = ets;
      etasb = 0.0;
      phisb = 0.0;
      for(Int_t kcell =0; kcell < nCell; kcell++)
	{
	  Int_t lcell = index[kcell];
	  if(lcell == jcell) continue; // cell itself
	  if(flagCell[lcell] != 0) continue; // cell used before
	  if(etCell[lcell] > etCell[jcell]) continue;  // can this happen
	  //calculate dr
	  deta = etaCell[lcell] - eta;
	  dphi = TMath::Abs(phiCell[lcell] - phi);
	  if (dphi > TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
	  dr = TMath::Sqrt(deta * deta + dphi * dphi);
	  if(dr <= rc){
	    // calculate offset from initiate cell
	    deta = etaCell[lcell] - eta0;
	    dphi = phiCell[lcell] - phi0;
	    if (dphi < -TMath::Pi()) dphi= dphi + 2.0 * TMath::Pi();
	    if (dphi > TMath::Pi()) dphi = dphi - 2.0 * TMath::Pi();
	    etas = etas + etCell[lcell]*deta;
	    phis = phis + etCell[lcell]*dphi;
	    ets = ets + etCell[lcell];
	    //new weighted eta and phi including this cell
	    eta = eta0 + etas/ets;
	    phi = phi0 + phis/ets;
	    // if cone does not move much, just go to next step
	    dphib = TMath::Abs(phi - phib);
	    if (dphib > TMath::Pi()) dphib = 2. * TMath::Pi() - dphib;
	    dr = TMath::Sqrt((eta-etab)*(eta-etab) + dphib * dphib);
	    if(dr <= minmove) break;
	    // cone should not move more than max_mov
	    dr = TMath::Sqrt((etas/ets)*(etas/ets) + (phis/ets)*(phis/ets));
	    if(dr > maxmove){
	      eta = etab;
	      phi = phib;
	      ets = etsb;
	      etas = etasb;
	      phis = phisb;
	    } else { // store this loop information
	      etab  = eta;
	      phib  = phi;
	      etsb  = ets;
	      etasb = etas;
	      phisb = phis;
	    }
	  } // inside cone
        }//end of cells loop looking centroide

        //avoid cones overloap (to be implemented in the future)

        //flag cells in Rc, estimate total energy in cone
      Float_t etCone = 0.0;
      Int_t   nCellIn = 0;
      Int_t   nCellOut = 0;
      rc = header->GetRadius();

      for(Int_t ncell =0; ncell < nCell; ncell++)
	{
	  if(flagCell[ncell] != 0) continue; // cell used before
	  //calculate dr
	  deta = etaCell[ncell] - eta;
	  //	    if(deta <= rc){ // Added to improve velocity -> to be tested
	  dphi = phiCell[ncell] - phi;
	  if (dphi < -TMath::Pi()) dphi= -dphi - 2.0 * TMath::Pi();
	  if (dphi > TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
	  //	      if(dphi <= rc){ // Added to improve velocity -> to be tested
	  dr = TMath::Sqrt(deta * deta + dphi * dphi);
	  if(dr <= rc){  // cell in cone
	    flagCell[ncell] = -1;
	    etCone+=etCell[ncell];
	    nCellIn++;
	  }
	  else nCellOut++;
	  //	} // end deta <= rc
	  //	    } // end dphi <= rc
        }

      // select jets with et > background
      // estimate max fluctuation of background in cone
      Double_t ncellin = (Double_t)nCellIn;
      Double_t ntcell  = (Double_t)nCell;
      Double_t etbmax = (etbgTotal + dEtTotal )*(ncellin/(ntcell));
      // min cone et
      Double_t etcmin = etCone ;  // could be used etCone - etmin !!
      //decisions !! etbmax < etcmin

      for(Int_t mcell =0; mcell < nCell; mcell++)
	{
	  if(flagCell[mcell] == -1){
	    if(etbmax < etcmin)
	      flagCell[mcell] = 1; //flag cell as used
	    else
	      flagCell[mcell] = 0; // leave it free
	  }
        }
      //store tmp jet info !!!
      if(etbmax < etcmin) 
	{
	  etaAlgoJet[nJets]    = eta;
	  phiAlgoJet[nJets]    = phi;
	  etAlgoJet[nJets]     = etCone;
	  ncellsAlgoJet[nJets] = nCellIn;
	  nJets++;
	}

    } // end of cells loop

  for(Int_t p = 0; p < nJets; p++)
    {
      etaJet[p]    = etaAlgoJet[p];
      phiJet[p]    = phiAlgoJet[p];
      etJet[p]     = etAlgoJet[p];
      etallJet[p]  = etAlgoJet[p];
      ncellsJet[p] = ncellsAlgoJet[p];
    }

  //delete
  delete[] index;

}

////////////////////////////////////////////////////////////////////////
void AliUA1JetFinderV2::RunAlgoritmC(Float_t etbgTotal, Double_t dEtTotal, Int_t& nJets,
		   Float_t* const etJet,Float_t* const etaJet, Float_t* const phiJet,
		   Float_t* const etallJet, Int_t* const ncellsJet)
{
  // Dump lego
  // Check enough space! *to be done*
  AliUA1JetHeaderV1* header = (AliUA1JetHeaderV1*) fHeader;
  Float_t etCell[60000];   //! Cell Energy
  Float_t etaCell[60000];  //! Cell eta
  Float_t phiCell[60000];  //! Cell phi
  Int_t   flagCell[60000]; //! Cell flag
  
  Int_t nCell = 0;
  TAxis* xaxis = fLego->GetXaxis();
  TAxis* yaxis = fLego->GetYaxis();
  Float_t e = 0.0;
  for (Int_t i = 1; i <= header->GetLegoNbinEta(); i++) 
    {
      for (Int_t j = 1; j <= header->GetLegoNbinPhi(); j++)
	{
	  e = fLego->GetBinContent(i,j);
	  if (e < 0.0) continue; // don't include this cells
	  Float_t eta  = xaxis->GetBinCenter(i);
	  Float_t phi  = yaxis->GetBinCenter(j);
	  etCell[nCell]  = e;
	  etaCell[nCell] = eta;
	  phiCell[nCell] = phi;
	  flagCell[nCell] = 0; //default
	  nCell++;
	}
    }

  // Parameters from header
  Float_t minmove = header->GetMinMove();
  Float_t maxmove = header->GetMaxMove();
  Float_t rc      = header->GetRadius();
  Float_t etseed  = header->GetEtSeed();

  // Tmp array of jets form algoritm
  Float_t etaAlgoJet[30]    = {0.};
  Float_t phiAlgoJet[30]    = {0.};
  Float_t etAlgoJet[30]     = {0.};
  Int_t   ncellsAlgoJet[30] = {0};

  // Run algorithm//

  // Sort cells by et
  Int_t * index  = new Int_t[nCell];
  TMath::Sort(nCell, etCell, index);
  // variable used in centroide loop
  Float_t eta   = 0.0;
  Float_t phi   = 0.0;
  Float_t eta0  = 0.0;
  Float_t phi0  = 0.0;
  Float_t etab  = 0.0;
  Float_t phib  = 0.0;
  Float_t etas  = 0.0;
  Float_t phis  = 0.0;
  Float_t ets   = 0.0;
  Float_t deta  = 0.0;
  Float_t dphi  = 0.0;
  Float_t dr    = 0.0;
  Float_t etsb  = 0.0;
  Float_t etasb = 0.0;
  Float_t phisb = 0.0;
  Float_t dphib = 0.0;

  for(Int_t icell = 0; icell < nCell; icell++)
    {
      Int_t jcell = index[icell];
      if(etCell[jcell] <= etseed) continue; // if cell energy is low et seed
      if(flagCell[jcell] != 0) continue; // if cell was used before
      
      eta  = etaCell[jcell];
      phi  = phiCell[jcell];
      eta0 = eta;
      phi0 = phi;
      etab = eta;
      phib = phi;
      ets  = etCell[jcell];
      etas = 0.0;
      phis = 0.0;
      etsb = ets;
      etasb = 0.0;
      phisb = 0.0;
      for(Int_t kcell =0; kcell < nCell; kcell++)
	{
	  Int_t lcell = index[kcell];
	  if(lcell == jcell) continue; // cell itself
	  if(flagCell[lcell] != 0) continue; // cell used before
	  if(etCell[lcell] > etCell[jcell]) continue; // can this happen
	  //calculate dr
	  deta = etaCell[lcell] - eta;
	  dphi = TMath::Abs(phiCell[lcell] - phi);
	  if (dphi > TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
	  dr = TMath::Sqrt(deta * deta + dphi * dphi);
	  if(dr <= rc)
	    {
	      // calculate offset from initiate cell
	      deta = etaCell[lcell] - eta0;
	      dphi = phiCell[lcell] - phi0;
	      if (dphi < -TMath::Pi()) dphi= dphi + 2.0 * TMath::Pi();
	      if (dphi > TMath::Pi()) dphi = dphi - 2.0 * TMath::Pi();
	      etas = etas + etCell[lcell]*deta;
	      phis = phis + etCell[lcell]*dphi;
	      ets = ets + etCell[lcell];
	      //new weighted eta and phi including this cell
	      eta = eta0 + etas/ets;
	      phi = phi0 + phis/ets;
	      // if cone does not move much, just go to next step
	      dphib = TMath::Abs(phi - phib);
	      if (dphib > TMath::Pi()) dphib = 2. * TMath::Pi() - dphib;
	      dr = TMath::Sqrt((eta-etab)*(eta-etab) + dphib * dphib);
	      if(dr <= minmove) break;
	      // cone should not move more than max_mov
	      dr = TMath::Sqrt((etas/ets)*(etas/ets) + (phis/ets)*(phis/ets));
	      if(dr > maxmove){
		eta = etab;
		phi = phib;
		ets = etsb;
		etas = etasb;
		phis = phisb;
	      } else { // store this loop information
		etab=eta;
		phib=phi;
		etsb = ets;
		etasb = etas;
		phisb = phis;
		}
	    } // inside cone
	}//end of cells loop looking centroide
      
      // Avoid cones overloap (to be implemented in the future)
      
      // Flag cells in Rc, estimate total energy in cone
      Float_t etCone   = 0.0;
      Int_t   nCellIn  = 0;
      Int_t   nCellOut = 0;
      rc = header->GetRadius();
      for(Int_t ncell =0; ncell < nCell; ncell++)
	{
	  if(flagCell[ncell] != 0) continue; // cell used before
	  //calculate dr
	  deta = etaCell[ncell] - eta;
	  dphi = phiCell[ncell] - phi;
	  if (dphi < -TMath::Pi()) dphi= -dphi - 2.0 * TMath::Pi();
	  if (dphi > TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
	  dr = TMath::Sqrt(deta * deta + dphi * dphi);
	  if(dr <= rc){  // cell in cone
	    flagCell[ncell] = -1;
	    etCone+=etCell[ncell];
	    nCellIn++;
	  }
	  else nCellOut++;
	}
      
      // Select jets with et > background
      // estimate max fluctuation of background in cone
      Double_t ncellin = (Double_t)nCellIn;
      Double_t ntcell  = (Double_t)nCell;
      Double_t etbmax = (etbgTotal + dEtTotal )*(ncellin/ntcell);
      // min cone et
      Double_t etcmin = etCone ;  // could be used etCone - etmin !!
      //decisions !! etbmax < etcmin
      
      for(Int_t mcell =0; mcell < nCell; mcell++){
	if(flagCell[mcell] == -1){
	  if(etbmax < etcmin)
	    flagCell[mcell] = 1; //flag cell as used
	  else
	    flagCell[mcell] = 0; // leave it free
	}
      }
      //store tmp jet info !!!
      
      if(etbmax < etcmin) {
	etaAlgoJet[nJets]    = eta;
	phiAlgoJet[nJets]    = phi;
	etAlgoJet[nJets]     = etCone;
	ncellsAlgoJet[nJets] = nCellIn;
	nJets++;
      }
      
    } // end of cells loop

  //reorder jets by et in cone
  //sort jets by energy
  Int_t * idx  = new Int_t[nJets];
  TMath::Sort(nJets, etAlgoJet, idx);
  for(Int_t p = 0; p < nJets; p++)
    {
      etaJet[p] = etaAlgoJet[idx[p]];
      phiJet[p] = phiAlgoJet[idx[p]];
      etJet[p] = etAlgoJet[idx[p]];
      etallJet[p] = etAlgoJet[idx[p]];
      ncellsJet[p] = ncellsAlgoJet[idx[p]];
    }
  
  //delete
  delete [] index;
  delete [] idx;

}

////////////////////////////////////////////////////////////////////////
void AliUA1JetFinderV2::SubtractBackg(const Int_t& nIn, const Int_t&nJ, Float_t& etbgTotalN, const Float_t* ptT, const Int_t* vectT, 
				      const Float_t* etaT, const Float_t* phiT, const Int_t* cFlagT, const Int_t* cFlag2T, 
				      const Int_t* sFlagT, Float_t* const etJet, const Float_t* etaJet, const Float_t* phiJet, 
				      Float_t* const etsigJet, Int_t* const multJet, Int_t* const injet)
{
  //
  // Background subtraction using cone method but without correction in dE/deta distribution
  // Cases to take into account the EMCal geometry are included
  //
  
  //calculate energy inside and outside cones
  AliUA1JetHeaderV1* header = (AliUA1JetHeaderV1*) fHeader;
  fOpt = fReader->GetReaderHeader()->GetDetector();
  Float_t rc= header->GetRadius();
  Float_t etIn[30] = {0.};
  Float_t etOut = 0;
    
  for(Int_t jpart = 0; jpart < nIn; jpart++){ // loop for all particles in array

    for(Int_t ijet=0; ijet<nJ; ijet++){

      Float_t deta = etaT[jpart] - etaJet[ijet];
      Float_t dphi = phiT[jpart] - phiJet[ijet];
      if (dphi < -TMath::Pi()) dphi= -dphi - 2.0 * TMath::Pi();
      if (dphi > TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
      
      Float_t dr = TMath::Sqrt(deta * deta + dphi * dphi);
      if(dr <= rc){ // particles inside this cone
	multJet[ijet]+=vectT[jpart];
	injet[jpart] = ijet;
	
	if(cFlagT[jpart] == 1 || cFlag2T[jpart] == 1){ // pt cut
	  etIn[ijet] += ptT[jpart];
	  if(sFlagT[jpart] == 1) etsigJet[ijet]+= ptT[jpart];
	}
	break;
      }
    }// end jets loop

    if(injet[jpart] == -1 && (cFlagT[jpart] == 1 || cFlag2T[jpart] == 1)){
      etOut += ptT[jpart]; // particle outside cones and pt cut
    }
  } //end particle loop

  //estimate jets and background areas
  // TPC case
  if(fOpt == 0 || fOpt == 1){
    Float_t areaJet[30];
    Float_t areaOut = 4*(header->GetLegoEtaMax())*TMath::Pi();

    for(Int_t k=0; k<nJ; k++){
      Float_t detamax = etaJet[k] + rc;
      Float_t detamin = etaJet[k] - rc;
      Float_t accmax = 0.0; Float_t accmin = 0.0;
      if(detamax > header->GetLegoEtaMax()){ // sector outside etamax
	Float_t h = header->GetLegoEtaMax() - etaJet[k];
	accmax = rc*rc*TMath::ACos(h/rc) - h*TMath::Sqrt(rc*rc - h*h);
      }
      if(detamin < header->GetLegoEtaMin()){ // sector outside etamin
	Float_t h = header->GetLegoEtaMax() + etaJet[k];
	accmin = rc*rc*TMath::ACos(h/rc) - h*TMath::Sqrt(rc*rc - h*h);
      }
      areaJet[k] = rc*rc*TMath::Pi() - accmax - accmin;
      areaOut = areaOut - areaJet[k];
    }
    //subtract background using area method
    for(Int_t ljet=0; ljet<nJ; ljet++){
      Float_t areaRatio = areaJet[ljet]/areaOut;
      etJet[ljet] = etIn[ljet]-etOut*areaRatio; // subtraction
    }

    // estimate new total background
    Float_t areaT = 4*(header->GetLegoEtaMax())*TMath::Pi();
    etbgTotalN = etOut*areaT/areaOut;
  }
  else { // If EMCal included
    Float_t areaJet[30];
    Float_t areaOut = 2*(header->GetLegoEtaMax())*(header->GetLegoPhiMax() - header->GetLegoPhiMin());
    for(Int_t k=0; k<nJ; k++){
      Float_t detamax = etaJet[k] + rc;
      Float_t detamin = etaJet[k] - rc;
      Float_t dphimax = phiJet[k] + rc;
      Float_t dphimin = phiJet[k] - rc;
      Float_t eMax = header->GetLegoEtaMax();
      Float_t eMin = header->GetLegoEtaMin();
      Float_t pMax = header->GetLegoPhiMax();
      Float_t pMin = header->GetLegoPhiMin();
      Float_t accetamax = 0.0; Float_t accetamin = 0.0;
      Float_t accphimax = 0.0; Float_t accphimin = 0.0;
      if((detamax > eMax && dphimax >= (pMin+2*rc) && dphimax <= pMax )||
	 (detamax > eMax && dphimin <= (pMax-2*rc) && dphimin >= pMin )){
	Float_t h = eMax - etaJet[k];
	accetamax = rc*rc*TMath::ACos(h/rc) - h*TMath::Sqrt(rc*rc - h*h);
      }
      if((detamin < eMin && dphimax >= (pMin+2*rc) && dphimax <= pMax )||
	 (detamin < eMin && dphimin <= (pMax-2*rc) && dphimin >= pMin )){
	Float_t h = eMax + etaJet[k];
	accetamin = rc*rc*TMath::ACos(h/rc) - h*TMath::Sqrt(rc*rc - h*h);
      }
      if((dphimax > pMax && detamax >= (eMin+2*rc) && detamax <= eMax )||
 	 (dphimax > pMax && detamin <= (eMax-2*rc) && detamin >= eMin )){
	Float_t h = pMax - phiJet[k];
	accphimax = rc*rc*TMath::ACos(h/rc) - h*TMath::Sqrt(rc*rc - h*h);
      }
      if((dphimin < eMin && detamax >= (eMin+2*rc) && detamax <= eMax )||
 	 (dphimin < eMin && detamin <= (eMax-2*rc) && detamin >= eMin )){
	Float_t h = phiJet[k] - pMin;
	accphimin = rc*rc*TMath::ACos(h/rc) - h*TMath::Sqrt(rc*rc - h*h);
      }
      
      if(detamax > eMax && dphimax > pMax ){
	Float_t he = eMax - etaJet[k];
	Float_t hp = pMax - phiJet[k];
	Float_t rlim = TMath::Sqrt(pow(he,2)+pow(hp,2));
	Float_t alphae = TMath::ACos(he/rc);
	Float_t alphap = TMath::ACos(hp/rc);
	Float_t alphad = (alphae+alphap)/2-TMath::Pi()/4;
	if(rlim <= rc){
	  accetamax = rc*rc*alphae - he*TMath::Sqrt(rc*rc - he*he);
	  accphimax = rc*rc*alphap - hp*TMath::Sqrt(rc*rc - hp*hp);
	}
	if(rlim > rc){
	  accetamax = rc*rc*alphae - he*TMath::Sqrt(rc*rc - he*he);
	  accphimax = rc*rc*alphap - hp*TMath::Sqrt(rc*rc - hp*hp)-
	    ((TMath::Sqrt(pow(rc,2)-pow(he,2))-hp)*(TMath::Sqrt(pow(rc,2)-pow(hp,2))-he))/2+
	    rc*rc*alphad - rc*rc*TMath::Sin(alphad)*TMath::Cos(alphad);
	}
      }

      if(detamax > eMax && dphimin < pMin ){
	Float_t he = eMax - etaJet[k];
	Float_t hp = phiJet[k] - pMin;
	Float_t rlim = TMath::Sqrt(pow(he,2)+pow(hp,2));
	Float_t alphae = TMath::ACos(he/rc);
	Float_t alphap = TMath::ACos(hp/rc);
	Float_t alphad = (alphae+alphap)/2-TMath::Pi()/4;
	if(rlim <= rc){
	  accetamax = rc*rc*alphae - he*TMath::Sqrt(rc*rc - he*he);
	  accphimin = rc*rc*alphap - hp*TMath::Sqrt(rc*rc - hp*hp);
	}
	if(rlim > rc){
	  accetamax = rc*rc*alphae - he*TMath::Sqrt(rc*rc - he*he);
	  accphimin = rc*rc*alphap - hp*TMath::Sqrt(rc*rc - hp*hp)-
	    ((TMath::Sqrt(pow(rc,2)-pow(he,2))-hp)*(TMath::Sqrt(pow(rc,2)-pow(hp,2))-he))/2+
	    rc*rc*alphad - rc*rc*TMath::Sin(alphad)*TMath::Cos(alphad);
	}
      }

      if(detamin < eMin && dphimax > pMax ){
	Float_t he = eMax + etaJet[k];
	Float_t hp = pMax - phiJet[k];
	Float_t rlim = TMath::Sqrt(pow(he,2)+pow(hp,2));
	Float_t alphae = TMath::ACos(he/rc);
	Float_t alphap = TMath::ACos(hp/rc);
	Float_t alphad = (alphae+alphap)/2-TMath::Pi()/4;
	if(rlim <= rc){
	  accetamin = rc*rc*alphae - he*TMath::Sqrt(rc*rc - he*he);
	  accphimax = rc*rc*alphap - hp*TMath::Sqrt(rc*rc - hp*hp);
	}
	if(rlim > rc){
	  accetamin = rc*rc*alphae - he*TMath::Sqrt(rc*rc - he*he);
	  accphimax = rc*rc*alphap - hp*TMath::Sqrt(rc*rc - hp*hp)-
	    ((TMath::Sqrt(pow(rc,2)-pow(he,2))-hp)*(TMath::Sqrt(pow(rc,2)-pow(hp,2))-he))/2+
	    rc*rc*alphad - rc*rc*TMath::Sin(alphad)*TMath::Cos(alphad);
	}
      }

      if(detamin < eMin && dphimin < pMin ){
	Float_t he = eMax + etaJet[k];
	Float_t hp = phiJet[k] - pMin;
	Float_t rlim = TMath::Sqrt(pow(he,2)+pow(hp,2));
	Float_t alphae = TMath::ACos(he/rc);
	Float_t alphap = TMath::ACos(hp/rc);
	Float_t alphad = (alphae+alphap)/2-TMath::Pi()/4;
	if(rlim <= rc){
	  accetamin = rc*rc*alphae - he*TMath::Sqrt(rc*rc - he*he);
	  accphimin = rc*rc*alphap - hp*TMath::Sqrt(rc*rc - hp*hp);
	}
	if(rlim > rc){
	  accetamin = rc*rc*alphae - he*TMath::Sqrt(rc*rc - he*he);
	  accphimin = rc*rc*alphap - hp*TMath::Sqrt(rc*rc - hp*hp)-
	    ((TMath::Sqrt(pow(rc,2)-pow(he,2))-hp)*(TMath::Sqrt(pow(rc,2)-pow(hp,2))-he))/2+
	    rc*rc*alphad - rc*rc*TMath::Sin(alphad)*TMath::Cos(alphad);
	}
      }
      areaJet[k] = rc*rc*TMath::Pi() - accetamax - accetamin - accphimax - accphimin;
      areaOut = areaOut - areaJet[k];
    } // end loop on jets
   
    //subtract background using area method
    for(Int_t ljet=0; ljet<nJ; ljet++){
      Float_t areaRatio = areaJet[ljet]/areaOut;
      etJet[ljet] = etIn[ljet]-etOut*areaRatio; // subtraction
    }

    // estimate new total background
    Float_t areaT = 2*(header->GetLegoEtaMax()*header->GetLegoPhiMax());
    etbgTotalN = etOut*areaT/areaOut;
  }
  
}

////////////////////////////////////////////////////////////////////////
void AliUA1JetFinderV2::SubtractBackgC(const Int_t& nIn, const Int_t&nJ, Float_t&etbgTotalN,
                      const Float_t* ptT, const Float_t* etaT, const Float_t* phiT,
                      Float_t* const etJet, const Float_t* etaJet, const Float_t* phiJet,
                      Float_t* const etsigJet,Int_t* const multJet, Int_t* const injet)
{
  //background subtraction using cone method but without correction in dE/deta distribution
  
  //calculate energy inside and outside cones
  AliUA1JetHeaderV1* header = (AliUA1JetHeaderV1*) fHeader;
  Float_t rc= header->GetRadius();
  Float_t etIn[30] = {0.};
  Float_t etOut = 0;
  for(Int_t jpart = 0; jpart < nIn; jpart++){ // loop for all particles in array
    // if((fReader->GetCutFlag(jpart)) != 1) continue; // pt cut
    for(Int_t ijet=0; ijet<nJ; ijet++){
      Float_t deta = etaT[jpart] - etaJet[ijet];
      Float_t dphi = phiT[jpart] - phiJet[ijet];
      if (dphi < -TMath::Pi()) dphi= -dphi - 2.0 * TMath::Pi();
      if (dphi > TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
      Float_t dr = TMath::Sqrt(deta * deta + dphi * dphi);
      if(dr <= rc){ // particles inside this cone
	multJet[ijet]++;
	injet[jpart] = ijet;
	if((fReader->GetCutFlag(jpart)) == 1){ // pt cut
	  etIn[ijet] += ptT[jpart];
	  if(fReader->GetSignalFlag(jpart) == 1) etsigJet[ijet]+= ptT[jpart];
	}
	break;
      }
    }// end jets loop
    if(injet[jpart] == -1 && fReader->GetCutFlag(jpart) == 1)
      etOut += ptT[jpart]; // particle outside cones and pt cut
  } //end particle loop
  
  //estimate jets and background areas
  Float_t areaJet[30];
  Float_t areaOut = 4*(header->GetLegoEtaMax())*TMath::Pi();
  for(Int_t k=0; k<nJ; k++){
    Float_t detamax = etaJet[k] + rc;
    Float_t detamin = etaJet[k] - rc;
    Float_t accmax = 0.0; Float_t accmin = 0.0;
    if(detamax > header->GetLegoEtaMax()){ // sector outside etamax
      Float_t h = header->GetLegoEtaMax() - etaJet[k];
      accmax = rc*rc*TMath::ACos(h/rc) - h*TMath::Sqrt(rc*rc - h*h);
    }
    if(detamin < header->GetLegoEtaMin()){ // sector outside etamin
      Float_t h = header->GetLegoEtaMax() + etaJet[k];
      accmin = rc*rc*TMath::ACos(h/rc) - h*TMath::Sqrt(rc*rc - h*h);
    }
    areaJet[k] = rc*rc*TMath::Pi() - accmax - accmin;
    areaOut = areaOut - areaJet[k];
  }
  //subtract background using area method
  for(Int_t ljet=0; ljet<nJ; ljet++){
    Float_t areaRatio = areaJet[ljet]/areaOut;
    etJet[ljet] = etIn[ljet]-etOut*areaRatio; // subtraction
  }
  
  // estimate new total background
  Float_t areaT = 4*(header->GetLegoEtaMax())*TMath::Pi();
  etbgTotalN = etOut*areaT/areaOut;
  
}


////////////////////////////////////////////////////////////////////////
void AliUA1JetFinderV2::SubtractBackgStat(const Int_t& nIn, const Int_t&nJ,Float_t&etbgTotalN,
					  const Float_t* ptT, const Float_t* etaT, const Float_t* phiT, const Int_t* cFlagT, 
					  const Int_t* sFlagT, Float_t* const etJet, const Float_t* etaJet, const Float_t* phiJet,
					  Float_t* const etsigJet, Int_t* const multJet, Int_t* const injet)
{

  //background subtraction using statistical method
  AliUA1JetHeaderV1* header = (AliUA1JetHeaderV1*) fHeader;
  Float_t etbgStat = header->GetBackgStat(); // pre-calculated background
  
  //calculate energy inside
  Float_t rc= header->GetRadius();
  Float_t etIn[30] = {0.};
  
  for(Int_t jpart = 0; jpart < nIn; jpart++)
    { // loop for all particles in array
      
      for(Int_t ijet=0; ijet<nJ; ijet++)
	{
	  Float_t deta = etaT[jpart] - etaJet[ijet];
	  Float_t dphi = phiT[jpart] - phiJet[ijet];
	  if (dphi < -TMath::Pi()) dphi= -dphi - 2.0 * TMath::Pi();
	  if (dphi > TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
	  Float_t dr = TMath::Sqrt(deta * deta + dphi * dphi);
	  if(dr <= rc){ // particles inside this cone
	    multJet[ijet]++;
	    injet[jpart] = ijet;
	    if(cFlagT[jpart] == 1){ // pt cut
	      etIn[ijet]+= ptT[jpart];
	      if(sFlagT[jpart] == 1) etsigJet[ijet] += ptT[jpart];
	    }
	    break;
	  }
	}// end jets loop
    } //end particle loop
  
  //calc jets areas
  Float_t areaJet[30];
  Float_t areaOut = 4*(header->GetLegoEtaMax())*TMath::Pi();
  for(Int_t k=0; k<nJ; k++)
    {
      Float_t detamax = etaJet[k] + rc;
      Float_t detamin = etaJet[k] - rc;
      Float_t accmax = 0.0; Float_t accmin = 0.0;
      if(detamax > header->GetLegoEtaMax()){ // sector outside etamax
	Float_t h = header->GetLegoEtaMax() - etaJet[k];
	accmax = rc*rc*TMath::ACos(h/rc) - h*TMath::Sqrt(rc*rc - h*h);
      }
      if(detamin < header->GetLegoEtaMin()){ // sector outside etamin
      Float_t h = header->GetLegoEtaMax() + etaJet[k];
      accmin = rc*rc*TMath::ACos(h/rc) - h*TMath::Sqrt(rc*rc - h*h);
      }
      areaJet[k] = rc*rc*TMath::Pi() - accmax - accmin;
    }

  //subtract background using area method
  for(Int_t ljet=0; ljet<nJ; ljet++){
    Float_t areaRatio = areaJet[ljet]/areaOut;
    etJet[ljet] = etIn[ljet]-etbgStat*areaRatio; // subtraction
  }
  
  etbgTotalN = etbgStat;
}

////////////////////////////////////////////////////////////////////////
void AliUA1JetFinderV2::SubtractBackgCone(const Int_t& nIn, const Int_t&nJ,Float_t& etbgTotalN,
                      const Float_t* ptT, const Float_t* etaT, const Float_t* phiT, const Int_t* cFlagT, const Int_t* sFlagT,
                      Float_t* const etJet, const Float_t* etaJet, const Float_t* phiJet,
                      Float_t* const etsigJet, Int_t* const multJet, Int_t* const injet)
{
  // Cone background subtraction method taking into acount dEt/deta distribution
  AliUA1JetHeaderV1* header = (AliUA1JetHeaderV1*) fHeader;
  //general
  Float_t rc= header->GetRadius();
  Float_t etamax = header->GetLegoEtaMax();
  Float_t etamin = header->GetLegoEtaMin();
  Int_t ndiv = 100;
  
  // jet energy and area arrays
  TH1F* hEtJet[30];
  TH1F* hAreaJet[30];
  for(Int_t mjet=0; mjet<nJ; mjet++){
    char hEtname[256]; char hAreaname[256];
    snprintf(hEtname, 256, "hEtJet%d", mjet); snprintf(hAreaname, 256, "hAreaJet%d", mjet);
    hEtJet[mjet] = new TH1F(hEtname,"et dist in eta ",ndiv,etamin,etamax);
    hAreaJet[mjet] = new TH1F(hAreaname,"area dist in eta ",ndiv,etamin,etamax);
  }
  // background energy and area
  TH1F* hEtBackg = new TH1F("hEtBackg"," backg et dist in eta ",ndiv,etamin,etamax);
  TH1F* hAreaBackg = new TH1F("hAreaBackg","backg area dist in eta ",ndiv,etamin,etamax);

  //fill energies
  for(Int_t jpart = 0; jpart < nIn; jpart++){ // loop for all particles in array
    for(Int_t ijet=0; ijet<nJ; ijet++){  // loop for all jets
      Float_t deta = etaT[jpart] - etaJet[ijet];
      Float_t dphi = phiT[jpart] - phiJet[ijet];
      if (dphi < -TMath::Pi()) dphi= -dphi - 2.0 * TMath::Pi();
      if (dphi > TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
      Float_t dr = TMath::Sqrt(deta * deta + dphi * dphi);
      if(dr <= rc){ // particles inside this cone
	injet[jpart] = ijet;
	multJet[ijet]++;
	if(cFlagT[jpart] == 1){// pt cut
	  hEtJet[ijet]->Fill(etaT[jpart],ptT[jpart]); //particle inside cone
	  if(sFlagT[jpart] == 1) etsigJet[ijet] += ptT[jpart];
	}
	break;
      }
    }// end jets loop

    if(injet[jpart] == -1  && cFlagT[jpart] == 1)
      hEtBackg->Fill(etaT[jpart],ptT[jpart]); // particle outside cones
  } //end particle loop

  //calc areas
  Float_t eta0 = etamin;
  Float_t etaw = (etamax - etamin)/((Float_t)ndiv);
  Float_t eta1 = eta0 + etaw;
  for(Int_t etabin = 0; etabin< ndiv; etabin++){ // loop for all eta bins
    Float_t etac = eta0 + etaw/2.0;
    Float_t areabg = etaw*2.0*TMath::Pi();
    for(Int_t ijet=0; ijet<nJ; ijet++){  // loop for all jets
      Float_t deta0 = TMath::Abs(eta0 - etaJet[ijet]);
      Float_t deta1 = TMath::Abs(eta1 - etaJet[ijet]);
      Float_t acc0 = 0.0; Float_t acc1 = 0.0;
      Float_t areaj = 0.0;
      if(deta0 > rc && deta1 < rc){
	acc1 = rc*rc*TMath::ACos(deta1/rc) - deta1*TMath::Sqrt(rc*rc - deta1*deta1);
	areaj = acc1;
      }
      if(deta0 < rc && deta1 > rc){
	acc0 = rc*rc*TMath::ACos(deta0/rc) - deta0*TMath::Sqrt(rc*rc - deta0*deta0);
	areaj = acc0;
      }
      if(deta0 < rc && deta1 < rc){
	acc0 = rc*rc*TMath::ACos(deta0/rc) - deta0*TMath::Sqrt(rc*rc - deta0*deta0);
	acc1 = rc*rc*TMath::ACos(deta1/rc) - deta1*TMath::Sqrt(rc*rc - deta1*deta1);
	if(eta1<etaJet[ijet]) areaj = acc1-acc0;  // case 1
	if((eta0 < etaJet[ijet]) && (etaJet[ijet]<eta1)) areaj = rc*rc*TMath::Pi() - acc1 -acc0; // case 2
	if(etaJet[ijet] < eta0) areaj = acc0 -acc1; // case 3
      }
      hAreaJet[ijet]->Fill(etac,areaj);
      areabg = areabg - areaj;
    } // end jets loop
    hAreaBackg->Fill(etac,areabg);
    eta0 = eta1;
    eta1 = eta1 + etaw;
  } // end loop for all eta bins

  //subtract background
  for(Int_t kjet=0; kjet<nJ; kjet++){
    etJet[kjet] = 0.0; // first  clear etJet for this jet
    for(Int_t bin = 0; bin< ndiv; bin++){
      if(hAreaJet[kjet]->GetBinContent(bin)){
	Float_t areab = hAreaBackg->GetBinContent(bin);
	Float_t etb = hEtBackg->GetBinContent(bin);
	Float_t areaR = (hAreaJet[kjet]->GetBinContent(bin))/areab;
	etJet[kjet] = etJet[kjet] + ((hEtJet[kjet]->GetBinContent(bin)) - etb*areaR); //subtraction
      }
    }
  }

  // calc background total
  Double_t etOut = hEtBackg->Integral();
  Double_t areaOut = hAreaBackg->Integral();
  Float_t areaT = 4*(header->GetLegoEtaMax())*TMath::Pi();
  etbgTotalN = etOut*areaT/areaOut;
  
  //delete
  for(Int_t ljet=0; ljet<nJ; ljet++){  // loop for all jets
    delete hEtJet[ljet];
    delete hAreaJet[ljet];
  }

  delete hEtBackg;
  delete hAreaBackg;
}

////////////////////////////////////////////////////////////////////////
void AliUA1JetFinderV2::SubtractBackgRatio(const Int_t& nIn, const Int_t&nJ,Float_t& etbgTotalN,
                      const Float_t* ptT, const Float_t* etaT, const Float_t* phiT, const Int_t* cFlagT, const Int_t* sFlagT,
                      Float_t* const etJet, const Float_t* etaJet, const Float_t* phiJet,
                      Float_t* const etsigJet, Int_t* const multJet, Int_t* const injet)
{
  // Ratio background subtraction method taking into acount dEt/deta distribution
  AliUA1JetHeaderV1* header = (AliUA1JetHeaderV1*) fHeader;
  //factor F calc before
  Float_t bgRatioCut = header->GetBackgCutRatio();
  
  //general
  Float_t rc= header->GetRadius();
  Float_t etamax = header->GetLegoEtaMax();
  Float_t etamin = header->GetLegoEtaMin();
  Int_t ndiv = 100;
  
  // jet energy and area arrays
  TH1F* hEtJet[30];
  TH1F* hAreaJet[30];
  for(Int_t mjet=0; mjet<nJ; mjet++){
    char hEtname[256]; char hAreaname[256];
    snprintf(hEtname, 256, "hEtJet%d", mjet); snprintf(hAreaname, 256, "hAreaJet%d", mjet);
    hEtJet[mjet] = new TH1F(hEtname,"et dist in eta ",ndiv,etamin,etamax);        // change range
    hAreaJet[mjet] = new TH1F(hAreaname,"area dist in eta ",ndiv,etamin,etamax);  // change range
  }
  // background energy and area
  TH1F* hEtBackg = new TH1F("hEtBackg"," backg et dist in eta ",ndiv,etamin,etamax);         // change range
  TH1F* hAreaBackg = new TH1F("hAreaBackg","backg area dist in eta ",ndiv,etamin,etamax);  // change range

  //fill energies
  for(Int_t jpart = 0; jpart < nIn; jpart++){ // loop for all particles in array
    for(Int_t ijet=0; ijet<nJ; ijet++){  // loop for all jets
      Float_t deta = etaT[jpart] - etaJet[ijet];
      Float_t dphi = phiT[jpart] - phiJet[ijet];
      if (dphi < -TMath::Pi()) dphi= -dphi - 2.0 * TMath::Pi();
      if (dphi > TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
      Float_t dr = TMath::Sqrt(deta * deta + dphi * dphi);
      if(dr <= rc){ // particles inside this cone
	multJet[ijet]++;
	injet[jpart] = ijet;
	if(cFlagT[jpart] == 1){ //pt cut
	  hEtJet[ijet]->Fill(etaT[jpart],ptT[jpart]); //particle inside cone and pt cut
	  if(sFlagT[jpart] == 1) etsigJet[ijet] += ptT[jpart];
	}
	break;
      }
    }// end jets loop
    if(injet[jpart] == -1) hEtBackg->Fill(etaT[jpart],ptT[jpart]); // particle outside cones
  } //end particle loop

  //calc areas
  Float_t eta0 = etamin;
  Float_t etaw = (etamax - etamin)/((Float_t)ndiv);
  Float_t eta1 = eta0 + etaw;
  for(Int_t etabin = 0; etabin< ndiv; etabin++){ // loop for all eta bins
    Float_t etac = eta0 + etaw/2.0;
    Float_t areabg = etaw*2.0*TMath::Pi();
    for(Int_t ijet=0; ijet<nJ; ijet++){  // loop for all jets
      Float_t deta0 = TMath::Abs(eta0 - etaJet[ijet]);
      Float_t deta1 = TMath::Abs(eta1 - etaJet[ijet]);
      Float_t acc0 = 0.0; Float_t acc1 = 0.0;
      Float_t areaj = 0.0;
      if(deta0 > rc && deta1 < rc){
	acc1 = rc*rc*TMath::ACos(deta1/rc) - deta1*TMath::Sqrt(rc*rc - deta1*deta1);
	areaj = acc1;
      }
      if(deta0 < rc && deta1 > rc){
	acc0 = rc*rc*TMath::ACos(deta0/rc) - deta0*TMath::Sqrt(rc*rc - deta0*deta0);
	areaj = acc0;
      }
      if(deta0 < rc && deta1 < rc){
	acc0 = rc*rc*TMath::ACos(deta0/rc) - deta0*TMath::Sqrt(rc*rc - deta0*deta0);
	acc1 = rc*rc*TMath::ACos(deta1/rc) - deta1*TMath::Sqrt(rc*rc - deta1*deta1);
	if(eta1<etaJet[ijet]) areaj = acc1-acc0;  // case 1
	if((eta0 < etaJet[ijet]) && (etaJet[ijet]<eta1)) areaj = rc*rc*TMath::Pi() - acc1 -acc0; // case 2
	if(etaJet[ijet] < eta0) areaj = acc0 -acc1; // case 3
      }
      hAreaJet[ijet]->Fill(etac,areaj);
      areabg = areabg - areaj;
    } // end jets loop
    hAreaBackg->Fill(etac,areabg);
    eta0 = eta1;
    eta1 = eta1 + etaw;
  } // end loop for all eta bins

  //subtract background
  for(Int_t kjet=0; kjet<nJ; kjet++){
    etJet[kjet] = 0.0; // first  clear etJet for this jet
    for(Int_t bin = 0; bin< ndiv; bin++){
      if(hAreaJet[kjet]->GetBinContent(bin)){
	Float_t areab = hAreaBackg->GetBinContent(bin);
	Float_t etb = hEtBackg->GetBinContent(bin);
	Float_t areaR = (hAreaJet[kjet]->GetBinContent(bin))/areab;
	etJet[kjet] = etJet[kjet] + ((hEtJet[kjet]->GetBinContent(bin)) - etb*areaR*bgRatioCut); //subtraction
      }
    }
  }

  // calc background total
  Double_t etOut = hEtBackg->Integral();
  Double_t areaOut = hAreaBackg->Integral();
  Float_t areaT = 4*(header->GetLegoEtaMax())*TMath::Pi();
  etbgTotalN = etOut*areaT/areaOut;
  
  //delete
  for(Int_t ljet=0; ljet<nJ; ljet++){  // loop for all jets
    delete hEtJet[ljet];
    delete hAreaJet[ljet];
  }
  
  delete hEtBackg;
  delete hAreaBackg;
}

////////////////////////////////////////////////////////////////////////
void AliUA1JetFinderV2::Reset()
{
  fLego->Reset();
  AliJetFinder::Reset();
}

////////////////////////////////////////////////////////////////////////
void AliUA1JetFinderV2::WriteJHeaderToFile() const
{
  AliUA1JetHeaderV1* header = (AliUA1JetHeaderV1*) fHeader;
  header->Write();
}

////////////////////////////////////////////////////////////////////////
void AliUA1JetFinderV2::InitTask(TChain* tree)
{

  // initializes some variables
  AliUA1JetHeaderV1* header = (AliUA1JetHeaderV1*) fHeader;
  // book lego
  fLego = new TH2F("legoH","eta-phi",
		   header->GetLegoNbinEta(), header->GetLegoEtaMin(),
		   header->GetLegoEtaMax(),  header->GetLegoNbinPhi(),
		   header->GetLegoPhiMin(),  header->GetLegoPhiMax());
  
  fDebug = fHeader->GetDebug();
  fOpt = fReader->GetReaderHeader()->GetDetector();
  
  // Tasks initialization
  if(fOpt>0)
    fReader->CreateTasks(tree);

}
