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

/* $Id */
 
//---------------------------------------------------------------------
// UA1 Cone Algorithm Jet finder for jet studies
// manages the search for jets using charged particle momentum and 
// neutral cell energy information
// Authors: Rafael.Diaz.Valdes@cern.ch 
//          magali.estienne@subatech.in2p3.fr
//          alexandre.shabetai@cern.ch
// ** 2011
// Modified accordingly to reader/finder splitting and new handling of neutral information
// Versions V1 and V2 merged
//---------------------------------------------------------------------

#include <TH2F.h>
#include <TMath.h>

#include "AliUA1JetFinder.h"
#include "AliUA1JetHeaderV1.h"
#include "AliJetCalTrk.h"
#include "AliJetBkg.h"
#include "AliAODJetEventBackground.h"
#include "AliAODJet.h"

ClassImp(AliUA1JetFinder)

////////////////////////////////////////////////////////////////////////

AliUA1JetFinder::AliUA1JetFinder():
  AliJetFinder(),
  fLego(0),  
  fJetBkg(new AliJetBkg())
{
  // Default constructor
}

//-----------------------------------------------------------------------
AliUA1JetFinder::~AliUA1JetFinder()
{
  // Destructor
  delete fLego;
  delete fJetBkg;

}

//-----------------------------------------------------------------------
void AliUA1JetFinder::FindJets()
{
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
  
  AliUA1JetHeader* header  = (AliUA1JetHeader*) fHeader;
  Int_t              nIn = fCalTrkEvent->GetNCalTrkTracks();
  fDebug   = fHeader->GetDebug();
  
  if (nIn <= 0) return;
  fJetBkg->SetHeader(fHeader);
  fJetBkg->SetCalTrkEvent(GetCalTrkEvent());
  fJetBkg->SetDebug(fDebug);
  // local arrays for input
  // ToDo: check memory fragmentation, maybe better to
  // define them globally and resize as needed
  // Fragmentation should be worse for low mult...
  Float_t* ptT      = new Float_t[nIn];
  Float_t* etaT     = new Float_t[nIn];
  Float_t* phiT     = new Float_t[nIn];
  Int_t*   injet    = new Int_t[nIn];
  Int_t*   injetOk  = new Int_t[nIn];

  memset(ptT,0,sizeof(Float_t)*nIn);
  memset(etaT,0,sizeof(Float_t)*nIn);
  memset(phiT,0,sizeof(Float_t)*nIn);
  memset(injet,0,sizeof(Int_t)*nIn);
  memset(injetOk,-1,sizeof(Int_t)*nIn);

  // load input vectors and calculate total energy in array

  // total energy in array
  Float_t etbgTotal = 0.;
  Float_t npart = 0.;
  Float_t etbg2 = 0.;

  for (Int_t i = 0; i < fCalTrkEvent->GetNCalTrkTracks(); i++){
    ptT[i]  = fCalTrkEvent->GetCalTrkTrack(i)->GetPt();
    etaT[i] = fCalTrkEvent->GetCalTrkTrack(i)->GetEta();
    phiT[i] = ((fCalTrkEvent->GetCalTrkTrack(i)->GetPhi() < 0) ? (fCalTrkEvent->GetCalTrkTrack(i)->GetPhi()) + 2 * TMath::Pi() : fCalTrkEvent->GetCalTrkTrack(i)->GetPhi());
    //fCalTrkEvent->GetCalTrkTrack(i)->Print(Form("%d",i));
    if (fCalTrkEvent->GetCalTrkTrack(i)->GetCutFlag() != 1) continue;
    fLego->Fill(etaT[i], phiT[i], ptT[i]);
    npart += 1;
    etbgTotal+= ptT[i];
    etbg2 += ptT[i]*ptT[i];
  }
  
  // calculate total energy and fluctuation in map
  Double_t meanpt = 0.;
  Double_t ptRMS = 0.;
  if(npart>0){
    meanpt = etbgTotal/npart;
    etbg2 = etbg2/npart;
    if(etbg2>(meanpt*meanpt)){// prenent NAN, should only happen due to numerical instabilities
      ptRMS = TMath::Sqrt(etbg2-meanpt*meanpt);
    }
  }
  Double_t dEtTotal = (TMath::Sqrt(npart))*TMath::Sqrt(meanpt * meanpt + ptRMS*ptRMS);
  
  // arrays to hold jets
  Float_t etaJet[kMaxJets];
  Float_t phiJet[kMaxJets];
  Float_t etJet[kMaxJets];
  Float_t etsigJet[kMaxJets]; //signal et in jet
  Float_t etallJet[kMaxJets];  // total et in jet (tmp variable)
  Int_t   ncellsJet[kMaxJets];
  Int_t   multJetT[kMaxJets];
  Int_t   multJetC[kMaxJets];
  Int_t   multJet[kMaxJets];
  Float_t *areaJet = new Float_t[kMaxJets];
  // Used for jet reordering at the end of the jet finding procedure
  Float_t etaJetOk[kMaxJets];
  Float_t phiJetOk[kMaxJets];
  Float_t etJetOk[kMaxJets];
  Float_t etsigJetOk[kMaxJets]; //signal et in jet
  Float_t etallJetOk[kMaxJets];  // total et in jet (tmp variable)
  Int_t   ncellsJetOk[kMaxJets];
  Int_t   multJetOk[kMaxJets];
  Float_t *areaJetOk = new Float_t[kMaxJets];
  Int_t   nJets; // to hold number of jets found by algorithm
  Int_t   nj;    // number of jets accepted
  Float_t prec  = header->GetPrecBg();
  Float_t bgprec = 1;

  while(bgprec > prec){
    //reset jet arrays in memory
    memset(etaJet,0,sizeof(Float_t)*kMaxJets);
    memset(phiJet,0,sizeof(Float_t)*kMaxJets);
    memset(etJet,0,sizeof(Float_t)*kMaxJets);
    memset(etallJet,0,sizeof(Float_t)*kMaxJets);
    memset(etsigJet,0,sizeof(Float_t)*kMaxJets);
    memset(ncellsJet,0,sizeof(Int_t)*kMaxJets);
    memset(multJetT,0,sizeof(Int_t)*kMaxJets);
    memset(multJetC,0,sizeof(Int_t)*kMaxJets);
    memset(multJet,0,sizeof(Int_t)*kMaxJets);
    memset(areaJet,0,sizeof(Float_t)*kMaxJets);
    memset(etaJetOk,0,sizeof(Float_t)*kMaxJets);
    memset(phiJetOk,0,sizeof(Float_t)*kMaxJets);
    memset(etJetOk,0,sizeof(Float_t)*kMaxJets);
    memset(etallJetOk,0,sizeof(Float_t)*kMaxJets);
    memset(etsigJetOk,0,sizeof(Float_t)*kMaxJets);
    memset(ncellsJetOk,0,sizeof(Int_t)*kMaxJets);
    memset(multJetOk,0,sizeof(Int_t)*kMaxJets);
    memset(areaJetOk,0,sizeof(Float_t)*kMaxJets);
    nJets = 0;
    nj = 0;
    // reset particles-jet array in memory
    memset(injet,-1,sizeof(Int_t)*nIn);
    //run cone algorithm finder
    RunAlgoritm(etbgTotal,dEtTotal,nJets,etJet,etaJet,phiJet,etallJet,ncellsJet);
    //run background subtraction
    if(nJets > header->GetNAcceptJets()) // limited number of accepted jets per event
      nj = header->GetNAcceptJets();
    else
      nj = nJets;

    //subtract background
    Float_t etbgTotalN = 0.0; //new background
    Float_t sigmaN = 0.0; //new background
    if(header->GetBackgMode() == 1) {// standard
      fJetBkg->SubtractBackg(nIn,nj,etbgTotalN,sigmaN,ptT,etaT,phiT,etJet,etaJet,phiJet,etsigJet,multJetT,multJetC,multJet,injet,areaJet);}
    if(header->GetBackgMode() == 2) //cone
      fJetBkg->SubtractBackgCone(nIn,nj,etbgTotalN,sigmaN,ptT,etaT,phiT,etJet,etaJet,phiJet,etsigJet,multJetT,multJetC,multJet,injet,areaJet);
    if(header->GetBackgMode() == 3) //ratio
      fJetBkg->SubtractBackgRatio(nIn,nj,etbgTotalN,sigmaN,ptT,etaT,phiT,etJet,etaJet,phiJet,etsigJet,multJetT,multJetC,multJet,injet,areaJet);
    if(header->GetBackgMode() == 4) //statistic
      fJetBkg->SubtractBackgStat(nIn,nj,etbgTotalN,sigmaN,ptT,etaT,phiT,etJet,etaJet,phiJet,etsigJet,multJetT,multJetC,multJet,injet,areaJet);

    //calc precision
    if(etbgTotalN != 0.0)
      bgprec = (etbgTotal - etbgTotalN)/etbgTotalN;
    else
      bgprec = 0;
    etbgTotal = etbgTotalN; // update with new background estimation

  } //end while

  // add tracks to the jet if it wasn't yet done
  if (header->GetBackgMode() == 0){
    Float_t rc= header->GetRadius();
    for(Int_t jpart = 0; jpart < nIn; jpart++){ // loop for all particles in array
      for(Int_t ijet=0; ijet<nJets; ijet++){
        Float_t deta = etaT[jpart] - etaJet[ijet];
        Float_t dphi = phiT[jpart] - phiJet[ijet];
        if (dphi < -TMath::Pi()) dphi= -dphi - 2.0 * TMath::Pi();
        if (dphi >  TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
        Float_t dr = TMath::Sqrt(deta * deta + dphi * dphi);
        if(dr <= rc){ // particles inside this cone
	  injet[jpart] = ijet;
          break;
        }
      }// end jets loop
    } //end particle loop
  }

  // add jets to list
  if (fDebug>1) printf("Found %d jets \n", nj);
  
  // Reorder jets by et in cone
  // Sort jets by energy
  Int_t idx[kMaxJets];
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
      areaJetOk[p]   = areaJet[idx[p]];
    }

  //////////////////////////

  Int_t nTracks = fCalTrkEvent->GetNCalTrkTracks();
  
  for(Int_t kj=0; kj<nj; kj++)
    {
      if ((etaJetOk[kj] > (header->GetJetEtaMax())) ||
 	  (etaJetOk[kj] < (header->GetJetEtaMin())) ||
 	  (phiJetOk[kj] > (header->GetJetPhiMax())) ||
 	  (phiJetOk[kj] < (header->GetJetPhiMin())) ||
 	  (etJetOk[kj] < header->GetMinJetEt())) continue; // acceptance eta range and etmin
      Float_t px=-999, py=-999 ,pz=-999 ,en=-999; // convert to 4-vector
      px = etJetOk[kj] * TMath::Cos(phiJetOk[kj]);
      py = etJetOk[kj] * TMath::Sin(phiJetOk[kj]);
      pz = etJetOk[kj] / TMath::Tan(2.0 * TMath::ATan(TMath::Exp(-etaJetOk[kj])));
      en = TMath::Sqrt(px * px + py * py + pz * pz);
      AliAODJet jet(px, py, pz, en);

      // Calc jet area if it wasn't yet done
      if (header->GetBackgMode() == 0){
        // calculate the area of the jet
        Float_t rc= header->GetRadius();
        areaJetOk[kj] = fJetBkg->CalcJetAreaEtaCut(rc,etaJetOk[kj]);
      }

      jet.SetEffArea(areaJetOk[kj],0.,0.,0.);
      jet.SetBgEnergy(etbgTotal,0.);
      if (fDebug>1) jet.Print(Form("%d",kj));
      
      for(Int_t jpart = 0; jpart < nTracks; jpart++) { // loop for all particles in array
        // Track to jet reordering
        if(injet[jpart] == idx[kj]){
          injetOk[jpart] = kj;
        }
        // Check if the particle belongs to the jet and add the ref
        if(injetOk[jpart] == kj && fCalTrkEvent->GetCalTrkTrack(jpart)->GetCutFlag() == 1) {
          jet.AddTrack(fCalTrkEvent->GetCalTrkTrack(jpart)->GetTrackObject());
	}
      }

      AddJet(jet);
      
    }

  //delete
  delete[] ptT;
  delete[] etaT;
  delete[] phiT;
  delete[] injet;
  delete[] injetOk;
  delete[] areaJet;
  delete[] areaJetOk;

}

//-----------------------------------------------------------------------
void AliUA1JetFinder::RunAlgoritm(Float_t etbgTotal, Double_t dEtTotal, Int_t& nJets,
				  Float_t* const etJet,Float_t* const etaJet, Float_t* const phiJet,
				  Float_t* const etallJet, Int_t* const ncellsJet)
{
  // Dump lego
  AliUA1JetHeader* header = (AliUA1JetHeader*) fHeader;
  const Int_t nBinsMax = 120000; // we use a fixed array not to fragment memory 
 
  const Int_t nBinEta = header->GetLegoNbinEta();
  const Int_t nBinPhi = header->GetLegoNbinPhi();
  if((nBinPhi*nBinEta)>nBinsMax){
    AliError("Too many bins of the ETA-PHI histogram");
  }

  Float_t etCell[nBinsMax] = {0.};   //! Cell Energy
  Float_t etaCell[nBinsMax] = {0.};  //! Cell eta
  Float_t phiCell[nBinsMax] = {0.};  //! Cell phi
  Short_t flagCell[nBinsMax] = {0}; //! Cell flag
  
  Int_t nCell = 0;
  TAxis* xaxis = fLego->GetXaxis();
  TAxis* yaxis = fLego->GetYaxis();
  Float_t e = 0.0;
  for (Int_t i = 1; i <= nBinEta; i++) {
    for (Int_t j = 1; j <= nBinPhi; j++) {
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
  Float_t etaAlgoJet[kMaxJets] = {0.0};
  Float_t phiAlgoJet[kMaxJets] = {0.0};
  Float_t etAlgoJet[kMaxJets] = {0.0};
  Int_t   ncellsAlgoJet[kMaxJets] = {0};

  // Run algorithm//
  
  // Sort cells by et
  Int_t index[nBinsMax];
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
        if(nJets<kMaxJets){
	  etaAlgoJet[nJets]    = eta;
	  phiAlgoJet[nJets]    = phi;
	  etAlgoJet[nJets]     = etCone;
	  ncellsAlgoJet[nJets] = nCellIn;
	  nJets++;
        }
	else{
	  AliError(Form("Too many jets (> %d) found by UA1JetFinder, adapt your cuts",kMaxJets));
	  break;
	}
      }      
    } // end of cells loop

  //reorder jets by et in cone
  //sort jets by energy
  for(Int_t p = 0; p < nJets; p++)
    {
      etaJet[p]    = etaAlgoJet[p];
      phiJet[p]    = phiAlgoJet[p];
      etJet[p]     = etAlgoJet[p];
      etallJet[p]  = etAlgoJet[p];
      ncellsJet[p] = ncellsAlgoJet[p];
    }

}

//-----------------------------------------------------------------------
void AliUA1JetFinder::Reset()
{
  fLego->Reset();
  AliJetFinder::Reset();

}

//-----------------------------------------------------------------------
void AliUA1JetFinder::WriteJHeaderToFile() const
{
  AliUA1JetHeader* header = (AliUA1JetHeader*) fHeader;
  header->Write();

}

//-----------------------------------------------------------------------
void AliUA1JetFinder::Init()
{

  // initializes some variables
  AliUA1JetHeader* header = (AliUA1JetHeader*) fHeader;
  // book lego
  fLego = new TH2F("legoH","eta-phi",
		   header->GetLegoNbinEta(), header->GetLegoEtaMin(),
		   header->GetLegoEtaMax(),  header->GetLegoNbinPhi(),
		   header->GetLegoPhiMin(),  header->GetLegoPhiMax());
  
}

