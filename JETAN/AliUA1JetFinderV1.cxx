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
// UA1 Cone Algorithm Jet finder
// manages the search for jets
// Author: Rafael.Diaz.Valdes@cern.ch
// (version in c++)
//---------------------------------------------------------------------

#include <TArrayF.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>

#include "AliUA1JetFinderV1.h"
#include "AliUA1JetHeaderV1.h"
#include "AliJetReaderHeader.h"
#include "AliJetReader.h"
#include "AliJetHeader.h"


#include "AliAODJet.h"
#include "AliLog.h"


ClassImp(AliUA1JetFinderV1)

/////////////////////////////////////////////////////////////////////

AliUA1JetFinderV1::AliUA1JetFinderV1() :
    AliJetFinder(),
  fLego(0),
  fhEtBackg(0),
  fhAreaBackg(0)
{
  // Constructor
  for(int i = 0;i < kMaxJets;i++){
    fhAreaJet[i] = fhEtJet[i] = 0;
  }
}

////////////////////////////////////////////////////////////////////////

AliUA1JetFinderV1::~AliUA1JetFinderV1()

{
  // destructor
  delete fLego;
  fLego = 0;
  if(fhEtBackg)delete  fhEtBackg;
  fhEtBackg = 0;
  if( fhAreaBackg) delete  fhAreaBackg;
  fhAreaBackg = 0;
  for(int i = 0;i < kMaxJets;i++){
    if(fhAreaJet[i])delete fhAreaJet[i];
    if(fhEtJet[i]) delete fhEtJet[i];
    fhAreaJet[i] = fhEtJet[i] = 0;
  }

}

////////////////////////////////////////////////////////////////////////


void AliUA1JetFinderV1::FindJets()

{
  //1) Fill cell map array
  //2) calculate total energy and fluctuation level
  //3) Run algorithm
  //   3.1) look centroides in cell map
  //   3.2) calculate total energy in cones
  //   3.3) flag as a possible jet
  //   3.4) reorder cones by energy
  //4) subtract backg in accepted jets
  //5) fill AliJet list

  // transform input to pt,eta,phi plus lego
    
  AliUA1JetHeaderV1* header = (AliUA1JetHeaderV1*) fHeader;
  TClonesArray *lvArray = fReader->GetMomentumArray();
  Int_t nIn =  lvArray->GetEntries();
  if (nIn <= 0) return;

  // local arrays for input
  // ToDo: check memory fragmentation, maybe better to 
  // define them globally and resize as needed
  // Fragementation should be worse for low mult...
  Float_t* ptT   = new Float_t[nIn];
  Float_t* etaT  = new Float_t[nIn];
  Float_t* phiT  = new Float_t[nIn];
  Int_t*   injet = new Int_t[nIn];

  memset(ptT,0,sizeof(Float_t)*nIn);
  memset(etaT,0,sizeof(Float_t)*nIn);
  memset(phiT,0,sizeof(Float_t)*nIn);


  // load input vectors and calculate total energy in array

  //total energy in array
  Float_t  etbgTotal = 0.0;
  Float_t npart = 0;
  Float_t etbg2 = 0;

  for (Int_t i = 0; i < nIn; i++){
    TLorentzVector *lv = (TLorentzVector*) lvArray->At(i);
    ptT[i]  = lv->Pt();
    etaT[i] = lv->Eta();
    phiT[i] = ((lv->Phi() < 0) ? (lv->Phi()) + 2 * TMath::Pi() : lv->Phi());
    if (fReader->GetCutFlag(i) != 1) continue;
    fLego ->Fill(etaT[i], phiT[i], ptT[i]);
    npart += 1;
    etbgTotal+= ptT[i];
    etbg2 += ptT[i]*ptT[i];
  }

  // calculate total energy and fluctuation in map
  Double_t meanpt = 0;
  Double_t ptRMS = 0;
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
  Int_t ncellsJet[kMaxJets];
  Int_t multJet[kMaxJets];
  Int_t nJets; // to hold number of jets found by algorithm
  Int_t nj;    // number of jets accepted
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
     memset(multJet,0,sizeof(Int_t)*kMaxJets);
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
     if(header->GetBackgMode() == 1) // standar
        SubtractBackg(nIn,nj,etbgTotalN,ptT,etaT,phiT,etJet,etaJet,phiJet,etsigJet,multJet,injet);
     if(header->GetBackgMode() == 2) //cone
        SubtractBackgCone(nIn,nj,etbgTotalN,ptT,etaT,phiT,etJet,etaJet,phiJet,etsigJet,multJet,injet);
     if(header->GetBackgMode() == 3) //ratio
        SubtractBackgRatio(nIn,nj,etbgTotalN,ptT,etaT,phiT,etJet,etaJet,phiJet,etsigJet,multJet,injet);
     if(header->GetBackgMode() == 4) //statistic
        SubtractBackgStat(nIn,nj,etbgTotalN,ptT,etaT,phiT,etJet,etaJet,phiJet,etsigJet,multJet,injet);
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
      for(Int_t ijet=0; ijet<nj; ijet++){
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
  Int_t idxjets[kMaxJets];
  Int_t nselectj = 0;

  TRefArray *refs = 0;
  Bool_t fromAod = !strcmp(fReader->ClassName(),"AliJetAODReader");
  if (fromAod) refs = fReader->GetReferences();
  Float_t rc= header->GetRadius();
  for(Int_t kj=0; kj<nj; kj++){
     if ((etaJet[kj] > (header->GetJetEtaMax())) ||
          (etaJet[kj] < (header->GetJetEtaMin())) ||
          (etJet[kj] < header->GetMinJetEt())) continue; // acceptance eta range and etmin
      Float_t px, py,pz,en; // convert to 4-vector
      px = etJet[kj] * TMath::Cos(phiJet[kj]);
      py = etJet[kj] * TMath::Sin(phiJet[kj]);
      pz = etJet[kj] / TMath::Tan(2.0 * TMath::ATan(TMath::Exp(-etaJet[kj])));
      en = TMath::Sqrt(px * px + py * py + pz * pz);
      
      AliAODJet jet(px, py, pz, en);

      if (fromAod){
        for(Int_t jpart = 0; jpart < nIn; jpart++) // loop for all particles in array
          if (injet[jpart] == kj && fReader->GetCutFlag(jpart) == 1)
		    jet.AddTrack(refs->At(jpart));  // check if the particle belongs to the jet and add the ref
      }
      
      //jet.Print("");
      
      // calculate the area of the jet
      Float_t detamax = etaJet[kj] + rc;
      Float_t detamin = etaJet[kj] - rc;
      Float_t accmax = 0.0; Float_t accmin = 0.0;
      if(detamax > header->GetLegoEtaMax()){ // sector outside etamax
         Float_t h = header->GetLegoEtaMax() - etaJet[kj];
         accmax = rc*rc*TMath::ACos(h/rc) - h*TMath::Sqrt(rc*rc - h*h);
      }
      if(detamin < header->GetLegoEtaMin()){ // sector outside etamin
         Float_t h = header->GetLegoEtaMax() + etaJet[kj];
         accmin = rc*rc*TMath::ACos(h/rc) - h*TMath::Sqrt(rc*rc - h*h);
      }
      Float_t areaJet = rc*rc*TMath::Pi() - accmax - accmin;
      // set both areas
      jet.SetEffArea(areaJet,areaJet);

      AddJet(jet);
      
      idxjets[nselectj] = kj;
      nselectj++;
  } //end particle loop

  //add signal percentage and total signal  in AliJets for analysis tool
  Float_t percentage[kMaxJets];
  Int_t ncells[kMaxJets];
  Int_t mult[kMaxJets];
  for(Int_t i = 0; i< nselectj; i++){
     percentage[i] = etsigJet[idxjets[i]]/etJet[idxjets[i]];
     ncells[i] = ncellsJet[idxjets[i]];
     mult[i] = multJet[idxjets[i]];
  }
   //add particle-injet relationship ///
   for(Int_t bj = 0; bj < nIn; bj++){
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
  delete [] etaT;
  delete [] phiT;
  delete [] injet;
}

////////////////////////////////////////////////////////////////////////

void AliUA1JetFinderV1::RunAlgoritm(Float_t etbgTotal, Double_t dEtTotal, Int_t& nJets,
                                  Float_t* etJet,Float_t* etaJet, Float_t* phiJet,
                                  Float_t* etallJet, Int_t* ncellsJet)
{

   //dump lego
  AliUA1JetHeaderV1* header = (AliUA1JetHeaderV1*) fHeader;
  const Int_t nBinsMax = 120000; // we use a fixed array not to fragment memory
  
  const Int_t nBinEta = header->GetLegoNbinEta();
  const Int_t nBinPhi = header->GetLegoNbinPhi();
  if((nBinPhi*nBinEta)>nBinsMax){
    AliError("Too many bins of the ETA-PHI histogram");
  }

  Float_t etCell[nBinsMax];   //! Cell Energy
  Float_t etaCell[nBinsMax];  //! Cell eta
  Float_t phiCell[nBinsMax];  //! Cell phi
  Short_t   flagCell[nBinsMax]; //! Cell flag

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
  //Float_t etmin   = header->GetMinJetEt();



  // tmp array of jets form algoritm
  Float_t etaAlgoJet[kMaxJets] = {0.0};
  Float_t phiAlgoJet[kMaxJets] = {0.0};
  Float_t etAlgoJet[kMaxJets] = {0.0};
  Int_t   ncellsAlgoJet[kMaxJets] = {0};

  //run algorithm//

  // sort cells by et
  Int_t  index[nBinsMax];
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
  

  for(Int_t icell = 0; icell < nCell; icell++){
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
        for(Int_t kcell =0; kcell < nCell; kcell++){
            Int_t lcell = index[kcell];
            if(lcell == jcell)                 continue; // cell itself
            if(flagCell[lcell] != 0)           continue; // cell used before
            if(etCell[lcell] > etCell[jcell])  continue; // can this happen
            //calculate dr
            deta = etaCell[lcell] - eta;
	    dphi = TMath::Abs(phiCell[lcell] - phi);
	    if (dphi > TMath::Pi()) dphi = 2. * TMath::Pi() - dphi;
	    dr = TMath::Sqrt(deta * deta + dphi * dphi);
            if(dr <= rc){
               // calculate offset from initiate cell
               deta = etaCell[lcell] - eta0;
               dphi = phiCell[lcell] - phi0;
               if (dphi < - TMath::Pi()) dphi=  dphi + 2.0 * TMath::Pi();
	       if (dphi >   TMath::Pi()) dphi = dphi - 2.0 * TMath::Pi();
	       
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
        rc = header->GetRadius();
        for(Int_t ncell =0; ncell < nCell; ncell++){
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

        // select jets with et > background
        // estimate max fluctuation of background in cone
        Double_t ncellin = (Double_t)nCellIn;
        Double_t ntcell  = (Double_t)nCell;
        Double_t etbmax = (etbgTotal + dEtTotal )*(ncellin/ntcell);
        // min cone et
        Double_t etcmin = etCone ;  // could be used etCone - etmin !!
        //desicions !! etbmax < etcmin
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
	   etaAlgoJet[nJets] = eta;
	   phiAlgoJet[nJets] = phi;
	   etAlgoJet[nJets] = etCone;
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
  Int_t idx[kMaxJets];
  TMath::Sort(nJets, etAlgoJet, idx); // sort only the found jets
  for(Int_t p = 0; p < nJets; p++){
     etaJet[p] = etaAlgoJet[idx[p]];
     phiJet[p] = phiAlgoJet[idx[p]];
     etJet[p] = etAlgoJet[idx[p]];
     etallJet[p] = etAlgoJet[idx[p]];
     ncellsJet[p] = ncellsAlgoJet[idx[p]];
  }

}
////////////////////////////////////////////////////////////////////////

void AliUA1JetFinderV1::SubtractBackg(const Int_t& nIn, const Int_t&nJ, Float_t&etbgTotalN,
				      const Float_t* ptT, const Float_t* etaT, const Float_t* phiT,
				      Float_t* etJet, const Float_t* etaJet, const Float_t* phiJet, Float_t* etsigJet,
				      Int_t* multJet, Int_t* injet)
{
  //background subtraction using cone method but without correction in dE/deta distribution

  //calculate energy inside and outside cones
  AliUA1JetHeaderV1* header = (AliUA1JetHeaderV1*) fHeader;
  Float_t rc= header->GetRadius();
  Float_t etIn[kMaxJets] = {0};
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
  Float_t areaJet[kMaxJets];
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

void AliUA1JetFinderV1::SubtractBackgStat(const Int_t& nIn, const Int_t&nJ,Float_t&etbgTotalN,
					  const Float_t* ptT, const Float_t* etaT, const Float_t* phiT,
					  Float_t* etJet, const Float_t* etaJet, const Float_t* phiJet, Float_t* etsigJet,
					  Int_t* multJet, Int_t* injet)
{

  //background subtraction using statistical method
  AliUA1JetHeaderV1* header = (AliUA1JetHeaderV1*) fHeader;
  Float_t etbgStat = header->GetBackgStat(); // pre-calculated background

  //calculate energy inside
  Float_t rc= header->GetRadius();
  Float_t etIn[kMaxJets] = {0.0};

  for(Int_t jpart = 0; jpart < nIn; jpart++){ // loop for all particles in array
     //if((fReader->GetCutFlag(jpart)) != 1) continue; // pt cut
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
                etIn[ijet]+= ptT[jpart];
                if(fReader->GetSignalFlag(jpart) == 1) etsigJet[ijet] += ptT[jpart];
             }
             break;
         }
     }// end jets loop
  } //end particle loop

  //calc jets areas
  Float_t areaJet[kMaxJets];
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
  }

  //subtract background using area method
  for(Int_t ljet=0; ljet<nJ; ljet++){
     Float_t areaRatio = areaJet[ljet]/areaOut;
     etJet[ljet] = etIn[ljet]-etbgStat*areaRatio; // subtraction
  }

  etbgTotalN = etbgStat;

}

////////////////////////////////////////////////////////////////////////

void AliUA1JetFinderV1::SubtractBackgCone(const Int_t& nIn, const Int_t&nJ,Float_t& etbgTotalN,
					  const Float_t* ptT, const Float_t* etaT, const Float_t* phiT,
					  Float_t* etJet, const Float_t* etaJet, const Float_t* phiJet, Float_t* etsigJet,
					  Int_t* multJet, Int_t* injet)
{
   // Cone background subtraction method taking into acount dEt/deta distribution
    AliUA1JetHeaderV1* header = (AliUA1JetHeaderV1*) fHeader;
   //general
   Float_t rc= header->GetRadius();
   Float_t etamax = header->GetLegoEtaMax();
   Float_t etamin = header->GetLegoEtaMin();
   Int_t ndiv = 100;

   // jet energy and area arrays
   Bool_t oldStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);

   for(Int_t mjet=0; mjet<nJ; mjet++){
     if(!fhEtJet[mjet]){ 
       fhEtJet[mjet] = new TH1F(Form("hEtJet%d", mjet),"et dist in eta ",ndiv,etamin,etamax);
     }
     if(!fhAreaJet[mjet]){        
       fhAreaJet[mjet] = new TH1F(Form("hEtJet%d", mjet),"area dist in eta ",ndiv,etamin,etamax);       
     }
     fhEtJet[mjet]->Reset();
     fhAreaJet[mjet]->Reset();
   }
   // background energy and area
   if(!fhEtBackg)fhEtBackg = new TH1F("hEtBackg"," backg et dist in eta ",ndiv,etamin,etamax);
   fhEtBackg->Reset();
   if(!fhAreaBackg) fhAreaBackg = new TH1F("hAreaBackg","backg area dist in eta ",ndiv,etamin,etamax);
   fhAreaBackg->Reset();
   TH1::AddDirectory(oldStatus);

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
             if((fReader->GetCutFlag(jpart)) == 1){// pt cut
                fhEtJet[ijet]->Fill(etaT[jpart],ptT[jpart]); //particle inside cone
                if(fReader->GetSignalFlag(jpart) == 1) etsigJet[ijet] += ptT[jpart];
             }
             break;
         }
     }// end jets loop
     if(injet[jpart] == -1  && fReader->GetCutFlag(jpart) == 1)
        fhEtBackg->Fill(etaT[jpart],ptT[jpart]); // particle outside cones
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
          fhAreaJet[ijet]->Fill(etac,areaj);
          areabg = areabg - areaj;
      } // end jets loop
      fhAreaBackg->Fill(etac,areabg);
      eta0 = eta1;
      eta1 = eta1 + etaw;
   } // end loop for all eta bins

   //subtract background
   for(Int_t kjet=0; kjet<nJ; kjet++){
       etJet[kjet] = 0.0; // first  clear etJet for this jet
       for(Int_t bin = 0; bin< ndiv; bin++){
           if(fhAreaJet[kjet]->GetBinContent(bin)){
              Float_t areab = fhAreaBackg->GetBinContent(bin);
              Float_t etb = fhEtBackg->GetBinContent(bin);
              Float_t areaR = (fhAreaJet[kjet]->GetBinContent(bin))/areab;
              etJet[kjet] = etJet[kjet] + ((fhEtJet[kjet]->GetBinContent(bin)) - etb*areaR); //subtraction
           }
       }
   }

   // calc background total
   Double_t etOut = fhEtBackg->Integral();
   Double_t areaOut = fhAreaBackg->Integral();
   Float_t areaT = 4*(header->GetLegoEtaMax())*TMath::Pi();
   etbgTotalN = etOut*areaT/areaOut;
}

////////////////////////////////////////////////////////////////////////


void AliUA1JetFinderV1::SubtractBackgRatio(const Int_t& nIn, const Int_t&nJ, Float_t& etbgTotalN,
					   const Float_t* ptT, const Float_t* etaT, const Float_t* phiT,
					   Float_t* etJet, const Float_t* etaJet, const Float_t* phiJet, Float_t* etsigJet,
					   Int_t* multJet, Int_t* injet)
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
   // jet energy and area arrays

   Bool_t oldStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);
   for(Int_t mjet=0; mjet<nJ; mjet++){
     if(!fhEtJet[mjet]){ 
       fhEtJet[mjet] = new TH1F(Form("hEtJet%d", mjet),"et dist in eta ",ndiv,etamin,etamax);
     }
     if(!fhAreaJet[mjet]){ 
       fhAreaJet[mjet] = new TH1F(Form("hAreaJet%d", mjet),"area dist in eta ",ndiv,etamin,etamax);       
     }
     fhEtJet[mjet]->Reset();
     fhAreaJet[mjet]->Reset();
   }
   // background energy and area
   if(!fhEtBackg)fhEtBackg = new TH1F("hEtBackg"," backg et dist in eta ",ndiv,etamin,etamax);
   fhEtBackg->Reset();
   if(!fhAreaBackg) fhAreaBackg = new TH1F("hAreaBackg","backg area dist in eta ",ndiv,etamin,etamax);
   fhAreaBackg->Reset();
   TH1::AddDirectory(oldStatus);

   //fill energies
   for(Int_t jpart = 0; jpart < nIn; jpart++){ // loop for all particles in array
     //if((fReader->GetCutFlag(jpart)) != 1) continue;
     for(Int_t ijet=0; ijet<nJ; ijet++){  // loop for all jets
         Float_t deta = etaT[jpart] - etaJet[ijet];
	      Float_t dphi = phiT[jpart] - phiJet[ijet];
         if (dphi < -TMath::Pi()) dphi= -dphi - 2.0 * TMath::Pi();
	      if (dphi > TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
	      Float_t dr = TMath::Sqrt(deta * deta + dphi * dphi);
         if(dr <= rc){ // particles inside this cone
            multJet[ijet]++;
            injet[jpart] = ijet;
            if((fReader->GetCutFlag(jpart)) == 1){ //pt cut
               fhEtJet[ijet]->Fill(etaT[jpart],ptT[jpart]); //particle inside cone and pt cut
               if(fReader->GetSignalFlag(jpart) == 1) etsigJet[ijet] += ptT[jpart];
            }
            break;
         }
     }// end jets loop
     if(injet[jpart] == -1) fhEtBackg->Fill(etaT[jpart],ptT[jpart]); // particle outside cones
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
          fhAreaJet[ijet]->Fill(etac,areaj);
          areabg = areabg - areaj;
      } // end jets loop
      fhAreaBackg->Fill(etac,areabg);
      eta0 = eta1;
      eta1 = eta1 + etaw;
   } // end loop for all eta bins

   //subtract background
   for(Int_t kjet=0; kjet<nJ; kjet++){
       etJet[kjet] = 0.0; // first  clear etJet for this jet
       for(Int_t bin = 0; bin< ndiv; bin++){
           if(fhAreaJet[kjet]->GetBinContent(bin)){
              Float_t areab = fhAreaBackg->GetBinContent(bin);
              Float_t etb = fhEtBackg->GetBinContent(bin);
              Float_t areaR = (fhAreaJet[kjet]->GetBinContent(bin))/areab;
              etJet[kjet] = etJet[kjet] + ((fhEtJet[kjet]->GetBinContent(bin)) - etb*areaR*bgRatioCut); //subtraction
           }
       }
   }

   // calc background total
   Double_t etOut = fhEtBackg->Integral();
   Double_t areaOut = fhAreaBackg->Integral();
   Float_t areaT = 4*(header->GetLegoEtaMax())*TMath::Pi();
   etbgTotalN = etOut*areaT/areaOut;
}

////////////////////////////////////////////////////////////////////////


void AliUA1JetFinderV1::Reset()
{
  fLego->Reset();
  AliJetFinder::Reset();
}

////////////////////////////////////////////////////////////////////////

void AliUA1JetFinderV1::WriteJHeaderToFile() const
{
  AliUA1JetHeaderV1* header = (AliUA1JetHeaderV1*) fHeader;
  header->Write();
}

////////////////////////////////////////////////////////////////////////

void AliUA1JetFinderV1::Init()
{
  // initializes some variables
  AliUA1JetHeaderV1* header = (AliUA1JetHeaderV1*) fHeader;
  fLego = new
    TH2F("legoH","eta-phi",
	 header->GetLegoNbinEta(), header->GetLegoEtaMin(),
	 header->GetLegoEtaMax(),  header->GetLegoNbinPhi(),
	 header->GetLegoPhiMin(),  header->GetLegoPhiMax());
  // Do not store in current dir
  fLego->SetDirectory(0);

}
