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
// UA1 Cone Algorithm Jet finder
// manages the search for jets
// Author: Rafael.Diaz.Valdes@cern.ch
// (version in c++)
//---------------------------------------------------------------------

#include <TLorentzVector.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TArrayF.h>
#include "AliUA1JetFinderV1.h"
#include "AliUA1JetHeaderV1.h"
#include "AliJetReaderHeader.h"
#include "AliJetReader.h"
#include "AliJet.h"


ClassImp(AliUA1JetFinderV1)

////////////////////////////////////////////////////////////////////////

AliUA1JetFinderV1::AliUA1JetFinderV1()

{
  // Constructor
  fHeader = 0x0;
  fLego   = 0x0;
}

////////////////////////////////////////////////////////////////////////

AliUA1JetFinderV1::~AliUA1JetFinderV1()

{
  // destructor
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
  TClonesArray *lvArray = fReader->GetMomentumArray();
  Int_t nIn =  lvArray->GetEntries();
  if (nIn == 0) return;

  // local arrays for input
  Float_t* ptT  = new Float_t[nIn];
  Float_t* etaT = new Float_t[nIn];
  Float_t* phiT = new Float_t[nIn];
  Int_t*   injet = new Int_t[nIn];

  //total energy in array
  Float_t  etbgTotal = 0.0;
  TH1F* hPtTotal = new TH1F("hPt","Pt distribution of all particles ",100,0.0,15.0);

  // load input vectors and calculate total energy in array
  for (Int_t i = 0; i < nIn; i++){
    TLorentzVector *lv = (TLorentzVector*) lvArray->At(i);
    ptT[i]  = lv->Pt();
    etaT[i] = lv->Eta();
    phiT[i] = ((lv->Phi() < 0) ? (lv->Phi()) + 2 * TMath::Pi() : lv->Phi());
    if (fReader->GetCutFlag(i) != 1) continue;
    fLego->Fill(etaT[i], phiT[i], ptT[i]);
    hPtTotal->Fill(ptT[i]);
    etbgTotal+= ptT[i];
  }
  fJets->SetNinput(nIn);

  // calculate total energy and fluctuation in map
  Double_t meanpt = hPtTotal->GetMean();
  Double_t ptRMS = hPtTotal->GetRMS();
  Double_t npart = hPtTotal->GetEntries();
  Double_t dEtTotal = (TMath::Sqrt(npart))*TMath::Sqrt(meanpt * meanpt + ptRMS*ptRMS);

  // arrays to hold jets
  Float_t* etaJet = new Float_t[30];
  Float_t* phiJet = new Float_t[30];
  Float_t* etJet  = new Float_t[30];
  Float_t* etsigJet  = new Float_t[30]; //signal et in jet
  Float_t* etallJet = new Float_t[30];  // total et in jet (tmp variable)
  Int_t* ncellsJet = new Int_t[30];
  Int_t* multJet  = new Int_t[30];
  Int_t nJets; // to hold number of jets found by algorithm
  Int_t nj;    // number of jets accepted
  Float_t prec  = fHeader->GetPrecBg();
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
     nJets = 0;
     nj = 0;
     // reset particles-jet array in memory
     memset(injet,-1,sizeof(Int_t)*nIn);
     //run cone algorithm finder
     RunAlgoritm(etbgTotal,dEtTotal,nJets,etJet,etaJet,phiJet,etallJet,ncellsJet);
     //run background subtraction
     if(nJets > fHeader->GetNAcceptJets()) // limited number of accepted jets per event
       nj = fHeader->GetNAcceptJets();
     else
       nj = nJets;
     //subtract background
     Float_t etbgTotalN = 0.0; //new background
     if(fHeader->GetBackgMode() == 1) // standar
        SubtractBackg(nIn,nj,etbgTotalN,ptT,etaT,phiT,etJet,etaJet,phiJet,etsigJet,multJet,injet);
     if(fHeader->GetBackgMode() == 2) //cone
        SubtractBackgCone(nIn,nj,etbgTotalN,ptT,etaT,phiT,etJet,etaJet,phiJet,etsigJet,multJet,injet);
     if(fHeader->GetBackgMode() == 3) //ratio
        SubtractBackgRatio(nIn,nj,etbgTotalN,ptT,etaT,phiT,etJet,etaJet,phiJet,etsigJet,multJet,injet);
     if(fHeader->GetBackgMode() == 4) //statistic
        SubtractBackgStat(nIn,nj,etbgTotalN,ptT,etaT,phiT,etJet,etaJet,phiJet,etsigJet,multJet,injet);
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
  for(Int_t kj=0; kj<nj; kj++){
     if ((etaJet[kj] > (fHeader->GetJetEtaMax())) ||
          (etaJet[kj] < (fHeader->GetJetEtaMin())) ||
          (etJet[kj] < fHeader->GetMinJetEt())) continue; // acceptance eta range and etmin
      Float_t px, py,pz,en; // convert to 4-vector
      px = etJet[kj] * TMath::Cos(phiJet[kj]);
      py = etJet[kj] * TMath::Sin(phiJet[kj]);
      pz = etJet[kj] / TMath::Tan(2.0 * TMath::ATan(TMath::Exp(-etaJet[kj])));
      en = TMath::Sqrt(px * px + py * py + pz * pz);
      fJets->AddJet(px, py, pz, en);
      idxjets[nselectj] = kj;
      nselectj++;
  }
  //add signal percentage and total signal  in AliJets for analysis tool
  Float_t* percentage  = new Float_t[nselectj];
  Int_t* ncells      = new Int_t[nselectj];
  Int_t* mult        = new Int_t[nselectj];
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
  fJets->SetNCells(ncells);
  fJets->SetPtFromSignal(percentage);
  fJets->SetMultiplicities(mult);
  fJets->SetInJet(injet);
  fJets->SetEtaIn(etaT);
  fJets->SetPhiIn(phiT);
  fJets->SetPtIn(ptT);
  fJets->SetEtAvg(etbgTotal/(4*(fHeader->GetLegoEtaMax())*TMath::Pi()));


  //delete
  delete ptT;
  delete etaT;
  delete phiT;
  delete injet;
  delete hPtTotal;
  delete etaJet;
  delete phiJet;
  delete etJet;
  delete etsigJet;
  delete etallJet;
  delete ncellsJet;
  delete multJet;
  delete idxjets;
  delete percentage;
  delete ncells;
  delete mult;


}

////////////////////////////////////////////////////////////////////////

void AliUA1JetFinderV1::RunAlgoritm(Float_t etbgTotal, Double_t dEtTotal, Int_t& nJets,
                                  Float_t* etJet,Float_t* etaJet, Float_t* phiJet,
                                  Float_t* etallJet, Int_t* ncellsJet)
{

   //dump lego
  // check enough space! *to be done*
  Float_t etCell[60000];   //! Cell Energy
  Float_t etaCell[60000];  //! Cell eta
  Float_t phiCell[60000];  //! Cell phi
  Int_t   flagCell[60000]; //! Cell flag

  Int_t nCell = 0;
  TAxis* xaxis = fLego->GetXaxis();
  TAxis* yaxis = fLego->GetYaxis();
  Float_t e = 0.0;
  for (Int_t i = 1; i <= fHeader->GetLegoNbinEta(); i++) {
      for (Int_t j = 1; j <= fHeader->GetLegoNbinPhi(); j++) {
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
  Float_t minmove = fHeader->GetMinMove();
  Float_t maxmove = fHeader->GetMaxMove();
  Float_t rc      = fHeader->GetRadius();
  Float_t etseed  = fHeader->GetEtSeed();
  //Float_t etmin   = fHeader->GetMinJetEt();



  // tmp array of jets form algoritm
  Float_t etaAlgoJet[30];
  Float_t phiAlgoJet[30];
  Float_t etAlgoJet[30];
  Int_t   ncellsAlgoJet[30];

  //run algorithm//

  // sort cells by et
  Int_t * index  = new Int_t[nCell];
  TMath::Sort(nCell, etCell, index);
  // variable used in centroide loop
  Float_t eta = 0.0;
  Float_t phi = 0.0;
  Float_t eta0 = 0.0;
  Float_t phi0 = 0.0;
  Float_t etab = 0.0;
  Float_t phib = 0.0;
  Float_t etas = 0.0;
  Float_t phis = 0.0;
  Float_t ets = 0.0;
  Float_t deta = 0.0;
  Float_t dphi = 0.0;
  Float_t dr = 0.0;
  Float_t etsb = 0.0;
  Float_t etasb = 0.0;
  Float_t phisb = 0.0;


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
            if(lcell == jcell) continue; // cell itself
            if(flagCell[lcell] != 0) continue; // cell used before
            if(etCell[lcell] > etCell[jcell]) continue;
            //calculate dr
            deta = etaCell[lcell] - eta;
	         dphi = phiCell[lcell] - phi;
	         if (dphi < -TMath::Pi()) dphi= -dphi - 2.0 * TMath::Pi();
	         if (dphi > TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
	         dr = TMath::Sqrt(deta * deta + dphi * dphi);
            if(dr <= rc){
               // calculate offset from initiate cell
               deta = etaCell[lcell] - eta0;
               dphi = phiCell[lcell] - phi0;
               if (dphi < -TMath::Pi()) dphi= -dphi - 2.0 * TMath::Pi();
	            if (dphi > TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
               etas = etas + etCell[lcell]*deta;
               phis = phis + etCell[lcell]*dphi;
               ets = ets + etCell[lcell];
               //new weighted eta and phi including this cell
               eta = eta0 + etas/ets;
               phi = phi0 + phis/ets;
               // if cone does not move much, just go to next step
               dr = TMath::Sqrt((eta-etab)*(eta-etab) + (phi-phib)*(phi-phib));
               if(dr <= minmove) break;
               // cone should not move more than max_mov
               dr = TMath::Sqrt((etas/ets)*(etas/ets) + (phis/ets)*(phis/ets));
               if(dr > maxmove){
                  eta = etab;
                  phi = phib;
                  ets = etsb;
                  etas = etasb;
                  phis = phisb;
               }else{ // store this loop information
                 etab=eta;
                 phib=phi;
                 etsb = ets;
                 etasb = etas;
                 phisb = phis;
               }
            }
        }//end of cells loop looking centroide

        //avoid cones overloap (to be implemented in the future)

        //flag cells in Rc, estimate total energy in cone
        Float_t etCone = 0.0;
        Int_t   nCellIn = 0;
        rc = fHeader->GetRadius();
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
             etaAlgoJet[nJets] = eta;
             phiAlgoJet[nJets] = phi;
             etAlgoJet[nJets] = etCone;
             ncellsAlgoJet[nJets] = nCellIn;
             nJets++;
        }

  } // end of cells loop

  //reorder jets by et in cone
  //sort jets by energy
  Int_t * idx  = new Int_t[nJets];
  TMath::Sort(nJets, etAlgoJet, idx);
  for(Int_t p = 0; p < nJets; p++){
     etaJet[p] = etaAlgoJet[idx[p]];
     phiJet[p] = phiAlgoJet[idx[p]];
     etJet[p] = etAlgoJet[idx[p]];
     etallJet[p] = etAlgoJet[idx[p]];
     ncellsJet[p] = ncellsAlgoJet[idx[p]];
  }


  //delete
  delete index;
  delete idx;

}
////////////////////////////////////////////////////////////////////////

void AliUA1JetFinderV1::SubtractBackg(Int_t& nIn, Int_t&nJ, Float_t&etbgTotalN,
                      Float_t* ptT, Float_t* etaT, Float_t* phiT,
                      Float_t* etJet,Float_t* etaJet, Float_t* phiJet, Float_t* etsigJet,
                      Int_t* multJet, Int_t* injet)
{
  //background subtraction using cone method but without correction in dE/deta distribution

  //calculate energy inside and outside cones
  Float_t rc= fHeader->GetRadius();
  Float_t etIn[30];
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
     if(injet[jpart] == -1 && fReader->GetSignalFlag(jpart) == 1)
        etOut += ptT[jpart]; // particle outside cones and pt cut
  } //end particle loop

  //estimate jets and background areas
  Float_t areaJet[30];
  Float_t areaOut = 4*(fHeader->GetLegoEtaMax())*TMath::Pi();
  for(Int_t k=0; k<nJ; k++){
      Float_t detamax = etaJet[k] + rc;
      Float_t detamin = etaJet[k] - rc;
      Float_t accmax = 0.0; Float_t accmin = 0.0;
      if(detamax > fHeader->GetLegoEtaMax()){ // sector outside etamax
         Float_t h = fHeader->GetLegoEtaMax() - etaJet[k];
         accmax = rc*rc*TMath::ACos(h/rc) - h*TMath::Sqrt(rc*rc - h*h);
      }
      if(detamin < fHeader->GetLegoEtaMin()){ // sector outside etamin
         Float_t h = fHeader->GetLegoEtaMax() + etaJet[k];
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
  Float_t areaT = 4*(fHeader->GetLegoEtaMax())*TMath::Pi();
  etbgTotalN = etOut*areaT/areaOut;


}

////////////////////////////////////////////////////////////////////////

void AliUA1JetFinderV1::SubtractBackgStat(Int_t& nIn, Int_t&nJ,Float_t&etbgTotalN,
                      Float_t* ptT, Float_t* etaT, Float_t* phiT,
                      Float_t* etJet,Float_t* etaJet, Float_t* phiJet, Float_t* etsigJet,
                      Int_t* multJet, Int_t* injet)
{

  //background subtraction using statistical method

  Float_t etbgStat = fHeader->GetBackgStat(); // pre-calculated background

  //calculate energy inside
  Float_t rc= fHeader->GetRadius();
  Float_t etIn[30];

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
  Float_t areaJet[30];
  Float_t areaOut = 4*(fHeader->GetLegoEtaMax())*TMath::Pi();
  for(Int_t k=0; k<nJ; k++){
      Float_t detamax = etaJet[k] + rc;
      Float_t detamin = etaJet[k] - rc;
      Float_t accmax = 0.0; Float_t accmin = 0.0;
      if(detamax > fHeader->GetLegoEtaMax()){ // sector outside etamax
         Float_t h = fHeader->GetLegoEtaMax() - etaJet[k];
         accmax = rc*rc*TMath::ACos(h/rc) - h*TMath::Sqrt(rc*rc - h*h);
      }
      if(detamin < fHeader->GetLegoEtaMin()){ // sector outside etamin
         Float_t h = fHeader->GetLegoEtaMax() + etaJet[k];
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

void AliUA1JetFinderV1::SubtractBackgCone(Int_t& nIn, Int_t&nJ,Float_t& etbgTotalN,
                      Float_t* ptT, Float_t* etaT, Float_t* phiT,
                      Float_t* etJet,Float_t* etaJet, Float_t* phiJet, Float_t* etsigJet,
                      Int_t* multJet, Int_t* injet)
{
   // Cone background subtraction method taking into acount dEt/deta distribution

   //general
   Float_t rc= fHeader->GetRadius();
   Float_t etamax = fHeader->GetLegoEtaMax();
   Float_t etamin = fHeader->GetLegoEtaMin();
   Int_t ndiv = 100;

   // jet energy and area arrays
   TH1F* hEtJet[30];
   TH1F* hAreaJet[30];
   for(Int_t mjet=0; mjet<nJ; mjet++){
     char hEtname[256]; char hAreaname[256];
     sprintf(hEtname, "hEtJet%d", mjet); sprintf(hAreaname, "hAreaJet%d", mjet);
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
             if((fReader->GetCutFlag(jpart)) == 1){// pt cut
                hEtJet[ijet]->Fill(etaT[jpart],ptT[jpart]); //particle inside cone
                if(fReader->GetSignalFlag(jpart) == 1) etsigJet[ijet] += ptT[jpart];
             }
             break;
         }
     }// end jets loop
     if(injet[jpart] == -1  && fReader->GetSignalFlag(jpart) == 1)
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
   Float_t areaT = 4*(fHeader->GetLegoEtaMax())*TMath::Pi();
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


void AliUA1JetFinderV1::SubtractBackgRatio(Int_t& nIn, Int_t&nJ,Float_t& etbgTotalN,
                      Float_t* ptT, Float_t* etaT, Float_t* phiT,
                      Float_t* etJet,Float_t* etaJet, Float_t* phiJet, Float_t* etsigJet,
                       Int_t* multJet, Int_t* injet)
{
   // Ratio background subtraction method taking into acount dEt/deta distribution

   //factor F calc before
    Float_t bgRatioCut = fHeader->GetBackgCutRatio();


   //general
   Float_t rc= fHeader->GetRadius();
   Float_t etamax = fHeader->GetLegoEtaMax();
   Float_t etamin = fHeader->GetLegoEtaMin();
   Int_t ndiv = 100;

   // jet energy and area arrays
   TH1F* hEtJet[30];
   TH1F* hAreaJet[30];
   for(Int_t mjet=0; mjet<nJ; mjet++){
     char hEtname[256]; char hAreaname[256];
     sprintf(hEtname, "hEtJet%d", mjet); sprintf(hAreaname, "hAreaJet%d", mjet);
     hEtJet[mjet] = new TH1F(hEtname,"et dist in eta ",ndiv,etamin,etamax);        // change range
     hAreaJet[mjet] = new TH1F(hAreaname,"area dist in eta ",ndiv,etamin,etamax);  // change range
  }
   // background energy and area
   TH1F* hEtBackg = new TH1F("hEtBackg"," backg et dist in eta ",ndiv,etamin,etamax);         // change range
   TH1F* hAreaBackg = new TH1F("hAreaBackg","backg area dist in eta ",ndiv,etamin,etamax);  // change range

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
               hEtJet[ijet]->Fill(etaT[jpart],ptT[jpart]); //particle inside cone and pt cut
               if(fReader->GetSignalFlag(jpart) == 1) etsigJet[ijet] += ptT[jpart];
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
   Float_t areaT = 4*(fHeader->GetLegoEtaMax())*TMath::Pi();
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


void AliUA1JetFinderV1::Reset()
{
  fLego->Reset();
  fJets->ClearJets();
}

////////////////////////////////////////////////////////////////////////

void AliUA1JetFinderV1::WriteJHeaderToFile()
{
  fOut->cd();
  fHeader->Write();
}

////////////////////////////////////////////////////////////////////////

void AliUA1JetFinderV1::Init()
{
  // initializes some variables

   // book lego
  fLego = new
    TH2F("legoH","eta-phi",
	 fHeader->GetLegoNbinEta(), fHeader->GetLegoEtaMin(),
	 fHeader->GetLegoEtaMax(),  fHeader->GetLegoNbinPhi(),
	 fHeader->GetLegoPhiMin(), fHeader->GetLegoPhiMax());

}
