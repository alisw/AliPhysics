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
// FastJet finder
// interface to FastJet algorithm
// Author: Rafael.Diaz.Valdes@cern.ch
// kt using NlnN 
//---------------------------------------------------------------------

#include <Riostream.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TArrayF.h>
#include <TRandom.h>
#include <TClonesArray.h>

#include "AliFastJetFinder.h"
#include "AliFastJetHeader.h"
#include "AliJetReaderHeader.h"
#include "AliJetReader.h"
#include "AliJet.h"
//for FastJet finder
#include "FjPseudoJet.hh"
#include "FjClusterSequence.hh"


ClassImp(AliFastJetFinder);


////////////////////////////////////////////////////////////////////////

AliFastJetFinder::AliFastJetFinder():
    AliJetFinder(),
    fHeader(0x0),
    fLego(0x0),
    fLegoSignal(0x0)
{
  // Constructor
}

////////////////////////////////////////////////////////////////////////

AliFastJetFinder::~AliFastJetFinder()

{
  // destructor
}

////////////////////////////////////////////////////////////////////////


void AliFastJetFinder::FindJets()

{
  //create cells and jets array
  // 1) transform input to pt,eta,phi plus lego
  TClonesArray *lvArray = fReader->GetMomentumArray();
  Int_t nIn =  lvArray->GetEntries();
  if (nIn == 0) return;

  // local arrays for particles input
  Float_t* enT  = new Float_t[nIn];
  Float_t* ptT  = new Float_t[nIn];
  Float_t* etaT = new Float_t[nIn];
  Float_t* phiT = new Float_t[nIn];
  Int_t*   injet = new Int_t[nIn];


  //total energy in array
  Float_t  EtTotal = 0.0;
  Float_t  meanptCell = 0.0;
  Float_t  sqptCell = 0.0;


  // load input vectors in fLego
  for (Int_t i = 0; i < nIn; i++){
    TLorentzVector *lv = (TLorentzVector*) lvArray->At(i);
    enT[i]  = lv->Energy();
    ptT[i]  = lv->Pt();
    etaT[i] = lv->Eta();
    phiT[i] = ((lv->Phi() < 0) ? (lv->Phi()) + 2 * TMath::Pi() : lv->Phi());
    if (fReader->GetCutFlag(i) == 1){
      fLego->Fill(etaT[i], phiT[i], ptT[i]);
      if(fReader->GetSignalFlag(i) == 1)
        fLegoSignal->Fill(etaT[i], phiT[i], ptT[i]);
      EtTotal= EtTotal+ptT[i];
    }
  }
  fJets->SetNinput(nIn);

  // add soft background fixed
  Int_t nsoft = (fHeader->GetLegoNbinEta())*(fHeader->GetLegoNbinPhi());
  Float_t* ptRndm =  new Float_t[nsoft];
  if(fHeader->AddSoftBackg()){
    gRandom->RndmArray(nsoft,ptRndm);
    for(Int_t isoft = 0; isoft < nsoft; isoft++){
        Float_t ptsoft  = 0.005*ptRndm[isoft];
        EtTotal= EtTotal+ptsoft;
    }
  }

  if(EtTotal == 0){
     delete enT;
     delete ptT;
     delete etaT;
     delete phiT;
     return;
  }

  //cell array
  Float_t* etCell = new Float_t[90000];   //! Cell Energy // check enough space! *to be done*
  Float_t* etaCell = new Float_t[90000];  //! Cell eta
  Float_t* phiCell = new Float_t[90000];  //! Cell phi
  Int_t*   jetflagCell = new Int_t[90000]; //! Cell flag for jets
  Float_t* etsigCell = new Float_t[90000];   // signal in this cell

  //jets array
  Float_t* etaJet = new Float_t[200];
  Float_t* phiJet = new Float_t[200];
  Float_t* etJet  = new Float_t[200];
  Float_t* etsigJet  = new Float_t[200]; //signal energy
  Float_t* etallJet = new Float_t[200];  //total energy in jet area
  Int_t*   ncellsJet = new Int_t[200];
  memset(etaJet,0,sizeof(Float_t)*200);
  memset(phiJet,0,sizeof(Float_t)*200);
  memset(etJet,0,sizeof(Float_t)*200);
  memset(etsigJet,0,sizeof(Float_t)*200);
  memset(etallJet,0,sizeof(Float_t)*200);
  memset(ncellsJet,0,sizeof(Int_t)*200);



  // load cells arrays
  Int_t nCell = 0;
  TAxis* xaxis = fLego->GetXaxis();
  TAxis* yaxis = fLego->GetYaxis();
  Float_t e = 0.0; Float_t esig = 0.0;

  for (Int_t ib = 1; ib <= fHeader->GetLegoNbinEta(); ib++) {
      for (Int_t jb = 1; jb <= fHeader->GetLegoNbinPhi(); jb++) {
	       e = fLego->GetBinContent(ib,jb);
	       if (e < 0.0) continue; // don't include this cells
	       Float_t eta  = xaxis->GetBinCenter(ib);
	       Float_t phi  = yaxis->GetBinCenter(jb);
          if(fHeader->AddSoftBackg())
             etCell[nCell]  = e + 0.005*ptRndm[nCell];
          else
             etCell[nCell]  = e;
          sqptCell = sqptCell + etCell[nCell]*etCell[nCell]; // xi^2 ////////
	       etaCell[nCell] = eta;
	       phiCell[nCell] = phi;
          jetflagCell[nCell] = -1; //default
          esig = fLegoSignal->GetBinContent(ib,jb);
          if(esig > 0.0)
             etsigCell[nCell] = esig;
          else
             etsigCell[nCell] = 0.0;
	       nCell++;
      }
  }

  meanptCell = EtTotal/(Float_t)nCell;
  sqptCell = sqptCell/(Float_t)nCell;

  Int_t nJets = 0;
  //call to FastJet Algorithm
  RunAlgorithm(nJets,etJet,etaJet,phiJet,etsigJet,etallJet,ncellsJet,
               nCell,etCell,etaCell,phiCell,etsigCell,jetflagCell);


  //subtract background
  SubtractBackg(nCell,jetflagCell,etCell,
                nJets,etJet,etallJet,ncellsJet,
                meanptCell,sqptCell,EtTotal);


  // add jets to list
  Int_t* index = new Int_t[nJets];
  Int_t nj = 0;
  for(Int_t kj=0; kj<nJets; kj++){
      if ((etaJet[kj] > (fHeader->GetJetEtaMax())) ||
          (etaJet[kj] < (fHeader->GetJetEtaMin())) ||
          (etJet[kj] < fHeader->GetMinJetEt())) continue; // acceptance eta range and etmin
      Float_t px, py,pz,en; // convert to 4-vector
      px = etJet[kj] * TMath::Cos(phiJet[kj]);
      py = etJet[kj] * TMath::Sin(phiJet[kj]);
      pz = etJet[kj] / TMath::Tan(2.0 * TMath::ATan(TMath::Exp(-etaJet[kj])));
      en = TMath::Sqrt(px * px + py * py + pz * pz);
      fJets->AddJet(px, py, pz, en);
      index[nj] = kj;
      nj++;
  }

  //add signal percentage and total signal  in AliJets for analysis tool
  Float_t* percentage  = new Float_t[nj];
  Int_t* ncells      = new Int_t[nj];
  Int_t* mult        = new Int_t[nj];

  for(Int_t i = 0; i< nj; i++){
     percentage[i] = etsigJet[index[i]]/etJet[index[i]];
     ncells[i] = ncellsJet[index[i]];
  }

   //reorder injet flags
  for(Int_t ipar = 0; ipar < nIn; ipar++){
     Float_t injetflag =0;
     Int_t iparCell = fLego->FindBin(etaT[ipar], phiT[ipar]);
     injet[ipar] = jetflagCell[iparCell];
     for(Int_t js = 0; js < nj; js++){
        if(injet[ipar] == index[js]){
          injet[ipar] = js;  // set the new jet id value
          mult[js]++; // add multiplicity in jet js
          injetflag = 1;
          break;
        }
     }
     if(injetflag == 0) injet[ipar] = -1; // set default value
  }


  fJets->SetNCells(ncells);
  fJets->SetPtFromSignal(percentage);
  fJets->SetMultiplicities(mult);
  fJets->SetInJet(injet);
  fJets->SetEtaIn(etaT);
  fJets->SetPhiIn(phiT);
  fJets->SetPtIn(ptT);
  fJets->SetEtAvg(meanptCell);

   //delete
  delete enT;
  delete ptT;
  delete etaT;
  delete phiT;
  delete injet;
  //cells
  delete etCell;
  delete etaCell;
  delete phiCell;
  delete jetflagCell;
  delete etsigCell;
  //jets
  delete etaJet;
  delete phiJet;
  delete etJet;
  delete etsigJet;
  delete etallJet;
  delete ncellsJet;

  delete index;
  delete percentage;
  delete ncells;
  delete mult;
  delete ptRndm;

}

////////////////////////////////////////////////////////////////////////
void AliFastJetFinder::RunAlgorithm(Int_t& nJets,Float_t* etJet,Float_t* etaJet,Float_t* phiJet,
                                    Float_t* etsigJet, Float_t* etallJet, Int_t* ncellsJet,
                                    Int_t& nCell,Float_t* etCell,Float_t* etaCell,Float_t* phiCell,
                                    Float_t* etsigCell, Int_t* jetflagCell)
{
   //FastJet objects
   vector<FjPseudoJet> input_cells; // create a vector
   for (Int_t i = 0; i < nCell; i++){
      if(etCell[i] == 0.0) continue; // not include cell empty
      Double_t px, py,pz,en; // convert to 4-vector
      px = etCell[i]*TMath::Cos(phiCell[i]);
      py = etCell[i]*TMath::Sin(phiCell[i]);
      pz = etCell[i]/TMath::Tan(2.0 * TMath::ATan(TMath::Exp(-etaCell[i])));
      en = TMath::Sqrt(px * px + py * py + pz * pz);
      FjPseudoJet input_cell(px,py,pz,en); // create FjPseudoJet object
      input_cell.set_user_index(i); //label the cell into Fastjet algortihm
      //push FjPseudoJet of (px,py,pz,en) onto back of the input_cells
      input_cells.push_back(input_cell);
   }

   //run the jet clustering with option R=1.0 and strategy= Best
   Double_t Rparam = fHeader->GetRadius(); // default 1.0;
   FjClusterSequence clust_seq(input_cells,Rparam);


   //vector to get clusters
   vector<FjPseudoJet> clusters;

   ///////////////////////////////////////////////////////////////////////////
   //extract the inclusive jets with pt> ptmin sorted by pt
   //Float_t areaT = 4*(fHeader->GetLegoEtaMax())*TMath::Pi();
  // Double_t ptbgRc = EtBackgT*(Rparam*Rparam*TMath::Pi()/areaT);
  // Double_t ptbgRcfluct = dEtTotal*Rparam*TMath::Sqrt(TMath::Pi()/areaT);
  // Double_t ptmin = ptbgRc + ptbgRcfluct;
   clusters = clust_seq.inclusive_jets(0);
   //////////////////////////////////////////////////////////////////////////

   /////////////////////////////////////////////////////////////////////////
   //extract the exclusive jets with dcut = 25 GeV**2 and sort them in order
   //of increasing pt
   //Float_t areaT = 4*(fHeader->GetLegoEtaMax())*TMath::Pi();
   //Double_t ptbgRc = EtBackgT*(Rparam*Rparam*TMath::Pi()/areaT);
   //Double_t ptbgRcfluct = dEtTotal*Rparam*TMath::Sqrt(TMath::Pi()/areaT);
   //Double_t ptmin = ptbgRc + ptbgRcfluct;
   //Double_t ktbackg = (fHeader->GetKfactor())*ptmin;
   //Double_t dcut = ktbackg*ktbackg;
   //clusters = sorted_by_pt(clust_seq.exclusive_jets(dcut));
   //clusters = sorted_by_pt(clust_seq.exclusive_jets(5));
   /////////////////////////////////////////////////////////////////////////


   //cout << " ktClusters " << clusters.size() << endl;
   nJets = (Int_t)clusters.size();
   ////////////////////////////////////////////////////////////////////////
   // add all clusters to jet arrays
   for(UInt_t ij = 0; ij < clusters.size(); ij++){
       //constituents
       vector<FjPseudoJet> constituents = clust_seq.constituents(clusters[ij]);
       //fill jet array info
       ncellsJet[ij] = (Int_t)constituents.size();
       phiJet[ij] = clusters[ij].phi();
       Float_t angle = TMath::ATan(clusters[ij].perp()/clusters[ij].pz());
       angle = ((angle < 0) ? angle + TMath::Pi() : angle);
       etaJet[ij] = - TMath::Log(TMath::Tan(angle/2.0));
       etJet[ij] = clusters[ij].perp();
       //get constituents cells
       for(UInt_t jc = 0; jc < constituents.size(); jc++){ // loop for all cells in ij cluster
           Int_t jcell = constituents[jc].user_index();
           jetflagCell[jcell] = ij; //flag this cell for jet
           etsigJet[ij] = etsigJet[ij] + etsigCell[jcell]; // add signal of this cell
           etallJet[ij] = etallJet[ij] + etCell[jcell];   // add total of this cell
       }
   }

}
////////////////////////////////////////////////////////////////////////
void AliFastJetFinder::SubtractBackg(Int_t& nCell, Int_t* jetflagCell, Float_t* etCell,
                                     Int_t& nJets, Float_t* etJet, Float_t* etallJet, Int_t* ncellsJet,
                                     Float_t& meanptCell, Float_t& sqptCell, Float_t& etBackg)
{
   // simplest method: subtract background from external region to jets

   //tmp array to flag jets
   Int_t flagJet[200];
   //tmp array to flag jet-cell status
   Int_t tmpjetflagCell[90000];

   Float_t etBackgOld = 0;
   Float_t prec  = fHeader->GetPrecBg();
   Float_t bgprec = 1;

   while(bgprec > prec){
        //clear tmpjetflagCell
        memset(tmpjetflagCell,-1,sizeof(Int_t)*90000); // init with -1 (all cells are background)
        //clear flagjet
        memset(flagJet,0,sizeof(Int_t)*200); // init with 0 (no flag jets)
        // select clusters > meantmpCell
        for(Int_t i = 0; i < nJets; i++){
            Float_t iptcell = etallJet[i]/(Float_t)ncellsJet[i];
            if(iptcell < meanptCell) continue; // cluster not selected
            // convert tmp cell background to jet cell
            for(Int_t ic = 0; ic < nCell; ic++){ //loop for all cells
                if(jetflagCell[ic] != i) continue; // other cells
                tmpjetflagCell[ic] = i; // convert to a jet cell
            }
            //load total energy in cluster
            etJet[i] = etallJet[i];
            flagJet[i] = 1; // flag jet
       }
       //subtract background
       for(Int_t j = 0; j < nCell; j++){ // loop for all cells
           Int_t idxjet = tmpjetflagCell[j];
           if(idxjet == -1) continue; // background cell
           if(idxjet == -2) continue; // background protected cell
           etJet[idxjet] = etJet[idxjet] - meanptCell;
       }
       // evaluate background fluctuations (rms value)
       Float_t rmsptCell = TMath::Sqrt(sqptCell - meanptCell*meanptCell);
       //fake jets
       for(Int_t k = 0; k < nJets; k++){
           if(flagJet[k] != 1) continue; // only flaged jets
           //if(etJet[k] > fHeader->GetMinJetEt()) continue;  // jet ok!!
           if(etJet[k] > rmsptCell*ncellsJet[k]) continue;  // jet ok!!
           //clear tmpjetflag in cells
           for(Int_t kc = 0; kc < nCell; kc++){ //loop for all cells
               if(tmpjetflagCell[kc] != k) continue; // other cells
               tmpjetflagCell[kc] = -1; // convert to background tmp cell
           }
           // clear all previous jet flags
           etJet[k] = 0;
           flagJet[k] = 0;
       }
       // recalculate background
       sqptCell = 0;
       etBackgOld = etBackg;
       etBackg = 0;
       Int_t nCellBackg = 0;
       for(Int_t l = 0; l < nCell; l++){ // loop for all cells
          if(tmpjetflagCell[l] != -1) continue; // cell included in some jet or protected
          nCellBackg++;
          etBackg = etBackg + etCell[l]; //add cell to background
          //calc sqptCell
          sqptCell = sqptCell + etCell[l]*etCell[l];
       }
       if(nCellBackg){
          meanptCell = etBackg/(Float_t)nCellBackg; // new pt cell mean value
          sqptCell = sqptCell/(Float_t)nCellBackg;
       }else{
          meanptCell = 0;
          sqptCell = 0;
       }
       // evaluate presicion values
       if(etBackg)
         bgprec = (etBackgOld - etBackg)/etBackg;
       else
         bgprec = 0;
   }

   // set etJet 0 for all clusters not flaged in order to
   for(Int_t m = 0; m < nJets; m++){
       if(flagJet[m] == 1) continue; // flaged jets
       etJet[m] = 0; //others clusters
   }


}

////////////////////////////////////////////////////////////////////////
void AliFastJetFinder::SubtractBackgArea(Int_t& nCell, Int_t* jetflagCell,
                                     Float_t* etCell,Int_t&nJets, Float_t* etJet, Float_t* etallJet)
{
   // area method: subtract background from external region to jets
   // using fixed area pi*Rc2

   // n cells contained in a cone Rc
   Double_t Rc = fHeader->GetRadius();
   Float_t nCellRc = Rc*Rc*TMath::Pi()/(0.015*0.015); // change in future !!!!
   //tmp array to flag fake jets
   Int_t fakeflagJet[100];
   memset(fakeflagJet,0,sizeof(Int_t)*100); // init with 0 (no fake jets)
   Int_t njfake = nJets;
   while(njfake > 0){
       //calc background per cell
       Int_t nCellBackg = 0;
       Float_t EtBackg = 0.0;
       for(Int_t i = 0; i < nCell; i++){ // loop for all cells
           if(jetflagCell[i] != -1) continue; // cell included in some jet
           nCellBackg++;
           EtBackg = EtBackg + etCell[i]; //add cell to background
       }
       //subtract background energy per jet
       for(Int_t l = 0; l < nJets; l++){
           if(fakeflagJet[l] == 1) continue; // fake jet
           etJet[l] = etallJet[l] - nCellRc*EtBackg/(Float_t)nCellBackg;
       }
       //fake jets analysis
       njfake = 0;
       for(Int_t k = 0; k < nJets; k++){
           if(fakeflagJet[k] == 1) continue; // fake jet
           if(etJet[k] < fHeader->GetMinJetEt()){  // a new fake jet
              //clear jet flag in cells
              for(Int_t kc = 0; kc < nCell; kc++){ //loop for all cells
                  Int_t kidx = jetflagCell[kc];
                  if(kidx != k) continue; // other cells
                  jetflagCell[kc] = -1; // convert to background cell
              }
              fakeflagJet[k] = 1; // mark as a fake jet
              njfake++; //count fakes in this loop
           }
       }
   }
}

////////////////////////////////////////////////////////////////////////



void AliFastJetFinder::Reset()
{
  fLego->Reset();
  fLegoSignal->Reset();
  fJets->ClearJets();

}
////////////////////////////////////////////////////////////////////////

void AliFastJetFinder::WriteJHeaderToFile()
{
  fOut->cd();
  fHeader->Write();
}

////////////////////////////////////////////////////////////////////////

void AliFastJetFinder::Init()
{
  // initializes some variables
   // book lego
  fLego = new
    TH2F("legoH","eta-phi",
	 fHeader->GetLegoNbinEta(), fHeader->GetLegoEtaMin(),
	 fHeader->GetLegoEtaMax(),  fHeader->GetLegoNbinPhi(),
	 fHeader->GetLegoPhiMin(), fHeader->GetLegoPhiMax());
  fLegoSignal = new
    TH2F("legoSignalH","eta-phi signal",
	 fHeader->GetLegoNbinEta(), fHeader->GetLegoEtaMin(),
	 fHeader->GetLegoEtaMax(),  fHeader->GetLegoNbinPhi(),
	 fHeader->GetLegoPhiMin(), fHeader->GetLegoPhiMax());

}
