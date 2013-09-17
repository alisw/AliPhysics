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

//--------------------------------------------------
// Method implementation for background studies and background subtraction with UA1 algorithms
//
// Author: magali.estienne@subatech.in2p3.fr
//-------------------------------------------------

#include <Riostream.h> 
#include <TList.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>

#include "AliAODJetEventBackground.h"
#include "AliUA1JetHeaderV1.h"
#include "AliJetCalTrk.h"
#include "AliJetBkg.h"


using namespace std;

ClassImp(AliJetBkg)

////////////////////////////////////////////////////////////////////////

AliJetBkg::AliJetBkg():
  TObject(),
  fEvent(0x0),
  fHeader(0x0),
  fDebug(0),
  fhEtBackg(0x0),
  fhAreaBackg(0x0)
{
  // Default constructor
  for(int i = 0;i < kMaxJets;i++){
    fhAreaJet[i] = fhEtJet[i] = 0;
  }
}

//----------------------------------------------------------------
AliJetBkg::AliJetBkg(const AliJetBkg& input):
  TObject(input),
  fEvent(input.fEvent),
  fHeader(input.fHeader),
  fDebug(input.fDebug),
  fhEtBackg(input.fhEtBackg),
  fhAreaBackg(input.fhAreaBackg)
{
  // copy constructor
  for(int i = 0;i < kMaxJets;i++){
    fhAreaJet[i] = input.fhAreaJet[i];
    fhEtJet[i] = input.fhEtJet[i];
  }

}

//----------------------------------------------------------------
AliJetBkg::~AliJetBkg()
{
  // Destructor
  if(fhEtBackg) delete  fhEtBackg;
  if(fhAreaBackg) delete  fhAreaBackg;
   for(int i = 0;i < kMaxJets;i++){
     if(fhAreaJet[i]) delete fhAreaJet[i];
     if(fhEtJet[i])  delete fhEtJet[i];
   }

}

//----------------------------------------------------------------
Bool_t AliJetBkg::PtCutPass(Int_t id, Int_t nTracks)
{
  // Check if track or cell passes the cut flag
  if(id < nTracks && fEvent->GetCalTrkTrack(id)->GetCutFlag() == 1) 
   return kTRUE;
  else return kFALSE;

}

//----------------------------------------------------------------
Bool_t AliJetBkg::SignalCutPass(Int_t id, Int_t nTracks)
{
  // Check if track or cell passes the cut flag
  if(id < nTracks && fEvent->GetCalTrkTrack(id)->GetSignalFlag() == 1)
    return kTRUE;
  else return kFALSE;

}

//----------------------------------------------------------------
Float_t AliJetBkg::CalcJetAreaEtaCut(Float_t radius, const Float_t etaJet)
{
  // Calculate jet area taking into account an acceptance cut in eta
  AliUA1JetHeader* header = (AliUA1JetHeader*) fHeader;
  Float_t detamax = etaJet + radius;
  Float_t detamin = etaJet - radius;
  Float_t accmax = 0.0; Float_t accmin = 0.0;
  if(detamax > header->GetLegoEtaMax()){ // sector outside etamax
    Float_t h = header->GetLegoEtaMax() - etaJet;
    accmax = radius*radius*TMath::ACos(h/radius) - h*TMath::Sqrt(radius*radius - h*h);
  }
  if(detamin < header->GetLegoEtaMin()){ // sector outside etamin
    Float_t h = header->GetLegoEtaMax() + etaJet;
    accmin = radius*radius*TMath::ACos(h/radius) - h*TMath::Sqrt(radius*radius - h*h);
  }
  
  return radius*radius*TMath::Pi() - accmax - accmin;

}

//----------------------------------------------------------------
void AliJetBkg::CalcJetAndBckgAreaEtaCut(Bool_t calcOutsideArea, Float_t radius, const Int_t nJets, const Float_t* etaJet, Float_t* &areaJet, Float_t &areaOut)
{
  // Calculate jet and bacground areas taking into account an acceptance cut in eta

  AliUA1JetHeader* header = (AliUA1JetHeader*) fHeader;
  areaOut = (header->GetLegoEtaMax()-header->GetLegoEtaMin())*(header->GetLegoPhiMax() - header->GetLegoPhiMin());
  for(Int_t k=0; k<nJets; k++){
    areaJet[k] = CalcJetAreaEtaCut(radius, etaJet[k]);
    if(calcOutsideArea) areaOut = areaOut - areaJet[k];
  }

}

//----------------------------------------------------------------
void AliJetBkg::SubtractBackg(const Int_t& nIn, const Int_t&nJ, Float_t&etbgTotalN, Float_t&sigmaN,
                              const Float_t* ptT, const Float_t* etaT, const Float_t* phiT,
                              Float_t* etJet, const Float_t* etaJet, const Float_t* phiJet,
                              Float_t* etsigJet, Int_t* multJetT, Int_t* multJetC, Int_t* multJet,
                              Int_t* injet, Float_t* &areaJet)
{
  //
  // Background subtraction using cone method but without correction in dE/deta distribution
  // Cases to take into account the EMCal geometry are included
  //
  
  AliUA1JetHeader* header = (AliUA1JetHeader*) fHeader;
  //calculate energy inside and outside cones
  fDebug = header->GetDebug();
  Float_t rc = header->GetRadius();
  Float_t etOut = 0;
  // Get number of tracks from EventCalTrk
  Int_t nTracks = fEvent->GetNCalTrkTracks();

  Float_t etIn[kMaxJets] = {0};
  Float_t areaOut = 0.;

  for(Int_t jpart = 0; jpart < nIn; jpart++){ // loop for all particles in array

    for(Int_t ijet=0; ijet<nJ; ijet++){
      Float_t deta = etaT[jpart] - etaJet[ijet];
      Float_t dphi = phiT[jpart] - phiJet[ijet];
      if (dphi < -TMath::Pi()) dphi= -dphi - 2.0 * TMath::Pi();
      if (dphi > TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
      Float_t dr = TMath::Sqrt(deta * deta + dphi * dphi);
      if(dr <= rc){ // particles inside this cone
        if(jpart < nTracks) multJetT[ijet]++;
        else multJetC[ijet]++;
        multJet[ijet]++;
        injet[jpart] = ijet;
        if(PtCutPass(jpart,nTracks)){ // pt cut
          etIn[ijet] += ptT[jpart];
          if(SignalCutPass(jpart,nTracks))
            etsigJet[ijet]+= ptT[jpart];
        }
        break;
      }
    }// end jets loop

    if((injet[jpart] == -1) &&
       (PtCutPass(jpart,nTracks))){
      etOut += ptT[jpart]; // particle outside cones and pt cut
    }
  } //end particle loop

  // Calculate jet and background areas
  Bool_t calcAreaOut = kTRUE;
  CalcJetAndBckgAreaEtaCut(calcAreaOut,rc, nJ, etaJet, areaJet, areaOut);

  //subtract background using area method
  for(Int_t ljet=0; ljet<nJ; ljet++){
    Float_t areaRatio = areaJet[ljet]/areaOut;
    etJet[ljet] = etIn[ljet]-etOut*areaRatio; // subtraction
  }

  // estimate new total background
  Float_t areaT = 0;
  areaT = (header->GetLegoEtaMax()-header->GetLegoEtaMin())*(header->GetLegoPhiMax()-header->GetLegoPhiMin());
  etbgTotalN = etOut*areaT/areaOut;

  // estimate standard deviation of background
  Int_t count = 0;
  for(Int_t jpart = 0; jpart < nIn; jpart++){ // loop for all particles in array
    if((injet[jpart] == -1) &&
       (PtCutPass(jpart,nTracks))){
      sigmaN += etbgTotalN/areaT - ptT[jpart];
      // To be checked (Division by jet area to obtain standard deviation of rho ?)

      count=count+1;
    }
  }
  if (count>0)
    sigmaN=TMath::Sqrt(TMath::Abs(sigmaN)/count);

}

//----------------------------------------------------------------
void AliJetBkg::SubtractBackgStat(const Int_t& nIn, const Int_t&nJ,Float_t&etbgTotalN, Float_t&sigmaN,
                                  const Float_t* ptT, const Float_t* etaT, const Float_t* phiT,
                                  Float_t* etJet, const Float_t* etaJet, const Float_t* phiJet,
                                  Float_t* etsigJet, Int_t* multJetT, Int_t* multJetC, Int_t* multJet,
                                  Int_t* injet, Float_t* &areaJet)
{
  //
  //background subtraction using statistical method
  // Cases to take into account the EMCal geometry are included
  //

  AliUA1JetHeader* header = (AliUA1JetHeader*) fHeader;
  Float_t etbgStat = header->GetBackgStat(); // pre-calculated background
  
  //calculate energy inside
  Float_t rc= header->GetRadius();
  Float_t etIn[kMaxJets] = {0.0};
  // Get number of tracks from EventCalTrk
  Int_t nTracks = fEvent->GetNCalTrkTracks();
  Float_t areaOut = 0.;

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
            if(jpart < nTracks) multJetT[ijet]++;
            else multJetC[ijet]++;
            multJet[ijet]++;
            injet[jpart] = ijet;

            if(PtCutPass(jpart,nTracks)){ // pt cut
              etIn[ijet] += ptT[jpart];
	      if(SignalCutPass(jpart,nTracks))
		etsigJet[ijet]+= ptT[jpart];
            }
            break;
          }
	}// end jets loop
    } //end particle loop
  
  // Calculate jet and background areas
  Bool_t calcAreaOut = kFALSE;
  CalcJetAndBckgAreaEtaCut(calcAreaOut,rc, nJ, etaJet, areaJet, areaOut);

  //subtract background using area method
  for(Int_t ljet=0; ljet<nJ; ljet++){
    Float_t areaRatio = areaJet[ljet]/areaOut;
    etJet[ljet] = etIn[ljet]-etbgStat*areaRatio; // subtraction
  }
  Int_t count=0;
  etbgTotalN = etbgStat;

  // estimate standard deviation of background
  Float_t areaT = 0;
  areaT = (header->GetLegoEtaMax()-header->GetLegoEtaMin())*(header->GetLegoPhiMax()-header->GetLegoPhiMin());
  for(Int_t jpart = 0; jpart < nIn; jpart++){ // loop for all particles in array
    if((injet[jpart] == -1) &&
       (PtCutPass(jpart,nTracks))){
      sigmaN += etbgTotalN/areaT - ptT[jpart];
      count=count+1;
    }
  }
  if(count>0)sigmaN=TMath::Sqrt(TMath::Abs(sigmaN)/count);

}

//----------------------------------------------------------------
void AliJetBkg::SubtractBackgCone(const Int_t& nIn, const Int_t&nJ,Float_t& etbgTotalN, Float_t&sigmaN,
                                  const Float_t* ptT, const Float_t* etaT, const Float_t* phiT, Float_t* etJet,
                                  const Float_t* etaJet, const Float_t* phiJet, Float_t* etsigJet,
                                  Int_t* multJetT, Int_t* multJetC, Int_t* multJet, Int_t* injet, Float_t* &/*areaJet*/)
{
  //
  // Cone background subtraction method taking into acount dEt/deta distribution
  // Cases to take into account the EMCal geometry are not included
  //

  AliUA1JetHeader* header = (AliUA1JetHeader*) fHeader;
  //general
  Float_t rc= header->GetRadius();
  Float_t etamax = header->GetLegoEtaMax();
  Float_t etamin = header->GetLegoEtaMin();
  Int_t ndiv = 100;
  // Get number of tracks from EventCalTrk
  Int_t nTracks = fEvent->GetNCalTrkTracks();
 
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
        if(jpart < nTracks) multJetT[ijet]++;
        else multJetC[ijet]++;
        multJet[ijet]++;
        injet[jpart] = ijet;

        if(PtCutPass(jpart,nTracks)){ // pt cut
          fhEtJet[ijet]->Fill(etaT[jpart],ptT[jpart]); //particle inside cone
          if(SignalCutPass(jpart,nTracks))
            etsigJet[ijet]+= ptT[jpart];
        }
        break;
      }
    }// end jets loop

    if((injet[jpart] == -1)  &&
       (PtCutPass(jpart,nTracks) == 1))
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
  Float_t areaT = (header->GetLegoEtaMax()-header->GetLegoEtaMin())*(header->GetLegoPhiMax()-header->GetLegoPhiMin());
  etbgTotalN = etOut*areaT/areaOut;

  Int_t count=0;
  for(Int_t jpart = 0; jpart < nIn; jpart++){ // loop for all particles in array
    if((injet[jpart] == -1) &&
       (PtCutPass(jpart,nTracks))){
      sigmaN += etbgTotalN/areaT - ptT[jpart];
      count=count+1;
    }
  }
  sigmaN=TMath::Sqrt(TMath::Abs(sigmaN)/count);
  
}

//----------------------------------------------------------------
void AliJetBkg::SubtractBackgRatio(const Int_t& nIn, const Int_t&nJ,Float_t& etbgTotalN, Float_t&sigmaN,
                                   const Float_t* ptT,const Float_t* etaT, const Float_t* phiT,
                                   Float_t* etJet, const Float_t* etaJet, const Float_t* phiJet,
                                   Float_t* etsigJet, Int_t* multJetT, Int_t* multJetC, Int_t* multJet,
                                   Int_t* injet,  Float_t* &/*areaJet*/)
{
  // Ratio background subtraction method taking into acount dEt/deta distribution
  // Cases to take into account the EMCal geometry are not included

  AliUA1JetHeader* header = (AliUA1JetHeader*) fHeader;
  //factor F calc before
  Float_t bgRatioCut = header->GetBackgCutRatio();
  
  //general
  Float_t rc= header->GetRadius();
  Float_t etamax = header->GetLegoEtaMax();
  Float_t etamin = header->GetLegoEtaMin();
  Int_t ndiv = 100;
  // Get number of tracks from EventCalTrk
  Int_t nTracks = fEvent->GetNCalTrkTracks();
  
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
    for(Int_t ijet=0; ijet<nJ; ijet++){  // loop for all jets
      Float_t deta = etaT[jpart] - etaJet[ijet];
      Float_t dphi = phiT[jpart] - phiJet[ijet];
      if (dphi < -TMath::Pi()) dphi= -dphi - 2.0 * TMath::Pi();
      if (dphi > TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
      Float_t dr = TMath::Sqrt(deta * deta + dphi * dphi);
      if(dr <= rc){ // particles inside this cone
        if(jpart < nTracks) multJetT[ijet]++;
        else multJetC[ijet]++;
        multJet[ijet]++;
        injet[jpart] = ijet;

        if(PtCutPass(jpart,nTracks)){ // pt cut
          fhEtJet[ijet]->Fill(etaT[jpart],ptT[jpart]); //particle inside cone
	  if(SignalCutPass(jpart,nTracks))
	    etsigJet[ijet]+= ptT[jpart];
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
  Float_t areaT = (header->GetLegoEtaMax()-header->GetLegoEtaMin())*(header->GetLegoPhiMax()-header->GetLegoPhiMin());
  etbgTotalN = etOut*areaT/areaOut;
 
  Int_t count=0;
    
  for(Int_t jpart = 0; jpart < nIn; jpart++){ // loop for all particles in array
    if((injet[jpart] == -1) &&
       (PtCutPass(jpart,nTracks))){
      sigmaN += etbgTotalN/areaT - ptT[jpart];
      count=count+1;
    }
  }
  sigmaN=TMath::Sqrt(TMath::Abs(sigmaN)/count);

}

