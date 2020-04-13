/**************************************************************************
 * Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////
//
// Base class for AOD reconstructed heavy-flavour 2-prong decay
//
// Author: A.Dainese, andrea.dainese@lnl.infn.it
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODMCParticle.h"

/// \cond CLASSIMP
ClassImp(AliAODRecoDecayHF2Prong);
/// \endcond

//--------------------------------------------------------------------------
AliAODRecoDecayHF2Prong::AliAODRecoDecayHF2Prong() :
  AliAODRecoDecayHF()
{
  //
  /// Default Constructor
  //
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF2Prong::AliAODRecoDecayHF2Prong(AliAODVertex *vtx2,
						 Double_t *px,Double_t *py,Double_t *pz,
						 Double_t *d0,Double_t *d0err,Float_t dca) :
  AliAODRecoDecayHF(vtx2,2,0,px,py,pz,d0,d0err)
{
  //
  /// Constructor with AliAODVertex for decay vertex
  //
  SetDCA(dca);
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF2Prong::AliAODRecoDecayHF2Prong(AliAODVertex *vtx2,
						 Double_t *d0,Double_t *d0err,Float_t dca) :
  AliAODRecoDecayHF(vtx2,2,0,d0,d0err)
{
  //
  /// Constructor with AliAODVertex for decay vertex and without prongs momenta
  //
  SetDCA(dca);
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF2Prong::AliAODRecoDecayHF2Prong(const AliAODRecoDecayHF2Prong &source) :
  AliAODRecoDecayHF(source)
{
  //
  /// Copy constructor
  //
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF2Prong &AliAODRecoDecayHF2Prong::operator=(const AliAODRecoDecayHF2Prong &source)
{
  //
  /// assignment operator
  //
  if(&source == this) return *this;

  AliAODRecoDecayHF::operator=(source);

  return *this;
}
//--------------------------------------------------------------------------
Bool_t AliAODRecoDecayHF2Prong::SelectD0(const Double_t *cuts,Int_t &okD0,Int_t &okD0bar) 
  const {
///
/// This function compares the D0 with a set of cuts:
///
/// cuts[0] = inv. mass half width [GeV]
/// cuts[1] = dca [cm]
/// cuts[2] = cosThetaStar
/// cuts[3] = pTK [GeV/c]
/// cuts[4] = pTPi [GeV/c]
/// cuts[5] = d0K [cm]   upper limit!
/// cuts[6] = d0Pi [cm]  upper limit!
/// cuts[7] = d0d0 [cm^2]
/// cuts[8] = cosThetaPoint
///
/// If the D0/D0bar doesn't pass the cuts it sets the weights to 0
/// If neither D0 nor D0bar pass the cuts return kFALSE
///
  Double_t mD0,mD0bar,ctsD0,ctsD0bar;
  okD0=1; okD0bar=1;

  Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();

  if(PtProng(1) < cuts[3] || PtProng(0) < cuts[4]) okD0 = 0;
  if(PtProng(0) < cuts[3] || PtProng(1) < cuts[4]) okD0bar = 0;
  if(!okD0 && !okD0bar) return kFALSE;

  if(TMath::Abs(Getd0Prong(1)) > cuts[5] || 
     TMath::Abs(Getd0Prong(0)) > cuts[6]) okD0 = 0;
  if(TMath::Abs(Getd0Prong(0)) > cuts[6] ||
     TMath::Abs(Getd0Prong(1)) > cuts[5]) okD0bar = 0;
  if(!okD0 && !okD0bar) return kFALSE;

  if(GetDCA() > cuts[1]) { okD0 = okD0bar = 0; return kFALSE; }

  InvMassD0(mD0,mD0bar);
  if(TMath::Abs(mD0-mD0PDG)    > cuts[0]) okD0 = 0;
  if(TMath::Abs(mD0bar-mD0PDG) > cuts[0]) okD0bar = 0;
  if(!okD0 && !okD0bar) return kFALSE;

  CosThetaStarD0(ctsD0,ctsD0bar);
  if(TMath::Abs(ctsD0)    > cuts[2]) okD0 = 0;
  if(TMath::Abs(ctsD0bar) > cuts[2]) okD0bar = 0;
  if(!okD0 && !okD0bar) return kFALSE;

  if(Prodd0d0() > cuts[7]) { okD0 = okD0bar = 0; return kFALSE; }

  if(CosPointingAngle()   < cuts[8]) { okD0 = okD0bar = 0; return kFALSE; }

  return kTRUE;
}
//-----------------------------------------------------------------------------
Bool_t AliAODRecoDecayHF2Prong::SelectBtoJPSI(const Double_t *cuts,Int_t &okB)
  const {
///
/// This function compares the Secondary JPSI candidates with a set of cuts:
///
/// cuts[0] = inv. mass half width [GeV]
/// cuts[1] = dca [cm]
/// cuts[2] = cosThetaStar (negative electron)
/// cuts[3] = pTP [GeV/c]
/// cuts[4] = pTN [GeV/c]
/// cuts[5] = d0P [cm]   upper limit!
/// cuts[6] = d0N [cm]  upper limit!
/// cuts[7] = d0d0 [cm^2]
/// cuts[8] = cosThetaPoint
///
/// If the candidate doesn't pass the cuts it sets the weight to 0
/// and return kFALSE
///
  Double_t mJPsi,ctsJPsi;
  okB=1;

  Double_t mJPSIPDG = TDatabasePDG::Instance()->GetParticle(443)->Mass();

  if(PtProng(1) < cuts[3] || PtProng(0) < cuts[4]) okB = 0;
  if(!okB) return kFALSE;

  if(TMath::Abs(Getd0Prong(1)) > cuts[5] ||
     TMath::Abs(Getd0Prong(0)) > cuts[6]) okB = 0;
  if(!okB) return kFALSE;

  if(GetDCA() > cuts[1]) { okB = 0; return kFALSE; }

  mJPsi=InvMassJPSIee();
  if(TMath::Abs(mJPsi-mJPSIPDG)    > cuts[0]) okB = 0;
  if(!okB) return kFALSE;

  ctsJPsi=CosThetaStarJPSI();
  if(TMath::Abs(ctsJPsi)    > cuts[2]) okB = 0;
  if(!okB) return kFALSE;

  if(Prodd0d0() > cuts[7]) { okB = 0; return kFALSE; }

  if(CosPointingAngle()   < cuts[8]) { okB = 0; return kFALSE; }

  return kTRUE;
}
//-----------------------------------------------------------------------------
Int_t AliAODRecoDecayHF2Prong::MatchToMCB2Prong(Int_t pdgabs,Int_t pdgabs2prong,Int_t *pdgDg,Int_t *pdgDg2prong,TClonesArray *mcArray) const
{
    //
    // Check if this candidate is matched to a MC signal
    // If no, return -1
    // If yes, return label (>=0) of the AliAODMCParticle
    // If yes, but rejected by momentum conservation check, return -1*label (to allow checks at task level)
    // 
    // NB: This function is only for non-resonant decays (for both 2 prongs)
    // NB: Loosened cut on mom. conserv. (needed because of small issue in ITS Upgrade productions)

    Int_t signlabMother = 1;
  
    // Check number of daughters. Candidate is AliAODRecoDecayHF2Prong, so only continue when 2 daughters
    Int_t ndg = GetNDaughters();
    if (!ndg) { AliError("HF2Prong: No daughters available"); return -1;}
    if (ndg != 2) { AliError(Form("HF2Prong: %d daughters instead of 2",ndg)); return -1;}
    
    // loop on daughters and write the labels
    Int_t dgLabels[2] = {-1};
    if(pdgabs2prong == 0){
        for(Int_t i=0; i<ndg; i++) {
            AliAODTrack *trk = (AliAODTrack*)GetDaughter(i);
            dgLabels[i] = trk->GetLabel();
        }
    } else {
        AliAODTrack *trk = (AliAODTrack*)GetDaughter(0);
        dgLabels[0] = trk->GetLabel();
        AliAODRecoDecayHF2Prong* daug2prong = (AliAODRecoDecayHF2Prong*)GetDaughter(1);
        //Daughter prong (independent on their decay time) will also be far from PV, so call this function again.
        Int_t pdgempty[2]={0,0};
        dgLabels[1] = daug2prong->MatchToMCB2Prong(pdgabs2prong, 0, pdgDg2prong, pdgempty, mcArray);
        if (dgLabels[1] < -1) signlabMother = -1;
    }
    if (dgLabels[0] == -1) return -1;
    if (dgLabels[1] == -1) return -1;
    
    Int_t labMom[2] = {0, 0};
    Int_t i, j, lab, labMother, pdgMother, pdgPart;
    AliAODMCParticle *part = 0;
    AliAODMCParticle *mother = 0;
    Double_t pxSumDgs = 0., pySumDgs = 0., pzSumDgs = 0.;
    Bool_t pdgUsed[2] = {kFALSE, kFALSE};
    
    // loop on daughter labels
    for (i = 0; i < ndg; i++){
        labMom[i] = -1;
        lab = TMath::Abs(dgLabels[i]);
        if (lab < 0){
            printf("daughter with negative label %d\n", lab);
            return -1;
        }
        part = (AliAODMCParticle*)mcArray->At(lab);
        if (!part){
            printf("no MC particle\n");
            return -1;
        }
        
        // check the PDG of the daughter, if requested
        pdgPart = TMath::Abs(part->GetPdgCode());
        for (j = 0; j < ndg; j++){
            if (!pdgUsed[j] && pdgPart == pdgDg[j]){
                pdgUsed[j] = kTRUE;
                break;
            }
        }
        
        mother = part;
        while (mother->GetMother() >= 0){
            labMother = mother->GetMother();
            mother = (AliAODMCParticle*)mcArray->At(labMother);
            if (!mother){
                printf("no MC mother particle\n");
                break;
            }
            pdgMother = TMath::Abs(mother->GetPdgCode());
            if (pdgMother == pdgabs){
                labMom[i] = labMother;
                // keep sum of daughters' momenta, to check for mom conservation
                pxSumDgs += part->Px();
                pySumDgs += part->Py();
                pzSumDgs += part->Pz();
                break;
            } else break;
        }
        if (labMom[i] == -1) return -1; // mother PDG not ok for this daughter
    } // end loop on daughters
    
    // check if the candidate is signal
    labMother = labMom[0];
    // all labels have to be the same and !=-1
    for (i = 0; i < ndg; i++){
        if (labMom[i] == -1)        return -1;
        if (labMom[i] != labMother) return -1;
    }
    
    // check that all daughter PDGs are matched
    for (i = 0; i < ndg; i++){
        if (pdgUsed[i] == kFALSE) return -1;
    }
    
    // check the number of daughters (we are not looking at resonant decay)
    if(mother->GetNDaughters() != 2) return -1;
  
    // Check for mom conservation
    mother = (AliAODMCParticle*)mcArray->At(labMother);
    Double_t pxMother = mother->Px();
    Double_t pyMother = mother->Py();
    Double_t pzMother = mother->Pz();
    // within 2% (fix for ITSUpgrade MC's: 0.00001 -> 0.02)
    if ((TMath::Abs(pxMother - pxSumDgs) / (TMath::Abs(pxMother) + 1.e-13)) > 0.02 &&
        (TMath::Abs(pyMother - pySumDgs) / (TMath::Abs(pyMother) + 1.e-13)) > 0.02 &&
        (TMath::Abs(pzMother - pzSumDgs) / (TMath::Abs(pzMother) + 1.e-13)) > 0.02)
    {
        AliWarning(Form("Mom. cons. not within 2.0%% perc for decay pdgabs = %d daughters = %d. Returning negative label!",pdgabs,(Int_t)mother->GetNDaughters()));
        signlabMother = -1;
    }
    labMother = signlabMother * labMother;
    return labMother;
}
//-----------------------------------------------------------------------------
Int_t AliAODRecoDecayHF2Prong::MatchToMCB3Prong(Int_t pdgabs, Int_t pdgabs3prong, Int_t *pdgBDg, Int_t *pdgDg3prong, TClonesArray *mcArray) const
{
  //
  // Check if this candidate is matched to a MC signal
  // If no, return -1
  // If yes, return label (>=0) of the AliAODMCParticle
  // If yes, but rejected by momentum conservation check, return -1*label (to allow checks at task level)
  //
  // NB: This function is only for non-resonant decays of the 2prong (3prong can decay resonant)
  // NB: Loosened cut on mom. conserv. (needed because of small issue in ITS Upgrade productions)
  //
  
  Int_t signlabMother = 1;
  
  // Check number of daughters. Candidate is AliAODRecoDecayHF2Prong, so only continue when 2 daughters
  Int_t ndg = GetNDaughters();
  if (!ndg) { AliWarning("HF2Prong: No daughters available"); return -1;}
  if (ndg != 2) { AliWarning(Form("HF2Prong: %d daughters instead of 2",ndg)); return -1;}
  
  // loop on daughters and write the labels
  Int_t dgLabels[2] = {-1};
  
  AliAODTrack *trk = (AliAODTrack*)GetDaughter(1);
  if(!trk) return -1;
  dgLabels[1] = trk->GetLabel();
  
  AliAODRecoDecayHF3Prong* daug3prong = (AliAODRecoDecayHF3Prong*)GetDaughter(0);
  if(!daug3prong) return -1;
  //Similar as daug3prong->MatchToMC(pdgabs3prong,mcArray,3,pdgDg3prong)
  //but w/o mom conservation check (fix for ITSUpgrade MC's)
  Int_t ndg3pr = daug3prong->GetNDaughters();
  if (!ndg3pr) { AliWarning("HF3Prong: No daughters available"); return -1;}
  if (ndg3pr != 3) { AliWarning(Form("HF3Prong: %d daughters instead of 3",ndg3pr)); return -1;}
  
  Int_t dgLabels3pr[3] = {-1};
  for(Int_t i=0; i<ndg3pr; i++) {
    AliAODTrack *trk = (AliAODTrack*)daug3prong->GetDaughter(i);
    dgLabels3pr[i] = trk->GetLabel();
  }
  
  Int_t labMom3pr[3]={0,0,0};
  Int_t lab3pr, labMother3pr, pdgMother3pr, pdgPart3pr;
  AliAODMCParticle *part3pr=0;
  AliAODMCParticle *mother3pr=0;
  Double_t pxSumDgs3pr=0., pySumDgs3pr=0., pzSumDgs3pr=0.;
  Bool_t pdgUsed3pr[3]={kFALSE,kFALSE,kFALSE};
  
  // loop on daughter labels
  for(int i=0; i<ndg3pr; i++) {
    labMom3pr[i]=-1;
    
    lab3pr = TMath::Abs(dgLabels3pr[i]);
    if(lab3pr<0) { AliWarning(Form("daughter with negative label %d",lab3pr)); return -1; }
    
    part3pr = (AliAODMCParticle*)mcArray->At(lab3pr);
    if(!part3pr) { AliWarning("no MC particle"); return -1; }
    
    // check the PDG of the daughter
    pdgPart3pr=TMath::Abs(part3pr->GetPdgCode());
    for(int j=0; j<ndg3pr; j++) {
      if(!pdgUsed3pr[j] && pdgPart3pr==pdgDg3prong[j]) {
        pdgUsed3pr[j]=kTRUE;
        break;
      }
    }
    
    mother3pr = part3pr;
    while(mother3pr->GetMother()>=0) {
      labMother3pr=mother3pr->GetMother();
      mother3pr = (AliAODMCParticle*)mcArray->At(labMother3pr);
      if(!mother3pr) { AliWarning("no MC mother particle"); break; }
      
      pdgMother3pr = TMath::Abs(mother3pr->GetPdgCode());
      if(pdgMother3pr==pdgabs3prong) {
        labMom3pr[i]=labMother3pr;
        // keep sum of daughters' momenta, to check for mom conservation
        pxSumDgs3pr += part3pr->Px();
        pySumDgs3pr += part3pr->Py();
        pzSumDgs3pr += part3pr->Pz();
        break;
      } else if(pdgMother3pr>pdgabs3prong || pdgMother3pr<10) {
        break;
      }
    }
    if(labMom3pr[i]==-1) return -1; // mother PDG not ok for this daughter
  } // end loop on daughters
  
  // check if the candidate is signal
  labMother3pr=labMom3pr[0];
  // all labels have to be the same and !=-1
  for(int i=0; i<ndg3pr; i++) {
    if(labMom3pr[i]==-1)           return -1;
    if(labMom3pr[i]!=labMother3pr) return -1;
  }
  
  // check that all daughter PDGs are matched
  for (int i = 0; i < ndg3pr; i++){
    if (pdgUsed3pr[i] == kFALSE) return -1;
  }
  
  mother3pr = (AliAODMCParticle*)mcArray->At(labMother3pr);
  Double_t pxMother3pr = mother3pr->Px();
  Double_t pyMother3pr = mother3pr->Py();
  Double_t pzMother3pr = mother3pr->Pz();
  // within 2% (fix for ITSUpgrade MC's: 0.00001 -> 0.02)
  if((TMath::Abs(pxMother3pr-pxSumDgs3pr)/(TMath::Abs(pxMother3pr)+1.e-13)) > 0.02 &&
     (TMath::Abs(pyMother3pr-pySumDgs3pr)/(TMath::Abs(pyMother3pr)+1.e-13)) > 0.02 &&
     (TMath::Abs(pzMother3pr-pzSumDgs3pr)/(TMath::Abs(pzMother3pr)+1.e-13)) > 0.02) {
    AliWarning(Form("Mom. cons. not within 2.0%% perc for decay pdgabs = %d daughters = %d. Returning negative label!",pdgabs3prong,(Int_t)mother3pr->GetNDaughters()));
    signlabMother = -1;
  }
  dgLabels[0] = labMother3pr;
  
  Int_t labMom[2] = {0, 0};
  Int_t lab, labMother, pdgMother, pdgPart;
  AliAODMCParticle *part = 0;
  AliAODMCParticle *mother = 0;
  Double_t pxSumDgs = 0., pySumDgs = 0., pzSumDgs = 0.;
  Bool_t pdgUsed[2] = {kFALSE, kFALSE};
  
  // loop on daughter labels
  for (int i = 0; i < ndg; i++){
    labMom[i] = -1;
    
    lab = TMath::Abs(dgLabels[i]);
    if(lab<0) { AliWarning(Form("daughter with negative label %d",lab)); return -1; }
    
    part = (AliAODMCParticle*)mcArray->At(lab);
    if(!part) { AliWarning("no MC particle"); return -1; }
    
    
    // check the PDG of the daughter, if requested
    pdgPart = TMath::Abs(part->GetPdgCode());
    for (int j = 0; j < ndg; j++){
      if (!pdgUsed[j] && pdgPart == pdgBDg[j]){
        pdgUsed[j] = kTRUE;
        break;
      }
    }
    
    mother = part;
    while (mother->GetMother() >= 0){
      labMother = mother->GetMother();
      mother = (AliAODMCParticle*)mcArray->At(labMother);
      if(!mother) { AliWarning("no MC mother particle"); break; }
      
      pdgMother = TMath::Abs(mother->GetPdgCode());
      if (pdgMother == pdgabs){
        labMom[i] = labMother;
        // keep sum of daughters' momenta, to check for mom conservation
        pxSumDgs += part->Px();
        pySumDgs += part->Py();
        pzSumDgs += part->Pz();
        break;
      } else break;
    }
    if (labMom[i] == -1) return -1; // mother PDG not ok for this daughter
  } // end loop on daughters
  
  // check if the candidate is signal
  labMother = labMom[0];
  // all labels have to be the same and !=-1
  for (int i = 0; i < ndg; i++){
    if (labMom[i] == -1)        return -1;
    if (labMom[i] != labMother) return -1;
  }
  
  // check that all daughter PDGs are matched
  for (int i = 0; i < ndg; i++){
    if (pdgUsed[i] == kFALSE) return -1;
  }
  
  // check the number of daughters (we are not looking at resonant decay)
  if(mother->GetNDaughters() != 2) return -1;
  
  // Check for mom conservation
  mother = (AliAODMCParticle*)mcArray->At(labMother);
  Double_t pxMother = mother->Px();
  Double_t pyMother = mother->Py();
  Double_t pzMother = mother->Pz();
  // within 2% (fix for ITSUpgrade MC's: 0.00001 -> 0.02)
  if((TMath::Abs(pxMother-pxSumDgs)/(TMath::Abs(pxMother)+1.e-13)) > 0.02 &&
     (TMath::Abs(pyMother-pySumDgs)/(TMath::Abs(pyMother)+1.e-13)) > 0.02 &&
     (TMath::Abs(pzMother-pzSumDgs)/(TMath::Abs(pzMother)+1.e-13)) > 0.02) {
    AliWarning(Form("Mom. cons. not within 2.0%% perc for decay pdgabs = %d daughters = %d. Returning negative label!",pdgabs,(Int_t)mother->GetNDaughters()));
    signlabMother = -1;
  }
  labMother = signlabMother * labMother;
  
  return labMother;
}
