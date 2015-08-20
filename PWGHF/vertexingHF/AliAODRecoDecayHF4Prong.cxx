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
// Base class for AOD reconstructed heavy-flavour 4-prong decay
//
// Authors: G.E.Bruno Giuseppe.Bruno@to.infn.it, R.Romita Rossella.Romita@ba.infn.it
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF4Prong.h"

ClassImp(AliAODRecoDecayHF4Prong)

//--------------------------------------------------------------------------
AliAODRecoDecayHF4Prong::AliAODRecoDecayHF4Prong() :
  AliAODRecoDecayHF(), 
  fDist12toPrim(0),
  fDist3toPrim(0),
  fDist4toPrim(0)
{
  //
  // Default Constructor
  //
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF4Prong::AliAODRecoDecayHF4Prong(AliAODVertex *vtx2,
						 Double_t *px,Double_t *py,Double_t *pz,
						 Double_t *d0,Double_t *d0err,
						 Double_t *dca, //Double_t sigvert,
						 Double_t dist12,Double_t dist3,
						 Double_t dist4,
						 Short_t charge) :
  AliAODRecoDecayHF(vtx2,4,charge,px,py,pz,d0,d0err),
  // fSigmaVert(sigvert),
  fDist12toPrim(dist12),
  fDist3toPrim(dist3),
  fDist4toPrim(dist4)
{
  //
  // Constructor with AliAODVertex for decay vertex
  //
  SetDCAs(6,dca);
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF4Prong::AliAODRecoDecayHF4Prong(AliAODVertex *vtx2,
						 Double_t *d0,Double_t *d0err,
						 Double_t *dca, //Double_t sigvert,
						 Double_t dist12,Double_t dist3, 
						 Double_t dist4, 
						 Short_t charge) :
  AliAODRecoDecayHF(vtx2,4,charge,d0,d0err),
  //fSigmaVert(sigvert),
  fDist12toPrim(dist12),
  fDist3toPrim(dist3),
  fDist4toPrim(dist4)
{
  //
  // Constructor with AliAODVertex for decay vertex and without prongs momenta
  //
  SetDCAs(6,dca);
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF4Prong::AliAODRecoDecayHF4Prong(const AliAODRecoDecayHF4Prong &source) :
  AliAODRecoDecayHF(source),
  //fSigmaVert(source.fSigmaVert),
  fDist12toPrim(source.fDist12toPrim),
  fDist3toPrim(source.fDist3toPrim),
  fDist4toPrim(source.fDist4toPrim)
{
  //
  // Copy constructor
  //
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF4Prong &AliAODRecoDecayHF4Prong::operator=(const AliAODRecoDecayHF4Prong &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;

  AliAODRecoDecayHF::operator=(source);

  fDist12toPrim= source.fDist12toPrim;
  fDist3toPrim= source.fDist3toPrim;
  fDist4toPrim= source.fDist4toPrim;
  //fSigmaVert= source.fSigmaVert;

  return *this;
}
//--------------------------------------------------------------------------
void AliAODRecoDecayHF4Prong::InvMassD0(Double_t mD0[2]) const {
  //
  // Mass for the two D0 hypotheses
  //
  UInt_t pdg[4];
  pdg[0]=211; pdg[1]=321; pdg[2]=211; pdg[3]=211;
  mD0[0]=InvMass(4,pdg);
  pdg[1]=211; pdg[3]=321;
  mD0[1]=InvMass(4,pdg);
  
  return;
}
//--------------------------------------------------------------------------
void AliAODRecoDecayHF4Prong::InvMassD0bar(Double_t mD0bar[2]) const {
  //
  // Mass for the two D0bar hypotheses
  //
  UInt_t pdg[4];
  pdg[0]=321; pdg[1]=211; pdg[2]=211; pdg[3]=211;
  mD0bar[0]=InvMass(4,pdg);
  pdg[0]=211; pdg[2]=321;
  mD0bar[1]=InvMass(4,pdg);
  
  return;
}
//--------------------------------------------------------------------------

Bool_t AliAODRecoDecayHF4Prong::SelectD0(const Double_t *cuts,Int_t &okD0,Int_t &okD0bar) const
{
  //
  // This function compares the D0 with a set of cuts:
  // 
  // cuts[0] = D0 invariant mass 
  // cuts[1] = DCA between opposite sign tracks 
  // cuts[2] = Distance between primary and two tracks vertex fDist12toPrim
  // cuts[3] = Distance between primary and three tracks vertex fDist3toPrim
  // cuts[4] = Distance between primary and two tracks vertex fDist4toPrim
  // cuts[5] = Cosinus of the pointing angle
  // cuts[6] = Transverse momentum of the D0 candidate
  // cuts[7] = Mass Pi+Pi- = mass of the rho0
  // cuts[8] = PID cut (one K in the quadruplet)
  //
  // If candidate D0 does not pass the cuts return kFALSE
  //
  
  okD0=0; okD0bar=0;
  Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  Double_t mD0[2];
  Double_t mD0bar[2];
  InvMassD0(mD0);
  InvMassD0bar(mD0bar);
  Bool_t goodMass=kFALSE;
  if(TMath::Abs(mD0[0]-mD0PDG)<=cuts[0]) {goodMass=kTRUE; okD0=1;}
  if(TMath::Abs(mD0[1]-mD0PDG)<=cuts[0]) {goodMass=kTRUE; okD0=1;}
  if(TMath::Abs(mD0bar[0]-mD0PDG)<=cuts[0]) {goodMass=kTRUE; okD0bar=1;}
  if(TMath::Abs(mD0bar[1]-mD0PDG)<=cuts[0]) {goodMass=kTRUE; okD0bar=1;}
  if(!goodMass) return kFALSE;
  
  //DCA btw opposite sign tracks
  if(cuts[1]>0.){
    if(GetDCA(0)>cuts[1]) return kFALSE;
    if(GetDCA(1)>cuts[1]) return kFALSE;
    if(GetDCA(2)>cuts[1]) return kFALSE;
    if(GetDCA(3)>cuts[1]) return kFALSE;
    if(GetDCA(4)>cuts[1]) return kFALSE;
    if(GetDCA(5)>cuts[1]) return kFALSE;
  }
  //2track cuts
  if(cuts[2]>0.){
    if(fDist12toPrim>10.)return kFALSE;
    if(fDist12toPrim<cuts[2])return kFALSE;
  }
  //3track cuts
  if(cuts[3]>0.){
    if(fDist3toPrim<cuts[3])return kFALSE;
  }
  //4track cuts
  if(cuts[4]>0.){
    if(fDist4toPrim<cuts[4])return kFALSE;
  }
  if(cuts[5]>-1.1){
    if(CosPointingAngle()<cuts[5])return kFALSE;
  }
  if(cuts[6]>0.){
    if(Pt()<cuts[6])return kFALSE;
  }
  if(cuts[7]>0.){
    Double_t massD0[2];
    Double_t massD0bar[2];
    Bool_t good=CutRhoMass(massD0,massD0bar,cuts[0],cuts[7]);
    if(!good) return kFALSE;
  }
  
  return kTRUE;
}
//----------------------------------------------------------------------------
Bool_t AliAODRecoDecayHF4Prong::CutRhoMass(Double_t massD0[2],Double_t massD0bar[2],Double_t cutMass,Double_t cutRho) const 
{
  //
  // Cut on rho->pipi mass for any of the pairs
  //
  Bool_t isGood=kFALSE;
  Int_t nprongs=4;
  for(Int_t i=0;i<2;i++){massD0[i]=0.;massD0bar[i]=0.;}
  Bool_t isRho=kFALSE;
  Bool_t isTrue=kFALSE;
  Double_t mPDG=TDatabasePDG::Instance()->GetParticle(421)->Mass();
  Double_t mPDGRho=TDatabasePDG::Instance()->GetParticle(113)->Mass();
  Double_t minv01=InvMassRho(0,1);
  if(TMath::Abs(minv01-mPDGRho)<cutRho) isRho=kTRUE;
  if(isRho){
    UInt_t pdg1[4]={211,211,321,211};
    Double_t mass1=InvMass(nprongs,pdg1);
    if(TMath::Abs(mass1-mPDG)<cutMass) {isTrue=kTRUE;isGood=kTRUE;}
    if(isTrue) massD0bar[1]=mass1;
    isTrue=kFALSE;
    UInt_t pdg2[4]={211,211,211,321};
    Double_t mass2=InvMass(4,pdg2);
    if(TMath::Abs(mass2-mPDG)<cutMass) {isTrue=kTRUE;isGood=kTRUE;}
    if(isTrue) massD0[1]=mass2;
    isTrue=kFALSE;
  }
  Double_t minv03=InvMassRho(0,3);
  if(TMath::Abs(minv03-mPDGRho)<cutRho) isRho=kTRUE;
  if(isRho){
    UInt_t pdg1[4]={211,211,321,211};
    Double_t mass1=InvMass(4,pdg1);
    if(TMath::Abs(mass1-mPDG)<cutMass) {isTrue=kTRUE;isGood=kTRUE;}
    if(isTrue) massD0bar[1]=mass1;
    isTrue=kFALSE;
    UInt_t pdg2[4]={211,321,211,211};
    Double_t mass2=InvMass(4,pdg2);
    if(TMath::Abs(mass2-mPDG)<cutMass) {isTrue=kTRUE;isGood=kTRUE;}
    if(isTrue) massD0[0]=mass2;
    isTrue=kFALSE;
  }
  Double_t minv12=InvMassRho(1,2);
  if(TMath::Abs(minv12-mPDGRho)<cutRho) isRho=kTRUE;
  if(isRho){
    UInt_t pdg1[4]={321,211,211,211};
    Double_t mass1=InvMass(4,pdg1);
    if(TMath::Abs(mass1-mPDG)<cutMass) {isTrue=kTRUE;isGood=kTRUE;}
    if(isTrue) massD0bar[0]=mass1;
    isTrue=kFALSE;
    UInt_t pdg2[4]={211,211,211,321};
    Double_t mass2=InvMass(4,pdg2);
    if(TMath::Abs(mass2-mPDG)<cutMass) {isTrue=kTRUE;isGood=kTRUE;}
    if(isTrue) massD0[1]=mass2;
    isTrue=kFALSE;
  }
  Double_t minv23=InvMassRho(2,3);
  if(TMath::Abs(minv23-mPDGRho)<cutRho) isRho=kTRUE;
  if(isRho){
    UInt_t pdg1[4]={321,211,211,211};
    Double_t mass1=InvMass(4,pdg1);
    if(TMath::Abs(mass1-mPDG)<cutMass) {isTrue=kTRUE;isGood=kTRUE;}
    if(isTrue) massD0bar[0]=mass1;
    isTrue=kFALSE;
    UInt_t pdg2[4]={211,321,211,211};
    Double_t mass2=InvMass(4,pdg2);
    if(TMath::Abs(mass2-mPDG)<cutMass) {isTrue=kTRUE;isGood=kTRUE;}
    if(isTrue) massD0[0]=mass2;
    isTrue=kFALSE;
  }
  return isGood;
}
