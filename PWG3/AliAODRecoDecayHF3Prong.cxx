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

/////////////////////////////////////////////////////////////
//
// Base class for AOD reconstructed heavy-flavour 3-prong decay
//
// Author: E.Bruna bruna@to.infn.it, F.Prino prino@to.infn.it
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF3Prong.h"

ClassImp(AliAODRecoDecayHF3Prong)

//--------------------------------------------------------------------------
AliAODRecoDecayHF3Prong::AliAODRecoDecayHF3Prong() :
  AliAODRecoDecayHF(), 
  fSigmaVert(0),
  fDist12toPrim(0),
  fDist23toPrim(0)
{
  //
  // Default Constructor
  //
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF3Prong::AliAODRecoDecayHF3Prong(AliAODVertex *vtx2,
						 Double_t *px,Double_t *py,Double_t *pz,
						 Double_t *d0,Double_t *d0err,
						 Double_t *dca, Double_t sigvert,
						 Double_t dist12,Double_t dist23,Short_t charge) :
  AliAODRecoDecayHF(vtx2,3,charge,px,py,pz,d0,d0err),
  fSigmaVert(sigvert),
  fDist12toPrim(dist12),
  fDist23toPrim(dist23)
{
  //
  // Constructor with AliAODVertex for decay vertex
  //
  Double_t dcafloat[3];
  for(Int_t i=0;i<3;i++)dcafloat[i]=dca[i];
  SetDCAs(3,dcafloat);
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF3Prong::AliAODRecoDecayHF3Prong(AliAODVertex *vtx2,
						 Double_t *d0,Double_t *d0err,
						 Double_t *dca, Double_t sigvert,
						 Double_t dist12,Double_t dist23, Short_t charge) :
  AliAODRecoDecayHF(vtx2,3,charge,d0,d0err),
  fSigmaVert(sigvert),
  fDist12toPrim(dist12),
  fDist23toPrim(dist23)
{
  //
  // Constructor with AliAODVertex for decay vertex and without prongs momenta
  //
  Double_t dcafloat[3];
  for(Int_t i=0;i<3;i++)dcafloat[i]=dca[i];
  SetDCAs(3,dcafloat);
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF3Prong::AliAODRecoDecayHF3Prong(const AliAODRecoDecayHF3Prong &source) :
  AliAODRecoDecayHF(source),
  fSigmaVert(source.fSigmaVert),
  fDist12toPrim(source.fDist12toPrim),
  fDist23toPrim(source.fDist23toPrim)
{
  //
  // Copy constructor
  //
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF3Prong &AliAODRecoDecayHF3Prong::operator=(const AliAODRecoDecayHF3Prong &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;
  fOwnPrimaryVtx = source.fOwnPrimaryVtx;
  fSecondaryVtx = source.fSecondaryVtx;
  fCharge = source.fCharge;
  fNProngs = source.fNProngs;
  fNDCA = source.fNDCA;
  fNPID = source.fNPID;
  fEventNumber = source.fEventNumber;
  fRunNumber = source.fRunNumber;
  fDist12toPrim= source.fDist12toPrim;
  fDist23toPrim= source.fDist23toPrim;
  fSigmaVert= source.fSigmaVert;
  if(source.GetNProngs()>0) {
    fd0 = new Double_t[GetNProngs()];
    fd0err = new Double_t[GetNProngs()];
    memcpy(fd0,source.fd0,GetNProngs()*sizeof(Double_t));
    memcpy(fd0err,source.fd0err,GetNProngs()*sizeof(Double_t));
    if(source.fPx) {
      fPx = new Double_t[GetNProngs()];
      fPy = new Double_t[GetNProngs()];
      fPz = new Double_t[GetNProngs()];
      memcpy(fPx,source.fPx,GetNProngs()*sizeof(Double_t));
      memcpy(fPy,source.fPy,GetNProngs()*sizeof(Double_t));
      memcpy(fPz,source.fPz,GetNProngs()*sizeof(Double_t));
    }
    if(source.fPID) {
      fPID = new Double_t[5*GetNProngs()];
      memcpy(fPID,source.fPID,GetNProngs()*sizeof(Double_t));
    }
    if(source.fDCA) {
      fDCA = new Double32_t[GetNProngs()*(GetNProngs()-1)/2];
      memcpy(fDCA,source.fDCA,(GetNProngs()*(GetNProngs()-1)/2)*sizeof(Float_t));
    }
  }
  return *this;
}
//--------------------------------------------------------------------------
Bool_t AliAODRecoDecayHF3Prong::SelectDplus(const Double_t *cuts)
  const {
//
// This function compares the Dplus with a set of cuts:
//
// cuts[0] = inv. mass half width [GeV]   
// cuts[1] = pTK [GeV/c]
// cuts[2] = pTPi [GeV/c]
// cuts[3] = d0K [cm]   lower limit!
// cuts[4] = d0Pi [cm]  lower limit!
// cuts[5] = dist12 (cm)
// cuts[6] = sigmavert (cm)
// cuts[7] = dist prim-sec (cm)
// cuts[8] = pM=Max{pT1,pT2,pT3} (GeV/c)
// cuts[9] = cosThetaPoint
// cuts[10] = Sum d0^2 (cm^2)
// cuts[11] = dca cut (cm)
//
// If candidate Dplus does not pass the cuts return kFALSE
//

  Double_t mDplusPDG = TDatabasePDG::Instance()->GetParticle(411)->Mass();
  Double_t mDplus=InvMassDplus();
  if(TMath::Abs(mDplus-mDplusPDG)>cuts[0])return kFALSE;
  //single track
  if(TMath::Abs(PtProng(1)) < cuts[1] || TMath::Abs(Getd0Prong(1))<cuts[3])return kFALSE;//Kaon
  if(TMath::Abs(PtProng(0)) < cuts[2] || TMath::Abs(Getd0Prong(0))<cuts[4])return kFALSE;//Pion1
  if(TMath::Abs(PtProng(2)) < cuts[2] || TMath::Abs(Getd0Prong(2))<cuts[4])return kFALSE;//Pion2

  //DCA
  for(Int_t i=0;i<3;i++) if(cuts[11]>0 && GetDCA(i)>cuts[11])return kFALSE;

  //2track cuts
  if(fDist12toPrim<cuts[5] || fDist23toPrim<cuts[5])return kFALSE;
  if(Getd0Prong(0)*Getd0Prong(1)<0. && Getd0Prong(2)*Getd0Prong(1)<0.)return kFALSE;

  //sec vert
  if(fSigmaVert>cuts[6])return kFALSE;

  if(DecayLength()<cuts[7])return kFALSE;

  if(TMath::Abs(PtProng(0))<cuts[8] && TMath::Abs(PtProng(1))<cuts[8] && TMath::Abs(PtProng(2))<cuts[8])return kFALSE;
  if(CosPointingAngle()   < cuts[9])return kFALSE;
  Double_t sum2=Getd0Prong(0)*Getd0Prong(0)+Getd0Prong(1)*Getd0Prong(1)+Getd0Prong(2)*Getd0Prong(2);
  if(sum2<cuts[10])return kFALSE;
  return kTRUE;
}
