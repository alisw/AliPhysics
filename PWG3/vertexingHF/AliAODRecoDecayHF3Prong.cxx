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
// Base class for AOD reconstructed heavy-flavour 3-prong decay
//
// Author: E.Bruna bruna@to.infn.it, F.Prino prino@to.infn.it
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODTrack.h"

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
  SetDCAs(3,dca);
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
  SetDCAs(3,dca);
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

  AliAODRecoDecayHF::operator=(source);

  fDist12toPrim= source.fDist12toPrim;
  fDist23toPrim= source.fDist23toPrim;
  fSigmaVert= source.fSigmaVert;

  return *this;
}
//-------------Temporary fix: overload AliAODRecoDecay::Charge()-----------
Short_t AliAODRecoDecayHF3Prong::Charge() const {
  Short_t charge=0;
  for(Int_t daught=0;daught<=2;daught++)charge=charge+((AliAODTrack*)GetDaughter(daught))->Charge();
  return charge;
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
  for(Int_t i=0;i<3;i++) if(GetDCA(i)>cuts[11])return kFALSE;

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
//--------------------------------------------------------------------------
Bool_t AliAODRecoDecayHF3Prong::SelectDs(const Double_t *cuts,Int_t &okDsKKpi,Int_t &okDspiKK, Int_t &okMassPhi, Int_t &okMassK0star)
  const {
//
// This function compares the Ds with a set of cuts 
// (same variables as D+, for now)
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
// cuts[12] = max. inv. mass difference(Mphi-MKK) [GeV] 
// cuts[13] = max. inv. mass difference(MK0*-MKpi) [GeV] 
//
// If candidate Ds does not pass the cuts return kFALSE
//
  Double_t mDsKKpi,mDspiKK;
  okDsKKpi=1; okDspiKK=1;
  okMassPhi=0; okMassK0star=0;

  Double_t mDsPDG = TDatabasePDG::Instance()->GetParticle(431)->Mass();

  mDsKKpi=InvMassDsKKpi();
  mDspiKK=InvMassDspiKK();

  if(TMath::Abs(mDsKKpi-mDsPDG)>cuts[0]) okDsKKpi = 0;
  if(TMath::Abs(mDspiKK-mDsPDG)>cuts[0]) okDspiKK = 0;
  if(!okDsKKpi && !okDspiKK) return kFALSE;

  //single track
  if(TMath::Abs(PtProng(0)) < cuts[1] || TMath::Abs(Getd0Prong(0))<cuts[3])return kFALSE;//Kaon1
  if(TMath::Abs(PtProng(1)) < cuts[1] || TMath::Abs(Getd0Prong(1))<cuts[3])return kFALSE;//Kaon2
  if(TMath::Abs(PtProng(2)) < cuts[2] || TMath::Abs(Getd0Prong(2))<cuts[4])return kFALSE;//Pion

  // cuts on resonant decays (via Phi or K0*)
  Double_t mPhiPDG = TDatabasePDG::Instance()->GetParticle(333)->Mass();
  Double_t mK0starPDG = TDatabasePDG::Instance()->GetParticle(313)->Mass();
  if(okDsKKpi){
    Double_t mass01phi=InvMass2Prongs(0,1,321,321);
    Double_t mass12K0s=InvMass2Prongs(1,2,321,211);
    if(TMath::Abs(mass01phi-mPhiPDG)<cuts[12]) okMassPhi=1;
    if(TMath::Abs(mass12K0s-mK0starPDG)<cuts[13]) okMassK0star = 1;
    if(!okMassPhi && !okMassK0star) okDsKKpi=kFALSE;
  }
  if(okDspiKK){
    Double_t mass01K0s=InvMass2Prongs(0,1,211,321);
    Double_t mass12phi=InvMass2Prongs(1,2,321,321);
    if(TMath::Abs(mass01K0s-mK0starPDG)<cuts[13]) okMassK0star = 1;
    if(TMath::Abs(mass12phi-mPhiPDG)<cuts[12]) okMassPhi=1;
    if(!okMassPhi && !okMassK0star) okDspiKK=kFALSE;
  }
  if(!okDsKKpi && !okDspiKK) return kFALSE;



  
  //DCA
  for(Int_t i=0;i<3;i++) if(GetDCA(i)>cuts[11])return kFALSE;

  //2track cuts
  if(fDist12toPrim<cuts[5] || fDist23toPrim<cuts[5])return kFALSE;

  //sec vert
  if(fSigmaVert>cuts[6])return kFALSE;

  if(DecayLength()<cuts[7])return kFALSE;

  if(TMath::Abs(PtProng(0))<cuts[8] && TMath::Abs(PtProng(1))<cuts[8] && TMath::Abs(PtProng(2))<cuts[8])return kFALSE;
  if(CosPointingAngle()   < cuts[9])return kFALSE;
  Double_t sum2=Getd0Prong(0)*Getd0Prong(0)+Getd0Prong(1)*Getd0Prong(1)+Getd0Prong(2)*Getd0Prong(2);
  if(sum2<cuts[10])return kFALSE;

  return kTRUE;
}
//--------------------------------------------------------------------------
Bool_t AliAODRecoDecayHF3Prong::SelectLc(const Double_t *cuts,Int_t &okLcpKpi,Int_t &okLcpiKp)
  const {
//
// This function compares the Lc with a set of cuts 
// (same variables as D+, for now)
//
// cuts[0] = inv. mass half width [GeV]   
// cuts[1] = pTP [GeV/c]
// cuts[2] = pTPi and pTK [GeV/c]
// cuts[3] = d0P [cm]   lower limit!
// cuts[4] = d0Pi and d0K [cm]  lower limit!
// cuts[5] = dist12 (cm)
// cuts[6] = sigmavert (cm)
// cuts[7] = dist prim-sec (cm)
// cuts[8] = pM=Max{pT1,pT2,pT3} (GeV/c)
// cuts[9] = cosThetaPoint
// cuts[10] = Sum d0^2 (cm^2)
// cuts[11] = dca cut (cm)
//
// If candidate Lc does not pass the cuts return kFALSE
//
  Double_t mLcpKpi,mLcpiKp;
  okLcpKpi=1; okLcpiKp=1;

  Double_t mLcPDG = TDatabasePDG::Instance()->GetParticle(4122)->Mass();

  mLcpKpi=InvMassLcpKpi();
  mLcpiKp=InvMassLcpiKp();

  if(TMath::Abs(mLcpKpi-mLcPDG)>cuts[0]) okLcpKpi = 0;
  if(TMath::Abs(mLcpiKp-mLcPDG)>cuts[0]) okLcpiKp = 0;
  if(!okLcpKpi && !okLcpiKp) return kFALSE;

  //single track
  if(TMath::Abs(PtProng(0)) < cuts[1] || TMath::Abs(Getd0Prong(0))<cuts[3])return kFALSE;//Proton
  if(TMath::Abs(PtProng(1)) < cuts[2] || TMath::Abs(Getd0Prong(1))<cuts[4])return kFALSE;//Kaon
  if(TMath::Abs(PtProng(2)) < cuts[2] || TMath::Abs(Getd0Prong(2))<cuts[4])return kFALSE;//Pion

  //DCA
  for(Int_t i=0;i<3;i++) if(GetDCA(i)>cuts[11])return kFALSE;

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
