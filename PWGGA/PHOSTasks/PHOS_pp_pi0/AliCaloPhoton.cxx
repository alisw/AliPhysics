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
 
//_________________________________________________________________________
// Minimal class to store photon infomation for pi0, tagged and isolated photon analysis
//
//-- Author: Dmitri Peressounko (RRC "KI")

#include "AliCaloPhoton.h"
ClassImp(AliCaloPhoton) 
//===============================================
AliCaloPhoton::AliCaloPhoton() :
  TLorentzVector(),
  fMomV2(),
  fDisp(0),
  fDisp2(0),
  fTof(0),
  fCpv(0),
  fCpv2(0),
  fPCA(0),
  fTrig(0),
  fIsTagged(0),
  fIsIsolated(0),
  fIsPhoton(0),
  fUnfolded(0),
  fModule(0),
  fBC(0),
  fBadDist(0),
  fBadDistfp(0.),
  fNCells(0),
  fFiducialArea(0),
  fPi0Decayflag(0),
  fPi0Id(0),
  fConverted(0),
  fConvertedPartner(0),
  fIsolationTag(0),
  fTagInfo(0),
  fPrimary(-1),
  fPrimaryAtVertex(-1),
  fX(0.),
  fY(0.),
  fZ(0.),
  fLambda0(0.),
  fLambda1(0.),
  fTime(0.),
  fPartnerPt(0),
  fWeight(1.),
  fNsigmaCPV(-1),
  fNsigmaFullDisp(-1),
  fNsigmaCoreDisp(-1),
  fTOFCutEfficiency(1.),
  fEmbEventID(-1),
  fCluster(0x0)
{
  ResetTagWeights() ;

} 
//===============================================
AliCaloPhoton::AliCaloPhoton(Double_t px,Double_t py,Double_t pz,Double_t energy):
  TLorentzVector(px,py,pz,energy),
  fMomV2(),
  fDisp(0),
  fDisp2(0),
  fTof(0),
  fCpv(0),
  fCpv2(0),
  fPCA(0),
  fTrig(0),
  fIsTagged(0),
  fIsIsolated(0),
  fIsPhoton(0),
  fUnfolded(0),
  fModule(0),
  fBC(0),
  fBadDist(0),
  fBadDistfp(0.),
  fNCells(0),
  fFiducialArea(0),
  fPi0Decayflag(0),
  fPi0Id(0),
  fConverted(0),
  fConvertedPartner(0),
  fIsolationTag(0),
  fTagInfo(0),
  fPrimary(-1),
  fPrimaryAtVertex(-1),
  fX(0.),
  fY(0.),
  fZ(0.),
  fLambda0(0.),
  fLambda1(0.),
  fTime(0.),
  fPartnerPt(0),
  fWeight(1.),
  fNsigmaCPV(-1),
  fNsigmaFullDisp(-1),
  fNsigmaCoreDisp(-1),
  fTOFCutEfficiency(1.),
  fEmbEventID(-1),
  fCluster(0x0)
{
    ResetTagWeights() ;
}
//===============================================
Bool_t AliCaloPhoton::IsPIDOK(Int_t ipid)const{
  // returns true if photon satisfies given PID criterium
  switch(ipid){
  case 0: return kTRUE ; //No PID at all
  case 1: return fPCA; //Overall PID calculated in AliPHOSPIDv1
  case 2: return fDisp ;   //only dispersion cut
  case 3: return fTof ;    //Only TOF cut
  case 4: return fCpv ;    //Only CPV cut
  case 5: return fDisp&&fTof ;  //Dispersion and TOF
  case 6: return fDisp&&fCpv ;  //Dispersion and CPV
  case 7: return fTof && fCpv;  //TOF and CPV
  case 8: return fDisp&&fTof&&fCpv ; // all 3 cuts
  default: return kFALSE ; //Not known combination
  }
}
