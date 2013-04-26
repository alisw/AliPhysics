/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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



//
// Class to handle the AOD tracks with good HMPID data 
// Author: Levente Molnar
// levente.molnar@cern.ch , March 2012
// 

#include "AliAODHMPIDrings.h"

ClassImp(AliAODHMPIDrings)

//________________________________________________________________________________________________________________________________________________________

AliAODHMPIDrings::AliAODHMPIDrings()
                   :TObject(),
                    fHmpidAODtrkId(0),
                    fHmpidAODqn(0),
                    fHmpidAODcluIdx(0),
                    fHmpidAODtrkTheta(0),
                    fHmpidAODtrkPhi(0),
                    fHmpidAODsignal(0),
                    fHmpidAODocc(0),
                    fHmpidAODchi2(0),
                    fHmpidAODtrkX(0),
                    fHmpidAODtrkY(0),
                    fHmpidAODmipX(0),
                    fHmpidAODmipY(0)

{
  //default ctor 
   for(Int_t isp = 0 ; isp <AliPID::kSPECIES; isp++)   fHmpidAODpid[isp] = 0;
       for ( Int_t ico = 0 ; ico < 3; ico++) fHMPIDmom[ico] = 0;
}

//________________________________________________________________________________________________________________________________________________________
AliAODHMPIDrings::AliAODHMPIDrings(
                    Int_t trkId,
                    Int_t qn, 
                    Int_t cluIdx,
                    Double_t  trkTheta,
                    Double_t trkPhi,
                    Double_t signal,
                    Double_t occ,
                    Double_t chi2,
                    Double_t trkX,
                    Double_t trkY,
                    Double_t mipX,
                    Double_t mipY,
                    Double_t *pid,
                    Double_t *p         ):
                    TObject(),
                    
                    fHmpidAODtrkId(trkId),
                    fHmpidAODqn(qn),
                    fHmpidAODcluIdx(cluIdx),
                    fHmpidAODtrkTheta(trkTheta),
                    fHmpidAODtrkPhi(trkPhi),
                    fHmpidAODsignal(signal),
                    fHmpidAODocc(occ),
                    fHmpidAODchi2(chi2),
                    fHmpidAODtrkX(trkX),
                    fHmpidAODtrkY(trkY),
                    fHmpidAODmipX(mipX),
                    fHmpidAODmipY(mipY)

                    
{
       //             
       for(Int_t isp = 0 ; isp <AliPID::kSPECIES; isp++)   fHmpidAODpid[isp] = pid[isp];
       for ( Int_t ico = 0 ; ico < 3; ico++) fHMPIDmom[ico] = p[ico];
                           
}
//________________________________________________________________________________________________________________________________________________________
AliAODHMPIDrings::AliAODHMPIDrings(const AliAODHMPIDrings& hmpidAOD):
    
                    TObject(hmpidAOD),
                    fHmpidAODtrkId(hmpidAOD.fHmpidAODtrkId),
                    fHmpidAODqn(hmpidAOD.fHmpidAODqn),
                    fHmpidAODcluIdx(hmpidAOD.fHmpidAODcluIdx),
                    fHmpidAODtrkTheta(hmpidAOD.fHmpidAODtrkTheta),
                    fHmpidAODtrkPhi(hmpidAOD.fHmpidAODtrkPhi),
                    fHmpidAODsignal(hmpidAOD.fHmpidAODsignal),
                    fHmpidAODocc(hmpidAOD.fHmpidAODocc),
                    fHmpidAODchi2(hmpidAOD.fHmpidAODchi2),
                    fHmpidAODtrkX(hmpidAOD.fHmpidAODtrkX),
                    fHmpidAODtrkY(hmpidAOD.fHmpidAODtrkY),
                    fHmpidAODmipX(hmpidAOD.fHmpidAODmipX),
                    fHmpidAODmipY(hmpidAOD.fHmpidAODmipY)

                    
{
       //             
       for(Int_t isp = 0 ; isp <AliPID::kSPECIES; isp++)   fHmpidAODpid[isp] = hmpidAOD.fHmpidAODpid[isp];
       for ( Int_t ico = 0 ; ico < 3; ico++) fHMPIDmom[ico] = hmpidAOD.fHMPIDmom[ico];
                           
}

//________________________________________________________________________________________________________________________________________________________
AliAODHMPIDrings& AliAODHMPIDrings::operator=(const AliAODHMPIDrings& hmpidAOD)
{
     if (this!=&hmpidAOD) {   
                    AliAODHMPIDrings::operator=(hmpidAOD);  
                    fHmpidAODtrkId = hmpidAOD.fHmpidAODtrkId;        
                    fHmpidAODqn = hmpidAOD.fHmpidAODqn;
                    fHmpidAODcluIdx = hmpidAOD.fHmpidAODcluIdx;
                    fHmpidAODtrkTheta = hmpidAOD.fHmpidAODtrkTheta;
                    fHmpidAODtrkPhi = hmpidAOD.fHmpidAODtrkPhi;
                    fHmpidAODsignal = hmpidAOD.fHmpidAODsignal;
                    fHmpidAODocc = hmpidAOD.fHmpidAODocc;
                    fHmpidAODchi2 = hmpidAOD.fHmpidAODchi2;
                    fHmpidAODtrkX = hmpidAOD.fHmpidAODtrkX;
                    fHmpidAODtrkY = hmpidAOD.fHmpidAODtrkY;
                    fHmpidAODmipX = hmpidAOD.fHmpidAODmipX;
                    fHmpidAODmipY = hmpidAOD.fHmpidAODmipY;
                    
                    for(Int_t isp = 0 ; isp <AliPID::kSPECIES; isp++)   fHmpidAODpid[isp] = hmpidAOD.fHmpidAODpid[isp];
                    for ( Int_t ico = 0 ; ico < 3; ico++) fHMPIDmom[ico] = hmpidAOD.fHMPIDmom[ico];

     }
     
     return *this;
                                  
}
//________________________________________________________________________________________________________________________________________________________
void AliAODHMPIDrings::GetHmpPidProbs(Double_t *pid) const
{
  // Gets probabilities of each particle type (in HMPID)
  for (Int_t i=0; i<AliPID::kSPECIES; i++) pid[i]=fHmpidAODpid[i];
}
//________________________________________________________________________________________________________________________________________________________
void  AliAODHMPIDrings::GetHmpMom(Double_t *mom) const
{
  for( Int_t ico = 0 ; ico < 3; ico++) mom[ico] = fHMPIDmom[ico];
}
//________________________________________________________________________________________________________________________________________________________
void AliAODHMPIDrings::SetHmpPidProbs(Double_t *pid)
{
  // Gets probabilities of each particle type (in HMPID)
  for (Int_t i=0; i<AliPID::kSPECIES; i++) fHmpidAODpid[i] = pid[i];
}

//________________________________________________________________________________________________________________________________________________________
void  AliAODHMPIDrings::SetHmpMom(Double_t *mom)
{
  for( Int_t ico = 0 ; ico < 3; ico++) fHMPIDmom[ico] = mom[ico];  
}
//________________________________________________________________________________________________________________________________________________________

