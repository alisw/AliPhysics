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

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliHMPIDPid                                                          //
//                                                                      //
// HMPID class to perfom particle identification                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "AliHMPIDPid.h"       //class header
#include "AliHMPIDParam.h"     //class header
#include <AliESDEvent.h>       //FindPid()

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDPid::AliHMPIDPid():TTask("HMPIDrec","HMPIDPid")
{
//..
//init of data members
//..
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDPid::FindPid(AliESDEvent *pESD)
{
// Calculates probability to be a electron-muon-pion-kaon-proton
// from the given Cerenkov angle and momentum assuming no initial particle composition
// (i.e. apriory probability to be the particle of the given sort is the same for all sorts)

  AliPID ppp; //needed
  Double_t pid[AliPID::kSPECIES],h[AliPID::kSPECIES];
   
  for(Int_t iTrk=0;iTrk<pESD->GetNumberOfTracks();iTrk++){//ESD tracks loop
    AliESDtrack *pTrk = pESD->GetTrack(iTrk);// get next reconstructed track
    if(pTrk->GetHMPIDsignal()<=0){//HMPID does not find anything reasonable for this track, assign 0.2 for all species
      for(Int_t iPart=0;iPart<AliPID::kSPECIES;iPart++) pid[iPart]=1.0/AliPID::kSPECIES;
      pTrk->SetHMPIDpid(pid);
      continue;
    } 
    Double_t pmod = pTrk->GetP();
    Double_t hTot=0;
    for(Int_t iPart=0;iPart<AliPID::kSPECIES;iPart++){
      Double_t mass = AliPID::ParticleMass(iPart);
      Double_t cosThetaTh = TMath::Sqrt(mass*mass+pmod*pmod)/(AliHMPIDParam::Instance()->MeanIdxRad()*pmod);
      if(cosThetaTh<1) //calculate the height of theoretical theta ckov on the gaus of experimental one
        h[iPart] =TMath::Gaus(TMath::ACos(cosThetaTh),pTrk->GetHMPIDsignal(),TMath::Sqrt(pTrk->GetHMPIDchi2()),kTRUE);      
      else             //beta < 1/ref. idx. => no light at all  
        h[iPart] =0 ;       
      hTot    +=h[iPart]; //total height of all theoretical heights for normalization
    }//species loop
     
    Double_t hMin=TMath::Gaus(pTrk->GetHMPIDsignal()-4*TMath::Sqrt(pTrk->GetHMPIDchi2()),pTrk->GetHMPIDsignal(),TMath::Sqrt(pTrk->GetHMPIDchi2()),kTRUE);//5 sigma protection
    
    for(Int_t iPart=0;iPart<AliPID::kSPECIES;iPart++)//species loop to assign probabilities
      if(hTot>hMin)  
        pid[iPart]=h[iPart]/hTot;
      else                               //all theoretical values are far away from experemental one
        pid[iPart]=1.0/AliPID::kSPECIES; 
    pTrk->SetHMPIDpid(pid);
  }//ESD tracks loop
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
