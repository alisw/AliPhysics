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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for RICH reconstruction                                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliRICHReconstructor.h"
#include "AliRICHClusterFinder.h"
#include "AliRICHHelix.h"
#include <AliRunLoader.h>
#include <AliRun.h>
#include <AliESD.h>

ClassImp(AliRICHReconstructor)

//__________________________________________________________________________________________________
void AliRICHReconstructor::Reconstruct(AliRunLoader* pAL) const
{
//Finds clusters out of digits
  AliDebug(1,"Start cluster finder.");AliRICHClusterFinder clus(GetRICH(pAL));  clus.Exec();
}
//__________________________________________________________________________________________________
void AliRICHReconstructor::FillESD(AliRunLoader* /*pAL*/, AliESD* pESD) const
{
//This methode fills AliESDtrack with information from RICH  
  AliDebug(1,Form("Start with %i tracks",pESD->GetNumberOfTracks()));
/*  const Double_t masses[5]={0.000511,0.105658,0.139567,0.493677,0.93828};//electron,muon,pion,kaon,proton
  const Double_t refIndex = 1.29052;

  Double_t thetaTh[5];
  Double_t sinThetaThNorm;
  Double_t sigmaThetaTh[5];
  Double_t height[5];
  Double_t totalHeight=0;
  Double_t x[3],p[3]; //tmp storage for track parameters
 
  for(Int_t iTrackN=0;iTrackN<pESD->GetNumberOfTracks();iTrackN++){//ESD tracks loop
    AliESDtrack *pTrack = pESD->GetTrack(iTrackN);// get next reconstructed track
//    if((pTrack->GetStatus()&AliESDtrack::kTOFout)==0) continue; //ignore tracks not recontructed by TOF
    pTrack->GetXYZ(x); pTrack->GetPxPyPz(p);          Double_t pmod=pTrack->GetP();//get running track parameters
    continue;  
//  AliRICHHelix helix(x[0],x[1],x[2],p[0],p[1],p[2]);      //construct helix from running track parameters
    
//    TVector rad(1,5); TVector pc(1,5);
//    helix.RichIntersect(GetRICH(pAL)->P(),rad,pc); //returns cross point of track with RICH PC in LRS
    
    for(Int_t iPart=4;iPart>=0;iPart--){
      Double_t cosThetaTh = TMath::Sqrt(masses[iPart]*masses[iPart]+pmod*pmod)/(refIndex*pmod);
      if(cosThetaTh>=1) {break;}
      thetaTh[iPart] = TMath::ACos(cosThetaTh);
      sinThetaThNorm = TMath::Sin(thetaTh[iPart])/TMath::Sqrt(1-1/(refIndex*refIndex));
      sigmaThetaTh[iPart] = (0.014*(1/sinThetaThNorm-1) + 0.0043)*1.25;
      height[iPart] = TMath::Gaus(thetaExp,thetaTh[iPart],sigmaThetaTh[iPart]);
      totalHeight +=height[iPart];
    }
    
    Double_t richPID[5];
    for(Int_t iPart=0;iPart<5;iPart++) richPID[iPart] = height[iPart]/totalHeight;    
    pTrack->SetRICHpid(richPID); 
  }//ESD tracks loop
  */
}//FillESD
//__________________________________________________________________________________________________
AliRICH* AliRICHReconstructor::GetRICH(AliRunLoader* pAL) const
{
// get the RICH detector

  if (!pAL->GetAliRun()) pAL->LoadgAlice();
  if (!pAL->GetAliRun()) {AliError("couldn't get AliRun object"); return NULL;  }
  AliRICH* pRich = (AliRICH*) pAL->GetAliRun()->GetDetector("RICH");
  if (!pRich) {AliError("couldn't get RICH detector");    return NULL;  }
  return pRich;
}
//__________________________________________________________________________________________________
