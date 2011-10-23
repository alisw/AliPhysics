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
#include "AliHMPIDRecon.h"     //class header
#include <AliESDtrack.h>       //FindPid()
#include <TRandom.h>           //Resolution()

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDPid::AliHMPIDPid():TNamed("HMPIDrec","HMPIDPid")
{
//..
//init of data members
//..
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDPid::FindPid(AliESDtrack *pTrk,Int_t nsp,Double_t *prob)
{
// Calculates probability to be a electron-muon-pion-kaon-proton with the "amplitude" method
// from the given Cerenkov angle and momentum assuming no initial particle composition
// (i.e. apriory probability to be the particle of the given sort is the same for all sorts)

  AliPID *pPid = new AliPID();
  
  Double_t thetaCerExp = -999.;                                                                           
  if(pTrk->GetHMPIDsignal()<=0) thetaCerExp = pTrk->GetHMPIDsignal();                                                                           
  else                          thetaCerExp = pTrk->GetHMPIDsignal() - (Int_t)pTrk->GetHMPIDsignal();     //  measured thetaCherenkov
  
  if(thetaCerExp<=0){                                                                                     //HMPID does not find anything reasonable for this track, assign 0.2 for all species
    for(Int_t iPart=0;iPart<nsp;iPart++) prob[iPart]=1.0/(Float_t)nsp;
    delete pPid ; pPid=0x0; return;
  } 
  
  Double_t p[3] = {0}, pmod = 0;
  if(pTrk->GetOuterHmpPxPyPz(p))  pmod = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);  // Momentum of the charged particle
  
  else {                                         
    for(Int_t iPart=0;iPart<nsp;iPart++) prob[iPart]=1.0/(Float_t)nsp;
    delete pPid ; pPid=0x0; return;
  } 
  
  Double_t hTot=0;                               // Initialize the total height of the amplitude method
  Double_t *h = new Double_t [nsp];              // number of charged particles to be considered

  Bool_t desert = kTRUE;                                                                                                     //  Flag to evaluate if ThetaC is far ("desert") from the given Gaussians
  
  for(Int_t iPart=0;iPart<nsp;iPart++){                                                                                      //  for each particle
    
    h[iPart] = 0;                                                                                                            //  reset the height
    Double_t mass = pPid->ParticleMass(iPart);                                                                             //  with the given mass
    Double_t cosThetaTh = TMath::Sqrt(mass*mass+pmod*pmod)/(AliHMPIDParam::Instance()->MeanIdxRad()*pmod);                   //  evaluate the theor. Theta Cherenkov
    if(cosThetaTh>1) continue;                                                                                               //  no light emitted, zero height
    Double_t thetaCerTh = TMath::ACos(cosThetaTh);                                                                           //  theoretical Theta Cherenkov
    Double_t sigmaRing = Resolution(thetaCerTh,pTrk);
    
    if(sigmaRing==0) {
      for(Int_t jPart=0;jPart<nsp;jPart++) prob[jPart]=1.0/(Float_t)nsp;
      delete pPid ; pPid=0x0; delete [] h; return;
    } 
    
    if(TMath::Abs(thetaCerExp-thetaCerTh)<4*sigmaRing) desert = kFALSE;                                                                //   
    h[iPart] =TMath::Gaus(thetaCerTh,thetaCerExp,sigmaRing,kTRUE);
    hTot    +=h[iPart]; //total height of all theoretical heights for normalization
    
  }//species loop

  for(Int_t iPart=0;iPart<nsp;iPart++) {//species loop to assign probabilities
    
    if(!desert) prob[iPart]=h[iPart]/hTot;
    else prob[iPart]=1.0/(Float_t)nsp;            //all theoretical values are far away from experemental one
    
  }
  
  delete [] h;
  delete pPid ; pPid=0x0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliHMPIDPid::Resolution(Double_t thetaCerTh, AliESDtrack *pTrk)
{
  AliHMPIDParam *pParam = AliHMPIDParam::Instance();
      
  AliHMPIDRecon rec;
  Float_t xRa,yRa,thRa,phRa;
  pTrk->GetHMPIDtrk(xRa,yRa,thRa,phRa);
  rec.SetTrack(xRa,yRa,thRa,phRa);
  Double_t thetaMax = TMath::ACos(1./pParam->MeanIdxRad());
  Int_t nPhots = (Int_t)(21.*TMath::Sin(thetaCerTh)*TMath::Sin(thetaCerTh)/(TMath::Sin(thetaMax)*TMath::Sin(thetaMax))+0.01);

  Double_t sigmatot = 0;
  Int_t nTrks = 20;
  for(Int_t iTrk=0;iTrk<nTrks;iTrk++) {
    Double_t invSigma = 0;
    Int_t nPhotsAcc = 0;
    for(Int_t j=0;j<nPhots;j++){
      Double_t phi = gRandom->Rndm()*TMath::TwoPi();
      TVector2 pos; pos=rec.TracePhot(thetaCerTh,phi);
      if(!pParam->IsInside(pos.X(),pos.Y())) continue;
      if(pParam->IsInDead(pos.X(),pos.Y())) continue;
      Double_t sigma2 = pParam->Sigma2(thRa,phRa,thetaCerTh,phi);//photon candidate sigma^2
      if(sigma2!=0) {
        invSigma += 1./sigma2;
        nPhotsAcc++;
      }
    }      
    if(invSigma!=0) sigmatot += 1./TMath::Sqrt(invSigma);  
  }
  return sigmatot/nTrks;
}
