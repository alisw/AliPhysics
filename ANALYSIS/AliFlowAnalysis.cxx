#include "AliFlowAnalysis.h"
//________________________________
///////////////////////////////////////////////////////////
//
// class AliFlowAnalysis
//
// Flow Analysis
//
//
// S.Radomski@gsi.de
// Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////
/*********************************************************/

#include <TString.h>
#include <TParticle.h>

#include <AliStack.h>
#include <AliAOD.h>
#include <AliVAODParticle.h>
#include <AliAODParticleCut.h>

#include <AliESDtrack.h>
#include <AliESD.h>

ClassImp(AliFlowAnalysis)

AliFlowAnalysis::AliFlowAnalysis():
 fPartCut(0x0)
{
 //ctor
}
/*********************************************************/

AliFlowAnalysis::~AliFlowAnalysis()
{
 //dtor
  delete fPartCut;
}
/*********************************************************/

Int_t AliFlowAnalysis::Init()
{
  //Initilizes anaysis
  Info("Init","");
  return 0; 
}
/*********************************************************/

Int_t AliFlowAnalysis::ProcessEvent(AliAOD* aodrec, AliAOD* aodsim)
{
 
  Info("ProcessEvent","Sim AOD address %#x",aodsim);
  Double_t psi = 0, v2 = 0;
  if (aodrec)
   {
     GetFlow(aodrec,v2,psi);
     Info("ProcessEvent","Reconstructed Event: Event plane is %f, V2 is %f",psi,v2);
   }  

  if (aodsim)
   {
     GetFlow(aodsim,v2,psi);
     Info("ProcessEvent","Simulated Event: Event plane is %f, V2 is %f",psi,v2);
   }  
  
  return 0;
  
}
/*********************************************************/

Int_t AliFlowAnalysis::Finish()
{
  //Finish analysis and writes results
  Info("Init","Finish");
  return 0;
}
/*********************************************************/

Double_t AliFlowAnalysis::GetEventPlane(AliAOD* aod)
{
  //returns event plane in degrees
  if (aod == 0x0)
   {
     Error("AliFlowAnalysis::GetFlow","Pointer to AOD is NULL");
     return -1.0;
   }

  Double_t psi;
  Int_t mult = aod->GetNumberOfParticles();
  
  Double_t ssin = 0, scos = 0;

  for (Int_t i=0; i<mult; i++) 
   {
     AliVAODParticle* aodtrack = aod->GetParticle(i);
     if (aodtrack == 0x0)
      {
        Error("AliFlowAnalysis::GetEventPlane","Can not get track %d", i);
        continue;
      }
     
     if (fPartCut)
      if (fPartCut->Rejected(aodtrack))
        continue;

     Double_t phi = TMath::Pi()+TMath::ATan2(-aodtrack->Py(),-aodtrack->Px()); 
     
     ssin += TMath::Sin( 2.0 * phi );
     scos += TMath::Cos( 2.0 * phi );
   }

  psi = atan2 (ssin, scos) / 2.0;
  psi = psi * 180. / TMath::Pi(); 
  
  return psi;

}
/*********************************************************/

void AliFlowAnalysis::GetFlow(AliAOD* aod,Double_t& v2,Double_t& psi)
{
//returns flow parameters: v2 and event plane
  if (aod == 0x0)
   {
     Error("AliFlowAnalysis::GetFlow","Pointer to AOD is NULL");
     return;
   }
   
  psi = GetEventPlane(aod);
  Int_t mult = aod->GetNumberOfParticles();
  
  Double_t ssin = 0, scos = 0;

  for (Int_t i=0; i<mult; i++) 
   {
     AliVAODParticle* aodtrack = aod->GetParticle(i);
     if (aodtrack == 0x0)
      {
        Error("AliFlowAnalysis::GetEventPlane","Can not get track %d", i);
        continue;
      }
     if (fPartCut)
      if (fPartCut->Rejected(aodtrack))
        continue;
      
     Double_t phi = TMath::Pi()+TMath::ATan2(-aodtrack->Py(),-aodtrack->Px()); 
     ssin += TMath::Sin( 2.0 * (phi - psi));
     scos += TMath::Cos( 2.0 * (phi - psi));
   }
   
  v2 = TMath::Hypot(ssin,scos);
}


/*********************************************************/

Double_t AliFlowAnalysis::GetEventPlane(AliESD* esd)
{
  //returns event plane
  if (esd == 0x0)
   {
     ::Error("AliFlowAnalysis::GetFlow","Pointer to ESD is NULL");
     return -1.0;
   }

  Double_t psi;
  Int_t mult = esd->GetNumberOfTracks();
  
  Double_t ssin = 0, scos = 0;

  for (Int_t i=0; i<mult; i++) 
   {
     AliESDtrack* esdtrack = esd->GetTrack(i);
     if (esdtrack == 0x0)
      {
        ::Error("AliFlowAnalysis::GetEventPlane","Can not get track %d", i);
        continue;
      }
      
     Double_t mom[3];//momentum
     esdtrack->GetPxPyPz(mom); 
     Double_t phi = TMath::Pi()+TMath::ATan2(-mom[1],-mom[0]); 
     
     ssin += TMath::Sin( 2.0 * phi );
     scos += TMath::Cos( 2.0 * phi );
   }

  psi = atan2 (ssin, scos) / 2.0;
  psi = psi * 180. / TMath::Pi(); 
  
  return psi;

}
/*********************************************************/

void AliFlowAnalysis::GetFlow(AliESD* esd,Double_t& v2,Double_t& psi)
{
//returns flow parameters: v2 and event plane
  if (esd == 0x0)
   {
     ::Error("AliFlowAnalysis::GetFlow","Pointer to ESD is NULL");
     return;
   }
   
  psi = GetEventPlane(esd);
  Int_t mult = esd->GetNumberOfTracks();
  
  Double_t ssin = 0, scos = 0;

  for (Int_t i=0; i<mult; i++) 
   {
     AliESDtrack* esdtrack = esd->GetTrack(i);
     if (esdtrack == 0x0)
      {
        ::Error("AliFlowAnalysis::GetEventPlane","Can not get track %d", i);
        continue;
      }
      
     Double_t mom[3];//momentum
     esdtrack->GetPxPyPz(mom); 
     Double_t phi = TMath::Pi()+TMath::ATan2(-mom[1],-mom[0]); 
     
     ssin += TMath::Sin( 2.0 * (phi - psi));
     scos += TMath::Cos( 2.0 * (phi - psi));
   }
   
  v2 = TMath::Hypot(ssin,scos);
}
/*********************************************************/
