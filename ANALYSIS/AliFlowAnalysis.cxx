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
#include <AliESDtrack.h>
#include <AliESD.h>

ClassImp(AliFlowAnalysis)

Int_t AliFlowAnalysis::Init()
{
  //Initilizes anaysis
  Info("Init","");
  return 0; 
}

/*********************************************************/

Int_t AliFlowAnalysis::ProcessEvent(AliESD* esd, AliStack* stack)
{
 
  Info("ProcessEvent","Stack address %#x",stack);
  Double_t psi = GetEventPlane(esd);
  Info("ProcessEvent","Event plane is %f",psi);
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

Double_t AliFlowAnalysis::GetEventPlane(AliESD* esd)
{
  //returns event plane
  Double_t psi;
  Int_t mult = esd->GetNumberOfTracks();
  
  Double_t ssin = 0, scos = 0;

  for (Int_t i=0; i<mult; i++) 
   {
     AliESDtrack* esdtrack = esd->GetTrack(i);
     if (esdtrack == 0x0)
      {
        ::Error("GetEventPlane","Can not get track %d", i);
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
