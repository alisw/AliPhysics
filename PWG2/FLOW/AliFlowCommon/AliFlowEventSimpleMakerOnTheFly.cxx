/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/************************************ 
 * Create an event and perform full *
 * flow analysis 'on the fly'.      * 
 *                                  * 
 * author: Ante Bilandzic           * 
 *         (abilandzic@gmail.com)   *
 ************************************/ 
  
#include "Riostream.h"
#include "TMath.h"
#include "TF1.h"
#include "TRandom3.h"
#include "AliFlowEventSimpleMakerOnTheFly.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"

ClassImp(AliFlowEventSimpleMakerOnTheFly)

//========================================================================================================================================

AliFlowEventSimpleMakerOnTheFly::AliFlowEventSimpleMakerOnTheFly(UInt_t uiSeed):
fCount(0),
fMinMult(0),
fMaxMult(0),  
fPtSpectra(NULL),
fMass(0.13957),
fTemperature(0.44),
fPhiDistribution(NULL),
fV1(0.),
fV2(0.05),
fV3(0.),
fV4(0.),
fUniformFluctuationsV2(kFALSE),
fMinV2(0.04),
fMaxV2(0.06),
fPtDependentV2(kFALSE),
fV2vsPtCutOff(2.0),
fV2vsPtMax(0.2),
fEtaMinA(-0.8),
fEtaMaxA(-0.5),
fEtaMinB(0.5),
fEtaMaxB(0.8),
fNTimes(1),
fUniformAcceptance(kTRUE),
fPhiMin1(0.),              
fPhiMax1(0.),             
fProbability1(0.),       
fPhiMin2(0.),   
fPhiMax2(0.),            
fProbability2(0.),     
fPi(TMath::Pi())
{
 // Constructor.
  
 // Determine seed for gRandom:
 delete gRandom;
 gRandom = new TRandom3(uiSeed); // if uiSeed is 0, the seed is determined uniquely in space and time via TUUID
  
} // end of AliFlowEventSimpleMakerOnTheFly::AliFlowEventSimpleMakerOnTheFly(UInt_t uiSeed):

//====================================================================================================================

AliFlowEventSimpleMakerOnTheFly::~AliFlowEventSimpleMakerOnTheFly()
{
 // Destructor.

 if(fPtSpectra){delete fPtSpectra;}
 if(fPhiDistribution){delete fPhiDistribution;}

} // end of AliFlowEventSimpleMakerOnTheFly::~AliFlowEventSimpleMakerOnTheFly()	

//====================================================================================================================

void AliFlowEventSimpleMakerOnTheFly::Init()
{
 // Book all objects in this method.
 
 // a) Define the pt spectra;
 // b) Define the phi distribution.

 // a) Define the pt spectra:
 Double_t dPtMin = 0.; 
 Double_t dPtMax = 10.; 
 fPtSpectra = new TF1("fPtSpectra","x*TMath::Exp(-pow([0]*[0]+x*x,0.5)/[1])",dPtMin,dPtMax); // hardwired is Boltzmann distribution  
 fPtSpectra->SetParName(0,"Mass");
 fPtSpectra->SetParameter(0,fMass);
 fPtSpectra->SetParName(1,"Temperature");
 fPtSpectra->SetParameter(1,fTemperature);
 fPtSpectra->SetTitle("Boltzmann Distribution: f(p_{t}) = p_{t}exp[-(m^{2}+p_{t}^{2})^{1/2}/T];p_{t};f(p_{t})");
 
 // b) Define the phi distribution:
 Double_t dPhiMin = 0.; 
 Double_t dPhiMax = TMath::TwoPi();
 fPhiDistribution = new TF1("fPhiDistribution","1+2.*[1]*TMath::Cos(x-[0])+2.*[2]*TMath::Cos(2.*(x-[0]))+2.*[3]*TMath::Cos(3.*(x-[0]))+2.*[4]*TMath::Cos(4.*(x-[0]))",dPhiMin,dPhiMax);
 fPhiDistribution->SetParName(0,"Reaction Plane");
 fPhiDistribution->SetParameter(0,0.);
 fPhiDistribution->SetParName(1,"Directed Flow (v1)"); 
 fPhiDistribution->SetParameter(1,fV1);
 fPhiDistribution->SetParName(2,"Elliptic Flow (v2)");
 fPhiDistribution->SetParameter(2,fV2);
 fPhiDistribution->SetParName(3,"Triangular Flow (v3)");
 fPhiDistribution->SetParameter(3,fV3);
 fPhiDistribution->SetParName(4,"Quadrangular Flow (v4)");
 fPhiDistribution->SetParameter(4,fV4);
    
} // end of void AliFlowEventSimpleMakerOnTheFly::Init()

//====================================================================================================================

Bool_t AliFlowEventSimpleMakerOnTheFly::AcceptOrNot(AliFlowTrackSimple *pTrack)
{
 // For the case of non-uniform acceptance determine in this method if particle is accepted or rejected.
 
 Bool_t bAccept = kTRUE;
 
 if((pTrack->Phi() >= fPhiMin1*fPi/180.) && (pTrack->Phi() < fPhiMax1*fPi/180.) && gRandom->Uniform(0,1) > fProbability1) 
 {
  bAccept = kFALSE; // particle is rejected in the first non-uniform sector
 } else if((pTrack->Phi() >= fPhiMin2*fPi/180.) && (pTrack->Phi() < fPhiMax2*fPi/180.) && gRandom->Uniform(0,1) > fProbability2) 
    {
     bAccept = kFALSE; // particle is rejected in the second non-uniform sector
    } 
  
 return bAccept;
 
} // end of Bool_t AliFlowEventSimpleMakerOnTheFly::AcceptOrNot(AliFlowTrackSimple *pTrack);

//====================================================================================================================

AliFlowEventSimple* AliFlowEventSimpleMakerOnTheFly::CreateEventOnTheFly(AliFlowTrackSimpleCuts *cutsRP, AliFlowTrackSimpleCuts *cutsPOI)
{
 // Method to create event 'on the fly'.
 
 // a) Determine the multiplicity of an event;
 // b) Determine the reaction plane of an event;
 // c) If v2 fluctuates uniformly event-by-event, sample its value from [fMinV2,fMaxV2];
 // d) Create event 'on the fly';
 // e) Cosmetics for the printout on the screen.
 
 // a) Determine the multiplicity of an event:
 Int_t iMult = (Int_t)gRandom->Uniform(fMinMult,fMaxMult);
 
 // b) Determine the reaction plane of an event:
 Double_t dReactionPlane = gRandom->Uniform(0.,TMath::TwoPi());
 fPhiDistribution->SetParameter(0,dReactionPlane);

 // c) If v2 fluctuates uniformly event-by-event, sample its value from [fMinV2,fMaxV2]:
 if(fUniformFluctuationsV2)
 {
  fPhiDistribution->SetParameter(2,gRandom->Uniform(fMinV2,fMaxV2));
 } 

 // d) Create event 'on the fly':
 AliFlowEventSimple *pEvent = new AliFlowEventSimple(iMult); 
 pEvent->SetReferenceMultiplicity(iMult);
 pEvent->SetMCReactionPlaneAngle(dReactionPlane); 
 Int_t nRPs = 0; // number of particles tagged RP in this event
 Int_t nPOIs = 0; // number of particles tagged POI in this event
 for(Int_t p=0;p<iMult;p++)
 {
  AliFlowTrackSimple *pTrack = new AliFlowTrackSimple();
  pTrack->SetPt(fPtSpectra->GetRandom()); 
  if(fPtDependentV2 && !fUniformFluctuationsV2)
  {
   // v2(pt): for pt < fV2vsPtCutOff v2 increases linearly, for pt >= fV2vsPtCutOff v2 = fV2vsPtMax
   (pTrack->Pt() < fV2vsPtCutOff ? 
    fPhiDistribution->SetParameter(2,pTrack->Pt()*fV2vsPtMax/fV2vsPtCutOff) :
    fPhiDistribution->SetParameter(2,fV2vsPtMax)
   );
  } // end of if(fPtDependentV2)  
  pTrack->SetPhi(fPhiDistribution->GetRandom());
  pTrack->SetEta(gRandom->Uniform(-1.,1.));
  pTrack->SetCharge((gRandom->Integer(2)>0.5 ? 1 : -1));
  // Check uniform acceptance:
  if(!fUniformAcceptance && !AcceptOrNot(pTrack)){continue;}
  // Checking the RP cuts:  	 
  if(cutsRP->PassesCuts(pTrack))
  {
   pTrack->TagRP(kTRUE); 
   nRPs++; 
  }
  // Checking the POI cuts:  	 
  if(cutsPOI->PassesCuts(pTrack))
  {
   pTrack->TagPOI(kTRUE); 
   nPOIs++;
  }
  // Assign particles to eta subevents (needed only for Scalar Product method):
  if(pTrack->Eta()>=fEtaMinA && pTrack->Eta()<fEtaMaxA) 
  {
   pTrack->SetForSubevent(0);
  }
  if(pTrack->Eta()>=fEtaMinB && pTrack->Eta()<fEtaMaxB) 
  {
   pTrack->SetForSubevent(1);
  }  
  pEvent->AddTrack(pTrack);
  // Simulating nonflow:
  if(fNTimes>1)
  {
   for(Int_t nt=1;nt<fNTimes;nt++)
   {
    pEvent->AddTrack(pTrack->Clone());  
   } 
  } // end of if(fNTimes>1)       
 } // end of for(Int_t p=0;p<iMult;p++)
 pEvent->SetNumberOfRPs(fNTimes*nRPs);
 
 // e) Cosmetics for the printout on the screen:
 Int_t cycle = (fPtDependentV2 ? 10 : 100);
 if((++fCount % cycle) == 0) 
 {
  if(TMath::Abs(dReactionPlane)>1.e-44) 
  {
   cout<<" MC Reaction Plane Angle = "<<dReactionPlane<<endl;
  } else 
    {
     cout<<" MC Reaction Plane Angle is unknown :'( "<< endl;
    }     
  cout<<" # of simulated tracks  = "<<fNTimes*iMult<<endl;
  cout<<" # of RP tagged tracks  = "<<fNTimes*nRPs<<endl;
  cout<<" # of POI tagged tracks = "<<fNTimes*nPOIs<<endl;  
  cout <<"  .... "<<fCount<< " events processed ...."<<endl;
 } // end of if((++fCount % cycle) == 0) 

 return pEvent;
    
} // end of CreateEventOnTheFly()
 
//====================================================================================================================

