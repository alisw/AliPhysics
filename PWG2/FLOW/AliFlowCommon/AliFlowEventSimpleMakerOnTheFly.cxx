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

/********************************** 
 * create an event and perform    *
 * flow analysis 'on the fly'     * 
 *                                * 
 * authors: Raimond Snellings     *
 *           (snelling@nikhef.nl) * 
 *          Ante Bilandzic        * 
 *           (anteb@nikhef.nl)    *
 *********************************/
  
#include "Riostream.h"
#include "TMath.h"
#include "TF1.h"
#include "TRandom3.h"

#include "AliFlowEventSimpleMakerOnTheFly.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"

ClassImp(AliFlowEventSimpleMakerOnTheFly)


//========================================================================


AliFlowEventSimpleMakerOnTheFly::AliFlowEventSimpleMakerOnTheFly(UInt_t iseed):
  fMultiplicityOfRP(0),
  fMultiplicitySpreadOfRP(0.),
  fTemperatureOfRP(0.),  
  fUseConstantHarmonics(kFALSE),
  fV1RP(0.), 
  fV1SpreadRP(0.), 
  fV2RP(0.), 
  fV2SpreadRP(0.), 
  fV4RP(0.), 
  fV4SpreadRP(0.), 
  fV2RPMax(0.), 
  fPtCutOff(0.), 
  fPtSpectra(NULL),
  fPhiDistribution(NULL),
  fMyTRandom3(NULL),
  fCount(0),
  fNoOfLoops(1)
 {
  // constructor
  fMyTRandom3 = new TRandom3(iseed);   
  gRandom->SetSeed(fMyTRandom3->Integer(65539));
 }


//========================================================================


AliFlowEventSimpleMakerOnTheFly::~AliFlowEventSimpleMakerOnTheFly()
{
 // destructor
 if (fPtSpectra) delete fPtSpectra;
 if (fPhiDistribution) delete fPhiDistribution;
 if (fMyTRandom3) delete fMyTRandom3;
}


//========================================================================


void AliFlowEventSimpleMakerOnTheFly::Init()
{
 // define the pt spectra and phi distribution

 // pt spectra of pions (Boltzman):   
 Double_t dPtMin = 0.; // to be improved (move this to the body of contstructor?)
 Double_t dPtMax = 10.; // to be improved (move this to the body of contstructor?) 
  
 fPtSpectra = new TF1("fPtSpectra","[0]*x*TMath::Exp(-pow(0.13957*0.13957+x*x,0.5)/[1])",dPtMin,dPtMax);  
 fPtSpectra->SetParName(0,"Multiplicity of RPs");  
 fPtSpectra->SetParName(1,"Temperature of RPs");
 
 // phi distribution:
 Double_t dPhiMin = 0.; // to be improved (move this to the body of contstructor?)
 Double_t dPhiMax = TMath::TwoPi(); // to be improved (move this to the body of contstructor?)
  
 fPhiDistribution = new TF1("fPhiDistribution","1+2.*[0]*TMath::Cos(x-[2])+2.*[1]*TMath::Cos(2*(x-[2]))+2.*[3]*TMath::Cos(4*(x-[2]))",dPhiMin,dPhiMax);
 fPhiDistribution->SetParName(0,"directed flow");
 fPhiDistribution->SetParName(1,"elliptic flow"); 
 fPhiDistribution->SetParName(2,"Reaction Plane");
 fPhiDistribution->SetParName(3,"harmonic 4"); // to be improved (name)
}

//========================================================================

AliFlowEventSimple* AliFlowEventSimpleMakerOnTheFly::CreateEventOnTheFly()
{
  // method to create event on the fly
  
  AliFlowEventSimple* pEvent = new AliFlowEventSimple(fMultiplicityOfRP);
   
  // sampling the multiplicity:
  Int_t iNewMultiplicityOfRP = fMultiplicityOfRP;
  if(fMultiplicitySpreadOfRP>0.0) iNewMultiplicityOfRP = (Int_t)fMyTRandom3->Gaus(fMultiplicityOfRP,fMultiplicitySpreadOfRP);
  fPtSpectra->SetParameter(0,iNewMultiplicityOfRP);
  
  // set the 'temperature' of RPs
  fPtSpectra->SetParameter(1,fTemperatureOfRP);  
  
  // sampling the reaction plane
  Double_t dMCReactionPlaneAngle = fMyTRandom3->Uniform(0.,TMath::TwoPi());
  fPhiDistribution->SetParameter(2,dMCReactionPlaneAngle);

  // sampling the V1:
  Double_t dNewV1RP=fV1RP;
  if(fV1SpreadRP>0.0) {dNewV1RP = fMyTRandom3->Gaus(fV1RP,fV1SpreadRP);}
  fPhiDistribution->SetParameter(0,dNewV1RP);
 
  // sampling the V2:
  if(fUseConstantHarmonics)
  {
   Double_t dNewV2RP = fV2RP;
   if(fV2SpreadRP>0.0) dNewV2RP = fMyTRandom3->Gaus(fV2RP,fV2SpreadRP);
   fPhiDistribution->SetParameter(1,dNewV2RP);
  }
  
  // sampling the V4:
  Double_t dNewV4RP = fV4RP;
  if(fV4SpreadRP>0.0) dNewV4RP = fMyTRandom3->Gaus(fV4RP,fV4SpreadRP);
  fPhiDistribution->SetParameter(3,dNewV4RP);
   
  // eta:
  Double_t dEtaMin = -1.; // to be improved 
  Double_t dEtaMax = 1.; // to be improved 
  
  Int_t iGoodTracks = 0;
  Int_t iSelParticlesRP = 0;
  Int_t iSelParticlesPOI = 0;
  Double_t dTmpPt = 0.;
  Double_t dTmpEta = 0.;
  Double_t dTmpV2 = 0.;
  Double_t dTmpPhi = 0.;
  for(Int_t i=0;i<iNewMultiplicityOfRP;i++) {
    dTmpPt = fPtSpectra->GetRandom();
    dTmpEta = fMyTRandom3->Uniform(dEtaMin,dEtaMax);
    // to be improved:
    if(!fUseConstantHarmonics) {
      if(dTmpPt >= fPtCutOff) { 
	dTmpV2 = fV2RPMax;
      } else {
	dTmpV2 = fV2RPMax*(dTmpPt/fPtCutOff);
      }  
      fPhiDistribution->SetParameter(1,dTmpV2);         
    }
    dTmpPhi = fPhiDistribution->GetRandom();
    // add the track to the event
    for(Int_t d=0;d<fNoOfLoops;d++) {
      AliFlowTrackSimple* pTrack = new AliFlowTrackSimple();
      pTrack->SetPt(dTmpPt);
      pTrack->SetEta(dTmpEta);
      pTrack->SetPhi(dTmpPhi);
      pTrack->SetForRPSelection(kTRUE);
      iSelParticlesRP++;
      pTrack->SetForPOISelection(kTRUE);
      iSelParticlesPOI++;
      pEvent->TrackCollection()->Add(pTrack);
      iGoodTracks++;
    }
  }
  
  // update the event quantities
  pEvent->SetEventNSelTracksRP(iSelParticlesRP);  
  pEvent->SetNumberOfTracks(iGoodTracks);//tracks used either for RP or for POI selection
  pEvent->SetMCReactionPlaneAngle(dMCReactionPlaneAngle);

  if (!dMCReactionPlaneAngle == 0) cout<<" MC Reaction Plane Angle = "<<  dMCReactionPlaneAngle << endl;
  else cout<<" MC Reaction Plane Angle = unknown "<< endl;

  cout<<" iGoodTracks = "<< iGoodTracks << endl;
  cout<<" # of RP selected tracks = "<<iSelParticlesRP<<endl;
  cout<<" # of POI selected tracks = "<<iSelParticlesPOI<<endl;  
  cout << "# " << ++fCount << " events processed" << endl;
  
  return pEvent;  
 
} // end of CreateEventOnTheFly()


 
