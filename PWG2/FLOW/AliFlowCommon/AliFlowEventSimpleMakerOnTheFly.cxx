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


AliFlowEventSimpleMakerOnTheFly::AliFlowEventSimpleMakerOnTheFly():
 fMultiplicityOfRP(0),
 fMultiplicitySpreadOfRP(0.),
 fPtFormula(NULL),
 fPhiFormula(NULL) 
 {
  // constructor
 }


//========================================================================


AliFlowEventSimpleMakerOnTheFly::~AliFlowEventSimpleMakerOnTheFly()
{
 // destructor
  if (fPtFormula) delete fPtFormula;
  if (fPhiFormula) delete fPhiFormula;
}


//========================================================================


AliFlowEventSimple* AliFlowEventSimpleMakerOnTheFly::CreateEventOnTheFly()
{
  // method to create event on the fly
  
  AliFlowEventSimple* pEvent = new AliFlowEventSimple(fMultiplicityOfRP);
  
  
  // pt:   
  Double_t dPtMin = 0.; // to be improved 
  Double_t dPtMax = 10.; // to be improved 
  
  fPtFormula = new TF1("PtFormula","[0]*x*TMath::Exp(-x*x)",dPtMin,dPtMax);
  
  fPtFormula->SetParName(0,"Multiplicity of RPs");
  fPtFormula->SetParameter(0,fMultiplicityOfRP);
  

  // phi:
  Double_t dPhiMin = 0.; // to be improved 
  Double_t dPhiMax = TMath::TwoPi(); // to be improved 
  
  fPhiFormula = new TF1("phiDistribution","(1)*(1+2.*[0]*TMath::Cos(2*x))",dPhiMin,dPhiMax);
 
  Double_t dV2 = 0.044; // to be improved
 
  fPhiFormula->SetParName(0,"elliptic flow");
  fPhiFormula->SetParameter(0,dV2);
 
  
  // eta:
  Double_t dEtaMin = -1.; // to be improved 
  Double_t dEtaMax = 1.; // to be improved 
  TRandom3 myTRandom3;

  //reaction plane
  Double_t fMCReactionPlaneAngle = TMath::TwoPi()*myTRandom3.Rndm();

  Int_t iGoodTracks = 0;
  Int_t iSelParticlesRP = 0;
  Int_t iSelParticlesPOI = 0;
  for(Int_t i=0;i<fMultiplicityOfRP;i++) {
    AliFlowTrackSimple* pTrack = new AliFlowTrackSimple();
    pTrack->SetPt(fPtFormula->GetRandom());
    pTrack->SetEta(myTRandom3.Uniform(dEtaMin,dEtaMax));
    pTrack->SetPhi(fPhiFormula->GetRandom()-fMCReactionPlaneAngle);
    pTrack->SetForRPSelection(kTRUE);
    iSelParticlesRP++;
    pTrack->SetForPOISelection(kTRUE);
    iSelParticlesPOI++;

    pEvent->TrackCollection()->Add(pTrack);
    iGoodTracks++;
  }
 
  pEvent->SetEventNSelTracksRP(iSelParticlesRP);  
  pEvent->SetNumberOfTracks(iGoodTracks);//tracks used either for RP or for POI selection
  pEvent->SetMCReactionPlaneAngle(fMCReactionPlaneAngle);


 return pEvent;  
 
} // end of CreateEventOnTheFly()


 
