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
  fPtFormula(NULL),
  fPhiFormula(NULL),
  fMyTRandom3(NULL),
  fCount(0)
 {
  // constructor
   fMyTRandom3 = new TRandom3(iseed); 
 }


//========================================================================


AliFlowEventSimpleMakerOnTheFly::~AliFlowEventSimpleMakerOnTheFly()
{
 // destructor
  if (fPtFormula) delete fPtFormula;
  if (fPhiFormula) delete fPhiFormula;
  if (fMyTRandom3) delete  fMyTRandom3;
}


//========================================================================


AliFlowEventSimple* AliFlowEventSimpleMakerOnTheFly::CreateEventOnTheFly()
{
  // method to create event on the fly
  
  AliFlowEventSimple* pEvent = new AliFlowEventSimple(fMultiplicityOfRP);
  
  //reaction plane
  Double_t fMCReactionPlaneAngle = fMyTRandom3->Uniform(0.,TMath::TwoPi());
  
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


  Int_t iGoodTracks = 0;
  Int_t iSelParticlesRP = 0;
  Int_t iSelParticlesPOI = 0;
  for(Int_t i=0;i<fMultiplicityOfRP;i++) {
    AliFlowTrackSimple* pTrack = new AliFlowTrackSimple();
    pTrack->SetPt(fPtFormula->GetRandom());
    pTrack->SetEta(fMyTRandom3->Uniform(dEtaMin,dEtaMax));
    pTrack->SetPhi(fPhiFormula->GetRandom()+fMCReactionPlaneAngle);
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

  if (!fMCReactionPlaneAngle == 0) cout<<" MC Reaction Plane Angle = "<<  fMCReactionPlaneAngle << endl;
  else cout<<" MC Reaction Plane Angle = unknown "<< endl;

  cout<<" iGoodTracks = "<< iGoodTracks << endl;
  cout<<" # of RP selected tracks = "<<iSelParticlesRP<<endl;
  cout<<" # of POI selected tracks = "<<iSelParticlesPOI<<endl;  
  cout << "# " << ++fCount << " events processed" << endl;

 return pEvent;  
 
} // end of CreateEventOnTheFly()


 
