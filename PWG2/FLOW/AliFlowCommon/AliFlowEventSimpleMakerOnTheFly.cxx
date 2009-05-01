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
  fV1RP(0.), 
  fV1SpreadRP(0.), 
  fV2RP(0.), 
  fV2SpreadRP(0.), 
  fPtSpectra(NULL),
  fPhiDistribution(NULL),
  fMyTRandom3(NULL),
  fCount(0),
  fNoOfLoops(1)
 {
  // constructor
   fMyTRandom3 = new TRandom3(iseed); 
 }


//========================================================================


AliFlowEventSimpleMakerOnTheFly::~AliFlowEventSimpleMakerOnTheFly()
{
 // destructor
  if (fPtSpectra) delete fPtSpectra;
  if (fPhiDistribution) delete fPhiDistribution;
  if (fMyTRandom3) delete  fMyTRandom3;
}


//========================================================================


AliFlowEventSimple* AliFlowEventSimpleMakerOnTheFly::CreateEventOnTheFly()
{
  // method to create event on the fly
  
  AliFlowEventSimple* pEvent = new AliFlowEventSimple(fMultiplicityOfRP);
    
  // pt:   
  Double_t dPtMin = 0.; // to be improved 
  Double_t dPtMax = 10.; // to be improved 
  
  fPtSpectra = new TF1("fPtSpectra","[0]*x*TMath::Exp(-x*x)",dPtMin,dPtMax);  
  fPtSpectra->SetParName(0,"Multiplicity of RPs");  
  // sampling the multiplicity:
  Int_t fNewMultiplicityOfRP = fMyTRandom3->Gaus(fMultiplicityOfRP,fMultiplicitySpreadOfRP);
  fPtSpectra->SetParameter(0,fNewMultiplicityOfRP);
  
  
  // phi:
  Double_t dPhiMin = 0.; // to be improved 
  Double_t dPhiMax = TMath::TwoPi(); // to be improved 
  
  fPhiDistribution = new TF1("fPhiDistribution","1+2.*[0]*TMath::Cos(x-[2])+2.*[1]*TMath::Cos(2*(x-[2]))",dPhiMin,dPhiMax);
  
  // samling the reaction plane
  Double_t fMCReactionPlaneAngle = fMyTRandom3->Uniform(0.,TMath::TwoPi());
  fPhiDistribution->SetParName(2,"Reaction Plane");
  fPhiDistribution->SetParameter(2,fMCReactionPlaneAngle);

  // sampling the V1:
  fPhiDistribution->SetParName(0,"directed flow");
  Double_t fNewV1RP=fV1RP;
  if(fV1SpreadRP>0.0) {fNewV1RP = fMyTRandom3->Gaus(fV1RP,fV1SpreadRP);}
  fPhiDistribution->SetParameter(0,fNewV1RP);
 
  // sampling the V2:
  fPhiDistribution->SetParName(1,"elliptic flow");
  Double_t fNewV2RP = fV2RP;
  if(fV2SpreadRP>0.0) fNewV2RP = fMyTRandom3->Gaus(fV2RP,fV2SpreadRP);
  fPhiDistribution->SetParameter(1,fNewV2RP);
   
  // eta:
  Double_t dEtaMin = -1.; // to be improved 
  Double_t dEtaMax = 1.; // to be improved 
  
  Int_t iGoodTracks = 0;
  Int_t iSelParticlesRP = 0;
  Int_t iSelParticlesPOI = 0;
  Double_t fTmpPt =0;
  Double_t fTmpEta =0;
  Double_t fTmpPhi =0;
  
  for(Int_t i=0;i<fNewMultiplicityOfRP;i++) {
    fTmpPt = fPtSpectra->GetRandom();
    fTmpEta = fMyTRandom3->Uniform(dEtaMin,dEtaMax);
    fTmpPhi = fPhiDistribution->GetRandom();
    for(Int_t d=0;d<fNoOfLoops;d++) {
      AliFlowTrackSimple* pTrack = new AliFlowTrackSimple();
      pTrack->SetPt(fTmpPt);
      pTrack->SetEta(fTmpEta);
      pTrack->SetPhi(fTmpPhi);
      pTrack->SetForRPSelection(kTRUE);
      iSelParticlesRP++;
      pTrack->SetForPOISelection(kTRUE);
      iSelParticlesPOI++;
      pEvent->TrackCollection()->Add(pTrack);
      iGoodTracks++;
    }
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

  delete fPhiDistribution;
  delete fPtSpectra;
  return pEvent;  
 
} // end of CreateEventOnTheFly()


 
