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

#include "Riostream.h"
#include "TObjArray.h"
#include "TFile.h"
#include "TList.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "AliFlowVector.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowEventSimple.h"

/**************************************
 * AliFlowEventSimple: A simple event *
 *  for flow analysis                 * 
 *                                    * 
 * authors: Naomi van der Kolk        *
 *           (kolk@nikhef.nl)         *  
 *          Ante Bilandzic            *
 *           (anteb@nikhef.nl)        * 
 * ***********************************/
 
ClassImp(AliFlowEventSimple)

//-----------------------------------------------------------------------

AliFlowEventSimple::AliFlowEventSimple():
  fTrackCollection(NULL),
  fNumberOfTracks(0),
  fEventNSelTracksIntFlow(0)
{
  cout << "AliFlowEventSimple: Default constructor to be used only by root for io" << endl;
}

//-----------------------------------------------------------------------

AliFlowEventSimple::AliFlowEventSimple(Int_t aLenght):
  fTrackCollection(NULL),
  fNumberOfTracks(0),
  fEventNSelTracksIntFlow(0)
{
  //constructor 
  fTrackCollection =  new TObjArray(aLenght) ;
}

//-----------------------------------------------------------------------

AliFlowEventSimple::AliFlowEventSimple(const AliFlowEventSimple& anEvent):
  TObject(),
  fTrackCollection(anEvent.fTrackCollection),
  fNumberOfTracks(anEvent.fNumberOfTracks),
  fEventNSelTracksIntFlow(anEvent.fEventNSelTracksIntFlow)
{
  //copy constructor 
}

//-----------------------------------------------------------------------

AliFlowEventSimple& AliFlowEventSimple::operator=(const AliFlowEventSimple& anEvent)
{
  *fTrackCollection = *anEvent.fTrackCollection ;
  fNumberOfTracks = anEvent.fNumberOfTracks;
  fEventNSelTracksIntFlow = anEvent.fEventNSelTracksIntFlow;
  
  return *this;
}

//----------------------------------------------------------------------- 

AliFlowEventSimple::~AliFlowEventSimple()
{
  //destructor
  if (fTrackCollection) {
    fTrackCollection->Delete() ; delete fTrackCollection ;
  }
  else { 
    cout << "AliFlowEventSimple: Warning trying to delete track collections NULL pointer" << endl; 
  }
}

//----------------------------------------------------------------------- 

AliFlowTrackSimple* AliFlowEventSimple::GetTrack(Int_t i)
{
  //get track i from collection
  AliFlowTrackSimple* pTrack = (AliFlowTrackSimple*)TrackCollection()->At(i) ;
  return pTrack;
}

//-----------------------------------------------------------------------   
AliFlowVector AliFlowEventSimple::GetQ(Int_t n, TList *weightsList, Bool_t usePhiWeights, Bool_t usePtWeights, Bool_t useEtaWeights) 
{
  // calculate Q-vector in harmonic n without weights (default harmonic n=2)  
  Double_t dQX = 0.;
  Double_t dQY = 0.;
  AliFlowVector vQ;
  vQ.Set(0.,0.);
  
  Int_t iOrder = n;
  Double_t iUsedTracks = 0;
  Double_t dPhi=0.;
  Double_t dPt=0.;
  Double_t dEta=0.;
  
  AliFlowTrackSimple* pTrack = NULL;
 
  Int_t nBinsPhi=0; 
  Double_t dBinWidthPt=0.;
  Double_t dPtMin=0.;
  Double_t dBinWidthEta=0.;
  Double_t dEtaMin=0.;
 
  Double_t wPhi=1.; // weight Phi  
  Double_t wPt=1.;  // weight Pt 
  Double_t wEta=1.; // weight Eta 
  
  TH1F *phiWeights = NULL;
  TH1D *ptWeights  = NULL;
  TH1D *etaWeights = NULL;

  Double_t dSumOfWeightsToPower2 = 0.; // sum_{i=1}^{n} pow((wPhi*wPt*wEta)_i, 2)
  Double_t dSumOfWeightsToPower3 = 0.; // sum_{i=1}^{n} pow((wPhi*wPt*wEta)_i, 3)
  Double_t dSumOfWeightsToPower4 = 0.; // sum_{i=1}^{n} pow((wPhi*wPt*wEta)_i, 4)
  Double_t dSumOfWeightsToPower5 = 0.; // sum_{i=1}^{n} pow((wPhi*wPt*wEta)_i, 5)
  Double_t dSumOfWeightsToPower6 = 0.; // sum_{i=1}^{n} pow((wPhi*wPt*wEta)_i, 6)
  Double_t dSumOfWeightsToPower7 = 0.; // sum_{i=1}^{n} pow((wPhi*wPt*wEta)_i, 7)
  Double_t dSumOfWeightsToPower8 = 0.; // sum_{i=1}^{n} pow((wPhi*wPt*wEta)_i, 8) 

  if(weightsList)
  {
   if(usePhiWeights)
   {
    phiWeights = dynamic_cast<TH1F *>(weightsList->FindObject("phi_weights"));
    if(phiWeights) nBinsPhi = phiWeights->GetNbinsX();
   }          
   if(usePtWeights)
   {
    ptWeights = dynamic_cast<TH1D *>(weightsList->FindObject("pt_weights"));
    if(ptWeights)
    {
     dBinWidthPt = ptWeights->GetBinWidth(1); // assuming that all bins have the same width
     dPtMin = (ptWeights->GetXaxis())->GetXmin();
    } 
   }       
   if(useEtaWeights)
   {
    etaWeights = dynamic_cast<TH1D *>(weightsList->FindObject("eta_weights"));
    if(etaWeights)
    {
     dBinWidthEta = etaWeights->GetBinWidth(1); // assuming that all bins have the same width
     dEtaMin = (etaWeights->GetXaxis())->GetXmin();
    } 
   }          
  } // end of if(weightsList)
  
  // loop over tracks    
  for(Int_t i=0;i<fNumberOfTracks;i++)                               
  {
   pTrack = (AliFlowTrackSimple*)TrackCollection()->At(i); 
   if(pTrack)
   {
    if(pTrack->UseForIntegratedFlow()) 
    {
     dPhi = pTrack->Phi();
     dPt  = pTrack->Pt();
     dEta = pTrack->Eta();
     
     // determine Phi weight: (to be improved, I should here only access it + the treatment of gaps in the if statement)
     if(phiWeights && nBinsPhi)
     {
      wPhi = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(dPhi*nBinsPhi/TMath::TwoPi())));
     }
     // determine v'(pt) weight:    
     if(ptWeights && dBinWidthPt)
     {
      wPt=ptWeights->GetBinContent(1+(Int_t)(TMath::Floor((dPt-dPtMin)/dBinWidthPt))); 
     }            
     // determine v'(eta) weight:    
     if(etaWeights && dBinWidthEta)
     {
      wEta=etaWeights->GetBinContent(1+(Int_t)(TMath::Floor((dEta-dEtaMin)/dBinWidthEta))); 
     } 

     // building up the weighted Q-vector:       
     dQX += wPhi*wPt*wEta*TMath::Cos(iOrder*dPhi);
     dQY += wPhi*wPt*wEta*TMath::Sin(iOrder*dPhi); 
    
     // weighted multiplicity:
     iUsedTracks+=wPhi*wPt*wEta;
    
     // weights raised to various powers are summed up:
     dSumOfWeightsToPower2+=pow(wPhi*wPt*wEta, 2); 
     dSumOfWeightsToPower3+=pow(wPhi*wPt*wEta, 3); 
     dSumOfWeightsToPower4+=pow(wPhi*wPt*wEta, 4); 
     dSumOfWeightsToPower5+=pow(wPhi*wPt*wEta, 5); 
     dSumOfWeightsToPower6+=pow(wPhi*wPt*wEta, 6); 
     dSumOfWeightsToPower7+=pow(wPhi*wPt*wEta, 7); 
     dSumOfWeightsToPower8+=pow(wPhi*wPt*wEta, 8); 
     
    } // end of if (pTrack->UseForIntegratedFlow())
   } // end of if (pTrack)
   else {cerr << "no particle!!!"<<endl;}
  } // loop over particles
    
  vQ.Set(dQX,dQY);
  vQ.SetMult(iUsedTracks);
  vQ.SetSumOfWeightsToPower2(dSumOfWeightsToPower2);
  vQ.SetSumOfWeightsToPower3(dSumOfWeightsToPower3);
  vQ.SetSumOfWeightsToPower4(dSumOfWeightsToPower4);
  vQ.SetSumOfWeightsToPower5(dSumOfWeightsToPower5);
  vQ.SetSumOfWeightsToPower6(dSumOfWeightsToPower6);
  vQ.SetSumOfWeightsToPower7(dSumOfWeightsToPower7);
  vQ.SetSumOfWeightsToPower8(dSumOfWeightsToPower8);

  return vQ;
  
}










