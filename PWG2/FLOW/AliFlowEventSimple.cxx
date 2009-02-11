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
  fTrackCollection->Delete() ; delete fTrackCollection ;
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
  //calculate Q-vector in harmonic n without weights (default harmonic n=2)  
  Double_t dQX = 0.;
  Double_t dQY = 0.;
  AliFlowVector vQ;
  vQ.Set(0.,0.);
  
  Int_t iOrder = n;
  Int_t iUsedTracks = 0;
  Double_t dPhi=0.;
  Double_t dPt=0;
  Double_t dEta=0;
  
  AliFlowTrackSimple* pTrack = NULL;
 
  Int_t nBinsPhi=0; 
  Double_t dBinWidthPt=0.;
  Double_t dNormPt=0.;
  Double_t dBinWidthEta=0.;
  Double_t dNormEta=0.;
 
  Double_t wPhi=1.; //weight Phi  
  Double_t wPt=1.;  //weight Pt 
  Double_t wEta=1.; //weight Eta 
  
  TH1F *phiWeights = NULL;
  TH1D *ptWeights = NULL;
  TH1D *etaWeights = NULL;

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
     dNormPt = ptWeights->Integral();
    } 
   }       
   if(useEtaWeights)
   {
    etaWeights = dynamic_cast<TH1D *>(weightsList->FindObject("eta_weights"));
    if(etaWeights)
    {
     dBinWidthEta = etaWeights->GetBinWidth(1); // assuming that all bins have the same width
     dNormEta = etaWeights->Integral();
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
    
     //determine Phi weight: (to be improved, I should here only access it + the treatment of gaps in the if statement)
     if(phiWeights && (nBinsPhi!=0) && (phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(dPhi*nBinsPhi/TMath::TwoPi())))!=0))
     {
      wPhi=pow(nBinsPhi*phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(dPhi*nBinsPhi/TMath::TwoPi()))),-1);
     }
     //determine v'(pt) weight:    
     if(ptWeights && dBinWidthPt && dNormPt)
     {
      wPt=ptWeights->GetBinContent(1+(Int_t)(TMath::Floor(dPt/dBinWidthPt)))/dNormPt; 
     }            
     //determine v'(eta) weight:    
     if(etaWeights && dBinWidthEta && dNormEta)
     {
      wEta=etaWeights->GetBinContent(1+(Int_t)(TMath::Floor(dEta/dBinWidthEta)))/dNormEta; 
     } 
       
     //building up the weighted Q-vector:       
     dQX += wPhi*wPt*wEta*TMath::Cos(iOrder*dPhi);
     dQY += wPhi*wPt*wEta*TMath::Sin(iOrder*dPhi); 
     iUsedTracks++;
     
    }//end of if (pTrack->UseForIntegratedFlow())
   }//end of if (pTrack)
   else {cerr << "no particle!!!"<<endl;}
  }//loop over particles
    
  vQ.Set(dQX,dQY);
  vQ.SetMult(iUsedTracks);

  return vQ;
  
}










