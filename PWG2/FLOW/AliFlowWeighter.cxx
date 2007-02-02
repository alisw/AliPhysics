	//////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Author: Emanuele Simili
//
//////////////////////////////////////////////////////////////////////
//_____________________________________________________________
//
// Description: 
//         the AliFlowWeighter class generates the phi-weights which 
// are later used to correct for azimuthal anisotropy in the reconstruction 
// efficiency of the ALICE TPC. It also fills an histogram of normalised 
// particle abundancies, which can be used as bayesian weights for particle Id. 
//
// - The method Init() generates the histograms and opens a new fPhiWgt file.
// - The method WeightEvent(AliFlowEvent*) fills phi and PId histograms. 
//   It must be inserted in a loop over the event sample.
// - The method Finish() calculates the weights, saves the histograms and 
//   closes the file. The AliFlowSelection object is saved as well.
//
//////////////////////////////////////////////////////////////////////

#ifndef ALIFLOWANALYSER_CXX
#define ALIFLOWANALYSER_CXX

// ROOT things
#include <TROOT.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TObject.h>
#include <TObjArray.h>
#include <TOrdCollection.h>
#include <TH1.h>
#include <TVector2.h>

// Flow things
#include "AliFlowEvent.h"
#include "AliFlowTrack.h"
#include "AliFlowV0.h"
#include "AliFlowConstants.h"
#include "AliFlowSelection.h"
#include "AliFlowWeighter.h"

// ANSI things
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std; //required for resolving the 'cout' symbol

ClassImp(AliFlowWeighter) ;
//-----------------------------------------------------------------------
AliFlowWeighter::AliFlowWeighter(const AliFlowSelection* flowSelect)
{
 // default constructor (selection given or default selection)

 if(flowSelect) { fFlowSelect = new AliFlowSelection(*flowSelect) ; }
 else 		{ fFlowSelect = new AliFlowSelection() ; }
 
 // output file (histograms)
 fWgtFileName = "flowPhiWgt.hist.root" ;
 fWgtFile     = 0 ;
}
//-----------------------------------------------------------------------
AliFlowWeighter::~AliFlowWeighter()
{
 // default distructor (no actions) 
}
//-----------------------------------------------------------------------
Bool_t AliFlowWeighter::Init() 
{
// sets some defaults for the analysis
 
 cout << "* FlowWeighter *  -  Init()" << endl ; cout << endl ; 
 
 // Open output files (->plots)
 fWgtFile = new TFile(fWgtFileName.Data(), "RECREATE");
 fWgtFile->cd() ;

 // counters and pointers to 0
 fEventNumber = 0 ;
 fTrackNumber = 0 ; 	
 fNumberOfTracks = 0 ;	
 fNumberOfV0s    = 0 ;	
 fFlowEvent  = 0 ;
 fFlowTrack  = 0 ;
 fFlowTracks = 0 ;
 //for(Int_t ii=0;ii<3;ii++) { fVertex[ii] = 0 ; }

 // Histogram settings
 fPhiBins = Flow::nPhiBins ;     
 fPhiMin  = 0.;
 fPhiMax  = 2*TMath::Pi() ; 

 TString* histTitle ;
 for(int k = 0; k < Flow::nSels; k++)
 {
  histTitle = new TString("Flow_BayPidMult_Sel");
  *histTitle += k+1;
  fHistFull[k].fHistBayPidMult = new TH1F(histTitle->Data(), histTitle->Data(),Flow::nPid,-0.5,((Float_t)Flow::nPid-0.5));
  fHistFull[k].fHistBayPidMult->Sumw2() ;
  fHistFull[k].fHistBayPidMult->SetXTitle("e+/-  ,  mu+/-  ,  pi+/-  ,  K+/-  ,  p+/-  ,  d+/- ");
  fHistFull[k].fHistBayPidMult->SetYTitle("Counts");
  delete histTitle;

  for(int j = 0; j < Flow::nHars; j++) 
  {
   // Tpc - Phi lab
   histTitle = new TString("Flow_Phi_TPC_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPhi = new TH1D(histTitle->Data(), histTitle->Data(), fPhiBins, fPhiMin, fPhiMax);
   fHistFull[k].fHistFullHar[j].fHistPhi->SetXTitle("Azimuthal Angles (rad)");
   fHistFull[k].fHistFullHar[j].fHistPhi->SetYTitle("Counts");
   delete histTitle;
   // Tpc (plus)
   histTitle = new TString("Flow_Phi_TPCplus_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPhiPlus = new TH1D(histTitle->Data(), histTitle->Data(), fPhiBins, fPhiMin, fPhiMax);
   fHistFull[k].fHistFullHar[j].fHistPhiPlus->SetXTitle("Azimuthal Angles (rad)");
   fHistFull[k].fHistFullHar[j].fHistPhiPlus->SetYTitle("Counts");
   delete histTitle;
   // Tpc (minus)
   histTitle = new TString("Flow_Phi_TPCminus_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPhiMinus = new TH1D(histTitle->Data(), histTitle->Data(), fPhiBins, fPhiMin, fPhiMax);
   fHistFull[k].fHistFullHar[j].fHistPhiMinus->SetXTitle("Azimuthal Angles (rad)");
   fHistFull[k].fHistFullHar[j].fHistPhiMinus->SetYTitle("Counts");
   delete histTitle;
   // Tpc (cross)
   histTitle = new TString("Flow_Phi_TPCcross_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPhiAll = new TH1D(histTitle->Data(), histTitle->Data(), fPhiBins, fPhiMin, fPhiMax);
   fHistFull[k].fHistFullHar[j].fHistPhiAll->SetXTitle("Azimuthal Angles (rad)");
   fHistFull[k].fHistFullHar[j].fHistPhiAll->SetYTitle("Counts");
   delete histTitle;

   // Tpc - Phi lab flattened
   histTitle = new TString("Flow_Phi_Flat_TPC_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPhiFlat = new TH1D(histTitle->Data(), histTitle->Data(), fPhiBins, fPhiMin, fPhiMax);
   fHistFull[k].fHistFullHar[j].fHistPhiFlat->SetXTitle("Azimuthal Angles (rad)");
   fHistFull[k].fHistFullHar[j].fHistPhiFlat->SetYTitle("Counts");
   delete histTitle;
   // Tpc Plus
   histTitle = new TString("Flow_Phi_Flat_TPCplus_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPhiFlatPlus = new TH1D(histTitle->Data(), histTitle->Data(), fPhiBins, fPhiMin, fPhiMax);
   fHistFull[k].fHistFullHar[j].fHistPhiFlatPlus->SetXTitle("Azimuthal Angles (rad)");
   fHistFull[k].fHistFullHar[j].fHistPhiFlatPlus->SetYTitle("Counts");
   delete histTitle;
   // Tpc Minus
   histTitle = new TString("Flow_Phi_Flat_TPCminus_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPhiFlatMinus = new TH1D(histTitle->Data(), histTitle->Data(), fPhiBins, fPhiMin, fPhiMax);
   fHistFull[k].fHistFullHar[j].fHistPhiFlatMinus->SetXTitle("Azimuthal Angles (rad)");
   fHistFull[k].fHistFullHar[j].fHistPhiFlatMinus->SetYTitle("Counts");
   delete histTitle;
   // Tpc cross
   histTitle = new TString("Flow_Phi_Flat_TPCcross_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPhiFlatAll = new TH1D(histTitle->Data(), histTitle->Data(), fPhiBins, fPhiMin, fPhiMax);
   fHistFull[k].fHistFullHar[j].fHistPhiFlatAll->SetXTitle("Azimuthal Angles (rad)");
   fHistFull[k].fHistFullHar[j].fHistPhiFlatAll->SetYTitle("Counts");
   delete histTitle;
  }
 }
 
 return kTRUE ;
} 
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Bool_t AliFlowWeighter::Finish() 
{
 // Calls the method to fill wgt histograms, then saves them 
 // on the fWgtFile and closes the file .
 
 cout << "* FlowWeighter *  -  Finish()" << endl ; cout << endl ;

 Weightening() ;
 fWgtFile->cd() ; 

 // Write PhiWgt histograms
 fPhiWgtHistList->Write();
 delete fPhiWgtHistList ;
 
 // Write Bayesian Weights for P.Id.
 for(int k=0;k<Flow::nSels;k++) { fHistFull[k].fHistBayPidMult->Write() ; }

 // Write the AliFlowSelection object
 fFlowSelect->Write();
 delete fFlowSelect ;

 fWgtFile->Close();

 cout << "    Finish()  -  Wgt file closed : " << fWgtFileName.Data() << endl ; cout << endl ; 

 return kTRUE ;
}
//-----------------------------------------------------------------------
// ###
//-----------------------------------------------------------------------
Bool_t AliFlowWeighter::WeightEvent(AliFlowEvent* fFlowEvent)         
{
 // Reads the AliFlowEvent (* fFlowEvent) and loops over the 
 // AliFlowTraks to fill phi histograms .
 
 cout << " AliFlowWeighter::WeightEvent(" << fFlowEvent << " )   -   " << fEventNumber << endl ;
 if(!fFlowEvent) { return kFALSE ; }

 if(fFlowSelect->Select(fFlowEvent))	 // event selected - here below the ANALYSIS FLAGS are setted -
 {
  fFlowTracks = fFlowEvent->TrackCollection() ; 
  fNumberOfTracks = fFlowTracks->GetEntries() ;
  cout << "       event ID = " << fFlowEvent->EventID() << " :  found " << fNumberOfTracks << " AliFlowTracks . " << endl ;  
  fFlowEvent->SetNoWgt() ;
  fFlowEvent->SetSelections(fFlowSelect) ; // does the selection of tracks for r.p. calculation (sets flags in AliFlowTrack)
  
  TracksLoop(fFlowTracks) ;  
 }
 else 
 {
  cout << "       * " << fEventNumber << " (event ID = " << fFlowEvent->EventID() << ") discarded . " << endl ; 
  delete fFlowEvent  ; fFlowEvent = 0 ; 
  return kFALSE ;
 }
 fEventNumber++ ;
 
 return kTRUE ;
}
//-----------------------------------------------------------------------
// ###
//-----------------------------------------------------------------------
void AliFlowWeighter::TracksLoop(TObjArray* fFlowTracks) 
{
 // fills phi and PId histograms

 cout << " Tracks Loop . " << endl ; 
  
 Float_t phi, eta, zFirstPoint ; // , zLastPoint ;
 Int_t  fPidId ;
 Char_t pid[10] = "0" ;
 for(fTrackNumber=0;fTrackNumber<fNumberOfTracks;fTrackNumber++) 
 {
  fFlowTrack = (AliFlowTrack*)fFlowTracks->At(fTrackNumber) ;
  // cout << "Track n. " << fTrackNumber << endl ; fFlowTrack->Dump() ; 
 
  phi = fFlowTrack->Phi() ;		
  eta = fFlowTrack->Eta() ;		
  zFirstPoint = fFlowTrack->ZFirstPoint() ; 
  //zLastPoint = fFlowTrack->ZLastPoint() ;
  strcpy(pid,fFlowTrack->Pid()) ; 

  fPidId = -1 ;
  if(strstr(pid,"e"))  	    { fPidId = 0 ; }
  else if(strstr(pid,"mu")) { fPidId = 1 ; }
  else if(strstr(pid,"pi")) { fPidId = 2 ; }
  else if(strstr(pid,"k"))  { fPidId = 3 ; }
  else if(strstr(pid,"pr")) { fPidId = 4 ; }
  else if(strstr(pid,"d"))  { fPidId = 5 ; }

  // Looping over Selections and Harmonics
  for (int k = 0; k < Flow::nSels; k++) 
  {
   fFlowSelect->SetSelection(k) ;
   for (int j = 0; j < Flow::nHars; j++) 
   {
    fFlowSelect->SetHarmonic(j);
    if(fFlowSelect->Select(fFlowTrack))
    {
     Bool_t kTpcPlus  = kFALSE ;
     Bool_t kTpcMinus = kFALSE ;
     Bool_t kTpcAll   = kFALSE ;

     // Set Tpc (+ and -)
     if(fFlowTrack->FitPtsTPC())
     {
      if(zFirstPoint >= 0. && eta > 0.)	     { kTpcPlus  = kTRUE ; } 
      else if(zFirstPoint <= 0. && eta < 0.) { kTpcMinus = kTRUE ; }
      else				     { kTpcAll   = kTRUE ; }	  
     }
    
     // Phi distribution (particle for R.P.)
     Float_t wt = 1. ; // TMath::Abs(fFlowEvent->Weight(k, j, fFlowTrack)) ; 
     if(kTpcPlus)       { fHistFull[k].fHistFullHar[j].fHistPhiPlus->Fill(phi,wt) ;  } 
     else if(kTpcMinus) { fHistFull[k].fHistFullHar[j].fHistPhiMinus->Fill(phi,wt) ; } 
     else if(kTpcAll)   { fHistFull[k].fHistFullHar[j].fHistPhiAll->Fill(phi,wt) ;   } 
     fHistFull[k].fHistFullHar[j].fHistPhi->Fill(phi,wt) ;

     // PID Multiplicities (particle for R.P.) - just once for each selection
     if(j==0) 		{ fHistFull[k].fHistBayPidMult->Fill(fPidId) ; }
    }
   }
  }
 }
  
 return ;
}
//-----------------------------------------------------------------------
Bool_t AliFlowWeighter::Weightening()
{
 // Calculates weights, and fills PhiWgt histograms .
 // This is called at the end of the event loop .
 
 cout << " AliFlowWeighter::Weightening() " << endl ; cout << endl ;
 
 // PhiWgt histogram collection
 fPhiWgtHistList = new TOrdCollection(4*Flow::nSels*Flow::nHars) ;
 
 // Creates PhiWgt Histograms
 TString* histTitle ;
 for(Int_t k = 0; k < Flow::nSels; k++)
 {
  for(Int_t j = 0; j < Flow::nHars; j++) 
  {
   // Tpc plus
   histTitle = new TString("Flow_Phi_Weight_TPCplus_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPhiWgtPlus = new TH1D(histTitle->Data(),histTitle->Data(), fPhiBins, fPhiMin, fPhiMax);
   fHistFull[k].fHistFullHar[j].fHistPhiWgtPlus->Sumw2();
   fHistFull[k].fHistFullHar[j].fHistPhiWgtPlus->SetXTitle("Azimuthal Angles (rad)");
   fHistFull[k].fHistFullHar[j].fHistPhiWgtPlus->SetYTitle("PhiWgt");
   delete histTitle;
   // Tpc minus
   histTitle = new TString("Flow_Phi_Weight_TPCminus_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPhiWgtMinus = new TH1D(histTitle->Data(),histTitle->Data(), fPhiBins, fPhiMin, fPhiMax);
   fHistFull[k].fHistFullHar[j].fHistPhiWgtMinus->Sumw2();
   fHistFull[k].fHistFullHar[j].fHistPhiWgtMinus->SetXTitle("Azimuthal Angles (rad)");
   fHistFull[k].fHistFullHar[j].fHistPhiWgtMinus->SetYTitle("PhiWgt");
   delete histTitle;
   // Tpc cross
   histTitle = new TString("Flow_Phi_Weight_TPCcross_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPhiWgtAll = new TH1D(histTitle->Data(),histTitle->Data(), fPhiBins, fPhiMin, fPhiMax);
   fHistFull[k].fHistFullHar[j].fHistPhiWgtAll->Sumw2();
   fHistFull[k].fHistFullHar[j].fHistPhiWgtAll->SetXTitle("Azimuthal Angles (rad)");
   fHistFull[k].fHistFullHar[j].fHistPhiWgtAll->SetYTitle("PhiWgt");
   delete histTitle;
   // Tpc
   histTitle = new TString("Flow_Phi_Weight_TPC_Sel");
   *histTitle += k+1;
   histTitle->Append("_Har");
   *histTitle += j+1;
   fHistFull[k].fHistFullHar[j].fHistPhiWgt = new TH1D(histTitle->Data(),histTitle->Data(), fPhiBins, fPhiMin, fPhiMax);
   fHistFull[k].fHistFullHar[j].fHistPhiWgt->Sumw2();
   fHistFull[k].fHistFullHar[j].fHistPhiWgt->SetXTitle("Azimuthal Angles (rad)");
   fHistFull[k].fHistFullHar[j].fHistPhiWgt->SetYTitle("PhiWgt");
   delete histTitle;

   // Calculate PhiWgt
   Double_t meanPlus  = fHistFull[k].fHistFullHar[j].fHistPhiPlus->Integral() / (Double_t)fPhiBins ;
   Double_t meanMinus = fHistFull[k].fHistFullHar[j].fHistPhiMinus->Integral() / (Double_t)fPhiBins ;
   Double_t meanCross = fHistFull[k].fHistFullHar[j].fHistPhiAll->Integral() / (Double_t)fPhiBins ;
   Double_t meanTPC = fHistFull[k].fHistFullHar[j].fHistPhi->Integral() / (Double_t)fPhiBins ;

   // Tpc
   for(Int_t i=0;i<fPhiBins;i++) 
   {
    fHistFull[k].fHistFullHar[j].fHistPhiWgtPlus->SetBinContent(i+1,meanPlus);
    fHistFull[k].fHistFullHar[j].fHistPhiWgtPlus->SetBinError(i+1, 0.);
    fHistFull[k].fHistFullHar[j].fHistPhiWgtMinus->SetBinContent(i+1,meanMinus);
    fHistFull[k].fHistFullHar[j].fHistPhiWgtMinus->SetBinError(i+1, 0.);
    fHistFull[k].fHistFullHar[j].fHistPhiWgtAll->SetBinContent(i+1,meanCross);
    fHistFull[k].fHistFullHar[j].fHistPhiWgtAll->SetBinError(i+1, 0.);
    fHistFull[k].fHistFullHar[j].fHistPhiWgt->SetBinContent(i+1,meanTPC);
    fHistFull[k].fHistFullHar[j].fHistPhiWgt->SetBinError(i+1, 0.);
   }
   
   if(meanTPC==0) { cout << " Sel." << k << " , Har." << j << " :  empty histogram ! " << endl ; }
   else 
   {
    fHistFull[k].fHistFullHar[j].fHistPhiWgtPlus->Divide(fHistFull[k].fHistFullHar[j].fHistPhiPlus);
    fHistFull[k].fHistFullHar[j].fHistPhiWgtMinus->Divide(fHistFull[k].fHistFullHar[j].fHistPhiMinus);
    fHistFull[k].fHistFullHar[j].fHistPhiWgtAll->Divide(fHistFull[k].fHistFullHar[j].fHistPhiAll);
    fHistFull[k].fHistFullHar[j].fHistPhiWgt->Divide(fHistFull[k].fHistFullHar[j].fHistPhi);
   }

   fPhiWgtHistList->AddLast(fHistFull[k].fHistFullHar[j].fHistPhiWgtPlus);
   fPhiWgtHistList->AddLast(fHistFull[k].fHistFullHar[j].fHistPhiWgtMinus);
   fPhiWgtHistList->AddLast(fHistFull[k].fHistFullHar[j].fHistPhiWgtAll);
   fPhiWgtHistList->AddLast(fHistFull[k].fHistFullHar[j].fHistPhiWgt);
  }
 }

 return kTRUE ;
}
//-----------------------------------------------------------------------
void AliFlowWeighter::PrintBayesian(Int_t selN)  
{
 // Prints the normalized particle abundance of all events (selection selN).
 // Call this at the end of the loop, just before Finish() .

 if(selN>Flow::nSels) { selN = 0 ; } 
 Char_t* names[Flow::nPid] = {"e","mu","pi","k","p","d"} ;
 Double_t bayes = 0. ;
 Double_t totCount = (fHistFull[selN].fHistBayPidMult)->GetSumOfWeights() ;
 if(totCount) 
 { 
  cout << " Sel." << selN << " particles normalized abundance (tot. " << totCount << " tracks) : " ;
  for(Int_t ii=0;ii<Flow::nPid;ii++)
  {
   bayes = (fHistFull[selN].fHistBayPidMult->GetBinContent(ii+1) / totCount) ; 
   cout << bayes << "_" << names[ii] << " ; " ;
  }
 }
 else {  cout << " Sel." << selN << " :  empty P.Id. histogram ! " << endl ; }
 cout << endl ;
 
 return ;
}
//-----------------------------------------------------------------------


#endif
