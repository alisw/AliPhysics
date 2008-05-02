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

#define AliFlowAnalysisWithMCEventPlane_cxx
 
#include "Riostream.h"  //needed as include
#include "TFile.h"      //needed as include
#include "TProfile.h"   //needed as include
#include "TComplex.h"   //needed as include

class TH1F;

#include "AliFlowCommonConstants.h"    //needed as include
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "AliFlowAnalysisWithMCEventPlane.h"

class AliFlowVector;

// AliFlowAnalysisWithMCEventPlane:
// Description: Maker to analyze Flow from the generated MC reaction plane.
//              This class is used to get the real value of the flow 
//              to compare the other methods to when analysing simulated events
// author: N. van der Kolk (kolk@nikhef.nl)

ClassImp(AliFlowAnalysisWithMCEventPlane)

  //-----------------------------------------------------------------------
 
 AliFlowAnalysisWithMCEventPlane::AliFlowAnalysisWithMCEventPlane():
   fQ2sum(0),
   fEventNumber(0),
   fMult(0),
   fNbins(0),
   fEvent(0x0),
   fTrack(0x0),
   fDebug(kFALSE),
   fHistFileName(0),
   fHistFile(0),
   fCommonHists(0),
   fCommonHistsRes(0),
   fHistProFlow(0),
   fHistRP(0)

{

  // Constructor.
  fQ.Set(0.,0.);           // flow vector
  fQsum.Set(0.,0.);        // flow vector sum
}

 
 //-----------------------------------------------------------------------


 AliFlowAnalysisWithMCEventPlane::~AliFlowAnalysisWithMCEventPlane() 
 {
   //destructor
   
 }
 

//-----------------------------------------------------------------------
void AliFlowAnalysisWithMCEventPlane::Init() {

  //Define all histograms
  cout<<"---Analysis with the real MC Event Plane---"<<endl;

  Int_t fNbinsPt = AliFlowCommonConstants::GetNbinsPt();
  Double_t  fPtMin = AliFlowCommonConstants::GetPtMin();	     
  Double_t  fPtMax = AliFlowCommonConstants::GetPtMax();

  // analysis file (output)
  fHistFile = new TFile(fHistFileName.Data(),"RECREATE") ;

  fCommonHists = new AliFlowCommonHist("MC");
  fCommonHistsRes = new AliFlowCommonHistResults("MC");

  fHistProFlow = new TProfile("FlowPro_VPt_MC","FlowPro_VPt_MC",fNbinsPt,fPtMin,fPtMax);
  fHistProFlow->SetXTitle("Pt");
  fHistProFlow->SetYTitle("v2 (%)");

  fHistRP = new TH1F("Flow_RP_MC","Flow_RP_MC",100,0.,3.14);
  fHistRP->SetXTitle("Reaction Plane Angle");
  fHistRP->SetYTitle("Counts");

 
  fEventNumber = 0;  //set number of events to zero
      
} 
 
//-----------------------------------------------------------------------
 
void AliFlowAnalysisWithMCEventPlane::Make(AliFlowEventSimple* fEvent, Double_t fRP) {

  //Calculate v2 from the MC reaction plane
  if (fEvent) {
         
    //fill control histograms     
    fCommonHists->FillControlHistograms(fEvent);

    //get the Q vector from the FlowEvent
    fQ = fEvent->GetQ(); 
    //cout<<"fQ.Mod() = " << fQ.Mod() << endl;
    //for chi calculation:
    fQsum += fQ;
    //cout<<"fQsum.Mod() = "<<fQsum.Mod()<<endl;
    fQ2sum += fQ.Mod2();
    cout<<"fQ2sum = "<<fQ2sum<<endl;
        
    fHistRP->Fill(fRP);   
              
    //calculate flow
    //loop over the tracks of the event
    Int_t fNumberOfTracks = fEvent->NumberOfTracks(); 
    for (Int_t i=0;i<fNumberOfTracks;i++) 
      {
	fTrack = fEvent->GetTrack(i) ; 
	if (fTrack){
	  if (fTrack->UseForDifferentialFlow()) {
	    Double_t fPhi = fTrack->Phi();
	    //if (fPhi<0.) fPhi+=2*TMath::Pi();
	    //calculate flow v2:
	    Double_t fv2 = TMath::Cos(2*(fPhi-fRP));
	    Double_t fPt = fTrack->Pt();
	    //fill histogram
	    fHistProFlow->Fill(fPt,100*fv2);  
	  }  
	}//track selected
      }//loop over tracks
	  
    fEventNumber++;
    cout<<"@@@@@ "<<fEventNumber<<" events processed"<<endl;
  }
}

  //--------------------------------------------------------------------    
void AliFlowAnalysisWithMCEventPlane::Finish() {
   
  //*************make histograms etc. 
  if (fDebug) cout<<"AliFlowAnalysisWithMCEventPlane::Terminate()"<<endl;
     
  Int_t fNbinsPt = AliFlowCommonConstants::GetNbinsPt();
    
  TH1F* fHistPtDiff = fCommonHists->GetfHistPtDiff();
  Double_t fV = 0.;
  Double_t fErrV = 0.;
  Double_t fSum = 0.;
  for(Int_t b=0;b<fNbinsPt;b++){
    Double_t fv2pro = 0.;
    Double_t fErrdifcomb = 0.; //in case error from profile is correct
    if(fHistProFlow) {
      fv2pro = fHistProFlow->GetBinContent(b);
      fErrdifcomb = fHistProFlow->GetBinError(b);
      //fill TH1D
      fCommonHistsRes->FillDifferentialFlow(b, fv2pro, fErrdifcomb); 
      if (fHistPtDiff){
	//integrated flow
	Double_t fYield = fHistPtDiff->GetBinContent(b);
	fV += fv2pro/100*fYield ;
	fSum += fYield;
	//error on integrated flow
	fErrV += fYield*fYield*(fErrdifcomb/100)*(fErrdifcomb/100);
      }
    }
  }
  fV /= fSum;  //because pt distribution should be normalised
  fErrV /= fSum*fSum;
  fErrV = TMath::Sqrt(fErrV);
  cout<<"fV is "<<fV<<" +- "<<fErrV<<endl;
  fCommonHistsRes->FillIntegratedFlow(fV,fErrV); 

  // write to file
  fHistFile->Write();
    	  
  cout<<".....finished"<<endl;
 }

 
