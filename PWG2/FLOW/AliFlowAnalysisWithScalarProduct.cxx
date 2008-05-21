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

#define AliFlowAnalysisWithScalarProduct_cxx
 
#include "Riostream.h"  //needed as include
#include "TFile.h"      //needed as include
#include "TMath.h"
#include "TProfile.h"
#include "TVector2.h"

class TH1F;

#include "AliFlowCommonConstants.h"    //needed as include
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowCommonHist.h"
//#include "AliFlowCommonHistResults.h"
#include "AliFlowAnalysisWithScalarProduct.h"

class AliFlowVector;

// AliFlowAnalysisWithScalarProduct:
// Description: 
// Maker to analyze Flow with the Scalar product method.
//
// author: N. van der Kolk (kolk@nikhef.nl)

ClassImp(AliFlowAnalysisWithScalarProduct)

  //-----------------------------------------------------------------------
 
 AliFlowAnalysisWithScalarProduct::AliFlowAnalysisWithScalarProduct():
   fEventNumber(0),
   fEvent(0x0),
   fTrack(0x0),
   fDebug(kFALSE),
   fHistFileName(0),
   fHistFile(0),
   fHistProUQ(0),
   fCommonHists(0)
{

  // Constructor.
  fQ.Set(0.,0.);           // flow vector
  fU.Set(0.,0.);           // particle unit vector

}
 //-----------------------------------------------------------------------


 AliFlowAnalysisWithScalarProduct::~AliFlowAnalysisWithScalarProduct() 
 {
   //destructor
   
 }
 

//-----------------------------------------------------------------------
void AliFlowAnalysisWithScalarProduct::Init() {

  //Define all histograms
  cout<<"---Analysis with the Scalar Product Method---"<<endl;

  Int_t fNbinsPt = AliFlowCommonConstants::GetNbinsPt();
  Double_t  fPtMin = AliFlowCommonConstants::GetPtMin();	     
  Double_t  fPtMax = AliFlowCommonConstants::GetPtMax();

  // analysis file (output)
  fHistFile = new TFile(fHistFileName.Data(),"RECREATE") ;

  fHistProUQ = new TProfile("Flow_UQ_SP","Flow_UQ_SP",fNbinsPt,fPtMin,fPtMax);
  fHistProUQ->SetXTitle("p_t (GeV)");
  fHistProUQ->SetYTitle("<uQ>");

  fCommonHists = new AliFlowCommonHist("SP");
  //fCommonHistsRes = new AliFlowCommonHistResults("SP");
 
  fEventNumber = 0;  //set number of events to zero
      
} 
 
//-----------------------------------------------------------------------
 
void AliFlowAnalysisWithScalarProduct::Make(AliFlowEventSimple* fEvent) {

  //Fill histogram
  if (fEvent) {

    //fill control histograms     
    fCommonHists->FillControlHistograms(fEvent);
         
    //get the Q vector from the FlowEvent
    fQ = fEvent->GetQ();
    //Double_t fMult =  fQ.GetMult();
            
    //loop over the tracks of the event
    Int_t fNumberOfTracks = fEvent->NumberOfTracks(); 
    for (Int_t i=0;i<fNumberOfTracks;i++) 
      {

	fTrack = fEvent->GetTrack(i) ; 
	if (fTrack){
	  if (fTrack->UseForDifferentialFlow()) {
	  Double_t fPhi = fTrack->Phi();

	  //calculate fU
	  Double_t fUX = TMath::Cos(2*fPhi);
	  Double_t fUY = TMath::Sin(2*fPhi);
	  fU.Set(fUX,fUY);
	  Double_t fModulus = fU.Mod();
	  if (fModulus!=0.) fU.Set(fUX/fModulus,fUY/fModulus);  // make length 1
	  else cerr<<"fModulus is zero!"<<endl;

	  TVector2 fQm = fQ;
	  //subtrackt particle from the flowvector if used to define it
	  if (fTrack->UseForIntegratedFlow()) {
	    Double_t fQmX = fQm.X() - fUX;
	    Double_t fQmY = fQm.Y() - fUY;
	    fQm.Set(fQmX,fQmY);
	  }

	  //Double_t fUQ = scalar product of fU and fQm
	  Double_t fUQ = fU*fQm;
	  Double_t fPt = fTrack->Pt();
	  //fill the profile histogram
	  fHistProUQ->Fill(fPt,fUQ); 
	  }  
	}//track selected
      }//loop over tracks
	  
    fEventNumber++;
    cout<<"@@@@@ "<<fEventNumber<<" events processed"<<endl;
  }
}

  //--------------------------------------------------------------------    
void AliFlowAnalysisWithScalarProduct::Finish() {
   
  //*************make histograms etc. 
  if (fDebug) cout<<"AliFlowAnalysisWithScalarProduct::Terminate()"<<endl;

  fHistProUQ->Draw();

  // write to file
  fHistFile->Write();
    	  
  cout<<".....finished"<<endl;
 }

 