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
#include "TList.h"
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
   fDebug(kFALSE),
   fHistList(NULL),
   fHistProUQ(NULL),
   fCommonHists(NULL)
{
  // Constructor.
  fHistList = new TList();
}
 //-----------------------------------------------------------------------


 AliFlowAnalysisWithScalarProduct::~AliFlowAnalysisWithScalarProduct() 
 {
   //destructor
   delete fHistList;
 }
 

//-----------------------------------------------------------------------

void AliFlowAnalysisWithScalarProduct::WriteHistograms(TString* outputFileName)
{
 //store the final results in output .root file

  TFile *output = new TFile(outputFileName->Data(),"RECREATE");
  output->WriteObject(fHistList, "cobjSP","SingleKey");
  delete output;
}

//-----------------------------------------------------------------------

void AliFlowAnalysisWithScalarProduct::WriteHistograms(TString outputFileName)
{
 //store the final results in output .root file

  TFile *output = new TFile(outputFileName.Data(),"RECREATE");
  output->WriteObject(fHistList, "cobjSP","SingleKey");
  delete output;
}

//-----------------------------------------------------------------------
void AliFlowAnalysisWithScalarProduct::Init() {

  //Define all histograms
  cout<<"---Analysis with the Scalar Product Method--- Init"<<endl;

  Int_t iNbinsPt = AliFlowCommonConstants::GetNbinsPt();
  Double_t  dPtMin = AliFlowCommonConstants::GetPtMin();	     
  Double_t  dPtMax = AliFlowCommonConstants::GetPtMax();

  fHistProUQ = new TProfile("Flow_UQ_SP","Flow_UQ_SP",iNbinsPt,dPtMin,dPtMax);
  fHistProUQ->SetXTitle("p_t (GeV)");
  fHistProUQ->SetYTitle("<uQ>");
  fHistList->Add(fHistProUQ);

  fCommonHists = new AliFlowCommonHist("AliFlowCommonHistSP");
  fHistList->Add(fCommonHists);
  
  //fCommonHistsRes = new AliFlowCommonHistResults("SP");
  
  fEventNumber = 0;  //set number of events to zero    
}

//-----------------------------------------------------------------------
 
void AliFlowAnalysisWithScalarProduct::Make(AliFlowEventSimple* anEvent) {

  //Fill histogram
  if (anEvent) {

    //fill control histograms     
    fCommonHists->FillControlHistograms(anEvent);
         
    //get the Q vector from the FlowEvent
    AliFlowVector vQ = anEvent->GetQ();
                
    //loop over the tracks of the event
    AliFlowTrackSimple*   pTrack = NULL; 
    Int_t iNumberOfTracks = anEvent->NumberOfTracks(); 
    for (Int_t i=0;i<iNumberOfTracks;i++) 
      {
	pTrack = anEvent->GetTrack(i) ; 
	if (pTrack){
	  if (pTrack->UseForDifferentialFlow()) {
	  Double_t dPhi = pTrack->Phi();

	  //calculate vU
	  TVector2 vU;
	  Double_t dUX = TMath::Cos(2*dPhi);
	  Double_t dUY = TMath::Sin(2*dPhi);
	  vU.Set(dUX,dUY);
	  Double_t dModulus = vU.Mod();
	  if (dModulus!=0.) vU.Set(dUX/dModulus,dUY/dModulus);  // make length 1
	  else cerr<<"dModulus is zero!"<<endl;

	  TVector2 vQm = vQ;
	  //subtrackt particle from the flowvector if used to define it
	  if (pTrack->UseForIntegratedFlow()) {
	    Double_t dQmX = vQm.X() - dUX;
	    Double_t dQmY = vQm.Y() - dUY;
	    vQm.Set(dQmX,dQmY);
	  }

	  //dUQ = scalar product of vU and vQm
	  Double_t dUQ = vU * vQm;
	  Double_t dPt = pTrack->Pt();
	  //fill the profile histogram
	  fHistProUQ->Fill(dPt,dUQ); 
	  }  
	}//track selected
      }//loop over tracks
	 
    fEventNumber++;
    //    cout<<"@@@@@ "<<fEventNumber<<" events processed"<<endl;
  }
}

  //--------------------------------------------------------------------    
void AliFlowAnalysisWithScalarProduct::Finish() {
   
  //*************make histograms etc. 
  if (fDebug) cout<<"AliFlowAnalysisWithScalarProduct::Terminate()"<<endl;

  //  fHistProUQ->Draw();
     	  
  cout<<".....finished"<<endl;
 }


