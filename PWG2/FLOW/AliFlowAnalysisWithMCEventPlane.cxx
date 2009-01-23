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
#include "TList.h"

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
   fQsum(NULL),
   fQ2sum(0),
   fEventNumber(0),
   fDebug(kFALSE),
   fHistList(NULL),
   fCommonHists(NULL),
   fCommonHistsRes(NULL),
   fHistProFlow(NULL),
   fHistRP(NULL),
   fHistProIntFlowRP(NULL),
   fHistProDiffFlowPtRP(NULL),
   fHistProDiffFlowEtaRP(NULL),
   fHistProDiffFlowPtPOI(NULL),
   fHistProDiffFlowEtaPOI(NULL)
{

  // Constructor.
  fHistList = new TList();

  fQsum = new TVector2;        // flow vector sum
}

 
 //-----------------------------------------------------------------------


 AliFlowAnalysisWithMCEventPlane::~AliFlowAnalysisWithMCEventPlane() 
 {
   //destructor
   delete fHistList;
   delete fQsum;
 }
 
//-----------------------------------------------------------------------

void AliFlowAnalysisWithMCEventPlane::WriteHistograms(TString* outputFileName)
{
 //store the final results in output .root file
 TFile *output = new TFile(outputFileName->Data(),"RECREATE");
 output->WriteObject(fHistList, "cobjMCEP","SingleKey");
 delete output;
}

//-----------------------------------------------------------------------
void AliFlowAnalysisWithMCEventPlane::Init() {

  //Define all histograms
  cout<<"---Analysis with the real MC Event Plane---"<<endl;

  Int_t iNbinsPt = AliFlowCommonConstants::GetNbinsPt();
  Double_t dPtMin = AliFlowCommonConstants::GetPtMin();	     
  Double_t dPtMax = AliFlowCommonConstants::GetPtMax();
  
  Int_t iNbinsEta = AliFlowCommonConstants::GetNbinsEta();
  Double_t dEtaMin = AliFlowCommonConstants::GetEtaMin();	     
  Double_t dEtaMax = AliFlowCommonConstants::GetEtaMax();  

  fCommonHists = new AliFlowCommonHist("AliFlowCommonHistMCEP");
  fHistList->Add(fCommonHists);
  fCommonHistsRes = new AliFlowCommonHistResults("AliFlowCommonHistResultsMCEP");
  fHistList->Add(fCommonHistsRes);
  
  fHistProFlow = new TProfile("FlowPro_VPt_MCEP","FlowPro_VPt_MCEP",iNbinsPt,dPtMin,dPtMax);
  fHistProFlow->SetXTitle("P_{t}");
  fHistProFlow->SetYTitle("v_{2}");
  fHistList->Add(fHistProFlow);

  fHistRP = new TH1F("Flow_RP_MCEP","Flow_RP_MCEP",100,0.,3.14);
  fHistRP->SetXTitle("Reaction Plane Angle");
  fHistRP->SetYTitle("Counts");
  fHistList->Add(fHistRP);
  
  fHistProIntFlowRP = new TProfile("fHistProIntFlowRP","fHistProIntFlowRP",1,0.,1.);
  fHistProIntFlowRP->SetLabelSize(0.06);
  (fHistProIntFlowRP->GetXaxis())->SetBinLabel(1,"v_{n}{2}");
  fHistProIntFlowRP->SetYTitle("");
  fHistList->Add(fHistProIntFlowRP);
 
  fHistProDiffFlowPtRP = new TProfile("fHistProDiffFlowPtRP","fHistProDiffFlowPtRP",iNbinsPt,dPtMin,dPtMax);
  fHistProDiffFlowPtRP->SetXTitle("P_{t}");
  fHistProDiffFlowPtRP->SetYTitle("");
  fHistList->Add(fHistProDiffFlowPtRP);  
  
  fHistProDiffFlowEtaRP = new TProfile("fHistProDiffFlowEtaRP","fHistProDiffFlowEtaRP",iNbinsEta,dEtaMin,dEtaMax);
  fHistProDiffFlowEtaRP->SetXTitle("#eta");
  fHistProDiffFlowEtaRP->SetYTitle("");
  fHistList->Add(fHistProDiffFlowEtaRP);
  
  fHistProDiffFlowPtPOI = new TProfile("fHistProDiffFlowPtPOI","fHistProDiffFlowPtPOI",iNbinsPt,dPtMin,dPtMax);
  fHistProDiffFlowPtPOI->SetXTitle("P_{t}");
  fHistProDiffFlowPtPOI->SetYTitle("");
  fHistList->Add(fHistProDiffFlowPtPOI);  
  
  fHistProDiffFlowEtaPOI = new TProfile("fHistProDiffFlowEtaPOI","fHistProDiffFlowEtaPOI",iNbinsEta,dEtaMin,dEtaMax);
  fHistProDiffFlowEtaPOI->SetXTitle("#eta");
  fHistProDiffFlowEtaPOI->SetYTitle("");
  fHistList->Add(fHistProDiffFlowEtaPOI);          
 
  fEventNumber = 0;  //set number of events to zero
      
} 
 
//-----------------------------------------------------------------------
 
void AliFlowAnalysisWithMCEventPlane::Make(AliFlowEventSimple* anEvent, Double_t aRP) {

  //Calculate v2 from the MC reaction plane
  if (anEvent) {
         
    //fill control histograms     
    fCommonHists->FillControlHistograms(anEvent);

    //get the Q vector from the FlowEvent
    AliFlowVector vQ = anEvent->GetQ(); 
    //cout<<"vQ.Mod() = " << vQ.Mod() << endl;
    //for chi calculation:
    *fQsum += vQ;
    //cout<<"fQsum.Mod() = "<<fQsum.Mod()<<endl;
    fQ2sum += vQ.Mod2();
    //cout<<"fQ2sum = "<<fQ2sum<<endl;
        
    fHistRP->Fill(aRP);   
    
    Double_t dPhi = 0.;
    Double_t dv2  = 0.;
    Double_t dPt  = 0.;
    Double_t dEta = 0.;
    //Double_t dPi = TMath::Pi();                   
                                                             
    //calculate flow
    //loop over the tracks of the event
    Int_t iNumberOfTracks = anEvent->NumberOfTracks(); 
    for (Int_t i=0;i<iNumberOfTracks;i++) 
      {
	AliFlowTrackSimple* pTrack = anEvent->GetTrack(i) ; 
	if (pTrack){
	  if (pTrack->UseForIntegratedFlow()){
            dPhi = pTrack->Phi();
            dv2 = TMath::Cos(2*(dPhi-aRP));
	    dPt = pTrack->Pt();
	    dEta = pTrack->Eta();
            //integrated flow (RP):
            fHistProIntFlowRP->Fill(0.,dv2);
            //differential flow (Pt, RP):
            fHistProDiffFlowPtRP->Fill(dPt,dv2,1.);
            //differential flow (Eta, RP):
            fHistProDiffFlowEtaRP->Fill(dEta,dv2,1.);
          }
	  if (pTrack->UseForDifferentialFlow()) {
	    dPhi = pTrack->Phi();
	    //if (dPhi<0.) dPhi+=2*TMath::Pi();
	    //calculate flow v2:
	    dv2 = TMath::Cos(2*(dPhi-aRP));
	    dPt = pTrack->Pt();
	    dEta = pTrack->Eta();
	    //fill histogram
	    fHistProFlow->Fill(dPt,dv2);//to be removed 
	    //differential flow (Pt, POI):
            fHistProDiffFlowPtPOI->Fill(dPt,dv2,1.);
            //differential flow (Eta, POI):
            fHistProDiffFlowEtaPOI->Fill(dEta,dv2,1.); 
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
   
  Int_t iNbinsPt  = AliFlowCommonConstants::GetNbinsPt();  
  Int_t iNbinsEta = AliFlowCommonConstants::GetNbinsEta(); 
         
  //RP:
  //integrated flow:
  Double_t dVRP = fHistProIntFlowRP->GetBinContent(1);  
  Double_t dErrVRP = fHistProIntFlowRP->GetBinError(1);  
  fCommonHistsRes->FillIntegratedFlowRP(dVRP,dErrVRP);
  //differential flow (Pt): 
  Double_t dvPtRP = 0.;           
  Double_t dErrvPtRP = 0.;
  for(Int_t b=0;b<iNbinsPt;b++)
  {
   dvPtRP    = fHistProDiffFlowPtRP->GetBinContent(b+1);
   dErrvPtRP = fHistProDiffFlowPtRP->GetBinError(b+1);
   fCommonHistsRes->FillDifferentialFlowPtRP(b, dvPtRP , dErrvPtRP);
  }
  //differential flow (Eta): 
  Double_t dvEtaRP = 0.;           
  Double_t dErrvEtaRP = 0.;
  for(Int_t b=0;b<iNbinsEta;b++)
  {
   dvEtaRP    = fHistProDiffFlowEtaRP->GetBinContent(b+1);
   dErrvEtaRP = fHistProDiffFlowEtaRP->GetBinError(b+1);
   fCommonHistsRes->FillDifferentialFlowEtaRP(b, dvEtaRP , dErrvEtaRP);
  }
                                                                                                                                   
  //POI:
  TH1F* fHistPtDiff = fCommonHists->GetHistPtDiff();
  Double_t dV = 0.;
  Double_t dErrV = 0.;
  Double_t dSum = 0.;
  Double_t dv2proPt = 0.;
  Double_t dErrdifcombPt = 0.; 
  Double_t dv2proEta = 0.;
  Double_t dErrdifcombEta = 0.;   
  //Pt:
  if(fHistProFlow && fHistProDiffFlowPtPOI) {//to be removed (fHistProFlow)
    for(Int_t b=0;b<iNbinsPt;b++){
      //dv2pro = fHistProFlow->GetBinContent(b);//to be removed
      //dErrdifcomb = fHistProFlow->GetBinError(b);//to be removed
      dv2proPt = fHistProDiffFlowPtPOI->GetBinContent(b);
      dErrdifcombPt = fHistProDiffFlowPtPOI->GetBinError(b); //in case error from profile is correct
      //fill TH1D
      fCommonHistsRes->FillDifferentialFlow(b, dv2proPt, dErrdifcombPt);//to be removed
      fCommonHistsRes->FillDifferentialFlowPtPOI(b, dv2proPt, dErrdifcombPt); 
      if (fHistPtDiff){
	//integrated flow
	Double_t dYield = fHistPtDiff->GetBinContent(b);
	dV += dv2proPt*dYield ;
	dSum += dYield;
	//error on integrated flow
	dErrV += dYield*dYield*dErrdifcombPt*dErrdifcombPt;
      }
    }  
  } else { cout<<"fHistProFlow is NULL"<<endl; }
  if (dSum != 0. ) {
    dV /= dSum;  //because pt distribution should be normalised
    dErrV /= dSum*dSum;
    dErrV = TMath::Sqrt(dErrV); }
  cout<<"dV is "<<dV<<" +- "<<dErrV<<endl;
  fCommonHistsRes->FillIntegratedFlow(dV,dErrV);//to be removed 
  fCommonHistsRes->FillIntegratedFlowPOI(dV,dErrV);
  
  
  //Eta:
  if(fHistProDiffFlowEtaPOI)
  {
   for(Int_t b=0;b<iNbinsEta;b++)
   {
    dv2proEta = fHistProDiffFlowEtaPOI->GetBinContent(b);
    dErrdifcombEta = fHistProDiffFlowEtaPOI->GetBinError(b); //in case error from profile is correct
    //fill common hist results:
    fCommonHistsRes->FillDifferentialFlowEtaPOI(b, dv2proEta, dErrdifcombEta); 
   }
  }
       	      	  
  cout<<".....finished"<<endl;
}

 
 
