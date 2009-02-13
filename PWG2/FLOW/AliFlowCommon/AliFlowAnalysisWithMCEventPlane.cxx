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
   fHistProIntFlow(NULL),
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

void AliFlowAnalysisWithMCEventPlane::WriteHistograms(TString outputFileName)
{
 //store the final results in output .root file
 TFile *output = new TFile(outputFileName.Data(),"RECREATE");
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
  
  fHistProIntFlow = new TProfile("fHistProIntFlow","fHistProIntFlow",1,0.,1.);
  fHistProIntFlow->SetLabelSize(0.06);
  (fHistProIntFlow->GetXaxis())->SetBinLabel(1,"v_{n}{2}");
  fHistProIntFlow->SetYTitle("");
  fHistList->Add(fHistProIntFlow);
 
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
            //no-name int. flow (to be improved):
            fHistProIntFlow->Fill(0.,dv2);
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
	    //differential flow (Pt, POI):
            fHistProDiffFlowPtPOI->Fill(dPt,dv2,1.);
            //differential flow (Eta, POI):
            fHistProDiffFlowEtaPOI->Fill(dEta,dv2,1.); 
	  }	      
	}//track selected
      }//loop over tracks
	  
    fEventNumber++;
    //    cout<<"@@@@@ "<<fEventNumber<<" events processed"<<endl;
  }
}

  //--------------------------------------------------------------------    
void AliFlowAnalysisWithMCEventPlane::Finish() {
   
  //*************make histograms etc. 
  if (fDebug) cout<<"AliFlowAnalysisWithMCEventPlane::Terminate()"<<endl;
   
  Int_t iNbinsPt  = AliFlowCommonConstants::GetNbinsPt();  
  Int_t iNbinsEta = AliFlowCommonConstants::GetNbinsEta(); 
         
  // no-name int. flow (to be improved):
  Double_t dV = fHistProIntFlow->GetBinContent(1);  
  Double_t dErrV = fHistProIntFlow->GetBinError(1); // to be improved (treatment of errors for non-Gaussian distribution needed!)  
  // fill no-name int. flow (to be improved):
  fCommonHistsRes->FillIntegratedFlow(dV,dErrV);
  cout<<"dV{MC} is       "<<dV<<" +- "<<dErrV<<endl;
  
  //RP:
  TH1F* fHistPtRP = fCommonHists->GetHistPtInt(); // to be improved (change "int" and "diff" to RP and POI in common control histos)
  Double_t dYieldPtRP = 0.;
  Double_t dVRP = 0.;
  Double_t dErrVRP = 0.;
  Double_t dSumRP = 0.;
  //differential flow (RP, Pt): 
  Double_t dvPtRP = 0.;           
  Double_t dErrvPtRP = 0.;
  for(Int_t b=1;b<iNbinsPt;b++)
  {
   dvPtRP    = fHistProDiffFlowPtRP->GetBinContent(b);
   dErrvPtRP = fHistProDiffFlowPtRP->GetBinError(b);//to be improved (treatment of errors for non-Gaussian distribution needed!)
   fCommonHistsRes->FillDifferentialFlowPtRP(b, dvPtRP, dErrvPtRP);
   if(fHistPtRP){
	//integrated flow (RP)
	dYieldPtRP = fHistPtRP->GetBinContent(b);
	dVRP += dvPtRP*dYieldPtRP;
	dSumRP += dYieldPtRP;
	//error on integrated flow
	dErrVRP += dYieldPtRP*dYieldPtRP*dErrvPtRP*dErrvPtRP;
      }
  }
  if (dSumRP != 0. ) {
    dVRP /= dSumRP;  //because pt distribution should be normalised
    dErrVRP /= (dSumRP*dSumRP);
    dErrVRP = TMath::Sqrt(dErrVRP); 
  }
  // fill integrated flow (RP):
  fCommonHistsRes->FillIntegratedFlowRP(dVRP,dErrVRP);
  cout<<"dV{MC} (RP) is  "<<dVRP<<" +- "<<dErrVRP<<endl;
  
  //differential flow (RP, Eta): 
  Double_t dvEtaRP = 0.;           
  Double_t dErrvEtaRP = 0.;
  for(Int_t b=1;b<iNbinsEta;b++)
  {
   dvEtaRP    = fHistProDiffFlowEtaRP->GetBinContent(b);
   dErrvEtaRP = fHistProDiffFlowEtaRP->GetBinError(b);//to be improved (treatment of errors for non-Gaussian distribution needed!)
   fCommonHistsRes->FillDifferentialFlowEtaRP(b, dvEtaRP, dErrvEtaRP);
  }
                                                                                                                                   
  //POI:
  TH1F* fHistPtPOI = fCommonHists->GetHistPtDiff(); // to be improved (change "int" and "diff" to RP and POI in common control histos)
  Double_t dYieldPtPOI = 0.;
  Double_t dVPOI = 0.;
  Double_t dErrVPOI = 0.;
  Double_t dSumPOI = 0.;
  Double_t dv2proPtPOI = 0.;
  Double_t dErrdifcombPtPOI = 0.; 
  Double_t dv2proEtaPOI = 0.;
  Double_t dErrdifcombEtaPOI = 0.;   
  //Pt:
  if(fHistProFlow && fHistProDiffFlowPtPOI) {//to be removed (fHistProFlow)
    for(Int_t b=1;b<iNbinsPt;b++){
      //dv2pro = fHistProFlow->GetBinContent(b);//to be removed
      //dErrdifcomb = fHistProFlow->GetBinError(b);//to be removed
      dv2proPtPOI = fHistProDiffFlowPtPOI->GetBinContent(b);
      dErrdifcombPtPOI = fHistProDiffFlowPtPOI->GetBinError(b);//to be improved (treatment of errors for non-Gaussian distribution needed!)
      //fill TH1D
      fCommonHistsRes->FillDifferentialFlowPtPOI(b, dv2proPtPOI, dErrdifcombPtPOI); 
      if (fHistPtPOI){
	//integrated flow (POI)
	dYieldPtPOI = fHistPtPOI->GetBinContent(b);
	dVPOI += dv2proPtPOI*dYieldPtPOI;
	dSumPOI += dYieldPtPOI;
	//error on integrated flow
	dErrVPOI += dYieldPtPOI*dYieldPtPOI*dErrdifcombPtPOI*dErrdifcombPtPOI;
      }
    }//end of for(Int_t b=0;b<iNbinsPt;b++)  
  } else { cout<<"fHistProFlow is NULL"<<endl; }
  if (dSumPOI != 0. ) {
    dVPOI /= dSumPOI;  //because pt distribution should be normalised
    dErrVPOI /= (dSumPOI*dSumPOI);
    dErrVPOI = TMath::Sqrt(dErrVPOI); 
  }
  cout<<"dV{MC} (POI) is "<<dVPOI<<" +- "<<dErrVPOI<<endl;

  fCommonHistsRes->FillIntegratedFlowPOI(dVPOI,dErrVPOI);
  
  //Eta:
  if(fHistProDiffFlowEtaPOI)
  {
   for(Int_t b=1;b<iNbinsEta;b++)
   {
    dv2proEtaPOI = fHistProDiffFlowEtaPOI->GetBinContent(b);
    dErrdifcombEtaPOI = fHistProDiffFlowEtaPOI->GetBinError(b);//to be improved (treatment of errors for non-Gaussian distribution needed!)
    //fill common hist results:
    fCommonHistsRes->FillDifferentialFlowEtaPOI(b, dv2proEtaPOI, dErrdifcombEtaPOI); 
   }
  }   
  
  cout<<endl;     	      	  
  cout<<".....finished"<<endl;
}

 
 
