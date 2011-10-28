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

// AliFlowAnalysisWithLYZEventPlane:
//
// Class to do flow analysis with the event plane from the LYZ method
//
// author: N. van der Kolk (kolk@nikhef.nl)

/*
$Log$
*/ 

//#define AliFlowAnalysisWithLYZEventPlane_cxx
 
#include "Riostream.h"  //needed as include
#include "TMath.h"   //needed as include
#include "TProfile.h"   //needed as include
#include "TH1F.h"
#include "TFile.h"
#include "TList.h"
#include "TVector2.h"

#include "AliFlowLYZConstants.h"    //needed as include
#include "AliFlowCommonConstants.h" //needed as include
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "AliFlowLYZEventPlane.h"
#include "AliFlowAnalysisWithLYZEventPlane.h"
#include "AliFlowVector.h"


ClassImp(AliFlowAnalysisWithLYZEventPlane)

  //-----------------------------------------------------------------------
 
 AliFlowAnalysisWithLYZEventPlane::AliFlowAnalysisWithLYZEventPlane():
   fHistList(NULL),
   fSecondRunList(NULL),
   fSecondReDtheta(NULL),
   fSecondImDtheta(NULL),
   fFirstr0theta(NULL),
   fHistProVetaRP(NULL),
   fHistProVetaPOI(NULL),
   fHistProVPtRP(NULL),
   fHistProVPtPOI(NULL),
   fHistProWr(NULL),
   fHistProWrCorr(NULL),
   fHistQsumforChi(NULL),
   fHistDeltaPhi(NULL),
   fHistDeltaPhi2(NULL),
   fHistDeltaPhihere(NULL),
   fHistPhiEP(NULL),
   fHistPhiEPhere(NULL),
   fHistPhiLYZ(NULL),
   fHistPhiLYZ2(NULL),
   fCommonHists(NULL),
   fCommonHistsRes(NULL),
   fEventNumber(0),
   fQsum(NULL),
   fQ2sum(0)
{

  // Constructor.
  fQsum = new TVector2();        // flow vector sum

  fHistList = new TList();
  fSecondRunList = new TList();
}

 

 //-----------------------------------------------------------------------


 AliFlowAnalysisWithLYZEventPlane::~AliFlowAnalysisWithLYZEventPlane() 
 {
   //destructor
   delete fQsum;
   delete fHistList;
   delete fSecondRunList;
 }
 

//-----------------------------------------------------------------------

void AliFlowAnalysisWithLYZEventPlane::WriteHistograms(TString* outputFileName)
{
 //store the final results in output .root file

  TFile *output = new TFile(outputFileName->Data(),"RECREATE");
  //output->WriteObject(fHistList, "cobjLYZEP","SingleKey");
  fHistList->SetName("cobjLYZEP");
  fHistList->SetOwner(kTRUE);
  fHistList->Write(fHistList->GetName(), TObject::kSingleKey);
  delete output;
}

//-----------------------------------------------------------------------

void AliFlowAnalysisWithLYZEventPlane::WriteHistograms(TString outputFileName)
{
 //store the final results in output .root file

  TFile *output = new TFile(outputFileName.Data(),"RECREATE");
  //output->WriteObject(fHistList, "cobjLYZEP","SingleKey");
  fHistList->SetName("cobjLYZEP");
  fHistList->SetOwner(kTRUE);
  fHistList->Write(fHistList->GetName(), TObject::kSingleKey);
  delete output;
}

//-----------------------------------------------------------------------

void AliFlowAnalysisWithLYZEventPlane::WriteHistograms(TDirectoryFile *outputFileName)
{
 //store the final results in output .root file
 fHistList->SetName("cobjLYZEP");
 fHistList->SetOwner(kTRUE);
 outputFileName->Add(fHistList);
 outputFileName->Write(outputFileName->GetName(), TObject::kSingleKey);
}

//-----------------------------------------------------------------------

void AliFlowAnalysisWithLYZEventPlane::Init() {

  //Initialise all histograms
  cout<<"---Analysis with Lee Yang Zeros Event Plane Method---"<<endl;

  //save old value and prevent histograms from being added to directory
  //to avoid name clashes in case multiple analaysis objects are used
  //in an analysis
  Bool_t oldHistAddStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
 
  //input histograms
  if (fSecondRunList) {
    fSecondReDtheta = (TProfile*)fSecondRunList->FindObject("Second_FlowPro_ReDtheta_LYZSUM");
    fHistList->Add(fSecondReDtheta);

    fSecondImDtheta = (TProfile*)fSecondRunList->FindObject("Second_FlowPro_ImDtheta_LYZSUM");
    fHistList->Add(fSecondImDtheta);
    
    fFirstr0theta = (TProfile*)fSecondRunList->FindObject("First_FlowPro_r0theta_LYZSUM");
    fHistList->Add(fFirstr0theta);

    //warnings
    if (!fSecondReDtheta) {cout<<"fSecondReDtheta is NULL!"<<endl; }
    if (!fSecondImDtheta) {cout<<"fSecondImDtheta is NULL!"<<endl; }
    if (!fFirstr0theta)   {cout<<"fFirstr0theta is NULL!"<<endl; }

  }

  fCommonHists = new AliFlowCommonHist("AliFlowCommonHistLYZEP");
  fHistList->Add(fCommonHists);
  
  fCommonHistsRes = new AliFlowCommonHistResults("AliFlowCommonHistResultsLYZEP"); 
  fHistList->Add(fCommonHistsRes); 
    
  Int_t iNbinsPt = AliFlowCommonConstants::GetMaster()->GetNbinsPt();
  Int_t iNbinsEta = AliFlowCommonConstants::GetMaster()->GetNbinsEta();
  Double_t  dPtMin  = AliFlowCommonConstants::GetMaster()->GetPtMin();	     
  Double_t  dPtMax  = AliFlowCommonConstants::GetMaster()->GetPtMax();
  Double_t  dEtaMin = AliFlowCommonConstants::GetMaster()->GetEtaMin();	     
  Double_t  dEtaMax = AliFlowCommonConstants::GetMaster()->GetEtaMax();

  fHistProVetaRP = new TProfile("FlowPro_VetaRP_LYZEP","FlowPro_VetaRP_LYZEP",iNbinsEta,dEtaMin,dEtaMax);
  fHistProVetaRP->SetXTitle("rapidity");
  fHistProVetaRP->SetYTitle("v_{2}(#eta) for RP selection");
  fHistList->Add(fHistProVetaRP);

  fHistProVetaPOI = new TProfile("FlowPro_VetaPOI_LYZEP","FlowPro_VetaPOI_LYZEP",iNbinsEta,dEtaMin,dEtaMax);
  fHistProVetaPOI->SetXTitle("rapidity");
  fHistProVetaPOI->SetYTitle("v_{2}(#eta) for POI selection");
  fHistList->Add(fHistProVetaPOI);

  fHistProVPtRP = new TProfile("FlowPro_VPtRP_LYZEP","FlowPro_VPtRP_LYZEP",iNbinsPt,dPtMin,dPtMax);
  fHistProVPtRP->SetXTitle("Pt");
  fHistProVPtRP->SetYTitle("v_{2}(p_{T}) for RP selection");
  fHistList->Add(fHistProVPtRP);

  fHistProVPtPOI = new TProfile("FlowPro_VPtPOI_LYZEP","FlowPro_VPtPOI_LYZEP",iNbinsPt,dPtMin,dPtMax);
  fHistProVPtPOI->SetXTitle("p_{T}");
  fHistProVPtPOI->SetYTitle("v_{2}(p_{T}) for POI selection");
  fHistList->Add(fHistProVPtPOI);

  fHistProWr = new TProfile("FlowPro_Wr_LYZEP","FlowPro_Wr_LYZEP",100,0.,0.25);
  fHistProWr->SetXTitle("Q");
  fHistProWr->SetYTitle("Wr");
  fHistList->Add(fHistProWr);

  fHistQsumforChi = new TH1F("Flow_QsumforChi_LYZEP","Flow_QsumforChi_LYZEP",3,-1.,2.);
  fHistQsumforChi->SetXTitle("Qsum.X , Qsum.Y, Q2sum");
  fHistQsumforChi->SetYTitle("value");
  fHistList->Add(fHistQsumforChi);

  fHistDeltaPhi = new TH1F("Flow_DeltaPhi_LYZEP","Flow_DeltaPhi_LYZEP",100,0.,3.14);
  fHistDeltaPhi->SetXTitle("DeltaPhi");
  fHistDeltaPhi->SetYTitle("Counts");
  fHistList->Add(fHistDeltaPhi);

  fHistPhiLYZ = new TH1F("Flow_PhiLYZ_LYZEP","Flow_PhiLYZ_LYZEP",100,0.,3.14);
  fHistPhiLYZ->SetXTitle("Phi from LYZ");
  fHistPhiLYZ->SetYTitle("Counts");
  fHistList->Add(fHistPhiLYZ);

  fHistPhiEP = new TH1F("Flow_PhiEP_LYZEP","Flow_PhiEP_LYZEP",100,0.,3.14);
  fHistPhiEP->SetXTitle("Phi from EP");
  fHistPhiEP->SetYTitle("Counts");
  fHistList->Add(fHistPhiEP);

  fEventNumber = 0;  //set number of events to zero
      
  //restore old status
  TH1::AddDirectory(oldHistAddStatus);
} 
 
//-----------------------------------------------------------------------
 
void AliFlowAnalysisWithLYZEventPlane::Make(AliFlowEventSimple* anEvent, AliFlowLYZEventPlane* aLYZEP) {
  
  //Get the event plane and weight for each event
  if (anEvent) {
         
    //fill control histograms     
    fCommonHists->FillControlHistograms(anEvent);

    //get the Q vector from the FlowEvent
    AliFlowVector vQ = anEvent->GetQ(); 
    //if (vQ.X()== 0. && vQ.Y()== 0. ) { cout<<"Q vector is NULL!"<<endl; } //coding violation
    //Weight with the multiplicity
    Double_t dQX = 0.;
    Double_t dQY = 0.;
    if (TMath::AreEqualAbs(vQ.GetMult(),0.0,1e-10)) {
      dQX = vQ.X()/vQ.GetMult();
      dQY = vQ.Y()/vQ.GetMult();
    } else {cerr<<"vQ.GetMult() is zero!"<<endl; }
    vQ.Set(dQX,dQY);
    //cout<<"vQ("<<dQX<<","<<dQY<<")"<<endl;

    //for chi calculation:
    *fQsum += vQ;
    fHistQsumforChi->SetBinContent(1,fQsum->X());
    fHistQsumforChi->SetBinContent(2,fQsum->Y());
    fQ2sum += vQ.Mod2();
    fHistQsumforChi->SetBinContent(3,fQ2sum);
    //cout<<"fQ2sum = "<<fQ2sum<<endl;

    //call AliFlowLYZEventPlane::CalculateRPandW() here!
    aLYZEP->CalculateRPandW(vQ);

    Double_t dWR = aLYZEP->GetWR();     
    Double_t dRP = aLYZEP->GetPsi();

    //fHistProWr->Fill(vQ.Mod(),dWR); //this weight is always positive
    fHistPhiLYZ->Fill(dRP);   
    
    //plot difference between event plane from EP-method and LYZ-method
    Double_t dRPEP = vQ.Phi()/2;                              //gives distribution from (0 to pi)
    //Double_t dRPEP = 0.5*TMath::ATan2(vQ.Y(),vQ.X());       //gives distribution from (-pi/2 to pi/2)
    //cout<<"dRPEP = "<< dRPEP <<endl;
    fHistPhiEP->Fill(dRPEP);

    Double_t dDeltaPhi = dRPEP - dRP;
    if (dDeltaPhi < 0.) { dDeltaPhi += TMath::Pi(); }        //to shift distribution from (-pi/2 to pi/2) to (0 to pi)
    //cout<<"dDeltaPhi = "<<dDeltaPhi<<endl;
    fHistDeltaPhi->Fill(dDeltaPhi); 

    //Flip sign of WR
    Double_t dLow = TMath::Pi()/4.;
    Double_t dHigh = 3.*(TMath::Pi()/4.);
    if ((dDeltaPhi > dLow) && (dDeltaPhi < dHigh)){
      dRP -= (TMath::Pi()/2);
      dWR = -dWR;
      cerr<<"*** dRP modified ***"<<endl;
    }
    fHistProWr->Fill(vQ.Mod(),dWR); //corrected weight
       
    //calculate flow for RP and POI selections
    //loop over the tracks of the event
    Int_t iNumberOfTracks = anEvent->NumberOfTracks(); 
    for (Int_t i=0;i<iNumberOfTracks;i++) 
      {
	AliFlowTrackSimple* pTrack = anEvent->GetTrack(i) ; 
	if (pTrack){
	  Double_t dPhi = pTrack->Phi();
	  //if (dPhi<0.) fPhi+=2*TMath::Pi();
	  Double_t dPt  = pTrack->Pt();
	  Double_t dEta = pTrack->Eta();
	  //calculate flow v2:
	  Double_t dv2 = dWR * TMath::Cos(2*(dPhi-dRP));
	  if (pTrack->InRPSelection()) {
	    //fill histograms for RP selection
	    fHistProVetaRP -> Fill(dEta,dv2); 
	    fHistProVPtRP  -> Fill(dPt,dv2); 
	  }
	  if (pTrack->InPOISelection()) {
	    //fill histograms for POI selection
	    fHistProVetaPOI -> Fill(dEta,dv2); 
	    fHistProVPtPOI  -> Fill(dPt,dv2); 
	  }  
	}//track 
      }//loop over tracks
	  
    fEventNumber++;
    //    cout<<"@@@@@ "<<fEventNumber<<" events processed"<<endl;
  }
}

  //--------------------------------------------------------------------    
void AliFlowAnalysisWithLYZEventPlane::GetOutputHistograms(TList *outputListHistos){
 //get pointers to all output histograms (called before Finish()) 
 if (outputListHistos) {
    //Get the common histograms from the output list
    AliFlowCommonHist *pCommonHist = dynamic_cast<AliFlowCommonHist*> 
      (outputListHistos->FindObject("AliFlowCommonHistLYZEP"));
    AliFlowCommonHistResults *pCommonHistResults = dynamic_cast<AliFlowCommonHistResults*> 
      (outputListHistos->FindObject("AliFlowCommonHistResultsLYZEP"));

    TProfile* pHistProR0theta = dynamic_cast<TProfile*> 
      (outputListHistos->FindObject("First_FlowPro_r0theta_LYZSUM"));

    TProfile* pHistProVetaRP = dynamic_cast<TProfile*> 
      (outputListHistos->FindObject("FlowPro_VetaRP_LYZEP"));
    TProfile* pHistProVetaPOI = dynamic_cast<TProfile*> 
      (outputListHistos->FindObject("FlowPro_VetaPOI_LYZEP"));
    TProfile* pHistProVPtRP = dynamic_cast<TProfile*> 
      (outputListHistos->FindObject("FlowPro_VPtRP_LYZEP"));
    TProfile* pHistProVPtPOI = dynamic_cast<TProfile*> 
      (outputListHistos->FindObject("FlowPro_VPtPOI_LYZEP"));

    TH1F* pHistQsumforChi = dynamic_cast<TH1F*> 
      (outputListHistos->FindObject("Flow_QsumforChi_LYZEP"));

    if (pCommonHist && pCommonHistResults && pHistProR0theta &&
	pHistProVetaRP && pHistProVetaPOI && pHistProVPtRP && 
	pHistProVPtPOI && pHistQsumforChi ) {
      this -> SetCommonHists(pCommonHist);
      this -> SetCommonHistsRes(pCommonHistResults);
      this -> SetFirstr0theta(pHistProR0theta);
      this -> SetHistProVetaRP(pHistProVetaRP);
      this -> SetHistProVetaPOI(pHistProVetaPOI);
      this -> SetHistProVPtRP(pHistProVPtRP);
      this -> SetHistProVPtPOI(pHistProVPtPOI);
      this -> SetHistQsumforChi(pHistQsumforChi);
     }  
  } else { 
      cout<<"WARNING: Histograms needed to run Finish() are not accessable!"<<endl; 
    }    
}

  //--------------------------------------------------------------------    
void AliFlowAnalysisWithLYZEventPlane::Finish() {
   
  //plot histograms etc. 
  cout<<"AliFlowAnalysisWithLYZEventPlane::Finish()"<<endl;
  
  //constants:
  Double_t  dJ01 = 2.405; 
  Int_t iNtheta   = AliFlowLYZConstants::GetMaster()->GetNtheta();
  Int_t iNbinsPt  = AliFlowCommonConstants::GetMaster()->GetNbinsPt();
  Int_t iNbinsEta = AliFlowCommonConstants::GetMaster()->GetNbinsEta();
  //set the event number
  if (fCommonHists) {
    SetEventNumber((int)fCommonHists->GetHistQ()->GetEntries());
    //cout<<"number of events processed is "<<fEventNumber<<endl;
  } else {
    cout<<"Commonhist pointer is NULL."<<endl;
    cout<<"Leaving LYZ Event plane analysis!"<<endl;
    return;
  }

  //set the sum of Q vectors
  fQsum->Set(fHistQsumforChi->GetBinContent(1),fHistQsumforChi->GetBinContent(2));
  SetQ2sum(fHistQsumforChi->GetBinContent(3));  

  //calculate dV the mean of dVtheta
  Double_t  dVtheta = 0; 
  Double_t  dV = 0; 
  for (Int_t theta=0;theta<iNtheta;theta++)	{
    Double_t dR0 = fFirstr0theta->GetBinContent(theta+1); 
    if (TMath::AreEqualAbs(dR0,0.0,1e-10)) { dVtheta = dJ01/dR0 ;}
    dV += dVtheta;
  }
  dV /= iNtheta;

  //calculate dChi 
  Double_t  dSigma2 = 0;
  Double_t  dChi= 0;
  if (fEventNumber!=0) {
    *fQsum /= fEventNumber;
    //cerr<<"fQsum->X() = "<<fQsum->X()<<endl;
    //cerr<<"fQsum->Y() = "<<fQsum->Y()<<endl;
    fQ2sum /= fEventNumber;
    //cerr<<"fEventNumber = "<<fEventNumber<<endl;
    //cerr<<"fQ2sum = "<<fQ2sum<<endl;
    dSigma2 = fQ2sum - TMath::Power(fQsum->X(),2.) - TMath::Power(fQsum->Y(),2.) - TMath::Power(dV,2.);  //BP eq. 62
    //cerr<<"dSigma2"<<dSigma2<<endl;
    if (dSigma2>0) dChi = dV/TMath::Sqrt(dSigma2);
    else dChi = -1.;
    fCommonHistsRes->FillChi(dChi);

    // recalculate statistical errors on integrated flow
    //combining 5 theta angles to 1 relative error BP eq. 89
    Double_t dRelErr2comb = 0.;
    Int_t iEvts = fEventNumber; 
    if (iEvts!=0) {
      for (Int_t theta=0;theta<iNtheta;theta++){
	Double_t dTheta = ((double)theta/iNtheta)*TMath::Pi(); 
	Double_t dApluscomb = TMath::Exp((dJ01*dJ01)/(2*dChi*dChi)*
				       TMath::Cos(dTheta));
	Double_t dAmincomb = TMath::Exp(-(dJ01*dJ01)/(2*dChi*dChi)*
				      TMath::Cos(dTheta));
	dRelErr2comb += (1/(2*iEvts*(dJ01*dJ01)*TMath::BesselJ1(dJ01)*
			  TMath::BesselJ1(dJ01)))*
	  (dApluscomb*TMath::BesselJ0(2*dJ01*TMath::Sin(dTheta/2)) + 
	   dAmincomb*TMath::BesselJ0(2*dJ01*TMath::Cos(dTheta/2)));
      }
      dRelErr2comb /= iNtheta;
    }
    Double_t dRelErrcomb = TMath::Sqrt(dRelErr2comb);
    Double_t dVErr = dV*dRelErrcomb ; 
    fCommonHistsRes->FillIntegratedFlow(dV, dVErr); 

    cout<<"*************************************"<<endl;
    cout<<"*************************************"<<endl;
    cout<<"      Integrated flow from           "<<endl;
    cout<<"  Lee-Yang Zeroes Event Plane        "<<endl;
    cout<<endl;
    cout<<"dChi = "<<dChi<<endl;
    cout<<"dV = "<<dV<<" +- "<<dVErr<<endl;
    cout<<endl;
        
  }
  
  //copy content of profile into TH1D, add error and fill the AliFlowCommonHistResults

  //v as a function of eta for RP selection
  for(Int_t b=0;b<iNbinsEta;b++) {
    Double_t dv2pro  = fHistProVetaRP->GetBinContent(b);
    Double_t dNprime = fCommonHists->GetEntriesInEtaBinRP(b);  
    Double_t dErrdifcomb = 0.;  //set error to zero
    Double_t dErr2difcomb = 0.; //set error to zero
    //calculate error
    if (TMath::AreEqualAbs(dNprime,0.0,1e-10)) { 
      for (Int_t theta=0;theta<iNtheta;theta++) {
	Double_t dTheta = ((double)theta/iNtheta)*TMath::Pi(); 
	Double_t dApluscomb = TMath::Exp((dJ01*dJ01)/(2*dChi*dChi)*
					   TMath::Cos(dTheta));
	Double_t dAmincomb = TMath::Exp(-(dJ01*dJ01)/(2*dChi*dChi)*
					  TMath::Cos(dTheta));
	dErr2difcomb += (TMath::Cos(dTheta)/(4*dNprime*TMath::BesselJ1(dJ01)*
						 TMath::BesselJ1(dJ01)))*
	  ((dApluscomb*TMath::BesselJ0(2*dJ01*TMath::Sin(dTheta/2))) - 
	   (dAmincomb*TMath::BesselJ0(2*dJ01*TMath::Cos(dTheta/2))));
      } //loop over theta
    } 
      
    if (TMath::AreEqualAbs(dErr2difcomb, 0.0, 1e-10)) {
      dErr2difcomb /= iNtheta;
      dErrdifcomb = TMath::Sqrt(dErr2difcomb);
    }
    else {dErrdifcomb = 0.;}
    //fill TH1D
    fCommonHistsRes->FillDifferentialFlowEtaRP(b, dv2pro, dErrdifcomb); 
  } //loop over bins b
  
  
  //v as a function of eta for POI selection
  for(Int_t b=0;b<iNbinsEta;b++) {
    Double_t dv2pro  = fHistProVetaPOI->GetBinContent(b);
    Double_t dNprime = fCommonHists->GetEntriesInEtaBinPOI(b);   
    Double_t dErrdifcomb = 0.;  //set error to zero
    Double_t dErr2difcomb = 0.; //set error to zero
    //calculate error
    if (dNprime!=0.) { 
      for (Int_t theta=0;theta<iNtheta;theta++) {
	Double_t dTheta = ((double)theta/iNtheta)*TMath::Pi(); 
	Double_t dApluscomb = TMath::Exp((dJ01*dJ01)/(2*dChi*dChi)*
					 TMath::Cos(dTheta));
	Double_t dAmincomb = TMath::Exp(-(dJ01*dJ01)/(2*dChi*dChi)*
					TMath::Cos(dTheta));
	dErr2difcomb += (TMath::Cos(dTheta)/(4*dNprime*TMath::BesselJ1(dJ01)*
					     TMath::BesselJ1(dJ01)))*
	  ((dApluscomb*TMath::BesselJ0(2*dJ01*TMath::Sin(dTheta/2))) - 
	   (dAmincomb*TMath::BesselJ0(2*dJ01*TMath::Cos(dTheta/2))));
      } //loop over theta
    } 
      
    if (dErr2difcomb!=0.) {
      dErr2difcomb /= iNtheta;
      dErrdifcomb = TMath::Sqrt(dErr2difcomb);
    }
    else {dErrdifcomb = 0.;}
    //fill TH1D
    fCommonHistsRes->FillDifferentialFlowEtaPOI(b, dv2pro, dErrdifcomb); 
  } //loop over bins b
    
  //v as a function of Pt for RP selection
  TH1F* fHistPtRP = fCommonHists->GetHistPtRP(); //for calculating integrated flow
  Double_t dVRP = 0.;
  Double_t dSum = 0.;
  Double_t dErrV =0.;

  for(Int_t b=0;b<iNbinsPt;b++) {
    Double_t dv2pro  = fHistProVPtRP->GetBinContent(b);
    Double_t dNprime = fCommonHists->GetEntriesInPtBinRP(b);   
    Double_t dErrdifcomb = 0.;  //set error to zero
    Double_t dErr2difcomb = 0.; //set error to zero
    //calculate error
    if (dNprime!=0.) { 
      for (Int_t theta=0;theta<iNtheta;theta++) {
	Double_t dTheta = ((double)theta/iNtheta)*TMath::Pi(); 
	Double_t dApluscomb = TMath::Exp((dJ01*dJ01)/(2*dChi*dChi)*
					   TMath::Cos(dTheta));
	Double_t dAmincomb = TMath::Exp(-(dJ01*dJ01)/(2*dChi*dChi)*
					  TMath::Cos(dTheta));
	dErr2difcomb += (TMath::Cos(dTheta)/(4*dNprime*TMath::BesselJ1(dJ01)*
					       TMath::BesselJ1(dJ01)))*
	  ((dApluscomb*TMath::BesselJ0(2*dJ01*TMath::Sin(dTheta/2))) - 
	   (dAmincomb*TMath::BesselJ0(2*dJ01*TMath::Cos(dTheta/2))));
      } //loop over theta
    } 
      
    if (dErr2difcomb!=0.) {
      dErr2difcomb /= iNtheta;
      dErrdifcomb = TMath::Sqrt(dErr2difcomb);
      //cerr<<"dErrdifcomb = "<<dErrdifcomb<<endl;
    }
    else {dErrdifcomb = 0.;}
      
    //fill TH1D
    fCommonHistsRes->FillDifferentialFlowPtRP(b, dv2pro, dErrdifcomb);
    //calculate integrated flow for RP selection
    if (fHistPtRP){
      Double_t dYieldPt = fHistPtRP->GetBinContent(b);
      dVRP += dv2pro*dYieldPt;
      dSum +=dYieldPt;
      dErrV += dYieldPt*dYieldPt*dErrdifcomb*dErrdifcomb;
    } else { cout<<"fHistPtRP is NULL"<<endl; }
 
  } //loop over bins b

  if (TMath::AreEqualAbs(dSum, 0.0, 1e-10)) {
    dVRP /= dSum; //the pt distribution should be normalised
    dErrV /= (dSum*dSum);
    dErrV = TMath::Sqrt(dErrV);
  }
  fCommonHistsRes->FillIntegratedFlowRP(dVRP,dErrV);

  cout<<"dV(RP) = "<<dVRP<<" +- "<<dErrV<<endl;
  //cout<<endl;

       
  //v as a function of Pt for POI selection 
  TH1F* fHistPtPOI = fCommonHists->GetHistPtPOI(); //for calculating integrated flow
  Double_t dVPOI = 0.;
  dSum = 0.;
  dErrV =0.;
  
  for(Int_t b=0;b<iNbinsPt;b++) {
    Double_t dv2pro = fHistProVPtPOI->GetBinContent(b);
    Double_t dNprime = fCommonHists->GetEntriesInPtBinPOI(b);    
    
    //cerr<<"dNprime = "<<dNprime<<endl;
    //Int_t iNprime = TMath::Nint(fHistProVPtPOI->GetBinEntries(b));
    //cerr<<"iNprime = "<<iNprime<<endl;

    Double_t dErrdifcomb = 0.;  //set error to zero
    Double_t dErr2difcomb = 0.; //set error to zero
    //calculate error
    if (dNprime!=0.) { 
      for (Int_t theta=0;theta<iNtheta;theta++) {
	Double_t dTheta = ((double)theta/iNtheta)*TMath::Pi(); 
	Double_t dApluscomb = TMath::Exp((dJ01*dJ01)/(2*dChi*dChi)*
					   TMath::Cos(dTheta));
	Double_t dAmincomb = TMath::Exp(-(dJ01*dJ01)/(2*dChi*dChi)*
					  TMath::Cos(dTheta));
	dErr2difcomb += (TMath::Cos(dTheta)/(4*dNprime*TMath::BesselJ1(dJ01)*
						 TMath::BesselJ1(dJ01)))*
	  ((dApluscomb*TMath::BesselJ0(2*dJ01*TMath::Sin(dTheta/2))) - 
	   (dAmincomb*TMath::BesselJ0(2*dJ01*TMath::Cos(dTheta/2))));
      } //loop over theta
    } 
      
    if (dErr2difcomb!=0.) {
      dErr2difcomb /= iNtheta;
      dErrdifcomb = TMath::Sqrt(dErr2difcomb);
      //cerr<<"dErrdifcomb = "<<dErrdifcomb<<endl;
    }
    else {dErrdifcomb = 0.;}
	  
    //fill TH1D
    fCommonHistsRes->FillDifferentialFlowPtPOI(b, dv2pro, dErrdifcomb); 

    //calculate integrated flow for POI selection
    if (fHistPtPOI){
      Double_t dYieldPt = fHistPtPOI->GetBinContent(b);
      dVPOI += dv2pro*dYieldPt;
      dSum +=dYieldPt;
      dErrV += dYieldPt*dYieldPt*dErrdifcomb*dErrdifcomb;
    } else { cout<<"fHistPtPOI is NULL"<<endl; }
  } //loop over bins b

  if (dSum != 0.) {
    dVPOI /= dSum; //the pt distribution should be normalised
    dErrV /= (dSum*dSum);
    dErrV = TMath::Sqrt(dErrV);
  }
  fCommonHistsRes->FillIntegratedFlowPOI(dVPOI,dErrV);

  cout<<"dV(POI) = "<<dVPOI<<" +- "<<dErrV<<endl;
  cout<<endl;
  cout<<"*************************************"<<endl;
  cout<<"*************************************"<<endl;
    
  //cout<<".....finished"<<endl;
 }

