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

/*
$Log$
*/ 

//#define AliFlowAnalysisWithLYZEventPlane_cxx
 
#include "Riostream.h"  //needed as include

#include "TFile.h"
#include "TComplex.h"   //needed as include
#include "TCanvas.h"   //needed as include
#include "TLegend.h"   //needed as include
#include "TProfile.h"  //needed as include

class TH1F;

#include "AliFlowLYZConstants.h"    //needed as include
#include "AliFlowCommonConstants.h" //needed as include
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "AliFlowLYZEventPlane.h"
#include "AliFlowAnalysisWithLYZEventPlane.h"

class AliFlowVector;

// AliFlowAnalysisWithLYZEventPlane:
//
// Class to do flow analysis with the event plane from the LYZ method
//
// author: N. van der Kolk (kolk@nikhef.nl)


ClassImp(AliFlowAnalysisWithLYZEventPlane)

  //-----------------------------------------------------------------------
 
 AliFlowAnalysisWithLYZEventPlane::AliFlowAnalysisWithLYZEventPlane():
   fOutFile(0), 
   fFirstRunFile(0),
   fSecondRunFile(0),
   fFirstRunFileName(0),
   fSecondRunFileName(0),
   fOutFileName(0),
   fSecondReDtheta(0),
   fSecondImDtheta(0),
   fFirstr0theta(0),
   fSecondVPt(0),
   fHistProFlow(0),
   fHistProFlow2(0),
   fHistProWr(0),
   fHistProWrCorr(0),
   fHistFlow(0),
   fHistDeltaPhi(0),
   fHistDeltaPhi2(0),
   fHistDeltaPhihere(0),
   fHistPhiEP(0),
   fHistPhiEPhere(0),
   fHistPhiLYZ(0),
   fHistPhiLYZ2(0),
   fHistProR0theta(0),
   fHistProReDtheta(0),
   fHistProImDtheta(0),
   fCommonHists(0),
   fCommonHistsRes(0),
   fEventNumber(0),
   fQX(0),
   fQY(0),
   fQ2sum(0),
   fQtheta(0),
   fEvent(0x0),
   fTrack(0x0),
   fLYZEP(0x0)
{

  // Constructor.
  fQ.Set(0.,0.);           // flow vector
  fQsum.Set(0.,0.);        // flow vector sum
  fLYZEP = new AliFlowLYZEventPlane();
}

 

 //-----------------------------------------------------------------------


 AliFlowAnalysisWithLYZEventPlane::~AliFlowAnalysisWithLYZEventPlane() 
 {
   //destructor
   delete fLYZEP;
 }
 

//-----------------------------------------------------------------------
void AliFlowAnalysisWithLYZEventPlane::Init() {

  //Initialise all histograms
  cout<<"---Analysis with Lee Yang Zeros Event Plane Method---"<<endl;

  //input histograms
  if (fSecondRunFile->IsZombie()){ //check if file exists
    cout << "Error opening file, run first with fFirstrun = kTRUE" << endl;
    exit(-1);
  } else if (fSecondRunFile->IsOpen()){
    cout<<"----secondRunFile is open----"<<endl;
    fSecondVPt = (TProfile*)fSecondRunFile->Get("Flow_Differential_Pt_LYZ"); //to compare to
  }

  if (fFirstRunFile->IsZombie()){ //check if file exists
    cout << "Error opening file, run first with fFirstrun = kTRUE" << endl;
    exit(-1);
  } else if (fFirstRunFile->IsOpen()){
    cout<<"----firstRunFile is open----"<<endl<<endl;
    fFirstr0theta = (TProfile*)fFirstRunFile->Get("First_FlowPro_r0theta_LYZ");
  }


  //output file
  // ********make output file 
  // analysis file (output)
  fOutFile = new TFile(fOutFileName.Data(),"RECREATE") ;

  fCommonHists = new AliFlowCommonHist("LYZEP");
  fCommonHistsRes = new AliFlowCommonHistResults("LYZEP");
  
  // output histograms
  Int_t fNtheta = AliFlowLYZConstants::kTheta;
  Int_t fNbinsPt = AliFlowCommonConstants::GetNbinsPt();
  Double_t  fPtMin = AliFlowCommonConstants::GetPtMin();	     
  Double_t  fPtMax = AliFlowCommonConstants::GetPtMax();


  fHistProFlow = new TProfile("FlowPro_VPt_LYZEP","FlowPro_VPt_LYZEP",fNbinsPt,fPtMin,fPtMax);
  fHistProFlow->SetXTitle("Pt");
  fHistProFlow->SetYTitle("v2 (%)");
  
  fHistProWr = new TProfile("FlowPro_Wr_LYZEP","FlowPro_Wr_LYZEP",100,0.,0.25);
  fHistProWr->SetXTitle("Q");
  fHistProWr->SetYTitle("Wr");

  fHistDeltaPhi = new TH1F("Flow_DeltaPhi_LYZEP","Flow_DeltaPhi_LYZEP",100,0.,3.14);
  fHistDeltaPhi->SetXTitle("DeltaPhi");
  fHistDeltaPhi->SetYTitle("Counts");

  fHistPhiLYZ = new TH1F("Flow_PhiLYZ_LYZEP","Flow_PhiLYZ_LYZEP",100,0.,3.14);
  fHistPhiLYZ->SetXTitle("Phi from LYZ");
  fHistPhiLYZ->SetYTitle("Counts");

  fHistPhiEP = new TH1F("Flow_PhiEP_LYZEP","Flow_PhiEP_LYZEP",100,0.,3.14);
  fHistPhiEP->SetXTitle("Phi from EP");
  fHistPhiEP->SetYTitle("Counts");

  
  fHistProR0theta  = new TProfile("FlowPro_r0theta_LYZEP","FlowPro_r0theta_LYZEP",fNtheta,-0.5,fNtheta-0.5);
  fHistProR0theta->SetXTitle("#theta");
  fHistProR0theta->SetYTitle("r_{0}^{#theta}");

  fHistProReDtheta = new TProfile("FlowPro_ReDtheta_LYZEP","FlowPro_ReDtheta_LYZEP",fNtheta, -0.5, fNtheta-0.5);
  fHistProReDtheta->SetXTitle("#theta");
  fHistProReDtheta->SetYTitle("Re(D^{#theta})");

  fHistProImDtheta = new TProfile("FlowPro_ImDtheta_LYZEP","FlowPro_ImDtheta_LYZEP",fNtheta, -0.5, fNtheta-0.5);
  fHistProImDtheta->SetXTitle("#theta");
  fHistProImDtheta->SetYTitle("Im(D^{#theta})");

  fEventNumber = 0;  //set number of events to zero

      
} 
 
//-----------------------------------------------------------------------
 
void AliFlowAnalysisWithLYZEventPlane::Make(AliFlowEventSimple* fEvent, AliFlowLYZEventPlane* fLYZEP) {

  //Get the event plane and weight for each event
  if (fEvent) {
         
    //fill control histograms     
    fCommonHists->FillControlHistograms(fEvent);

    //get the Q vector from the FlowEvent
    fQ = fEvent->GetQ(); 
    //Weight with the multiplicity
    Double_t fQX = 0.;
    Double_t fQY = 0.;
    if (fQ.GetMult()!=0.) {
      fQX = fQ.X()/fQ.GetMult();
      fQY = fQ.Y()/fQ.GetMult();
    } else {cerr<<"fQ.GetMult() is zero!"<<endl; }
    fQ.Set(fQX,fQY);
    //cout<<"fQ.Mod() = " << fQ.Mod() << endl;
    //for chi calculation:
    fQsum += fQ;
    //cout<<"fQsum.Mod() = "<<fQsum.Mod()<<endl;
    fQ2sum += fQ.Mod2();
    cout<<"fQ2sum = "<<fQ2sum<<endl;

    //call AliFlowLYZEventPlane::CalculateRPandW() here!
    fLYZEP->CalculateRPandW(fQ);

    Double_t fWR = fLYZEP->GetWR();     
    Double_t fRP = fLYZEP->GetPsi();

    //fHistProWr->Fill(fQ.Mod(),fWR); //this weight is always positive
    fHistPhiLYZ->Fill(fRP);   
    
    //plot difference between event plane from EP-method and LYZ-method
    Double_t fRPEP = fQ.Phi()/2;                              //gives distribution from (0 to pi)
    //Float_t fRPEP = 0.5*TMath::ATan2(fQ.Y(),fQ.X());       //gives distribution from (-pi/2 to pi/2)
    //cout<<"fRPEP = "<< fRPEP <<endl;
    fHistPhiEP->Fill(fRPEP);

    Double_t fDeltaPhi = fRPEP - fRP;
    if (fDeltaPhi < 0.) { fDeltaPhi += TMath::Pi(); }        //to shift distribution from (-pi/2 to pi/2) to (0 to pi)
    //cout<<"fDeltaPhi = "<<fDeltaPhi<<endl;
    fHistDeltaPhi->Fill(fDeltaPhi); 

    //Flip sign of WR
    Double_t low = TMath::Pi()/4.;
    Double_t high = 3.*(TMath::Pi()/4.);
    if ((fDeltaPhi > low) && (fDeltaPhi < high)){
      fRP -= (TMath::Pi()/2);
      fWR = -fWR;
      cerr<<"*** fRP modified ***"<<endl;
    }
    fHistProWr->Fill(fQ.Mod(),fWR); //corrected weight
       
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
	    Double_t fv2 = fWR * TMath::Cos(2*(fPhi-fRP));
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
void AliFlowAnalysisWithLYZEventPlane::Finish() {
   
  //plot histograms etc. 
  cout<<"AliFlowAnalysisWithLYZEventPlane::Terminate()"<<endl;
  //constands:
  Double_t  fJ01 = 2.405; 
  Int_t fNtheta = AliFlowLYZConstants::kTheta;
  Int_t fNbinsPt = AliFlowCommonConstants::GetNbinsPt();
  //Double_t fErr2difcomb = 0.;
  //Double_t fErrdifcomb = 0.;
  
  //calculate fV the mean of fVtheta
  Double_t  fVtheta = 0; 
  Double_t  fV = 0; 
  for (Int_t theta=0;theta<fNtheta;theta++)	{
    Double_t fR0 = fFirstr0theta->GetBinContent(theta+1); 
    if (fR0!=0.) { fVtheta = fJ01/fR0 ;}
    fV += fVtheta;
  }
  fV /= fNtheta;

  //calculate fChi 
  Double_t  fSigma2 = 0;
  Double_t  fChi= 0;
  if (fEventNumber!=0) {
    fQsum /= fEventNumber;
    //cerr<<"fQsum.X() = "<<fQsum.X()<<endl;
    //cerr<<"fQsum.Y() = "<<fQsum.Y()<<endl;
    fQ2sum /= fEventNumber;
    cerr<<"fEventNumber = "<<fEventNumber<<endl;
    cerr<<"fQ2sum = "<<fQ2sum<<endl;
    fSigma2 = fQ2sum - TMath::Power(fQsum.X(),2.) - TMath::Power(fQsum.Y(),2.) - TMath::Power(fV,2.);  //BP eq. 62
    if (fSigma2>0) fChi = fV/TMath::Sqrt(fSigma2);
    else fChi = -1.;
    fCommonHistsRes->FillChi(fChi);
    cerr<<"fV = "<<fV<<" and chi = "<<fChi<<endl;
  }
  
  for(Int_t b=0;b<fNbinsPt;b++){
    Double_t fv2pro = 0.;
    Double_t fErr2difcomb = 0.;   
    Double_t fErrdifcomb = 0.;
    if(fHistProFlow) {
      fv2pro = fHistProFlow->GetBinContent(b);
      //calculate error
      for (Int_t theta=0;theta<fNtheta;theta++) {
	Double_t fTheta = ((double)theta/fNtheta)*TMath::Pi(); 
	Int_t fNprime = TMath::Nint(fHistProFlow->GetBinEntries(b));
	//cerr<<"fNprime = "<<fNprime<<endl;
	if (fNprime!=0.) { 
	  Double_t fApluscomb = TMath::Exp((fJ01*fJ01)/(2*fChi*fChi)*
					   TMath::Cos(fTheta));
	  Double_t fAmincomb = TMath::Exp(-(fJ01*fJ01)/(2*fChi*fChi)*
					  TMath::Cos(fTheta));
	  fErr2difcomb += (TMath::Cos(fTheta)/(4*fNprime*TMath::BesselJ1(fJ01)*
						 TMath::BesselJ1(fJ01)))*
	    ((fApluscomb*TMath::BesselJ0(2*fJ01*TMath::Sin(fTheta/2))) - 
	     (fAmincomb*TMath::BesselJ0(2*fJ01*TMath::Cos(fTheta/2))));
	} //if !=0
	//else { cout<<"fNprime = 0."<<endl; }
      } //loop over theta
      
      if (fErr2difcomb!=0.) {
	fErr2difcomb /= fNtheta;
	fErrdifcomb = TMath::Sqrt(fErr2difcomb)*100;
	//cerr<<"fErrdifcomb = "<<fErrdifcomb<<endl;
      }
      else {fErrdifcomb = 0.; }

      //fill TH1D
      fCommonHistsRes->FillDifferentialFlow(b, fv2pro, fErrdifcomb); 
    
    } //if fHistProFLow
    else  {
      cout << "Profile Hist missing" << endl;
      break;
    }
    
  } //loop over b

  // write to file
  fOutFile->Write();

  cout<<"Making some plots to check the results:"<<endl<<endl;

  TCanvas *canvas = new TCanvas("canvas","compare v2 vs pt",800,800);
  canvas->cd();
  fSecondVPt->SetLineColor(3);
  fSecondVPt->SetLineWidth(2);
  fSecondVPt->Draw();
  fHistFlow = fCommonHistsRes->GetfHistDiffFlow();
  fHistFlow->Draw("SAME");
  // draw the legend
  TLegend *legend2 = new TLegend(0.6,0.65,0.88,0.85);
  legend2->SetTextFont(72);
  legend2->SetTextSize(0.04);
  legend2->AddEntry(fSecondVPt,"stand. LYZ","lpe");
  legend2->AddEntry(fHistFlow,"new LYZ with calculated errors","lpe");
  legend2->Draw();
   
  TCanvas *canvas5 = new TCanvas("canvas5","phi and Delta Phi",800,600);
  canvas5->Divide(3,1);
  canvas5->cd(1); 
  fHistDeltaPhi->Draw();
  canvas5->cd(2); 
  fHistPhiLYZ->Draw();
  canvas5->cd(3);
  fHistPhiEP->Draw();
  	  
  cout<<".....finished"<<endl;
 }
