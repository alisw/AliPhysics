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
  fQsum(NULL),
  fQ2sum(0)
{

  // Constructor.
  fQsum = new TVector2();        // flow vector sum
}

 

 //-----------------------------------------------------------------------


 AliFlowAnalysisWithLYZEventPlane::~AliFlowAnalysisWithLYZEventPlane() 
 {
   //destructor
   delete fQsum;
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
  Int_t iNtheta = AliFlowLYZConstants::kTheta;
  Int_t iNbinsPt = AliFlowCommonConstants::GetNbinsPt();
  Double_t  dPtMin = AliFlowCommonConstants::GetPtMin();	     
  Double_t  dPtMax = AliFlowCommonConstants::GetPtMax();


  fHistProFlow = new TProfile("FlowPro_VPt_LYZEP","FlowPro_VPt_LYZEP",iNbinsPt,dPtMin,dPtMax);
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

  
  fHistProR0theta  = new TProfile("FlowPro_r0theta_LYZEP","FlowPro_r0theta_LYZEP",iNtheta,-0.5,iNtheta-0.5);
  fHistProR0theta->SetXTitle("#theta");
  fHistProR0theta->SetYTitle("r_{0}^{#theta}");

  fHistProReDtheta = new TProfile("FlowPro_ReDtheta_LYZEP","FlowPro_ReDtheta_LYZEP",iNtheta, -0.5, iNtheta-0.5);
  fHistProReDtheta->SetXTitle("#theta");
  fHistProReDtheta->SetYTitle("Re(D^{#theta})");

  fHistProImDtheta = new TProfile("FlowPro_ImDtheta_LYZEP","FlowPro_ImDtheta_LYZEP",iNtheta, -0.5, iNtheta-0.5);
  fHistProImDtheta->SetXTitle("#theta");
  fHistProImDtheta->SetYTitle("Im(D^{#theta})");

  fEventNumber = 0;  //set number of events to zero

      
} 
 
//-----------------------------------------------------------------------
 
void AliFlowAnalysisWithLYZEventPlane::Make(AliFlowEventSimple* anEvent, AliFlowLYZEventPlane* aLYZEP) {
  
  //Get the event plane and weight for each event
  if (anEvent) {
         
    //fill control histograms     
    fCommonHists->FillControlHistograms(anEvent);

    //get the Q vector from the FlowEvent
    AliFlowVector vQ = anEvent->GetQ(); 
    //Weight with the multiplicity
    Double_t dQX = 0.;
    Double_t dQY = 0.;
    if (vQ.GetMult()!=0.) {
      dQX = vQ.X()/vQ.GetMult();
      dQY = vQ.Y()/vQ.GetMult();
    } else {cerr<<"vQ.GetMult() is zero!"<<endl; }
    vQ.Set(dQX,dQY);
    //for chi calculation:
    *fQsum += vQ;
    fQ2sum += vQ.Mod2();
    cout<<"fQ2sum = "<<fQ2sum<<endl;

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
       
    //calculate flow
    //loop over the tracks of the event
    Int_t iNumberOfTracks = anEvent->NumberOfTracks(); 
    for (Int_t i=0;i<iNumberOfTracks;i++) 
      {
	AliFlowTrackSimple* pTrack = anEvent->GetTrack(i) ; 
	if (pTrack){
	  if (pTrack->UseForDifferentialFlow()) {
	    Double_t dPhi = pTrack->Phi();
	    //if (dPhi<0.) fPhi+=2*TMath::Pi();
	    //calculate flow v2:
	    Double_t dv2 = dWR * TMath::Cos(2*(dPhi-dRP));
	    Double_t dPt = pTrack->Pt();
	    //fill histogram
	    fHistProFlow->Fill(dPt,100*dv2);  
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
  Double_t  dJ01 = 2.405; 
  Int_t iNtheta = AliFlowLYZConstants::kTheta;
  Int_t iNbinsPt = AliFlowCommonConstants::GetNbinsPt();
    
  //calculate dV the mean of dVtheta
  Double_t  dVtheta = 0; 
  Double_t  dV = 0; 
  for (Int_t theta=0;theta<iNtheta;theta++)	{
    Double_t dR0 = fFirstr0theta->GetBinContent(theta+1); 
    if (dR0!=0.) { dVtheta = dJ01/dR0 ;}
    dV += dVtheta;
  }
  dV /= iNtheta;

  //calculate dChi 
  Double_t  dSigma2 = 0;
  Double_t  dChi= 0;
  if (fEventNumber!=0) {
    *fQsum /= fEventNumber;
    //cerr<<"fQsum.X() = "<<fQsum.X()<<endl;
    //cerr<<"fQsum.Y() = "<<fQsum.Y()<<endl;
    fQ2sum /= fEventNumber;
    cerr<<"fEventNumber = "<<fEventNumber<<endl;
    cerr<<"fQ2sum = "<<fQ2sum<<endl;
    dSigma2 = fQ2sum - TMath::Power(fQsum->X(),2.) - TMath::Power(fQsum->Y(),2.) - TMath::Power(dV,2.);  //BP eq. 62
    if (dSigma2>0) dChi = dV/TMath::Sqrt(dSigma2);
    else dChi = -1.;
    fCommonHistsRes->FillChi(dChi);
    cerr<<"dV = "<<dV<<" and chi = "<<dChi<<endl;
  }
  
  for(Int_t b=0;b<iNbinsPt;b++){
    Double_t dv2pro = 0.;
    Double_t dErr2difcomb = 0.;   
    Double_t dErrdifcomb = 0.;
    if(fHistProFlow) {
      dv2pro = fHistProFlow->GetBinContent(b);
      //calculate error
      for (Int_t theta=0;theta<iNtheta;theta++) {
	Double_t dTheta = ((double)theta/iNtheta)*TMath::Pi(); 
	Int_t iNprime = TMath::Nint(fHistProFlow->GetBinEntries(b));
	//cerr<<"iNprime = "<<iNprime<<endl;
	if (iNprime!=0) { 
	  Double_t dApluscomb = TMath::Exp((dJ01*dJ01)/(2*dChi*dChi)*
					   TMath::Cos(dTheta));
	  Double_t dAmincomb = TMath::Exp(-(dJ01*dJ01)/(2*dChi*dChi)*
					  TMath::Cos(dTheta));
	  dErr2difcomb += (TMath::Cos(dTheta)/(4*iNprime*TMath::BesselJ1(dJ01)*
						 TMath::BesselJ1(dJ01)))*
	    ((dApluscomb*TMath::BesselJ0(2*dJ01*TMath::Sin(dTheta/2))) - 
	     (dAmincomb*TMath::BesselJ0(2*dJ01*TMath::Cos(dTheta/2))));
	} //if !=0
	//else { cout<<"iNprime = 0"<<endl; }
      } //loop over theta
      
      if (dErr2difcomb!=0.) {
	dErr2difcomb /= iNtheta;
	dErrdifcomb = TMath::Sqrt(dErr2difcomb)*100;
	//cerr<<"dErrdifcomb = "<<dErrdifcomb<<endl;
      }
      else {dErrdifcomb = 0.; }

      //fill TH1D
      fCommonHistsRes->FillDifferentialFlow(b, dv2pro, dErrdifcomb); 
    
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

