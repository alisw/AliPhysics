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


#include "Riostream.h"                 //needed as include
#include "TObject.h"                   //needed as include
#include "AliFlowCommonConstants.h"    //needed as include
#include "AliFlowLYZConstants.h"       //needed as include
#include "AliFlowAnalysisWithLeeYangZeros.h"
#include "AliFlowLYZHist1.h"           //needed as include
#include "AliFlowLYZHist2.h"           //needed as include
#include "AliFlowCommonHist.h"         //needed as include
#include "AliFlowCommonHistResults.h"  //needed as include
#include "AliFlowEventSimple.h"        //needed as include
#include "AliFlowTrackSimple.h"        //needed as include

class AliFlowVector;

#include "TMath.h" //needed as include
#include "TFile.h" //needed as include
#include "TList.h"

class TComplex;
class TProfile;
class TH1F;
class TH1D;


//Description: Maker to analyze Flow using the LeeYangZeros method  
//             Equation numbers are from Big Paper (BP): Nucl. Phys. A 727, 373 (2003)
//             Practical Guide (PG):    J. Phys. G: Nucl. Part. Phys. 30, S1213 (2004)  
//             Adapted from StFlowLeeYangZerosMaker.cxx           
//             by Markus Oldenberg and Art Poskanzer, LBNL        
//             with advice from Jean-Yves Ollitrault and Nicolas Borghini   
//
//Author: Naomi van der Kolk (kolk@nikhef.nl)



ClassImp(AliFlowAnalysisWithLeeYangZeros)

  //-----------------------------------------------------------------------
 
  AliFlowAnalysisWithLeeYangZeros::AliFlowAnalysisWithLeeYangZeros():
    fQsum(NULL),
    fQ2sum(0),
    fEventNumber(0),
    fFirstRun(kTRUE),
    fUseSum(kTRUE),
    fDoubleLoop(kFALSE),
    fDebug(kFALSE),
    fHistList(NULL),
    fFirstRunList(NULL),
    fHistProVtheta(NULL),
    fHistProVetaRP(NULL),  
    fHistProVetaPOI(NULL), 
    fHistProVPtRP(NULL),   
    fHistProVPtPOI(NULL),  
    fHistProR0theta(NULL),
    fHistProReDenom(NULL),
    fHistProImDenom(NULL),
    fHistProReDtheta(NULL),
    fHistProImDtheta(NULL),
    fHistQsumforChi(NULL),
    fCommonHists(NULL),
    fCommonHistsRes(NULL)
  
{
  //default constructor
  if (fDebug) cout<<"****AliFlowAnalysisWithLeeYangZeros::AliFlowAnalysisWithLeeYangZeros default constructor****"<<endl;

  fHistList = new TList();
  fFirstRunList = new TList();

  for(Int_t i = 0;i<5;i++)
    {
      fHist1[i]=0;
      fHist2RP[i]=0;
      fHist2POI[i]=0;
    }

  fQsum = new TVector2();

}

//-----------------------------------------------------------------------


 AliFlowAnalysisWithLeeYangZeros::~AliFlowAnalysisWithLeeYangZeros() 
 {
   //default destructor
   if (fDebug) cout<<"****~AliFlowAnalysisWithLeeYangZeros****"<<endl;
   delete fQsum;
   delete fHistList;
   delete fFirstRunList;
   
 }
 
//-----------------------------------------------------------------------

void AliFlowAnalysisWithLeeYangZeros::WriteHistograms(TString* outputFileName)
{
 //store the final results in output .root file

  TFile *output = new TFile(outputFileName->Data(),"RECREATE");
  if (GetFirstRun()) {
    output->WriteObject(fHistList, "cobjLYZ1","SingleKey");
  }
  else {
    output->WriteObject(fHistList, "cobjLYZ2","SingleKey");
  }
  delete output;
}

//-----------------------------------------------------------------------

Bool_t AliFlowAnalysisWithLeeYangZeros::Init() 
{
  //init method 
  if (fDebug) cout<<"****AliFlowAnalysisWithLeeYangZeros::Init()****"<<endl;

  // Book histograms
  Int_t iNtheta = AliFlowLYZConstants::kTheta;
  Int_t iNbinsPt = AliFlowCommonConstants::GetNbinsPt();
  Int_t iNbinsEta = AliFlowCommonConstants::GetNbinsEta();

  Double_t  dPtMin = AliFlowCommonConstants::GetPtMin();	     
  Double_t  dPtMax = AliFlowCommonConstants::GetPtMax();
  Double_t  dEtaMin = AliFlowCommonConstants::GetEtaMin();	     
  Double_t  dEtaMax = AliFlowCommonConstants::GetEtaMax();
    
  //for control histograms
  if (fFirstRun){ fCommonHists = new AliFlowCommonHist("AliFlowCommonHistLYZ1");
  fHistList->Add(fCommonHists);
  fCommonHistsRes = new AliFlowCommonHistResults("AliFlowCommonHistResultsLYZ1");
  fHistList->Add(fCommonHistsRes);
  }
  else { fCommonHists = new AliFlowCommonHist("AliFlowCommonHistLYZ2");
  fHistList->Add(fCommonHists);
  fCommonHistsRes = new AliFlowCommonHistResults("AliFlowCommonHistResultsLYZ2");
  fHistList->Add(fCommonHistsRes); 
  }
    
  fHistQsumforChi = new TH1F("Flow_QsumforChi_LYZ","Flow_QsumforChi_LYZ",3,-1.,2.);
  fHistQsumforChi->SetXTitle("Qsum.X , Qsum.Y, Q2sum");
  fHistQsumforChi->SetYTitle("value");
  fHistList->Add(fHistQsumforChi);

  //for first loop over events 
  if (fFirstRun){
    fHistProR0theta  = new TProfile("First_FlowPro_r0theta_LYZ","First_FlowPro_r0theta_LYZ",iNtheta,-0.5,iNtheta-0.5);
    fHistProR0theta->SetXTitle("#theta");
    fHistProR0theta->SetYTitle("r_{0}^{#theta}");
    fHistList->Add(fHistProR0theta);

    fHistProVtheta  = new TProfile("First_FlowPro_Vtheta_LYZ","First_FlowPro_Vtheta_LYZ",iNtheta,-0.5,iNtheta-0.5);
    fHistProVtheta->SetXTitle("#theta");
    fHistProVtheta->SetYTitle("V_{n}^{#theta}");        
    fHistList->Add(fHistProVtheta);

    //class AliFlowLYZHist1 defines the histograms: fHistProGtheta, fHistProReGtheta, fHistProImGtheta, fHistProR0theta
    for (Int_t theta=0;theta<iNtheta;theta++) {  
      TString name = "AliFlowLYZHist1_";
      name += theta;
      fHist1[theta]=new AliFlowLYZHist1(theta, name);
      fHistList->Add(fHist1[theta]);
    }
         
  }
  //for second loop over events 
  else {
    fHistProReDenom = new TProfile("Second_FlowPro_ReDenom_LYZ","Second_FlowPro_ReDenom_LYZ" , iNtheta, -0.5, iNtheta-0.5);
    fHistProReDenom->SetXTitle("#theta");
    fHistProReDenom->SetYTitle("Re(Q^{#theta}e^{ir_{0}^{#theta}Q^{#theta}})");
    fHistList->Add(fHistProReDenom);

    fHistProImDenom = new TProfile("Second_FlowPro_ImDenom_LYZ","Second_FlowPro_ImDenom_LYZ" , iNtheta, -0.5, iNtheta-0.5);
    fHistProImDenom->SetXTitle("#theta");
    fHistProImDenom->SetYTitle("Im(Q^{#theta}e^{ir_{0}^{#theta}Q^{#theta}})");
    fHistList->Add(fHistProImDenom);

    fHistProVetaRP = new TProfile("Second_FlowPro_VetaRP_LYZ","Second_FlowPro_VetaRP_LYZ",iNbinsEta,dEtaMin,dEtaMax);
    fHistProVetaRP->SetXTitle("rapidity");
    fHistProVetaRP->SetYTitle("v_{2}(#eta) for RP selection");
    fHistList->Add(fHistProVetaRP);

    fHistProVetaPOI = new TProfile("Second_FlowPro_VetaPOI_LYZ","Second_FlowPro_VetaPOI_LYZ",iNbinsEta,dEtaMin,dEtaMax);
    fHistProVetaPOI->SetXTitle("rapidity");
    fHistProVetaPOI->SetYTitle("v_{2}(#eta) for POI selection");
    fHistList->Add(fHistProVetaPOI);

    fHistProVPtRP = new TProfile("Second_FlowPro_VPtRP_LYZ","Second_FlowPro_VPtRP_LYZ",iNbinsPt,dPtMin,dPtMax);
    fHistProVPtRP->SetXTitle("Pt");
    fHistProVPtRP->SetYTitle("v_{2}(p_{T}) for RP selection");
    fHistList->Add(fHistProVPtRP);

    fHistProVPtPOI = new TProfile("Second_FlowPro_VPtPOI_LYZ","Second_FlowPro_VPtPOI_LYZ",iNbinsPt,dPtMin,dPtMax);
    fHistProVPtPOI->SetXTitle("p_{T}");
    fHistProVPtPOI->SetYTitle("v_{2}(p_{T}) for POI selection");
    fHistList->Add(fHistProVPtPOI);

    fHistProReDtheta = new TProfile("Second_FlowPro_ReDtheta_LYZ","Second_FlowPro_ReDtheta_LYZ",iNtheta, -0.5, iNtheta-0.5);
    fHistProReDtheta->SetXTitle("#theta");
    fHistProReDtheta->SetYTitle("Re(D^{#theta})");
    fHistList->Add(fHistProReDtheta);

    fHistProImDtheta = new TProfile("Second_FlowPro_ImDtheta_LYZ","Second_FlowPro_ImDtheta_LYZ",iNtheta, -0.5, iNtheta-0.5);
    fHistProImDtheta->SetXTitle("#theta");
    fHistProImDtheta->SetYTitle("Im(D^{#theta})");
    fHistList->Add(fHistProImDtheta);

    //class AliFlowLYZHist2 defines the histograms: 
    for (Int_t theta=0;theta<iNtheta;theta++)  {  
      TString nameRP = "AliFlowLYZHist2RP_";
      nameRP += theta;
      fHist2RP[theta]=new AliFlowLYZHist2(theta,nameRP);
      fHistList->Add(fHist2RP[theta]);

      TString namePOI = "AliFlowLYZHist2POI_";
      namePOI += theta;
      fHist2POI[theta]=new AliFlowLYZHist2(theta,namePOI);
      fHistList->Add(fHist2POI[theta]);
    }
     
    //read histogram fHistProR0theta from the first run list
    if (fFirstRunList) {
      fHistProR0theta  = (TProfile*)fFirstRunList->FindObject("First_FlowPro_r0theta_LYZ");
      if (!fHistProR0theta) {cout<<"fHistProR0theta has a NULL pointer!"<<endl;}
      fHistList->Add(fHistProR0theta);
    } else { cout<<"list is NULL pointer!"<<endl; }

  }
   

  if (fDebug) cout<<"****Histograms initialised****"<<endl;
    
  fEventNumber = 0; //set event counter to zero
 
  return kTRUE; 
}
 
 //-----------------------------------------------------------------------
 
Bool_t AliFlowAnalysisWithLeeYangZeros::Make(AliFlowEventSimple* anEvent) 
{
  //make method
  if (fDebug) cout<<"****AliFlowAnalysisWithLeeYangZeros::Make()****"<<endl;
        
  //get tracks from event
  if (anEvent) {
    if (fFirstRun){
      fCommonHists->FillControlHistograms(anEvent);
      FillFromFlowEvent(anEvent);
    }
    else {
      fCommonHists->FillControlHistograms(anEvent);
      SecondFillFromFlowEvent(anEvent);
    }
  }
 
  else {
    cout<<"##### FlowLeeYangZero: Stack pointer null"<<endl;
    return kFALSE;
  }
  cout<<"^^^^read event "<<fEventNumber<<endl;
  fEventNumber++;
  
     
  return kTRUE; 
}


  //-----------------------------------------------------------------------     
 Bool_t AliFlowAnalysisWithLeeYangZeros::Finish() 
{
  //finish method
  if (fDebug) cout<<"****AliFlowAnalysisWithLeeYangZeros::Finish()****"<<endl; 
  
  //define variables for both runs
  Double_t  dJ01 = 2.405; 
  Int_t iNtheta = AliFlowLYZConstants::kTheta;
  //set the event number
  SetEventNumber((int)fCommonHists->GetHistMultOrig()->GetEntries());
  //cout<<"number of events processed is "<<fEventNumber<<endl; 

  //set the sum of Q vectors
  fQsum->Set(fHistQsumforChi->GetBinContent(1),fHistQsumforChi->GetBinContent(2));
  SetQ2sum(fHistQsumforChi->GetBinContent(3));  
    
  if (fFirstRun){
    Double_t  dR0 = 0;
    Double_t  dVtheta = 0; 
    Double_t  dv = 0;
    Double_t  dV = 0; 
    for (Int_t theta=0;theta<iNtheta;theta++)
      {
	//get the first minimum r0
	fHist1[theta]->FillGtheta();
	dR0 = fHist1[theta]->GetR0();
	    	   	   
	//calculate integrated flow
	if (dR0!=0) dVtheta = dJ01/dR0;
	//for estimating systematic error resulting from d0
	Double_t dBinsize = (AliFlowLYZConstants::fgMax)/(AliFlowLYZConstants::kNbins);
	Double_t dVplus = dJ01/(dR0+dBinsize);
	Double_t dVmin = dJ01/(dR0-dBinsize);
	dv = dVtheta;
	Double_t dvplus = dVplus;
	Double_t dvmin = dVmin;
	if (fDebug) cout<<"dv = "<<dv<<" and dvplus = "<<dvplus<< " and dvmin = "<<dvmin<<endl;
	     
	//fill the histograms
	fHistProR0theta->Fill(theta,dR0);
	fHistProVtheta->Fill(theta,dVtheta); 
	//get average value of fVtheta = fV
	dV += dVtheta;
      } //end of loop over theta

    //get average value of fVtheta = fV
    dV /=iNtheta;
       
    //sigma2 and chi 
    Double_t  dSigma2 = 0;
    Double_t  dChi= 0;
    if (fEventNumber!=0) {
      *fQsum /= fEventNumber;
      fQ2sum /= fEventNumber;
      dSigma2 = fQ2sum - TMath::Power(fQsum->X(),2.) - TMath::Power(fQsum->Y(),2.) - TMath::Power(dV,2.);  //BP eq. 62
      if (dSigma2>0) dChi = dV/TMath::Sqrt(dSigma2);
      else dChi = -1.;
      fCommonHistsRes->FillChiRP(dChi);
  
      cout<<"*************************************"<<endl;
      cout<<"*************************************"<<endl;
      cout<<"      Integrated flow from           "<<endl;
      cout<<"        Lee-Yang Zeroes              "<<endl;
      cout<<endl;
      cout<<"dV(RP) = "<<dV<<" and chi = "<<dChi<<endl;
      cout<<endl;
    }
	   
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

    //copy content of profile into TH1D and add error
    Double_t dv2pro = dV;   //in the case that fv is equal to fV
    Double_t dv2Err = dv2pro*dRelErrcomb ; 
    cout<<"dV(RP) = "<<dv2pro<<" +- "<<dv2Err<<endl;
    cout<<endl;
    cout<<"*************************************"<<endl;
    cout<<"*************************************"<<endl;
    fCommonHistsRes->FillIntegratedFlowRP(dv2pro, dv2Err);  


    if (fDebug) cout<<"****histograms filled****"<<endl;  
    
    return kTRUE;
    fEventNumber =0; //set to zero for second round over events
  }  //firstrun
   
 
  else {  //second run
   
    //calculate differential flow
    //declare variables
    TComplex cDenom, cNumerRP, cNumerPOI, cDtheta;
    Int_t m = 1;
    TComplex i = TComplex::I();
    Double_t dBesselRatio[3] = {1., 1.202, 2.69};
    Int_t iNbinsPt = AliFlowCommonConstants::GetNbinsPt();
    Int_t iNbinsEta = AliFlowCommonConstants::GetNbinsEta();

    Double_t dEtaRP, dPtRP, dReRatioRP, dVetaRP, dVPtRP, dEtaPOI, dPtPOI, dReRatioPOI, dVetaPOI, dVPtPOI;
         
    Double_t dR0 = 0.; 
    Double_t dVtheta = 0.;
    Double_t dV = 0.;
    Double_t dReDenom = 0.;
    Double_t dImDenom = 0.; 
    for (Int_t theta=0;theta<iNtheta;theta++)  { //loop over theta
      if (!fHistProR0theta) {
	cout << "Hist pointer R0theta in file does not exist" <<endl;
      }	else {
	dR0 = fHistProR0theta->GetBinContent(theta+1); //histogram starts at bin 1
	if (fDebug) cerr<<"dR0 = "<<dR0<<endl;
	if (dR0!=0) dVtheta = dJ01/dR0;
	dV += dVtheta; 
	// BP Eq. 9  -> Vn^theta = j01/r0^theta
	if (!fHistProReDenom && !fHistProImDenom) {
	  cout << "Hist pointer fDenom in file does not exist" <<endl;
	} else {
	  dReDenom = fHistProReDenom->GetBinContent(theta+1);
	  dImDenom = fHistProImDenom->GetBinContent(theta+1);
	  cDenom(dReDenom,dImDenom);
      
	  //for new method and use by others (only with the sum generating function):
	  if (fUseSum) {
	    cDtheta = dR0*cDenom/dJ01;
	    fHistProReDtheta->Fill(theta,cDtheta.Re());
	    fHistProImDtheta->Fill(theta,cDtheta.Im());
	  }
	 
	  cDenom *= TComplex::Power(i, m-1);
	  //cerr<<"TComplex::Power(i, m-1) = "<<TComplex::Power(i, m-1).Rho()<<endl; //checked ok
	 
	  //v as a function of eta for RP selection
	  for (Int_t be=1;be<=iNbinsEta;be++)  {
	    dEtaRP = fHist2RP[theta]->GetBinCenter(be);
	    cNumerRP = fHist2RP[theta]->GetNumerEta(be);
	    if (cNumerRP.Rho()==0) {
	      if (fDebug) cerr<<"WARNING: modulus of cNumerRP is zero in Finish()"<<endl;
	      dReRatioRP = 0;
	    }
	    else if (cDenom.Rho()==0) {
	      if (fDebug) cerr<<"WARNING: modulus of cDenom is zero"<<endl;
	      dReRatioRP = 0;
	    }
	    else {
	      dReRatioRP = (cNumerRP/cDenom).Re();
	    }
	    dVetaRP = dBesselRatio[m-1]*dReRatioRP*dVtheta; //BP eq. 12
	    fHistProVetaRP->Fill(dEtaRP,dVetaRP);
	  } //loop over bins be
	 

	  //v as a function of eta for POI selection
	  for (Int_t be=1;be<=iNbinsEta;be++)  {
	    dEtaPOI = fHist2POI[theta]->GetBinCenter(be);
	    cNumerPOI = fHist2POI[theta]->GetNumerEta(be);
	    if (cNumerPOI.Rho()==0) {
	      if (fDebug) cerr<<"WARNING: modulus of cNumerPOI is zero in Finish()"<<endl;
	      dReRatioPOI = 0;
	    }
	    else if (cDenom.Rho()==0) {
	      if (fDebug) cerr<<"WARNING: modulus of cDenom is zero"<<endl;
	      dReRatioPOI = 0;
	    }
	    else {
	      dReRatioPOI = (cNumerPOI/cDenom).Re();
	    }
	    dVetaPOI = dBesselRatio[m-1]*dReRatioPOI*dVtheta; //BP eq. 12
	    fHistProVetaPOI->Fill(dEtaPOI,dVetaPOI);
	  } //loop over bins be


	  //v as a function of Pt for RP selection
	  for (Int_t bp=1;bp<=iNbinsPt;bp++)  {
	    dPtRP = fHist2RP[theta]->GetBinCenterPt(bp);
	    cNumerRP = fHist2RP[theta]->GetNumerPt(bp);
	    if (cNumerRP.Rho()==0) {
	      if (fDebug) cerr<<"WARNING: modulus of cNumerRP is zero"<<endl;
	      dReRatioRP = 0;
	    }
	    else if (cDenom.Rho()==0) {
	      if (fDebug) cerr<<"WARNING: modulus of cDenom is zero"<<endl;
	      dReRatioRP = 0;
	    }
	    else {
	      dReRatioRP = (cNumerRP/cDenom).Re();
	    }
	    dVPtRP = dBesselRatio[m-1]*dReRatioRP*dVtheta; //BP eq. 12
	    fHistProVPtRP->Fill(dPtRP,dVPtRP);
	  } //loop over bins bp



	  //v as a function of Pt for POI selection
	  for (Int_t bp=1;bp<=iNbinsPt;bp++)  {
	    dPtPOI = fHist2POI[theta]->GetBinCenterPt(bp);
	    cNumerPOI = fHist2POI[theta]->GetNumerPt(bp);
	    if (cNumerPOI.Rho()==0) {
	      if (fDebug) cerr<<"WARNING: modulus of cNumerPOI is zero"<<endl;
	      dReRatioPOI = 0;
	    }
	    else if (cDenom.Rho()==0) {
	      if (fDebug) cerr<<"WARNING: modulus of cDenom is zero"<<endl;
	      dReRatioPOI = 0;
	    }
	    else {
	      dReRatioPOI = (cNumerPOI/cDenom).Re();
	    }
	    dVPtPOI = dBesselRatio[m-1]*dReRatioPOI*dVtheta; //BP eq. 12
	    fHistProVPtPOI->Fill(dPtPOI,dVPtPOI);
	  } //loop over bins bp

	}
      }

    }//end of loop over theta

    //calculate the average of fVtheta = fV
    dV /= iNtheta;

    //sigma2 and chi (for statistical error calculations)
    Double_t  dSigma2 = 0;
    Double_t  dChi= 0;
    if (fEventNumber!=0) {
      *fQsum /= fEventNumber;
      //cerr<<"fQsum.X() = "<<fQsum.X()<<endl;
      //cerr<<"fQsum.Y() = "<<fQsum.Y()<<endl;
      fQ2sum /= fEventNumber;
      //cout<<"fQ2sum = "<<fQ2sum<<endl;
      dSigma2 = fQ2sum - TMath::Power(fQsum->X(),2.) - TMath::Power(fQsum->Y(),2.) - TMath::Power(dV,2.);  //BP eq. 62
      if (dSigma2>0) dChi = dV/TMath::Sqrt(dSigma2);
      else dChi = -1.;
      fCommonHistsRes->FillChiRP(dChi);

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
      fCommonHistsRes->FillIntegratedFlowRP(dV,dVErr); 

      cout<<"*************************************"<<endl;
      cout<<"*************************************"<<endl;
      cout<<"      Integrated flow from           "<<endl;
      cout<<"        Lee-Yang Zeroes              "<<endl;
      cout<<endl;
      cout<<"dV(RP) = "<<dV<<" +- "<<dVErr<<" and chi = "<<dChi<<endl;
    }
	     
    //copy content of profile into TH1D and add error and fill the AliFlowCommonHistResults

    //v as a function of eta for RP selection
    for(Int_t b=0;b<iNbinsEta;b++) {
      Double_t dv2pro  = fHistProVetaRP->GetBinContent(b);
      Double_t dNprime = fCommonHists->GetEntriesInEtaBinRP(b);  
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
    } //loop over bins b
          
    //v as a function of Pt for POI selection 
    TH1F* fHistPtDiff = fCommonHists->GetHistPtDiff(); //for calculating integrated flow
    Double_t dVPOI = 0.;
    Double_t dSum = 0.;
    Double_t dErrV =0.;

    for(Int_t b=0;b<iNbinsPt;b++) {
      Double_t dv2pro = fHistProVPtPOI->GetBinContent(b);
      Double_t dNprime = fCommonHists->GetEntriesInPtBinPOI(b);   
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
      if (fHistPtDiff){
	Double_t dYieldPt = fHistPtDiff->GetBinContent(b);
	dVPOI += dv2pro*dYieldPt;
	dSum +=dYieldPt;
	dErrV += dYieldPt*dYieldPt*dErrdifcomb*dErrdifcomb;
      } else { cout<<"fHistPtDiff is NULL"<<endl; }
    } //loop over bins b

    if (dSum != 0.) {
      dVPOI /= dSum; //the pt distribution should be normalised
      dErrV /= (dSum*dSum);
      dErrV = TMath::Sqrt(dErrV);
    }
    cout<<"dV(POI) is "<<dVPOI<<" +- "<<dErrV<<endl;
    cout<<endl;
    cout<<"*************************************"<<endl;
    cout<<"*************************************"<<endl;
    fCommonHistsRes->FillIntegratedFlowPOI(dVPOI,dErrV);

  } //secondrun
   
  cout<<"----LYZ analysis finished....----"<<endl<<endl;

  return kTRUE;
}


//-----------------------------------------------------------------------
 
 Bool_t AliFlowAnalysisWithLeeYangZeros::FillFromFlowEvent(AliFlowEventSimple* anEvent) 
{ 
  // Get event quantities from AliFlowEvent for all particles

  if (fDebug) cout<<"****AliFlowAnalysisWithLeeYangZeros::FillFromFlowEvent()****"<<endl;
   
  if (!anEvent){
    cout<<"##### FlowLeeYangZero: FlowEvent pointer null"<<endl;
    return kFALSE;
  }
   
  //define variables
  TComplex cExpo, cGtheta, cGthetaNew, cZ;
  Int_t iNtheta = AliFlowLYZConstants::kTheta;
  Int_t iNbins = AliFlowLYZConstants::kNbins;
   

  //calculate flow
  Double_t dOrder = 2.;
      
  //get the Q vector 
  AliFlowVector vQ = anEvent->GetQ();
  //weight by the multiplicity
  Double_t dQX = 0;
  Double_t dQY = 0;
  if (vQ.GetMult() != 0) {
    dQX = vQ.X()/vQ.GetMult(); 
    dQY = vQ.Y()/vQ.GetMult();
  }
  vQ.Set(dQX,dQY);

  //for chi calculation:
  *fQsum += vQ;
  fHistQsumforChi->SetBinContent(1,fQsum->X());
  fHistQsumforChi->SetBinContent(2,fQsum->Y());
  fQ2sum += vQ.Mod2();
  fHistQsumforChi->SetBinContent(3,fQ2sum);
  //cerr<<"fQ2sum = "<<fQ2sum<<endl;

  for (Int_t theta=0;theta<iNtheta;theta++)
    {
      Double_t dTheta = ((double)theta/iNtheta)*TMath::Pi()/dOrder; 
	  
      //calculate dQtheta = cos(dOrder*(fPhi-dTheta);the projection of the Q vector on the reference direction dTheta
      Double_t dQtheta = GetQtheta(vQ, dTheta);
	     	   
      for (Int_t bin=1;bin<=iNbins;bin++)
	{
	  Double_t dR = fHist1[theta]->GetBinCenter(bin); //bincentre of bins in histogram  //FIXED???
	  //if (theta == 0) cerr<<"cerr::dR = "<<dR<<endl;
	  if (fUseSum)
	    {
	      //calculate the sum generating function
	      cExpo(0.,dR*dQtheta);                       //Re=0 ; Im=dR*dQtheta
	      cGtheta = TComplex::Exp(cExpo);
	    }
	  else
	    {
	      //calculate the product generating function
	      cGtheta = GetGrtheta(anEvent, dR, dTheta);  
	      if (cGtheta.Rho2() > 100.) break;
	    }

	  fHist1[theta]->Fill(dR,cGtheta);              //fill real and imaginary part of cGtheta

	} //loop over bins
    } //loop over theta 
 
   
  return kTRUE;

          
}

 //-----------------------------------------------------------------------   
 Bool_t AliFlowAnalysisWithLeeYangZeros::SecondFillFromFlowEvent(AliFlowEventSimple* anEvent) 
{ 
  //for differential flow

  if (fDebug) cout<<"****AliFlowAnalysisWithLeeYangZeros::SecondFillFromFlowEvent()****"<<endl;
    
  if (!anEvent){
    cout<<"##### FlowLeeYangZero: FlowEvent pointer null"<<endl;
    return kFALSE;
  }
   
  //define variables
  TComplex cExpo, cDenom, cNumerRP, cNumerPOI, cCosTermComplex;
  Double_t dPhi, dEta, dPt;
  Double_t dR0 = 0.;
  Double_t dCosTermRP = 0.;
  Double_t dCosTermPOI = 0.;
  Double_t m = 1.;
  Double_t dOrder = 2.;
  Int_t iNtheta = AliFlowLYZConstants::kTheta;
  
   
  //get the Q vector 
  AliFlowVector vQ = anEvent->GetQ();
  //weight by the multiplicity
  Double_t dQX = 0.;
  Double_t dQY = 0.;
  if (vQ.GetMult() != 0) {
    dQX = vQ.X()/vQ.GetMult(); 
    dQY = vQ.Y()/vQ.GetMult();
  }
  vQ.Set(dQX,dQY); 
              
  //for chi calculation:
  *fQsum += vQ;
  fHistQsumforChi->SetBinContent(1,fQsum->X());
  fHistQsumforChi->SetBinContent(2,fQsum->Y());
  fQ2sum += vQ.Mod2();
  fHistQsumforChi->SetBinContent(3,fQ2sum);

  for (Int_t theta=0;theta<iNtheta;theta++)
    {
      Double_t dTheta = ((double)theta/iNtheta)*TMath::Pi()/dOrder;   

      //calculate dQtheta = cos(dOrder*(dPhi-dTheta);the projection of the Q vector on the reference direction dTheta	  
      Double_t dQtheta = GetQtheta(vQ, dTheta);
      //cerr<<"dQtheta for fdenom = "<<dQtheta<<endl;
  	 
      //denominator for differential v
      if (fHistProR0theta) {
	dR0 = fHistProR0theta->GetBinContent(theta+1);
      }
      else {
	cout <<"pointer fHistProR0theta does not exist" << endl;
      }
      //cerr<<"dR0 = "<<dR0 <<endl;

      if (fUseSum) //sum generating function
	{
	  cExpo(0.,dR0*dQtheta);
	  cDenom = dQtheta*(TComplex::Exp(cExpo)); //BP eq 12
	  //loop over tracks in event
	  Int_t iNumberOfTracks = anEvent->NumberOfTracks();
	  for (Int_t i=0;i<iNumberOfTracks;i++)  {
	    AliFlowTrackSimple*  pTrack = anEvent->GetTrack(i);
	    if (pTrack) {
	      dEta = pTrack->Eta();
	      dPt = pTrack->Pt();
	      dPhi = pTrack->Phi();
	      if (pTrack->UseForIntegratedFlow()) { // RP selection
		dCosTermRP = cos(m*dOrder*(dPhi-dTheta));
		cNumerRP = dCosTermRP*(TComplex::Exp(cExpo));
		if (cNumerRP.Rho()==0) {cerr<<"WARNING: modulus of cNumerRP is zero in SecondFillFromFlowEvent"<<endl;}
		if (fDebug) cerr<<"modulus of cNumerRP is "<<cNumerRP.Rho()<<endl;
		if (fHist2RP[theta]) {
		  fHist2RP[theta]->Fill(dEta,dPt,cNumerRP); }
	      }
	      if (pTrack->UseForDifferentialFlow()) { //POI selection
		dCosTermPOI = cos(m*dOrder*(dPhi-dTheta));
		cNumerPOI = dCosTermPOI*(TComplex::Exp(cExpo));
		if (cNumerPOI.Rho()==0) {cerr<<"WARNING: modulus of cNumerPOI is zero in SecondFillFromFlowEvent"<<endl;}
		if (fDebug) cerr<<"modulus of cNumerPOI is "<<cNumerPOI.Rho()<<endl;
		if (fHist2POI[theta]) {
		  fHist2POI[theta]->Fill(dEta,dPt,cNumerPOI); }
	      }
	      else {
		cout << "fHist2 pointer mising" <<endl;
	      }
	    } //if track
	    else {cerr << "no particle!!!"<<endl;}
	  } //loop over tracks
		  
	} //sum
      else        //product generating function
	{
	  cDenom = GetDiffFlow(anEvent, dR0, theta); 
		   
	}//product
      if (fHistProReDenom && fHistProImDenom) {
	fHistProReDenom->Fill(theta,cDenom.Re());               //fill the real part of fDenom
	fHistProImDenom->Fill(theta,cDenom.Im());               //fill the imaginary part of fDenom
      }
      else {
	cout << "Pointers to cDenom  mising" << endl;
      }
     	       
    }//end of loop over theta
  
  return kTRUE;
 
  
}
 //-----------------------------------------------------------------------   
 Double_t AliFlowAnalysisWithLeeYangZeros::GetQtheta(AliFlowVector aQ, Double_t aTheta) 
{
  //calculate Qtheta. Qtheta is the sum over all particles of cos(dOrder*(dPhi-dTheta)) BP eq. 3
  if (fDebug) cout<<"****AliFlowAnalysisWithLeeYangZeros::GetQtheta()****"<<endl;

  Double_t dQtheta = 0.;
  Double_t dOrder = 2.;
  
  dQtheta = aQ.X()*cos(dOrder*aTheta)+aQ.Y()*sin(dOrder*aTheta);

  return dQtheta;
 
}
 

//-----------------------------------------------------------------------   
TComplex AliFlowAnalysisWithLeeYangZeros::GetGrtheta(AliFlowEventSimple* anEvent, Double_t aR, Double_t aTheta) 
{
  // Product Generating Function for LeeYangZeros method
  // PG Eq. 3 (J. Phys. G Nucl. Part. Phys 30 S1213 (2004))
  
  if (fDebug) cout<<"****AliFlowAnalysisWithLeeYangZeros::GetGrtheta()****"<<endl;
  
  
  TComplex cG = TComplex::One();
  Double_t dOrder =  2.;
  Double_t dWgt = 1.;
  
  Int_t iNumberOfTracks = anEvent->NumberOfTracks();
  
  for (Int_t i=0;i<iNumberOfTracks;i++) //loop over tracks in event
    {
      AliFlowTrackSimple* pTrack = anEvent->GetTrack(i) ; 
      if (pTrack){
	if (pTrack->UseForIntegratedFlow()) {
	  Double_t dPhi = pTrack->Phi();
	  Double_t dGIm = aR * dWgt*cos(dOrder*(dPhi - aTheta));
	  TComplex cGi(1., dGIm);
	  cG *= cGi;     //product over all tracks
	}
      }
      else {cerr << "no particle pointer !!!"<<endl;}
    }//loop over tracks
  
  return cG;
  
} 


//-----------------------------------------------------------------------   
TComplex AliFlowAnalysisWithLeeYangZeros::GetDiffFlow(AliFlowEventSimple* anEvent, Double_t aR0, Int_t theta) 
{
  // Sum for the denominator for diff. flow for the Product Generating Function for LeeYangZeros method
  // PG Eq. 9 (J. Phys. G Nucl. Part. Phys 30 S1213 (2004))
  // Also for v1 mixed harmonics: DF Eq. 5
  // It is the deriverative of Grtheta at r0 divided by Grtheta at r0
  
  if (fDebug) cout<<"****AliFlowAnalysisWithLeeYangZeros::GetGrtheta()****"<<endl;
  
  TComplex cG = TComplex::One();
  TComplex cdGr0(0.,0.);
  Double_t dOrder =  2.;
  Double_t dWgt = 1.;
  
  Int_t iNumberOfTracks = anEvent->NumberOfTracks();
  
  Int_t iNtheta = AliFlowLYZConstants::kTheta;
  Double_t dTheta = ((double)theta/iNtheta)*TMath::Pi()/dOrder;
  
  //for the denominator (use all RP selected particles)
  for (Int_t i=0;i<iNumberOfTracks;i++) //loop over tracks in event
    {
      AliFlowTrackSimple* pTrack = anEvent->GetTrack(i) ;  
      if (pTrack){
	if (pTrack->UseForIntegratedFlow()) {
	  Double_t dPhi = pTrack->Phi();
	  Double_t dCosTerm = dWgt*cos(dOrder*(dPhi - dTheta));
	  //GetGr0theta
	  Double_t dGIm = aR0 * dCosTerm;
	  TComplex cGi(1., dGIm);
	  TComplex cCosTermComplex(1., aR0*dCosTerm);
	  cG *= cGi;     //product over all tracks
	  //GetdGr0theta
	  cdGr0 +=(dCosTerm / cCosTermComplex);  //sum over all tracks
	}
      } //if particle
      else {cerr << "no particle!!!"<<endl;}
    }//loop over tracks
  
  //for the numerator
  for (Int_t i=0;i<iNumberOfTracks;i++) 
    {
      AliFlowTrackSimple* pTrack = anEvent->GetTrack(i) ;  
      if (pTrack){
	Double_t dEta = pTrack->Eta();
	Double_t dPt = pTrack->Pt();
	Double_t dPhi = pTrack->Phi();
	Double_t dCosTerm = cos(dOrder*(dPhi-dTheta));
	TComplex cCosTermComplex(1.,aR0*dCosTerm);
	//RP selection
	if (pTrack->UseForIntegratedFlow()) {
	  TComplex cNumerRP = cG*dCosTerm/cCosTermComplex;  //PG Eq. 9
	  fHist2RP[theta]->Fill(dEta,dPt,cNumerRP);  
	}
	//POI selection
	if (pTrack->UseForDifferentialFlow()) {
	  TComplex cNumerPOI = cG*dCosTerm/cCosTermComplex;  //PG Eq. 9
	  fHist2POI[theta]->Fill(dEta,dPt,cNumerPOI);  
	}
      } //if particle
      else {cerr << "no particle pointer!!!"<<endl;}
    }//loop over tracks
  
  TComplex cDenom = cG*cdGr0;  
  return cDenom;
  
} 

//----------------------------------------------------------------------- 

