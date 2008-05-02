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

/*
$Log$
*/ 

#include "Riostream.h"                 //needed as include
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
    fQ2sum(0),
    fQtheta(0),
    fEventNumber(0),
    fMult(0),
    fNbins(0),
    fTheta(0),
    fEvent(0),
    fTrack(0),
    fFirstRun(kTRUE),
    fUseSum(kTRUE),
    fDoubleLoop(kFALSE),
    fDebug(kFALSE),
    fHistFileName(0),
    fHistFile(0),
    fSummaryFile(0),
    firstRunFileName(0),
    firstRunFile(0),
    fHistProVtheta(0),
    fHistProVeta(0),
    fHistProVPt(0),
    fHistProR0theta(0),
    fHistProReDenom(0),
    fHistProImDenom(0),
    fHistProReDtheta(0),
    fHistProImDtheta(0),
    fCommonHists(0),
    fCommonHistsRes(0)
  
{
  //default constructor
  if (fDebug) cout<<"****AliFlowAnalysisWithLeeYangZeros::AliFlowAnalysisWithLeeYangZeros default constructor****"<<endl;

  // output file (histograms)
  TString fHistFileName = "outputFromLYZAnalysis.root" ;

  for(Int_t i = 0;i<5;i++)
    {
      fHist1[i]=0;
      fHist2[i]=0;
    }

  fQ.Set(0.,0.);
  fQsum.Set(0.,0.);

}


//-----------------------------------------------------------------------


 AliFlowAnalysisWithLeeYangZeros::~AliFlowAnalysisWithLeeYangZeros() 
 {
   //default destructor
   if (fDebug) cout<<"****~AliFlowAnalysisWithLeeYangZeros****"<<endl;
   delete fHistFile;
 }
 
 //-----------------------------------------------------------------------

Bool_t AliFlowAnalysisWithLeeYangZeros::Init() 
{
  //init method 
  if (fDebug) cout<<"****AliFlowAnalysisWithLeeYangZeros::Init()****"<<endl;

  // Open output files (->plots)
  fHistFile = new TFile(fHistFileName.Data(), "RECREATE");
  
    
  // Book histograms
  Int_t fNtheta = AliFlowLYZConstants::kTheta;
  Int_t fNbinsPt = AliFlowCommonConstants::GetNbinsPt();
  Int_t fNbinsEta = AliFlowCommonConstants::GetNbinsEta();

  Double_t  fPtMin = AliFlowCommonConstants::GetPtMin();	     
  Double_t  fPtMax = AliFlowCommonConstants::GetPtMax();
  Double_t  fEtaMin = AliFlowCommonConstants::GetEtaMin();	     
  Double_t  fEtaMax = AliFlowCommonConstants::GetEtaMax();
  
  
  //for control histograms
  fCommonHists = new AliFlowCommonHist("LYZ");
  fCommonHistsRes = new AliFlowCommonHistResults("LYZ");

  //for first loop over events 
  if (fFirstRun){
    fHistProR0theta  = new TProfile("First_FlowPro_r0theta_LYZ","First_FlowPro_r0theta_LYZ",fNtheta,-0.5,fNtheta-0.5);
    fHistProR0theta->SetXTitle("#theta");
    fHistProR0theta->SetYTitle("r_{0}^{#theta}");

    fHistProVtheta  = new TProfile("First_FlowPro_Vtheta_LYZ","First_FlowPro_Vtheta_LYZ",fNtheta,-0.5,fNtheta-0.5);
    fHistProVtheta->SetXTitle("#theta");
    fHistProVtheta->SetYTitle("V_{n}^{#theta}");        

    //class AliFlowLYZHist1 defines the histograms: fHistProGtheta, fHistProReGtheta, fHistProImGtheta, fHistProR0theta
    for (Int_t theta=0;theta<fNtheta;theta++) {  
      fHist1[theta]=new AliFlowLYZHist1(theta);
    }
     
  }
  //for second loop over events 
  else {
    fHistProReDenom = new TProfile("Second_FlowPro_ReDenom_LYZ","Second_FlowPro_ReDenom_LYZ" , fNtheta, -0.5, fNtheta-0.5);
    fHistProReDenom->SetXTitle("#theta");
    fHistProReDenom->SetYTitle("Re(Q^{#theta}e^{ir_{0}^{#theta}Q^{#theta}})");

    fHistProImDenom = new TProfile("Second_FlowPro_ImDenom_LYZ","Second_FlowPro_ImDenom_LYZ" , fNtheta, -0.5, fNtheta-0.5);
    fHistProImDenom->SetXTitle("#theta");
    fHistProImDenom->SetYTitle("Im(Q^{#theta}e^{ir_{0}^{#theta}Q^{#theta}})");

    fHistProVeta = new TProfile("Second_FlowPro_Veta_LYZ","Second_FlowPro_Veta_LYZ",fNbinsEta,fEtaMin,fEtaMax);
    fHistProVeta->SetXTitle("rapidity");
    fHistProVeta->SetYTitle("v (%)");

    fHistProVPt = new TProfile("Second_FlowPro_VPt_LYZ","Second_FlowPro_VPt_LYZ",fNbinsPt,fPtMin,fPtMax);
    fHistProVPt->SetXTitle("Pt");
    fHistProVPt->SetYTitle("v (%)");

    fHistProReDtheta = new TProfile("Second_FlowPro_ReDtheta_LYZ","Second_FlowPro_ReDtheta_LYZ",fNtheta, -0.5, fNtheta-0.5);
    fHistProReDtheta->SetXTitle("#theta");
    fHistProReDtheta->SetYTitle("Re(D^{#theta})");

    fHistProImDtheta = new TProfile("Second_FlowPro_ImDtheta_LYZ","Second_FlowPro_ImDtheta_LYZ",fNtheta, -0.5, fNtheta-0.5);
    fHistProImDtheta->SetXTitle("#theta");
    fHistProImDtheta->SetYTitle("Im(D^{#theta})");

    //class AliFlowLYZHist2 defines the histograms: 
    for (Int_t theta=0;theta<fNtheta;theta++)  {  
      fHist2[theta]=new AliFlowLYZHist2(theta);
    }
      
    //read hists from first run file
    //firstRunFile = new TFile("fof_flowLYZAnal_firstrun.root","READ");  //default is read
    if (firstRunFile->IsZombie()){ //check if file exists
      cout << "Error opening file, run first with fFirstrun = kTRUE" << endl;
      exit(-1);
    } else if (firstRunFile->IsOpen()){
      cout<<"----firstRunFile is open----"<<endl<<endl;
      fHistProR0theta  = (TProfile*)firstRunFile->Get("First_FlowPro_r0theta_LYZ");
    }    
  }
   

  if (fDebug) cout<<"****Histograms initialised****"<<endl;
    
  fEventNumber = 0; //set event counter to zero
 
  return kTRUE; 
}
 
 //-----------------------------------------------------------------------
 
Bool_t AliFlowAnalysisWithLeeYangZeros::Make(AliFlowEventSimple* fEvent) 
{
  //make method
  if (fDebug) cout<<"****AliFlowAnalysisWithLeeYangZeros::Make()****"<<endl;
        
  //get tracks from event
  if (fEvent) {
    if (fFirstRun){
      fCommonHists->FillControlHistograms(fEvent);
      FillFromFlowEvent(fEvent);
    }
    else {
      fCommonHists->FillControlHistograms(fEvent);
      SecondFillFromFlowEvent(fEvent);
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
  Double_t  fJ01 = 2.405; 
  Int_t fNtheta = AliFlowLYZConstants::kTheta;
  
  if (fFirstRun){
    Double_t  fR0 = 0;
    Double_t  fVtheta = 0; 
    Double_t  fv = 0;
    Double_t  fV = 0; 
    for (Int_t theta=0;theta<fNtheta;theta++)
      {
	//get the first minimum r0
	fHist1[theta]->FillGtheta();
	fR0 = fHist1[theta]->GetR0();
	    	   	   
	//calculate integrated flow
	if (fR0!=0) fVtheta = fJ01/fR0;
	//for estimating systematic error resulting from r0
	Double_t binsize = (AliFlowLYZConstants::fgMax)/(AliFlowLYZConstants::kNbins);
	Double_t fVplus = fJ01/(fR0+binsize);
	Double_t fVmin = fJ01/(fR0-binsize);
	fv = fVtheta;
	Double_t fvplus = fVplus;
	Double_t fvmin = fVmin;
	cout<<"fv = "<<fv<<" and fvplus = "<<fvplus<< " and fvmin = "<<fvmin<<endl;
	     
	//fill the histograms
	fHistProR0theta->Fill(theta,fR0);
	fHistProVtheta->Fill(theta,fVtheta); 
	//get average value of fVtheta = fV
	fV += fVtheta;
      } //end of loop over theta

    //get average value of fVtheta = fV
    fV /=fNtheta;
       
    //sigma2 and chi 
    Double_t  fSigma2 = 0;
    Double_t  fChi= 0;
    if (fEventNumber!=0) {
      fQsum /= fEventNumber;
      fQ2sum /= fEventNumber;
      fSigma2 = fQ2sum - TMath::Power(fQsum.X(),2.) - TMath::Power(fQsum.Y(),2.) - TMath::Power(fV,2.);  //BP eq. 62
      if (fSigma2>0) fChi = fV/TMath::Sqrt(fSigma2);
      else fChi = -1.;
      fCommonHistsRes->FillChi(fChi);
      cout<<"fV = "<<fV<<" and chi = "<<fChi<<endl;
    }
	   
    // recalculate statistical errors on integrated flow
    //combining 5 theta angles to 1 relative error BP eq. 89
    Double_t fRelErr2comb = 0.;
    Int_t nEvts = fEventNumber;
    for (Int_t theta=0;theta<fNtheta;theta++){
      Double_t fTheta = ((double)theta/fNtheta)*TMath::Pi(); 
      Double_t fApluscomb = TMath::Exp((fJ01*fJ01)/(2*fChi*fChi)*
				       TMath::Cos(fTheta));
      Double_t fAmincomb = TMath::Exp(-(fJ01*fJ01)/(2*fChi*fChi)*
				      TMath::Cos(fTheta));
      fRelErr2comb += (1/(2*nEvts*(fJ01*fJ01)*TMath::BesselJ1(fJ01)*
			  TMath::BesselJ1(fJ01)))*
	(fApluscomb*TMath::BesselJ0(2*fJ01*TMath::Sin(fTheta/2)) + 
	 fAmincomb*TMath::BesselJ0(2*fJ01*TMath::Cos(fTheta/2)));
    }
    fRelErr2comb /= fNtheta;
    Double_t fRelErrcomb = TMath::Sqrt(fRelErr2comb);
    cout<<"fRelErrcomb = "<<fRelErrcomb<<endl;         

    //copy content of profile into TH1D and add error
    Double_t fv2pro = fV;   //in the case that fv is equal to fV
    Double_t fv2Err = fv2pro*fRelErrcomb ; 
    cout<<"fv2pro +- fv2Err = "<<fv2pro<<" +- "<<fv2Err<<endl;
    fCommonHistsRes->FillIntegratedFlow(fv2pro, fv2Err); 


    if (fDebug) cout<<"****histograms filled****"<<endl;  
    
    //save histograms in file //temp for testing selector
    fHistFile->cd();
    fHistFile->Write();

    return kTRUE;
    fEventNumber =0; //set to zero for second round over events
  }  //firstrun
   
 
  else {  //second run
   
    //calculate differential flow
    //declare variables
    TComplex fDenom, fNumer, fDtheta;
    Int_t m = 1;
    TComplex i = TComplex::I();
    Double_t fBesselRatio[3] = {1., 1.202, 2.69};
    Int_t fNbinsPt = AliFlowCommonConstants::GetNbinsPt();
    Int_t fNbinsEta = AliFlowCommonConstants::GetNbinsEta();

    Double_t fEta, fPt, fReRatio, fVeta, fVPt;
        
     
    Double_t fR0 = 0.; 
    Double_t fVtheta = 0.;
    Double_t fV = 0.;
    Double_t fReDenom = 0.;
    Double_t fImDenom = 0.; 
    for (Int_t theta=0;theta<fNtheta;theta++)  { //loop over theta
      if (fHistProR0theta) {
	fR0 = fHistProR0theta->GetBinContent(theta+1);
	if (fDebug) cerr<<"fR0 = "<<fR0<<endl;
	if (fR0!=0) fVtheta = fJ01/fR0;
	fV += fVtheta; 
	// BP Eq. 9  -> Vn^theta = j01/r0^theta
	if (fHistProReDenom && fHistProImDenom) {
	  fReDenom = fHistProReDenom->GetBinContent(theta+1);
	  fImDenom = fHistProImDenom->GetBinContent(theta+1);
	}
	else {
	  cout << "Hist pointer fDenom in file does not exist" <<endl;
	}
	  
      }
      else {
	cout << "Hist pointer R0theta in file does not exist" <<endl;
      }
      //} //loop over theta
      
      fDenom(fReDenom,fImDenom);
	 

      //for new method and use by others (only with the sum generating function):
      if (fUseSum) {
	fR0 = fHistProR0theta->GetBinContent(theta+1); 
	fDtheta = fR0*fDenom/fJ01;
	fHistProReDtheta->Fill(theta,fDtheta.Re());
	fHistProImDtheta->Fill(theta,fDtheta.Im());
      }
	 
      fDenom *= TComplex::Power(i, m-1);
      //cerr<<"TComplex::Power(i, m-1) = "<<TComplex::Power(i, m-1).Rho()<<endl; //checked ok
	 
      //v as a function of eta
      for (Int_t be=1;be<=fNbinsEta;be++)  {
	fEta = fHist2[theta]->GetBinCenter(be);
	fNumer = fHist2[theta]->GetfNumer(be);
	if (fNumer.Rho()==0) {
	  if (fDebug) cerr<<"WARNING: modulus of fNumer is zero in Finish()"<<endl;
	  fReRatio = 0;
	}
	else if (fDenom.Rho()==0) {
	  if (fDebug) cerr<<"WARNING: modulus of fDenom is zero"<<endl;
	  fReRatio = 0;
	}
	else {
	  //if ( j==1 && theta==0) cerr<<"modulus of fNumer = "<<fNumer.Rho() <<endl; //always a number smaller than 1, or 0.
	  fReRatio = (fNumer/fDenom).Re();
	}

	fVeta = fBesselRatio[m-1]*fReRatio*fVtheta*100.; //BP eq. 12
	//if ( j==1 && theta==0) cerr<<"eta = "<<fEta<<" cerr::fReRatio for eta = "<<fReRatio<<" cerr::fVeta for eta = "<<fVeta<<endl;
	   
	fHistProVeta->Fill(fEta,fVeta);
      } //loop over bins be
	 
      //v as a function of Pt
      Int_t fNbinsPt = AliFlowCommonConstants::GetNbinsPt();
      for (Int_t bp=1;bp<=fNbinsPt;bp++)  {
	fPt = fHist2[theta]->GetBinCenterPt(bp);
	fNumer = fHist2[theta]->GetfNumerPt(bp);
	if (fNumer.Rho()==0) {
	  if (fDebug) cerr<<"modulus of fNumer is zero"<<endl;
	  fReRatio = 0;
	}
	else if (fDenom.Rho()==0) {
	  if (fDebug) cerr<<"modulus of fDenom is zero"<<endl;
	  fReRatio = 0;
	}
	else {
	  //if ( j==1 && theta==0) cerr<<"modulus of fNumer = "<<fNumer.Rho() <<endl; //always a number smaller than 1, or 0.
	  fReRatio = (fNumer/fDenom).Re();
	}
   
	fVPt = fBesselRatio[m-1]*fReRatio*fVtheta*100.; //BP eq. 12
	      
	fHistProVPt->Fill(fPt,fVPt);
      } //loop over bins bp

    }//end of loop over theta

    //calculate the average of fVtheta = fV
    fV /= fNtheta;

    //sigma2 and chi (for statistical error calculations)
    Double_t  fSigma2 = 0;
    Double_t  fChi= 0;
    if (fEventNumber!=0) {
      fQsum /= fEventNumber;
      //cerr<<"fQsum.X() = "<<fQsum.X()<<endl;
      //cerr<<"fQsum.Y() = "<<fQsum.Y()<<endl;
      fQ2sum /= fEventNumber;
      //cout<<"fQ2sum = "<<fQ2sum<<endl;
      fSigma2 = fQ2sum - TMath::Power(fQsum.X(),2.) - TMath::Power(fQsum.Y(),2.) - TMath::Power(fV,2.);  //BP eq. 62
      if (fSigma2>0) fChi = fV/TMath::Sqrt(fSigma2);
      else fChi = -1.;
      fCommonHistsRes->FillChi(fChi);
      cout<<"fV = "<<fV<<" and chi = "<<fChi<<endl;
    }
	     
    //copy content of profile into TH1D and add error
    for(Int_t b=0;b<fNbinsPt;b++) {
      Double_t fv2pro = fHistProVPt->GetBinContent(b);
      Double_t fNprime = fCommonHists->GetEntriesInPtBin(b);
      Double_t fErrdifcomb = 0.;  //set error to zero
      Double_t fErr2difcomb = 0.; //set error to zero
      //calculate error
      if (fNprime!=0.) { 
	for (Int_t theta=0;theta<fNtheta;theta++) {
	  Double_t fTheta = ((double)theta/fNtheta)*TMath::Pi(); 
	  Double_t fApluscomb = TMath::Exp((fJ01*fJ01)/(2*fChi*fChi)*
					   TMath::Cos(fTheta));
	  Double_t fAmincomb = TMath::Exp(-(fJ01*fJ01)/(2*fChi*fChi)*
					  TMath::Cos(fTheta));
	  fErr2difcomb += (TMath::Cos(fTheta)/(4*fNprime*TMath::BesselJ1(fJ01)*
						 TMath::BesselJ1(fJ01)))*
	    ((fApluscomb*TMath::BesselJ0(2*fJ01*TMath::Sin(fTheta/2))) - 
	     (fAmincomb*TMath::BesselJ0(2*fJ01*TMath::Cos(fTheta/2))));
	} //loop over theta
      } 
      
      if (fErr2difcomb!=0.) {
	fErr2difcomb /= fNtheta;
	fErrdifcomb = TMath::Sqrt(fErr2difcomb)*100;
	//cerr<<"fErrdifcomb = "<<fErrdifcomb<<endl;
      }
      else {fErrdifcomb = 0.;}
	  
      //fill TH1D
      fCommonHistsRes->FillDifferentialFlow(b, fv2pro, fErrdifcomb); 
    } //loop over bins b
 
    
    //save histograms in file
    fHistFile->cd();
    fHistFile->Write();
    fHistFile->Close();
    //Note that when the file is closed, all histograms and Trees in memory associated with this file are deleted
    if (fDebug) cout<<"****Histograms saved and fHistFile closed, all histograms deleted****"<<endl;
     
    //close the first run file 
    firstRunFile->Close();

     
  } //secondrun
   
  cout<<"----LYZ analysis finished....----"<<endl<<endl;

  return kTRUE;
}


//-----------------------------------------------------------------------
 
 Bool_t AliFlowAnalysisWithLeeYangZeros::FillFromFlowEvent(AliFlowEventSimple* fEvent) 
{ 
  // Get event quantities from AliFlowEvent for all particles

  if (fDebug) cout<<"****AliFlowAnalysisWithLeeYangZeros::FillFromFlowEvent()****"<<endl;
   
  if (!fEvent){
    cout<<"##### FlowLeeYangZero: FlowEvent pointer null"<<endl;
    return kFALSE;
  }
   
  //define variables
  TComplex fExpo, fGtheta, fGthetaNew, fZ;
  Int_t fNtheta = AliFlowLYZConstants::kTheta;
  Int_t fNbins = AliFlowLYZConstants::kNbins;
   

  //calculate flow
  Double_t fOrder = 2.;
      
  //get the Q vector 
  fQ = fEvent->GetQ();
  //weight by the multiplicity
  Double_t fQX = fQ.X()/fQ.GetMult(); 
  Double_t fQY = fQ.Y()/fQ.GetMult();
  fQ.Set(fQX,fQY);
  //for chi calculation:
  fQsum += fQ;
  fQ2sum += fQ.Mod2();
  //cerr<<"fQ2sum = "<<fQ2sum<<endl;

  for (Int_t theta=0;theta<fNtheta;theta++)
    {
      fTheta = ((double)theta/fNtheta)*TMath::Pi()/fOrder; 
	  
      //calculate fQtheta = cos(fOrder*(fPhi-fTheta);the projection of the Q vector on the reference direction fTheta
      fQtheta = GetQtheta(fQ, fTheta);
	     	   
      for (Int_t bin=1;bin<=fNbins;bin++)
	{
	  Double_t fR = fHist1[theta]->GetBinCenter(bin); //bincentre of bins in histogram  //FIXED???
	  //if (theta == 0) cerr<<"cerr::fR = "<<fR<<endl;
	  if (fUseSum)
	    {
	      //calculate the sum generating function
	      fExpo(0.,fR*fQtheta);                           //Re=0 ; Im=fR*fQtheta
	      fGtheta = TComplex::Exp(fExpo);
	    }
	  else
	    {
	      //calculate the product generating function
	      fGtheta = GetGrtheta(fEvent, fR, fTheta);  //make this function
	      if (fGtheta.Rho2() > 100.) break;
	    }

	  fHist1[theta]->Fill(fR,fGtheta);              //fill real and imaginary part of fGtheta

	} //loop over bins
    } //loop over theta 
 
   
  return kTRUE;

          
}

 //-----------------------------------------------------------------------   
 Bool_t AliFlowAnalysisWithLeeYangZeros::SecondFillFromFlowEvent(AliFlowEventSimple* fEvent) 
{ 
  //for differential flow

  if (fDebug) cout<<"****AliFlowAnalysisWithLeeYangZeros::SecondFillFromFlowEvent()****"<<endl;
    
  if (!fEvent){
    cout<<"##### FlowLeeYangZero: FlowEvent pointer null"<<endl;
    return kFALSE;
  }
   
  //define variables
  //TVector2 fQ;
  TComplex fExpo, fDenom, fNumer,fCosTermComplex;
  Double_t  fOrder, fPhi, fEta, fPt;
  Double_t fR0 = 0.;
  Double_t fCosTerm;
  Double_t m = 1.;
  Int_t fNtheta = AliFlowLYZConstants::kTheta;
     
  //calculate flow
  fOrder = 2.;
  
  //get the Q vector 
  fQ = fEvent->GetQ();
  //weight by the multiplicity
  Double_t fQX = fQ.X()/fQ.GetMult(); 
  Double_t fQY = fQ.Y()/fQ.GetMult();
  fQ.Set(fQX,fQY);               
  //for chi calculation:
  fQsum += fQ;
  fQ2sum += fQ.Mod2();

  for (Int_t theta=0;theta<fNtheta;theta++)
    {
      fTheta = ((double)theta/fNtheta)*TMath::Pi()/fOrder;   

      //calculate fQtheta = cos(fOrder*(fPhi-fTheta);the projection of the Q vector on the reference direction fTheta	  
      fQtheta = GetQtheta(fQ, fTheta);
      //cerr<<"fQtheta for fdenom = "<<fQtheta<<endl;
  	 
      //denominator for differential v
      if (fHistProR0theta) {
	fR0 = fHistProR0theta->GetBinContent(theta+1);
      }
      else {
	cout <<"pointer fHistProR0Theta does not exist" << endl;
      }
      //cerr<<"fR0 = "<<fR0 <<endl;

      if (fUseSum)                                                    //sum generating function
	{
	  fExpo(0.,fR0*fQtheta);
	  fDenom = fQtheta*(TComplex::Exp(fExpo)); //BP eq 12
	  //loop over tracks in event
	  Int_t fNumberOfTracks = fEvent->NumberOfTracks();
	  for (Int_t i=0;i<fNumberOfTracks;i++)  {
	    fTrack = fEvent->GetTrack(i);
	    if (fTrack) {
	      if (fTrack->UseForDifferentialFlow()) {
		fEta = fTrack->Eta();
		fPt = fTrack->Pt();
		fPhi = fTrack->Phi();
		fCosTerm = cos(m*fOrder*(fPhi-fTheta));
		//cerr<<"fCosTerm = "<<fCosTerm <<endl;
		fNumer = fCosTerm*(TComplex::Exp(fExpo));
		if (fNumer.Rho()==0) {cerr<<"WARNING: modulus of fNumer is zero in SecondFillFromFlowEvent"<<endl;}
		if (fDebug) cerr<<"modulus of fNumer is "<<fNumer.Rho()<<endl;
		if (fHist2[theta]) {
		  fHist2[theta]->Fill(fEta,fPt,fNumer);
		}
		else {
		  cout << "fHist2 pointer mising" <<endl;
		}
	      }
	    } //if particle
	    else {cerr << "no particle!!!"<<endl;}
	  } //loop over tracks
		  
	} //sum
      else                                                        //product generating function
	{
	  fDenom = GetDiffFlow(fEvent, fR0, theta); 
		   
	}//product
      if (fHistProReDenom && fHistProImDenom) {
	fHistProReDenom->Fill(theta,fDenom.Re());               //fill the real part of fDenom
	fHistProImDenom->Fill(theta,fDenom.Im());               //fill the imaginary part of fDenom
      }
      else {
	cout << "Pointers to fDenom  mising" << endl;
      }
     	       
    }//end of loop over theta
  
  return kTRUE;
 
  
}
 //-----------------------------------------------------------------------   
 Double_t AliFlowAnalysisWithLeeYangZeros::GetQtheta(TVector2 fQ, Double_t fTheta) 
{
  //calculate Qtheta. Qtheta is the sum over all particles of cos(fOrder*(fPhi-fTheta)) BP eq. 3
  if (fDebug) cout<<"****AliFlowAnalysisWithLeeYangZeros::GetQtheta()****"<<endl;

  Double_t fQtheta = 0.;
  Double_t fOrder = 2.;
  
  fQtheta = fQ.X()*cos(fOrder*fTheta)+fQ.Y()*sin(fOrder*fTheta);

  return fQtheta;
 
}
 

//-----------------------------------------------------------------------   
TComplex AliFlowAnalysisWithLeeYangZeros::GetGrtheta(AliFlowEventSimple* fEvent, Double_t fR, Double_t fTheta) 
{
  // Product Generating Function for LeeYangZeros method
  // PG Eq. 3 (J. Phys. G Nucl. Part. Phys 30 S1213 (2004))
  
  if (fDebug) cout<<"****AliFlowAnalysisWithLeeYangZeros::GetGrtheta()****"<<endl;
  
  
  TComplex fG = TComplex::One();
  Double_t fOrder =  2.;
  Double_t fWgt = 1.;
  
  Int_t fNumberOfTracks = fEvent->NumberOfTracks();
  
  for (Int_t i=0;i<fNumberOfTracks;i++) //loop over tracks in event
    {
      fTrack = fEvent->GetTrack(i) ; 
      if (fTrack){
	if (fTrack->UseForIntegratedFlow()) {
	  Double_t fPhi = fTrack->Phi();
	  Double_t fGIm = fR * fWgt*cos(fOrder*(fPhi - fTheta));
	  TComplex fGi(1., fGIm);
	  fG *= fGi;     //product over all tracks
	}
      }
      else {cerr << "no particle pointer !!!"<<endl;}
    }//loop over tracks
  
  return fG;
  
} 


//-----------------------------------------------------------------------   
TComplex AliFlowAnalysisWithLeeYangZeros::GetDiffFlow(AliFlowEventSimple* fEvent, Double_t fR0, Int_t theta) 
{
  // Sum for the denominator for diff. flow for the Product Generating Function for LeeYangZeros method
  // PG Eq. 9 (J. Phys. G Nucl. Part. Phys 30 S1213 (2004))
  // Also for v1 mixed harmonics: DF Eq. 5
  // It is the deriverative of Grtheta at r0 divided by Grtheta at r0
  
  if (fDebug) cout<<"****AliFlowAnalysisWithLeeYangZeros::GetGrtheta()****"<<endl;
  
  TComplex fG = TComplex::One();
  TComplex fdGr0(0.,0.);
  Double_t fOrder =  2.;
  Double_t fWgt = 1.;
  
  Int_t fNumberOfTracks = fEvent->NumberOfTracks();
  
  Int_t fNtheta = AliFlowLYZConstants::kTheta;
  Double_t fTheta = ((double)theta/fNtheta)*TMath::Pi()/fOrder;
  
  for (Int_t i=0;i<fNumberOfTracks;i++) //loop over tracks in event
    {
      fTrack = fEvent->GetTrack(i) ;  
      if (fTrack){
	if (fTrack->UseForDifferentialFlow()) {
	  Double_t fPhi = fTrack->Phi();
	  Double_t fCosTerm = fWgt*cos(fOrder*(fPhi - fTheta));
	  //GetGr0theta
	  Double_t fGIm = fR0 * fCosTerm;
	  TComplex fGi(1., fGIm);
	  fG *= fGi;     //product over all tracks
	  //GetdGr0theta
	  TComplex fCosTermComplex(1., fR0*fCosTerm);
	  fdGr0 +=(fCosTerm / fCosTermComplex);  //sum over all tracks
	}
      } //if particle
      else {cerr << "no particle!!!"<<endl;}
    }//loop over tracks
  
  
  for (Int_t i=0;i<fNumberOfTracks;i++) 
    {
      fTrack = fEvent->GetTrack(i) ;  
      if (fTrack){
	if (fTrack->UseForDifferentialFlow()) {
	  Double_t fEta = fTrack->Eta();
	  Double_t fPt = fTrack->Pt();
	  Double_t fPhi = fTrack->Phi();
	  Double_t fCosTerm = cos(fOrder*(fPhi-fTheta));
	  TComplex fCosTermComplex(1.,fR0*fCosTerm);
	  TComplex fNumer = fG*fCosTerm/fCosTermComplex;  //PG Eq. 9
	  fHist2[theta]->Fill(fEta,fPt,fNumer);
	}
      } //if particle
      else {cerr << "no particle pointer!!!"<<endl;}
    }//loop over tracks
  
  TComplex fDenom = fG*fdGr0;
  return fDenom;
  
} 

//----------------------------------------------------------------------- 

