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

#include "Riostream.h"
#include "AliFlowLeeYangZerosMaker.h"
#include "AliFlowEvent.h"
#include "AliFlowSelection.h"
#include "AliFlowConstants.h"     //??
#include "AliFlowLYZHist1.h"
#include "AliFlowLYZHist2.h"
#include "AliFlowLYZConstants.h"  //??
//#include "AliFlowLYZSummary.h"

#include "TMath.h"
#include "TComplex.h" 
#include "TProfile.h"
#include "TObjArray.h"
#include "TFile.h"
#include "TVector.h"
#include "TVector2.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

class TTree;
class TH1F;
class TH1D;
class TVector3;
class TProfile2D;
class TObject;

//class Riostream; //does not compile
//class TMath;     //does not compile
//class TVector;   //does not compile

//Description: Maker to analyze Flow using the LeeYangZeros method  
//             Equation numbers are from Big Paper (BP): Nucl. Phys. A 727, 373 (2003)
//             Practical Guide (PG):    J. Phys. G: Nucl. Part. Phys. 30, S1213 (2004)  
//             Adapted from StFlowLeeYangZerosMaker.cxx           
//             by Markus Oldenberg and Art Poskanzer, LBNL        
//             with advice from Jean-Yves Ollitrault and Nicolas Borghini   
//
//Author: Naomi van der Kolk (kolk@nikhef.nl)



ClassImp(AliFlowLeeYangZerosMaker)

 //-----------------------------------------------------------------------
 
  AliFlowLeeYangZerosMaker::AliFlowLeeYangZerosMaker():
    fFirstRun(kTRUE),
    fUseSum(kTRUE),
    fDebug(kFALSE),
    fHistFileName(0),
    fHistFile(0),
    fSummaryFile(0),
    firstRunFile(0)

{
  //default constructor
  if (fDebug) cout<<"****AliFlowLeeYangZerosMaker::AliFlowLeeYangZerosMaker default constructor****"<<endl;

  fFlowSelect = new AliFlowSelection();
  if (fDebug) { cerr<<"****fFlowSelect in constructor AliFlowLeeYangZerosMaker ("<<fFlowSelect<<")****"<<endl;}
  // output file (histograms)
  TString fHistFileName = "flowLYZAnalysPlot.root" ;
}
 

AliFlowLeeYangZerosMaker::AliFlowLeeYangZerosMaker(const AliFlowSelection* flowSelect):
  fFirstRun(kTRUE),
  fUseSum(kTRUE),
  fDebug(kFALSE),
  fHistFileName(0),
  fHistFile(0),
  fSummaryFile(0),
  firstRunFile(0)
{
  //custum constructor
  if (fDebug) cout<<"****AliFlowLeeYangZerosMaker::AliFlowLeeYangZerosMaker custum constructor****"<<endl;

  if(flowSelect) { fFlowSelect = new AliFlowSelection(*flowSelect); }
  else 	  { 
    fFlowSelect = new AliFlowSelection() ; 
    if (fDebug) cerr<<"****fFlowSelect in constructor AliFlowLeeYangZerosMaker ("<<fFlowSelect<<")****"<<endl;
  }
  // output file (histograms)
  TString fHistFileName = "flowLYZAnalysPlot.root" ;
}
 
 //-----------------------------------------------------------------------


 AliFlowLeeYangZerosMaker::~AliFlowLeeYangZerosMaker() 
 {
   //default destructor
   if (fDebug) cout<<"****~AliFlowLeeYangZerosMaker****"<<endl;
   delete fHistFile;
 }
 
 //-----------------------------------------------------------------------

Bool_t AliFlowLeeYangZerosMaker::Init() 
{
  //init method 
  if (fDebug) cout<<"****AliFlowLeeYangZerosMaker::Init()****"<<endl;

  // Open output files (->plots)
  fHistFile = new TFile(fHistFileName.Data(), "RECREATE");
  //fHistFile->cd() ;  //all histograms will be saved in this file

  //for each harmonic ???
  fQsum.Set(0.,0.);
  fQ2sum = 0.;
  
  // Book histograms
  Int_t fNtheta = AliFlowLYZConstants::kTheta;
  
  
  //for control histograms
  fHistOrigMult = new TH1F("Control_FlowLYZ_OrigMult", "Control_FlowLYZ_OrigMult",1000, 0., 10000.);
  fHistOrigMult->SetXTitle("Original Multiplicity");
  fHistOrigMult->SetYTitle("Counts");

  fHistMult = new TH1F("Control_FlowLYZ_Mult", "Control_FlowLYZ_Mult",1000, 0., 10000.);
  fHistMult->SetXTitle("Multiplicity from selection");
  fHistMult->SetYTitle("Counts");

  fHistQ = new TH1F("Control_FlowLYZ_Q","Control_FlowLYZ_Q",500, 0., 10.);
  fHistQ->SetXTitle("Qvector");
  fHistQ->SetYTitle("Counts");

  fHistPt = new TH1F("Control_FlowLYZ_Pt","Control_FlowLYZ_Pt",200, 0., 10.);
  fHistPt->SetXTitle("Pt (GeV/c)");
  fHistPt->SetYTitle("Counts");

  fHistPhi = new TH1F("Control_FlowLYZ_Phi","Control_FlowLYZ_Phi",70, 0., 7.);
  fHistPhi->SetXTitle("Phi");
  fHistPhi->SetYTitle("Counts");

  fHistEta = new TH1F("Control_FlowLYZ_Eta","Control_FlowLYZ_Eta",40, 0., 2.);
  fHistEta->SetXTitle("Eta");
  fHistEta->SetYTitle("Counts");

  fHistQtheta = new TH1F("Control_FlowLYZ_Qtheta","Control_FlowLYZ_Qtheta",50,-1000.,1000.);
  fHistQtheta->SetXTitle("Qtheta");
  fHistQtheta->SetYTitle("Counts");
  
  if (fFirstRun){
    //for first loop over events
    fHistProR0thetaHar1  = new TProfile("First_FlowProLYZ_r0theta_Har1","First_FlowProLYZ_r0theta_Har1",fNtheta,-0.5,fNtheta-0.5);
    fHistProR0thetaHar1->SetXTitle("#theta");
    fHistProR0thetaHar1->SetYTitle("r_{0}^{#theta}");
  
    fHistProR0thetaHar2  = new TProfile("First_FlowProLYZ_r0theta_Har2","First_FlowProLYZ_r0theta_Har2",fNtheta,-0.5,fNtheta-0.5);
    fHistProR0thetaHar2->SetXTitle("#theta");
    fHistProR0thetaHar2->SetYTitle("r_{0}^{#theta}");

    fHistProVthetaHar1  = new TProfile("First_FlowProLYZ_Vtheta_Har1","First_FlowProLYZ_Vtheta_Har1",fNtheta,-0.5,fNtheta-0.5);
    fHistProVthetaHar1->SetXTitle("#theta");
    fHistProVthetaHar1->SetYTitle("V_{n}^{#theta}");

    fHistProVthetaHar2  = new TProfile("First_FlowProLYZ_Vtheta_Har2","First_FlowProLYZ_Vtheta_Har2",fNtheta,-0.5,fNtheta-0.5);
    fHistProVthetaHar2->SetXTitle("#theta");
    fHistProVthetaHar2->SetYTitle("V_{n}^{#theta}");

    fHistProVR0 = new TProfile("First_FlowProLYZ_vR0","First_FlowProLYZ_vR0",2,0.5,2.5,-100.,100.);
    fHistProVR0->SetXTitle("Harmonic");
    fHistProVR0->SetYTitle("v integrated from r0 (%)");

    fHistVR0 = new TH1D("First_FlowLYZ_vR0","First_FlowLYZ_vR0",2,0.5,2.5);
    fHistVR0->SetXTitle("Harmonic");
    fHistVR0->SetYTitle("v integrated from r0 (%)");

    fHistProV = new TProfile("First_FlowProLYZ_V","First_FlowProLYZ_V",2,0.5,2.5,-1000.,1000.);
    fHistProV->SetXTitle("Harmonic");
    fHistProV->SetYTitle("v integrated");
  
    //class AliFlowLYZHist1 defines the histograms: fHistProGtheta, fHistProReGtheta, fHistProImGtheta, fHistProR0theta
    for (Int_t j=0;j<AliFlowConstants::kHars;j++)  
      {
	for (Int_t theta=0;theta<fNtheta;theta++)
	  {  
	    fHist1[j][theta]=new AliFlowLYZHist1(theta,j+1);
	  }
      }
  }
  else {
    //for second loop over events
    fHistProReDenomHar1 = new TProfile("Second_FlowProLYZ_ReDenom_Har1","Second_FlowProLYZ_ReDenom_Har1" , fNtheta, -0.5, fNtheta-0.5);
    fHistProReDenomHar1->SetXTitle("#theta");
    fHistProReDenomHar1->SetYTitle("Re(Q^{#theta}e^{ir_{0}^{#theta}Q^{#theta}})");

    fHistProReDenomHar2 = new TProfile("Second_FlowProLYZ_ReDenom_Har2","Second_FlowProLYZ_ReDenom_Har2" , fNtheta, -0.5, fNtheta-0.5);
    fHistProReDenomHar2->SetXTitle("#theta");
    fHistProReDenomHar2->SetYTitle("Re(Q^{#theta}e^{ir_{0}^{#theta}Q^{#theta}})");

    fHistProImDenomHar1 = new TProfile("Second_FlowProLYZ_ImDenom_Har1","Second_FlowProLYZ_ImDenom_Har1" , fNtheta, -0.5, fNtheta-0.5);
    fHistProImDenomHar1->SetXTitle("#theta");
    fHistProImDenomHar1->SetYTitle("Im(Q^{#theta}e^{ir_{0}^{#theta}Q^{#theta}})");

    fHistProImDenomHar2 = new TProfile("Second_FlowProLYZ_ImDenom_Har2","Second_FlowProLYZ_ImDenom_Har2" , fNtheta, -0.5, fNtheta-0.5);
    fHistProImDenomHar2->SetXTitle("#theta");
    fHistProImDenomHar2->SetYTitle("Im(Q^{#theta}e^{ir_{0}^{#theta}Q^{#theta}})");

    fHistProVetaHar1 = new TProfile("Second_FlowProLYZ_Veta_Har1","Second_FlowProLYZ_Veta_Har1",40,-2.,2.);
    fHistProVetaHar1->SetXTitle("rapidity");
    fHistProVetaHar1->SetYTitle("v (%)");

    fHistProVetaHar2 = new TProfile("Second_FlowProLYZ_Veta_Har2","Second_FlowProLYZ_Veta_Har2",40,-2.,2.);
    fHistProVetaHar2->SetXTitle("rapidity");
    fHistProVetaHar2->SetYTitle("v (%)");

    fHistProVPtHar1 = new TProfile("Second_FlowProLYZ_VPt_Har1","Second_FlowProLYZ_VPt_Har1",100,0.,10.);
    fHistProVPtHar1->SetXTitle("Pt");
    fHistProVPtHar1->SetYTitle("v (%)");

    fHistProVPtHar2 = new TProfile("Second_FlowProLYZ_VPt_Har2","Second_FlowProLYZ_VPt_Har2",100,0.,10.);
    fHistProVPtHar2->SetXTitle("Pt");
    fHistProVPtHar2->SetYTitle("v (%)");

    fHistVPtHar2 = new TH1D("Second_FlowLYZ_VPt_Har2","Second_FlowLYZ_VPt_Har2",100,0.,10.);
    fHistVPtHar2->SetXTitle("Pt");
    fHistVPtHar2->SetYTitle("v (%)");

    fHistProReDtheta = new TProfile("Second_FlowProLYZ_ReDtheta_Har2","Second_FlowProLYZ_ReDtheta_Har2",fNtheta, -0.5, fNtheta-0.5);
    fHistProReDtheta->SetXTitle("#theta");
    fHistProReDtheta->SetYTitle("Re(D^{#theta})");

    fHistProImDtheta = new TProfile("Second_FlowProLYZ_ImDtheta_Har2","Second_FlowProLYZ_ImDtheta_Har2",fNtheta, -0.5, fNtheta-0.5);
    fHistProImDtheta->SetXTitle("#theta");
    fHistProImDtheta->SetYTitle("Im(D^{#theta})");


    //class AliFlowLYZHist2 defines the histograms: 
    for (Int_t j=0;j<AliFlowConstants::kHars;j++) 
      {
	for (Int_t theta=0;theta<fNtheta;theta++)
	  {  
	    fHist2[j][theta]=new AliFlowLYZHist2(theta,j+1);
	  }
      }

    //read hists from first run file
    //firstRunFile = new TFile("fof_flowLYZAnal_firstrun.root","READ");  //default is read
    if (firstRunFile->IsZombie()){ //check if file exists
      cout << "Error opening file, run first with fFirstrun = kTRUE" << endl;
      exit(-1);
    } else if (firstRunFile->IsOpen()){
      cout<<"----firstRunFile is open----"<<endl<<endl;
      fHistProVthetaHar1  = (TProfile*)firstRunFile->Get("First_FlowProLYZ_Vtheta_Har1");
      fHistProVthetaHar2  = (TProfile*)firstRunFile->Get("First_FlowProLYZ_Vtheta_Har2");
      fHistProR0thetaHar2  = (TProfile*)firstRunFile->Get("First_FlowProLYZ_r0theta_Har2");
      fHistProV = (TProfile*)firstRunFile->Get("First_FlowProLYZ_V");
    }    
  }
   

  if (fDebug) cout<<"****Histograms initialised****"<<endl;
  if (fDebug) cout<<"****fFlowSelect in Init() "<<fFlowSelect<<"****"<<endl;
   
  fEventNumber = 0; //set event counter to zero
  /*
  if (fUseSum)
    {
      //initialize LYZ summary class
      fLYZSummary = new AliFlowLYZSummary(); 
      fSummaryFile = new TFile("fFlowSummary.root","RECREATE","Flow LYZ summary file");
      fSummaryFile->SetCompressionLevel(1);
      fFlowTree = new TTree("FlowTree", "Flow Summary Tree");
      fFlowTree->SetAutoSave(1000000);  // autosave when 1 Mbyte written
      fFlowTree->Branch("fLYZSummary","AliFlowLYZSummary",&fLYZSummary,25000,99);
    }
  */
  return kTRUE; 
}
 
 //-----------------------------------------------------------------------
 
Bool_t AliFlowLeeYangZerosMaker::Make(AliFlowEvent* fFlowEvent) 
{
  //make method
  if (fDebug) cout<<"****AliFlowLeeYangZerosMaker::Make()****"<<endl;
        
  //get tracks from event
  if (fFlowEvent) {
    fFlowTracks = fFlowEvent->TrackCollection();
    if (fDebug) cout<<"****fFlowSelect in Make() "<<fFlowSelect<<"****"<<endl;
    if (fDebug) fFlowSelect->PrintSelectionList() ;
    if (fDebug) fFlowSelect->PrintList() ;
     
    if (fFlowSelect && fFlowSelect->Select(fFlowEvent))   // check if event is selected
      {  
	fFlowTracks = fFlowEvent->TrackCollection();     //get tracks from event
	fFlowEvent->SetSelections(fFlowSelect) ;		// does the selection of tracks for r.p. calculation (sets flags in AliFlowTrack)
	if (fFirstRun){
	  MakeControlHistograms(fFlowEvent);
	  FillFromFlowEvent(fFlowEvent);
	}
	else {
	  MakeControlHistograms(fFlowEvent);
	  SecondFillFromFlowEvent(fFlowEvent);
	}
      }
  }
  else {
    cout<<"##### FlowLeeYangZero: FlowEvent pointer null"<<endl;
    return kFALSE;
  }

  fEventNumber++;
     
  return kTRUE; 
}


  //-----------------------------------------------------------------------     
 Bool_t AliFlowLeeYangZerosMaker::Finish() 
{
  //finish method
  if (fDebug) cout<<"****AliFlowLeeYangZerosMaker::Finish()****"<<endl; 
  
 //define variables for both runs
  Float_t  fR0 = 0;
  Float_t  fv = 0;
  Float_t  fVtheta = 0; 
  Float_t  fSigma2 = 0;
  Float_t  fChi= 0;
  Float_t  fJ01 = 2.405; 
  Int_t fNtheta = AliFlowLYZConstants::kTheta;
  Float_t  fV = 0; 
  
  if (fFirstRun){
    for (Int_t j=0;j<AliFlowConstants::kHars;j++)    
      //loop over harmonics j=0 ->first harmonic, j=1 ->second harmonic
      {
	Float_t  fMeanMult = fHistMult->GetMean();
	//cerr<<"fMeanMult = "<<fMeanMult<<endl;

	for (Int_t theta=0;theta<fNtheta;theta++)
	  {
	    //get the first minimum r0
	    fHist1[j][theta]->FillGtheta();
	    fR0 = fHist1[j][theta]->GetR0();
	    //cerr<<"fR0 = "<<fR0<<endl;
	   	   
	    //calculate integrated flow
	    if (fR0!=0) fVtheta = fJ01/fR0;
	    if (fMeanMult!=0.) 
	      {
		fv = fVtheta/fMeanMult;
		cerr<<"fv = "<<fv<<endl;
	      }
	   

	    //fill the histograms
	    if (j==0) 
	      {
		fHistProR0thetaHar1->Fill(theta,fR0);
		fHistProVthetaHar1->Fill(theta,fVtheta);
		fHistProV->Fill(j+1,fVtheta);    //profile takes care of the averaging over theta.
	      }
	    if (j==1) 
	      {
		fHistProR0thetaHar2->Fill(theta,fR0);
		fHistProVthetaHar2->Fill(theta,fVtheta);
		fHistProV->Fill(j+1,fVtheta);   //profile takes care of the averaging over theta.
	      }
	     
	    fHistProVR0->Fill(j+1,fv*100);    //*100 to get %     //profile takes care of the averaging over theta.              
	  
	  } //end of loop over theta
       
	//sigma2 and chi 
	if (j==1)   //second harmonic only temporarily
	  { 
	    if (fEventNumber!=0) {
	      fQsum /= fEventNumber;
	      //cerr<<"fQsum.X() = "<<fQsum.X()<<endl;
	      //cerr<<"fQsum.Y() = "<<fQsum.Y()<<endl;
	      fQ2sum /= fEventNumber;
	      //cerr<<"fQ2sum = "<<fQ2sum<<endl;
	      fV = fHistProV->GetBinContent(j+1);
	      fSigma2 = fQ2sum - TMath::Power(fQsum.X(),2.) - TMath::Power(fQsum.Y(),2.) - TMath::Power(fV,2.);  //BP eq. 62
	      if (fSigma2>0) fChi = fV/TMath::Sqrt(fSigma2);
	      else fChi = -1.;
	      cerr<<"fV = "<<fV<<" and chi = "<<fChi<<endl;
	    }
	   
	  } //j==1

      } //end of loop over harmonics

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
	 fAmincomb*TMath::BesselJ0(2*fJ01)*TMath::Cos(fTheta/2));
    }
    fRelErr2comb /= fNtheta;
    Double_t fRelErrcomb = TMath::Sqrt(fRelErr2comb);
    cerr<<"fRelErrcomb = "<<fRelErrcomb<<endl;         

    //copy content of profile into TH1D and add error
    for(Int_t b=0;b<2;b++){
      Double_t fv2pro = fHistProVR0->GetBinContent(b+1);
      fHistVR0->SetBinContent(b+1,fv2pro);
      Double_t fv2Err = fv2pro*fRelErrcomb ; 
      cerr<<"fv2pro +- fv2Err = "<<fv2pro<<" +- "<<fv2Err<<" for bin "<<b+1<<endl;
      fHistVR0->SetBinError(b+1,fv2Err);
    }
   

    if (fDebug) cout<<"****histograms filled****"<<endl;  
    /*
      if (fUseSum)
      {
      if (fSummaryFile->IsOpen())
      {
      fSummaryFile->Write(0,TObject::kOverwrite);
      fSummaryFile->Close();
      }
      }
    */
    
    //save histograms in file //temp for testing selector
    fHistFile->cd();
    fHistFile->Write();

    return kTRUE;
    fEventNumber =0; //set to zero for second round over events
  }  //firstrun
   
 
  else {  //second run
   
    //calculate differential flow
    //declare variables
    Float_t fEta, fPt, fReRatio, fVeta, fVPt;
    Float_t fReDenom = 0;
    Float_t fImDenom = 0; 
    Double_t fR0 = 0;
    TComplex fDenom, fNumer, fDtheta;
    Int_t m = 1;
    TComplex i = TComplex::I();
    Float_t fBesselRatio[3] = {1., 1.202, 2.69};
    Double_t fErrdifcomb = 0.;  //set error to zero
    Double_t fErr2difcomb = 0.; //set error to zero
 
    /*
    //v1 integrated
    for (Int_t j=0;j<AliFlowConstants::kHars;j++)   //loop over harmonics j=0 ->first harmonic, j=1 ->second harmonic
    {
    for (Int_t theta=0;theta<fNtheta;theta++)
    {
    //get the first minimum r0
    //fR0 = GetR0(fHist1[j][theta]->FillGtheta()); 

    fReDenom = fHistProReDenomHar2->GetBinContent(theta+1);
    fImDenom = fHistProImDenomHar2->GetBinContent(theta+1);
    TComplex fDenom(fReDenom,fImDenom);
	  
    //complete!!

    }//end of loop over theta
    }//end of loop over harmonics
    */

    //differential flow
    //for (Int_t j=0;j<AliFlowConstants::kHars;j++)    //loop over harmonics j=0 ->first harmonic, j=1 ->second harmonic
    //{
    Int_t j=1; //temp only harm 2
    //fFlowSelect->SetHarmonic(j);  //not needed here?
     
    for (Int_t theta=0;theta<fNtheta;theta++)  { //loop over theta
      if (j==0) {
	fR0 = fHistProR0thetaHar1->GetBinContent(theta+1);
	fVtheta = 2.405/fR0;
	//fVtheta =  fHistProVthetaHar1->GetBinContent(theta+1);  // BP Eq. 9  -> Vn^theta = j01/r0^theta
	fReDenom = fHistProReDenomHar1->GetBinContent(theta+1);
	fImDenom = fHistProImDenomHar1->GetBinContent(theta+1);
      }
      if (j==1) {
	fR0 = fHistProR0thetaHar2->GetBinContent(theta+1);
	if (fDebug) cerr<<"fR0 = "<<fR0<<endl;
	fVtheta = 2.405/fR0;                                    // BP Eq. 9  -> Vn^theta = j01/r0^theta
	fReDenom = fHistProReDenomHar2->GetBinContent(theta+1);
	fImDenom = fHistProImDenomHar2->GetBinContent(theta+1);
      }
	 
      fDenom(fReDenom,fImDenom);
	 

      //for new method and use by others (only with the sum generating function):
      if (fUseSum) {
	fR0 = fHistProR0thetaHar2->GetBinContent(theta+1); 
	fDtheta = fR0*fDenom;
	fHistProReDtheta->Fill(theta,fDtheta.Re());
	fHistProImDtheta->Fill(theta,fDtheta.Im());
      }
	 
      fDenom *= TComplex::Power(i, m-1);
      //cerr<<"TComplex::Power(i, m-1) = "<<TComplex::Power(i, m-1).Rho()<<endl; //checked ok
	 
      //v as a function of eta	  
      Int_t fEtaBins = AliFlowLYZConstants::kEtaBins;
      for (Int_t be=1;be<=fEtaBins;be++)  {
	fEta = fHist2[j][theta]->GetBinCenter(be);
	fNumer = fHist2[j][theta]->GetfNumer(be);
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
	      	      
	if (j==0) fHistProVetaHar1->Fill(fEta,fVeta);
	if (j==1) fHistProVetaHar2->Fill(fEta,fVeta);
      } //loop over bins be
	 
	//v as a function of Pt
      Int_t fPtBins = AliFlowLYZConstants::kPtBins;
      for (Int_t bp=1;bp<=fPtBins;bp++)  {
	fPt = fHist2[j][theta]->GetBinCenterPt(bp);
	fNumer = fHist2[j][theta]->GetfNumerPt(bp);
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
	//cerr<<"fBesselRatio[m-1] = "<<fBesselRatio[m-1]<<endl;   //checked ok
	//if ( j==1 && theta==0) cerr<<"pt = "<<fPt<<" cerr::fReRatio for pt = "<<fReRatio<<" cerr::fVPt for pt = "<<fVeta<<endl;
	      
	if (j==0) fHistProVPtHar1->Fill(fPt,fVPt);
	if (j==1) fHistProVPtHar2->Fill(fPt,fVPt);
      } //loop over bins bp

    }//end of loop over theta


    //sigma2 and chi (for statistical error calculations)
    if (j==1) {  //second harmonic only temporarily
      
      if (fEventNumber!=0) {
	fQsum /= fEventNumber;
	//cerr<<"fQsum.X() = "<<fQsum.X()<<endl;
	//cerr<<"fQsum.Y() = "<<fQsum.Y()<<endl;
	fQ2sum /= fEventNumber;
	cerr<<"fQ2sum = "<<fQ2sum<<endl;
	fV = fHistProV->GetBinContent(j+1);
	fSigma2 = fQ2sum - TMath::Power(fQsum.X(),2.) - TMath::Power(fQsum.Y(),2.) - TMath::Power(fV,2.);  //BP eq. 62
	if (fSigma2>0) fChi = fV/TMath::Sqrt(fSigma2);
	else fChi = -1.;
	cerr<<"fV = "<<fV<<" and chi = "<<fChi<<endl;
      }
	   
    } //j==1

    
    //copy content of profile into TH1D and add error
    for(Int_t b=0;b<100;b++){
      Double_t fv2pro = fHistProVPtHar2->GetBinContent(b);
      //calculate error
      for (Int_t theta=0;theta<fNtheta;theta++) {
	Double_t fTheta = ((double)theta/fNtheta)*TMath::Pi(); 
	Int_t Nprime = fHist2[j][theta]->GetNprimePt(b);
	//cerr<<"Nprime = "<<Nprime<<endl;
	if (Nprime!=0.) { 
	  Double_t fApluscomb = TMath::Exp((fJ01*fJ01)/(2*fChi*fChi)*
					   TMath::Cos(fTheta));
	  Double_t fAmincomb = TMath::Exp(-(fJ01*fJ01)/(2*fChi*fChi)*
					  TMath::Cos(fTheta));
	  fErr2difcomb += (TMath::Cos(fTheta)/(4*Nprime*TMath::BesselJ1(fJ01)*
						 TMath::BesselJ1(fJ01)))*
	    ((fApluscomb*TMath::BesselJ0(2*fJ01*TMath::Sin(fTheta/2))) - 
	     (fAmincomb*TMath::BesselJ0(2*fJ01*TMath::Cos(fTheta/2))));
	}
      } //loop over theta
      
      if (fErr2difcomb!=0.) {
	fErr2difcomb /= fNtheta;
	fErrdifcomb = TMath::Sqrt(fErr2difcomb)*100;
	//cerr<<"fErrdifcomb = "<<fErrdifcomb<<endl;
      }
      else {fErrdifcomb = 0.;}
	  
      //fill TH1D
      if (j==1) {
	fHistVPtHar2->SetBinContent(b,fv2pro);
	fHistVPtHar2->SetBinError(b,fErrdifcomb);
      }
      //check that error is set
      //if (j==1) { cout<<"difference between calculated error and error in hostogram: "<< fErrdifcomb - fHistVPtHar2->GetBinError(b)<<endl; }
	    
    } //loop over bins b

	    
    //} //end of loop over harmonics         //temporarily out

    //save histograms in file
    fHistFile->cd();
    fHistFile->Write();
    fHistFile->Close();
    //Note that when the file is closed, all histograms and Trees in memory associated with this file are deleted
    if (fDebug) cout<<"****Histograms saved and fHistFile closed, all histograms deleted****"<<endl;
     
    //close the first run file 
    firstRunFile->Close();

     
  } //secondrun
   
  delete fFlowSelect;
  cout<<"----LYZ analysis finished....----"<<endl<<endl;

  return kTRUE;
}


//-----------------------------------------------------------------------
  Bool_t AliFlowLeeYangZerosMaker::MakeControlHistograms(AliFlowEvent* fFlowEvent) 
{ 
  //contol histograms of pt, eta, phi, Q, mult
  if (fDebug) cout<<"****AliFlowLeeYangZerosMaker::MakeControlHistograms()****"<<endl; 
   
  //define variables
  TVector2 fQ;
  Float_t fPt, fPhi, fEta, fQmult;

  if (!fFlowEvent){
    cout<<"##### FlowLeeYangZero: FlowEvent pointer null"<<endl;
    return kFALSE;
  }

  if (!fFlowSelect){
    cout<<"##### FlowLeeYangZero: FlowSelect pointer null"<<endl;
    return kFALSE;
  }

  //set selection and harmonic
  fFlowSelect->SetSelection(1); 
  fFlowSelect->SetHarmonic(1); //second harmonic
     
  //cerr<<"selection in MakeControlHistograms()"<<endl;
  //fFlowSelect->PrintSelectionList() ;
  //fFlowSelect->PrintList() ;

  fFlowTracks = fFlowEvent->TrackCollection();
  Int_t fNumberOfTracks = fFlowTracks->GetEntries();
  fHistOrigMult->Fill(fNumberOfTracks);
  //cerr<<"fNumberOfTracks = "<<fNumberOfTracks<<endl;
  Int_t fMult = fFlowEvent->Mult(fFlowSelect);  // Multiplicity of tracks in the specified Selection
  fHistMult->Fill(fMult);   
  //cerr<<"Mult = "<<fMult<<endl;
  
  fQ = fFlowEvent ->Q(fFlowSelect);
  fQmult = fQ.Mod()/TMath::Sqrt(fNumberOfTracks);
  fHistQ->Fill(fQmult);
   
  Int_t tempmult = 0; //for testing
  for (Int_t i=0;i<fNumberOfTracks;i++) 
    {
      fFlowTrack = (AliFlowTrack*)fFlowTracks->At(i) ;   //get object at array position i
      //if (fFlowSelect->SelectPart(fFlowTrack))           //if track is selected
      if (fFlowSelect->Select(fFlowTrack))  
	{
	  fPt = fFlowTrack->Pt();
	  fHistPt->Fill(fPt);
	  tempmult++;
	  fPhi = fFlowTrack->Phi();
	  if (fPhi<0.) fPhi+=2*TMath::Pi();
	  fHistPhi->Fill(fPhi);
	  fEta = fFlowTrack->Eta();
	  fHistEta->Fill(fEta);
	}
    }
  if (fMult!=tempmult){cerr<<"ERROR: Mult() is not tempmult! "<<fMult<<" :: "<<tempmult<<endl<<endl;}
  //else {cerr<<"Mult()= tempmult "<<fMult<<" :: "<<tempmult<<endl<<endl;}

  return kTRUE; 


}



//-----------------------------------------------------------------------
 
 Bool_t AliFlowLeeYangZerosMaker::FillFromFlowEvent(AliFlowEvent* fFlowEvent) 
{ 
  // Get event quantities from AliFlowEvent for all particles

  if (fDebug) cout<<"****AliFlowLeeYangZerosMaker::FillFromFlowEvent()****"<<endl;
   
  if (!fFlowEvent){
    cout<<"##### FlowLeeYangZero: FlowEvent pointer null"<<endl;
    return kFALSE;
  }
   
  if (!fFlowSelect){
    cout<<"##### FlowLeeYangZero: FlowSelect pointer null"<<endl;
    return kFALSE;
  }


  //define variables
  TComplex fExpo, fGtheta, fGthetaNew, fZ;
  //Int_t m;   
  Int_t fNtheta = AliFlowLYZConstants::kTheta;
  Int_t fNbins = AliFlowLYZConstants::kNbins;
   

  //calculate flow
  fFlowSelect->SetSelection(1); 
      
  for (Int_t j=0;j<AliFlowConstants::kHars;j++)   //loop over harmonics j=0 ->first harmonic, j=1 ->second harmonic
    {
      Int_t m=1;
      fFlowSelect->SetHarmonic(j);
      Float_t fOrder = (double)((j+1)/m);
      //cerr<<"fOrder = "<<fOrder<<endl;

      //get the Q vector from the FlowEvent
      TVector2 fQ = fFlowEvent->Q(fFlowSelect); 
      //for chi calculation:
      if (j==1)  //second harmonic only temporarily
	{
	  fQsum += fQ;
	  fQ2sum += fQ.Mod2();
	  //cerr<<"fQ2sum = "<<fQ2sum<<endl;
	}

      fFlowTracks = fFlowEvent->TrackCollection();

      for (Int_t theta=0;theta<fNtheta;theta++)
	{
	  fTheta = ((float)theta/fNtheta)*TMath::Pi()/fOrder; 
	  //cerr<<"fTheta = "<<fTheta<<endl;
	   
	  //calculate fQtheta = cos(fOrder*(fPhi-fTheta);the projection of the Q vector on the reference direction fTheta
	  fQtheta = GetQtheta(fFlowSelect,fFlowTracks,fTheta);
	  //save fQtheta in AliFlowLYZSummary class

	  //something
	  //AliFlowLYZSummary::SetQtheta(theta,fQtheta);

	  if (j==1 && theta==0) fHistQtheta->Fill(fQtheta);
	  //cerr<<"fQtheta = "<<fQtheta<<endl;
	   
	  for (Int_t bin=1;bin<=fNbins;bin++)
	    {
	      Float_t fR = fHist1[j][theta]->GetBinCenter(bin); //bincentre of bins in histogram
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
		  fGtheta = GetGrtheta(fFlowSelect,fR,fTheta);  //make this function
		  if (fGtheta.Rho2() > 100.) break;
		}

	      fHist1[j][theta]->Fill(fR,fGtheta);              //fill real and imaginary part of fGtheta

	    } //loop over bins
	} //loop over theta 
    } //loop over harmonics
   
  return kTRUE;

          
}

 //-----------------------------------------------------------------------   
 Bool_t AliFlowLeeYangZerosMaker::SecondFillFromFlowEvent(AliFlowEvent* fFlowEvent) 
{ 
  //for differential flow

  if (fDebug) cout<<"****AliFlowLeeYangZerosMaker::SecondFillFromFlowEvent()****"<<endl;
    
  if (!fFlowEvent){
    cout<<"##### FlowLeeYangZero: FlowEvent pointer null"<<endl;
    return kFALSE;
  }
   
  if (!fFlowSelect){
    cout<<"##### FlowLeeYangZero: FlowSelect pointer null"<<endl;
    return kFALSE;
  }

  //define variables
  //TVector2 fQ;
  TComplex fExpo, fDenom, fNumer,fCosTermComplex;
  Float_t  fOrder, fR0, fPhi, fEta, fPt;
  Double_t fCosTerm;
  Double_t m = 1.;
  Int_t fNtheta = AliFlowLYZConstants::kTheta;
   
  //calculate flow
  fFlowSelect->SetSelection(1);

  //cerr<<"selection in SecondFillFromFlowEvent()"<<endl;
  //fFlowSelect->PrintSelectionList() ;
  //fFlowSelect->PrintList() ;

      
  for (Int_t j=0;j<AliFlowConstants::kHars;j++)   //loop over harmonics j=0 ->first harmonic, j=1 ->second harmonic
    {
      m=1;
      fFlowSelect->SetHarmonic(j);
      fOrder = (double)((j+1)/m);

      //get the Q vector from the FlowEvent
      TVector2 fQ = fFlowEvent->Q(fFlowSelect); 
      //for chi calculation:
      if (j==1)  //second harmonic only temporarily
	{
	  fQsum += fQ;
	  fQ2sum += fQ.Mod2();
	  //cerr<<"fQ2sum = "<<fQ2sum<<endl;
	}

      fFlowTracks = fFlowEvent->TrackCollection();
       
      for (Int_t theta=0;theta<fNtheta;theta++)
	{
	  fTheta = ((float)theta/fNtheta)*TMath::Pi()/fOrder;   

	  //calculate fQtheta = cos(fOrder*(fPhi-fTheta);the projection of the Q vector on the reference direction fTheta
	  fQtheta = GetQtheta(fFlowSelect,fFlowTracks,fTheta);
	  //cerr<<"fQtheta for fdenom = "<<fQtheta<<endl;
	   
	  /*
	    if (j==0)  //first harmonic
	    {
	    //denominator for differential v
	    fR0 = fHistProR0thetaHar1->GetBinContent(theta+1);
	    fExpo(0.,fR0*fQtheta);
	    fDenom = fQtheta*(TComplex::Exp(fExpo)); //BP eq 12
	    //cerr<<"fDenom.Re() = "<<fDenom.Re()<<endl;
	    //cerr<<"fDenom.Im() = "<<fDenom.Im()<<endl;

	    //denominator for differential v
	    // ****put in product generating function!!

	       
	    fHistProReDenomHar1->Fill(theta,fDenom.Re());               //fill the real part of fDenom
	    fHistProImDenomHar1->Fill(theta,fDenom.Im());               //fill the imaginary part of fDenom
	    }
	  */

	  if (j==1)  //second harmonic
	    {
	      //denominator for differential v
	      fR0 = fHistProR0thetaHar2->GetBinContent(theta+1);
	      //cerr<<"fR0 = "<<fR0 <<endl;

	      if (fUseSum)                                                    //sum generating function
		{
		  fExpo(0.,fR0*fQtheta);
		  fDenom = fQtheta*(TComplex::Exp(fExpo)); //BP eq 12
		  //loop over tracks in event
		  fFlowTracks = fFlowEvent->TrackCollection();
		  Int_t fNumberOfTracks = fFlowTracks->GetEntries();
		  for (Int_t i=0;i<fNumberOfTracks;i++) 
		    {
		      fFlowTrack = (AliFlowTrack*)fFlowTracks->At(i) ;                //get object at array position i
		      if (fFlowSelect->SelectPart(fFlowTrack))                        //if track is selected
			{
			  fPhi = fFlowTrack->Phi();
			  fCosTerm = cos(m*fOrder*(fPhi-fTheta));
			  //cerr<<"fCosTerm = "<<fCosTerm <<endl;
			  fNumer = fCosTerm*(TComplex::Exp(fExpo));
			  if (fNumer.Rho()==0) {cerr<<"WARNING: modulus of fNumer is zero in SecondFillFromFlowEvent"<<endl;}
			  if (fDebug) cerr<<"modulus of fNumer is "<<fNumer.Rho()<<endl;
			  fEta = fFlowTrack->Eta();  //rapidity
			  fPt = fFlowTrack->Pt();
			  fHist2[j][theta]->Fill(fEta,fPt,fNumer);
			} //if
		    } //loop over tracks
		   
		} //sum
	      else                                                        //product generating function
		{
		  fDenom = GetDiffFlow(fFlowSelect,fR0,theta); 
		   
		}//product
	      
	      fHistProReDenomHar2->Fill(theta,fDenom.Re());               //fill the real part of fDenom
	      fHistProImDenomHar2->Fill(theta,fDenom.Im());               //fill the imaginary part of fDenom
	    }//j==1
   	   	   
	   	       
	}//end of loop over theta

    }//loop over harmonics

  
  return kTRUE;
 
  
 
}
 //-----------------------------------------------------------------------   
 Double_t AliFlowLeeYangZerosMaker::GetQtheta(AliFlowSelection* fFlowSelect, TObjArray* fFlowTracks, Float_t fTheta) 
{
  //calculate Qtheta. Qtheta is the sum over all particles of cos(fOrder*(fPhi-fTheta)) BP eq. 3
  if (fDebug) cout<<"****AliFlowLeeYangZerosMaker::GetQtheta()****"<<endl;

  if (!fFlowSelect){
     cout<<"##### FlowLeeYangZero: FlowSelect pointer null"<<endl;
     return kFALSE;
   }


  Double_t fQtheta = 0.;
  Int_t fHarmonic = fFlowSelect->Har(); 
  Double_t fOrder = (double)(fHarmonic+1);
       
  Int_t fNumberOfTracks = fFlowTracks->GetEntries();
  //cerr<<"GetQtheta::fNumberOfTracks = "<<fNumberOfTracks<<endl;
  
  for (Int_t i=0;i<fNumberOfTracks;i++)                  //loop over tracks in event
    {
      fFlowTrack = (AliFlowTrack*)fFlowTracks->At(i) ;   //get object at array position i
      if(fFlowSelect->Select(fFlowTrack))                //if track is selected  //gives the same number of particles as Mult(fFlowSelect) method
	{
     	  Float_t fPhi = fFlowTrack->Phi();
	  fQtheta += cos(fOrder*(fPhi-fTheta));
	}
     
    }//loop over tracks
  
  return fQtheta;
 
}
 
//-----------------------------------------------------------------------   
TComplex AliFlowLeeYangZerosMaker::GetGrtheta(AliFlowSelection* fFlowSelect, Float_t fR, Float_t fTheta) 
{
  // Product Generating Function for LeeYangZeros method
  // PG Eq. 3 (J. Phys. G Nucl. Part. Phys 30 S1213 (2004))

  if (fDebug) cout<<"****AliFlowLeeYangZerosMaker::GetGrtheta()****"<<endl;
  
  if (!fFlowSelect){
    cout<<"##### FlowLeeYangZero: FlowSelect pointer null"<<endl;
    return kFALSE;
  }
  
  TComplex fG = TComplex::One();
  Int_t fHarmonic = fFlowSelect->Har(); 
  Double_t fOrder = (double)(fHarmonic+1);
  Double_t fWgt = 1.;

  Int_t fNumberOfTracks = fFlowTracks->GetEntries();
  //cerr<<"GetGrtheta::fNumberOfTracks = "<<fNumberOfTracks<<endl; 

  for (Int_t i=0;i<fNumberOfTracks;i++) //loop over tracks in event
    {
      fFlowTrack = (AliFlowTrack*)fFlowTracks->At(i) ;                //get object at array position i
      if (fFlowSelect->SelectPart(fFlowTrack))                        //if track is selected
	{
	  Float_t fPhi = fFlowTrack->Phi();
	  Double_t fGIm = fR * fWgt*cos(fOrder*(fPhi - fTheta));
	  TComplex fGi(1., fGIm);
	  fG *= fGi;     //product over all tracks
	}//if
    }//loop over tracks

  return fG;

  } 


//-----------------------------------------------------------------------   
TComplex AliFlowLeeYangZerosMaker::GetDiffFlow(AliFlowSelection* fFlowSelect, Float_t fR0, Int_t theta) 
{
  // Sum for the denominator for diff. flow for the Product Generating Function for LeeYangZeros method
  // PG Eq. 9 (J. Phys. G Nucl. Part. Phys 30 S1213 (2004))
  // Also for v1 mixed harmonics: DF Eq. 5
  // It is the deriverative of Grtheta at r0 divided by Grtheta at r0

  if (fDebug) cout<<"****AliFlowLeeYangZerosMaker::GetGrtheta()****"<<endl;

  if (!fFlowSelect){
    cout<<"##### FlowLeeYangZero: FlowSelect pointer null"<<endl;
    return kFALSE;
  }


  TComplex fG = TComplex::One();
  TComplex fdGr0(0.,0.);
  Int_t fHarmonic = fFlowSelect->Har(); 
  Double_t fOrder = (double)(fHarmonic+1);
  Double_t fWgt = 1.;

  Int_t fNumberOfTracks = fFlowTracks->GetEntries();
  //cerr<<"GetGrtheta::fNumberOfTracks = "<<fNumberOfTracks<<endl; 

  Int_t fNtheta = AliFlowLYZConstants::kTheta;
  Float_t fTheta = ((float)theta/fNtheta)*TMath::Pi()/fOrder;

  for (Int_t i=0;i<fNumberOfTracks;i++) //loop over tracks in event
    {
      fFlowTrack = (AliFlowTrack*)fFlowTracks->At(i) ;                //get object at array position i
      if (fFlowSelect->SelectPart(fFlowTrack))                        //if track is selected
	{
	  Float_t fPhi = fFlowTrack->Phi();
	  Double_t fCosTerm = fWgt*cos(fOrder*(fPhi - fTheta));
	  //GetGr0theta
	  Double_t fGIm = fR0 * fCosTerm;
	  TComplex fGi(1., fGIm);
	  fG *= fGi;     //product over all tracks
	  //GetdGr0theta
	  TComplex fCosTermComplex(1., fR0*fCosTerm);
	  fdGr0 +=(fCosTerm / fCosTermComplex);  //sum over all tracks
	}//if
    }//loop over tracks

  for (Int_t i=0;i<fNumberOfTracks;i++) 
    {
      fFlowTrack = (AliFlowTrack*)fFlowTracks->At(i) ;                //get object at array position i
      if (fFlowSelect->SelectPart(fFlowTrack))                        //if track is selected
	{
	  Float_t fPhi = fFlowTrack->Phi();
	  Double_t fCosTerm = cos(fOrder*(fPhi-fTheta));
	  TComplex fCosTermComplex(1.,fR0*fCosTerm);
	  
	  TComplex fNumer = fG*fCosTerm/fCosTermComplex;  //PG Eq. 9
	  Float_t fEta = fFlowTrack->Eta();  //rapidity
	  Float_t fPt = fFlowTrack->Pt();
	  fHist2[1][theta]->Fill(fEta,fPt,fNumer);
	}//if
    }//loop over tracks

  TComplex fDenom = fG*fdGr0;
  return fDenom;

  } 

//-------------------------------------------------
