#define AliAnalysisTaskRLLYZNewMethod_cxx
 
#include <iostream>


#include "TChain.h"
#include "TH1.h"
#include "TList.h"
#include "TFile.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TParticle.h"
#include "TComplex.h"
#include "TCanvas.h"

#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliStack.h"
#include <AliHeader.h>
#include <AliGenEventHeader.h>

#include "AliAnalysisTaskRL.h"
#include "AliAnalysisTaskRLLYZNewMethod.h"
#include "AliFlowConstants.h"
#include "AliFlowLYZConstants.h"
#include "AliFlowSelection.h"
#include "AliFlowEvent.h"
#include "AliFlowMaker.h"
#include "AliFlowTrack.h"
//#include "AliFlowLeeYangZerosMaker.h"

//#include "TObjectTable.h"



ClassImp(AliAnalysisTaskRLLYZNewMethod)

 //-----------------------------------------------------------------------
 
 AliAnalysisTaskRLLYZNewMethod::AliAnalysisTaskRLLYZNewMethod(const char *name) :
   AliAnalysisTaskRL(name,""),
   fESD(0)
   {

  // Constructor.
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  DefineInput(1, TList::Class()); 
  DefineInput(2, TList::Class()); 
  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class());  
}

 
 //-----------------------------------------------------------------------


 AliAnalysisTaskRLLYZNewMethod::~AliAnalysisTaskRLLYZNewMethod() 
 {
   //destructor
   
 }
 
//-----------------------------------------------------------------------


void AliAnalysisTaskRLLYZNewMethod::ConnectInputData(Option_t *) {
  // Initialize branches.
  printf("   ConnectInputData of task %s\n", GetName());
  //cerr<<"fESD ("<<fESD<<")"<<endl;
  if (!fESD) {
    //cerr<<"no fESD"<<endl;
    char ** address = (char **)GetBranchAddress(0, "ESD");
    if (address) fESD = (AliESD*)(*address);
    if (!fESD) {
      //cerr<<"still no fESD"<<endl;
      fESD = new AliESD();
      SetBranchAddress(0, "ESD", &fESD);
      cerr<<"new fESD"<<endl;
    }
  }
}

//-----------------------------------------------------------------------
void AliAnalysisTaskRLLYZNewMethod::CreateOutputObjects() {

  
  fFlowMaker = new AliFlowMaker();
  cerr<<"create fFlowMaker ("<<fFlowMaker<<")"<<endl;
  fFlowMaker->SetNHitsCut(1);
  fFlowMaker->SetECut(0.01,100.);
  fFlowMaker->PrintCutList();

  fFlowSelect = new AliFlowSelection();
  cerr<<"create fFlowSelect ("<<fFlowSelect<<")"<<endl;
  // Event Cuts
  fFlowSelect->SetCentralityCut(-1) ;
  fFlowSelect->SetRunIdCut(-1) ;
  // R.P. calculation cuts
  for(int j=0;j<AliFlowConstants::kHars;j++)
    {
      fFlowSelect->SetEtaCut(0., 2., j, 1) ;
      fFlowSelect->SetPtCut(0.1, 10. , j, 1);  
    }
  fFlowSelect->SetConstrainCut(kTRUE) ;
  fFlowSelect->SetDcaGlobalCut(0.,0.1);
  // Correlation analysis cuts (not all of them)
  fFlowSelect->SetEtaPart(-1.1,1.1);
  fFlowSelect->SetPtPart(0.1,10.);   
  fFlowSelect->SetConstrainablePart(kTRUE);
  fFlowSelect->SetDcaGlobalPart(0.,0.1);
  // V0 analysis cuts (not all of them ... they are useless anyway)
  fFlowSelect->SetV0Mass(0.4875,0.5078) ;	 // Mk0 = 0.49765
  fFlowSelect->SetV0SideBands(0.1) ;
  fFlowSelect->SetV0Pt(0.1,10.) ;
  fFlowSelect->SetV0Eta(-2.1,2.1) ;
  // print list :
  //cout << " . Selection for R.P. calculation: " << endl ;
  fFlowSelect->PrintSelectionList() ;
  //cout << " . Selection for correlation analysis: " << endl ;
  fFlowSelect->PrintList() ;
  //cout << " . Selection for V0 analysis: " << endl ;
  fFlowSelect->PrintV0List() ;

  // Get data from input slots 
  fFirstRunFile = (TFile*)GetInputData(1);
  cerr<<"fFirstRunFile ("<<fFirstRunFile<<")"<<endl;
  cerr<<"fFirstRunFile -> IsOpen() = "<<fFirstRunFile -> IsOpen()<<endl;

  fSecondRunFile = (TFile*)GetInputData(2);
  cerr<<"fSecondRunFile ("<<fSecondRunFile<<")"<<endl;
  cerr<<"fSecondRunFile -> IsOpen() = "<<fSecondRunFile -> IsOpen()<<endl;  
     
  //input histograms
  h1 = ( TProfile*)fSecondRunFile->Get("Second_FlowProLYZ_ReDtheta_Har2");
  h2 = ( TProfile*)fSecondRunFile->Get("Second_FlowProLYZ_ImDtheta_Har2");
  h3 = ( TH1F*)fSecondRunFile->Get("Control_FlowLYZ_Qtheta"); //only for debugging
  p1 = (TProfile*)fFirstRunFile->Get("First_FlowProLYZ_r0theta_Har2");
  p2 = (TProfile*)fSecondRunFile->Get("Second_FlowProLYZ_VPt_Har2");
  p3 = (TProfile*)fFirstRunFile->Get("First_FlowProLYZ_r0theta_Har2"); //double??
  p4 = (TProfile*)fSecondRunFile->Get("Second_FlowProLYZ_ReDtheta_Har2");
  p5 = (TProfile*)fSecondRunFile->Get("Second_FlowProLYZ_ImDtheta_Har2");


  //output file
  // ********make output file 
  // analysis file (output)
  fOutfile = new TFile("outtestNewMethodTask.root","RECREATE") ;
  

  // output histograms
  fHistProFlow = new TProfile("NewSecond_FlowProLYZ_VPt_Har2","NewSecond_FlowProLYZ_VPt_Har2",100,0.,10.);
  fHistProFlow->SetXTitle("Pt");
  fHistProFlow->SetYTitle("v (%)");

  fHistQtheta = new TH1F("NewControl_FlowLYZ_Qtheta","NewControl_FlowLYZ_Qtheta",50,-1000.,1000.);
  fHistQtheta->SetXTitle("Qtheta");
  fHistQtheta->SetYTitle("Counts");

  fHistProR0thetaHar2  = new TProfile("NewFirst_FlowProLYZ_r0theta_Har2","NewFirst_FlowProLYZ_r0theta_Har2",5,-0.5,4.5);
  fHistProR0thetaHar2->SetXTitle("#theta");
  fHistProR0thetaHar2->SetYTitle("r_{0}^{#theta}");

  fHistProReDtheta = new TProfile("NewSecond_FlowProLYZ_ReDtheta_Har2","NewSecond_FlowProLYZ_ReDtheta_Har2",5, -0.5, 4.5);
  fHistProReDtheta->SetXTitle("#theta");
  fHistProReDtheta->SetYTitle("Re(D^{#theta})");

  fHistProImDtheta = new TProfile("NewSecond_FlowProLYZ_ImDtheta_Har2","NewSecond_FlowProLYZ_ImDtheta_Har2",5, -0.5, 4.5);
  fHistProImDtheta->SetXTitle("#theta");
  fHistProImDtheta->SetYTitle("Im(D^{#theta})");

  // Needed otherwise pointer to fOutfile not available in Exec task
  PostData(0, fOutfile);        
} 
 
//-----------------------------------------------------------------------
 
void AliAnalysisTaskRLLYZNewMethod::Exec(Option_t *) {

  
  // Get data from input slot 0
  TTree *tinput = (TTree*)GetInputData(0);
  Long64_t ientry = tinput->GetReadEntry();
  if (AliAnalysisTaskRL::GetEntry(ientry) == kFALSE) {
    printf("Couldn't get event from the runLoader\n");
    return;
  }
  
  if (!fESD) {
    cout << "No ESD branch available" << endl;
    return;
  }
  
  cerr<<"fESD ("<<fESD <<") in begin Exec"<<endl;
  cerr<<"number of tracks: "<<fESD->GetNumberOfTracks()<<endl;

  AliStack* stack = GetStack();
  if (!stack) {
    AliDebug(AliLog::kError, "Stack not available");
    // return kFALSE;
  }

  cerr<<"fFlowMaker ("<<fFlowMaker<<")"<<endl;
  cerr<<"fFlowSelect ("<<fFlowSelect<<")"<<endl;
  
  fFlowEvent = new AliFlowEvent() ;
  cerr<<"create fFlowEvent ("<<fFlowEvent<<")"<<endl;
  
  if (!fFlowMaker){cerr<<"no fFlowMaker: nullpointer"<<endl;}
  else { 
    if (!fESD) { cerr<<"no fESD: NULL pointer"<<endl;}
    else {
      cerr<<"fFlowMaker ("<<fFlowMaker<<"), fFlowEvent ("<<fFlowEvent<<") and fESD ("<<fESD<<") available"<<endl;
      fFlowEvent = fFlowMaker->FillFlowEvent(fESD);
      
      if (!fFlowEvent){ cerr<<"no fFlowEvent: NULL pointer"<<endl; }
      else {

	// Analysis
	cout<<"---New Lee Yang Zeros Flow Analysis---"<<endl;

	//declare variables
	Float_t cosTerm, sinTerm;
	TComplex fDtheta;
	Int_t fNtheta = AliFlowLYZConstants::kTheta;
      
	fFlowSelect->SetSelection(1);
	fFlowSelect->SetHarmonic(1); //second harmonic

	if(fFlowSelect->Select(fFlowEvent)){	 // event selected 
	 
	  cerr<<"event selected"<<endl;
	  fFlowEvent->SetSelections(fFlowSelect) ;
	  TObjArray* fFlowTracks = fFlowEvent->TrackCollection();
	  Int_t fNumberOfTracks = fFlowTracks->GetEntries();
	  cosTerm = 0.;
	  sinTerm = 0.;
	  
	  for (Int_t theta=0;theta<fNtheta;theta++)	{
	    Float_t fTheta = ((float)theta/fNtheta)*TMath::Pi()/2;  
	    //Calculate Qtheta 
	    Double_t fQtheta = 0.;
	    for (Int_t i=0;i<fNumberOfTracks;i++)  {    //loop over tracks in event
	      fFlowTrack = (AliFlowTrack*)fFlowTracks->At(i) ;   //get object at array position i
	      if(fFlowSelect->Select(fFlowTrack)) {
		Float_t fPhi = fFlowTrack->Phi();
		fQtheta += cos(2*(fPhi-fTheta));
	      }//track selected
     	    }//loop over tracks

	    if (theta==0) fHistQtheta->Fill(fQtheta);  //is correct!
	    cerr<<"fQtheta = "<<fQtheta<<endl;

	    //get R0
	    Double_t fR0 = p1->GetBinContent(theta+1); 
	    fHistProR0thetaHar2->Fill(theta,fR0); //is correct!
	    cerr<<"fR0 = "<<fR0<<endl;

	    //get Dtheta
	    Double_t ReDtheta = h1->GetBinContent(theta+1);
	    Double_t ImDtheta = h2->GetBinContent(theta+1);
	    fDtheta(ReDtheta,ImDtheta);

	    fHistProReDtheta->Fill(theta,ReDtheta); //is correct!
	    fHistProImDtheta->Fill(theta,ImDtheta); //is correct!
	    cerr<<"Dtheta stored"<<endl;

	    TComplex fExpo(0.,fR0*fQtheta);                  //Complex number: 0 +(i r0 Qtheta)
	    TComplex ratio =(TComplex::Exp(fExpo))/fDtheta;  //(e^(i r0 Qtheta))/Dtheta

	    //sum over theta
	    cosTerm += ratio.Re() * TMath::Cos(2*fTheta);    //Re{(e^(i r0 Qtheta))/Dtheta } cos(2 theta)  
	    sinTerm += ratio.Re() * TMath::Sin(2*fTheta);    //Re{(e^(i r0 Qtheta))/Dtheta } sin(2 theta)

	  }//loop over theta

	  //average over theta
	  cosTerm /= fNtheta;  
	  sinTerm /= fNtheta;
	  cerr<<"cosTerm and sinTerm: "<<cosTerm<<" "<<sinTerm<<endl;

	  //calculate fWR
	  Float_t fWR = TMath::Sqrt(cosTerm*cosTerm + sinTerm*sinTerm);
	  cerr<<"fWR = "<<fWR<<endl;

	  //calculate fRP
	  Float_t fRP = 0.5*TMath::ATan2(sinTerm,cosTerm);
	  cerr<<"fRP = "<<fRP<<endl;

	  //loop over the tracks of the event
	  for (Int_t i=0;i<fNumberOfTracks;i++) 
	    {
	      fFlowTrack = (AliFlowTrack*)fFlowTracks->At(i) ;   //get object at array position i
	      if (fFlowSelect->SelectPart(fFlowTrack))           //if track is selected
		{
		  Float_t fPhi = fFlowTrack->Phi();
		  if (fPhi<0.) fPhi+=2*TMath::Pi();
		  //calculate flow v2:
		  Float_t fv2 = fWR * TMath::Cos(2*(fPhi-fRP));
		  Float_t fPt = fFlowTrack->Pt();
		  //fill histogram
		  fHistProFlow->Fill(fPt,2.405*100*fv2);   //2.405 = j01, was missing from method, do not know where
		}//track selected
	    }//loop over tracks
	  
	}//event selected

      }
    }
  }
  
  
  delete fFlowEvent;
  cerr<<"delete fFlowEvent ("<<fFlowEvent<<")"<<endl;
  
  //PostData(0, fOutFile); 
}

  //--------------------------------------------------------------------    
void AliAnalysisTaskRLLYZNewMethod::Terminate(Option_t *) {
   
  //*************make histograms etc. 
 
  fOutfile->Write();
  PostData(0, fOutfile); 
  //save histograms in file //temp for testing selector
  //fOutfile->cd();
  //fOutfile->Write();
  
  TCanvas *canvas = new TCanvas("canvas","compare v2 vs pt",800,800);
  //canvas->Divide(2,1);
  canvas->cd(1); 
  p2->SetLineColor(3);
  p2->SetLineWidth(2);
  p2->Draw();
  fHistProFlow->Draw("SAME");

  TCanvas *canvas2 = new TCanvas("canvas2","compare Qtheta",800,600);
  canvas2->cd(); 
  h3->SetLineColor(3);
  h3->SetLineWidth(2);
  h3->Draw();
  fHistQtheta->Draw("SAME");


  TCanvas *canvas3 = new TCanvas("canvas3","compare r0",800,600);
  canvas3->cd(); 
  p3->SetLineColor(3);
  p3->SetLineWidth(2);
  p3->Draw();
  fHistProR0thetaHar2->Draw("SAME");

  TCanvas *canvas4 = new TCanvas("canvas4","compare Dtheta",800,600);
  canvas4->Divide(2,1);
  canvas4->cd(1); 
  h1->SetLineColor(3);
  h1->SetLineWidth(2);
  h1->Draw();
  fHistProReDtheta->Draw("SAME");
  canvas4->cd(2); 
  h2->SetLineColor(3);
  h2->SetLineWidth(2);
  h2->Draw();
  fHistProImDtheta->Draw("SAME");
 

  delete fFlowMaker;
  delete fFlowSelect;

  cout<<".....finished"<<endl;
 }

 

