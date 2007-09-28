/* AliSelectorLYZNewMethod.cxx, v1.0 30/07/2007 kolk Exp */
/* derived from AliSelectorFoF.cxx, v1.1 01/02/2007 esimili Exp */
/* derived from AliSelector.cxx,v 1.17 2006/08/31 jgrosseo Exp */

// The class definition in esdV0.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called everytime a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("AliSelector.C")
// Root > T->Process("AliSelector.C","some options")
// Root > T->Process("AliSelector.C+")
//

#include "AliSelectorLYZNewMethod.h"

#include <TSystem.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TTimeStamp.h>
#include <TMath.h>
#include <TVector2.h>
#include <TComplex.h>
#include <TObject.h>
#include <TObjArray.h>
#include <TSelector.h>

#include "AliLog.h"		  
#include "AliESD.h"		  
#include "AliESDtrack.h"  
	   
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std; //required for resolving the 'cout' symbol

#include "AliFlowConstants.h"
#include "AliFlowLYZConstants.h"
#include "AliFlowEvent.h"
#include "AliFlowTrack.h"
#include "AliFlowSelection.h"
#include "AliFlowMaker.h"


ClassImp(AliSelectorLYZNewMethod)

//-----------------------------------------------------------------------

AliSelectorLYZNewMethod::AliSelectorLYZNewMethod() :
  TSelector(),
  fTree(0),
  fESD(0),
  fCountFiles(0),
  fKineFile(0)
 {
  //
  // Constructor. Initialization of pointers
  //
     
}

//-----------------------------------------------------------------------

AliSelectorLYZNewMethod::~AliSelectorLYZNewMethod()
{
  //
  // Destructor
  //

 if (fTree) { fTree->ResetBranchAddresses() ; }

 if (fESD)
 {
   delete fESD;
   fESD = 0;
 }
}

//-----------------------------------------------------------------------

void AliSelectorLYZNewMethod::CheckOptions()
{
  // checks the option string for the debug flag
  
  AliLog::SetClassDebugLevel(ClassName(), AliLog::kInfo);

  TString option = GetOption();

  if (option.Contains("moredebug"))
  {
    printf("Enabling verbose debug mode for %s\n", ClassName());
    AliLog::SetClassDebugLevel(ClassName(), AliLog::kDebug+1);
    AliInfo(Form("Called with option %s.", option.Data()));
  }
  else if (option.Contains("debug"))
  {
    printf("Enabling debug mode for %s\n", ClassName());
    AliLog::SetClassDebugLevel(ClassName(), AliLog::kDebug);
    AliInfo(Form("Called with option %s.", option.Data()));
  }
}

//-----------------------------------------------------------------------

void AliSelectorLYZNewMethod::Begin(TTree*)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  cerr << " HERE I begin !!! " << endl ; cout << endl ;
  

  CheckOptions();

  AliDebug(AliLog::kDebug, "============BEGIN===========");

  // Maker part :
 fFlowMaker = new AliFlowMaker() ;
 // ESD Cuts
 fFlowMaker->SetNHitsCut(1) ;
 fFlowMaker->SetECut(0.01,100.) ; 
 fFlowMaker->PrintCutList() ;


  // AliFlowSelection...
  fFlowSelect = new AliFlowSelection() ;
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
    
  // read input files (read)
  TString firstRunFileName  = "fof_flowLYZAnal_firstrun.root"  ;
  TString secondRunFileName = "fof_flowLYZAnal_secondrun.root" ;

  fFirstRunFile = new TFile(firstRunFileName.Data(),"READ");
  if(!fFirstRunFile || fFirstRunFile->IsZombie()) { cerr << " ERROR: NO first Run file... " << endl ; }
  fSecondRunFile = new TFile(secondRunFileName.Data(),"READ");
  if(!fSecondRunFile || fSecondRunFile->IsZombie()) { cerr << " ERROR: NO Second Run file... " << endl ; }
  cerr<<"....input files read"<<endl;
  
  //input histograms
  h1 = ( TProfile*)fSecondRunFile->Get("Second_FlowProLYZ_ReDtheta_Har2");
  h2 = ( TProfile*)fSecondRunFile->Get("Second_FlowProLYZ_ImDtheta_Har2");
  h3 = ( TH1F*)fSecondRunFile->Get("Control_FlowLYZ_Qtheta"); //only for debugging
  p1 = (TProfile*)fFirstRunFile->Get("First_FlowProLYZ_r0theta_Har2");
  p2 = (TProfile*)fSecondRunFile->Get("Second_FlowProLYZ_VPt_Har2");
  p3 = (TProfile*)fFirstRunFile->Get("First_FlowProLYZ_r0theta_Har2"); //double??
  p4 = (TProfile*)fSecondRunFile->Get("Second_FlowProLYZ_ReDtheta_Har2");
  p5 = (TProfile*)fSecondRunFile->Get("Second_FlowProLYZ_ImDtheta_Har2");

  
  // analysis file (output)
  fOutfile = new TFile("outtestNewMethod.root","RECREATE") ;
  fOutfile->cd() ; 
  
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

 cout << "                ... init is done " << endl ; 
}

//-----------------------------------------------------------------------

void AliSelectorLYZNewMethod::SlaveBegin(TTree* tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  CheckOptions();

  AliDebug(AliLog::kDebug, "=======SLAVEBEGIN========");
  AliDebug(AliLog::kDebug, Form("Hostname: %s", gSystem->HostName()));
  AliDebug(AliLog::kDebug, Form("Time: %s", gSystem->Now().AsString()));

  if (tree != 0) { Init(tree) ; }
}

//-----------------------------------------------------------------------

void AliSelectorLYZNewMethod::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses of the tree
  // will be set. It is normaly not necessary to make changes to the
  // generated code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running with PROOF.

  AliDebug(AliLog::kDebug, "=========Init==========");

  fTree = tree;

  if (fTree == 0)
  {
   AliDebug(AliLog::kError, "ERROR: tree argument is 0.");
   return;
  }

  // Set branch address
  fTree->SetBranchAddress("ESD", &fESD);
  if (fESD != 0) { AliDebug(AliLog::kInfo, "INFO: Found ESD branch in chain.") ; }
}

//-----------------------------------------------------------------------

Bool_t AliSelectorLYZNewMethod::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. Typically here the branch pointers
  // will be retrieved. It is normaly not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed.

  AliDebug(AliLog::kDebug, "=========NOTIFY==========");
  AliDebug(AliLog::kDebug, Form("Hostname: %s", gSystem->HostName()));
  AliDebug(AliLog::kDebug, Form("Time: %s", TTimeStamp(time(0)).AsString()));

  ++fCountFiles;
  if (fTree)
  {
    TFile *f = fTree->GetCurrentFile();
    AliDebug(AliLog::kInfo, Form("Processing %d. file %s", fCountFiles, f->GetName()));
  }
  else
  {
    AliDebug(AliLog::kError, "fTree not available");
  }

  DeleteKinematicsFile();

  return kTRUE;
}

//-----------------------------------------------------------------------

Bool_t AliSelectorLYZNewMethod::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either TTree::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.

  // WARNING when a selector is used with a TChain, you must use
  //  the pointer to the current TTree to call GetEntry(entry).
  //  The entry is always the local entry number in the current tree.
  //  Assuming that fTree is the pointer to the TChain being processed,
  //  use fTree->GetTree()->GetEntry(entry).

  AliDebug(AliLog::kDebug, Form("=========PROCESS========== Entry %lld", entry));

  if(!fTree)
  {
   AliDebug(AliLog::kError, "ERROR: fTree is 0.") ;
   return kFALSE ;
  }

  fEventNumber = entry ;
  fTree->GetTree()->GetEntry(fEventNumber) ;

  if(fESD) { AliDebug(AliLog::kDebug, Form("ESD: We have %d tracks.", fESD->GetNumberOfTracks())); }
  cout << " event !!! " << entry << endl ;

  fRunID = fESD->GetRunNumber() ;
  //  fEventNumber = fESD->GetEventNumber() ;
  fEventNumber = -1 ;
  fNumberOfTracks = fESD->GetNumberOfTracks() ;
  fNumberOfV0s = fESD->GetNumberOfV0s() ;

  cout << " *evt n. " << fEventNumber << " (run " << fRunID << ") " << endl ;
  cout << "  tracks: " << fNumberOfTracks << " ,   v0s " << fNumberOfV0s << endl ;

  // Instantiate a new AliFlowEvent
  cout << " filling the flow event :| " << endl ;
  fFlowEvent = fFlowMaker->FillFlowEvent(fESD) ;
  if(!fFlowEvent) { cout << "! something bad occurred !" << endl ; return kFALSE ; }
  else 		  { cout << "# event done :) " << entry << "     # ok ! #" << endl ; }  
     

   
  // Analysis
  cout<<"---New Lee Yang Zeros Flow Analysis---"<<endl;

  //declare variables
  Float_t cosTerm, sinTerm;
  TComplex fDtheta;
  Int_t fNtheta = AliFlowLYZConstants::kTheta;
      
  fFlowSelect->SetSelection(1);
  fFlowSelect->SetHarmonic(1); //second harmonic

  if(fFlowSelect->Select(fFlowEvent))	 // event selected 
    {
	  	  
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



  delete fFlowEvent; 
  
   
  return kTRUE;
}

//-----------------------------------------------------------------------

void AliSelectorLYZNewMethod::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  AliDebug(AliLog::kDebug, "=======SLAVETERMINATE=======");

  DeleteKinematicsFile();
}

//-----------------------------------------------------------------------

void AliSelectorLYZNewMethod::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliDebug(AliLog::kDebug, "=========TERMINATE==========");

  cout << " Finished ... " << endl ;

  //save histograms in file //temp for testing selector
  fOutfile->cd();
  fOutfile->Write();

  //plot some histograms or store them in a file
  //
  TCanvas *canvas = new TCanvas("canvas","compare v2 vs pt",800,800);
  //canvas->Divide(2,1);
  canvas->cd(1); 
  p2->SetLineColor(3);
  p2->SetLineWidth(2);
  p2->Draw();
  fHistProFlow->Draw("SAME");
  //canvas->cd(2); 
  //fHistProDivide = new TProfile("Divide_VPt_Har2","Divide_VPt_Har2",20,0.,2.);
  //Float_t integrated1 = 0.;
  //Float_t integrated2 = 0.;
  //for (Int_t b=1;b<=20;b++)
  //{
  //    Float_t pt = fHistProFlow->GetBinCenter(b);
  //    Float_t v1 = fHistProFlow->GetBinContent(b);
  //    Float_t v2 = p2->GetBinContent(b);
  //    if (v1==0.) Float_t v3 = 0.;
  //    else Float_t v3 = v2/v1;
  //    fHistProDivide->Fill(pt,v3);
  //    integrated1 += v1;
  //    integrated2 += v2;
  //  }
  //fHistProDivide ->Draw();
  //cout<<"new method gives: "<<integrated1<<" and LYZ gives: "<<integrated2<<". Ratio between them: "<<integrated2/integrated1<<endl;

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
 

  delete fFlowSelect ;
  
  //cout << endl ; 
  return ;
} 

//-----------------------------------------------------------------------

TTree* AliSelectorLYZNewMethod::GetKinematics()
{
  // Returns kinematics tree corresponding to current ESD active in fTree
  // Loads the kinematics from the kinematics file, the file is identified by replacing "AliESDs" to
  // "Kinematics" in the file path of the ESD file. This is a hack, to be changed!

  if (!fKineFile)
  {
    if(!fTree->GetCurrentFile()) { return 0 ; }

    TString fileName(fTree->GetCurrentFile()->GetName());
    fileName.ReplaceAll("AliESDs", "Kinematics");

    // temporary workaround for PROOF bug #18505
    fileName.ReplaceAll("#Kinematics.root#Kinematics.root", "#Kinematics.root");

    AliDebug(AliLog::kInfo, Form("Opening %s", fileName.Data()));

    fKineFile = TFile::Open(fileName);
    if(!fKineFile) { return 0 ; }
  }

  return dynamic_cast<TTree*> (fKineFile->Get(Form("Event%d/TreeK", fTree->GetTree()->GetReadEntry())));
}

//-----------------------------------------------------------------------

void AliSelectorLYZNewMethod::DeleteKinematicsFile()
{
  //
  // Closes the kinematics file and deletes the pointer.
  //

  if (fKineFile)
  {
    fKineFile->Close();
    delete fKineFile;
    fKineFile = 0;
  }
}

//-----------------------------------------------------------------------


