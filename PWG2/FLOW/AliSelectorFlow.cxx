/* $Id$ */
/* derived from AliSelector.cxx,v 1.17 2006/08/31 jgrosseo Exp $ */

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

#include "AliSelectorFlow.h"

#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TRegexp.h>
#include <TTime.h>
#include <TFriendElement.h>
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TTimeStamp.h>
#include <TMath.h>

#include "/localstore/alice/alice_new/AliRoot_Head/include/AliLog.h"      //  
#include "/localstore/alice/alice_new/AliRoot_Head/include/AliESD.h"      //  
#include "/localstore/alice/alice_new/AliRoot_Head/include/AliESDtrack.h" //  
#include "/localstore/alice/alice_new/AliRoot_Head/include/AliESDv0.h"    //  

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std; //required for resolving the 'cout' symbol

#include "AliFlowEvent.h"
#include "AliFlowTrack.h"
#include "AliFlowV0.h"
#include "AliFlowConstants.h"

ClassImp(AliSelectorFlow)

//-----------------------------------------------------------------------

AliSelectorFlow::AliSelectorFlow() :
  TSelector(),
  fTree(0),
  fESD(0),
  fCountFiles(0),
  fKineFile(0)
{
  //
  // Constructor. Initialization of pointers
  //
  
  fFlowEventFileName = "flowEvtS.root" ;    

  fEtrkLow = 0.01  ;
  fEtrkHig = 100. ;
  fHitsTrk = 1 ;
}

//-----------------------------------------------------------------------

AliSelectorFlow::~AliSelectorFlow()
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

void AliSelectorFlow::CheckOptions()
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

void AliSelectorFlow::Begin(TTree*)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  cout << " HERE I begin !!! " << endl ; cout << endl ;

  CheckOptions();

  AliDebug(AliLog::kDebug, "============BEGIN===========");

 // Opens OutputFile
 pFlowfile = new TFile(fFlowEventFileName.Data(),"RECREATE") ;
 pFlowfile->cd() ; 
 cout << " OUTput File Name : " << fFlowEventFileName.Data() << endl ; 
  
 // Resets counters
 fGoodTracks = 0  ;  
 fGoodV0s    = 0  ;	
 fGoodTracksEta = 0  ;
 fPosiTracks = 0  ;  
 fNegaTracks = 0  ;  
 for(Int_t bb=0;bb<5;bb++) { fBayesianAll[bb] = 0 ; } ;
 fSumAll = 0 ;
 fCutEvts = 0  ;
 fCutTrks = 0  ;
 fCutV0s  = 0  ; 
 
 fLoopV0 = kFALSE ;     // kTRUE ;     // something is wrong with them ... too many!?
 fDoNothing = kFALSE ;  // kFALSE ;
 fOnFlyAnalysis = kFALSE ;  
}

//-----------------------------------------------------------------------

void AliSelectorFlow::SlaveBegin(TTree* tree)
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

void AliSelectorFlow::Init(TTree *tree)
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

Bool_t AliSelectorFlow::Notify()
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

Bool_t AliSelectorFlow::Process(Long64_t entry)
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

  if (!fTree)
  {
    AliDebug(AliLog::kError, "ERROR: fTree is 0.");
    return kFALSE;
  }

  fEventNumber = entry ;
  fTree->GetTree()->GetEntry(fEventNumber) ;

  if(fESD) { AliDebug(AliLog::kDebug, Form("ESD: We have %d tracks.", fESD->GetNumberOfTracks())); }
  // cout << " event !!! " << entry << endl ;

  fRunID = fESD->GetRunNumber() ;
  fEventNumber = fESD->GetEventNumber() ;
  fNumberOfTracks = fESD->GetNumberOfTracks() ;
  fNumberOfV0s = fESD->GetNumberOfV0s() ;

  cout << " *evt n. " << fEventNumber << " (run " << fRunID << ") " << endl ;
  cout << "  tracks: " << fNumberOfTracks << " ,   v0s " << fNumberOfV0s << endl ;

 // Dummy Loop (for testing)
  if(fDoNothing) 
  { 
   pFlowEvent = new AliFlowEvent() ;  
   pFlowEvent->SetRunID(fRunID) ;  		
   pFlowEvent->SetEventID(fEventNumber) ;  	
   pFlowEvent->SetOrigMult((UInt_t)fNumberOfTracks) ;
   return pFlowEvent ; 
  }

 // Evt cuts (dummy) :
  // if(pFlowCuts->CheckEvent(fESD)) { fCutEvts++ ; ...  }

 // Instantiate a new AliFlowEvent
  pFlowEvent = FillFlowEvent(fESD) ;

   if(!pFlowEvent) 
   { 
    cout << "!something bad occurred" << endl ;
    // fCutEvts++ ; continue ; 
    return kFALSE ;
   }
   else {  cout << "!event filled " << endl ; }

  cout << " save event ... " << endl ;

  TString strID = "" ; strID += entry ; 
  pFlowfile->cd() ; 
  pFlowEvent->Write(strID.Data()) ;

  cout << " saved  " << strID.Data() << "                # thanks gesus ! #       " << endl ;
  cout << endl ;  

  if(fOnFlyAnalysis)  { cout << "seeeeeeeeeeeeeeeeeeeeeee :D  " << endl ;} 
   
  return kTRUE;
}

//-----------------------------------------------------------------------

void AliSelectorFlow::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  AliDebug(AliLog::kDebug, "=======SLAVETERMINATE=======");

  DeleteKinematicsFile();
}

//-----------------------------------------------------------------------

void AliSelectorFlow::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliDebug(AliLog::kDebug, "=========TERMINATE==========");

 cout << " file closed . " << endl ;
 pFlowfile->Close() ; 
}

//-----------------------------------------------------------------------

TTree* AliSelectorFlow::GetKinematics()
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

void AliSelectorFlow::DeleteKinematicsFile()
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


//----------------------------------------------------------------------
//*** FLOW SPECIFIC METHODS (to fill the flowEvents) ***
//----------------------------------------------------------------------

AliFlowEvent* AliSelectorFlow::FillFlowEvent(AliESD* fESD)
{
 // Fills the AliFlowEvent . 
 // Calls the methods to fills track & v0 arrays . 
 // Re-Shuffles the tracks (if fShuffle) . 

 fEventNumber = fESD->GetEventNumber() ; 
 fNumberOfTracks = fESD->GetNumberOfTracks() ; 
 fNumberOfV0s = fESD->GetNumberOfV0s() ; 
 fRunID = fESD->GetRunNumber() ;

 fTrackNumber = 0 ;
 fV0Number = 0 ;

 cout << "FillFlowEvent()   " << fEventNumber << "   (run " << fRunID << ") " << endl ; 
   
 // Instantiate a new AliFlowEvent
 pFlowEvent = new AliFlowEvent() ; if(!pFlowEvent) { return 0 ; }

 // Event id 
 pFlowEvent->SetRunID(fRunID) ;  		
 pFlowEvent->SetEventID(fEventNumber) ;  	
 pFlowEvent->SetOrigMult((UInt_t)fNumberOfTracks) ;
						
 // Run information (fixed - ???)
 fMagField = fESD->GetMagneticField() ;
 pFlowEvent->SetMagneticField(fMagField) ;	
 pFlowEvent->SetCenterOfMassEnergy(5500) ;	
 pFlowEvent->SetBeamMassNumberEast(208) ; 	
 pFlowEvent->SetBeamMassNumberWest(208) ; 	

 // Trigger information (now is: ULon64_t - some trigger mask)
 pFlowEvent->SetL0TriggerWord((Int_t)fESD->GetTriggerMask()); 
 
 // Get primary vertex position
 AliESDVertex* pVertex = (AliESDVertex*)fESD->GetVertex() ;
 Double_t position[3] ; 
 pVertex->GetXYZ(position) ;
 pFlowEvent->SetVertexPos((Float_t)position[0],(Float_t)position[1],(Float_t)position[2]) ; 

 // Zero Degree Calorimeter information
 Int_t zdcp = fESD->GetZDCParticipants() ; 
 Float_t zdce[3] ; 
 zdce[0] = fESD->GetZDCN1Energy() + fESD->GetZDCN2Energy(); 
 zdce[1] = fESD->GetZDCP1Energy() + fESD->GetZDCP2Energy() ; 
 zdce[2] = fESD->GetZDCEMEnergy() ;
 pFlowEvent->SetZDCpart(zdcp);  			
 pFlowEvent->SetZDCenergy(zdce[0],zdce[1],zdce[2]);	
 
  cout << "  looping tracks ... (real)  " << fNumberOfTracks << endl ;

  Int_t badTrks = 0 ;
  for(Int_t itr=0;itr<fNumberOfTracks;itr++)   // Loop over tracks
  {
   fTrackNumber = itr ;
   fTrack = fESD->GetTrack(fTrackNumber) ;
   
  // Trk cuts (temporary) :
   Int_t idXt[180] ;  // used for Cluster Map ( see AliESDtrack::GetTPCclusters() )
   Int_t   nHits = fTrack->GetTPCclusters(idXt) ;
   Float_t E	 = fTrack->GetP() ; 
   if((nHits<fHitsTrk) || ((E<fEtrkLow) || (E>fEtrkHig)))  { fCutTrks++ ; continue ; }
   // -or-
   // if(pFlowCuts->CheckTrack(pTrack)) { fCutTracks++ ; ...  }  // Apply track Cuts

  // Instantiate a new AliFlowTrack and fills it
   pFlowTrack = FillFlowTrack(fTrack) ;   
   if(!pFlowTrack) { badTrks++ ; continue ; }

   pFlowEvent->TrackCollection()->Add(pFlowTrack) ;   
   fGoodTracks++ ;				      
  }
  cout << " " << fCutTrks+badTrks << " tracks from event n." << fEventNumber << " have been trown away . " << endl ; 

  if(fLoopV0)
  {
   cout << "  looping v0s ...  " << fNumberOfV0s << endl ;

   Int_t badV0s = 0 ;
   for(Int_t v0n=0;v0n<fNumberOfV0s;v0n++)   // Loop over v0 in AliESD
   {
    fV0Number = v0n ;
    fV0 = fESD->GetV0(fV0Number) ;			

   // Instantiate a new AliFlowV0 and fills it
    pFlowV0 =  FillFlowV0(fV0) ;
    if(!pFlowV0) {  badV0s++ ; continue ; }
    
    pFlowEvent->V0Collection()->Add(pFlowV0) ;
    fGoodV0s++ ; 			
   }
   cout << " " << fCutV0s+badV0s << " v0s from event n." << fEventNumber << " have been trown away . " << endl ; 
  }
 
 // Evt setting stuff
 pFlowEvent->SetCentrality();	
 				
 return pFlowEvent ;
}

//----------------------------------------------------------------------

AliFlowTrack* AliSelectorFlow::FillFlowTrack(AliESDtrack* fTrack)
{

 // cout << "FillFlowTrack() " << endl ; 
 // fTrack->Dump() ;

 TString ntra = "" ; ntra += fTrackNumber ;
 pFlowTrack = new AliFlowTrack(ntra.Data()) ;	

 // ESD particle label (link: KineTree-ESD)
 Int_t label = TMath::Abs(fTrack->GetLabel());
 pFlowTrack->SetLabel(label) ; 			

 //cout  << "      - track " << fTrackNumber << " ... pTrack->GetLabel() = " << label << endl ;

 // signed DCA from ESDtrack
 Float_t xy = 0 ; Float_t z = 0 ; 
 fTrack->GetImpactParameters(xy,z) ; 				// this returns (0,0) !!!
// if(xy == 0 && z == 0)  					// if 0 get it from the ITS method
// {
//   cout << " ------------- no dca !!!" << endl ; 
//   AliITStrackV2 *itsTrack = 0 ;	
//   UInt_t status = fTrack->GetStatus();
//   if((status&AliESDtrack::kITSrefit)==0) 
//   { 
//    xy = 0 ; 	
//   }
//   else  		   		 
//   {
//    itsTrack = new AliITStrackV2(*fTrack) ;
//    if(itsTrack) { xy = itsTrack->GetD() ; }  		// signed DCA from ITStrack
//    delete itsTrack ;
//   }  
// }
 pFlowTrack->SetDcaSigned(xy,z) ; 		    

 // UnConstrained (global) first
 Double_t gD[3] ; 				
 fTrack->GetPxPyPz(gD) ;			
 // -
 Float_t phiGl = (Float_t)Phi(gD) ;  
 if(phiGl<0) { phiGl += 2*TMath::Pi() ; }
 pFlowTrack->SetPhiGlobal(phiGl) ;		
 Float_t ptGl = (Float_t)Pt(gD) ; 
 pFlowTrack->SetPtGlobal(ptGl) ;		
 Float_t etaGl = (Float_t)Eta(gD) ; 
 pFlowTrack->SetEtaGlobal(etaGl) ;		

 // Constrained (NEW)
 Double_t cD[3] ;
 Double_t par1 ; Double_t par2 ;  Double_t par3[3] ;
 if(fTrack->GetConstrainedExternalParameters(par1,par2,par3))
 {
  fTrack->GetConstrainedPxPyPz(cD) ;     
 }
 else { for(Int_t iii=0;iii<3;iii++) { cD[iii] =0 ; } }

 if(Norm(cD)!=0.)   // ConstrainedPxPyPz != 0 if ConstrainedChi2 < something ...
 {  							
  Float_t phi = (Float_t)Phi(cD) ; 
  if(phi<0) { phi += 2*TMath::Pi() ; }
  pFlowTrack->SetPhi(phi) ;                 		
  Float_t pt = (Float_t)Pt(cD) ; 
  pFlowTrack->SetPt(pt) ;          			
  Float_t eta = (Float_t)Eta(cD) ; 
  pFlowTrack->SetEta(eta) ; 				
 
  // number of constrainable tracks with |eta| < Flow::fEtaGood (0.9)
  if(TMath::Abs(eta) < Flow::fEtaGood)  { fGoodTracksEta++ ; }
 }
 else  // in case Constriction impossible for track, fill the UnConstrained (global)
 {
  fUnconstrained++ ; 	
  pFlowTrack->SetPhi(0.) ;
  pFlowTrack->SetPt(0.) ;   
  pFlowTrack->SetEta(0.) ; 
 }	     

 // positive - negative tracks
 Int_t trk_sign = (Int_t)fTrack->GetSign() ; 
 pFlowTrack->SetCharge(trk_sign) ;		
 if(trk_sign>0) 	{ fPosiTracks++ ; }
 else if(trk_sign<0) 	{ fNegaTracks++ ; }
 else 			{ return 0 ; }

 // Tracking parameters (fit , TPC , ITS , dE/dx)
 pFlowTrack->SetChi2(fTrack->GetConstrainedChi2()) ;		
 pFlowTrack->SetTrackLength(fTrack->GetIntegratedLength()) ;	
 // -
 Int_t idXt[180] ; // used for Cluster Map ( see AliESDtrack::GetTPCclusters() )    // old:    Int
 Int_t idX[6] ;    // used for Cluster Map ( see AliESDtrack::GetITSclusters() )    // old:    UInt
 Int_t idxr[130] ; // used for Cluster Map ( see AliESDtrack::GetTRDclusters() )    // old:    UInt
 Int_t nClus = 0 ;	
 Int_t fNFound = 0 ;  							// *!* fNFoundable (in AliTPCtrack) ... added by M.Ianov 
 Double_t pVecAt[3] ; Bool_t boh ; Float_t pAt = (Float_t)Norm(gD) ; 	// to get p at each det.
 // -
 boh = fTrack->GetInnerParam()->GetPxPyPzAt(Flow::fTPCx, fMagField, pVecAt) ;
 pAt = (Float_t)Norm(pVecAt) ;
 nClus = fTrack->GetTPCclusters(idXt) ;
 fNFound = fTrack->GetTPCNclsF() ;  // was 160  		// *!* fNFoundable (in AliTPCtrack) ... added by M.Ianov 
 pFlowTrack->SetMaxPtsTPC(fNFound) ; 	 			
 pFlowTrack->SetFitPtsTPC(nClus) ;                		
 pFlowTrack->SetDedxTPC(fTrack->GetTPCsignal()) ; 		
 pFlowTrack->SetChi2TPC((Float_t)(fTrack->GetTPCchi2())) ; 	
 pFlowTrack->SetPatTPC(pAt) ; 					
 // -
 boh = fTrack->GetInnerParam()->GetPxPyPzAt(Flow::fITSx, fMagField, pVecAt) ;
 pAt = (Float_t)Norm(pVecAt) ;
 nClus = fTrack->GetITSclusters(idX) ;
 fNFound = 6 ;
 pFlowTrack->SetMaxPtsITS(fNFound) ; 	 			
 pFlowTrack->SetFitPtsITS(nClus) ;				
 pFlowTrack->SetDedxITS(fTrack->GetITSsignal()) ;		
 pFlowTrack->SetChi2ITS((Float_t)(fTrack->GetITSchi2())) ; 	
 pFlowTrack->SetPatITS(pAt) ; 					
 // -
 boh = fTrack->GetInnerParam()->GetPxPyPzAt(Flow::fITSx, fMagField, pVecAt) ;
 pAt = (Float_t)Norm(pVecAt) ;
 fNFound = fTrack->GetTRDncls() ;  // was 130
 nClus = fTrack->GetTRDclusters(idxr) ;
 pFlowTrack->SetMaxPtsTRD(fNFound) ;	 			
 pFlowTrack->SetNhitsTRD(nClus) ;				
 pFlowTrack->SetSigTRD(fTrack->GetTRDsignal()) ;		
 pFlowTrack->SetChi2TRD((Float_t)fTrack->GetTRDchi2()) ; 	
 pFlowTrack->SetPatTRD(pAt) ; 					
 // -
 boh = fTrack->GetInnerParam()->GetPxPyPzAt(Flow::fITSx, fMagField, pVecAt) ;
 pAt = (Float_t)Norm(pVecAt) ;
 fNFound = 0 ; if(fTrack->GetTOFCalChannel() > 0) { fNFound = 1 ; }
 nClus = fTrack->GetTOFcluster() ;
 pFlowTrack->SetMaxPtsTOF(fNFound) ;				
 pFlowTrack->SetNhitsTOF(nClus) ;				
 pFlowTrack->SetTofTOF(fTrack->GetTOFsignal()) ;		
 pFlowTrack->SetChi2TOF(fTrack->GetTOFchi2()) ; 		
 pFlowTrack->SetPatTOF(pAt) ; 					
 // -
 Double_t rIn[3]  ; rIn[0] = 0.  ;  rIn[1] = 0. ;  rIn[2] = 0. ;
 Double_t rOut[3] ; rOut[0] = 0. ; rOut[1] = 0. ; rOut[2] = 0. ;
 // -
 fTrack->GetInnerXYZ(rIn) ;			 		
 pFlowTrack->SetZFirstPoint(rIn[2]) ; 
 //fTrack->GetXYZAt(Flow::fTPCx,fMagField,rOut) ;		
 fTrack->GetOuterXYZ(rOut) ;			 		
 pFlowTrack->SetZLastPoint(rOut[2]) ; 
 
 // ESD-P.Id. = 5-vector of Best detectors probabilities for [e , mu , pi , K , p] 
 Double_t trkPid[5] ; fTrack->GetESDpid(trkPid) ;		
 Double_t trkPid6[Flow::nPid] ; 
 for(Int_t bb=0;bb<5;bb++) { trkPid6[bb] = trkPid[bb] ; } 
 trkPid6[5] = 0. ;						// *!* implement P.Id. for Deuterim

 // Bayesian P.Id. method (weighted probabilities for [e , mu , pi , K , p , d])
 Double_t Bsum = 0 ; 
 Double_t bayePid[Flow::nPid] ;    // normalized P.id
 Double_t storedPid[Flow::nPid] ;  // stored P.id
 for(Int_t nB=0;nB<Flow::nPid;nB++)  { Bsum += trkPid6[nB]*Flow::fBayesian[nB] ; }
 if(Bsum)
 {
  for(Int_t nB=0;nB<Flow::nPid;nB++) 
  { 
   bayePid[nB] = trkPid6[nB]*Flow::fBayesian[nB] / Bsum ; 
   storedPid[nB] = trkPid6[nB] ; 
  }
 }
 else { cout << " ERROR - Empty Bayesian Vector !!! " << endl ; }
 
 pFlowTrack->SetElectronPositronProb(storedPid[0]);               
 pFlowTrack->SetMuonPlusMinusProb(storedPid[1]);
 pFlowTrack->SetPionPlusMinusProb(storedPid[2]);
 pFlowTrack->SetKaonPlusMinusProb(storedPid[3]);
 pFlowTrack->SetProtonPbarProb(storedPid[4]);
 pFlowTrack->SetDeuteriumAntiDeuteriumProb(storedPid[5]); 	// *!* implement P.Id. for Deuterim

 // P.id. label given via the weighted prob.
 const Int_t code[]   =  {11,13,211,321,2212,10010020} ;
 Int_t kkk = 2 ; 			// if No id. -> then is a Pi
 Float_t pid_max = bayePid[2] ; 	// (if all equal, Pi probability get's the advantage to be the first)
 for(Int_t iii=0; iii<5; iii++) 
 {
  if(bayePid[iii]>pid_max) { kkk = iii ; pid_max = bayePid[iii] ; }  // !!! Bayesian as well !!!
 }
 fBayesianAll[kkk]++ ; fSumAll++ ; 	// goes on filling the vector of observed abundance 
 //-
 Int_t pdg_code = trk_sign*code[kkk] ;
 pFlowTrack->SetMostLikelihoodPID(pdg_code);
 
 return pFlowTrack ; 
}

//-----------------------------------------------------------------------

AliFlowV0* AliSelectorFlow::FillFlowV0(AliESDv0* fV0)
{
 // Fills the AliFlowV0 .
 // Sets the index of the 2 daughter tracks in the track array . 
  
 // cout << "FillFlowV0()" << endl ; 
 // fV0->Dump() ;

 TString ntra = "" ; ntra += fV0Number ;
 pFlowV0 = new AliFlowV0(ntra.Data()) ;

 Double_t Pxyz[3] ; 		// reconstructed momentum of the V0
 fV0->GetPxPyPz(Pxyz[0],Pxyz[1],Pxyz[2]) ;

 Float_t phi = (Float_t)Phi(Pxyz) ; if(phi<0) { phi += 2*TMath::Pi() ; }
 pFlowV0->SetPhi(phi) ;		
 Float_t pt = (Float_t)Pt(Pxyz) ; 
 pFlowV0->SetPt(pt) ;		
 Float_t eta = (Float_t)Eta(Pxyz) ; 
 pFlowV0->SetEta(eta) ; 	

 Double_t xyz[3] ; 		// reconstructed position of the V0 
 fV0->GetXYZ(xyz[0],xyz[1],xyz[2]) ;	        
 pFlowV0->SetCrossPoint(xyz[0],xyz[1],xyz[2]) ;

 // V0's impact parameter & error
 pFlowV0->SetDca((Float_t)fV0->GetD()) ;	
 //pFlowV0->SetSigma(1.) ;			

 // chi2 of the V0 
 //pFlowV0->SetChi2((Float_t)fV0->GetChi2()) ;	  // AliRoot v4-04-Release
 //pFlowV0->SetChi2((Float_t)fV0->GetChi2V0()) ;  // AliRoot v4-04-Release (but a different one!)	
 // ...when they'll stop annoying me with stupid name changes (see above) I may enable these lines
 pFlowV0->SetChi2(1.) ;

 // P.id. 
 Int_t pdg_code = fV0->GetPdgCode() ;
 pFlowV0->SetMostLikelihoodPID(pdg_code); 	

 // mass 
 pFlowV0->SetVmass((Float_t)fV0->GetEffMass()) ; 
 
 //Int_t pN = fV0->GetPindex() ;
 //Int_t nN = fV0->GetNindex() ;
 AliFlowTrack* pos = 0 ;
 AliFlowTrack* neg = 0 ;
 //if(fMovedTr[pN]>=0) { pos = (AliFlowTrack*)pFlowEvent->TrackCollection()->At(fMovedTr[pN]) ; }
 //if(fMovedTr[nN]>=0) { neg = (AliFlowTrack*)pFlowEvent->TrackCollection()->At(fMovedTr[nN]) ; }
 pFlowV0->SetDaughters(pos,neg) ;

 return pFlowV0 ;
}

//----------------------------------------------------------------------
//*** USEFULL METHODS (to use 3-array) ***
//-----------------------------------------------------------------------

Double_t AliSelectorFlow::Norm(Double_t nu[3])
{ 
 // returns the norm of a double[3] 

 Double_t norm2 = nu[0]*nu[0] + nu[1]*nu[1] + nu[2]*nu[2] ;
 return TMath::Sqrt(norm2) ; 
}

//-----------------------------------------------------------------------

Double_t AliSelectorFlow::Phi(Double_t nu[3])
{
 // returns the azimuthal angle of a double[3] 

 if(nu[0]==0 && nu[1]==0) { return 0. ; }
 else 			  { return TMath::ATan2(nu[1],nu[0]) ; }
}

//-----------------------------------------------------------------------

Double_t AliSelectorFlow::Pt(Double_t nu[3])
{
 // returns the transvers momentum of a double[3] 

 Double_t trans = nu[0]*nu[0] + nu[1]*nu[1] ;
 return TMath::Sqrt(trans) ; 
}

//-----------------------------------------------------------------------

Double_t AliSelectorFlow::Eta(Double_t nu[3])
{
 // returns the PseudoRapidity of a double[3] 
 // if transvers momentum = 0 --> returns +/- 1.000

 Double_t m = Norm(nu) ;
 if(nu[0]!=0 || nu[1]!=0) { return 0.5*TMath::Log((m+nu[2])/(m-nu[2])) ; }
 else     	 	  { return TMath::Sign((Double_t)1000.,nu[2]) ; }
}

//-----------------------------------------------------------------------



