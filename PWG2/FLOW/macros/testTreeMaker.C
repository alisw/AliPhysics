/////////////////////////////////////////////////////////////
//
// $Id$
//
// Author: Emanuele Simili
//
/////////////////////////////////////////////////////////////
//
// Description: 
//    AliRoot macro to make an AliFlowEvent tree from AliESDs 
//
/////////////////////////////////////////////////////////////

#include <vector>
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TFile.h"
#include "TObjArray"
#include "TTree"
#include "TStopwatch.h"

using namespace std; //required for resolving the 'cout' symbol

void testTreeMaker(TString output = "flowEvtsTree.root")
{
 cout << " . Here the new flow maker (TTree) ... " << endl ;
 cout << endl ;

 bool kOne = kFALSE ;

 TStopwatch timer;
 timer.Start();

 gSystem->Load("libPWG2flow.so");

 // esd chain //

 TString fESDfileName = "AliESDs.root" ;
 TString fESDtree = "esdTree" ;
 TString fESDbranch = "ESD" ;

 TChain * pESDchain = new TChain(fESDtree.Data()) ;

 if(kOne)
 {
  pESDchain->Add(fESDfileName.Data()) ; // nFiles++ ;
 }
 else
 {
  pESDchain->Add("1/AliESDs.root") ; // nFiles++ ;
  pESDchain->Add("2/AliESDs.root") ; // nFiles++ ;
  pESDchain->Add("3/AliESDs.root") ; // nFiles++ ;
  pESDchain->Add("4/AliESDs.root") ; // nFiles++ ;
  pESDchain->Add("5/AliESDs.root") ; // nFiles++ ;
  pESDchain->Add("6/AliESDs.root") ; // nFiles++ ;
  pESDchain->Add("7/AliESDs.root") ; // nFiles++ ;
  pESDchain->Add("8/AliESDs.root") ; // nFiles++ ;
  pESDchain->Add("9/AliESDs.root") ; // nFiles++ ;
 }

 Int_t fNumberOfEvents = (Int_t)pESDchain->GetEntries() ;
 cout << " tot. " << fNumberOfEvents << " events in the TChain ... " << endl ; cout << endl ;

 AliESD * pEsd = 0 ;
 pESDchain->SetBranchAddress(fESDbranch.Data(),&pEsd) ;

 // flow maker //

 AliFlowMaker * flowMaker = new  AliFlowMaker() ;
 // cuts, etc.
 flowMaker->SetNHitsCut(1) ;
 flowMaker->SetECut(0.01,100.) ;
 //flowMaker->SetLabelCut(..,..) ;
 flowMaker->PrintCutList() ;

 // output file //

 TFile * flowFile = new TFile(output.Data(),"RECREATE") ;
 flowFile->cd() ; 

 // AliFlowEvents tree //

 TTree * flowTree = new TTree("flowTree","flowTree",2) ;
 AliFlowEvent * flowEvt = new AliFlowEvent() ;
 TBranch * flowEventBranch = flowTree->Branch("flowBranch","AliFlowEvent",&flowEvt,32000,2);

 // loop //

 Int_t evtN = 0 ;
 for(evtN=0;evtN<fNumberOfEvents;evtN++)
 {
  pESDchain->GetEntry(evtN,1) ;

  Int_t evtNN = pEsd->GetEventNumber() ;
  Int_t nTrk = pEsd->GetNumberOfTracks() ;
  Int_t nV0s = pEsd->GetNumberOfV0s() ;
  cout << endl ; cout << " Event " << evtN << "  ( " << evtNN << " )  : " << nTrk << " tracks  &  " << nV0s << " v0s ." << endl ;

  flowEvt = flowMaker->FillFlowEvent(pEsd) ;
  TString evtID = "" ; evtID += evtN ; 
  flowEvt->SetName(evtID.Data()) ;
  
  cout << " Event filled " << flowEvt << " (" << flowEvt->GetName() << ")  ... " << endl ;
  cout << endl ; cout << " trks : " << flowEvt->TrackCollection()->GetEntries() << endl ;
  // flowEvt->Dump() ; cout << endl ;

  flowTree->Fill() ;
  cout <<  " Event " << evtN << "  ( " << evtID.Data() << " )  -  on the Tree (" << flowTree->GetName() << ") ." << endl;

  //delete flowEvt ;
 }
 flowFile->Write() ; 

 flowFile->Close() ; 
 cout <<  " File : " << flowFile->GetName() << "  Closed successfully  (" << evtN << " events saved) ." << endl ;
 cout << endl ; 
 
 //delete flowEvt ;
 //delete flowTree ;
 //delete flowFile ;

 cout << " Finished ... " << endl ;
 cout << "  nTracks:  " << flowMaker->GetNgoodTracks() << endl ;   
 cout << "  nV0s:  " << flowMaker->GetNgoodV0s()  << endl ;  	     
 cout << "  nTracks (|eta|<0.5):  " << flowMaker->GetNgoodTracksEta() << endl ; 
 cout << "  nTracks+:  " << flowMaker->GetNposiTracks() << endl ; 	     
 cout << "  nTracks-:  " << flowMaker->GetNnegaTracks() << endl ; 	     
 cout << "  nTracks unconstrained:  " << flowMaker->GetNunconstrained() << endl ; 	 
 cout << "  Bayesian :  " ; 
 for(int ii=0;ii<5;ii++) { cout << flowMaker->GetBayesianNorm(ii) << "   " ; } 
 cout << " . " << endl ; 

 timer.Stop() ;
 cout << endl ;
 timer.Print() ;
 cout << " . here it was (treeMaker) ... " << endl ;  //juice!
 cout << endl ;

 // break ;

}
