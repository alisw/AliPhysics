/////////////////////////////////////////////////////////////
//
// $Id$
//
// Author: Emanuele Simili
//
/////////////////////////////////////////////////////////////
//
// Description: AliRoot macro to make AliFlowEvents from AliESDs (new way) 
//
/////////////////////////////////////////////////////////////

#include <vector>
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TFile.h"
#include "TObjArray"
#include "TStopwatch.h"

using namespace std; //required for resolving the 'cout' symbol

void testMaker(TString output = "pFlowEvts.root")
{
 cout << " . Here the new flow maker (2007) ... " << endl ;
 cout << endl ;

 bool kOne = kFALSE ;

 TStopwatch timer;
 timer.Start();

 gSystem->Load("libPWG2flow.so");

 // output file //

 TFile * fFlowfile = new TFile(output.Data(),"RECREATE") ;
 //fFlowfile->cd() ; 

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

 // loop //

 Int_t evtN = 0 ;
 AliFlowEvent * flowEvt = 0 ;
 for(evtN=0;evtN<fNumberOfEvents;evtN++)
 {
  pESDchain->GetEntry(evtN,1) ;

  Int_t evtNN = pEsd->GetEventNumber() ;
  Int_t nTrk = pEsd->GetNumberOfTracks() ;
  Int_t nV0s = pEsd->GetNumberOfV0s() ;
  cout << endl ; cout << " Event " << evtN << "  ( " << evtNN << " )  : " << nTrk << " tracks  &  " << nV0s << " v0s ." << endl ;

  flowEvt = flowMaker->FillFlowEvent(pEsd) ;
  cout << " Event filled " << flowEvt << " ... " << endl ;

  TString evtID = "" ; evtID += evtN ; 
  fFlowfile->cd() ; 
  flowEvt->Write(evtID.Data()) ;
  cout <<  " Event " << evtN << "  ( " << evtID.Data() << " )  -  written on disk (" << output << ") ." << endl;
  delete flowEvt ;
 }
 
 fFlowfile->Close() ; 

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
 cout << " . here it was (maker) ... " << endl ;  //juice!
 cout << endl ;

 // break ;

}
