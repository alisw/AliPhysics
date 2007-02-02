/////////////////////////////////////////////////////////////
//
// $Id$
//
// Author: Emanuele Simili
//
/////////////////////////////////////////////////////////////
//
// Description: ROOT macro to perform the fill phi weights from AliFlowEvents (new way)
//
/////////////////////////////////////////////////////////////

#include <vector>
#include <iostream>
#include <fstream>
#include "TVector.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TMath.h"
#include "TObjArray"
#include "TStopwatch.h"

using namespace std; //required for resolving the 'cout' symbol

void testWeight(TString wgts = "pFlowPhiWgt.root")
{
 cout << " . Here the new flow Weighter (2007) ... " << endl ;
 cout << endl ;

 bool kOne = kFALSE ;

 TStopwatch timer;
 timer.Start();

 gSystem->Load("libPhysics.so");
 gSystem->Load("libPWG2.so");
 gSystem->Load("libPWG2flow.so");

 TStopwatch timer;
 timer.Start();


 // selection object (cuts) //

 AliFlowSelection* select = new AliFlowSelection() ;
 
 // Event Cuts
 select->SetCentralityCut(-1) ;
 select->SetRunIdCut(-1) ;

 // R.P. calculation cuts
 for(int j=0;j<Flow::nHars;j++)
 {
  select->SetEtaCut(0., 1.1, j, 1) ;
  select->SetPtCut(0.1, 10. , j, 1);
 }
 select->SetConstrainCut(kTRUE) ;
 select->SetDcaGlobalCut(0.,0.1);
 //select->SetNhitsCut(3., 1) ;
 //select->SetPidCut("pi");

 // couts
 cout << " . Selection for R.P. calculation: " << endl ;
 select->PrintSelectionList() ;

 // flow weighter //

 AliFlowWeighter* fwgt = new AliFlowWeighter(select) ;

 // output file //

 TFile * wgtFile = new TFile(wgts.Data(),"RECREATE") ;
 fwgt->SetWgtFile(wgtFile) ;
 cout << " . Writing Histograms on  : " << fwgt->GetWgtFileName()  << "  . " << endl ;

 // init (make histograms, start the analysis)
 cout << endl ;
 kOne = fwgt->Init() ;  if(kOne) { cout << "   ok! " << endl ;}
 
 // flowEvents chain (imput) //

 TString input = "pFlowEvts.root" ;
 TFile* flowEventsFile = new TFile(input.Data(), "READ") ; // flowEventsFile->ls() ;
 Int_t nEvts = flowEventsFile->GetNkeys() ; 
 cout << " . Found  " << nEvts << " AliFlowEvents in file " << input.Data() << endl ;
 TList* flowEventsList = (TList*)flowEventsFile->GetListOfKeys() ; 

 // event loop
 cout << endl ;
 AliFlowEvent* flowEvt = new AliFlowEvent() ;
 for(Int_t ie=0;ie<nEvts;ie++)
 {
  TString evtName = flowEventsList->At(ie)->GetName() ;
  flowEventsFile->GetObject(evtName.Data(),flowEvt) ;
  // cout << "dumping event " << ie << " : " << endl ; flowEvt->Dump() ; cout << endl ; 
  Bool_t succ = fwgt->WeightEvent(flowEvt) ;  
  if(succ) { cout << ie << " done ... " << endl ; } 
 }
 
 // p.Id
 cout << endl ;
 cout << "Particles normalized abundance:" << endl;  // shows the bayesian vector
 for(Int_t k=0;k<2;k++)  { fwgt->PrintBayesian(k) ; }
 cout << endl ;

 // finish (save hists, close file)
 cout << endl ;
 kOne = fwgt->Finish() ; if(kOne) { cout << "   ok! " << endl ;}

 timer.Stop();
 cout << endl ;
 cout << " . " ; timer.Print();
 cout << " . here it was (weighter) ... " << endl ;  //juice!
 cout << endl ;

 // break ;
 
}
