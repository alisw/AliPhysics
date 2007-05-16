/////////////////////////////////////////////////////////////
//
// $Id$
//
// Author: Emanuele Simili
//
/////////////////////////////////////////////////////////////
//
// Description: ROOT macro to perform the Flow analysis on AliFlowEvents (new way)
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
#include "TKey.h"
#include "TList.h"

using namespace std; //required for resolving the 'cout' symbol

int testAnal(int cen = -1)
{
 cout << " . Here the new flow analysis (2007 rs) ... " << endl ;
 cout << endl ;

 bool kOne = kFALSE ;
 int limit = -1 ;

 if(limit>0)
 {
  cout << " . limited to " << limit << "events . " << endl ;
  cout << endl ;
 }
 
 // histogram file (output) //

 TString output = "flowEvtsAnal" ; 
 if(cen >= 0) { output += cen ; }
 output += ".root" ;
 
 // flowEvents chain (imput) //

 TString fFlowFileName = "flowEvts" ;
 if(cen >= 0) { fFlowFileName += cen ; }
 fFlowFileName += ".root" ;

 // Wgt file (optional) //

 TString wgtFileName = "flowPhiWgt" ;
 if(cen >= 0) { wgtFileName += cen ; }
 wgtFileName += ".root" ;

 // start, load libs //

 TStopwatch timer;
 timer.Start();

 //gSystem->Load("libPhysics.so");
 gSystem->Load("libPWG2flow.so");

 // open output file //

 TFile * fFlowfile = new TFile(output.Data(),"RECREATE") ;
 //fFlowfile->cd() ; 

 // selection object (cuts) // 

 AliFlowSelection* select = new AliFlowSelection() ;
 
 // Event Cuts
 select->SetCentralityCut(-1) ;  // cen
 select->SetRunIdCut(-1) ;

 // R.P. calculation cuts
 for(int j=0;j<AliFlowConstants::kHars;j++)
 {
  select->SetEtaCut(0., 0.9, j, 1) ;
  select->SetPtCut(0.1, 10. , j, 1);
 }
 select->SetConstrainCut(kTRUE) ;
 select->SetDcaGlobalCut(0.,0.1);
 //select->SetNhitsCut(0, 1) ;
 //select->SetPidCut("pi");

 // Tracks Correlation analysis cuts
 select->SetEtaPart(-0.9,0.9);
 select->SetPtPart(0.1,10.);
 select->SetConstrainablePart(kTRUE);
 select->SetDcaGlobalPart(0.,0.1);
 //select->SetPPart(0.1,10.);
 //select->SetEtaAbsPart(0.5,1.1);
 //select->SetPidPart("pi");
 //select->SetYPart(..,..);
 //select->SetPidProbPart(0.5,1.);
 //select->SetFitPtsPart(3, 480);
 //select->SetFitOverMaxPtsPart(0.1,3.);
 //select->SetDedxPtsPart(1.,0.);
 //select->SetChiSqPart(0.,100.);

 // V0 Correlation analysis cuts
 //select->SetV0DcaCross(0.,0.1) ;
 select->SetV0Mass(0.48,0.52) ;	    // Mk0 = 0.49765
 select->SetV0SideBands(0.10) ;
 //select->SetV0P(0.,10.) ;
 select->SetV0Pt(0.1,10.) ;
 select->SetV0Eta(-2.1,2.1) ;
 //select->SetV0EtaAbs(0.,0.9) ;
 //select->SetV0Pid("0") ;
 //select->SetV0Y(0.,10.) ;
 //select->SetV0Lenght(0.,100.) ;
 //select->SetV0LenghtOverSigma(0.,5.) ;
 //select->SetV0ChiSqPart(0.,100.) ;

 // couts
 cout << " . Selection for R.P. calculation: " << endl ;
 select->PrintSelectionList() ;
 cout << " . Selection for correlation analysis: " << endl ;
 select->PrintList() ;
 cout << " . Selection for V0 analysis: " << endl ;
 select->PrintV0List() ;

 // flow analyser //

 AliFlowAnalyser* flow = new AliFlowAnalyser(select) ;

 // output file
 flow->SetHistFileName(output.Data()) ;
 cout << " . Writing Histograms on  : " << flow->GetHistFileName()   << "  . " << endl ;

 // Wgt s
 TFile* wgtFile = new TFile(wgtFileName.Data(),"READ");
 if(!wgtFile || wgtFile->IsZombie()) { cout << " . NO phi Weights . " << endl ;}
 else { flow->FillWgtArrays(wgtFile) ; cout << " . Weights from  : " << flow->GetWgtFileName() << "  . " << endl ; }

 // Analysis settings
 flow->SetFlowForV0(kFALSE) ;   // default kTRUE.
 flow->SetSub(1) ; //eta // flow->SetSub(0) ; //rnd // flow->SetSub(-1) ; //charged
 //flow->SetV1Ep1Ep2() ;        // default kFALSE.
 //flow->SetRedoWgt();      	// default kFALSE. recalculates phiWgt (even if phiWgt file is already there)
 flow->SetUsePhiWgt(kFALSE) ;   // default kTRUE if phiWgt file is there, kFALSE if not (& phiWgt file is created)
 flow->SetUseOnePhiWgt() ; // or // flow->SetUseFirstLastPhiWgt() ; // uses 1 or 3 wgt histograms (default is 1)
 flow->SetUseBayWgt(kFALSE) ;   // default kFALSE. uses bayesian weights in P.id.
 //flow->SetUsePtWgt();		// default kFALSE. uses pT as a weight for RP determination
 //flow->SetUseEtaWgt();	// default kFALSE. uses eta as a weight for RP determination
 //flow->SetCustomRespFunc();	// default kFALSE. uses the combined response function from the ESD

 // init (make histograms, start the analysis)
 cout << endl ;
 kOne = flow->Init() ;  if(kOne) { cout << "   ok! " << endl ;}

 // organizing event loop
 TFile* flowEventsFile = new TFile(fFlowFileName.Data(), "READ") ; // flowEventsFile->ls() ;
 Int_t nEvts = flowEventsFile->GetNkeys() ; 
 cout << " . Found  " << nEvts << " AliFlowEvents in file " << fFlowFileName.Data() << endl ;
 TList* flowEventsList = (TList*)flowEventsFile->GetListOfKeys() ; 
 AliFlowEvent* flowEvt = new AliFlowEvent() ;
 TIter next(flowEventsList); 
 TKey* key ;

 cout << endl ;

 // Loop over the events
 Int_t count = 0 ;
 while( key=(TKey *)next() ) 
 {
  // TString evtName = flowEventsList->At(ie)->GetName() ;
  // flowEventsFile->GetObject(evtName.Data(),flowEvt) ;
  // cout << "dumping event " << ie << " : " << endl ; flowEvt->Dump() ; cout << endl ; 

  flowEvt = (AliFlowEvent *)key->ReadObj();
  if(!flowEvt) break;

  // Process event .......
  Bool_t succ = flow->Analyze(flowEvt) ;  
  flow->PrintEventQuantities() ;
  if(succ) { cout << count << " done ... " << endl ; }  
  
  delete flowEvt ;

  if(count == limit) { break ; }  
  count++ ;
 }

 // p.Id
 cout << endl ;
 cout << "Particles normalized abundance:" << endl;  // shows the bayesian vector
 flow->PrintRunBayesian() ; cout << endl;
 
 // resolution
 cout << endl ;
 flow->Resolution() ;

 // saves the wgt file
 cout << endl ;
 flow->Weightening() ;

 // finish (save hists, close file)
 cout << endl ;
 kOne = flow->Finish() ; if(kOne) { cout << "   ok! " << endl ;}

 timer.Stop();
 cout << endl ;
 cout << " . " ; timer.Print();
 cout << " . here it was (analyser rs) ... " << endl ;  //juice!
 cout << endl ;
 
 // cout << endl ; cout << " Memory Check (from Paul)" << endl ; 
 // gObjectTable->Print();
 // cout << endl ; cout << endl ;

 return cen ;
}
