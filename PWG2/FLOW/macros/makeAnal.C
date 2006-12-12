/////////////////////////////////////////////////////////////
//
// $Id$
//
// Author: Emanuele Simili
//
/////////////////////////////////////////////////////////////
//
// Description:  ROOT macro to perform the AliFlowAnalysis 
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

void makeAnal(TString wch = "ESD")
{
 cout << " . Here the newly compiled flow analysis . " << endl ;
 cout << " . AliFlowEvent(s) --> root histograms  " << endl ;
 cout << endl ;

 gSystem->Load("libPhysics.so");
 gSystem->Load("libPWG2.so");
 gSystem->Load("libPWG2flow.so");

 //char *loadFrom = gSystem->ExpandPathName("$ALICE_ROOT/PWG2/FLOW/AliFlow_Anal.so") ;
 //cout << " . *** Loading libs from   ' " << loadFrom << " ' ." << endl ;
 //gSystem->Load(loadFrom);
 cout << " .                          ...libs  loaded" << endl ;

 int cc = -1 ;  	  // centrality class
 TString ver = "h_" ;     // names ... according to the FlowMaker 
 int selNumber = 0 ;      // ... 

 TStopwatch timer;
 timer.Start();

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

// Correlation analysis cuts
 select->SetEtaPart(-1.1,1.1);
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

// V0 analysis cuts
 //select->SetV0DcaCross(0.,0.1) ;
 select->SetV0Mass(0.4875,0.5078) ;	// Mk0 = 0.49765
 select->SetV0SideBands(0.08) ;
 //select->SetV0P(0.,10.) ;
 select->SetV0Pt(0.1,10.) ;
 select->SetV0Eta(-2.1,2.1) ;
 //select->SetV0EtaAbs(0.,0.9) ;
 //select->SetV0Pid("0") ;
 //select->SetV0Y(0.,10.) ;
 //select->SetV0Lenght(0.,100.) ;
 //select->SetV0LenghtOverSigma(0.,5.) ;
 //select->SetV0ChiSqPart(0.,100.) ;

 cout << " . Selection for R.P. calculation: " << endl ;
 select->PrintSelectionList() ;
 cout << " . Selection for correlation analysis: " << endl ;
 select->PrintList() ;
 cout << " . Selection for V0 analysis: " << endl ;
 select->PrintV0List() ;

 AliFlowAnalysisMaker* flow = new AliFlowAnalysisMaker("FlowAnal",*select) ;

// Input/Output/Wgt files : auto-naming
 TString inpt = ver ; inpt += wch ; inpt += "flow" ;
 TString outp = inpt ; outp += "AnalPlot" ;
 if(selNumber) { outp += selNumber ; }
 if(cc>=0) { outp += "_" ; outp += cc ; }
 TString wgt  = "flowPhiWgt." ; wgt += wch ;
 inpt += ".root" ; outp += ".root" ; wgt += ".root" ;

// tmporary setting ---
 inpt = "flowEvtS.root" ;
 outp = "flowEvtSanalPlot.root" ;
// tmporary setting ---

//OUTput file
 flow->SetHistFileName(outp.Data());       // __MC/ESDflowAnalPlot.root
 cout << " . Writing Histograms on  : " << flow->GetHistFileName()   << "  . " << endl ;

// auto file name
 flow->SetOneInputFile(kTRUE) ;
 flow->SetInputFileName(inpt.Data()) ;    // __MC/ESDflow.root
 cout << " . Reading FlowEvents from  : " << flow->GetInputFileName()  << "  . " << endl ;
//INput: one flowEvent file
 //flow->SetOneInputFile(kTRUE) ;
 //flow->SetInputFileName("neo_ESDflow.root") ;
//INput: more than one file
//  //flow->SetInputFileNames("...") ;
//  cout << " . Reading FlowEvents from  : " << "  tot files  " << endl ;

//WGT file
 flow->SetPhiWgtFileName(wgt.Data());      // flowPhiWgt.MC/ESD.root
 cout << "# Weights File (if there)  : " << flow->GetPhiWgtFileName() << "  . " << endl ;
 cout << endl ;

// Analysis settings
 flow->SetFlowForV0() ;         // default kTRUE.
 // flow->SetEtaSub() ;          // default is Disabled - DO NOT USE (kFALSE) AS AN ARGUMENT !!!
 //flow->SetV1Ep1Ep2() ;        // default kFALSE.
 flow->SetShuffle() ;           // default kFALSE. shuffles track array
 //flow->SetRedoWgt();      	// default kFALSE. recalculates phiWgt (even if phiWgt file is already there)
 flow->SetUsePhiWgt() ;         // default kTRUE if phiWgt file is there, kFALSE if not (& phiWgt file is created)
 //flow->SetUseBayWgt() ;   	// default kFALSE. uses bayesian weights in P.id.
 //flow->SetUsePtWgt();		// default kFALSE. uses pT as a weight for RP determination
 //flow->SetUseEtaWgt();		// default kFALSE. uses eta as a weight for RP determination
 flow->SetUseOnePhiWgt() ; // or // flow->SetUseFirstLastPhiWgt() ; // uses 1 or 3 wgt histograms (default is 1)
 //
 flow->SetMakeAll() ;           // default kTRUE. this should speed up the execution time
 flow->SetMaxLabel(1000) ;	// THERE WAS A CRASH WHEN LOOPING OVER MORE FLOWEVENT FILES !!!
 flow->SetFillLabels(kFALSE) ;   // ... ~A CHANGE WAS MADE IN THE ANALYSIS TO FIX IT (// *temp*)
 //
 flow->SetDebugg(1) ;		// more couts
 //Flow::mDebug = kTRUE ;	// more couts from singleton

// init (make histograms, start the analysis)
 Int_t step = 0 ;
 step += flow->Init() ;
 if(step<1) { cout << "#! UNEXPECTED TERMINATION !# " << step << endl ; break ; }

// make (loop)
 step += flow->Make() ;
 if(step<2) { cout << "#! UNEXPECTED TERMINATION !# " << step << endl ; break ; }
 cout << endl;
 cout << "Particles normalized abundance:" << endl;  // shows the bayesian vector
 flow->PrintRunBayesian() ; cout << endl;

// finish (save hists, close files)
 step += flow->Finish() ;
 if(step<3) { cout << "#! UNEXPECTED TERMINATION !# " << step << endl ; break ; }

 timer.Stop();
 cout << endl ;
 cout << " . " ; timer.Print();
 cout << " . here it was (analysis) ... " << endl ;  //juice!
 cout << endl ;
}
