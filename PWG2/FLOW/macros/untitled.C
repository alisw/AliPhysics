// exampleMacro.C

#include <vector>
#include <iostream>
#include <fstream>
#include "TVector.h"
#include "TVector2.h"
#include "TMath.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TStopwatch.h"
#include "TList.h"

//#include "AliFlowConstants.h"

using namespace std; //required for resolving the 'cout' symbol

bool untitled(int something=0)
{
 cout << " . Example of how to use the new flow package (2007) ...  " << something << endl ;
 cout << endl ;

 TStopwatch timer;
 timer.Start();

 gSystem->Load("libPhysics.so");
 gSystem->Load("AliFlow_All.so");

 // allocate space for event quantities //
 
  Float_t  fQ[2][2];		// flow vector
  Float_t  fPsi[2][2];		// event plane angle
  Int_t    fMult[2][2];     	// multiplicity
  Float_t  fQnorm[2][2];    	// Q/Sqrt(Mult)

 // output file //

 // ...

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
 select->SetNhitsCut(1., 1) ;

 // couts
 cout << " . Selection for R.P. calculation: " << endl ;
 select->PrintSelectionList() ;
 cout << " . Selection for correlation analysis: " << endl ;
 select->PrintList() ;
 cout << " . Selection for V0 analysis: " << endl ;
 select->PrintV0List() ;

 // Wgt file // 

 // ...

 // flowEvents chain (imput) //

 // organise the loop     ... in a very clever way :)
 TString input = "fFlowEvts.root" ;
 TFile* flowEventsFile = new TFile(input.Data(), "READ") ; // flowEventsFile->ls() ;
 Int_t nEvts = flowEventsFile->GetNkeys() ; 
 cout << " . Found  " << nEvts << " AliFlowEvents in file " << input.Data() << endl ;
 TList* flowEventsList = (TList*)flowEventsFile->GetListOfKeys() ; 
 cout << endl ;

 // event loop
 AliFlowEvent* flowEvt = new AliFlowEvent() ;
 for(Int_t ie=0;ie<nEvts;ie++)
 {
  TString evtName = flowEventsList->At(ie)->GetName() ;
  flowEventsFile->GetObject(evtName.Data(),flowEvt) ;
  //flowEvt->Dump() ;
  cout << "* event ID = " << flowEvt->EventID() << " :  found " << flowEvt->FlowEventMult() << " AliFlowTracks, and " << flowEvt->V0Mult() << " AliFlowV0s . " << endl ;  

  if(select->Select(flowEvt))		// event selected
  {
   cout << "  selected ... " <<  endl ; 
   
   flowEvt->SetSelections(select) ;	// applies the selection of tracks for RP calculation
   flowEvt->MakeAll() ; 		// calculates Q, Psi, Mult for RP analysis [nSel]*[nHar] times

   int selCheck = 0 ;
   for (int k = 0; k < Flow::nSels; k++) 
   {
    select->SetSelection(k) ;		// set the Sel. of interest
    //cout << "  selN : " << k << " = " << select->Sel() <<  endl ; 
    for (int j = 0; j < Flow::nHars; j++) 
    {
     select->SetHarmonic(j) ;		// set the Har. of interest 
     //cout << "  harN : " << j << " = " << select->Har() << endl ; 
     select->SetSubevent(-1) ;		// set the SubEvt of interest (-1 for the full event) 
     //cout << "  subN : " << -1 << " = " << select->Sub() << endl ; 

     cout << "  calculating  ...  " << k << "," << j <<  endl ; 

     fQ[k][j]	  = flowEvt->Q(select) ;	   //cout << "  Q ..."      << fQ[k][j]	  << endl ; 
     fPsi[k][j]   = flowEvt->Psi(select) ;	   //cout << "  Psi  ..."   << fPsi[k][j]   << endl ; 
     fQnorm[k][j] = flowEvt->NormQ(select).Mod() ; //cout << "  NormQ  ..." << fQnorm[k][j] << endl ; 
     fMult[k][j]  = flowEvt->Mult(select) ;	   //cout << "  Mult  ..."  << fMult[k][j]  << endl ; 

     selCheck += fMult[k][j] ; 
    }
   }
  }
  if(selCheck==0) { cout << " Bad Event : " << ie << endl ; }
  else 		  { cout << " Event : " << ie << " :  psi[1][1] = " << fPsi[1][1] << " (" << fMult[1][1] << " traks)" << endl ; }
 }
 
 return kTRUE ;
}
