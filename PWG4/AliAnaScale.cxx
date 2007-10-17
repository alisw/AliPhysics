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

/* $Id: */

//_________________________________________________________________________
// A basic analysis task to analyse photon detected by PHOS
//
//*-- Yves Schutz 
//////////////////////////////////////////////////////////////////////////////

#include <TFile.h> 
#include <TH1.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TCanvas.h> 

#include "AliAnaScale.h" 
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "Riostream.h"

//______________________________________________________________________________
AliAnaScale::AliAnaScale() : 
  fDebug(0),
  fScale(1.0),
  fInputList(0x0), 
  fOutputList(0x0) 
{
  //Default constructor
}
//______________________________________________________________________________
AliAnaScale::AliAnaScale(const char *name) : 
  AliAnalysisTask(name,""),  
  fDebug(0),
  fScale(1.0), 
  fInputList(0x0), 
  fOutputList(0x0) 
{
  // Constructor.
  // Called only after the event loop
  SetPostEventLoop(kTRUE);
  // Input slot #0 
  DefineInput(0,  TList::Class()) ; 
  // Output slot 
  DefineOutput(0,  TList::Class()) ; 
}

//______________________________________________________________________________
AliAnaScale::~AliAnaScale()
{
  // dtor
  
}


//______________________________________________________________________________
void AliAnaScale::ConnectInputData(const Option_t*)
{
  // Initialisation of branch container and histograms 
    
  AliInfo(Form("*** Initialization of %s", GetName())) ; 
  fInputList     = dynamic_cast<TList*>(GetInputData(0)) ;  
}
//________________________________________________________________________
void AliAnaScale::CreateOutputObjects()
{  
  // Create the outputs containers
  // Is created in Exec(), because the input must be available

}

//______________________________________________________________________________
void AliAnaScale::Exec(Option_t *) 
{
  // Do the Scaling
    
  fOutputList = new TList() ; 
  fOutputList->SetName(GetName()) ; 
  TIter next(fInputList) ; 	
  TObject * h ; 
  while ( (h = next()) ) { 
    if(h){
      if ( strcmp(h->ClassName(),"TNtuple") ) {
      char name[20] ; 
      sprintf(name, "%sScaled", h->GetName()) ; 
      TH1 * hout = dynamic_cast<TH1*> (h->Clone(name)) ; 
      hout->Scale(fScale) ;  
      fOutputList->Add(hout) ; 
      } 
      else  fOutputList->Add(h) ; 
    }
  }
  cout<<"end"<<endl;
  PostData(0, fOutputList);
}


//______________________________________________________________________________
void AliAnaScale::Init()
{
  // Intialisation of parameters
  AliInfo("Doing initialisation") ;
  // nothing to be done
}

//______________________________________________________________________________
void AliAnaScale::Terminate(Option_t *)
{
  // Processing when the event loop is ended
  

}
