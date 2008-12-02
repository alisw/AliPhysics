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

//______________________________________________________________________________
AliAnaScale::AliAnaScale() : 
  fDebug(0),
  fScale(1.0),
  fInputList(0x0), 
  fOutputList(0x0),
  fSumw2(0),
  fhCount() 
{
  //Default constructor
}
//______________________________________________________________________________
AliAnaScale::AliAnaScale(const char *name) : 
  AliAnalysisTask(name,""),  
  fDebug(0),
  fScale(1.0), 
  fInputList(0x0), 
  fOutputList(0x0), 
  fSumw2(0),
  fhCount(0) 
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
    
  if(fDebug > 1) printf("*** Initialization of %s \n", GetName()) ; 
  fInputList     = dynamic_cast<TList*>(GetInputData(0)) ;  
}
//________________________________________________________________________
void AliAnaScale::CreateOutputObjects()
{  
  // Create the outputs containers

  fOutputList = new TList() ; 
  fOutputList->SetName(GetName()) ; 

  fhCount =new TH1F("hCount","count files",1,0,1);  
  fOutputList->Add(fhCount);

}

//______________________________________________________________________________
void AliAnaScale::Exec(Option_t *) 
{
  // Do the Scaling

  if(fDebug > 0 ) printf(">>>>> Scaling factor %e, do Sumw2 %d <<<<< \n",fScale,fSumw2) ;

  TIter next(fInputList) ; 	
  TObject * h ; 
  while ( (h = next()) ) { 
    if(h){
      if ( strcmp(h->ClassName(),"TNtuple") ) {
      char name[128] ; 
      sprintf(name, "%sScaled", h->GetName()) ; 
      TH1 * hout = dynamic_cast<TH1*> (h->Clone(name)) ; 
     
      if(fSumw2) hout->Sumw2();
      hout->Scale(fScale) ;  
      fOutputList->Add(hout) ; 
      } 
      else  fOutputList->Add(h) ; 
    }
  }
  // number of files

  //File scaled, needed for file merging on grid
  fhCount->Fill(0);
 
  PostData(0, fOutputList);
}


//______________________________________________________________________________
void AliAnaScale::Init()
{
  // Intialisation of parameters
  if(fDebug > 0 )printf("No initialization in scale class \n") ;

}

//______________________________________________________________________________
//void AliAnaScale::Terminate(Option_t *)
//{
//  // Processing when the event loop is ended
//  
//
//}
