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

//______________________________________________________________________________
AliAnaScale::AliAnaScale() : 
  fDebug(0),
  fScale(1.0),
  fInputList(0x0), 
  fOutputList(0x0), 
  fhInPHOSEnergy(0),
  fhOuPHOSEnergy(0)
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
  fhInPHOSEnergy(0),
  fhOuPHOSEnergy(0)
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
  fhInPHOSEnergy = dynamic_cast<TH1D*>(fInputList->At(2));
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
    
  fhOuPHOSEnergy = static_cast<TH1D*>(fhInPHOSEnergy->Clone("PHOSEnergyScaled")) ;
  
  // create output container
  
  fOutputList = new TList() ; 
  fOutputList->SetName(GetName()) ; 

  fOutputList->AddAt(fhOuPHOSEnergy,          0) ; 
  fhOuPHOSEnergy->Scale(fScale) ; 
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
  
  AliInfo(Form(" *** %s Report:", GetName())) ; 
  printf("        PHOS Energy Integral In         : %5.3e \n", fhInPHOSEnergy->Integral() ) ;
  printf("        PHOS Energy Integral Ou         : %5.3e \n", fhOuPHOSEnergy->Integral() ) ;

  TCanvas  * cPHOS = new TCanvas("cPHOS", "PHOS ESD Test", 400, 10, 600, 700) ;

  cPHOS->Divide(2, 1);
  cPHOS->cd(1) ; 
  if ( fhInPHOSEnergy->GetMaximum() > 0. ) 
    gPad->SetLogy();
  fhInPHOSEnergy->SetAxisRange(0, 25.);
  fhInPHOSEnergy->SetLineColor(2);
  fhInPHOSEnergy->Draw();

  cPHOS->cd(2) ; 
  if ( fhOuPHOSEnergy->GetMaximum() > 0. ) 
    gPad->SetLogy();
  fhOuPHOSEnergy->SetAxisRange(0, 25.);
  fhOuPHOSEnergy->SetLineColor(2);
  fhOuPHOSEnergy->Draw();
}
