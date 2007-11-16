/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                              *
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

/* $Id$ */
//
// An analysis task to check the ESD VZERO data in simulated data
// An analysis task to check the ESD VZERO data in simulated data
// An analysis task to check the ESD VZERO data in simulated data
// An analysis task to check the ESD VZERO data in simulated data
// An analysis task to check the ESD VZERO data in simulated data
//
//////////////////////////////////////////////////////////////////////////////

 
#include <TROOT.h>
#include <TChain.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TFile.h> 
#include <TString.h> 

#include "AliVZEROQATask.h" 
#include "AliESD.h" 
#include "AliESDVZERO.h"
#include "AliLog.h"

//______________________________________________________________________________
AliVZEROQATask::AliVZEROQATask(const char *name) : 
  AliAnalysisTask(name,""),  
  fChain(0),
  fESD(0),
  fOutputContainer(0), 
  fhVZERONbPMA(0),
  fhVZERONbPMC(0),
  fhVZEROMultA(0),
  fhVZEROMultC(0)
     
{
  // Constructor.
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0,  TObjArray::Class()) ; 
}

AliVZEROQATask::AliVZEROQATask(const AliVZEROQATask& ta) : 
AliAnalysisTask(ta.GetName(),""),  
fChain(ta.fChain),
fESD(ta.fESD),
fOutputContainer(ta.fOutputContainer), 
fhVZERONbPMA(ta.fhVZERONbPMA),
fhVZERONbPMC(ta.fhVZERONbPMC),
fhVZEROMultA(ta.fhVZEROMultA),
fhVZEROMultC(ta.fhVZEROMultC)

{
	// copy constructor
}

//_____________________________________________________________________________
AliVZEROQATask& AliVZEROQATask::operator = (const AliVZEROQATask& ap)
{
	// assignment operator
	
	this->~AliVZEROQATask();
	new(this) AliVZEROQATask(ap);
	return *this;
}

//______________________________________________________________________________
AliVZEROQATask::~AliVZEROQATask()
{
  // dtor
  
  fOutputContainer->Clear(); 
  delete fOutputContainer; 
  
  delete fhVZERONbPMA;
  delete fhVZERONbPMC; 
  delete fhVZEROMultA;
  delete fhVZEROMultC;
}

//______________________________________________________________________________
void AliVZEROQATask::ConnectInputData(const Option_t*)
{
  // Initialises branch container and histograms 
    
  AliInfo(Form("*** Initialization of %s", GetName())) ; 
  
  // Gets input data
  fChain = dynamic_cast<TChain *>(GetInputData(0)) ;
  if (!fChain) {
    AliError(Form("Input 0 for %s not found\n", GetName()));
    return ;
  }
  
  // One should first check if the branch address was taken by some other task
  char ** address = (char **)GetBranchAddress(0, "ESD");
  if (address) {
    fESD = (AliESD*)(*address);
  } else {
    fESD = new AliESD();
    SetBranchAddress(0, "ESD", &fESD);
  }
}

//________________________________________________________________________
void AliVZEROQATask::CreateOutputObjects()
{  
  // Creates histograms 
     
  OpenFile(0) ; 

  fhVZERONbPMA  = new TH1I("Nb of fired PMs in V0A", "VZERONbPMA" ,100 ,0 ,99);
  fhVZERONbPMC  = new TH1I("Nb of fired PMs in V0C", "VZERONbPMC" ,100 ,0 ,99);
  fhVZEROMultA  = new TH1I("Multiplicity in V0A", "VZEROMultA" ,50 ,0 ,49);
  fhVZEROMultC  = new TH1I("Multiplicity in V0C", "VZEROMultC" ,50 ,0 ,49);
  
  // Creates output container
  
  fOutputContainer = new TObjArray(4); 
  fOutputContainer->SetName(GetName()) ; 
  fOutputContainer->AddAt(fhVZERONbPMA, 0); 
  fOutputContainer->AddAt(fhVZERONbPMC, 1); 
  fOutputContainer->AddAt(fhVZEROMultA, 2); 
  fOutputContainer->AddAt(fhVZEROMultC, 3); 
   
}

//______________________________________________________________________________
void AliVZEROQATask::Exec(Option_t *) 
{
  // Processing of one event
  Long64_t entry = fChain->GetReadEntry() ;

  if (!fESD) {
    AliError("fESD is not connected to the input!") ; 
    return ; 
  }
  
  if ( !((entry-1)%100) ) 
    AliInfo(Form("%s ----> Processing event # %lld",  (dynamic_cast<TChain *>(fChain))->GetFile()->GetName(), entry)) ; 
 
  AliESDVZERO *esdVZERO=fESD->GetVZEROData();
   
  if (esdVZERO) { 
    fhVZERONbPMA->Fill(esdVZERO->GetNbPMV0A());
    fhVZERONbPMC->Fill(esdVZERO->GetNbPMV0C());  
    fhVZEROMultA->Fill(esdVZERO->GetMTotV0A());
    fhVZEROMultC->Fill(esdVZERO->GetMTotV0C());  
  }
  PostData(0, fOutputContainer);
  
}

//______________________________________________________________________________
void AliVZEROQATask::Terminate(Option_t *)
{
  // Processed when the event loop is ended
  
  fOutputContainer = (TObjArray*)GetOutputData(0);
  fhVZERONbPMA     = (TH1I*)fOutputContainer->At(0);
  fhVZERONbPMC     = (TH1I*)fOutputContainer->At(1);
  fhVZEROMultA     = (TH1I*)fOutputContainer->At(2);
  fhVZEROMultC     = (TH1I*)fOutputContainer->At(3);
  
  Bool_t problem = kFALSE ; 
  AliInfo(Form(" *** %s Report:", GetName())) ; 
  printf("        V0A Multiplicity Mean : %5.3f , RMS : %5.3f \n",fhVZEROMultA->GetMean(),fhVZEROMultA->GetRMS());
  printf("        V0C Multiplicity Mean : %5.3f , RMS : %5.3f \n",fhVZEROMultC->GetMean(),fhVZEROMultC->GetRMS());

  TCanvas  * c1 = new TCanvas("Number of PMs fired in V0A", "Number of PMs fired in V0A", 1);
  fhVZERONbPMA->SetAxisRange(0, 99);
  fhVZERONbPMA->SetLineColor(2);
  fhVZERONbPMA->Draw("SAME");
  c1->Update();
 
  TCanvas  * c2 = new TCanvas("Number of PMs fired in V0C", "Number of PMs fired in V0C", 1);
  fhVZERONbPMC->SetAxisRange(0,99);
  fhVZERONbPMC->SetLineColor(2);
  fhVZERONbPMC->Draw("SAME");
  c2->Update();

  TCanvas  * c3 = new TCanvas("Multiplicity in V0A", "Multiplicity in V0A", 1);
  fhVZEROMultA->SetAxisRange(0, 49);
  fhVZEROMultA->SetLineColor(2);
  fhVZEROMultA->Draw("SAME");
  c3->Update();
 
  TCanvas  * c4 = new TCanvas("Multiplicity in V0C", "Multiplicity in V0C", 1);
  fhVZEROMultC->SetAxisRange(0,49);
  fhVZEROMultC->SetLineColor(2);
  fhVZEROMultC->Draw("SAME");
  c4->Update();

  c1->Print("V0AMultiplicity.eps");
  c2->Print("V0CMultiplicity.eps");
  c3->Print("NumberV0APMs.eps");
  c4->Print("NumberV0CPMs.eps");
  
  char line[1024] ; 
  sprintf(line, ".!tar -zcf %s.tar.gz *.eps", GetName()) ; 
  gROOT->ProcessLine(line);
  sprintf(line, ".!rm -fR *.eps"); 
  gROOT->ProcessLine(line);
 
  AliInfo(Form("!!! All the eps files are in %s.tar.gz !!! ", GetName())) ;
  
  TString report ; 
  if(problem)
    report="Problems found, please check!!!";  
  else 
    report="OK";
  
  AliInfo(Form("*** %s Summary Report: %s\n",GetName(), report.Data())) ; 
  
}
