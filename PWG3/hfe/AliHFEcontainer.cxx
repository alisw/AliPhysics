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
//
// HFE correction framework container
// Contains many single containers
// Extra fuctionality like appending added
//
// Author:
//   Markus Fasel <M.Fasel@gsi.de>
//
#include <iostream>
#include <TAxis.h>
#include <TClass.h>
#include <TCollection.h>
#include <THashList.h>
#include <THnSparse.h>
#include <TList.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>

#include "AliCFContainer.h"
#include "AliHFEcontainer.h"
#include "AliHFEtools.h"

ClassImp(AliHFEcontainer)
ClassImp(AliHFEcontainer::AliHFEvarInfo)

//__________________________________________________________________
AliHFEcontainer::AliHFEcontainer():
  TNamed("HFEcontainer", ""),
  fContainers(NULL),
  fCorrelationMatrices(NULL),
  fVariables(NULL),
  fNVars(0),
  fNEvents(0)
{
  //
  // Default constructor
  //
}

//__________________________________________________________________
AliHFEcontainer::AliHFEcontainer(const Char_t *name):
  TNamed(name, ""),
  fContainers(NULL),
  fCorrelationMatrices(NULL),
  fVariables(NULL),
  fNVars(0),
  fNEvents(0)
{
  //
  // Default constructor
  //
  fContainers = new THashList();
  fContainers->SetOwner();
}

//__________________________________________________________________
AliHFEcontainer::AliHFEcontainer(const Char_t *name, UInt_t nVar):
  TNamed(name, ""),
  fContainers(NULL),
  fCorrelationMatrices(NULL),
  fVariables(NULL),
  fNVars(0),
  fNEvents(0)
{
  //
  // Constructor
  // Setting Number of Variables too
  //
  fContainers = new THashList();
  fContainers->SetOwner();
  SetNumberOfVariables(nVar);
}

//__________________________________________________________________
AliHFEcontainer::AliHFEcontainer(const AliHFEcontainer &ref):
  TNamed(ref),
  fContainers(NULL),
  fCorrelationMatrices(NULL),
  fVariables(NULL),
  fNVars(ref.fNVars),
  fNEvents(ref.fNEvents)
{
  //
  // Copy constructor
  // creates a new object with new (empty) containers
  //
  if(fNVars){
    fVariables = new TObjArray(fNVars);
    AliHFEvarInfo *vtmp = NULL;
    for(UInt_t ivar = 0; ivar < fNVars; ivar++){
      vtmp = static_cast<AliHFEvarInfo *>(ref.fVariables->UncheckedAt(ivar));
      fVariables->AddAt(new AliHFEvarInfo(*vtmp), ivar);
    }
  }
  fContainers = new THashList;
  fContainers->SetOwner();
  AliCFContainer *ctmp = NULL;
  for(Int_t ien = 0; ien < ref.fContainers->GetEntries(); ien++){
    ctmp = static_cast<AliCFContainer *>(ref.fContainers->At(ien));
    CreateContainer(ctmp->GetName(), ctmp->GetTitle(), ctmp->GetNStep());
  }
  // Copy also correlation matrices
  if(ref.fCorrelationMatrices){
    THnSparseF *htmp = NULL;
    fCorrelationMatrices = new THashList;
    fCorrelationMatrices->SetOwner();
    for(Int_t ien = 0; ien < ref.fCorrelationMatrices->GetEntries(); ien++){
      htmp = static_cast<THnSparseF *>(ref.fCorrelationMatrices->At(ien));
      CreateCorrelationMatrix(htmp->GetName(), htmp->GetTitle());
    }
  }
}

//__________________________________________________________________
AliHFEcontainer &AliHFEcontainer::operator=(const AliHFEcontainer &ref){
  //
  // Assignment operator
  // Cleanup old object, create a new one with new containers inside
  //
  this->~AliHFEcontainer(); // cleanup old object before creating the new onwe
  TNamed::operator=(ref);
  fContainers = new THashList();
  fCorrelationMatrices = NULL;
  fNVars = ref.fNVars;
  if(fNVars){
    fVariables = new TObjArray(fNVars);
    AliHFEvarInfo *vtmp = NULL;
    for(UInt_t ivar = 0; ivar < fNVars; ivar++){
      vtmp = static_cast<AliHFEvarInfo *>(ref.fVariables->UncheckedAt(ivar));
      fVariables->AddAt(new AliHFEvarInfo(*vtmp), ivar);
    }
  } else {
    // No varible defined, do not try to copy anything
    fVariables = NULL;
    return *this;
  }

  // Reference contains content, try copying also the containers and the correlation matrices
  fContainers = new THashList();
  AliCFContainer *ctmp = NULL;
  for(Int_t ien = 0; ien < ref.fContainers->GetEntries(); ien++){
    ctmp = static_cast<AliCFContainer *>(ref.fContainers->At(ien));
    fContainers->Add(new AliCFContainer(*ctmp));
  }
  // Copy also correlation matrices
  if(ref.fCorrelationMatrices){
    THnSparseF *htmp = NULL;
    fCorrelationMatrices = new THashList;
    fCorrelationMatrices->SetOwner();
    for(Int_t ien = 0; ien < ref.fCorrelationMatrices->GetEntries(); ien++){
      htmp = static_cast<THnSparseF *>(ref.fCorrelationMatrices->At(ien));
      CreateCorrelationMatrix(htmp->GetName(), htmp->GetTitle());
    }
  }
  return *this;
}

//__________________________________________________________________
AliHFEcontainer::~AliHFEcontainer(){
  //
  // Destructor
  //
  delete fContainers;
  if(fCorrelationMatrices) delete fCorrelationMatrices;
  if(fVariables){
    fVariables->Delete();
    delete fVariables;
  }
}

//__________________________________________________________________
Long64_t AliHFEcontainer::Merge(TCollection *coll){
  //
  // Merge Container
  //
  if(!coll)
    return 0;
  if(coll->IsEmpty())
    return 1;

  TIter iter(coll);
  TObject *o = NULL;
  Long64_t count = 0;
  while((o = iter())){
    AliHFEcontainer *cont = dynamic_cast<AliHFEcontainer *>(o);
    if(!cont) continue;

    // Merge the two TObjArrays
    TList containers;
    containers.Add(cont->fContainers);
    fContainers->Merge(&containers);

    if(fCorrelationMatrices && cont->fCorrelationMatrices){
      containers.Clear();
      containers.Add(cont->fCorrelationMatrices);
      fCorrelationMatrices->Merge(&containers);
    }

    fNEvents += cont->GetNumberOfEvents();
    count++;
  }
  return count + 1;
}

//__________________________________________________________________
void AliHFEcontainer::SetNumberOfVariables(UInt_t nVar){
  //
  // Define the number of variables 
  // Initialize containers for the variable informations
  //
  if(fNVars) return;

  fNVars = nVar;
  fVariables = new TObjArray(nVar);
  for(UInt_t ivar = 0; ivar < nVar; ivar++)
    fVariables->AddAt(new AliHFEvarInfo, ivar);
}

//__________________________________________________________________
void AliHFEcontainer::CreateContainer(const Char_t *name, const Char_t *title, UInt_t nStep){
  //
  // Create a new Correction Framework Container and store it 
  //
  if(fContainers->FindObject(name)){
    AliError(Form("Container %s already exists. Cannot replace it!", name));
    return;
  }
  
  Int_t *nBins = new Int_t[fNVars];
  AliHFEvarInfo *var = NULL;
  for(UInt_t ivar = 0; ivar < fNVars; ivar++){ 
    var = dynamic_cast<AliHFEvarInfo *>(fVariables->UncheckedAt(ivar));
    nBins[ivar] = var ? var->GetNumberOfBins() : 0;
  }
  AliCFContainer *cont = new AliCFContainer(name, title, nStep, fNVars, nBins);
  for(UInt_t ivar = 0; ivar < fNVars; ivar++){
    var = dynamic_cast<AliHFEvarInfo *>(fVariables->UncheckedAt(ivar));
    if(var){
      cont->SetBinLimits(ivar, var->GetBinning());
      cont->SetVarTitle(ivar, var->GetVarName()->Data());
    }
  }
  delete[] nBins;
  fContainers->Add(cont);
  AliInfo(Form("Container %s created with %d cut steps", name, nStep));
}

//__________________________________________________________________
void AliHFEcontainer::CreateCorrelationMatrix(const Char_t *name, const Char_t *title){
  //
  // Create Correlation Matrix
  //
  if(!fCorrelationMatrices){
    fCorrelationMatrices = new THashList;
    fCorrelationMatrices->SetName("fCorrelationMatrices");
    fCorrelationMatrices->SetOwner();
  }

  Int_t *nBins = new Int_t[2*fNVars];
  AliHFEvarInfo *var = NULL;
  for(UInt_t ivar = 0; ivar < fNVars; ivar++){
    var = dynamic_cast<AliHFEvarInfo *>(fVariables->UncheckedAt(ivar));
    if(var){
      nBins[ivar] = var->GetNumberOfBins();
      nBins[ivar+fNVars] = var->GetNumberOfBins();
    }
  }

  THnSparseF * hTmp = new THnSparseF(name, title, 2*fNVars, nBins);
  for(UInt_t ivar = 0; ivar < fNVars; ivar++){
    var = dynamic_cast<AliHFEvarInfo *>(fVariables->UncheckedAt(ivar));
    if(var){
      hTmp->SetBinEdges(ivar,var->GetBinning());
      //hTmp->GetAxis(ivar)->Set(var->GetNumberOfBins(), var->GetBinning());
      hTmp->GetAxis(ivar)->SetTitle(var->GetVarName()->Data());
      //hTmp->GetAxis(ivar + fNVars)->Set(var->GetNumberOfBins(), var->GetBinning());
      hTmp->GetAxis(ivar + fNVars)->SetTitle(Form("%s_{MC}", var->GetVarName()->Data()));
      hTmp->SetBinEdges(ivar+fNVars,var->GetBinning());
    }
  }
  hTmp->Sumw2();
  fCorrelationMatrices->AddLast(hTmp);
}

//__________________________________________________________________
AliCFContainer *AliHFEcontainer::GetCFContainer(const Char_t *name) const{
  //
  // Find a given container 
  //
  return dynamic_cast<AliCFContainer *>(fContainers->FindObject(name));
}

//__________________________________________________________________
THnSparseF *AliHFEcontainer::GetCorrelationMatrix(const Char_t *name) const{
  //
  // Find Correlation Matrix
  //
  if(fCorrelationMatrices) return dynamic_cast<THnSparseF *>(fCorrelationMatrices->FindObject(name));
  else return 0x0;

}

//__________________________________________________________________
void AliHFEcontainer::FillCFContainer(const Char_t *name, UInt_t step, const Double_t * const content, Double_t weight) const {
  //
  // Fill container
  //
  AliCFContainer *cont = GetCFContainer(name);
  if(!cont) return;
  cont->Fill(content, step, weight);
}

//__________________________________________________________________
void AliHFEcontainer::FillCFContainerStepname(const Char_t *name, const Char_t *steptitle, const Double_t * const content, Double_t weight)const{
  //
  // Fill container
  //
  AliCFContainer *cont = GetCFContainer(name);
  if(!cont) return;
  // find the matching step title
  Int_t mystep = -1;
  for(Int_t istep = 0; istep < cont->GetNStep(); istep++){
    TString tstept = cont->GetStepTitle(istep);
    if(!tstept.CompareTo(steptitle)){
      mystep = istep;
      break;
    }
  }
  if(mystep < 0){
    // step not found
    AliDebug(1, Form("Step %s not found in container %s", steptitle, name));
    return;
  }
  AliDebug(1, Form("Filling step %s(%d) for container %s", steptitle, mystep, name));
  cont->Fill(content, mystep, weight);
}

//__________________________________________________________________
AliCFContainer *AliHFEcontainer::MakeMergedCFContainer(const Char_t *name, const Char_t *title, const Char_t* contnames) const {
  //
  // Merge CF Container out of several containers 
  // Container names are separated by :
  // returns a new object which has to be taken care of by the user
  //

  TObjArray *containers = TString(contnames).Tokenize(":");
  // we first need the size of the container to be merged
  Int_t nStepMerged = 0;
  AliCFContainer *ctemp = NULL;
  TObjString *cname = NULL;
  for(Int_t icont = 0; icont < containers->GetEntries(); icont++){
    cname = dynamic_cast<TObjString *>(containers->At(icont));
    ctemp = dynamic_cast<AliCFContainer *>(fContainers->FindObject(cname->String().Data()));
    if(!ctemp){
      AliWarning(Form("Container %s not found. It will be unprocessed", cname->String().Data()));
      continue;
    }
    nStepMerged += ctemp->GetNStep(); 
  }
  AliInfo("Please Ignore the messgae comming from AliCFContainer!");
  Int_t *dummyBinning = new Int_t[fNVars];
  for(UInt_t ibin = 0; ibin < fNVars; ibin++) dummyBinning[ibin] = 1;
  AliCFContainer *cmerged = new AliCFContainer(name, title, nStepMerged, fNVars, dummyBinning);
  delete[] dummyBinning;
  // Fill container with content
  AliInfo("Filling new container");
  Int_t cstep = 0;
  for(Int_t icont = 0; icont < containers->GetEntries(); icont++){
    cname = dynamic_cast<TObjString *>(containers->At(icont));
    ctemp = dynamic_cast<AliCFContainer *>(fContainers->FindObject(cname->String().Data()));
    if(!ctemp) continue;
    for(Int_t istep = 0; istep < ctemp->GetNStep(); istep++)
      cmerged->SetGrid(cstep++, new AliCFGridSparse(*ctemp->GetGrid(istep)));
  }
  return cmerged;
}

//__________________________________________________________________
void AliHFEcontainer::SetStepTitle(const Char_t *contname, const Char_t *steptitle, UInt_t step){
  //
  // Set title for given analysis step in container with name contname
  //
  AliCFContainer *cont = GetCFContainer(contname);
  if(!cont) return;
  if(step >= static_cast<UInt_t>(cont->GetNStep())) return;
  cont->SetStepTitle(step, steptitle);
}

//__________________________________________________________________
void AliHFEcontainer::MakeLinearBinning(UInt_t var, UInt_t nBins, Double_t begin, Double_t end){
  //
  // Set Linear binning for the given container
  //
  AliHFEvarInfo *myvar = dynamic_cast<AliHFEvarInfo *>(fVariables->UncheckedAt(var));
  if(myvar) myvar->SetBinning(nBins, AliHFEtools::MakeLinearBinning(nBins, begin, end));
}

//__________________________________________________________________
void AliHFEcontainer::MakeLogarithmicBinning(UInt_t var, UInt_t nBins, Double_t begin, Double_t end){
  //
  // Set Logarithmic binning for the given container
  //
  AliHFEvarInfo *myvar = dynamic_cast<AliHFEvarInfo *>(fVariables->UncheckedAt(var));
  if(myvar) myvar->SetBinning(nBins, AliHFEtools::MakeLogarithmicBinning(nBins, begin, end));
}

//__________________________________________________________________
void AliHFEcontainer::MakeUserDefinedBinning(UInt_t var, UInt_t nBins, const Double_t *binning){
  //
  // Set User defined binning
  //
  AliHFEvarInfo *myvar = dynamic_cast<AliHFEvarInfo *>(fVariables->UncheckedAt(var));
  if(myvar) myvar->SetBinning(nBins, binning);
}

//__________________________________________________________________
void AliHFEcontainer::SetVariableName(UInt_t var, const Char_t *varname){
  //
  // Variable name
  // 
  AliHFEvarInfo *myvar = dynamic_cast<AliHFEvarInfo *>(fVariables->UncheckedAt(var));
  if(myvar) myvar->SetVarName(varname);
}

//__________________________________________________________________
Int_t AliHFEcontainer::GetNumberOfCFContainers() const{
  //
  // Get the number of entries
  //
  return fContainers->GetEntries();
}

//__________________________________________________________________
void AliHFEcontainer::Print(const Option_t *)const{
  //
  // Print Container Status
  //
  std::cout << "Container status: " << std::endl;
  std::cout << "=====================================================\n";
  std::cout << "Number of variables: " << fNVars << std::endl;
  if(fNVars){
    UInt_t nVars = fVariables ? fVariables->GetEntriesFast() : 0;
    if(nVars != fNVars)
      std::cout << "Inconsistency in number of Variables [" << fNVars << "|" << nVars << "]" << std::endl;
    AliHFEvarInfo *var = NULL;
    if(fVariables){
      for(UInt_t ivar = 0; ivar < fNVars; ivar++){
        var = dynamic_cast<AliHFEvarInfo *>(fVariables->UncheckedAt(ivar));
        if(var)
          std::cout << "Variable " << ivar << ": Name: " << var->GetVarName()->Data() << ", Number of Bins: " << var->GetNumberOfBins() << std::endl;
      }
    }
  }
  std::cout << std::endl;

  // Print CF Containers:
  if(fContainers){
    std::cout << "Containers[" << fContainers->GetEntries() << "]: "<< std::endl;
    std::cout << "=====================================================\n";
    for(Int_t icont = 0; icont < fContainers->GetEntries(); icont++){
      AliCFContainer *c = dynamic_cast<AliCFContainer *>(fContainers->At(icont));
      if(c){
        std::cout << "Name: " << c->GetName() << ", Title: "  << c->GetTitle() << std::endl;
        for(Int_t istep = 0; istep < c->GetNStep(); istep++)
          std::cout << "Step " << istep << ": Title " << c->GetStepTitle(istep) << std::endl;
      }
      std::cout << "------------------------------------------------------\n";
    }
  }
  std::cout << "Number of Events: " << fNEvents << std::endl;
}

//------------------------------------ Content of class AliHFEvarInfo -----------------------------------
//__________________________________________________________________
AliHFEcontainer::AliHFEvarInfo::AliHFEvarInfo():
  TObject(),
  fVarName(NULL),
  fBinning(NULL)
{
  // Default constructor
  fBinning = new TArrayD;
  fVarName = new TString;
}

//__________________________________________________________________
AliHFEcontainer::AliHFEvarInfo::AliHFEvarInfo(const Char_t *name):
  TObject(),
  fVarName(NULL),
  fBinning(NULL)
{
  fBinning = new TArrayD;
  fVarName = new TString(name);
}

//__________________________________________________________________
AliHFEcontainer::AliHFEvarInfo::AliHFEvarInfo(const AliHFEvarInfo &ref):
  TObject(ref),
  fVarName(NULL),
  fBinning(NULL)
{
  //
  // Constructor
  //
  fVarName = new TString(*(ref.fVarName));
  fBinning = new TArrayD(*(ref.fBinning));
}

//__________________________________________________________________
AliHFEcontainer::AliHFEvarInfo &AliHFEcontainer::AliHFEvarInfo::operator=(const AliHFEvarInfo &ref){
  //
  // Assignment operator
  //
  TObject::operator=(ref);
  *fVarName = *(ref.fVarName);
  *fBinning = *(ref.fBinning);
  return *this;
}

//__________________________________________________________________
AliHFEcontainer::AliHFEvarInfo::~AliHFEvarInfo(){
  //
  // Destructor
  //
  delete fVarName;
  delete fBinning;
}

//__________________________________________________________________
void AliHFEcontainer::AliHFEvarInfo::SetVarName(const Char_t *name){
  //
  // Setter for var name
  //
  *fVarName = name;
}

//__________________________________________________________________
void AliHFEcontainer::AliHFEvarInfo::SetBinning(UInt_t nBins, const Double_t *content){
  // Setter for binning
  //
  fBinning->Set(nBins + 1, content);
}

