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
  fVariables(NULL),
  fNVars(0),
  fNEvents(0)
{
  //
  // Default constructor
  //
  fContainers = new THashList();
}

//__________________________________________________________________
AliHFEcontainer::AliHFEcontainer(const Char_t *name):
  TNamed(name, ""),
  fContainers(NULL),
  fVariables(NULL),
  fNVars(0),
  fNEvents(0)
{
  //
  // Default constructor
  //
  fContainers = new THashList();
}

//__________________________________________________________________
AliHFEcontainer::AliHFEcontainer(const Char_t *name, UInt_t nVar):
  TNamed(name, ""),
  fContainers(NULL),
  fVariables(NULL),
  fNVars(0),
  fNEvents(0)
{
  //
  // Constructor
  // Setting Number of Variables too
  //
  fContainers = new THashList();
  SetNumberOfVariables(nVar);
}

//__________________________________________________________________
AliHFEcontainer::AliHFEcontainer(const AliHFEcontainer &ref):
  TNamed(ref),
  fContainers(NULL),
  fVariables(NULL),
  fNVars(ref.fNVars),
  fNEvents(ref.fNEvents)
{
  //
  // Copy constructor
  // creates a new object with new containers
  //
  fContainers = new THashList;
  for(Int_t ien = 0; ien < ref.fContainers->GetEntries(); ien++)
    fContainers->Add(new AliCFContainer(*dynamic_cast<AliCFContainer *>(ref.fContainers->At(ien))));
  if(fNVars){
    fVariables = new TObjArray(fNVars);
    for(UInt_t ivar = 0; ivar < fNVars; ivar++)
      fVariables->AddAt(new AliHFEvarInfo(*dynamic_cast<AliHFEvarInfo *>(ref.fVariables->UncheckedAt(ivar))), ivar);
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
  fNVars = ref.fNVars;
  for(Int_t ien = 0; ien < ref.fContainers->GetEntries(); ien++)
    fContainers->Add(new AliCFContainer(*dynamic_cast<AliCFContainer *>(ref.fContainers->At(ien))));
  if(fNVars){
    fVariables = new TObjArray(fNVars);
    for(UInt_t ivar = 0; ivar < fNVars; ivar++)
      fVariables->AddAt(new AliHFEvarInfo(*dynamic_cast<AliHFEvarInfo *>(ref.fVariables->UncheckedAt(ivar))), ivar);
  } else {
    fVariables = NULL;
  }

  return *this;
}

//__________________________________________________________________
AliHFEcontainer::~AliHFEcontainer(){
  //
  // Destructor
  //
  fContainers->Delete();
  delete fContainers;
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

  TIterator *iter = coll->MakeIterator();
  TObject *o = NULL;
  Long64_t count = 0;
  while((o = iter->Next())){
    AliHFEcontainer *cont = dynamic_cast<AliHFEcontainer *>(o);
    if(!cont) continue;

    // Merge the two TObjArrays
    TList containers;
    containers.Add(cont->fContainers);
    fContainers->Merge(&containers);

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
  for(UInt_t ivar = 0; ivar < fNVars; ivar++) nBins[ivar] = (dynamic_cast<AliHFEvarInfo *>(fVariables->UncheckedAt(ivar)))->GetNumberOfBins();
  AliHFEvarInfo *var = NULL;
  AliCFContainer *cont = new AliCFContainer(name, title, nStep, fNVars, nBins);
  for(UInt_t ivar = 0; ivar < fNVars; ivar++){
    var = dynamic_cast<AliHFEvarInfo *>(fVariables->UncheckedAt(ivar));
    cont->SetBinLimits(ivar, var->GetBinning());
    cont->SetVarTitle(ivar, var->GetVarName()->Data());
  }
  delete[] nBins;
  fContainers->Add(cont);
  AliInfo(Form("Container %s created with %d cut steps", name, nStep));
}

//__________________________________________________________________
AliCFContainer *AliHFEcontainer::GetCFContainer(const Char_t *name){
  //
  // Find a given container 
  //
  return dynamic_cast<AliCFContainer *>(fContainers->FindObject(name));
}

//__________________________________________________________________
void AliHFEcontainer::FillCFContainer(const Char_t *name, UInt_t step, Double_t *content){
  //
  // Fill container
  //
  AliCFContainer *cont = GetCFContainer(name);
  if(!cont) return;
  cont->Fill(content, step);
}

//__________________________________________________________________
AliCFContainer *AliHFEcontainer::MakeMergedCFContainer(const Char_t *name, const Char_t *title, const Char_t* contnames){
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
  delete dummyBinning;
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
void AliHFEcontainer::MakeLinearBinning(UInt_t var, UInt_t nBins, Double_t begin, Double_t end){
  //
  // Set Linear binning for the given container
  //
  (dynamic_cast<AliHFEvarInfo *>(fVariables->UncheckedAt(var)))->SetBinning(nBins, AliHFEtools::MakeLinearBinning(nBins, begin, end));
}

//__________________________________________________________________
void AliHFEcontainer::MakeLogarithmicBinning(UInt_t var, UInt_t nBins, Double_t begin, Double_t end){
  //
  // Set Logarithmic binning for the given container
  //
  (dynamic_cast<AliHFEvarInfo *>(fVariables->UncheckedAt(var)))->SetBinning(nBins, AliHFEtools::MakeLogarithmicBinning(nBins, begin, end));
}

//__________________________________________________________________
void AliHFEcontainer::SetVariableName(UInt_t var, const Char_t *varname){
  //
  // Variable name
  // 
  (dynamic_cast<AliHFEvarInfo *>(fVariables->UncheckedAt(var)))->SetVarName(varname);
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
    for(UInt_t ivar = 0; ivar < fNVars; ivar++){
      var = dynamic_cast<AliHFEvarInfo *>(fVariables->UncheckedAt(ivar));
      std::cout << "Variable " << ivar << ": Name: " << var->GetVarName()->Data() << ", Number of Bins: " << var->GetNumberOfBins() << std::endl;
    }
  }
  std::cout << std::endl;

  // Print CF Containers:
  std::cout << "Containers[" << fContainers->GetEntries() << "]: "<< std::endl;
  std::cout << "=====================================================\n";
  for(Int_t icont = 0; icont < fContainers->GetEntries(); icont++){
    AliCFContainer *c = dynamic_cast<AliCFContainer *>(fContainers->At(icont));
    std::cout << "Name: " << c->GetName() << ", Title: "  << c->GetTitle() << std::endl;
    for(Int_t istep = 0; istep < c->GetNStep(); istep++)
      std::cout << "Step " << istep << ": Title " << c->GetStepTitle(istep) << std::endl;
    std::cout << "------------------------------------------------------\n";
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
void AliHFEcontainer::AliHFEvarInfo::SetBinning(UInt_t nBins, Double_t *content){
  // Setter for binning
  //
  fBinning->Set(nBins + 1, content);
}

