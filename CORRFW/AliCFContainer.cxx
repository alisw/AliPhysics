/* $Id$ */
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
//--------------------------------------------------------------------//
//                                                                    //
// AliCFContainer Class                                               //
// Class to accumulate data on an N-dimensional grids, at different   //
// selection stages. To be used as an input to get corrections for    //
// Reconstruction & Trigger efficiency                                // 
//                                                                    //
// -- Author : S.Arcelli                                              //
//--------------------------------------------------------------------//
//
//
#include <AliLog.h>
#include "AliCFGridSparse.h"
#include "AliCFContainer.h"
#include "TAxis.h"
//____________________________________________________________________
ClassImp(AliCFContainer)

//____________________________________________________________________
AliCFContainer::AliCFContainer() : 
  AliCFFrame(),
  fNStep(0),
  fGrid(0x0)
{
  //
  // default constructor
  //
}
//____________________________________________________________________
AliCFContainer::AliCFContainer(const Char_t* name, const Char_t* title) : 
  AliCFFrame(name,title),
  fNStep(0),
  fGrid(0x0)
{
  // default constructor
}

//____________________________________________________________________
AliCFContainer::AliCFContainer(const Char_t* name, const Char_t* title, const Int_t nSelSteps, const Int_t nVarIn, const Int_t* nBinIn) :  
  AliCFFrame(name,title),
  fNStep(nSelSteps),
  fGrid(0x0)
{
  //
  // main constructor
  //

  // The grids 
  fGrid = new AliCFGridSparse*[fNStep]; //the grids at the various selection steps
  char gname[100];
  for (Int_t istep=0; istep<fNStep; istep++) {
    sprintf(gname,"%s%s%i",GetName(),"_SelStep", istep);
    fGrid[istep] = new AliCFGridSparse(gname,title,nVarIn,nBinIn);
    fGrid[istep]->SumW2();
  }
  AliInfo(Form("Grids created for %d steps required  \n =>  Don't forget to set the bin limits !!",fNStep));
}
//____________________________________________________________________
AliCFContainer::AliCFContainer(const AliCFContainer& c) :
  AliCFFrame(c.fName,c.fTitle),
  fNStep(0),
  fGrid(0x0)
{
  //
  // copy constructor
  //
  c.Copy(*this);
}
//____________________________________________________________________
AliCFContainer::~AliCFContainer()
{
  //
  // destructor
  //
  if (fGrid) delete [] fGrid;
}
//____________________________________________________________________
AliCFContainer &AliCFContainer::operator=(const AliCFContainer &c)
{
  //
  // assigment operator
  //
  if (this != &c) c.Copy(*this);
  return *this;
} 

//____________________________________________________________________
void AliCFContainer::Copy(TObject& c) const
{
  //
  // copy function
  //
  AliCFFrame::Copy(c);
  AliCFContainer& target = (AliCFContainer &) c;
  target.fNStep = fNStep;
  target.fGrid  = new AliCFGridSparse*[fNStep];
  for (Int_t iStep=0; iStep<fNStep; iStep++) {
    if (fGrid[iStep])  target.fGrid[iStep] = new AliCFGridSparse(*(fGrid[iStep]));
  }
}

//____________________________________________________________________
void AliCFContainer::Fill(const Double_t *var, Int_t istep, Double_t weight)
{
  //
  // Fills the grid at selection step istep for a set of values of the 
  // input variables, with a given weight (by default w=1)
  //
  if(istep >= fNStep || istep < 0){
    AliError("Non-existent selection step, grid was not filled");
    return;
  }
  fGrid[istep]->Fill(var,weight);
}
//___________________________________________________________________
TH1D *AliCFContainer::Project(Int_t ivar, Int_t istep) const
{
  //
  // returns 1-D projection along variable ivar at selection step istep
  //
  if(istep >= fNStep || istep < 0){
    AliError("Non-existent selection step, return NULL");
    return 0x0;
  }
  if(ivar >= GetNVar() || ivar < 0){
    AliError("Non-existent variable, return NULL");
    return 0x0;
  }

  return fGrid[istep]->Project(ivar);
}
//___________________________________________________________________
TH2D *AliCFContainer::Project(Int_t ivar1, Int_t ivar2, Int_t istep) const
{
  //
  // returns 2-D projection along variables ivar1,ivar2 at selection step istep
  //
  if(istep >= fNStep || istep < 0){
    AliError("Non-existent selection step, return NULL");
    return 0x0;
  }
  if(ivar1 >= GetNVar() || ivar1 < 0 || ivar2 >= GetNVar() || ivar2 < 0){
    AliError("Non-existent variable, return NULL");
    return 0x0;
  }

  return fGrid[istep]->Project(ivar1,ivar2);
}
//___________________________________________________________________
TH3D *AliCFContainer::Project(Int_t ivar1, Int_t ivar2, Int_t ivar3, Int_t istep) const
{
  //
  // returns 3-D projection along variables ivar1,ivar2,ivar3 
  // at selection step istep
  //
  if(istep >= fNStep || istep < 0){
    AliError("Non-existent selection step, return NULL");
    return 0x0;
  }
  if(ivar1 >= GetNVar() || ivar1 < 0 || 
     ivar2 >= GetNVar() || ivar2 < 0 ||
     ivar3 >= GetNVar() || ivar3 < 0 ) {
    AliError("Non-existent variable, return NULL");
    return 0x0;
  }

  return fGrid[istep]->Project(ivar1,ivar2,ivar3);
}
//___________________________________________________________________
TH1D *AliCFContainer::ShowSlice(Int_t ivar, const Double_t *varMin, const Double_t* varMax, Int_t istep, Bool_t useBins) const
{
  //
  // Make a slice along variable ivar at selection level istep in range [varMin,varMax]
  // If useBins=true, varMin and varMax are taken as bin numbers
  //
  if(istep >= fNStep || istep < 0){
    AliError("Non-existent selection step, return NULL");
    return 0x0;
  }
  if (ivar >= GetNVar() || ivar < 0) {
    AliError("Non-existent variable, return NULL");
    return 0x0;
  }
  return (TH1D*)fGrid[istep]->Slice(ivar,varMin,varMax,useBins);
}
//___________________________________________________________________
TH2D *AliCFContainer::ShowSlice(Int_t ivar1, Int_t ivar2, const Double_t *varMin, const Double_t* varMax, Int_t istep, Bool_t useBins) const
{
  //
  // Make a slice along variables ivar1 and ivar2 at selection level istep in range [varMin,varMax]
  // If useBins=true, varMin and varMax are taken as bin numbers
  //
  if(istep >= fNStep || istep < 0){
    AliError("Non-existent selection step, return NULL");
    return 0x0;
  }
  if (ivar1 >= GetNVar() || ivar1 < 0 || ivar2 >= GetNVar() || ivar2 < 0) {
    AliError("Non-existent variable, return NULL");
    return 0x0;
  }
  return (TH2D*)fGrid[istep]->Slice(ivar1,ivar2,varMin,varMax,useBins);
}
//___________________________________________________________________
TH3D *AliCFContainer::ShowSlice(Int_t ivar1, Int_t ivar2, Int_t ivar3, const Double_t *varMin, const Double_t* varMax, Int_t istep, Bool_t useBins) const
{
  //
  // Make a slice along variables ivar1, ivar2and ivar3 at selection level istep in range [varMin,varMax]
  // If useBins=true, varMin and varMax are taken as bin numbers
  //
  if(istep >= fNStep || istep < 0){
    AliError("Non-existent selection step, return NULL");
    return 0x0;
  }
  if (ivar1 >= GetNVar() || ivar1 < 0 || ivar2 >= GetNVar() || ivar2 < 0 || ivar3 >= GetNVar() || ivar3 < 0) {
    AliError("Non-existent variable, return NULL");
    return 0x0;
  }
  return (TH3D*)fGrid[istep]->Slice(ivar1,ivar2,ivar3,varMin,varMax,useBins);
}
//____________________________________________________________________
AliCFContainer* AliCFContainer::MakeSlice(Int_t nVars, const Int_t* vars, const Double_t* varMin, const Double_t* varMax, Bool_t useBins) const
{
  //
  // Makes a slice along the "nVars" variables defined in the array "vars[nVars]" for all the container steps.
  // The ranges of ALL the container variables must be defined in the array varMin[GetNVar()] and varMax[GetNVar()]
  // The function returns a new container of nVars variables.
  // If useBins=true, varMin and varMax are taken as bin numbers
  //
  Int_t* steps = new Int_t[fNStep];
  for (Int_t iStep=0;iStep<fNStep;iStep++) steps[iStep]=iStep;
  AliCFContainer* out = MakeSlice(nVars,vars,varMin,varMax,fNStep,steps,useBins);
  delete [] steps ;
  return out;
}

//____________________________________________________________________
AliCFContainer* AliCFContainer::MakeSlice(Int_t nVars, const Int_t* vars, const Double_t* varMin, const Double_t* varMax, Int_t nSteps, const Int_t* steps, Bool_t useBins) const
{
  //
  // Makes a slice along the "nVars" variables defined in the array "vars[nVars]" for the given "nSteps" defined in "steps[nSteps]".
  // The ranges of ALL the container variables must be defined in the array varMin[GetNVar()] and varMax[GetNVar()]
  // The function returns a new container of nVars variables.
  // If useBins=true, varMin and varMax are taken as bin numbers
  //

  if (nVars < 1 || nVars > GetNVar())   AliError("Bad number of dimensions required for the slice");
  if (nSteps< 1 || nSteps> fNStep)  AliError("Bad number of steps required for the slice");

  AliInfo(Form("Making a slice in %d dimension(s)",nVars));

  // create the output grids
  AliCFGridSparse** grids = new AliCFGridSparse*[nSteps] ;
  for (Int_t iStep=0; iStep<nSteps; iStep++) grids[iStep] = fGrid[steps[iStep]]->Project(nVars,vars,varMin,varMax,useBins);
  
  TAxis ** axis = new TAxis*[nVars];
  for (Int_t iVar=0; iVar<nVars; iVar++) axis[iVar] = ((AliCFGridSparse*)grids[0])->GetGrid()->GetAxis(iVar); //same axis for every grid

  //define new binning for new container
  Int_t* bins=new Int_t[nVars];
  for (Int_t iVar=0; iVar<nVars; iVar++) bins[iVar] = axis[iVar]->GetNbins();

  AliCFContainer* out = new AliCFContainer(fName,fTitle,nSteps,nVars,bins);

  //set the bin limits
  for (Int_t iVar=0; iVar<nVars; iVar++) {
    Int_t nBins = bins[iVar];
    Double_t *array = new Double_t[nBins+1];
    for (Int_t iBin=1; iBin<=nBins; iBin++) {
      array[iBin-1] = axis[iVar]->GetBinLowEdge(iBin);
    }
    array[nBins] = axis[iVar]->GetBinUpEdge(nBins);
    out->SetBinLimits(iVar,array);
    delete [] array;
  }

  //set grid for the given steps
  for (Int_t iStep=0; iStep<nSteps; iStep++) out->SetGrid(iStep,grids[iStep]);

  delete [] bins;
  delete [] axis ;
  return out;
}

//____________________________________________________________________
Long64_t AliCFContainer::Merge(TCollection* list)
{
  // Merge a list of AliCorrection objects with this (needed for
  // PROOF). 
  // Returns the number of merged objects (including this).

  if (!list)
    return 0;
  
  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;
  
  Int_t count = 0;
  while ((obj = iter->Next())) {
    AliCFContainer* entry = dynamic_cast<AliCFContainer*> (obj);
    if (entry == 0) 
      continue;
    this->Add(entry);
    count++;
  }

  return count+1;
}

//____________________________________________________________________
void AliCFContainer::Add(const AliCFContainer* aContainerToAdd, Double_t c)
{
  //
  //add the content of container aContainerToAdd to the current one
  //
  if ((aContainerToAdd->GetNStep()      != fNStep)          ||
      (aContainerToAdd->GetNVar()       != GetNVar())       ||
      (aContainerToAdd->GetNBinsTotal() != GetNBinsTotal()))
    {
      AliError("Different number of steps/sensitive variables/grid elements: cannot add the containers");
      return;
    }
  for (Int_t istep=0; istep<fNStep; istep++) {
    fGrid[istep]->Add(aContainerToAdd->GetGrid(istep),c);
  }
}
//____________________________________________________________________
Float_t AliCFContainer::GetOverFlows( Int_t ivar, Int_t istep, Bool_t exclusive) const {
  //
  // Get overflows in variable var at selection level istep
  // Set 'exclusive' to true for an exclusive check on variable ivar
  //
  if(istep >= fNStep || istep < 0){
    AliError("Non-existent selection step, return -1");
    return -1.;
  }
  return fGrid[istep]->GetOverFlows(ivar,exclusive);
} 
//____________________________________________________________________
Float_t AliCFContainer::GetUnderFlows( Int_t ivar, Int_t istep, Bool_t exclusive) const {
  //
  // Get underflows in variable var at selection level istep
  // Set 'exclusive' to true for an exclusive check on variable ivar
  //
  if(istep >= fNStep || istep < 0){
    AliError("Non-existent selection step, return -1");
    return -1.;
  }
  return fGrid[istep]->GetUnderFlows(ivar,exclusive);
} 
//____________________________________________________________________
Float_t AliCFContainer::GetEntries(Int_t istep) const {
  //
  // Get total entries in variable var at selection level istep
  //
  if(istep >= fNStep || istep < 0){
    AliError("Non-existent selection step, return -1");
    return -1.;
  }
  return fGrid[istep]->GetEntries();
} 
// //____________________________________________________________________
// Int_t AliCFContainer::GetEmptyBins(Int_t istep) const {
//   //
//   // Get empty bins in variable var at selection level istep
//   //
//   if(istep >= fNStep || istep < 0){
//     AliError("Non-existent selection step, return -1");
//     return -1;
//   }
//   return fGrid[istep]->GetEmptyBins();
// } 
// //____________________________________________________________________
// Int_t AliCFContainer::GetEmptyBins( Int_t istep, Double_t *varMin, Double_t* varMax) const {
//   //
//   // Get empty bins in a range in variable var at selection level istep
//   //
//   if(istep >= fNStep || istep < 0){
//     AliError("Non-existent selection step, return -1");
//     return -1;
//   }
//   return fGrid[istep]->GetEmptyBins(varMin,varMax);
// } 
//_____________________________________________________________________
Double_t AliCFContainer::GetIntegral( Int_t istep) const 
{
  //
  // Get Integral over the grid at selection level istep
  //
  if(istep >= fNStep || istep < 0){
    AliError("Non-existent selection step, return -1");
    return -1.;
  }
  return fGrid[istep]->GetIntegral();
}
//_____________________________________________________________________
// Double_t AliCFContainer::GetIntegral( Int_t istep, Double_t *varMin, Double_t* varMax ) const 
// {
//   //
//   // Get Integral over the grid in a range at selection level istep
//   //
//   if(istep >= fNStep || istep < 0){
//     AliError("Non-existent selection step, return -1");
//     return -1.;
//   }
//   return fGrid[istep]->GetIntegral(varMin,varMax);
// }
//_____________________________________________________________________
void AliCFContainer::SetRangeUser(Int_t ivar, Double_t varMin, Double_t varMax, Int_t istep) 
{
  //
  // set axis range at step istep
  //
  if (istep >= fNStep || istep < 0){
    AliError("Non-existent selection step");
    return ;
  }
  if (ivar >= GetNVar() || ivar < 0){
    AliError("Non-existent selection var");
    return ;
  }
  fGrid[istep]->SetRangeUser(ivar,varMin,varMax);
}

//_____________________________________________________________________
void AliCFContainer::SetRangeUser(Double_t* varMin, Double_t* varMax, Int_t istep) 
{
  //
  // set all axis ranges at step istep according to varMin and varMax values
  //
  if (istep >= fNStep || istep < 0){
    AliError("Non-existent selection step");
    return ;
  }
  ((AliCFGridSparse*)fGrid[istep])->SetRangeUser(varMin,varMax);
}
