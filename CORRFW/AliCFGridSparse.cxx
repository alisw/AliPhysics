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
// AliCFGridSparse Class                                              //
// Class to accumulate data on an N-dimensional grid, to be used      //
// as input to get corrections for Reconstruction & Trigger efficiency// 
// Based on root THnSparse                                            //
// -- Author : S.Arcelli                                              //
// Still to be done:                                                  //
// --Interpolate among bins in a range                                // 
//--------------------------------------------------------------------//
//
//
#include "AliCFGridSparse.h"
#include "THnSparse.h"
#include "AliLog.h"
#include "TMath.h"
#include "TROOT.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TAxis.h"
#include "AliCFUnfolding.h"

//____________________________________________________________________
ClassImp(AliCFGridSparse)

//____________________________________________________________________
AliCFGridSparse::AliCFGridSparse() : 
  AliCFFrame(),
  fSumW2(kFALSE),
  fData(0x0)
{
  // default constructor
}
//____________________________________________________________________
AliCFGridSparse::AliCFGridSparse(const Char_t* name, const Char_t* title) : 
  AliCFFrame(name,title),
  fSumW2(kFALSE),
  fData(0x0)
{
  // default constructor
}
//____________________________________________________________________
AliCFGridSparse::AliCFGridSparse(const Char_t* name, const Char_t* title, Int_t nVarIn, const Int_t * nBinIn) :  
  AliCFFrame(name,title),
  fSumW2(kFALSE),
  fData(0x0)
{
  //
  // main constructor
  //

  fData = new THnSparseF(name,title,nVarIn,nBinIn);
}

//____________________________________________________________________
AliCFGridSparse::~AliCFGridSparse()
{
  //
  // destructor
  //
  if (fData) delete fData;
}

//____________________________________________________________________
AliCFGridSparse::AliCFGridSparse(const AliCFGridSparse& c) :
  AliCFFrame(c),
  fSumW2(kFALSE),
  fData(0x0)
{
  //
  // copy constructor
  //
  ((AliCFGridSparse &)c).Copy(*this);
}

//____________________________________________________________________
AliCFGridSparse& AliCFGridSparse::operator=(const AliCFGridSparse &c)
{
  //
  // assigment operator
  //
  if (this != &c) c.Copy(*this);
  return *this;
} 

//____________________________________________________________________
void AliCFGridSparse::SetBinLimits(Int_t ivar, Double_t min, Double_t max)
{
  //
  // set a uniform binning for variable ivar
  //
  Int_t nBins = GetNBins(ivar);
  Double_t * array = new Double_t[nBins+1];
  for (Int_t iEdge=0; iEdge<=nBins; iEdge++) array[iEdge] = min + iEdge * (max-min)/nBins ;
  fData->SetBinEdges(ivar, array);
  delete [] array ;
} 

//____________________________________________________________________
void AliCFGridSparse::SetBinLimits(Int_t ivar, const Double_t *array)
{
  //
  // setting the arrays containing the bin limits 
  //
  fData->SetBinEdges(ivar, array);
} 

//____________________________________________________________________
void AliCFGridSparse::Fill(const Double_t *var, Double_t weight)
{
  //
  // Fill the grid,
  // given a set of values of the input variable, 
  // with weight (by default w=1)
  //
  fData->Fill(var,weight);
}

//___________________________________________________________________
AliCFGridSparse* AliCFGridSparse::MakeSlice(Int_t nVars, const Int_t* vars, const Double_t* varMin, const Double_t* varMax, Bool_t useBins) const
{
  //
  // projects the grid on the nVars dimensions defined in vars.
  // axis ranges can be defined in arrays varMin, varMax
  // If useBins=true, varMin and varMax are taken as bin numbers
  //

  // binning for new grid
  Int_t* bins = new Int_t[nVars];
  for (Int_t iVar=0; iVar<nVars; iVar++) {
    bins[iVar] = GetNBins(vars[iVar]);
  }
  
  // create new grid sparse
  AliCFGridSparse* out = new AliCFGridSparse(fName,fTitle,nVars,bins);

  //set the range in the THnSparse to project
  THnSparse* clone = ((THnSparse*)fData->Clone());
  if (varMin && varMax) {
    for (Int_t iAxis=0; iAxis<GetNVar(); iAxis++) {
      SetAxisRange(clone->GetAxis(iAxis),varMin[iAxis],varMax[iAxis],useBins);
    }
  }
  else AliInfo("Keeping same axis ranges");

  out->SetGrid(clone->Projection(nVars,vars));
  delete [] bins;
  delete clone;
  return out;
}


//____________________________________________________________________
Float_t AliCFGridSparse::GetBinCenter(Int_t ivar, Int_t ibin) const
{
  //
  // Returns the center of specified bin for variable axis ivar
  // 
  
  return (Float_t) fData->GetAxis(ivar)->GetBinCenter(ibin);
}

//____________________________________________________________________
Float_t AliCFGridSparse::GetBinSize(Int_t ivar, Int_t ibin) const
{
  //
  // Returns the size of specified bin for variable axis ivar
  // 
  
  return (Float_t) fData->GetAxis(ivar)->GetBinUpEdge(ibin) - fData->GetAxis(ivar)->GetBinLowEdge(ibin);
}

//____________________________________________________________________
Float_t AliCFGridSparse::GetEntries() const
{
  //
  // total entries (including overflows and underflows)
  //

  return fData->GetEntries();
}

//____________________________________________________________________
Float_t AliCFGridSparse::GetElement(Long_t index) const
{
  //
  // Returns content of grid element index 
  //
  
  return fData->GetBinContent(index);
}
//____________________________________________________________________
Float_t AliCFGridSparse::GetElement(const Int_t *bin) const
{
  //
  // Get the content in a bin corresponding to a set of bin indexes
  //
  return fData->GetBinContent(bin);

}  
//____________________________________________________________________
Float_t AliCFGridSparse::GetElement(const Double_t *var) const
{
  //
  // Get the content in a bin corresponding to a set of input variables
  //

  Long_t index = fData->GetBin(var,kFALSE);
  if (index<0) return 0.;
  return fData->GetBinContent(index);
} 

//____________________________________________________________________
Float_t AliCFGridSparse::GetElementError(Long_t index) const
{
  //
  // Returns the error on the content 
  //

  return fData->GetBinContent(index);
}
//____________________________________________________________________
Float_t AliCFGridSparse::GetElementError(const Int_t *bin) const
{
 //
  // Get the error in a bin corresponding to a set of bin indexes
  //
  return fData->GetBinError(bin);

}  
//____________________________________________________________________
Float_t AliCFGridSparse::GetElementError(const Double_t *var) const
{
  //
  // Get the error in a bin corresponding to a set of input variables
  //

  Long_t index=fData->GetBin(var,kFALSE); //this is the THnSparse index (do not allocate new cells if content is empy)
  if (index<0) return 0.;
  return fData->GetBinError(index);
} 

//____________________________________________________________________
void AliCFGridSparse::SetElement(Long_t index, Float_t val)
{
  //
  // Sets grid element value
  //
  Int_t* bin = new Int_t[GetNVar()];
  fData->GetBinContent(index,bin); //affects the bin coordinates
  SetElement(bin,val);
  delete [] bin ;
}

//____________________________________________________________________
void AliCFGridSparse::SetElement(const Int_t *bin, Float_t val)
{
  //
  // Sets grid element of bin indeces bin to val
  //
  fData->SetBinContent(bin,val);
}
//____________________________________________________________________
void AliCFGridSparse::SetElement(const Double_t *var, Float_t val) 
{
  //
  // Set the content in a bin to value val corresponding to a set of input variables
  //
  Long_t index=fData->GetBin(var,kTRUE); //THnSparse index: allocate the cell
  Int_t *bin = new Int_t[GetNVar()];
  fData->GetBinContent(index,bin); //trick to access the array of bins
  SetElement(bin,val);
  delete [] bin;
}

//____________________________________________________________________
void AliCFGridSparse::SetElementError(Long_t index, Float_t val)
{
  //
  // Sets grid element iel error to val (linear indexing) in AliCFFrame
  //
  Int_t *bin = new Int_t[GetNVar()];
  fData->GetBinContent(index,bin);
  SetElementError(bin,val);
  delete [] bin;
}

//____________________________________________________________________
void AliCFGridSparse::SetElementError(const Int_t *bin, Float_t val)
{
  //
  // Sets grid element error of bin indeces bin to val
  //
  fData->SetBinError(bin,val);
}
//____________________________________________________________________
void AliCFGridSparse::SetElementError(const Double_t *var, Float_t val) 
{
  //
  // Set the error in a bin to value val corresponding to a set of input variables
  //
  Long_t index=fData->GetBin(var); //THnSparse index
  Int_t *bin = new Int_t[GetNVar()];
  fData->GetBinContent(index,bin); //trick to access the array of bins
  SetElementError(bin,val);
  delete [] bin;
}

//____________________________________________________________________
void AliCFGridSparse::SumW2()
{
  //
  //set calculation of the squared sum of the weighted entries
  //
  if(!fSumW2){
    fData->CalculateErrors(kTRUE); 
  }
  fSumW2=kTRUE;
}

//____________________________________________________________________
void AliCFGridSparse::Add(const AliCFGridSparse* aGrid, Double_t c)
{
  //
  //add aGrid to the current one
  //

  if (aGrid->GetNVar() != GetNVar()){
    AliError("Different number of variables, cannot add the grids");
    return;
  } 
  
  if (!fSumW2  && aGrid->GetSumW2()) SumW2();
  fData->Add(aGrid->GetGrid(),c);
}

//____________________________________________________________________
void AliCFGridSparse::Add(const AliCFGridSparse* aGrid1, const AliCFGridSparse* aGrid2, Double_t c1,Double_t c2)
{
  //
  //Add aGrid1 and aGrid2 and deposit the result into the current one
  //

  if (GetNVar() != aGrid1->GetNVar() || GetNVar() != aGrid2->GetNVar()) {
    AliInfo("Different number of variables, cannot add the grids");
    return;
  } 
  
  if (!fSumW2  && (aGrid1->GetSumW2() || aGrid2->GetSumW2())) SumW2();

  fData->Reset();
  fData->Add(aGrid1->GetGrid(),c1);
  fData->Add(aGrid2->GetGrid(),c2);
}

//____________________________________________________________________
void AliCFGridSparse::Multiply(const AliCFGridSparse* aGrid, Double_t c)
{
  //
  // Multiply aGrid to the current one
  //

  if (aGrid->GetNVar() != GetNVar()) {
    AliError("Different number of variables, cannot multiply the grids");
    return;
  } 
  
  if(!fSumW2  && aGrid->GetSumW2()) SumW2();
  THnSparse *h = aGrid->GetGrid();
  fData->Multiply(h);
  fData->Scale(c);
}

//____________________________________________________________________
void AliCFGridSparse::Multiply(const AliCFGridSparse* aGrid1, const AliCFGridSparse* aGrid2, Double_t c1,Double_t c2)
{
  //
  //Multiply aGrid1 and aGrid2 and deposit the result into the current one
  //

  if (GetNVar() != aGrid1->GetNVar() || GetNVar() != aGrid2->GetNVar()) {
    AliError("Different number of variables, cannot multiply the grids");
    return;
  }
  
  if(!fSumW2  && (aGrid1->GetSumW2() || aGrid2->GetSumW2())) SumW2();

  fData->Reset();
  THnSparse *h1 = aGrid1->GetGrid();
  THnSparse *h2 = aGrid2->GetGrid();
  h2->Multiply(h1);
  h2->Scale(c1*c2);
  fData->Add(h2);
}

//____________________________________________________________________
void AliCFGridSparse::Divide(const AliCFGridSparse* aGrid, Double_t c)
{
  //
  // Divide aGrid to the current one
  //

  if (aGrid->GetNVar() != GetNVar()) {
    AliError("Different number of variables, cannot divide the grids");
    return;
  } 
  
  if (!fSumW2  && aGrid->GetSumW2()) SumW2();

  THnSparse *h1 = aGrid->GetGrid();
  THnSparse *h2 = (THnSparse*)fData->Clone();
  fData->Divide(h2,h1);
  fData->Scale(c);
}

//____________________________________________________________________
void AliCFGridSparse::Divide(const AliCFGridSparse* aGrid1, const AliCFGridSparse* aGrid2, Double_t c1,Double_t c2, Option_t *option)
{
  //
  //Divide aGrid1 and aGrid2 and deposit the result into the current one
  //binomial errors are supported
  //

  if (GetNVar() != aGrid1->GetNVar() || GetNVar() != aGrid2->GetNVar()) {
    AliError("Different number of variables, cannot divide the grids");
    return;
  } 
  
  if (!fSumW2  && (aGrid1->GetSumW2() || aGrid2->GetSumW2())) SumW2();

  THnSparse *h1= aGrid1->GetGrid();
  THnSparse *h2= aGrid2->GetGrid();
  fData->Divide(h1,h2,c1,c2,option);
}


//____________________________________________________________________
void AliCFGridSparse::Rebin(const Int_t* group)
{
  //
  // rebin the grid according to Rebin() as in THnSparse
  // Please notice that the original number of bins on
  // a given axis has to be divisible by the rebin group.
  //

  for(Int_t i=0;i<GetNVar();i++){
    if (group[i]!=1) AliInfo(Form(" merging bins along dimension %i in groups of %i bins", i,group[i]));
  }

  THnSparse *rebinned =fData->Rebin(group);
  fData->Reset();
  fData = rebinned;
}
//____________________________________________________________________
void AliCFGridSparse::Scale(Long_t index, const Double_t *fact)
{
  //
  //scale content of a certain cell by (positive) fact (with error)
  //

  if (GetElement(index)==0 || fact[0]==0) return;

  Double_t in[2], out[2];
  in[0]=GetElement(index);
  in[1]=GetElementError(index);
  GetScaledValues(fact,in,out);
  SetElement(index,out[0]);
  if (fSumW2) SetElementError(index,out[1]);
}
//____________________________________________________________________
void AliCFGridSparse::Scale(const Int_t *bin, const Double_t *fact)
{
  //
  //scale content of a certain cell by (positive) fact (with error)
  //
  if(GetElement(bin)==0 || fact[0]==0)return;

  Double_t in[2], out[2];
  in[0]=GetElement(bin);
  in[1]=GetElementError(bin);
  GetScaledValues(fact,in,out);
  SetElement(bin,out[0]);
  if(fSumW2)SetElementError(bin,out[1]);
  
}
//____________________________________________________________________
void AliCFGridSparse::Scale(const Double_t *var, const Double_t *fact) 
{
  //
  //scale content of a certain cell by (positive) fact (with error)
  //
  if(GetElement(var)==0 || fact[0]==0)return;

  Double_t in[2], out[2];
  in[0]=GetElement(var);
  in[1]=GetElementError(var);
  GetScaledValues(fact,in,out);
  SetElement(var,out[0]);
  if(fSumW2)SetElementError(var,out[1]);
  
}
//____________________________________________________________________
void AliCFGridSparse::Scale(const Double_t* fact)
{
  //
  //scale contents of the whole grid by fact
  //

  for (Long_t iel=0; iel<GetNFilledBins(); iel++) {
    Scale(iel,fact);
  }
}
//____________________________________________________________________
Long_t AliCFGridSparse::GetEmptyBins() const {
  //
  // Get empty bins 
  //

  return (GetNBinsTotal() - GetNFilledBins()) ;
} 

//____________________________________________________________________
Int_t AliCFGridSparse::CheckStats(Double_t thr) const
{
  //
  // Count the cells below a certain threshold
  //
  Int_t ncellsLow=0;
  for (Int_t i=0; i<GetNBinsTotal(); i++) {
    if (GetElement(i)<thr) ncellsLow++;
  }
  return ncellsLow;
}

//_____________________________________________________________________
Double_t AliCFGridSparse::GetIntegral() const 
{
  //
  // Get full Integral
  //
  return fData->ComputeIntegral();  
} 

//____________________________________________________________________
Long64_t AliCFGridSparse::Merge(TCollection* list)
{
  //
  // Merge a list of AliCFGridSparse with this (needed for PROOF). 
  // Returns the number of merged objects (including this).
  //

  if (!list)
    return 0;
  
  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;
  
  Int_t count = 0;
  while ((obj = iter->Next())) {
    AliCFGridSparse* entry = dynamic_cast<AliCFGridSparse*> (obj);
    if (entry == 0) 
      continue;
    this->Add(entry);
    count++;
  }

  return count+1;
}

//____________________________________________________________________
void AliCFGridSparse::GetScaledValues(const Double_t *fact, const Double_t *in, Double_t *out) const{
  //
  // scale input *in and its error by (positive) fact (with error)
  // and erite it to *out
  //
  out[0]=in[0]*fact[0];
  out[1]=TMath::Sqrt(in[1]*in[1]/in[0]/in[0]
		     +fact[1]*fact[1]/fact[0]/fact[0])*out[0];
    
}

//____________________________________________________________________
void AliCFGridSparse::Copy(TObject& c) const
{
  //
  // copy function
  //
  AliCFFrame::Copy(c);
  AliCFGridSparse& target = (AliCFGridSparse &) c;
  target.fSumW2 = fSumW2 ;
  if (fData) {
    target.fData = (THnSparse*)fData->Clone();
  }
}

//____________________________________________________________________
TH1* AliCFGridSparse::Slice(Int_t iVar1, Int_t iVar2, Int_t iVar3, const Double_t *varMin, const Double_t *varMax, Bool_t useBins) const
{
  //
  // return a slice on variables iVar1 (and optionnally iVar2 (and iVar3)) while axis ranges are defined with varMin,varMax
  // arrays varMin and varMax contain the min and max values of each variable.
  // therefore varMin and varMax must have their dimensions equal to GetNVar()
  // If useBins=true, varMin and varMax are taken as bin numbers
  // if varmin or varmax point to null, all the range is taken, including over- and underflows

  THnSparse* clone = (THnSparse*)fData->Clone();
  if (varMin != 0x0 && varMax != 0x0) {
    for (Int_t iAxis=0; iAxis<GetNVar(); iAxis++) SetAxisRange(clone->GetAxis(iAxis),varMin[iAxis],varMax[iAxis],useBins);
  }

  TH1* projection = 0x0 ;
  TString name,title;
  GetProjectionName (name ,iVar1,iVar2,iVar3);
  GetProjectionTitle(title,iVar1,iVar2,iVar3);

  if (iVar3<0) {
    if (iVar2<0) {
      if (iVar1 >= GetNVar() || iVar1 < 0 ) {
	AliError("Non-existent variable, return NULL");
	return 0x0;
      }
      projection = (TH1D*)clone->Projection(iVar1); 
      projection->SetTitle(Form("%s_proj-%s",GetTitle(),GetVarTitle(iVar1)));
      for (Int_t iBin=1; iBin<=projection->GetNbinsX(); iBin++) {
        Int_t origBin = GetAxis(iVar1)->GetFirst()+iBin-1;
	TString binLabel = GetAxis(iVar1)->GetBinLabel(origBin) ;
	if (binLabel.CompareTo("") != 0) projection->GetXaxis()->SetBinLabel(iBin,binLabel);
      }
    }
    else {
      if (iVar1 >= GetNVar() || iVar1 < 0 ||
	  iVar2 >= GetNVar() || iVar2 < 0 ) {
	AliError("Non-existent variable, return NULL");
	return 0x0;
      }
      projection = (TH2D*)clone->Projection(iVar2,iVar1); 
      for (Int_t iBin=1; iBin<=projection->GetNbinsX(); iBin++) {
        Int_t origBin = GetAxis(iVar1)->GetFirst()+iBin-1;
	TString binLabel = GetAxis(iVar1)->GetBinLabel(origBin) ;
	if (binLabel.CompareTo("") != 0) projection->GetXaxis()->SetBinLabel(iBin,binLabel);
      }
      for (Int_t iBin=1; iBin<=projection->GetNbinsY(); iBin++) {
        Int_t origBin = GetAxis(iVar2)->GetFirst()+iBin-1;
	TString binLabel = GetAxis(iVar2)->GetBinLabel(origBin) ;
	if (binLabel.CompareTo("") != 0) projection->GetYaxis()->SetBinLabel(iBin,binLabel);
      }
    }
  }
  else {
    if (iVar1 >= GetNVar() || iVar1 < 0 ||
	iVar2 >= GetNVar() || iVar2 < 0 ||
	iVar3 >= GetNVar() || iVar3 < 0 ) {
      AliError("Non-existent variable, return NULL");
      return 0x0;
    }
    projection = (TH3D*)clone->Projection(iVar1,iVar2,iVar3); 
    for (Int_t iBin=1; iBin<=projection->GetNbinsX(); iBin++) {
      Int_t origBin = GetAxis(iVar1)->GetFirst()+iBin-1;
      TString binLabel = GetAxis(iVar1)->GetBinLabel(origBin) ;
      if (binLabel.CompareTo("") != 0) projection->GetXaxis()->SetBinLabel(iBin,binLabel);
    }
    for (Int_t iBin=1; iBin<=projection->GetNbinsY(); iBin++) {
      Int_t origBin = GetAxis(iVar2)->GetFirst()+iBin-1;
      TString binLabel = GetAxis(iVar2)->GetBinLabel(origBin) ;
      if (binLabel.CompareTo("") != 0) projection->GetYaxis()->SetBinLabel(iBin,binLabel);
    }
    for (Int_t iBin=1; iBin<=projection->GetNbinsZ(); iBin++) {
      Int_t origBin = GetAxis(iVar3)->GetFirst()+iBin-1;
      TString binLabel = GetAxis(iVar3)->GetBinLabel(origBin) ;
      if (binLabel.CompareTo("") != 0) projection->GetZaxis()->SetBinLabel(iBin,binLabel);
    }
  }
  
  projection->SetName (name .Data());
  projection->SetTitle(title.Data());

  delete clone;
  return projection ;
}

//____________________________________________________________________
void AliCFGridSparse::SetAxisRange(TAxis* axis, Double_t min, Double_t max, Bool_t useBins) const {
  //
  // sets axis range, and forces bit TAxis::kAxisRange
  //
  
  if (useBins) axis->SetRange    ((Int_t)min,(Int_t)max);
  else         axis->SetRangeUser(       min,       max);
  //axis->SetBit(TAxis::kAxisRange); // uncomment when ROOT TAxis is fixed
}

//____________________________________________________________________
void AliCFGridSparse::SetRangeUser(Int_t iVar, Double_t varMin, Double_t varMax, Bool_t useBins) const {
  //
  // set range of axis iVar. 
  //
  SetAxisRange(fData->GetAxis(iVar),varMin,varMax,useBins);
	//AliInfo(Form("AliCFGridSparse axis %d range has been modified",iVar));
	TAxis* currAxis = fData->GetAxis(iVar);
  TString outString = Form("%s new range: %.1f < %s < %.1f", GetName(), currAxis->GetBinLowEdge(currAxis->GetFirst()), currAxis->GetTitle(), currAxis->GetBinUpEdge(currAxis->GetLast()));
  TString binLabel = currAxis->GetBinLabel(currAxis->GetFirst());
  if ( ! binLabel.IsNull() ) {
    outString += " ( ";
    for ( Int_t ibin = currAxis->GetFirst(); ibin <= currAxis->GetLast(); ibin++ ) {
      outString += Form("%s ", currAxis->GetBinLabel(ibin));
    }
    outString += ")";
  }
  AliWarning(outString.Data());

}

//____________________________________________________________________
void AliCFGridSparse::SetRangeUser(const Double_t *varMin, const Double_t *varMax, Bool_t useBins) const {
  //
  // set range of every axis. varMin and varMax must be of dimension GetNVar()
  //
  for (Int_t iAxis=0; iAxis<GetNVar() ; iAxis++) { // set new range for every axis
    SetRangeUser(iAxis,varMin[iAxis],varMax[iAxis], useBins);
  }
  AliInfo("AliCFGridSparse axes ranges have been modified");
}

//____________________________________________________________________
Float_t AliCFGridSparse::GetOverFlows(Int_t ivar, Bool_t exclusive) const
{
  //
  // Returns overflows in variable ivar
  // Set 'exclusive' to true for an exclusive check on variable ivar
  //
  Int_t* bin = new Int_t[GetNVar()];
  memset(bin, 0, sizeof(Int_t) * GetNVar());
  Float_t ovfl=0.;
  for (Long64_t i = 0; i < fData->GetNbins(); i++) {
    Double_t v = fData->GetBinContent(i, bin);
    Bool_t add=kTRUE;
    if (exclusive) {
      for(Int_t j=0;j<GetNVar();j++){
	if(ivar==j)continue;
	if((bin[j]==0) || (bin[j]==GetNBins(j)+1))add=kFALSE;
      }
    }
    if(bin[ivar]==GetNBins(ivar)+1 && add) ovfl+=v;
  }

  delete[] bin;
  return ovfl;
}

//____________________________________________________________________
Float_t AliCFGridSparse::GetUnderFlows(Int_t ivar, Bool_t exclusive) const
{
  //
  // Returns exclusive overflows in variable ivar
  // Set 'exclusive' to true for an exclusive check on variable ivar
  //
  Int_t* bin = new Int_t[GetNVar()];
  memset(bin, 0, sizeof(Int_t) * GetNVar());
  Float_t unfl=0.;
  for (Long64_t i = 0; i < fData->GetNbins(); i++) {
    Double_t v = fData->GetBinContent(i, bin);
    Bool_t add=kTRUE;
    if (exclusive) {
      for(Int_t j=0;j<GetNVar();j++){
	if(ivar==j)continue;
	if((bin[j]==0) || (bin[j]==GetNBins(j)+1))add=kFALSE;
      }
    }
    if(bin[ivar]==0 && add) unfl+=v;
  }

  delete[] bin;
  return unfl;
}


//____________________________________________________________________
void AliCFGridSparse::Smooth() {
  //
  // smoothing function: TO USE WITH CARE
  //

  AliInfo("Your GridSparse is going to be smoothed");
  AliInfo(Form("N TOTAL  BINS : %li",GetNBinsTotal()));
  AliInfo(Form("N FILLED BINS : %li",GetNFilledBins()));
  AliCFUnfolding::SmoothUsingNeighbours(fData);
}
