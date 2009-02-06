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

//____________________________________________________________________
ClassImp(AliCFGridSparse)

//____________________________________________________________________
AliCFGridSparse::AliCFGridSparse() : 
  AliCFVGrid(),
  fData(0x0)
{
  // default constructor
}
//____________________________________________________________________
AliCFGridSparse::AliCFGridSparse(const Char_t* name, const Char_t* title) : 
  AliCFVGrid(name,title),
  fData(0x0)
{
  // default constructor
}
//____________________________________________________________________
AliCFGridSparse::AliCFGridSparse(const Char_t* name, const Char_t* title, const Int_t nVarIn, const Int_t * nBinIn, const Double_t *binLimitsIn) :  
  AliCFVGrid(name,title,nVarIn,nBinIn,binLimitsIn),
  fData(0x0)
{
  //
  // main constructor
  //

  fData=new THnSparseF(name,title,fNVar,fNVarBins);
  
  if(binLimitsIn){
    for(Int_t ivar=0;ivar<fNVar;ivar++){    
      Int_t nbins=fNVarBins[ivar]+1;
      Double_t *array= new Double_t[nbins];
      for(Int_t i=0;i<nbins;i++){
	array[i]=fVarBinLimits[fOffset[ivar]+i];
      } 
      fData->SetBinEdges(ivar, array);
      delete array;
    }
  }
}
//____________________________________________________________________
AliCFGridSparse::AliCFGridSparse(const AliCFGridSparse& c) : 
  AliCFVGrid(c),
  fData(c.fData)
{
  //
  // copy constructor
  //
  ((AliCFGridSparse &)c).Copy(*this);
}

//____________________________________________________________________
AliCFGridSparse::~AliCFGridSparse()
{
  //
  // destructor
  //
  if(fData) delete fData;
}

//____________________________________________________________________
AliCFGridSparse &AliCFGridSparse::operator=(const AliCFGridSparse &c)
{
  //
  // assigment operator
  //
  if (this != &c)
    ((AliCFGridSparse &) c).Copy(*this);
  return *this;
} 

//____________________________________________________________________
void AliCFGridSparse::SetBinLimits(Int_t ivar, Double_t *array)
{
  //
  // setting the arrays containing the bin limits 
  //
  fData->SetBinEdges(ivar, array);
  //then fill the appropriate array in ALICFFrame, to be able to use 
  //the getter, in case....
  Int_t nbins=fNVarBins[ivar]+1;
  for(Int_t i=0;i<nbins;i++){
    fVarBinLimits[fOffset[ivar]+i] =array[i];
  } 
} 

//____________________________________________________________________
void AliCFGridSparse::Fill(Double_t *var, Double_t weight)
{
  //
  // Fill the grid,
  // given a set of values of the input variable, 
  // with weight (by default w=1)
  //
  fData->Fill(var,weight);
}

//___________________________________________________________________
TH1D *AliCFGridSparse::Project(Int_t ivar) const
{
  //
  // Make a 1D projection along variable ivar 
  //

  TH1D *hist=fData->Projection(ivar);
  return hist;
}
//___________________________________________________________________
TH2D *AliCFGridSparse::Project(Int_t ivar1, Int_t ivar2) const
{
  //
  // Make a 2D projection along variables ivar1 & ivar2 
  //

  TH2D *hist=fData->Projection(ivar2,ivar1); //notice inverted axis (THnSparse uses TH3 2d-projection convention...)
  return hist;

}
//___________________________________________________________________
TH3D *AliCFGridSparse::Project(Int_t ivar1, Int_t ivar2, Int_t ivar3) const
{
  //
  // Make a 3D projection along variables ivar1 & ivar2 & ivar3 
  //

  TH3D *hist=fData->Projection(ivar1,ivar2,ivar3); 
  return hist;

}

//___________________________________________________________________
AliCFGridSparse* AliCFGridSparse::Project(Int_t nVars, Int_t* vars, Double_t* varMin, Double_t* varMax) const
{

  // binning for new grid
  Int_t* bins = new Int_t[nVars];
  for (Int_t iVar=0; iVar<nVars; iVar++) {
    bins[iVar] = fNVarBins[vars[iVar]];
  }
  
  // create new grid sparse
  AliCFGridSparse* out = new AliCFGridSparse(fName,fTitle,nVars,bins);

  //set the range in the THnSparse to project
  THnSparse* clone = ((THnSparse*)fData->Clone());
  for (Int_t iAxis=0; iAxis<fNVar; iAxis++) {
    clone->GetAxis(iAxis)->SetRangeUser(varMin[iAxis],varMax[iAxis]);
  }
  out->SetGrid(clone->Projection(nVars,vars));
  return out;
}

//____________________________________________________________________
Float_t AliCFGridSparse::GetOverFlows(Int_t ivar) const
{
  //
  // Returns exclusive overflows in variable ivar
  //
  Int_t* bin = new Int_t[fNVar];
  memset(bin, 0, sizeof(Int_t) * fNVar);
  Float_t ovfl=0.;
  for (Long64_t i = 0; i < fData->GetNbins(); ++i) {
    Double_t v = fData->GetBinContent(i, bin);
    Bool_t add=kTRUE;
    for(Int_t j=0;j<fNVar;j++){
      if(ivar==j)continue;
      if((bin[j]==0) || (bin[j]==fNVarBins[j]+1))add=kFALSE;
    }
    if(bin[ivar]==fNVarBins[ivar]+1 && add) ovfl+=v;
  }

  delete[] bin;
  return ovfl;
}

//____________________________________________________________________
Float_t AliCFGridSparse::GetUnderFlows(Int_t ivar) const
{
  //
  // Returns exclusive overflows in variable ivar
  //
  Int_t* bin = new Int_t[fNVar];
  memset(bin, 0, sizeof(Int_t) * fNVar);
  Float_t unfl=0.;
  for (Long64_t i = 0; i < fData->GetNbins(); ++i) {
    Double_t v = fData->GetBinContent(i, bin);
    Bool_t add=kTRUE;
    for(Int_t j=0;j<fNVar;j++){
      if(ivar==j)continue;
      if((bin[j]==0) || (bin[j]==fNVarBins[j]+1))add=kFALSE;
    }
    if(bin[ivar]==0 && add) unfl+=v;
  }

  delete[] bin;
  return unfl;
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
Float_t AliCFGridSparse::GetElement(Int_t index) const
{
  //
  // Returns content of grid element index according to the
  // linear indexing in AliCFFrame
  //
  Int_t *bin = new Int_t[fNVar];
  GetBinIndex(index, bin);
  for(Int_t i=0;i<fNVar;i++)fIndex[i]=bin[i]+1; //consistency with AliCFGrid
  Float_t val=GetElement(fIndex);
  delete [] bin;
  return val; 
  
}
//____________________________________________________________________
Float_t AliCFGridSparse::GetElement(Int_t *bin) const
{
  //
  // Get the content in a bin corresponding to a set of bin indexes
  //
  return fData->GetBinContent(bin);

}  
//____________________________________________________________________
Float_t AliCFGridSparse::GetElement(Double_t *var) const
{
  //
  // Get the content in a bin corresponding to a set of input variables
  //

  Long_t index=fData->GetBin(var,kFALSE); //this is the THnSparse index (do not allocate new cells if content is empty)
  if(index<0){
    return 0.;
  }else{
    return fData->GetBinContent(index);
  }
} 

//____________________________________________________________________
Float_t AliCFGridSparse::GetElementError(Int_t index) const
{
  //
  // Returns the error on the content of a bin according to a linear
  // indexing in AliCFFrame
  //

  Int_t *bin = new Int_t[fNVar];
  GetBinIndex(index, bin);
  for(Int_t i=0;i<fNVar;i++)fIndex[i]=bin[i]+1; //consistency with AliCFGrid
  Float_t val=GetElementError(fIndex);
  delete [] bin;
  return val;

}
//____________________________________________________________________
Float_t AliCFGridSparse::GetElementError(Int_t *bin) const
{
 //
  // Get the error in a bin corresponding to a set of bin indexes
  //
  return fData->GetBinError(bin);

}  
//____________________________________________________________________
Float_t AliCFGridSparse::GetElementError(Double_t *var) const
{
  //
  // Get the error in a bin corresponding to a set of input variables
  //

  Long_t index=fData->GetBin(var,kFALSE); //this is the THnSparse index (do not allocate new cells if content is empy)
  if(index<0){
    return 0.;
  }else{
    return fData->GetBinError(index);
  }
} 


//____________________________________________________________________
void AliCFGridSparse::SetElement(Int_t index, Float_t val)
{
  //
  // Sets grid element iel to val (linear indexing) in AliCFFrame
  //
  Int_t *bin = new Int_t[fNVar];
  GetBinIndex(index, bin);
  for(Int_t i=0;i<fNVar;i++)fIndex[i]=bin[i]+1;
  SetElement(fIndex,val);
  delete [] bin;
}
//____________________________________________________________________
void AliCFGridSparse::SetElement(Int_t *bin, Float_t val)
{
  //
  // Sets grid element of bin indeces bin to val
  //
  fData->SetBinContent(bin,val);
}
//____________________________________________________________________
void AliCFGridSparse::SetElement(Double_t *var, Float_t val) 
{
  //
  // Set the content in a bin to value val corresponding to a set of input variables
  //
  Long_t index=fData->GetBin(var); //THnSparse index: allocate the cell
  Int_t *bin = new Int_t[fNVar];
  fData->GetBinContent(index,bin); //trick to access the array of bins
  fData->SetBinContent(bin,val);
  delete [] bin;

}

//____________________________________________________________________
void AliCFGridSparse::SetElementError(Int_t index, Float_t val)
{
  //
  // Sets grid element iel error to val (linear indexing) in AliCFFrame
  //
  Int_t *bin = new Int_t[fNVar];
  GetBinIndex(index, bin);
  for(Int_t i=0;i<fNVar;i++)fIndex[i]=bin[i]+1;
  SetElementError(fIndex,val);
  delete [] bin;
}
//____________________________________________________________________
void AliCFGridSparse::SetElementError(Int_t *bin, Float_t val)
{
  //
  // Sets grid element error of bin indeces bin to val
  //
  fData->SetBinError(bin,val);
}
//____________________________________________________________________
void AliCFGridSparse::SetElementError(Double_t *var, Float_t val) 
{
  //
  // Set the error in a bin to value val corresponding to a set of input variables
  //
  Long_t index=fData->GetBin(var); //THnSparse index
  Int_t *bin = new Int_t[fNVar];
  fData->GetBinContent(index,bin); //trick to access the array of bins
  fData->SetBinError(bin,val);
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
void AliCFGridSparse::Add(AliCFVGrid* aGrid, Double_t c)
{
  //
  //add aGrid to the current one
  //


  if(aGrid->GetNVar()!=fNVar){
    AliInfo("Different number of variables, cannot add the grids");
    return;
  } 
  if(aGrid->GetNDim()!=fNDim){
    AliInfo("Different number of dimensions, cannot add the grids!");
    return;
  } 
  
  if(!fSumW2  && aGrid->GetSumW2())SumW2();


  fData->Add(((AliCFGridSparse*)aGrid)->GetGrid(),c);

}

//____________________________________________________________________
void AliCFGridSparse::Add(AliCFVGrid* aGrid1, AliCFVGrid* aGrid2, Double_t c1,Double_t c2)
{
  //
  //Add aGrid1 and aGrid2 and deposit the result into the current one
  //

  if(fNVar!=aGrid1->GetNVar()|| fNVar!=aGrid2->GetNVar()){
    AliInfo("Different number of variables, cannot add the grids");
    return;
  } 
  if(fNDim!=aGrid1->GetNDim()|| fNDim!=aGrid2->GetNDim()){
    AliInfo("Different number of dimensions, cannot add the grids!");
    return;
  } 
  
  if(!fSumW2  && (aGrid1->GetSumW2() || aGrid2->GetSumW2()))SumW2();


  fData->Reset();
  fData->Add(((AliCFGridSparse*)aGrid1)->GetGrid(),c1);
  fData->Add(((AliCFGridSparse*)aGrid2)->GetGrid(),c2);

}

//____________________________________________________________________
void AliCFGridSparse::Multiply(AliCFVGrid* aGrid, Double_t c)
{
  //
  // Multiply aGrid to the current one
  //


  if(aGrid->GetNVar()!=fNVar){
    AliInfo("Different number of variables, cannot multiply the grids");
    return;
  } 
  if(aGrid->GetNDim()!=fNDim){
    AliInfo("Different number of dimensions, cannot multiply the grids!");
    return;
  } 
  
  if(!fSumW2  && aGrid->GetSumW2())SumW2();

  THnSparse *h= ((AliCFGridSparse*)aGrid)->GetGrid();

  fData->Multiply(h);
  fData->Scale(c);

}

//____________________________________________________________________
void AliCFGridSparse::Multiply(AliCFVGrid* aGrid1, AliCFVGrid* aGrid2, Double_t c1,Double_t c2)
{
  //
  //Multiply aGrid1 and aGrid2 and deposit the result into the current one
  //

  if(fNVar!=aGrid1->GetNVar()|| fNVar!=aGrid2->GetNVar()){
    AliInfo("Different number of variables, cannot multiply the grids");
    return;
  } 
  if(fNDim!=aGrid1->GetNDim()|| fNDim!=aGrid2->GetNDim()){
    AliInfo("Different number of dimensions, cannot multiply the grids!");
    return;
  } 
  
  if(!fSumW2  && (aGrid1->GetSumW2() || aGrid2->GetSumW2()))SumW2();


  fData->Reset();
  THnSparse *h1= ((AliCFGridSparse*)aGrid1)->GetGrid();
  THnSparse *h2= ((AliCFGridSparse*)aGrid2)->GetGrid();
  h2->Multiply(h1);
  h2->Scale(c1*c2);
  fData->Add(h2);
}



//____________________________________________________________________
void AliCFGridSparse::Divide(AliCFVGrid* aGrid, Double_t c)
{
  //
  // Divide aGrid to the current one
  //


  if(aGrid->GetNVar()!=fNVar){
    AliInfo("Different number of variables, cannot divide the grids");
    return;
  } 
  if(aGrid->GetNDim()!=fNDim){
    AliInfo("Different number of dimensions, cannot divide the grids!");
    return;
  } 
  
  if(!fSumW2  && aGrid->GetSumW2())SumW2();

  THnSparse *h= ((AliCFGridSparse*)aGrid)->GetGrid();

  fData->Divide(h);
  fData->Scale(c);

}

//____________________________________________________________________
void AliCFGridSparse::Divide(AliCFVGrid* aGrid1, AliCFVGrid* aGrid2, Double_t c1,Double_t c2, Option_t *option)
{
  //
  //Divide aGrid1 and aGrid2 and deposit the result into the current one
  //bynomial errors are supported
  //

  if(fNVar!=aGrid1->GetNVar()|| fNVar!=aGrid2->GetNVar()){
    AliInfo("Different number of variables, cannot divide the grids");
    return;
  } 
  if(fNDim!=aGrid1->GetNDim()|| fNDim!=aGrid2->GetNDim()){
    AliInfo("Different number of dimensions, cannot divide the grids!");
    return;
  } 
  
  if(!fSumW2  && (aGrid1->GetSumW2() || aGrid2->GetSumW2()))SumW2();


  THnSparse *h1= ((AliCFGridSparse*)aGrid1)->GetGrid();
  THnSparse *h2= ((AliCFGridSparse*)aGrid2)->GetGrid();
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

  for(Int_t i=0;i<fNVar;i++){
    if(group[i]!=1)AliInfo(Form(" merging bins along dimension %i in groups of %i bins", i,group[i]));
  }

  THnSparse *rebinned =fData->Rebin(group);
  fData->Reset();
  fData = rebinned;

  //redefine the needed stuff

  Int_t ndimTot=1;
  Int_t nbinTot=0;

  //number of bins in each dimension, auxiliary variables

  for(Int_t ivar=0;ivar<fNVar;ivar++){
    Int_t nbins = fData->GetAxis(ivar)->GetNbins();
    fNVarBins[ivar]=nbins;
    ndimTot*=fNVarBins[ivar];
    nbinTot+=(fNVarBins[ivar]+1);
    Int_t offset=0;
    for(Int_t i =0;i<ivar;i++)offset+=(fNVarBins[i]+1);      
    fOffset[ivar]=offset;
    Int_t prod=1;
    for(Int_t i=0;i<ivar;i++)prod*=fNVarBins[i];
    fProduct[ivar]=prod;
  }

  fNDim=ndimTot;

  //now the array of bin limits

  delete fVarBinLimits;
  fNVarBinLimits=nbinTot;
  fVarBinLimits=new Double_t[fNVarBinLimits];

  for(Int_t ivar=0;ivar<fNVar;ivar++){
    Double_t low = fData->GetAxis(ivar)->GetXmin();
    Double_t high = fData->GetAxis(ivar)->GetXmax();    
    const TArrayD *xbins = fData->GetAxis(ivar)->GetXbins();
    if (xbins->fN == 0){
      for(Int_t ibin=0;ibin<=fNVarBins[ivar];ibin++){
	fVarBinLimits[ibin+fOffset[ivar]] = low + ibin*(high-low)/((Double_t) fNVarBins[ivar]);
      }
    }
    else{
      
      for(Int_t ibin=0;ibin<=fNVarBins[ivar];ibin++) {
	fVarBinLimits[ibin+fOffset[ivar]] = xbins->At(ibin);
      }
    }
  }   
  
}
//____________________________________________________________________
void AliCFGridSparse::Copy(TObject& c) const
{
  //
  // copy function
  //
  AliCFGridSparse& target = (AliCFGridSparse &) c;

  if(fData)target.fData = fData;
}

//____________________________________________________________________
TH1D* AliCFGridSparse::Slice(Int_t iVar, Double_t *varMin, Double_t *varMax) const
{
  //
  // return a slice (1D-projection) on variable iVar while axis ranges are defined with varMin,varMax
  // arrays varMin and varMax contain the min and max values of each variable.
  // therefore varMin and varMax must have their dimensions equal to fNVar
  //
  
  THnSparse* clone = (THnSparse*)fData->Clone();
  for (Int_t iAxis=0; iAxis<fNVar; iAxis++) {
    clone->GetAxis(iAxis)->SetRangeUser(varMin[iAxis],varMax[iAxis]);
  }
  return clone->Projection(iVar);
}

//____________________________________________________________________
TH2D* AliCFGridSparse::Slice(Int_t iVar1, Int_t iVar2, Double_t *varMin, Double_t *varMax) const
{
  //
  // return a slice (2D-projection) on variables iVar1 and iVar2 while axis ranges are defined with varMin,varMax
  // arrays varMin and varMax contain the min and max values of each variable.
  // therefore varMin and varMax must have their dimensions equal to fNVar
  //
  
  THnSparse* clone = (THnSparse*)fData->Clone();
  for (Int_t iAxis=0; iAxis<fNVar; iAxis++) {
    clone->GetAxis(iAxis)->SetRangeUser(varMin[iAxis],varMax[iAxis]);
  }
  return clone->Projection(iVar1,iVar2);
}

//____________________________________________________________________
TH3D* AliCFGridSparse::Slice(Int_t iVar1, Int_t iVar2, Int_t iVar3, Double_t *varMin, Double_t *varMax) const
{
  //
  // return a slice (3D-projection) on variables iVar1, iVar2 and iVar3 while axis ranges are defined with varMin,varMax
  // arrays varMin and varMax contain the min and max values of each variable.
  // therefore varMin and varMax must have their dimensions equal to fNVar
  //

  THnSparse* clone = (THnSparse*)fData->Clone();
  for (Int_t iAxis=0; iAxis<fNVar; iAxis++) {
    clone->GetAxis(iAxis)->SetRangeUser(varMin[iAxis],varMax[iAxis]);
  }
  return clone->Projection(iVar1,iVar2,iVar3);
}

//____________________________________________________________________
void AliCFGridSparse::SetRangeUser(Int_t iVar, Double_t varMin, Double_t varMax) {
  //
  // set range of axis iVar. 
  //
  fData->GetAxis(iVar)->SetRangeUser(varMin,varMax);
  AliWarning(Form("THnSparse axis %d range has been modified",iVar));
}

//____________________________________________________________________
void AliCFGridSparse::SetRangeUser(Double_t *varMin, Double_t *varMax) {
  //
  // set range of every axis. varMin and varMax must be of dimension fNVar
  //
  for (Int_t iAxis=0; iAxis<fNVar ; iAxis++) { // set new range for every axis
    SetRangeUser(iAxis,varMin[iAxis],varMax[iAxis]);
  }
  AliWarning("THnSparse axes ranges have been modified");
}

//____________________________________________________________________
void AliCFGridSparse::UseAxisRange(Bool_t b) const {
  for (Int_t iAxis=0; iAxis<fNVar; iAxis++) fData->GetAxis(iAxis)->SetBit(TAxis::kAxisRange,b);
}
