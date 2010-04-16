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
TH1D *AliCFGridSparse::Project(Int_t ivar) const
{
  //
  // Make a 1D projection along variable ivar 
  //

  TH1D *hist=fData->Projection(ivar);
  hist->SetXTitle(GetVarTitle(ivar));
  hist->SetName(Form("%s_proj-%s",GetName(),GetVarTitle(ivar)));
  hist->SetTitle(Form("%s: projection on %s",GetTitle(),GetVarTitle(ivar)));
  for (Int_t iBin=1; iBin<=GetNBins(ivar); iBin++) {
    TString binLabel = GetAxis(ivar)->GetBinLabel(iBin) ;
    if (binLabel.CompareTo("") != 0) hist->GetXaxis()->SetBinLabel(iBin,binLabel);
  }
  return hist;
}
//___________________________________________________________________
TH2D *AliCFGridSparse::Project(Int_t ivar1, Int_t ivar2) const
{
  //
  // Make a 2D projection along variables ivar1 & ivar2 
  //

  TH2D *hist=fData->Projection(ivar2,ivar1); //notice inverted axis (THnSparse uses TH3 2d-projection convention...)
  hist->SetXTitle(GetVarTitle(ivar1));
  hist->SetYTitle(GetVarTitle(ivar2));
  hist->SetName(Form("%s_proj-%s,%s",GetName(),GetVarTitle(ivar1),GetVarTitle(ivar2)));
  hist->SetTitle(Form("%s: projection on %s-%s",GetTitle(),GetVarTitle(ivar1),GetVarTitle(ivar2)));
  for (Int_t iBin=1; iBin<=GetNBins(ivar1); iBin++) {
    TString binLabel = GetAxis(ivar1)->GetBinLabel(iBin) ;
    if (binLabel.CompareTo("") != 0) hist->GetXaxis()->SetBinLabel(iBin,binLabel);
  }
  for (Int_t iBin=1; iBin<=GetNBins(ivar2); iBin++) {
    TString binLabel = GetAxis(ivar2)->GetBinLabel(iBin) ;
    if (binLabel.CompareTo("") != 0) hist->GetYaxis()->SetBinLabel(iBin,binLabel);
  }
  return hist;
}
//___________________________________________________________________
TH3D *AliCFGridSparse::Project(Int_t ivar1, Int_t ivar2, Int_t ivar3) const
{
  //
  // Make a 3D projection along variables ivar1 & ivar2 & ivar3 
  //

  TH3D *hist=fData->Projection(ivar1,ivar2,ivar3); 
  hist->SetXTitle(GetVarTitle(ivar1));
  hist->SetYTitle(GetVarTitle(ivar2));
  hist->SetZTitle(GetVarTitle(ivar3));
  hist->SetName(Form("%s_proj-%s,%s,%s",GetName(),GetVarTitle(ivar1),GetVarTitle(ivar2),GetVarTitle(ivar3)));
  hist->SetTitle(Form("%s: projection on %s-%s-%s",GetTitle(),GetVarTitle(ivar1),GetVarTitle(ivar2),GetVarTitle(ivar3)));
  for (Int_t iBin=1; iBin<=GetNBins(ivar1); iBin++) {
    TString binLabel = GetAxis(ivar1)->GetBinLabel(iBin) ;
    if (binLabel.CompareTo("") != 0) hist->GetXaxis()->SetBinLabel(iBin,binLabel);
  }
  for (Int_t iBin=1; iBin<=GetNBins(ivar2); iBin++) {
    TString binLabel = GetAxis(ivar2)->GetBinLabel(iBin) ;
    if (binLabel.CompareTo("") != 0) hist->GetYaxis()->SetBinLabel(iBin,binLabel);
  }
  for (Int_t iBin=1; iBin<=GetNBins(ivar3); iBin++) {
    TString binLabel = GetAxis(ivar3)->GetBinLabel(iBin) ;
    if (binLabel.CompareTo("") != 0) hist->GetZaxis()->SetBinLabel(iBin,binLabel);
  }
  return hist;
}

//___________________________________________________________________
AliCFGridSparse* AliCFGridSparse::Project(Int_t nVars, const Int_t* vars, const Double_t* varMin, const Double_t* varMax, Bool_t useBins) const
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
      if (useBins)  clone->GetAxis(iAxis)->SetRange((Int_t)varMin[iAxis],(Int_t)varMax[iAxis]);
      else          clone->GetAxis(iAxis)->SetRangeUser(varMin[iAxis],varMax[iAxis]);
    }
  }
  else AliInfo("Keeping same axis ranges");

  out->SetGrid(clone->Projection(nVars,vars));
  delete [] bins;
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

//_____________________________________________________________________
// Int_t AliCFGridSparse::GetEmptyBins( Double_t *varMin, Double_t* varMax ) const 
// {
//   //
//   // Get empty bins in a range specified by varMin and varMax
//   //

//   Int_t *indexMin = new Int_t[GetNVar()];
//   Int_t *indexMax = new Int_t[GetNVar()];

//   //Find out the min and max bins

//   for (Int_t i=0;i<GetNVar();i++) {
//     Double_t xmin=varMin[i]; // the min values  
//     Double_t xmax=varMax[i]; // the min values  
//     Int_t nbins = fData->GetAxis(i)->GetNbins()+1;
//     Double_t *bins=new Double_t[nbins];
//     for(Int_t ibin =0;ibin<nbins;ibin++){
//      bins[ibin] = fVarBinLimits[ibin+fOffset[i]];
//     }
//     indexMin[i] = TMath::BinarySearch(nbins,bins,xmin);
//     indexMax[i] = TMath::BinarySearch(nbins,bins,xmax);
//     delete [] bins;
//   }

//   Int_t val=0;
//   for(Int_t i=0;i<fNDim;i++){
//     for (Int_t j=0;j<GetNVar();j++)fIndex[j]=GetBinIndex(j,i);
//     Bool_t isIn=kTRUE;
//     for (Int_t j=0;j<GetNVar();j++){
//       if(!(fIndex[j]>=indexMin[j] && fIndex[j]<=indexMax[j]))isIn=kFALSE;   
//     }
//     if(isIn && GetElement(i)<=0)val++;     
//   }
//   AliInfo(Form(" the empty bins = %i ",val)); 

//   delete [] indexMin;
//   delete [] indexMax;
//   return val;
// } 

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

//_____________________________________________________________________
// Double_t AliCFGridSparse::GetIntegral(Int_t *binMin, Int_t* binMax ) const 
// {
//   //
//   // Get Integral in a range of bin indeces (extremes included)
//   //

//   Double_t val=0;

//   for(Int_t i=0;i<GetNVar();i++){
//     if(binMin[i]<1)binMin[i]=1;
//     if(binMax[i]>fNVarBins[i])binMax[i]=fNVarBins[i];
//     if((binMin[i]>binMax[i])){
//       AliInfo(Form(" Bin indices in variable %i in reverse order, please check!", i));
//       return val;
//     }
//   }
//   val=GetSum(0,binMin,binMax);

//   return val;
// } 

//_____________________________________________________________________
// Double_t AliCFGridSparse::GetIntegral(Double_t *varMin, Double_t* varMax ) const 
// {
//   //
//   // Get Integral in a range (extremes included)
//   //

//   Int_t *indexMin=new Int_t[GetNVar()];
//   Int_t *indexMax=new Int_t[GetNVar()];

//   //Find out the min and max bins

//   for(Int_t i=0;i<GetNVar();i++){
//     Double_t xmin=varMin[i]; // the min values  
//     Double_t xmax=varMax[i]; // the min values  
//     Int_t nbins=fNVarBins[i]+1;
//     Double_t *bins=new Double_t[nbins];
//     for(Int_t ibin =0;ibin<nbins;ibin++){
//      bins[ibin] = fVarBinLimits[ibin+fOffset[i]];
//     }
//     indexMin[i] = TMath::BinarySearch(nbins,bins,xmin);
//     indexMax[i] = TMath::BinarySearch(nbins,bins,xmax);
//     delete [] bins;
//   }

//   //move to the TH/THnSparse convention in N-dim bin numbering 
//   for(Int_t i=0;i<GetNVar();i++){
//     indexMin[i]+=1;
//     indexMax[i]+=1;
//   }

//   Double_t val=GetIntegral(indexMin,indexMax);

//   delete [] indexMin;
//   delete [] indexMax;

//   return val;
// } 

// //_____________________________________________________________________
// Double_t AliCFGridSparse::GetSum(Int_t ivar, Int_t *binMin, Int_t* binMax) const 
// {
//   //
//   // recursively add over nested loops.... 
//   //

//   static Double_t val;
//   if (ivar==0) val=0.;
//   for(Int_t ibin=binMin[ivar]-1 ; ibin<=binMax[ivar]-1 ; ibin++){
//     //-1 is to move from TH/ThnSparse N-dim bin convention to one in AliCFFrame
//     fIndex[ivar]=ibin;
//     if(ivar<GetNVar()-1) {
//       val=GetSum(ivar+1,binMin,binMax);
//     }
//     else {
//       Int_t iel=GetBinIndex(fIndex);
//       val+=GetElement(iel);
//     }
//   }
  
//   return val;
// }

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
TH1D* AliCFGridSparse::Slice(Int_t iVar, const Double_t *varMin, const Double_t *varMax, Bool_t useBins) const
{
  //
  // return a slice (1D-projection) on variable iVar while axis ranges are defined with varMin,varMax
  // arrays varMin and varMax contain the min and max values of each variable.
  // therefore varMin and varMax must have their dimensions equal to GetNVar()
  // If useBins=true, varMin and varMax are taken as bin numbers
  
  THnSparse* clone = (THnSparse*)fData->Clone();
  for (Int_t iAxis=0; iAxis<GetNVar(); iAxis++) {
    if (useBins)  clone->GetAxis(iAxis)->SetRange((Int_t)varMin[iAxis],(Int_t)varMax[iAxis]);
    else          clone->GetAxis(iAxis)->SetRangeUser(varMin[iAxis],varMax[iAxis]);
  }

  TH1D* projection = 0x0 ;

  TObject* objInMem = gROOT->FindObject(Form("%s_proj_%d",clone->GetName(),iVar)) ;
  if (objInMem) {
    TString name(objInMem->ClassName());
    if (name.CompareTo("TH1D")==0) projection = (TH1D*) objInMem ;
    else projection = clone->Projection(iVar);
  }
  else projection = clone->Projection(iVar); 
  delete clone;
  for (Int_t iBin=1; iBin<=GetNBins(iVar); iBin++) {
    TString binLabel = GetAxis(iVar)->GetBinLabel(iBin) ;
    if (binLabel.CompareTo("") != 0) projection->GetXaxis()->SetBinLabel(iBin,binLabel);
  }
  return projection ;
}

//____________________________________________________________________
TH2D* AliCFGridSparse::Slice(Int_t iVar1, Int_t iVar2, const Double_t *varMin, const Double_t *varMax, Bool_t useBins) const
{
  //
  // return a slice (2D-projection) on variables iVar1 and iVar2 while axis ranges are defined with varMin,varMax
  // arrays varMin and varMax contain the min and max values of each variable.
  // therefore varMin and varMax must have their dimensions equal to GetNVar()
  // If useBins=true, varMin and varMax are taken as bin numbers
  
  THnSparse* clone = (THnSparse*)fData->Clone();
  for (Int_t iAxis=0; iAxis<GetNVar(); iAxis++) {
    if (useBins)  clone->GetAxis(iAxis)->SetRange((Int_t)varMin[iAxis],(Int_t)varMax[iAxis]);
    else          clone->GetAxis(iAxis)->SetRangeUser(varMin[iAxis],varMax[iAxis]);
  }

  TH2D* projection = 0x0 ;

  TObject* objInMem = gROOT->FindObject(Form("%s_proj_%d_%d",clone->GetName(),iVar2,iVar1)) ;
  if (objInMem) {
    TString name(objInMem->ClassName());
    if (name.CompareTo("TH2D")==0) projection = (TH2D*) objInMem ;
    else projection = clone->Projection(iVar1,iVar2);
  }
  else projection = clone->Projection(iVar1,iVar2); 
  delete clone;
  for (Int_t iBin=1; iBin<=GetNBins(iVar1); iBin++) {
    TString binLabel = GetAxis(iVar1)->GetBinLabel(iBin) ;
    if (binLabel.CompareTo("") != 0) projection->GetXaxis()->SetBinLabel(iBin,binLabel);
  }
  for (Int_t iBin=1; iBin<=GetNBins(iVar2); iBin++) {
    TString binLabel = GetAxis(iVar2)->GetBinLabel(iBin) ;
    if (binLabel.CompareTo("") != 0) projection->GetYaxis()->SetBinLabel(iBin,binLabel);
  }
  return projection ;
}

//____________________________________________________________________
TH3D* AliCFGridSparse::Slice(Int_t iVar1, Int_t iVar2, Int_t iVar3, const Double_t *varMin, const Double_t *varMax, Bool_t useBins) const
{
  //
  // return a slice (3D-projection) on variables iVar1, iVar2 and iVar3 while axis ranges are defined with varMin,varMax
  // arrays varMin and varMax contain the min and max values of each variable.
  // therefore varMin and varMax must have their dimensions equal to GetNVar()
  // If useBins=true, varMin and varMax are taken as bin numbers

  THnSparse* clone = (THnSparse*)fData->Clone();
  for (Int_t iAxis=0; iAxis<GetNVar(); iAxis++) {
    if (useBins)  clone->GetAxis(iAxis)->SetRange((Int_t)varMin[iAxis],(Int_t)varMax[iAxis]);
    else          clone->GetAxis(iAxis)->SetRangeUser(varMin[iAxis],varMax[iAxis]);
  }

  TH3D* projection = 0x0 ;

  TObject* objInMem = gROOT->FindObject(Form("%s_proj_%d_%d_%d",clone->GetName(),iVar1,iVar2,iVar3)) ;
  if (objInMem) {
    TString name(objInMem->ClassName());
    if (name.CompareTo("TH3D")==0) projection = (TH3D*) objInMem ;
    else projection = clone->Projection(iVar1,iVar2,iVar3);
  }
  else projection = clone->Projection(iVar1,iVar2,iVar3); 
  delete clone;
  for (Int_t iBin=1; iBin<=GetNBins(iVar1); iBin++) {
    TString binLabel = GetAxis(iVar1)->GetBinLabel(iBin) ;
    if (binLabel.CompareTo("") != 0) projection->GetXaxis()->SetBinLabel(iBin,binLabel);
  }
  for (Int_t iBin=1; iBin<=GetNBins(iVar2); iBin++) {
    TString binLabel = GetAxis(iVar2)->GetBinLabel(iBin) ;
    if (binLabel.CompareTo("") != 0) projection->GetYaxis()->SetBinLabel(iBin,binLabel);
  }
  for (Int_t iBin=1; iBin<=GetNBins(iVar3); iBin++) {
    TString binLabel = GetAxis(iVar3)->GetBinLabel(iBin) ;
    if (binLabel.CompareTo("") != 0) projection->GetZaxis()->SetBinLabel(iBin,binLabel);
  }
  return projection ;
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
void AliCFGridSparse::SetRangeUser(const Double_t *varMin, const Double_t *varMax) {
  //
  // set range of every axis. varMin and varMax must be of dimension GetNVar()
  //
  for (Int_t iAxis=0; iAxis<GetNVar() ; iAxis++) { // set new range for every axis
    SetRangeUser(iAxis,varMin[iAxis],varMax[iAxis]);
  }
  AliWarning("THnSparse axes ranges have been modified");
}

//____________________________________________________________________
void AliCFGridSparse::UseAxisRange(Bool_t b) const {
  for (Int_t iAxis=0; iAxis<GetNVar(); iAxis++) fData->GetAxis(iAxis)->SetBit(TAxis::kAxisRange,b);
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

