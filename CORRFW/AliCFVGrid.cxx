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
// AliCFVGrid Class                                                   //
// Just an interface to handle both AliCFGrid and AliCFGridSparse     //
// implementations                                                    //
//                                                                    //
// -- Author : S.Arcelli                                              //
//--------------------------------------------------------------------//
//
//
#include "AliCFVGrid.h"
#include "TMath.h"
#include <AliLog.h>

//____________________________________________________________________
ClassImp(AliCFVGrid)

//____________________________________________________________________
AliCFVGrid::AliCFVGrid() : 
  AliCFFrame(),
  fSumW2(kFALSE)
{
  // default constructor
}
//____________________________________________________________________
AliCFVGrid::AliCFVGrid(const Char_t* name, const Char_t* title) : 
  AliCFFrame(name,title),
  fSumW2(kFALSE)
{
  // default constructor
}

//____________________________________________________________________
AliCFVGrid::AliCFVGrid(const Char_t* name, const Char_t* title, const Int_t nVarIn, const Int_t * nBinIn, const Double_t *binLimitsIn) :  
  AliCFFrame(name,title,nVarIn,nBinIn,binLimitsIn),
  fSumW2(kFALSE)
{
  //
  // main constructor
  //

}

//____________________________________________________________________
AliCFVGrid::AliCFVGrid(const AliCFVGrid& c) : 
  AliCFFrame(c),
  fSumW2(c.fSumW2)
{
  //
  // copy constructor
  //
 ((AliCFVGrid &)c).Copy(*this);

}
//____________________________________________________________________
AliCFVGrid::~AliCFVGrid()
{
  //
  // destructor
  //
}

//____________________________________________________________________
AliCFVGrid &AliCFVGrid::operator=(const AliCFVGrid &c)
{
  //
  // assigment operator
  //
  if (this != &c)
    AliCFFrame::operator=(c);
  return *this;
} 

//____________________________________________________________________
void AliCFVGrid::Copy(TObject& c) const
{
  //
  // copy function
  //
  AliCFVGrid& target = (AliCFVGrid &) c;

  target.fSumW2=fSumW2;
}

//____________________________________________________________________
void AliCFVGrid::Scale(Int_t index, Double_t *fact)
{
  //
  //scale content of a certain cell by (positive) fact (with error)
  //
  if(GetElement(index)==0 || fact[0]==0)return;

  Double_t in[2], out[2];
  in[0]=GetElement(index);
  in[1]=GetElementError(index);
  GetScaledValues(fact,in,out);
  SetElement(index,out[0]);
  if(fSumW2)SetElementError(index,out[1]);
}
//____________________________________________________________________
void AliCFVGrid::Scale(Int_t *bin, Double_t *fact)
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
void AliCFVGrid::Scale(Double_t *var, Double_t *fact) 
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
void AliCFVGrid::Scale( Double_t *fact) 
{
  //
  //scale contents of the whole grid by fact
  //

  for(Int_t iel=0;iel<fNDim;iel++){
    Scale(iel,fact);
  }
}

//____________________________________________________________________
Int_t AliCFVGrid::GetEmptyBins() const {
  //
  // Get empty bins 
  //
  Int_t val=0;
  for(Int_t i=0;i<fNDim;i++){
    if(GetElement(i)<=0)val++;     
  }
  return val;
} 
//_____________________________________________________________________
Int_t AliCFVGrid::GetEmptyBins( Double_t *varMin, Double_t* varMax ) const 
{
  //
  // Get empty bins in a range
  //

  Int_t *indexMin=new Int_t[fNVar];
  Int_t *indexMax=new Int_t[fNVar];

  //Find out the min and max bins

  for(Int_t i=0;i<fNVar;i++){
    Double_t xmin=varMin[i]; // the min values  
    Double_t xmax=varMax[i]; // the min values  
    Int_t nbins=fNVarBins[i]+1;
    Double_t *bins=new Double_t[nbins];
    for(Int_t ibin =0;ibin<nbins;ibin++){
     bins[ibin] = fVarBinLimits[ibin+fOffset[i]];
    }
    indexMin[i] = TMath::BinarySearch(nbins,bins,xmin);
    indexMax[i] = TMath::BinarySearch(nbins,bins,xmax);
    delete [] bins;
  }

  Int_t val=0;
  for(Int_t i=0;i<fNDim;i++){
    for (Int_t j=0;j<fNVar;j++)fIndex[j]=GetBinIndex(j,i);
    Bool_t isIn=kTRUE;
    for (Int_t j=0;j<fNVar;j++){
      if(!(fIndex[j]>=indexMin[j] && fIndex[j]<=indexMax[j]))isIn=kFALSE;   
    }
    if(isIn && GetElement(i)<=0)val++;     
  }
  AliInfo(Form(" the empty bins = %i ",val)); 

  delete [] indexMin;
  delete [] indexMax;
  return val;
} 
//____________________________________________________________________
Int_t AliCFVGrid::CheckStats(Double_t thr) const
{
  //
  // Count the cells below a certain threshold
  //
  Int_t ncellsLow=0;
  for(Int_t i=0;i<fNDim;i++){
    if(GetElement(i)<thr)ncellsLow++;
  }
  return ncellsLow;
}
//_____________________________________________________________________
Double_t AliCFVGrid::GetIntegral() const 
{
  //
  // Get full Integral
  //
  Double_t val=0;
  for(Int_t i=0;i<fNDim;i++){
    val+=GetElement(i);     
  }
  return val;  
} 
//_____________________________________________________________________
Double_t AliCFVGrid::GetIntegral(Int_t *binMin, Int_t* binMax ) const 
{
  //
  // Get Integral in a range of bin indeces (extremes included)
  //

  Double_t val=0;

  for(Int_t i=0;i<fNVar;i++){
    if(binMin[i]<1)binMin[i]=1;
    if(binMax[i]>fNVarBins[i])binMax[i]=fNVarBins[i];
    if((binMin[i]>binMax[i])){
      AliInfo(Form(" Bin indeces in variable %i in reverse order, please check!", i));
      return val;
    }
  }
  val=GetSum(0,binMin,binMax);

  return val;

} 
//_____________________________________________________________________
Double_t AliCFVGrid::GetIntegral(Double_t *varMin, Double_t* varMax ) const 
{
  //
  // Get Integral in a range (extremes included)
  //

  Int_t *indexMin=new Int_t[fNVar];
  Int_t *indexMax=new Int_t[fNVar];

  //Find out the min and max bins

  for(Int_t i=0;i<fNVar;i++){
    Double_t xmin=varMin[i]; // the min values  
    Double_t xmax=varMax[i]; // the min values  
    Int_t nbins=fNVarBins[i]+1;
    Double_t *bins=new Double_t[nbins];
    for(Int_t ibin =0;ibin<nbins;ibin++){
     bins[ibin] = fVarBinLimits[ibin+fOffset[i]];
    }
    indexMin[i] = TMath::BinarySearch(nbins,bins,xmin);
    indexMax[i] = TMath::BinarySearch(nbins,bins,xmax);
    delete [] bins;
  }

  //move to the TH/THnSparse convention in N-dim bin numbering 
  for(Int_t i=0;i<fNVar;i++){
    indexMin[i]+=1;
    indexMax[i]+=1;
  }

  Double_t val=GetIntegral(indexMin,indexMax);

  delete [] indexMin;
  delete [] indexMax;

  return val;
} 

//_____________________________________________________________________
Double_t AliCFVGrid::GetSum(Int_t ivar, Int_t *binMin, Int_t* binMax) const 
{
  //
  // recursively add over nested loops.... 
  //
  static Double_t val;
  if(ivar==0)val=0.;
  for(Int_t ibin=binMin[ivar]-1;ibin<=binMax[ivar]-1;ibin++){
  //-1 is to move from TH/ThnSparse N-dim bin convention to one in AliCFFrame
    fIndex[ivar]=ibin;
    if(ivar<fNVar-1) {
      val=GetSum(ivar+1,binMin,binMax);
    }
    else {
      Int_t iel=GetBinIndex(fIndex);
      val+=GetElement(iel);
    }
  }

  return val;
}

//____________________________________________________________________
Long64_t AliCFVGrid::Merge(TCollection* list)
{
  // Merge a list of AliCF grids with this (needed for
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
    AliCFVGrid* entry = dynamic_cast<AliCFVGrid*> (obj);
    if (entry == 0) 
      continue;
    this->Add(entry);
    count++;
  }

  return count+1;
}

//____________________________________________________________________
void AliCFVGrid::GetScaledValues(Double_t *fact, Double_t *in, Double_t *out) const{
  //
  // scale input *in and it error by (positive) fact (with error)
  // and erite it to *out
  //
    out[0]=in[0]*fact[0];
    out[1]=TMath::Sqrt(in[1]*in[1]/in[0]/in[0]
		       +fact[1]*fact[1]/fact[0]/fact[0])*out[0];
    
}
