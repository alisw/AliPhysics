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
// AliCFFrame Class                                                 //
// Class to accumulate data on an N-dimensional grid, to be used      //
// as input to get corrections for Reconstruction & Trigger efficiency// 
//                                                                    //
// -- Author : S.Arcelli                                              //
// Still to be done:                                                  //
// --Implement methods to merge cells                                 //
// --Interpolate among bins in a range                                // 
//--------------------------------------------------------------------//
//
//

#include "TSystem.h"
#include <TFile.h>
#include <AliLog.h>
#include "AliCFFrame.h"

//____________________________________________________________________
ClassImp(AliCFFrame)

//____________________________________________________________________
AliCFFrame::AliCFFrame() : 
  TNamed(),
  fNVar(0),
  fNDim(0),
  fNVarBinLimits(0),
  fNVarBins(0x0),
  fIndex(0x0),
  fProduct(0x0),
  fOffset(0x0),
  fVarBinLimits(0x0)
{
  // default constructor
}
//____________________________________________________________________
AliCFFrame::AliCFFrame(const Char_t* name, const Char_t* title) : 
  TNamed(name,title),
  fNVar(0),
  fNDim(0),
  fNVarBinLimits(0),
  fNVarBins(0x0),
  fIndex(0x0),
  fProduct(0x0),
  fOffset(0x0),
  fVarBinLimits(0x0)
{
  // named constructor
}
//____________________________________________________________________
AliCFFrame::AliCFFrame(const Char_t* name, const Char_t* title, const Int_t nVarIn, const Int_t * nBinIn, const Double_t *binLimitsIn) :  
  TNamed(name,title),
  fNVar(0),
  fNDim(0),
  fNVarBinLimits(0),
  fNVarBins(0x0),
  fIndex(0x0),
  fProduct(0x0),
  fOffset(0x0),
  fVarBinLimits(0x0)
{
  //
  // main constructor
  //

  //the number of variables on the grid
  fNVar=nVarIn;

  // the binning in each variable
  fNVarBins= new Int_t[fNVar];

  fIndex= new Int_t[fNVar];

  Int_t ndimTot=1;
  Int_t nbinTot=0;


  //Calculate total number of elements and bins...

  for(Int_t i=0;i<fNVar;i++){
    fNVarBins[i]=nBinIn[i];
    ndimTot*=nBinIn[i];
    nbinTot+=(nBinIn[i]+1);
  } 


  //total number of elements

  fNDim=ndimTot;

  fOffset= new Int_t[fNVar];
  fProduct= new Int_t[fNVar];

  for(Int_t ivar=0;ivar<fNVar;ivar++){
    Int_t offset=0;
    for(Int_t i =0;i<ivar;i++)offset+=(fNVarBins[i]+1);      
    fOffset[ivar]=offset;
    Int_t prod=1;
    for(Int_t i=0;i<ivar;i++)prod*=fNVarBins[i];
    fProduct[ivar]=prod;
  }

  //The bin limits

  fNVarBinLimits=nbinTot;
  fVarBinLimits=new Double_t[fNVarBinLimits];
  if(binLimitsIn){
    for(Int_t i=0;i<fNVarBinLimits;i++){
      fVarBinLimits[i] =binLimitsIn[i];
    } 
  }
}
//____________________________________________________________________
AliCFFrame::AliCFFrame(const AliCFFrame& c) : 
  TNamed(),
  fNVar(0),
  fNDim(0),
  fNVarBinLimits(0),
  fNVarBins(0x0),
  fIndex(0x0),
  fProduct(0x0),
  fOffset(0x0),
  fVarBinLimits(0x0)
{
  //
  // copy constructor
  //
  ((AliCFFrame &)c).Copy(*this);
}
//____________________________________________________________________
AliCFFrame::~AliCFFrame()
{
  //
  // destructor
  //
  if(fNVarBins)delete [] fNVarBins;
  if(fVarBinLimits)delete [] fVarBinLimits;
  if(fIndex)delete [] fIndex;
  if(fProduct)delete [] fProduct;
  if(fOffset)delete [] fOffset;

}
//____________________________________________________________________
AliCFFrame &AliCFFrame::operator=(const AliCFFrame &c)
{
  //
  // assigment operator
  //
  if (this != &c)
    ((AliCFFrame &) c).Copy(*this);
  return *this;
} 
//____________________________________________________________________
void AliCFFrame::SetBinLimits(Int_t ivar, Double_t *array)
{
  //
  // setting the arrays containing the bin limits 
  //
  Int_t nbins=fNVarBins[ivar]+1;
  for(Int_t i=0;i<nbins;i++){
    fVarBinLimits[fOffset[ivar]+i] =array[i];
  } 
} 
//____________________________________________________________________
void AliCFFrame::GetBinLimits(Int_t ivar, Double_t *array) const
{
  //
  // getting the arrays containing the bin limits 
  //
  Int_t nbins=fNVarBins[ivar]+1;
  for(Int_t i=0;i<nbins;i++){
    array[i]=fVarBinLimits[fOffset[ivar]+i];
  } 
} 
//____________________________________________________________________
Int_t AliCFFrame::GetBinIndex(Int_t *ibin) const
{
  //
  // getting the element number on the grid for a given set of bin indeces 
  //
  Int_t ind=ibin[fNVar-1];
  for(Int_t i=fNVar-1;i>0;i--){
    ind=ibin[i-1]+ind*fNVarBins[i-1]; 
  }
  return ind;
} 
//____________________________________________________________________
Int_t AliCFFrame::GetBinIndex(Int_t ivar, Int_t ind) const
{
  //
  // getting the bin index on a given dim. ivar for a given element ind on the grid
  //
  Int_t index=-1;
  for(Int_t i=fNVar-1;i>=0;i--){
    index=ind/fProduct[i];
    fIndex[i]=index;
    ind-=index*fProduct[i];
  }
  
  index=fIndex[ivar];
  return index;
}
//____________________________________________________________________
void AliCFFrame::GetBinIndex(Int_t ind, Int_t *bins ) const
{
  //
  // getting the set of bin indeces for a given element ind on the grid
  //
  Int_t index=-1;
  for(Int_t i=fNVar-1;i>=0;i--){
    index=ind/fProduct[i];
    bins[i]=index;
    ind-=index*fProduct[i];
  }
  return;
}
//____________________________________________________________________
Double_t AliCFFrame::GetBinCenter(Int_t ivar, Int_t ibin) const
{
  //
  // getting the bin center of a given bin ibin along variable ivar
  //
 
  Double_t binMin=fVarBinLimits[fOffset[ivar]+ibin];
  Double_t binMax=fVarBinLimits[fOffset[ivar]+ibin+1];
  Double_t val=0.5*(binMin+binMax);
  return val;
}
//____________________________________________________________________
void AliCFFrame::GetBinCenters(Int_t *ibin, Double_t *binCenter) const
{
  //
  // gives the centers of the N-dim bin identified by a set of bin indeces 
  // invar
  //
  for(Int_t i=0;i<fNVar;i++){
    binCenter[i]=GetBinCenter(i,ibin[i]);
  }
} //____________________________________________________________________
Double_t AliCFFrame::GetBinSize(Int_t ivar, Int_t ibin) const
{
  //
  // getting the bin size on axis ivar of a given bin ibin 
  //

  if(ibin>=fNVarBins[ivar]){
    AliWarning(Form("bin index out of range, number of bins in variable %i is %i", ivar,fNVarBins[ivar]));
    return -1;
  } 

  Double_t binMin=fVarBinLimits[fOffset[ivar]+ibin];
  Double_t binMax=fVarBinLimits[fOffset[ivar]+ibin+1];
  Double_t val=binMax-binMin;
  return val;
}
//____________________________________________________________________
void AliCFFrame::GetBinSizes(Int_t *ibin, Double_t *binSizes) const
{
  //
  // gives the sizes of the N-dim bin identified by a set of bin indeces 
  // ibin
  //
  for(Int_t i=0;i<fNVar;i++){
    binSizes[i]=GetBinSize(i,ibin[i]);
  }
} 

//____________________________________________________________________
void AliCFFrame::Copy(TObject& c) const
{
  //
  // copy function
  //
  AliCFFrame& target = (AliCFFrame &) c;

  target.fNVar=fNVar;
  target.fNDim=fNDim;
  target.fNVarBinLimits=fNVarBinLimits;
  if (fNVarBins)
    target.fNVarBins = fNVarBins;
  if (fVarBinLimits)
    target.fVarBinLimits = fVarBinLimits;
  if (fProduct)
    target.fProduct = fProduct;
  if (fOffset)
    target.fOffset = fOffset;  
}
//____________________________________________________________________
void AliCFFrame::PrintBinLimits()
{
  //
  // printing the array containing the bin limits 
  //
  for(Int_t i=0;i<fNVarBinLimits;i++){
    AliInfo(Form("bin limit index %i is: %f",i, fVarBinLimits[i]));
  } 
} 
//____________________________________________________________________
void AliCFFrame::PrintNBins()
{
  //
  // printing the array containing the # of bins  
  //
  for(Int_t i=0;i<fNVar;i++){
    AliInfo(Form("bins in var %i are: %i",i, fNVarBins[i]));
  } 
} 
//____________________________________________________________________
void AliCFFrame::Save(const Char_t *outfile) const
{
  //
  // Save the grid to a root file
  //

  const char *dirname = "./";
  TString filename = outfile;
  TFile *file=0x0;
  if((gSystem->FindFile(dirname,filename))!=NULL){
      file = new TFile( outfile,"UPDATE");
  }
  else{
    file = new TFile( outfile,"RECREATE");
  } 
  file->cd();
  //write the object to a file
  this->Write(GetName(),TObject::kSingleKey);
  file->Close();
  delete file;
}
