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

/* $Id: AliRun.cxx 30859 2009-02-02 16:24:37Z fca $ */

#include <string.h>
#include "AliVectorSparse.h"
#include <TString.h>

/**********************************************************************************************/
/* Sparse vector class, used as row of the AliMatrixSparse class                              */
/*                                                                                            */ 
/* Author: ruben.shahoyan@cern.ch                                                             */
/*                                                                                            */ 
/**********************************************************************************************/

ClassImp(AliVectorSparse)

//___________________________________________________________
AliVectorSparse::AliVectorSparse()
  : fNElems(0),fIndex(0),fElems(0) {}

//___________________________________________________________
AliVectorSparse::AliVectorSparse(const AliVectorSparse& src)
  : TObject(src),fNElems(src.fNElems),fIndex(0),fElems(0)
{
  fIndex = new UShort_t[fNElems];
  fElems = new Double_t[fNElems];
  memcpy(fIndex,src.fIndex,fNElems*sizeof(UShort_t));
  memcpy(fElems,src.fElems,fNElems*sizeof(Double_t));
}

//___________________________________________________________
void AliVectorSparse::Clear(Option_t*)
{
  delete[] fIndex; fIndex = 0;
  delete[] fElems; fElems = 0;
  fNElems = 0;
}

//___________________________________________________________
AliVectorSparse& AliVectorSparse::operator=(const AliVectorSparse& src)
{
  if (&src==this) return *this;
  Clear();
  TObject::operator=(src);
  fNElems = src.fNElems;
  fIndex = new UShort_t[fNElems];
  fElems = new Double_t[fNElems];
  memcpy(fIndex,src.fIndex,fNElems*sizeof(UShort_t));
  memcpy(fElems,src.fElems,fNElems*sizeof(Double_t));
  //
  return *this;
}

//___________________________________________________________
Double_t AliVectorSparse::FindIndex(Int_t ind) const
{
  // return an element with given index
  //printf("V: findindex\n");
  int first = 0;
  int last = fNElems-1;
  while (first<=last) {
    int mid = (first+last)>>1;
    if (ind>fIndex[mid]) first = mid+1;
    else if (ind<fIndex[mid]) last = mid-1;
    else return fElems[mid];
  }
  return 0.0;
}

//___________________________________________________________
void AliVectorSparse::SetToZero(Int_t ind)
{
  // set element to 0 if it was already defined
  int first = 0;
  int last = fNElems-1;
  while (first<=last) {
    int mid = (first+last)>>1;
    if (ind>fIndex[mid]) first = mid+1;
    else if (ind<fIndex[mid]) last = mid-1;
    else {fElems[mid] = 0.; return;}
  }
}

//___________________________________________________________
Double_t& AliVectorSparse::FindIndexAdd(Int_t ind)
{
  // increment an element with given index
  //printf("V: findindexAdd\n");
  int first = 0;
  int last = fNElems-1;
  while (first<=last) {
    int mid = (first+last)>>1;
    if (ind>fIndex[mid]) first = mid+1;
    else if (ind<fIndex[mid]) last = mid-1;
    else return fElems[mid];
  }
  // need to insert a new element
  UShort_t *arrI = new UShort_t[fNElems+1];
  memcpy(arrI,fIndex,first*sizeof(UShort_t));
  arrI[first] = ind;
  memcpy(arrI+first+1,fIndex+first,(fNElems-first)*sizeof(UShort_t));
  delete[] fIndex;
  fIndex = arrI;
  //
  Double_t   *arrE = new Double_t[fNElems+1];
  memcpy(arrE,fElems,first*sizeof(Double_t));
  arrE[first] = 0;
  memcpy(arrE+first+1,fElems+first,(fNElems-first)*sizeof(Double_t));
  delete[] fElems;
  fElems = arrE;
  //
  fNElems++;
  return fElems[first];
  //
}

//__________________________________________________________
void AliVectorSparse::ReSize(Int_t sz,Bool_t copy)
{
  if (sz<1) {Clear(); return;}
    // need to insert a new element
  UShort_t *arrI = new UShort_t[sz];
  Double_t *arrE = new Double_t[sz];
  memset(arrI,0,sz*sizeof(UShort_t));
  memset(arrE,0,sz*sizeof(Double_t));
  //
  if (copy && fIndex) {
    int cpsz = TMath::Min(fNElems,sz);
    memcpy(arrI,fIndex,cpsz*sizeof(UShort_t));
    memcpy(arrE,fElems,cpsz*sizeof(Double_t));
  }
  delete[] fIndex;
  delete[] fElems;
  fIndex = arrI;
  fElems = arrE;
  fNElems = sz;
  //
}

//__________________________________________________________
void AliVectorSparse::SortIndices(Bool_t valuesToo)
{
  // sort indices in increasing order. Used to fix the row after ILUk decomposition
  for (int i=fNElems;i--;) for (int j=i;j--;) if (fIndex[i]<fIndex[j]) { //swap
	UShort_t tmpI = fIndex[i]; fIndex[i] = fIndex[j]; fIndex[j]=tmpI;
	if (valuesToo) {Double_t tmpV = fElems[i];fElems[i]=fElems[j];fElems[j]=tmpV;}
      }
}

//__________________________________________________________
void AliVectorSparse::Print(Option_t* opt)  const
{
  TString sopt = opt; sopt.ToLower();
  int ndig = sopt.Atoi();
  if (ndig<=1) ndig = 2;
  sopt = "%2d:%+.";
  sopt += ndig;
  sopt += "e |";
  printf("|");
  for (int i=0;i<fNElems;i++) printf(sopt.Data(),fIndex[i],fElems[i]);
  printf("\n");
}

//___________________________________________________________
void AliVectorSparse::Add(Double_t *valc,Int_t *indc,Int_t n)
{
  // add indiced array to row. Indices must be in increasing order
  int indx;
  int nadd = 0;
  //
  int last = fNElems-1;
  int mid = 0;
  for (int i=n;i--;) {
    // if the element with this index is already defined, just add the value
    int first = 0;
    Bool_t toAdd = kTRUE;
    indx = indc[i];
    while (first<=last) {
      mid = (first+last)>>1;
      if (indx>fIndex[mid]) first = mid+1;
      else if (indx<fIndex[mid]) last = mid-1;
      else {
	fElems[mid] += valc[i];
	indc[i] = -1;
	toAdd = kFALSE;
	last = mid-1;   // profit from the indices being ordered
	break;
      }
    }
    if (toAdd) nadd++;
  }
  //
  if (nadd<1) return; // nothing to do anymore
  // 
  // need to expand the row
  UShort_t *arrI = new UShort_t[fNElems+nadd];
  Double_t *arrE = new Double_t[fNElems+nadd];
  // copy old elems embedding the new ones
  int inew=0,iold=0;
  for (int i=0;i<n;i++) {  
    if ( (indx=indc[i])<0) continue;
    while (iold<fNElems && fIndex[iold]<indx) {
      arrI[inew]   = fIndex[iold];
      arrE[inew++] = fElems[iold++];
    }
    arrI[inew] = indx;
    arrE[inew++] = valc[i];
  }
  // copy the rest
  while (iold<fNElems) {
    arrI[inew]   = fIndex[iold];
    arrE[inew++] = fElems[iold++];
  }
  //
  delete[] fIndex;
  delete[] fElems;
  fIndex = arrI;
  fElems = arrE;
  //
  fNElems += nadd;
  //
}

