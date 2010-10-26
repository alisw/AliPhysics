/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------
//           Implementation of the Event Time class
//           for the Event Data Summary class
//           This class contains the Event Time
//           as estimated by the TOF combinatorial algorithm
// Origin: A.De Caro, decaro@lsa.infn.it
//-----------------------------------------------------------------

//---- standard headers ----
#include "Riostream.h"
//---- Root headers --------
#include "TArrayF.h"
#include "TArrayI.h"
//---- AliRoot headers -----
#include "AliTOFHeader.h"


ClassImp(AliTOFHeader)

//--------------------------------------------------------------------------
AliTOFHeader::AliTOFHeader() :
  TObject(),
  fDefaultEventTimeValue(0.),
  fDefaultEventTimeRes(0.),
  fNbins(0),
  fEventTimeValues(new TArrayF(1)),
  fEventTimeRes(new TArrayF(1)),
  fNvalues(new TArrayI(1)),
  fTOFtimeResolution(0.),
  fT0spread(0.)
{
  //
  // Default Constructor
  //

}

//--------------------------------------------------------------------------
AliTOFHeader::AliTOFHeader(Float_t defEvTime, Float_t defResEvTime,
			   Int_t nDifPbins, Float_t *times, Float_t *res,
			   Int_t *nPbin, Float_t tofTimeRes, Float_t t0spread) :
  TObject(),
  fDefaultEventTimeValue(defEvTime),
  fDefaultEventTimeRes(defResEvTime),
  fNbins(nDifPbins),
  fEventTimeValues(new TArrayF(nDifPbins)),
  fEventTimeRes(new TArrayF(nDifPbins)),
  fNvalues(new TArrayI(nDifPbins)),
  fTOFtimeResolution(tofTimeRes),
  fT0spread(t0spread)
{
  //
  // Constructor for vertex Z from pixels
  //

  for (Int_t ii=0; ii<fNbins; ii++) {
    fEventTimeValues->SetAt(times[ii],ii);
    fEventTimeRes->SetAt(res[ii],ii);
    fNvalues->SetAt(nPbin[ii],ii);
  }

}

//--------------------------------------------------------------------------
AliTOFHeader::AliTOFHeader(const AliTOFHeader &source):
  TObject(source),
  fDefaultEventTimeValue(source.fDefaultEventTimeValue),
  fDefaultEventTimeRes(source.fDefaultEventTimeRes),
  fNbins(source.fNbins),
  fEventTimeValues(new TArrayF(fNbins)),
  fEventTimeRes(new TArrayF(fNbins)),
  fNvalues(new TArrayI(fNbins)),
  fTOFtimeResolution(source.fTOFtimeResolution),
  fT0spread(source.fT0spread)
{
  //
  // Copy constructor
  //


  for(Int_t i=0;i<fNbins;i++) {
    fEventTimeValues->SetAt(source.fEventTimeValues->GetAt(i),i);
    fEventTimeRes->SetAt(source.fEventTimeRes->GetAt(i),i);
    fNvalues->SetAt(source.fNvalues->GetAt(i),i);
  }

}
//--------------------------------------------------------------------------
AliTOFHeader &AliTOFHeader::operator=(const AliTOFHeader &source){
  //
  // assignment operator
  //
  if(&source != this){
    TObject::operator=(source);

    fDefaultEventTimeValue=source.fDefaultEventTimeValue;
    fDefaultEventTimeRes=source.fDefaultEventTimeRes;
    fNbins=source.fNbins;
    fEventTimeValues=new TArrayF(fNbins);
    fEventTimeRes=new TArrayF(fNbins);
    fNvalues=new TArrayI(fNbins);
    fTOFtimeResolution=source.fTOFtimeResolution;
    fT0spread=source.fT0spread;
    for(Int_t i=0;i<fNbins;i++) {
      fEventTimeValues->SetAt(source.fEventTimeValues->GetAt(i),i);
      fEventTimeRes->SetAt(source.fEventTimeRes->GetAt(i),i);
      fNvalues->SetAt(source.fNvalues->GetAt(i),i);
    }
  }
  return *this;
}
//--------------------------------------------------------------------------
void AliTOFHeader::Copy(TObject &obj) const {
  
  // this overwrites the virtual TOBject::Copy()
  // to allow run time copying without casting
  // in AliESDEvent

  if(this==&obj)return;
  AliTOFHeader *robj = dynamic_cast<AliTOFHeader*>(&obj);
  if(!robj)return; // not an AliTOFHeader
  *robj = *this;

}

//--------------------------------------------------------------------------
AliTOFHeader::~AliTOFHeader()
{

  fEventTimeValues->Reset();
  fEventTimeRes->Reset();
  fNvalues->Reset();

  delete fEventTimeValues;
  delete fEventTimeRes;
  delete fNvalues;

}

//--------------------------------------------------------------------------
void AliTOFHeader::SetNbins(Int_t nbins)
{
  //
  //
  //

  fNbins=nbins;
  fEventTimeValues->Set(nbins);
  fEventTimeRes->Set(nbins);
  fNvalues->Set(nbins);

}
//--------------------------------------------------------------------------
void AliTOFHeader::SetEventTimeValues(TArrayF *arr)
{
  //
  //
  //

  fNbins=arr->GetSize();
  fEventTimeValues->Set(arr->GetSize());
  for (Int_t ii=0; ii<arr->GetSize(); ii++)
    fEventTimeValues->SetAt(arr->GetAt(ii),ii);

}
//--------------------------------------------------------------------------
void AliTOFHeader::SetEventTimeRes(TArrayF *arr)
{
  //
  //
  //

  fNbins=arr->GetSize();
  fEventTimeRes->Set(arr->GetSize());
  for (Int_t ii=0; ii<arr->GetSize(); ii++)
    fEventTimeRes->SetAt(arr->GetAt(ii),ii);

}
//--------------------------------------------------------------------------
void AliTOFHeader::SetNvalues(TArrayI *arr)
{
  //
  //
  //

  fNbins=arr->GetSize();
  fNvalues->Set(arr->GetSize());
  for (Int_t ii=0; ii<arr->GetSize(); ii++)
    fNvalues->SetAt(arr->GetAt(ii),ii);

}
