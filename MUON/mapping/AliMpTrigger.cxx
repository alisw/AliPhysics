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
* about the suitability of this software for any purpeateose. It is      *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id$
// $MpId: AliMpTrigger.cxx,v 1.2 2006/03/02 16:35:31 ivana Exp $

#include "AliMpTrigger.h"

#include "AliLog.h"
#include "AliMpSlat.h"

#include "Riostream.h"
#include "TArrayI.h"
#include "TObjArray.h"

ClassImp(AliMpTrigger)

namespace
{
  Bool_t IsEqual(Double_t a, Double_t b, Double_t precision)
{
    if (b)
    {
      Double_t diff = TMath::Abs(b-a)/TMath::Abs(b);
      if ( diff < precision ) 
      {
        return kTRUE;
      }
    }
    else
    {
      if ( !a ) return kTRUE;
    }
    return kFALSE;
}
}

//_____________________________________________________________________________
AliMpTrigger::AliMpTrigger()
: TObject(), fId(""), fPlaneType(kNonBendingPlane), 
fMaxNofPadsY(0), fDX(0), fDY(0)
{
}

//_____________________________________________________________________________
AliMpTrigger::AliMpTrigger(const char* slatType, AliMpPlaneType bendingOrNot)
:  TObject(), fId(slatType), fPlaneType(bendingOrNot), 
fMaxNofPadsY(0), fDX(0), fDY(0)
{
}

//_____________________________________________________________________________
AliMpTrigger::~AliMpTrigger()
{
  AliDebug(1,Form("this=%p before fSlats.Delete()",this));			
  fSlats.Delete();
  AliDebug(1,Form("this=%p after fSlats.Delete()",this));			
}

//_____________________________________________________________________________
Bool_t
AliMpTrigger::AdoptLayer(AliMpSlat* slat)
{
  AliDebug(1,Form("%s is adopting %s :\n",
                  GetID(),slat->GetID()));

  // Check that we keep our size constant.
  
  const Double_t precision = 1E-3;
  
  if ( GetSize() > 0 && 
       ( !::IsEqual(slat->DX(),fDX,precision) || 
         !::IsEqual(slat->DY(),fDY,precision) )
     )
  {
    AliError(Form("In %s trying to add a layer (%e,%e) of a different size than "
             "mine (%e,%e)\n",GetID(),slat->DX(),slat->DY(),
                  fDX,fDY));
    return kFALSE;
  }
  fSlats.Add(slat);
  fMaxNofPadsY = std::max(slat->GetMaxNofPadsY(),fMaxNofPadsY);
  fDX = std::max(fDX,slat->DX());
  fDY = std::max(fDY,slat->DY());
  return kTRUE;
}

//_____________________________________________________________________________
TVector2
AliMpTrigger::Dimensions() const
{
  return TVector2(DX(),DY());
}

//_____________________________________________________________________________
Double_t
AliMpTrigger::DX() const
{
  return fDX;
}

//_____________________________________________________________________________
Double_t
AliMpTrigger::DY() const
{
  return fDY;
}

//_____________________________________________________________________________
void 
AliMpTrigger::GetAllLocalBoardNumbers(TArrayI& lbn) const
{
  Int_t n(0);
  for ( Int_t i = 0; i < GetSize(); ++i )
  {
    n += GetLayer(i)->GetNofElectronicCards();
  }
  
  lbn.Set(n);

  Int_t index(0);
  
  for ( Int_t i = 0; i < GetSize(); ++i )
  {
    TArrayI slbn;
    GetLayer(i)->GetAllMotifPositionsIDs(slbn);
    for ( Int_t j = 0; j < slbn.GetSize(); ++j )
    {
      lbn[index] = slbn[j];
      ++index;
    }
  }
}

//_____________________________________________________________________________
const char*
AliMpTrigger::GetID() const
{
  return fId.Data();
}

//_____________________________________________________________________________
const char*
AliMpTrigger::GetName() const
{
  TString name(GetID());
  if ( fPlaneType == kBendingPlane )
  {
    name += ".Bending";
  }
  else if ( fPlaneType == kNonBendingPlane )
  {
    name += ".NonBending";
  }
  else
  {
    name += ".Invalid";
  }
  return name.Data();
}

//_____________________________________________________________________________
AliMpSlat*
AliMpTrigger::GetLayer(int layer) const
{
  if ( IsLayerValid(layer) )
  {
    return (AliMpSlat*)fSlats.At(layer);
  }
  return 0;
}

//_____________________________________________________________________________
Int_t
AliMpTrigger::GetNofPadsX() const
{
  if ( !GetSize() ) return -1;
  if ( GetLayer(0) )
  {
    return GetLayer(0)->GetNofPadsX();
  }
  return -1;
}

//_____________________________________________________________________________
Int_t
AliMpTrigger::GetMaxNofPadsY() const
{
  return fMaxNofPadsY;
}

//_____________________________________________________________________________
Int_t
AliMpTrigger::GetSize() const
{
  return fSlats.GetEntriesFast();
}

//_____________________________________________________________________________
Bool_t
AliMpTrigger::IsLayerValid(int layer) const
{
  if ( layer >= 0 && layer < GetSize() )
  {
    return kTRUE;
  }
  return kFALSE;
}

//_____________________________________________________________________________
AliMpPlaneType
AliMpTrigger::PlaneType() const
{
  return fPlaneType;
}

//_____________________________________________________________________________
TVector2
AliMpTrigger::Position() const
{
  return TVector2(DX(),DY());
}

//_____________________________________________________________________________
void
AliMpTrigger::Print(Option_t* opt) const
{
  cout << "AliMpTrigger::" << GetID();
  if ( GetSize() == 0 )
  {
    cout << " Empty";
  }
  else if ( GetSize() > 1 )
  {
    cout << " Number of layers : " << GetSize();
  }
  else 
  {
    cout << " One layer";
  }
  cout << endl;
  for ( Int_t i = 0; i < GetSize(); ++i ) 
  {
    cout << "   ";
    GetLayer(i)->Print(opt);
  }
}

//_____________________________________________________________________________
//_____________________________________________________________________________
//_____________________________________________________________________________
//_____________________________________________________________________________
