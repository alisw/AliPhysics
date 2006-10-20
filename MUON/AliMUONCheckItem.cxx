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

// $Id$

#include "AliMUONCheckItem.h"

#include "AliLog.h"
#include "AliMpExMap.h"
#include "Riostream.h"
#include "AliMUONCheckItemIterator.h"

/// \class AliMUONCheckItem
///
/// A structure used to gather information at different levels (ch,manu,de,chamber)
///
/// Used by AliMUON2DStoreValidator to present results in a concise way
///
/// 
/// \author Laurent Aphecetche

/// \cond CLASSIMP
ClassImp(AliMUONCheckItem)
/// \endcond

//_____________________________________________________________________________
AliMUONCheckItem::AliMUONCheckItem(Int_t id, Int_t maxNumber, const char* name) : 
TNamed(name,name),
fID(id),
fDead(-1),
fMaximum(maxNumber),
fMissing(new AliMpExMap(kTRUE))
{
  /// ctor. id is the number of that item, maxNumber is the maximum number
  /// of sub-item it can contains, and name is a label, e.g. de, chamber, manu.
  /// Note that name="manu" has a special influence on the IsDead() method.
  
  fMissing->SetSize(fMaximum);
  AliDebug(1,Form("ID %d maxNumber %d name %s",id,maxNumber,name));
}

//_____________________________________________________________________________
AliMUONCheckItem::~AliMUONCheckItem()
{
  /// dtor
  delete fMissing;
}

//_____________________________________________________________________________
Bool_t AliMUONCheckItem::AddItem(Int_t id, TObject* item)
{
  /// Add an item, if possible.
  
  if ( IsFull() ) 
  {
    AliError("I'm already full!");
    return kFALSE;
  }
  
  TObject* test = GetItem(id);
  if (test) 
  {
    AliError(Form("id %d is already there !",id));
    return kFALSE;
  }
  else
  {
    fMissing->Add(id,item);
    fDead=-1;
  }  
  return kTRUE;
}

//_____________________________________________________________________________
void
AliMUONCheckItem::ComputeDead() const
{
  /// Decide whether this item is completely dead, which is determined by
  /// the fact that all its sub-items are dead, or for name="manu", by
  /// the fact that all channels are missing, i.e. IsFull()==kTRUE
  
  TString name(GetName());
  name.ToLower();

  if ( name.Contains("manu") )
  {
    if ( IsFull() ) 
    {
      fDead=1;
    }
    else
    {
      fDead=0;
    }
  }
  else
  {
    AliMUONCheckItemIterator it(*this);
    AliMUONCheckItem* item;
    it.First();
    Int_t ndead(0);
    fDead=0;
    while ( ( item = dynamic_cast<AliMUONCheckItem*>(it.Next()) ) )
    {
      if ( item->IsDead() ) ++ndead;
    }
    if ( ndead == fMaximum ) fDead = 1;
  }
}

//_____________________________________________________________________________
TObject* 
AliMUONCheckItem::GetItem(Int_t id) const
{
  /// Return item of a given id
  return fMissing->GetValue(id);
}

//_____________________________________________________________________________
Bool_t 
AliMUONCheckItem::IsDead() const
{
  /// Return (and compute it first if not done already) dead status
  if ( fDead == -1 )
  {
    ComputeDead();
  }
  return (fDead==1);
}

//_____________________________________________________________________________
Bool_t 
AliMUONCheckItem::IsFull() const
{
  /// Whether we have as many sub-items as possible
  return (fMissing->GetSize() == fMaximum);
}
  
//_____________________________________________________________________________
void
AliMUONCheckItem::Print(Option_t* opt) const
{
  /// output to screen
  cout << Form("<AliMUONCheckItem> %s ID %d has %d items over %d max. Dead %d",
               GetName(),fID,fMissing->GetSize(),fMaximum,IsDead()) << endl;  
  TString sopt(opt);
  sopt.ToLower();
  if (sopt.Contains("all") )
  {
    TObject* object(0x0);
  
    AliMUONCheckItemIterator it(*this);
  
    it.First();
  
    while ( ( object = it.Next() ) )
    {
      object->Print(opt);
    }
  }
}
