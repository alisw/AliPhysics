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

#include "AliMUONPainterRegistry.h"

#include "AliMpManuIterator.h"
#include "AliMUON2DMap.h"
#include "AliMUONCalibParamND.h"
#include "AliMUONPainterMatrix.h"
#include "AliLog.h"
#include <TGMenu.h>
#include <TGWindow.h>
#include <THashList.h>
#include <TObjArray.h>
#include <TString.h>
#include <Riostream.h>

///\class AliMUONPainterRegistry
///
/// Registry for AliMUONVPainter related stuff : painter data sources
/// and painter matrices
///
///\author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONPainterRegistry)
///\endcond

AliMUONPainterRegistry* AliMUONPainterRegistry::fgInstance(0x0);

//_____________________________________________________________________________
AliMUONPainterRegistry::AliMUONPainterRegistry() : TObject(), TQObject(),
fPainterMatrices(new TObjArray),
fHistoryMenu(0x0),
fMenuBar(0x0),
fHistoryCounter(0)
{
  /// ctor
  fPainterMatrices->SetOwner(kTRUE);
}

//_____________________________________________________________________________
AliMUONPainterRegistry::~AliMUONPainterRegistry()
{
  /// dtor
  delete fPainterMatrices;
}

//_____________________________________________________________________________
Int_t 
AliMUONPainterRegistry::FindIndexOf(AliMUONPainterMatrix* group) const
{
  /// Get the index of a given painterMatrix
  return fPainterMatrices->IndexOf(group);
}

//_____________________________________________________________________________
AliMUONPainterMatrix*
AliMUONPainterRegistry::PainterMatrix(const char* name) const
{
  /// Get a painterMatrix by name
  return static_cast<AliMUONPainterMatrix*>(fPainterMatrices->FindObject(name));
}

//_____________________________________________________________________________
void
AliMUONPainterRegistry::HistoryMenuActivated(Int_t i)
{
  /// A painterMatrix was chosen from the history menu
  
  AliDebug(1,Form("i=%d",i));
  
  TGMenuEntry* entry = fHistoryMenu->GetEntry(i);
  
  AliMUONPainterMatrix* group = reinterpret_cast<AliMUONPainterMatrix*>(entry->GetUserData());
  
  PainterMatrixWantToShow(group);
}

//_____________________________________________________________________________
AliMUONPainterRegistry*
AliMUONPainterRegistry::Instance()
{
  /// Get unique instance of this class
  if ( !fgInstance ) fgInstance = new AliMUONPainterRegistry;
  return fgInstance;
}

//_____________________________________________________________________________
AliMUONPainterMatrix* 
AliMUONPainterRegistry::PainterMatrix(Int_t i) const
{
  /// Get one painter matrix
  if ( i >= 0 && i <= fPainterMatrices->GetLast() )
  {
    return static_cast<AliMUONPainterMatrix*>(fPainterMatrices->At(i));
  }
  else
  {
    AliError(Form("Index out of bounds : %d / %d",i,fPainterMatrices->GetLast()+1));
    return 0x0;
  }
}

//_____________________________________________________________________________
void 
AliMUONPainterRegistry::PainterMatrixWantToShow(const AliMUONPainterMatrix* group)
{
  /// A given paintermatrix want to appear on screen
  Long_t param[] = { (Long_t)group };
  
  Emit("PainterMatrixWantToShow(AliMUONPainterMatrix*)",param);  
}

//_____________________________________________________________________________
void
AliMUONPainterRegistry::AddToHistory(AliMUONPainterMatrix* group)
{
  /// Add a matrix to the history
  
  if ( !fHistoryMenu && fMenuBar ) 
  {
    fHistoryMenu =  new TGPopupMenu(gClient->GetRoot());
    TGPopupMenu* before = 0x0; //FIXME: could try to find a place where to put it (e.g. before Help ?)
    
    fMenuBar->AddPopup("&History",fHistoryMenu, new TGLayoutHints(kLHintsNormal),before);
    
    fHistoryMenu->Connect("Activated(Int_t)",
                          "AliMUONPainterRegistry",this,
                          "HistoryMenuActivated(Int_t)");
    
    AliDebug(1,Form("HistoryMenu create at %x",fHistoryMenu));
  }
  
  if ( fHistoryMenu ) 
  {
    TIter next(fHistoryMenu->GetListOfEntries());
    TGMenuEntry* e(0x0);
    
    while ( ( e = static_cast<TGMenuEntry*>(next()) ) )
    {
      if ( e->GetUserData() == group ) 
      {
        fHistoryMenu->DeleteEntry(e);
        break;
      }
    }
    
    e = static_cast<TGMenuEntry*>(fHistoryMenu->GetListOfEntries()->First());
    
    fHistoryMenu->AddEntry(group->GetName(),++fHistoryCounter,(void*)group,0x0,e);
  }
  else
  {
    AliError("fHistoryMenu is null. We probably did not find the relevant menu entry ?");
  }
}

//_____________________________________________________________________________
void 
AliMUONPainterRegistry::PainterMatrixWasRegistered(const AliMUONPainterMatrix* group)
{
  /// A new painter matrix was registered
  Long_t param[] = { (Long_t)group };
  
  Emit("PainterMatrixWasRegistered(AliMUONPainterMatrix*)",param);
}

//_____________________________________________________________________________
void 
AliMUONPainterRegistry::PainterMatrixWasUnregistered(const AliMUONPainterMatrix* group)
{
  /// A painter matrix was unregistered
  Long_t param[] = { (Long_t)group };
  
  Emit("PainterMatrixWasUnregistered(AliMUONPainterMatrix*)",param);
}

//_____________________________________________________________________________
void 
AliMUONPainterRegistry::Print(Option_t* opt) const
{
  /// Printout
  TString sopt(opt);
  sopt.ToUpper();
  
  cout << "Number of painter matrices = " << NumberOfPainterMatrices() << endl;
  
  if ( sopt.Contains("FULL") || sopt.Contains("MATRIX") )
  {
    TIter next(fPainterMatrices);
    AliMUONPainterMatrix* matrix;
    
    while ( ( matrix = static_cast<AliMUONPainterMatrix*>(next()) ) )
    {
      matrix->Print();
    }
  }
  
}

//_____________________________________________________________________________
Int_t
AliMUONPainterRegistry::Register(AliMUONPainterMatrix* group)
{
  /// group is adopted, i.e. the registry becomes the owner of it.
  fPainterMatrices->AddLast(group);
  
  PainterMatrixWasRegistered(group);
  
  return fPainterMatrices->IndexOf(group);
}

//_____________________________________________________________________________
Int_t 
AliMUONPainterRegistry::NumberOfPainterMatrices() const
{
  /// The number of painter matrices we handle
  return fPainterMatrices->GetLast()+1;
}

//_____________________________________________________________________________
Bool_t 
AliMUONPainterRegistry::Unregister(AliMUONPainterMatrix* group)
{
  /// Unregister some matrix
  
  if (!group) return kFALSE;
  
  PainterMatrixWasUnregistered(group);
  
  TObject* o = fPainterMatrices->Remove(group);
  if ( o ) 
  {
    delete o;
  }
  else
  {
    AliError(Form("Could not unregister group named %s",group->GetName()));
  }
  return ( o != 0x0 );
}
