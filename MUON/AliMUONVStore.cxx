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

//-----------------------------------------------------------------------------
/// \class AliMUONVStore
///
/// A store is a container, which can be searched for (using FindObject methods),
/// iterated upon (using CreateIterator() method), and into which you can 
/// add objects (using Add)
///
/// In addition, a store can be connected to a TTree.
///
/// The general way of dealing with I/O for MUON is a two stage process :
/// 
/// 1) first get a TTree pointer using the AliLoader mechanism (this is AliRoot
///    general)
///
/// 2) connect that TTree to a MUON (virtual) data container using
///    the container's Connect(TTree&) method (this is MUON specific)
///
/// Example for reading digits for nevents
///
/// \code
///
/// AliMUONVDigitStore* digitStore(0x0);
///
/// AliLoader* loader = ... ( get loader from somewhere, e.g. AliRunLoader::GetDetectorLoader());
/// loader->LoadDigits("READ"); // load digits
///
/// for ( Int_t i = 0; i < nevents; ++i ) 
/// {
///   TTree* treeD = loader->TreeD(); // get the tree
///   if (!digitStore) digitStore = static_cast<AliMUONVDigitStore*>
///     (AliMUONVDigitStore::CreateStore(*treeD,"Digit")); // creates a container for digits 
///   (concrete class is given by the tree itself)
///   digitStore->Connect(*treeD);
///   treeD->GetEvent(0); // actual reading of the data
///
///   the digitStore is now filled and ready to be used
///   ....
///
///   digitStore->Clear(); // reset once used
/// }
///
/// \endcode
///
/// Please note that for reading data, you do *not* need to know the concrete
/// container class, as it is given to you by the TTree
///
/// Example for writing digits
///
/// \code
///
/// get the loader and do a loader->LoadDigits("RECREATE" or "UPDATE") 
/// (generally done by the framework itself)
///
/// AliMUONVDigitStore* digitStore = new AliMUONDigitStoreV1
/// // for writing, must decide on the concrete store class to use
///
/// for ( Int_t i = 0; i < nevents; ++i ) 
/// {
///   TTree* treeD = loader->TreeD();
///   digitStore->Connect(*treeD);
///   
///   ... put some digits in the digitStore
///
///   treeD->Fill();
///
///   loader->WriteDigits("OVERWRITE");
/// }
///
/// delete digitStore;
///
/// loader->UnloadDigits();
///
/// \endcode
///
/// In the write case, one *must* specify a concrete class for the container
///
//-----------------------------------------------------------------------------

#include "AliMUONVStore.h"

#include "AliMUONTreeManager.h"
#include "AliLog.h"

#include <TRegexp.h>
#include <TClass.h>

/// \cond CLASSIMP
ClassImp(AliMUONVStore)
/// \endcond

//_____________________________________________________________________________
AliMUONVStore::AliMUONVStore() : TObject()
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONVStore::~AliMUONVStore()
{
  /// dtor
}

//_____________________________________________________________________________
Bool_t
AliMUONVStore::Connect(TTree&, Bool_t) const
{
  /// Connect to a Ttree
  AliError("Not implemented");
  return kFALSE;
}

//_____________________________________________________________________________
AliMUONVStore* 
AliMUONVStore::Create(TTree& tree, const char* what)
{
  /// Create a store from a tree. Forwarded to AliMUONTreeManager::CreateStore
  AliMUONTreeManager tman;

  TObject* o = tman.CreateObject(tree,what);;
  if (o)
  {
    AliMUONVStore* c = dynamic_cast<AliMUONVStore*>(o);
    if (!c)
    {
      AliErrorClass(Form("Object of class %s cannot be cast to an AliMUONVStore",
                    o->ClassName()));
    }
    return c;
  }
  return 0x0;
}

//_____________________________________________________________________________
TObject*
AliMUONVStore::FindObject(Int_t, Int_t) const
{
  /// Find an object using 2 identifiers
  AliError("(Int_t,Int_t) : Not implemented");
  return 0;
}

//______________________________________________________________________________
TObject*
AliMUONVStore::FindObject(const char *name) const
{
  // Find an object in this collection using its name. Requires a sequential
  // scan till the object has been found. Returns 0 if object with specified
  // name is not found.
  
  TIter next(CreateIterator());
  TObject *obj;
  
  while ((obj = next()))
    if (!strcmp(name, obj->GetName())) return obj;
  return 0;
}

//______________________________________________________________________________
TObject*
AliMUONVStore::FindObject(const TObject *obj) const
{
  // Find an object in this store using the object's IsEqual()
  // member function. Requires a sequential scan till the object has
  // been found. Returns 0 if object is not found.
  // Typically this function is overridden by a more efficient version
  // in concrete collection classes.
  
  TIter next(CreateIterator());
  TObject *ob;
  
  while ((ob = next()))
    if (ob->IsEqual(obj)) return ob;
  return 0;
}

//_____________________________________________________________________________
TObject* 
AliMUONVStore::FindObject(UInt_t uniqueID) const
{
  /// Generic find method. Should be overriden by derived class if it can
  /// be made more efficient there.
  
  AliDebug(1,Form("uniqueID=%u",uniqueID));
  
  TIter next(CreateIterator());
  TObject* o;
  while ( ( o = next() ) ) 
  {
    if ( o->GetUniqueID() == uniqueID ) return o;
  }
  return 0x0;
}

//______________________________________________________________________________
Int_t
AliMUONVStore::GetSize(Int_t /*i*/) const
{
  /// Number of objects store for "i", whatever that can means
  AliError("Not implemented");
  return 0;
}

//______________________________________________________________________________
void 
AliMUONVStore::Print(Option_t *wildcard) const
{
  // Print all objects in this store.
  // Wildcarding is supported, e.g. wildcard="xxx*" prints only objects
  // with names matching xxx*.
  
  if (!wildcard) wildcard = "";
  TRegexp re(wildcard, kTRUE);
  Int_t nch = strlen(wildcard);
  TIter next(CreateIterator());
  TObject *object;
  
  while ((object = next())) {
    TString s = object->GetName();
    if (nch && s != wildcard && s.Index(re) == kNPOS) continue;
    object->Print();
  }
}

//_____________________________________________________________________________
void 
AliMUONVStore::Print(Option_t *wildcard, Option_t *option) const
{
  // Print all objects in this store, passing option to the
  // objects Print() method.
  // Wildcarding is supported, e.g. wildcard="xxx*" prints only objects
  // with names matching xxx*.
  
  if (!wildcard) wildcard = "";
  TRegexp re(wildcard, kTRUE);
  Int_t nch = strlen(wildcard);
  TIter next(CreateIterator());
  TObject *object;
  
  while ((object = next())) {
    TString s = object->GetName();
    if (nch && s != wildcard && s.Index(re) == kNPOS) continue;
    object->Print(option);
  }
}

