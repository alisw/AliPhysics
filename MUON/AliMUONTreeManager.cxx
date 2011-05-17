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
/// \class AliMUONTreeManager
///
/// Helper class to handle the relationships TTree<->MUON data containers
/// 
/// The general way of dealing with I/O for MUON is a two stage process :
/// 
/// 1) first get a TTree pointer using the AliLoader mechanism (this is AliRoot
///    general)
///
/// 2) connect that TTree to a MUON (virtual) data container using
///    the container's Connect(TTree&) method (this is MUON specific)
///
/// This class helps implementing stage 2 in the relevant store implementations
///
/// \see AliMUONVStore
///
/// Basically, the relationship Tree<->Store is possible because
/// the TreeManager, when creating branches, uses the UserInfo part of the
/// TTree to store the relationship branch -> MUON store classname.
///
/// \author Laurent Aphecetche, Subatech
//-----------------------------------------------------------------------------

#include "AliMUONTreeManager.h"

#include "AliLog.h"
#include "AliMUONObjectPair.h"
#include <TList.h>
#include <TObjString.h>
#include <TTree.h>
#include <TBranch.h>
#include <Riostream.h>

/// \cond CLASSIMP 
ClassImp(AliMUONTreeManager)
/// \endcond

//_____________________________________________________________________________
AliMUONTreeManager::AliMUONTreeManager() : TObject()
{
  /// Default ctor
}

//_____________________________________________________________________________
AliMUONTreeManager::~AliMUONTreeManager()
{
  /// Dtor
}

//_____________________________________________________________________________
void
AliMUONTreeManager::GetEvent(TTree& tree, Int_t event) const
{
  /// Equivalent to tree.GetEvent(event) is NDEBUG is not defined,
  /// otherwise insure that selected branches have a non-zero address
  /// (the contrary indicating we'll get a memory leak when reading the
  /// tree).
  
#ifndef NDEBUG
    TObjArray* branches = tree.GetListOfBranches();
    TIter next(branches);
    TBranch* branch;
    Bool_t error(kFALSE);
    
    while ( ( branch = static_cast<TBranch*>(next()) ) )
    {
      TString bname(branch->GetName());
      if ( branch->GetAddress() == 0 &&
           tree.GetBranchStatus(bname.Data()) == 1 ) 
      {
        AliError(Form("Branch %s has status 1 and no address",bname.Data()));
        error = kTRUE;
      }
    }
    
    if ( error ) 
    {
      ShowStatus(tree);
    }
#endif
  
  tree.GetEvent(event);
}

//_____________________________________________________________________________
void
AliMUONTreeManager::ShowStatus(TTree& tree) const
{
  /// Show main branches status and address
  TObjArray* branches = tree.GetListOfBranches();
  TIter next(branches);
  TBranch* branch;
  
  while ( ( branch = static_cast<TBranch*>(next()) ) )
  {
    TString bname(branch->GetName());
    Int_t status = tree.GetBranchStatus(bname.Data());
    cout << Form("%50s Status %d Address %p",bname.Data(),status,branch->GetAddress()) << endl;
  }
}

//_____________________________________________________________________________
void 
AliMUONTreeManager::AddClassName(TTree& tree, const char* pattern, 
                                 const char* className) const
{
  /// Adds a association (pattern,className) to the UserInfo() of tree
  /// It is mandatory to use this method in MakeBranches(), as this is the key
  /// to get an automatic container creation from a tree.
  
  TString test(GetClassName(tree,pattern,kFALSE));
  if ( test.Length() == 0 )
  {
    // not already there
    TList* userInfo = tree.GetUserInfo();
    userInfo->Add(new AliMUONObjectPair(new TObjString(pattern),
                                        new TObjString(className),
                                        kTRUE,kTRUE));
  }
}

//_____________________________________________________________________________
Bool_t
AliMUONTreeManager::MakeBranch(TTree& tree, const char* storeClassName,
                               const char* branchClassName,
                               const char* branchName, 
                               void* address,
                               Int_t bufferSize, Int_t splitLevel) const
{
  /// Create a branch in the tree
  AddClassName(tree,branchName,storeClassName);
  TBranch* branch = tree.Branch(branchName,branchClassName,address,bufferSize,splitLevel);
  return ( branch != 0x0 );
}

//_____________________________________________________________________________
Bool_t 
AliMUONTreeManager::SetAddress(TTree& tree, const char* branchName, 
                               void* address) const
{
  /// Set the address for one branch
  TBranch* branch = tree.GetBranch(branchName);
  if (branch)
  {
    branch->SetAddress(address);
    return kTRUE;
  }
  AliError(Form("Did not find branch %s in tree %s",branchName,tree.GetName()));
  return kFALSE;
}

//_____________________________________________________________________________
TObject* 
AliMUONTreeManager::CreateObject(const TTree& tree, const char* detail) const
{
  /// Object creation from tree, using
  /// the (pattern,className) pairs stored in the UserInfo() of TTree.
  TString className(GetClassName(tree,detail,kTRUE));
  TClass* classPtr = TClass::GetClass(className.Data());
  if (!classPtr)
  {
    AliError(Form("Could not get class %s",className.Data()));
    return 0x0;
  }
  return reinterpret_cast<TObject*>(classPtr->New());
}

//_____________________________________________________________________________
void 
AliMUONTreeManager::UpdateBranchStatuses(TTree& tree, const char* pattern) const
{
  /// Loop over tree branches and set their status if their name matches
  /// the pattern : 
  /// - zero if branch address is null and 
  /// - one otherwise
  /// This will avoid memory leak, if we set the address to only part
  /// of the branches before doing a tree.GetEvent(i)
  ///
  /// WARNING : this is a time consuming operation, as we loop over branches
  /// at least twice (one here, and the TTree::SetBranchStatus is itself
  /// a loop). So use only when necessary.
  ///
  
  TIter next(tree.GetListOfBranches());
  TBranch* branch;
  
  while ( ( branch = static_cast<TBranch*>(next()) ) )
  {
    TString bname(branch->GetName());
    if ( bname.Contains(pattern) )
    {
      tree.SetBranchStatus(Form("%s*",branch->GetName()),1);
    }
    else
    {
      if ( !branch->GetAddress() )
      {
        tree.SetBranchStatus(Form("%s*",branch->GetName()),0);
      }
    }
  }

}

//_____________________________________________________________________________
const char* 
AliMUONTreeManager::GetClassName(const TTree& tree, const char* pattern,
                                 Bool_t makeDefault) const
{
  /// Find out, using the TTree::UserInfo, the classname corresponding to 
  /// pattern.
  /// If makeDefault=true and we cannot find the pattern in the UserInfo,
  /// we return DefaultClassName(pattern)
  ///
  
  TTree& vtree = const_cast<TTree&>(tree); // not pretty, but the GetUserInfo is not const...
  
  TIter next(vtree.GetUserInfo());
  
  TObject* object;
  
  while ( ( object = next() ) )
  {
    AliMUONObjectPair* pair = static_cast<AliMUONObjectPair*>(object);
    TString key = (static_cast<TObjString*>(pair->First()))->String();
    TString value = (static_cast<TObjString*>(pair->Second()))->String();
    if ( key.Contains(pattern,TString::kIgnoreCase) ) 
    {
      return value.Data();
    }
  }
  
  if ( makeDefault ) return DefaultClassName(tree.GetName(),pattern);
  
  return "";
}

//_____________________________________________________________________________
const char* 
AliMUONTreeManager::DefaultClassName(const char* treeName, const char* pattern) const
{
  /// For backward compatibility only. Decides, based on the tree name and a 
  /// pattern, which store class should be used.
  
  TString name(treeName);
  TString spattern(pattern);
  spattern.ToUpper();
  
  if ( name == "TreeH" ) 
  {
    return "AliMUONHitStoreV1";
  }
  
  if ( name == "TreeD" || name == "TreeS" ) 
  {
    if ( spattern.Contains("TRIGGER") ) return "AliMUONTriggerStoreV1";
    if ( spattern.Contains("DIGIT") ) return "AliMUONDigitStoreV1";
  }
  
  if ( name == "TreeR" ) 
  {
    if ( spattern.Contains("CLUSTER") ) return "AliMUONClusterStoreV1";
    if ( spattern.Contains("TRIGGER") ) return "AliMUONTriggerStoreV1";
  }
  
  if ( name == "TreeT" ) 
  {
    if ( spattern.Contains("TRIGGER" ) ) return "AliMUONTriggerTrackStoreV1";
    if ( spattern.Contains("TRACK") ) return "AliMUONTrackStoreV1";
  }
  
  AliError(Form("Do not know how to create default class for tree %s pattern %s",
                treeName,pattern));
  return "";
}







