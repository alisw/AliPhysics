#ifndef ALIMUONSEGMENTTREE_H
#define ALIMUONSEGMENTTREE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup geometry
/// \class AliMUONSegmentTree
/// \brief Implementation of a segment tree
/// 
// author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

#ifndef ROOT_TObjArray
#  include "TObjArray.h"
#endif

class TArrayD;
class AliMUONNode;

class AliMUONSegmentTree : public TObject
{
public:
  AliMUONSegmentTree(const TArrayD& values);
  virtual ~AliMUONSegmentTree();
  
  AliMUONNode* Build(const TArrayD& values, Int_t i, Int_t j);
  
  void Print(Option_t* opt="") const;
  
  /// Get the stack
  const TObjArray& Stack() const { return fStack; }
  
  /// Reset the stack
  void ResetStack() { fStack.Clear(); }
  
  void Contribution(double b, double e);
  
  void InsertInterval(double b, double e);
  
  void DeleteInterval(double d, double e);
  
private:
  /// not implemented
  AliMUONSegmentTree(const AliMUONSegmentTree& rhs);
  /// not implemented
  AliMUONSegmentTree& operator=(const AliMUONSegmentTree& rhs);
  
  AliMUONNode* fRoot; ///< root of the tree
  TObjArray fStack; ///< array of AliMUONSegment objects
  
  ClassDef(AliMUONSegmentTree,1) // Implementation of a segment tree
};

#endif
