#ifndef ALIMUONNODE_H
#define ALIMUONNODE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup geometry
/// \class AliMUONNode
/// \brief A node of a segment tree
/// 
// author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class TObjArray;

class AliMUONNode : public TObject
{
public:
  AliMUONNode(Double_t a, Double_t b, Double_t midpoInt_t);
  virtual ~AliMUONNode();
  
  void Print(const char* opt="") const;
  
  void Contribution(Double_t b, Double_t e, TObjArray& stack);
  
  void InsertInterval(Double_t b, Double_t e, TObjArray& stack);
  
  void DeleteInterval(Double_t b, Double_t e, TObjArray& stack);
  
  Bool_t IsFullyContained(Double_t b, Double_t e) const;
  
  void Update();
  
  void Demote();
  
  void Promote();
  
  /// Get cardinality
  Int_t C() const { return fC; }
  
  /// Increase cardinality
  void C(Int_t v) { fC += v; }
  
  /// Get potent state
  Int_t P() const { return fP; }
  
  /// Set left node
  void LeftNode(AliMUONNode* n) { fLeftNode = n; }
  
  /// Set right node
  void RightNode(AliMUONNode* n) { fRightNode = n; }
  
private:
  
  /// not implemented
  AliMUONNode(const AliMUONNode& node); 
  /// not implemented
  AliMUONNode& operator=(const AliMUONNode& node);  
  AliMUONNode* fLeftNode; ///< left node
  AliMUONNode* fRightNode; ///< right node
  
  Double_t fMin; ///< Min
  Double_t fMax; ///< Max
  Double_t fMidPoint; ///< (Min+Max)/2
  
  Int_t fC; ///< cardinality
  Int_t fP; ///< potent state
  
  ClassDef(AliMUONNode,0) // A node of a segment tree
};

#endif
