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

///
/// \class AliMUONSegmentTree
///
/// Implementation of a segment tree, which is used to make contour
/// merging (see AliMUONContourMaker)
///
/// \author Laurent Aphecetche, Subatech
///

#include "AliMUONSegmentTree.h"

#include "AliLog.h"
#include "AliMUONNode.h"
#include "Riostream.h"
#include "TArrayD.h"
#include "TMath.h"

/// \cond CLASSIMP
ClassImp(AliMUONSegmentTree)
/// \endcond

//_____________________________________________________________________________
AliMUONSegmentTree::AliMUONSegmentTree(const TArrayD& values)
: fRoot(0x0), fStack()
{
  /// Values should be sorted and have at least 2 elements.
  
  fStack.SetOwner(kTRUE);
  
  if ( values.GetSize() < 2 ) 
  {
    AliError("cannot build a segmenttree with less than 2 values !");
    TObject* forceACrash(0x0);
    forceACrash->Print();
  }
  
  fRoot = Build(values,0,values.GetSize()-1);
}

//_____________________________________________________________________________
AliMUONSegmentTree::~AliMUONSegmentTree()
{
  /// dtor
  delete fRoot;
}

//_____________________________________________________________________________
AliMUONNode* 
AliMUONSegmentTree::Build(const TArrayD& values, Int_t i, Int_t j)
{
  /// Build the segment tree from a list of values
  
  double midpoint(TMath::Sqrt(-1.0));
  Int_t mid((i+j)/2);
  
  if ( mid != i && mid != j ) midpoint = values[mid];
  
  AliMUONNode* node = new AliMUONNode(values[i],values[j],midpoint);
  
  if ( j - i == 1 ) return node;
  
  node->LeftNode(Build(values,i,(i+j)/2));
  node->RightNode(Build(values,(i+j)/2,j));
  
  return node;
}

//_____________________________________________________________________________
void 
AliMUONSegmentTree::Contribution(double b, double e)
{
  /// Compute the contribution of edge (b,e)
  fRoot->Contribution(b,e,fStack);
}

//_____________________________________________________________________________
void 
AliMUONSegmentTree::InsertInterval(double b, double e)
{
  /// Insert interval (b,e)
  fRoot->InsertInterval(b,e,fStack);
}

//_____________________________________________________________________________
void 
AliMUONSegmentTree::DeleteInterval(double b, double e)
{
  /// Delete interval (b,e)
  fRoot->DeleteInterval(b,e,fStack);
}

//_____________________________________________________________________________
void 
AliMUONSegmentTree::Print(Option_t*) const
{
  /// Printout
  if (fRoot) 
    fRoot->Print(); 
  else 
    cout << "Empty binary tree" << endl; 
}
