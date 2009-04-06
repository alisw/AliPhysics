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

/// \class AliMUONNode
/// 
/// A node of a segment tree
///
/// For the details of the meaning of cardinality and potent data
/// members, please see Diane L. Souvaine and Iliana Bjorling-Sachs,
/// Proceedings of the IEEE, Vol. 80, No. 9, September 1992, p. 1449
///
/// 
/// \Author Laurent Aphecetche, Subatech

#include "AliMUONNode.h"

#include "AliLog.h"
#include "AliMUONSegment.h"
#include "Riostream.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TString.h"

///\cond CLASSIMP
ClassImp(AliMUONNode)
///\endcond

//_____________________________________________________________________________
AliMUONNode::AliMUONNode(Double_t a, Double_t b, Double_t midpoint)
: fLeftNode(0x0), fRightNode(0x0), fMin(a), fMax(b), fMidPoint(midpoint), fC(0), fP(0)
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONNode::~AliMUONNode()
{
  /// dtor
  delete fLeftNode;
  delete fRightNode;
}

//_____________________________________________________________________________
void 
AliMUONNode::Print(const char* opt) const
{
  /// Printout
  cout << opt << Form("[%7.2f,%7.2f]",fMin,fMax);
  if ( !TMath::IsNaN(fMidPoint) ) cout << Form(" (%7.2f)",fMidPoint);
  cout << endl;
  
  TString sopt(opt);
  sopt += "   ";
  
  if ( fLeftNode ) 
  {
    fLeftNode->Print(sopt.Data());
  }
  if ( fRightNode ) 
  {
    fRightNode->Print(sopt.Data());
  }
}

//_____________________________________________________________________________
void 
AliMUONNode::Contribution(Double_t b, Double_t e, TObjArray& stack)
{
  /// Contribution of an edge (b,e) to the final contour
  if ( fMax < fMin ) 
  {
    AliError(Form("fMax(%10.5f) < fMin(%10.5f",fMax,fMin));
  }
  
  if ( fC == 0 ) 
  {
    if ( IsFullyContained(b,e) && fP == 0 ) 
    {
      AliMUONSegment* back = static_cast<AliMUONSegment*>(stack.Last());

      if ( back && AliMUONSegment::AreEqual(back->EndY(),fMin) )
      {
        // merge to existing segment
        Double_t y(back->StartY());
        back->Set(0.0,y,0.0,fMax);
      }
      else
      {
        // add a new segment
        stack.Add(new AliMUONSegment(0.0,fMin,0.0,fMax));
      }
    }
    else
    {
      if ( b < fMidPoint ) 
      {
        fLeftNode->Contribution(b,e,stack);
      }
      if ( fMidPoint < e ) 
      {
        fRightNode->Contribution(b,e,stack);
      }
    }
  }
}

//_____________________________________________________________________________
Bool_t 
AliMUONNode::IsFullyContained(Double_t b, Double_t e) const
{
  /// Whether this node's interval is fully contained into [b,e]
  
  return ( ( b < fMin || AliMUONSegment::AreEqual(b,fMin) ) && ( fMax < e || AliMUONSegment::AreEqual(e,fMax)) );
}

//_____________________________________________________________________________
void 
AliMUONNode::InsertInterval(Double_t b, Double_t e, TObjArray& stack)
{
  /// Insert an interval
  if ( IsFullyContained(b,e) ) 
  {
    C(1);
  }
  else
  {
    if ( b < fMidPoint ) 
    {
      fLeftNode->InsertInterval(b,e,stack);
    }
    if ( fMidPoint <  e ) 
    {
      fRightNode->InsertInterval(b,e,stack);
    }
  }
  Update();
}

//_____________________________________________________________________________
void 
AliMUONNode::DeleteInterval(Double_t b, Double_t e, TObjArray& stack)
{
  /// Delete an interval
  if ( IsFullyContained(b,e) ) 
  {
    C(-1);
  }
  else
  {
    if ( fC > 0 ) Demote();
    if ( b < fMidPoint )
    {
      fLeftNode->DeleteInterval(b,e,stack);
    }
    
    if ( fMidPoint < e ) 
    {
      fRightNode->DeleteInterval(b,e,stack);
    }
  }
  Update();
}

//_____________________________________________________________________________
void 
AliMUONNode::Update()
{
  /// Update internal values
  if ( !fLeftNode ) 
  {
    fP = 0;
  }
  else
  {
    if (fLeftNode->C() > 0 && fRightNode->C() > 0 )
    {
      Promote();
    }
    if (fLeftNode->C()==0 && fRightNode->C()==0 && fLeftNode->P()==0 && fRightNode->P()==0 ) 
    {
      fP = 0;
    }
    else
    {
      fP = 1;
    }
  }
}

//_____________________________________________________________________________
void 
AliMUONNode::Promote()
{
  /// Promote node
  fLeftNode->C(-1);
  fRightNode->C(-1);
  C(+1);
}

//_____________________________________________________________________________
void 
AliMUONNode::Demote()
{
  /// Demote node
  fLeftNode->C(+1);
  fRightNode->C(+1);
  C(-1);
  fP = 1;
}

