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

#include "AliMUONPainterGroup.h"

#include "AliMUONVPainter.h"
#include "AliMUONVTrackerData.h"
#include "AliLog.h"
#include <Riostream.h>
#include <TObjArray.h>
#include <float.h>

///\class AliMUONPainterGroup
///
/// A group of AliMUONVPainter
///
///\author Laurent Aphecetche, Subatech

using std::cout;
using std::endl;
///\cond CLASSIMP
ClassImp(AliMUONPainterGroup)
///\endcond

//_____________________________________________________________________________
AliMUONPainterGroup::AliMUONPainterGroup()
: TObject(),
fType(""),
fIsResponder(kFALSE),
fIsVisible(kTRUE),
fData(0x0),
fDataIndex(-1),
fDataMin(FLT_MAX),
fDataMax(-FLT_MAX),
fPainters(0x0),
fDepth(-1),
fIsOutlined(kTRUE)
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONPainterGroup::AliMUONPainterGroup(const char* type, Int_t depth)
: TObject(),
 fType(type),
 fIsResponder(kFALSE),
 fIsVisible(kTRUE),
 fData(0x0),
 fDataIndex(-1),
 fDataMin(FLT_MAX),
 fDataMax(-FLT_MAX),
 fPainters(0x0),
 fDepth(depth),
 fIsOutlined(kTRUE)
{
   /// ctor
   if ( fType == "" || fDepth < 0 ) 
   {
     AliFatal("Sorry guy.");
   }
}

//_____________________________________________________________________________
AliMUONPainterGroup::~AliMUONPainterGroup()
{
  /// dtor
  delete fPainters;
}

//_____________________________________________________________________________
Bool_t
AliMUONPainterGroup::Add(AliMUONVPainter* painter)
{
  /// Add a painter to this group (must be of the correct type)
  
   if ( fType != painter->Type() ) 
   {
     AliError(Form("Cannot add painter of type %s to this = %s",
                   painter->Type(), fType.Data()));
     return kFALSE;
   }
  
  if ( fDepth != painter->Depth() )
  {
    AliError(Form("Cannot add painter of depth %d to this = %d",
                  painter->Depth(), fDepth));
    return kFALSE;
  }
  
  if (!fPainters)
  {
    fPainters = new TObjArray;
  }
  
  painter->SetMotherGroup(this);
  
  fPainters->Add(painter);

  return kTRUE;
}

//_____________________________________________________________________________
Int_t 
AliMUONPainterGroup::Compare(const TObject* obj) const
{
  /// Compare with another group (based on type)
  
  const AliMUONPainterGroup* group = static_cast<const AliMUONPainterGroup*>(obj);
  return fType.CompareTo(group->Type());
}

//_____________________________________________________________________________
void
AliMUONPainterGroup::ComputeDataRange(Double_t& dataMin, Double_t& dataMax)
{
  /// Compute the data range spanned by this group
  dataMin = FLT_MAX;
  dataMax = -FLT_MAX;
  
  if ( !fData || fDataIndex < 0 ) return;

  TIter next(fPainters);
  AliMUONVPainter* p;
  
  while ( ( p = static_cast<AliMUONVPainter*>(next()) ) )
  {
    Double_t min, max;
    p->ComputeDataRange(*fData,fDataIndex,min,max);
    dataMin = TMath::Min(min,dataMin);
    dataMax = TMath::Max(max,dataMax);
  }
}

//_____________________________________________________________________________
void
AliMUONPainterGroup::Draw(Option_t* opt)
{
  /// Draw our painters
  TIter next(fPainters);
  TObject* o;
  while ( ( o = next() ) )
  {
    o->Draw(opt);
  }
}  

//_____________________________________________________________________________
AliMUONVPainter* 
AliMUONPainterGroup::First() const
{
  /// Get the first painter in group
  if ( fPainters ) 
  {
    return static_cast<AliMUONVPainter*>(fPainters->First());
  }
  return 0x0;
}

//_____________________________________________________________________________
Int_t
AliMUONPainterGroup::GetLineColor() const
{
  /// Get line color of this group's painters
  if ( fPainters ) 
  {
    return static_cast<AliMUONVPainter*>(fPainters->First())->GetLineColor();
  }
  return 1;
}

//_____________________________________________________________________________
Int_t
AliMUONPainterGroup::GetLineWidth() const
{
  /// Get line width of this group's painters
  if ( fPainters ) 
  {
    return static_cast<AliMUONVPainter*>(fPainters->First())->GetLineWidth();
  }
  return 1;
}

//_____________________________________________________________________________
Bool_t
AliMUONPainterGroup::Matches(const char* pattern) const
{
  /// Whether our type matches "pattern"
  TString spattern(pattern);
  
  if ( spattern == "*" || fType.Contains(pattern) )
  {
    return kTRUE;
  }
  return kFALSE;
}

//_____________________________________________________________________________
void
AliMUONPainterGroup::Print(Option_t* opt) const
{
  /// Printout
  cout << "Type " << fType.Data() << " Depth " << fDepth;
  if ( IsResponder() ) cout << " is responder ";
  if ( IsVisible() ) cout << " is visible ";
  if ( IsPlotter() ) 
  {
    cout << Form(" is plotter for data %p %s dimension %d %s plot range = %e, %e",
                 fData,(fData ? fData->Name() : ""),
                 fDataIndex,( (fData && fDataIndex>=0 ) ? 
                              fData->DimensionName(fDataIndex).Data() : ""),
                 DataMin(),DataMax());
  }
  if ( IsOutlined() ) 
  {
    cout << " is outlined";
  }
  if ( fPainters ) 
  {
    cout << " contains " << fPainters->GetLast()+1 << " painters";
  }
  
  cout << endl;
  
  TString sopt(opt);
  sopt.ToUpper();
  if ( sopt == "FULL" ) 
  {
    TIter next(fPainters);
    AliMUONVPainter* painter;
    while ( ( painter = static_cast<AliMUONVPainter*>(next()) ) )
    {
      cout << "    ";
      painter->Print();
    }
  }
}

//_____________________________________________________________________________
void 
AliMUONPainterGroup::SetData(AliMUONVTrackerData* data, Int_t dataIndex)
{ 
  /// Set the data to be plotted
  fData = data; 
  fDataIndex = dataIndex; 
  fDataMax = -FLT_MAX;
  fDataMin = FLT_MAX;
}

//_____________________________________________________________________________
void
AliMUONPainterGroup::SetLine(Int_t lineColor, Int_t lineWidth)
{
  /// Set our outline attributes
  TIter next(fPainters);
  AliMUONVPainter* painter;
  while ( ( painter = static_cast<AliMUONVPainter*>(next()) ) )
  {
    painter->SetLineColor(lineColor);
    painter->SetLineWidth(lineWidth);
  }
}  

