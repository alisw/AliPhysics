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

#include "AliMUONPainterMatrix.h"

#include "AliMUONPainterGroup.h"
#include "AliMUONVPainter.h"
#include "AliLog.h"
#include <Riostream.h>
#include <TObjArray.h>
#include <float.h>
#include <TObjString.h>

///\class AliMUONPainterMatrix
///
/// Matrix of AliMUONVPainter
///
///\author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONPainterMatrix)
///\endcond

//_____________________________________________________________________________
AliMUONPainterMatrix::AliMUONPainterMatrix(const char* name, Int_t nx, Int_t ny)
: TObject(),
  fBasename(name),
  fName(""),
  fNx(nx),
  fNy(ny),
  fPainters(new TObjArray(fNx*fNy)),
  fAttributes()
{
    /// ctor

    fPainters->SetOwner(kTRUE);
    if ( fNx*fNy > 1 ) 
    {
      fAttributes.SetSingle(kFALSE);
    }
    
    fName = NameIt(name,fAttributes);
}

//_____________________________________________________________________________
AliMUONPainterMatrix::~AliMUONPainterMatrix()
{
  /// dtor
  delete fPainters;
}

//_____________________________________________________________________________
void 
AliMUONPainterMatrix::Adopt(AliMUONVPainter* painter)
{
  /// Adopt a given painter
  fPainters->AddLast(painter);
  UpdateAttributes();
}

//_____________________________________________________________________________
void
AliMUONPainterMatrix::UpdateAttributes()
{
  /// Update our attributes (using our painters' attributes)
  
  Bool_t cathode0(kFALSE);
  Bool_t cathode1(kFALSE);
  Bool_t bending(kFALSE);
  Bool_t nonbending(kFALSE);
  Bool_t front(kFALSE);
  Bool_t back(kFALSE);
  Bool_t cathplaneexclusive(kFALSE);
  Bool_t cathplanedisabled(kFALSE);
  
  for ( Int_t i = 0; i < Size(); ++i )
  {
    AliMUONAttPainter att = Painter(i)->Attributes();
    
    if ( att.IsCathodeDefined() ) 
    {
      if ( att.IsCathode0() ) cathode0 = kTRUE;
      if ( att.IsCathode1() ) cathode1 = kTRUE;
    }

    if ( att.IsPlaneDefined() ) 
    {
      if ( att.IsBendingPlane() ) bending = kTRUE;
      if ( att.IsNonBendingPlane() ) nonbending = kTRUE;
    }
    
    if ( att.IsFrontView() ) front = kTRUE;
    if ( att.IsBackView() ) back = kTRUE;
    
    if ( att.IsCathodeAndPlaneMutuallyExclusive() ) cathplaneexclusive = kTRUE;
    
    if ( att.IsCathodeAndPlaneDisabled() ) cathplanedisabled = kTRUE;
  }
  
  fAttributes.SetCathode(cathode0,cathode1);
  fAttributes.SetPlane(bending,nonbending);
  fAttributes.SetViewPoint(front,back);
  fAttributes.SetCathodeAndPlaneMutuallyExclusive(cathplaneexclusive);
  fAttributes.SetCathodeAndPlaneDisabled(cathplanedisabled);
  
  fName = NameIt(fBasename,fAttributes);
}

//_____________________________________________________________________________
TString
AliMUONPainterMatrix::NameIt(const TString& basename, const AliMUONAttPainter& att)
{
  /// Build a name 
  TString name(basename);
  
  name += "-";
  name += att.Name();
  
  return name;
}

//_____________________________________________________________________________
void 
AliMUONPainterMatrix::ComputeDataRange()
{
  /// Compute the data range spanned by the painters in this matrix
  
  Double_t dataMin(FLT_MAX);
  Double_t dataMax(-FLT_MAX);
  Bool_t atLeastOnePlotter(kFALSE);
  
  for ( Int_t i = 0; i < Size(); ++i ) 
  {
    AliMUONVPainter* p = Painter(i);
    AliMUONPainterGroup* g = p->PlotterGroup();

    Double_t min(FLT_MAX);
    Double_t max(-FLT_MAX);

    if ( g ) 
    {
      atLeastOnePlotter = kTRUE;
      g->ComputeDataRange(min,max);
      if ( min <= max ) 
      {
        dataMin = TMath::Min(min,dataMin);
        dataMax = TMath::Max(max,dataMax);
      }
    }

    AliDebug(1,Form("painter %s group %s min %e max %e dataMin,Max=%7.3f,%7.3f",
                    p->GetName(),
                    g ? g->Type() : "none",
                    min,max,
                    dataMin,dataMax));
  }

  if ( dataMin > dataMax && atLeastOnePlotter ) 
  {
    AliError(Form("data min %e > max %e : setting both to 0.0",
                    dataMin,dataMax));
    dataMin = dataMax = 0.0;
  }
  
  AliDebug(1,Form("Final dataMin,Max=%7.3f,%7.3f",dataMin,dataMax));
  
  SetDataRange(dataMin,dataMax);
}

//_____________________________________________________________________________
void 
AliMUONPainterMatrix::Connect(const char* sourceMethod, const char* destClassName, 
                              void* destObject, const char* destMethod)
{
  /// Connect our painters
  
  for ( Int_t i = 0; i < Size(); ++i )
  {
    Painter(i)->Connect(sourceMethod,destClassName,destObject,destMethod);
  }
}

//_____________________________________________________________________________
void
AliMUONPainterMatrix::GetDataRange(Double_t& dataMin, Double_t& dataMax) const
{
  /// Get the data range spanned by the painters in this matrix
  
  dataMin=FLT_MAX;
  dataMax=-FLT_MAX;
  
  for ( Int_t i = 0; i < Size(); ++i ) 
  {
    AliMUONVPainter* p = Painter(i);
    if ( p )
    {
      AliMUONPainterGroup* g = p->PlotterGroup();
      if ( g ) 
      {
        dataMin = TMath::Min(dataMin,g->DataMin());
        dataMax = TMath::Max(dataMax,g->DataMax());
      }
    }
  }
}

//_____________________________________________________________________________
void 
AliMUONPainterMatrix::GetTypes(TObjArray& types) const
{
  /// Get the types of the painters in this matrix
  
  types.SetOwner(kTRUE);
  types.Clear();
  
  for ( Int_t i = 0; i < Size(); ++i ) 
  {
    AliMUONVPainter* p = Painter(i);
    TObjArray ptypes;
    p->GetTypes(ptypes);
    TIter next(&ptypes);
    TObject* o;
    while ( ( o = next() ) )
    {
      if ( ! types.FindObject(o) )
      {
        types.AddLast(o->Clone());
      }
    }
  }  
}

//_____________________________________________________________________________
AliMUONVPainter* 
AliMUONPainterMatrix::Painter(Int_t index) const
{
  /// Get a given painter
  
  if ( index <= fPainters->GetLast() ) 
  {
    return static_cast<AliMUONVPainter*>(fPainters->At(index));
  }
  return 0x0;
}

//_____________________________________________________________________________
AliMUONVTrackerData* 
AliMUONPainterMatrix::Data() const
{
  /// Return our data
  AliMUONPainterGroup* group = Painter(0)->PlotterGroup();
  return ( group ? group->Data() : 0x0 );
}

//_____________________________________________________________________________
TString 
AliMUONPainterMatrix::DataPattern() const
{
  /// Return our data pattern
  AliMUONPainterGroup* group = Painter(0)->PlotterGroup();
  return ( group ? group->Type() : "" );
}

//_____________________________________________________________________________
Int_t 
AliMUONPainterMatrix::DataIndex() const
{
  /// Return our data index
  AliMUONPainterGroup* group = Painter(0)->PlotterGroup();
  return ( group ? group->DataIndex() : -1 );
}

//_____________________________________________________________________________
void
AliMUONPainterMatrix::SetData(const char* pattern, AliMUONVTrackerData* d,
                              Int_t indexInData)
{
  /// Set the data to be plotted
  
  for ( Int_t i = 0; i < Size(); ++i )
  {
    AliMUONVPainter* painter = Painter(i);
    painter->SetData(pattern,d,indexInData);
  }
}

//_____________________________________________________________________________
void
AliMUONPainterMatrix::SetDataRange(Double_t dataMin, Double_t dataMax)
{
  /// Set the data range
  
  for ( Int_t i = 0; i < Size(); ++i ) 
  {
    AliMUONVPainter* p = Painter(i);
    AliMUONPainterGroup* g = p->PlotterGroup();
    if ( g ) 
    {
      g->SetDataRange(dataMin,dataMax);
    }
  }       
}

//_____________________________________________________________________________
Int_t 
AliMUONPainterMatrix::Size() const
{
  /// Return the number of painters we actually handle
  return fPainters->GetLast()+1;
}

//_____________________________________________________________________________
void
AliMUONPainterMatrix::Print(Option_t*) const
{
  /// Printout
  cout << "Basename=" << fBasename.Data() << " Name=" << fName.Data() 
  << " Nx=" << fNx << " Ny=" << fNy << " Att=" << fAttributes.GetName() << endl;
}

//_____________________________________________________________________________
//void 
//AliMUONPainterMatrix::ChangeAttributes(const AliMUONAttPainter& attributes)
//{
//  /// Change painters' attributes
//  
//  AliWarning("Implement me !");
//  
//  //  for ( Int_t i = 0; i < Size(); ++i ) 
//  //  {
//  //    Painter(i)->SetAttributes(attributes);
//  //  }
//}

//_____________________________________________________________________________
AliMUONPainterMatrix*
AliMUONPainterMatrix::Clone(const AliMUONAttPainter& attributes) const
{
  /// Clone with given attributes
  
  AliMUONPainterMatrix* clone = new AliMUONPainterMatrix(Basename().Data(),Nx(),Ny());

  for ( Int_t i = 0; i < Size(); ++i ) 
  {
    AliMUONVPainter* oldPainter = Painter(i);
    
    AliMUONVPainter* newPainter(0x0);
    
    newPainter = AliMUONVPainter::CreatePainter(oldPainter->ClassName(),
                                                attributes,
                                                oldPainter->ID0(),
                                                oldPainter->ID1());
    
    if (newPainter)
    {
      newPainter->UpdateGroupsFrom(*(oldPainter->Master()));
      clone->Adopt(newPainter);
    }
    else
    {
      AliError(Form("Failed to create painter of class %s ID0 %d ID1 %d",
                    oldPainter->ClassName(),
                    oldPainter->ID0(),
                    oldPainter->ID1()));
    }
  }
  
  return clone;
}

//_____________________________________________________________________________
void
AliMUONPainterMatrix::SetOutlined(const char* pattern, Bool_t value)
{
  /// Calls SetOutlined for all our painters

  for ( Int_t i = 0; i < Size(); ++i ) 
  {
    Painter(i)->SetOutlined(pattern,value);
  }
}

//_____________________________________________________________________________
void
AliMUONPainterMatrix::SetResponder(const char* pattern)
{
  /// Calls SetResponder for all our painters
  for ( Int_t i = 0; i < Size(); ++i ) 
  {
    Painter(i)->SetResponder(pattern);
  }
}

//_____________________________________________________________________________
AliMUONAttPainter
AliMUONPainterMatrix::Validate(const AliMUONAttPainter& att) const
{
  /// Normalize attributes

  AliMUONAttPainter a;
  
  for ( Int_t i = 0; i < Size() && a.IsValid(); ++i ) 
  {
    a = Painter(i)->Validate(att);
  }
  return a;
}


