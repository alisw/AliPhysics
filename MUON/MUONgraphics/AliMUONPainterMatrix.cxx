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

#include "AliLog.h"
#include "AliMUONPainterGroup.h"
#include "AliMUONPainterHelper.h"
#include "AliMUONVPainter.h"
#include "AliMUONVTrackerData.h"
#include "TCanvas.h"
#include "TGClient.h"
#include "TPaveLabel.h"
#include <Riostream.h>
#include <TBox.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TVirtualPad.h>
#include <float.h>
#include "AliMUONPainterEnv.h"
#include <cassert>

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
  fWhatname(""),
  fNx(nx),
  fNy(ny),
  fPainters(new TObjArray(fNx*fNy)),
  fAttributes(),
  fName()
{
  /// ctor
  
  fPainters->SetOwner(kTRUE);
  if ( fNx*fNy > 1 ) 
  {
    fAttributes.SetSingle(kFALSE);
  }
  SetName();
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
AliMUONPainterMatrix*
AliMUONPainterMatrix::Clone(const AliMUONAttPainter& attributes) const
{
  /// Clone with given attributes
  
  AliMUONPainterMatrix* clone = new AliMUONPainterMatrix(Basename(),Nx(),Ny());

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
  }

  if ( dataMin > dataMax && atLeastOnePlotter ) 
  {
    AliError(Form("data min %e > max %e : setting both to 0.0",
                    dataMin,dataMax));
    dataMin = dataMax = 0.0;
  }
  
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
TCanvas*
AliMUONPainterMatrix::CreateCanvas(Int_t x, Int_t y, Int_t w, Int_t h)
{
  /// Generate a canvas to show the painter matrix.
  ///
  /// Layout is the following :
  ///
  /// ----------------------------------------------------
  /// |    title describing what is plotted              |
  /// ----------------------------------------------------
  /// |                                        |         |
  /// |                                        |         |
  /// |                                        |         |
  /// |                                        |         |
  /// |                                        |         |
  /// |             painter themselves         | color   |
  /// |                                        | range   |
  /// |                                        |         |
  /// |                                        |         |
  /// ----------------------------------------------------
  ///
  
  Int_t mw = ( w <= 0 ? TMath::Nint(gClient->GetDisplayWidth()*0.9) : w );
  Int_t mh = ( h <= 0 ? TMath::Nint(gClient->GetDisplayHeight()*0.9) : h );
  
  TString name(GetName());
  
  TCanvas* d = new TCanvas(name.Data(),name.Data(),x,y,mw,mh);

  TVirtualPad* pTitle = new TPad(Form("%s-title",name.Data()),Form("%s-title",name.Data()),0,0.9,1.0,0.99);
  
  pTitle->Draw();
  
  pTitle->cd();
  
  TPaveLabel* text = new TPaveLabel(0,0,1,1,"");
  text->SetFillStyle(0);
  text->SetFillColor(0);
  text->SetTextColor(4);
  text->SetBorderSize(0);
  
  text->SetLabel(name.Data());
  
  text->Draw();
  
  d->cd();
  
  TVirtualPad* pMatrix = new TPad(Form("%s-matrix",name.Data()),Form("%s-matrix",name.Data()),0,0,0.9,0.89);
  
  pMatrix->Draw();
  pMatrix->cd();
  
  Draw();
  
  d->cd();
  
  TVirtualPad* pColor = new TPad(Form("%s-color",name.Data()),Form("%s-color",name.Data()),0.91,0.01,0.99,0.89);
    
  pColor->Range(0,0,1,1);

  pColor->Draw();
  
  pColor->cd();
  
  Int_t ndivisions(20);
  
  Double_t rangeXmin(0.1);
  Double_t rangeXmax(0.9);
  
  Double_t ymin, ymax;
  
  GetDataRange(ymin,ymax);
    
  Double_t min(0.0);
  Double_t max(1.0);
  
  Double_t step = (max-min)/ndivisions;

  Double_t hsize = 1.0/(ndivisions+2);

  Double_t ypos = 1.0;
  
  for ( Int_t i = -1; i < ndivisions+1; ++i ) 
  {
    Double_t value = max - (min + step*i);
    
    Int_t color = AliMUONPainterHelper::Instance()->ColorFromValue(value,min,max);
    
    Bool_t limit(kFALSE);
    
    TString label;
    TString sign;
    
    Double_t yvalue(0.0);
    
    if ( i == -1 )
    {
      yvalue = ymax;
      limit = kTRUE;
      sign = ">";
    }
    else if ( i == ndivisions )
    {
      yvalue = ymin;
      limit = kTRUE;
      sign = "<=";
    }
    
    if (limit)
    {
      if ( TMath::Abs(yvalue) < 1E5 ) 
      {
        label = Form("%s %7.2f",sign.Data(),yvalue);    
      }
      else
      {
        label = Form("%s %e",sign.Data(),yvalue);
      }
    }

    TPaveLabel* box = new TPaveLabel(rangeXmin,TMath::Max(0.001,ypos-hsize),rangeXmax,ypos,label.Data(),"");    
    
    ypos -= hsize;
    
    box->SetFillColor(color);
    box->SetTextColor( i == -1 ? 0 : 1 );
    box->SetBorderSize(1);
    box->SetLineColor(1);
    box->Draw();
  }  
  
  d->SetEditable(kFALSE);
  
  return d;
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
AliMUONPainterMatrix::Draw(Option_t*)
{
  /// Append our painters to the current pad

  if (!gPad)
  {
    gROOT->MakeDefCanvas();
  }

  TVirtualPad* pad = gPad;

  gPad->Divide(Nx(),Ny());

  for ( Int_t i = 0; i < Size(); ++i )
  {
    AliMUONVPainter* painter = Painter(i);
    pad->cd(i+1);
    painter->Draw("R");
  }

  AppendPad("");
}

//_____________________________________________________________________________
std::string
AliMUONPainterMatrix::NameIt(const char* whatname, const char* basename, const AliMUONAttPainter& att)
{
  /// Build a name
  if ( strlen(whatname) > 0 )
  {
    return Form("%s-%s-%s",whatname,basename,att.GetName());
  }
  else
  {
    return Form("nodata-%s-%s",basename,att.GetName());
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
void
AliMUONPainterMatrix::Print(Option_t*) const
{
  /// Printout
  std::cout << "Whatname=" << fWhatname.Data() << " Basename=" << fBasename.Data()
  << " Nx=" << fNx << " Ny=" << fNy << " Att=" << fAttributes.GetName() << std::endl;
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
  
  if ( d ) 
  {
    fWhatname = Form("%s-%s",d->GetName(),d->DimensionName(indexInData).Data());
  }
  else
  {
    fWhatname = "";
  }
  
  SetName();
  
  // get the data source range, if any is given
  
  if ( d )
  {
    AliMUONPainterEnv* env = AliMUONPainterHelper::Instance()->Env();
    
    TString desc = env->DataSourceDescriptor(d->GetName());
    TString dimensionName = d->DimensionName(DataIndex());
    
    Double_t xmin,xmax;

    Bool_t ok = env->Ranges2DimensionRange(env->Descriptor2Ranges(desc),dimensionName.Data(),xmin,xmax);
    
    if (ok)
    {
      SetDataRange(xmin,xmax);
      env->Save();
    }
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
  
  AliMUONVTrackerData* data = Data();
  
  if ( data )
  {
    AliMUONPainterEnv* env = AliMUONPainterHelper::Instance()->Env();
    
    TString dimensionName = data->DimensionName(DataIndex());
    env->SetDimensionRange(data->GetName(),dimensionName,dataMin,dataMax);
  }
  
}

//_____________________________________________________________________________
void AliMUONPainterMatrix::SetName()
{
	/// Build our name
//	fName = NameIt(fWhatname.Data(),fBasename.Data(),fAttributes);
	fName = "nodata";

	if ( fWhatname.Length() > 0 )
	{
		fName = fWhatname;
	}
	fName += "-";
	fName += fBasename;
	fName += "-";
	fName += fAttributes.GetName();
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
Int_t
AliMUONPainterMatrix::Size() const
{
  /// Return the number of painters we actually handle
  return fPainters->GetLast()+1;
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

  SetName();
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
