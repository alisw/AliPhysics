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

#include "AliMUONVPainter.h"

#include "AliLog.h"
#include "AliMUON2DMap.h"
#include "AliMUONCalibParamND.h"
#include "AliMUONContour.h"
#include "AliMUONContourPainter.h"
#include "AliMUONObjectPair.h"
#include "AliMUONPainterGroup.h"
#include "AliMUONPainterHelper.h"
#include "AliMUONPainterDataRegistry.h"
#include "AliMUONTrackerDataHistogrammer.h"
#include "AliMUONVTrackerData.h"
#include "AliMpManuUID.h"
#include <Riostream.h>
#include <TCanvas.h>
#include <TClass.h>
#include <TClassMenuItem.h>
#include <TH1.h>
#include <TList.h>
#include <TMap.h>
#include <TMath.h>
#include <TMethodCall.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TVirtualPad.h>
#include <cassert>
#include <float.h>

/// \class AliMUONVPainter
///
/// Base class for a graphical object representing some part of the
/// MUON tracking system.
/// 
/// A painter is a graphical representation of some part (e.g. detection element,
/// full chamber, one manu, etc...) of the MUON tracking system.
///
/// A painter is a double fold hierarchical structure. 
///
/// First, a painter is part of a tree (mother->childrens), that describe
/// the natural organization of the spectrometer. For instance, a chamber
/// painter has children that are the detection element, the detection elements
/// themselves contain manus, which in turn contain channels.
///
/// Second, a painter contains a number of "painter groups" (see AliMUONPainterGroup). 
/// A group gather all the painters of the same type, 
/// where the type is a string identifying which part of system we're dealing
/// with (chamber, DE, manu, etc...)
///
/// The groups are there to ease the manipulation of similar painters, e.g. if
/// we want to hide all detection elements, we hide the "detection element group"
/// Some special groups are the responder and the plotter groups. The responder
/// group is the group which is currently responding to mouse events. 
/// The plotter group is the group which is supposed to represent some data.
/// 
/// There are two ways to represent the painter on screen. In any case, we can
/// outline the painter (i.e. draw its borders) (see AliMUONVPainter::PaintOutline).
/// In the cases where the painter is attached to some data source (i.e. it is
/// used to represent some data about its type, e.g. the mean charge on some manu),
/// we can draw the full area of the contour, using some color (see
/// AliMUONVPainter::PaintArea).
///
/// Note that you can outline several types of painters (aka groups) at the same
/// time, but you cannot plot several groups at the same time.
///
/// Painters are TQObject so they can emit signals.
///
/// Currently emitted signal are : 
///
/// void Clicked(AliMUONVPainter* painter, Double_t*);
/// DoubleClicked(AliMUONVPainter* painter, Double_t*);
///
/// to know which and where a painter was (double-) clicked.
///
/// \author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONVPainter)
///\endcond

//_____________________________________________________________________________
AliMUONVPainter::AliMUONVPainter(TRootIOCtor*) : TObject(), 
TQObject(),
fHistogram(0x0),
fPainterGroups(0x0),
fResponderGroup(0x0),
fName(""),
fPathName(""),
fType(""),
fMother(0x0),
fGroup(0x0),
fContour(0x0),
fChildren(0x0),
fPlotterGroup(0x0),
fBorderFactor(1.1),
fPad(0x0),
fAttributes(),
fLineColor(1),
fLineWidth(1),
fIsValid(kTRUE)
{
  /// streamer ctor
}

//_____________________________________________________________________________
AliMUONVPainter::AliMUONVPainter(const char* type)
: TObject(), 
  TQObject(),
  fHistogram(0x0),
  fPainterGroups(0x0),
  fResponderGroup(0x0),
  fName(""),
  fPathName(""),
  fType(type),
  fMother(0x0),
  fGroup(0x0),
  fContour(0x0),
  fChildren(0x0),
  fPlotterGroup(0x0),
  fBorderFactor(1.1),
  fPad(0x0),
  fAttributes(),
  fLineColor(1),
  fLineWidth(1),
  fIsValid(kTRUE)
{
    /// ctor
    SetID(-1,-1);
}

//_____________________________________________________________________________
AliMUONVPainter::AliMUONVPainter(const AliMUONVPainter& rhs)
: TObject(rhs),
TQObject(),
fHistogram(0x0),
fPainterGroups(0x0),
fResponderGroup(0x0),
fName(""),
fPathName(""),
fType(""),
fMother(0x0),
fGroup(0x0),
fContour(0x0),
fChildren(0x0),
fPlotterGroup(0x0),
fBorderFactor(1.0),
fPad(0x0),
fAttributes(),
fLineColor(-1),
fLineWidth(-1),
fIsValid(kTRUE)
{
  /// copy ctor
  rhs.Copy(*this);
}

//_____________________________________________________________________________
AliMUONVPainter& 
AliMUONVPainter::operator=(const AliMUONVPainter& rhs)
{
  /// assignment operator
  if ( this != &rhs ) 
  {
    rhs.Copy(*this);
  }
  return *this;
}

//_____________________________________________________________________________
AliMUONVPainter::~AliMUONVPainter()
{
  /// dtor
  delete fChildren;
  delete fHistogram;
}

//_____________________________________________________________________________
AliMpArea
AliMUONVPainter::Area() const
{
  /// Return the area covered by this painter
  if ( fContour ) 
  {
    return fContour->Area();
  }
  else
  {
    AliWarning("Returning an invalid area, as contour is not defined");
    return AliMpArea();
  }
}

//_____________________________________________________________________________
void
AliMUONVPainter::Add(AliMUONVPainter* painter)
{
  /// Add a child painter
  if (!fChildren) fChildren = new TObjArray;
  assert(painter->Mother()==0x0);
  fChildren->Add(painter);
  painter->SetMother(this);
}

//_____________________________________________________________________________
TCollection*
AliMUONVPainter::Children() const
{
  /// Return the list of childrens
  return fChildren;
}

//_____________________________________________________________________________
void
AliMUONVPainter::Clicked(AliMUONVPainter* painter, Double_t* values)
{
  /// Let our mother emit the signal as clients are probably connected to
  /// our (grand)mother, not to us

  if ( Mother() ) 
  {
    Mother()->Clicked(painter,values);
  }
  else
  {
    Long_t param[] = { (Long_t)painter,(Long_t)values };
  
    Emit("Clicked(AliMUONVPainter*,Double_t*)",param);
  }
}

//_____________________________________________________________________________
void
AliMUONVPainter::ShiftClicked(AliMUONVPainter* painter, Double_t* values)
{
  /// Let our mother emit the signal as clients are probably connected to
  /// our (grand)mother, not to us
  
  if ( Mother() ) 
  {
    Mother()->ShiftClicked(painter,values);
  }
  else
  {
    Long_t param[] = { (Long_t)painter,(Long_t)values };
    
    Emit("ShiftClicked(AliMUONVPainter*,Double_t*)",param);
  }
}

//_____________________________________________________________________________
void
AliMUONVPainter::ComputeDataRange(const AliMUONVTrackerData&, Int_t,
                                  Double_t&, Double_t&) const
{
  /// Should compute the min and max of a given data source
  AliError("Not implemented. Please fixe me");
}

//_____________________________________________________________________________
TString
AliMUONVPainter::ContourName() const
{
  /// Default implementation of the contour name.
  
  TString name(PathName());

  name += "-";
  name += fAttributes.Name();
  
  return name;
}

//_____________________________________________________________________________
void
AliMUONVPainter::Copy(TObject& object) const
{
  /// Copy this to object.
  
  TObject::Copy(object);

  AliMUONVPainter& painter = static_cast<AliMUONVPainter&>(object);

  painter.fType = fType;
  painter.fName = fName;
  painter.fPathName = fPathName;
  
  painter.fMother = 0x0;
  painter.fContour = fContour;
  
  painter.fGroup = 0x0;
  painter.fResponderGroup = 0x0;
  painter.fPlotterGroup = 0x0;
  
  painter.fBorderFactor = fBorderFactor;

  painter.fAttributes = fAttributes;
  
  painter.fAttributes.SetCathodeAndPlaneDisabled(kFALSE);
  
  painter.fPad = fPad;
  
  painter.fLineColor = fLineColor;
  painter.fLineWidth = fLineWidth;
  
  painter.fIsValid = fIsValid;
  
  delete painter.fChildren;
  painter.fChildren = 0x0;
  
  painter.fID[0] = fID[0];
  painter.fID[1] = fID[1];
  
  delete painter.fHistogram;
  painter.fHistogram = 0x0;
  
  TIter next(fChildren);
  AliMUONVPainter* p;
  
  while ( ( p = static_cast<AliMUONVPainter*>(next()) ) )
  {
    painter.Add(static_cast<AliMUONVPainter*>(p->Clone()));
  }
    
  painter.UpdateGroupsFrom(*this);
  
  object.ResetBit(kCanDelete);
}

//_____________________________________________________________________________
AliMUONPainterGroup*
AliMUONVPainter::CreateGroup(const char* type, Int_t depth)
{
  /// Create a painter group at a given depth
  
  if (!fPainterGroups) fPainterGroups = new TMap;
  TObject* o = fPainterGroups->GetValue(type);
  if (o)
  {
    AliError(Form("Group %s is already there ! Check this",type));
    return 0x0;
  }
  AliMUONPainterGroup* group = new AliMUONPainterGroup(type,depth);
  fPainterGroups->Add(new TObjString(type),group);
  return group;
}

//_____________________________________________________________________________
void
AliMUONVPainter::CreateGroups()
{
  /// Groups our children into groups
  
  if ( Mother() ) 
  {
    AliFatal("Not supposed to create groups for a children");
  }
  
  TList list;
  FlatList(list);
  
  TIter next(&list);
  AliMUONVPainter* painter;
  
  while ( ( painter = static_cast<AliMUONVPainter*>(next()) ) )
  {
    AliMUONPainterGroup* group = Group(painter->Type());
    if (!group) 
    {
      group = CreateGroup(painter->Type(),painter->Depth());
    }
    group->Add(painter);
  }
}

//_____________________________________________________________________________
AliMUONVPainter*
AliMUONVPainter::Detach() const
{
  /// Make this a new top painter (i.e. a master)
  
  AliDebug(1,Form("Detaching %s",GetName()));
           
  AliMUONVPainter* p = static_cast<AliMUONVPainter*>(Clone());
  
  AliMUONVPainter* master = Master();
  
  if ( master )
  {
    AliDebug(1,Form("UpdatingGroups of the detached painter %s from its master %s",
                    p->GetName(),master->GetName()));
    p->UpdateGroupsFrom(*master);
  }
  
  return p;
}

//_____________________________________________________________________________
Int_t
AliMUONVPainter::Depth() const
{
  /// Return our depth in the hierarchy
  
  if ( Mother() ) 
  {
    return Mother()->Depth() + 1;
  }
  else
  {
    return 0;
  }
}

//_____________________________________________________________________________
Int_t	
AliMUONVPainter::DistancetoPrimitive(Int_t px, Int_t py)
{
  /// See TObject::DistancetoPrimitive
  
  static const Int_t kBigValue = 999999;

  if (!gPad) return kBigValue;
  
  Double_t x,y;
  
  AliMUONVPainter* painter = GetPainter(px,py,x,y);
  
  x=y=0.0; // to avoid compiler warning
  
  if ( painter == this) return 0;
  
  return kBigValue;
}

//_____________________________________________________________________________
void
AliMUONVPainter::DoubleClicked(AliMUONVPainter*, Double_t*)
{
  /// Should emit the DoubleClicked signal (if I knew how to detect those events...)
  
  AliWarning("Please implement me !");

  //  if ( fMother )
//  {
//    // let our top mother emit the signal as clients are probably connected to
//    // our mother, not to us
//    Top()->DoubleClicked(painter,values);
//  }
//  else
//  {
//    Long_t param[] = { (Long_t)painter,(Long_t)values };
//    
//    Emit("DoubleClicked(AliMUONVPainter*,Double_t*)",param);
//  }
}

//_____________________________________________________________________________
void
AliMUONVPainter::Draw(Option_t* opt)
{
  /// Append ourselves to the current pad
  
  if (!gPad) 
  {
    gROOT->MakeDefCanvas();
  }
  
  Bool_t kMustSetRange(kFALSE);
 
  TString sopt(opt);
  sopt.ToUpper();
 
  if (sopt.Contains("R") ) kMustSetRange=kTRUE;
  
  if (kMustSetRange)
  {
    Double_t x1,y1,x2,y2;
    GetBoundingBox(x1,y1,x2,y2);
    if ( gPad) gPad->Range(x1,y1,x2,y2);
  }
 
  if ( !fMother && !fPainterGroups ) 
  {
    CreateGroups();
  }
  
  TIter next(fChildren);
  AliMUONVPainter* painter;
  while ( ( painter = static_cast<AliMUONVPainter*>(next()) ) )
  {
    painter->Draw();
  }  
  
  AppendPad(opt);
  
  fPad = gPad;
}

//_____________________________________________________________________________
void 
AliMUONVPainter::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
  /// Handle graphics events
  
  Double_t x,y;
    
  AliMUONVPainter* painter = GetPainter(px,py,x,y);

  if ( painter == this ) 
  {
    Double_t values[] = { x,y };
  
    switch (event)
    {
      case kButton2Up:
        ShiftClicked(this,values);
        break;
      case kButton1Up:
        Clicked(this,values);      
        break;
      case kButton1Double:
        //the following statement is required against other loop executions before returning (depending on the time between the clicks) 
        gPad->GetCanvas()->HandleInput((EEventType)-1,0,0); 
        DoubleClicked(this,values);
        break;
    }
  }
}

//_____________________________________________________________________________
void
AliMUONVPainter::FlatList(TList& list)
{
  /// Make a flat list of our children (and ourselves)
  
  TIter next(fChildren);
  AliMUONVPainter* painter;
  while ( ( painter = static_cast<AliMUONVPainter*>(next())))
  {
    painter->FlatList(list);
  }
  
  list.Add(this);
}

//_____________________________________________________________________________
void 
AliMUONVPainter::GetBoundingBox(Double_t& x1, Double_t& y1, 
                                Double_t& x2, Double_t& y2) const
{
  /// Get the bounding box = our area
  AliMpArea area(Area().GetPositionX(),
                 Area().GetPositionY(),
                 Area().GetDimensionX()*fBorderFactor,
                 Area().GetDimensionY()*fBorderFactor);

  x1 = area.LeftBorder();
  y1 = area.DownBorder();
  x2 = area.RightBorder();
  y2 = area.UpBorder();
}

//_____________________________________________________________________________
char*	
AliMUONVPainter::GetObjectInfo(Int_t, Int_t) const
{
  /// See TObject::GetObjectInfo
  return const_cast<char*>(GetName());
}

//_____________________________________________________________________________
AliMUONVPainter* 
AliMUONVPainter::GetPainter(Int_t px, Int_t py, Double_t& x, Double_t& y) const
{
  /// Get the responder painter at integer position (px,py), and get back its
  /// absolute position (x,y)
  
  PixelToPad(px,py,x,y);
  
  if ( !IsInside(x,y) ) return 0x0;
  
  if ( fGroup->IsResponder() ) return const_cast<AliMUONVPainter*>(this);
  
  if (fChildren)
  {
    TIter next(fChildren);
    AliMUONVPainter* painter;
    
    while ( ( painter = static_cast<AliMUONVPainter*>(next()) ) )
    {
      AliMUONVPainter* p = painter->GetPainter(px,py,x,y);
      if (p) return p;
    }
  }  
  
  return 0x0;
}

//_____________________________________________________________________________
void
AliMUONVPainter::GetTypes(TObjArray& types) const
{
  /// Get the list of types (as a TObjArray of TObjString) 
  /// of our hierarchy, sorted alphabetically
  
  types.SetOwner(kTRUE);
  types.Clear();

  TObjArray tmp;
  tmp.SetOwner(kFALSE);
  
  TIter next(fPainterGroups);
  
  TObjString* str;
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    AliMUONPainterGroup* group = Group(str->String().Data());
    tmp.AddLast(group);
  }
  
  tmp.Sort();
  
  Int_t n = tmp.GetLast()+1;
  
  Int_t* index = new Int_t[n];
  
  Int_t* a = new Int_t[n];
  
  for ( Int_t i = 0; i < n; ++i )
  {
    AliMUONPainterGroup* group = static_cast<AliMUONPainterGroup*>(tmp.At(i));
    a[i] = group->Depth();
  }
  
  TMath::Sort(n,a,index,kFALSE);
  
  for ( Int_t i = 0; i < n; ++i ) 
  {
    AliMUONPainterGroup* group = static_cast<AliMUONPainterGroup*>(tmp.At(index[i]));
    types.AddLast(new TObjString(group->Type()));
  }
  
  delete[] index;
  delete[] a;
}

//_____________________________________________________________________________
AliMUONPainterGroup*
AliMUONVPainter::Group(const char* type) const
{
  /// Returns a group of a given type
  if (!fPainterGroups) return 0x0;
  return static_cast<AliMUONPainterGroup*>(fPainterGroups->GetValue(type));
}

//_____________________________________________________________________________
AliMUONPainterGroup*
AliMUONVPainter::Group(Int_t depth) const
{
  /// Returns a group of a given depth
  if (!fPainterGroups) return 0x0;
  TIter next(fPainterGroups);
  TObjString* groupName;
  while ( ( groupName = static_cast<TObjString*>(next()) ) )
  {
    AliMUONPainterGroup* group = static_cast<AliMUONPainterGroup*>
    (fPainterGroups->GetValue(groupName->String().Data()));
    if ( group->Depth() == depth ) 
    {
      return group;
    }
  }
  return 0x0;
}

//_____________________________________________________________________________
Bool_t 
AliMUONVPainter::IsInside(Double_t x, Double_t y) const
{
  /// Whether point (x,y) is inside our contour
  if (!fContour) return kFALSE;
  return fContour->IsInside(x,y);
}

//_____________________________________________________________________________
Bool_t 
AliMUONVPainter::IsResponder() const
{
  /// Whether we're responding to mouse events
  return MotherGroup()->IsResponder();
}

//_____________________________________________________________________________
AliMUONVPainter*
AliMUONVPainter::Master() const
{
  /// Return the top of the hierarchy
  
  /// if we get no mother, we are the master
  
  if ( Mother() == 0x0 ) return const_cast<AliMUONVPainter*>(this);
  
  AliMUONVPainter* p = Mother();
  
  while ( p->Mother() )
  {
    p = p->Mother();
  }
  
  return p;
}

//_____________________________________________________________________________
void
AliMUONVPainter::Paint(Option_t*)
{
  /// Paint ourselves of screen
  /// If we have some data (i.e. we're belonging to the plotter group)
  /// we use PaintArea.
  /// And if must be outlined, then we do that too.
  
  if ( !MotherGroup()->IsVisible() ) return;

  if ( MotherGroup()->IsPlotter() ) 
  {
    PaintArea(*(MotherGroup()->Data()),
              MotherGroup()->DataIndex(),
              MotherGroup()->DataMin(),
              MotherGroup()->DataMax());
  }
  
  if ( MotherGroup()->IsOutlined() )
  {
    PaintOutline();
  }
  
  if ( IsExcluded() )
  {
    AliMUONContourPainter::Paint(*fContour,1,1,2); // red fill with black thin outline
  }
}

//_____________________________________________________________________________
TString
AliMUONVPainter::Describe(const AliMUONVTrackerData&, Int_t, Double_t, Double_t)
{
  /// Default implementation (must be overriden)
  AliError(Form("%s : implement me",GetName()));
  return "";
}

//_____________________________________________________________________________
void
AliMUONVPainter::PaintArea(const AliMUONVTrackerData&, Int_t, Double_t, Double_t)
{
  /// Default implementation (must be overriden)
  AliError(Form("%s : implement me",GetName()));
  return;
}

//_____________________________________________________________________________
void
AliMUONVPainter::PaintOutline(Int_t color, Int_t width, Double_t /*x*/, Double_t /*y*/)
{
  /// Default implementation is simply a drawing of the contour lines,
  /// not using the optional (x,y)
  Int_t c = color >= 0 ? color : GetLineColor();
  Int_t w = width >= 0 ? width : GetLineWidth();
  
  AliMUONContourPainter::Paint(*fContour,c,w);
}

//_____________________________________________________________________________
void 
AliMUONVPainter::PixelToPad(Int_t px, Int_t py, Double_t& x, Double_t& y)
{
  /// convert (px,py) into pad position (x,y)
  
  x = gPad->PadtoX(gPad->AbsPixeltoX(px));
  y = gPad->PadtoY(gPad->AbsPixeltoY(py));
}

//_____________________________________________________________________________
void
AliMUONVPainter::Print(Option_t* opt) const
{
  /// Printout
  for ( Int_t i = 0; i < Depth()*4; ++i ) 
  {
    cout << " ";
  }
  
  if ( !IsValid() ) cout << "!!!INVALID!!!" << endl;
  
  cout << Form("%p Name %s Depth %d ContourName %s ID=(%d,%d)",
               this,GetName(),Depth(),ContourName().Data(),ID0(),ID1());
  
  if ( fResponderGroup )
  {
    cout << Form(" Responder group %p %s",fResponderGroup,fResponderGroup->Type());
  }
  if ( fPlotterGroup )
  {
    cout << Form(" Plotter group %p %s",fPlotterGroup,fPlotterGroup->Type());
  }
  if ( Mother() )
  {
    cout << Form(" Mother %p %s",Mother(),Mother()->GetName());
  }
  if ( MotherGroup() )
  {
    cout << Form(" Group %p %s ",MotherGroup(),MotherGroup()->Type());
  }
  
  if ( fChildren ) 
  {
    cout << Form(" %d children",fChildren->GetLast()+1);
  }
  
  cout << endl;
  
  TString sopt(opt);
  sopt.ToUpper();
  
  if ( fChildren && ( sopt == "FULL" || sopt == "CHILD" ) ) 
  {
    TIter next(fChildren);
    AliMUONVPainter* painter;
    while ( ( painter = static_cast<AliMUONVPainter*>(next()) ) )
    {
      painter->Print(opt);
    }
  }
  
  if ( fPainterGroups && ( sopt == "FULL" || sopt == "GROUP" ) )
  {
    TIter next(fPainterGroups);
    TObjString* groupName;
    while ( ( groupName = static_cast<TObjString*>(next()) ) )
    {
      AliMUONPainterGroup* group = Group(groupName->String().Data());
      group->Print(opt);
    }
  }
}

//_____________________________________________________________________________
void 
AliMUONVPainter::SetAttributes(const AliMUONAttPainter& attributes)
{
  /// Set our attributes  
  fAttributes = attributes;
}

//_____________________________________________________________________________
void 
AliMUONVPainter::SetContour(AliMUONContour* contour)
{
  /// Set out contour
  if (!contour)
  {
    AliError(Form("Setting a null contour for painter %s : bad idea !",PathName().Data()));
  }
  fContour = contour;
}

//_____________________________________________________________________________
void
AliMUONVPainter::SetData(const char* pattern, AliMUONVTrackerData* data,
                         Int_t dataIndex)
{
  /// Tell all painters which type matches pattern that they should
  /// monitor a given data source
  
  if ( !fPainterGroups ) 
  {
    CreateGroups();
  }
  
  if ( data ) 
  {
    data->Connect("Destroyed()",ClassName(),this,Form("SetData(=\"%s\",0x0,-1)",pattern));
  }
  
  TIter next(fPainterGroups);
  TObjString* str;
  
  fPlotterGroup = 0x0;
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    AliMUONPainterGroup* group = static_cast<AliMUONPainterGroup*>(fPainterGroups->GetValue(str));
        
    if ( group->Matches(pattern) )
    {
      group->SetData(data,dataIndex);
      if ( data ) 
      {        
        fPlotterGroup = group;
      }
    }
    else
    {
      group->SetData(0x0,-1);
    }
  }
  
  // Update context menus
  TList list;
  FlatList(list);
  
  TIter pnext(&list);
  AliMUONVPainter* p;
  
  AliMUONPainterGroup* group = Master()->PlotterGroup();
  
  while ( ( p = static_cast<AliMUONVPainter*>(pnext()) ) )
  {
    TList* l = p->IsA()->GetMenuList();
  
    l->Delete();
  
    TClassMenuItem* n(0x0);
    
    l->Add(new TClassMenuItem(TClassMenuItem::kPopupUserFunction,p->IsA(),
                              "Include","Include",p,"",-1,kTRUE));
    l->Add(new TClassMenuItem(TClassMenuItem::kPopupUserFunction,p->IsA(),
                              "Exclude","Exclude",p,"",-1,kTRUE));    
    
    if ( group )  
    {
      if ( data && data->IsHistogrammed(0) ) 
      {
        // Add histo drawing to the popup menu
        TString name("Draw histogram of ");
        
        name += data->ExternalDimensionName(0);
        
        n = new TClassMenuItem(TClassMenuItem::kPopupUserFunction,p->IsA(),
                               name.Data(),"DrawHistogram0",p,"",-1,kTRUE);
        l->Add(n);
        
        name += " clone";
        
        n = new TClassMenuItem(TClassMenuItem::kPopupUserFunction,p->IsA(),
                               name.Data(),"DrawHistogramClone0",p,"",-1,kTRUE);
        l->Add(n);
      }
      
      Int_t nd = data->IsSingleEvent() ? data->ExternalDimension() : data->ExternalDimension()*2;
      
      for ( Int_t i = 0; i < nd; ++i ) 
      {
        n = new TClassMenuItem(TClassMenuItem::kPopupUserFunction,p->IsA(),
                               Form("Draw %s clone",data->DimensionName(i).Data()),
                               Form("DrawInternalHistogramClone%d",i),p,"",-1,kTRUE);
        l->Add(n);
      } 
    }
  }
}

//_____________________________________________________________________________
void
AliMUONVPainter::DrawInternalHistogram(Int_t dim) const
{
  /// Draw histogram (and delete the previous one)
  
  delete fHistogram;
  fHistogram = 0x0;
  
  DrawInternalHistogramClone(dim);
}

//_____________________________________________________________________________
void
AliMUONVPainter::DrawInternalHistogramClone(Int_t dim) const
{
  /// Draw histogram 
  
  fHistogram = AliMUONTrackerDataHistogrammer::CreateHisto(*this,-1,dim);
  
  if (fHistogram) 
  {
    new TCanvas();
    fHistogram->Draw();
  }
}

//_____________________________________________________________________________
void
AliMUONVPainter::DrawHistogram(Double_t* values) const
{
  /// Draw histogram (and delete the previous one)

  delete fHistogram;
  fHistogram = 0x0;
  
  DrawHistogramClone(values);
}

//_____________________________________________________________________________
void
AliMUONVPainter::DrawHistogramClone(Double_t*) const
{
  /// Draw histogram 
  
  fHistogram = AliMUONTrackerDataHistogrammer::CreateHisto(*this,0,-1);
  
  if (fHistogram) 
  {
    new TCanvas();
    fHistogram->Draw();
  }
}

//_____________________________________________________________________________
void
AliMUONVPainter::FillManuList(TObjArray& manuList) const
{
  /// Append to manulist
  /// This is the default implementation, which just calls the FillManuList
  /// of all our children.
  /// Some derived class might need to override this in order to exclude
  /// some children from the fill.
  
  TIter next(Children());
  
  AliMUONVPainter* p;
  
  while ( ( p = static_cast<AliMUONVPainter*>(next()) ) )
  {
    p->FillManuList(manuList);
  }
}

//_____________________________________________________________________________
void
AliMUONVPainter::SetLine(Int_t depth, Int_t lineColor, Int_t lineWidth)
{
  /// Set the line attributes of painters at a given depth
  AliMUONPainterGroup* group = Group(depth);
  if ( group )
  {
    group->SetLine(lineColor,lineWidth);
  }
}

//_____________________________________________________________________________
void
AliMUONVPainter::SetMother(AliMUONVPainter* painter)
{
  /// Set our mother
  fMother = painter;
}

//_____________________________________________________________________________
void
AliMUONVPainter::SetOutlined(const char* pattern, Bool_t flag)
{
  /// Decide whether or not painters which type matches pattern 
  /// should be outlined
  
  AliDebug(1,Form("pattern=%s flag=%d",pattern,flag));
  
  if (!fPainterGroups)
  {
    CreateGroups();
  }
  
  TIter next(fPainterGroups);
  TObjString* str;
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    AliMUONPainterGroup* group = static_cast<AliMUONPainterGroup*>(fPainterGroups->GetValue(str));
    if ( group->Matches(pattern) )
    {
      group->SetOutlined(flag);
    }
  }
}

//_____________________________________________________________________________
void
AliMUONVPainter::SetResponder(const char* pattern)
{
  /// Set the painters matching pattern to be the responder
  
  AliDebug(1,Form("pattern=%s",pattern));
  
  if (!fPainterGroups)
  {
    CreateGroups();
  }
  
  TIter next(fPainterGroups);
  TObjString* str;
  
  fResponderGroup = 0x0;
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    AliMUONPainterGroup* group = static_cast<AliMUONPainterGroup*>(fPainterGroups->GetValue(str));
    if ( group->Matches(pattern) )
    {
      AliDebug(1,Form("group %s is matching pattern %s : setting to responder",
                      group->Type(),pattern));
      group->SetResponder(kTRUE);
      fResponderGroup = group;
    }
    else
    {
      group->SetResponder(kFALSE);
    }
  }
}

//_____________________________________________________________________________
void
AliMUONVPainter::SetResponder(Int_t depth)
{
  /// Select as responder the *first* group that has a given depth
  
  AliDebug(1,Form("depth=%d",depth));
  
  if (!fPainterGroups)
  {
    CreateGroups();
  }
  
  TIter next(fPainterGroups);
  TObjString* str;
  
  fResponderGroup = 0x0;
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    AliMUONPainterGroup* group = static_cast<AliMUONPainterGroup*>(fPainterGroups->GetValue(str));
    if ( group->Depth() == depth ) 
    {
      AliDebug(1,Form("group %s has correct depth = %d, using as responder",
                      group->Type(),depth));
      group->SetResponder(kTRUE);
      fResponderGroup = group;
      break;
    }
    else
    {
      group->SetResponder(kFALSE);
    }
  }
}

//_____________________________________________________________________________
void
AliMUONVPainter::SetVisible(const char* pattern, Bool_t flag)
{
  /// Decide whether the painters matching pattern should be visible or not
  AliDebug(1,Form("pattern=%s flag=%d",pattern,flag));
  
  if (!fPainterGroups)
  {
    CreateGroups();
  }
  
  TIter next(fPainterGroups);
  TObjString* str;
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    AliMUONPainterGroup* group = static_cast<AliMUONPainterGroup*>(fPainterGroups->GetValue(str));
    if ( group->Matches(pattern) )
    {
      group->SetVisible(flag);
    }
  }
}

//_____________________________________________________________________________
void
AliMUONVPainter::UpdateGroupsFrom(const AliMUONVPainter& painter)
{
  /// (re)Create groups
  delete fPainterGroups;
  fPainterGroups = 0x0;
  
  CreateGroups();
  
  // and copy the status of responder, plotter and visible
  if ( painter.ResponderGroup() ) 
  {
    SetResponder(painter.ResponderGroup()->Type());
  }
  
  if ( painter.PlotterGroup() ) 
  {
    SetData(painter.PlotterGroup()->Type(),
            painter.PlotterGroup()->Data(),
            painter.PlotterGroup()->DataIndex());
    PlotterGroup()->SetDataRange(painter.PlotterGroup()->DataMin(),
                                 painter.PlotterGroup()->DataMax());
  }
  
  TObjArray types;
  painter.GetTypes(types);
  TIter next(&types);
  TObjString* groupName;
  
  while ( ( groupName = static_cast<TObjString*>(next()) ) )
  {
    AliMUONPainterGroup* group = painter.Group(groupName->String().Data());      
    if ( group->IsVisible() ) 
    {
      SetVisible(group->Type(),kTRUE);
    }
    else
    {
      SetVisible(group->Type(),kFALSE);
    }

    if ( group->IsOutlined() ) 
    {
      SetOutlined(group->Type(),kTRUE);
    }
    else
    {
      SetOutlined(group->Type(),kFALSE);
    }
    
    SetLine(group->Depth(),group->GetLineColor(),group->GetLineWidth());
  }
  
}

//_____________________________________________________________________________
void
AliMUONVPainter::Include()
{
  /// Include this painter
  AliInfo(GetName());
  
  /// Update the global interactive read out configuration  
  WriteIROC(1);
}

//_____________________________________________________________________________
void
AliMUONVPainter::GetIROCManuList(TObjArray& manuList)
{
  /// Get the list of manus spanned by this painter AND by its dual

  FillManuList(manuList);
  
  // get our dual
  AliMUONAttPainter att(Attributes());
  
  att.Invert();
  
  att.SetCathodeAndPlaneDisabled(kTRUE);
  
  AliMUONVPainter* p = AliMUONVPainter::CreatePainter(ClassName(),att,ID0(),ID1());
  
  if (p)
  {
    p->FillManuList(manuList);
  }
  
  delete p;
}

//_____________________________________________________________________________
void
AliMUONVPainter::WriteIROC(Double_t value)
{
  /// Update the interactive readout configuration
  
  TObjArray manuList;
  GetIROCManuList(manuList);
  
  AliMpManuUID* muid;
  TIter nextm(&manuList);
  AliMUON2DMap store(true);
  
  while ((muid=static_cast<AliMpManuUID*>(nextm())))
  {
    AliMUONVCalibParam* param = new AliMUONCalibParamND(1,64,
                                                        muid->DetElemId(),
                                                        muid->ManuId(),value);
    store.Add(param);
  }
  
  InteractiveReadOutConfig()->Replace(store);
}

//_____________________________________________________________________________
void
AliMUONVPainter::Exclude()
{
  /// Exclude this painter
  AliInfo(GetName());
  
  /// Update the global interactive read out configuration
  WriteIROC(0.0);
}

//_____________________________________________________________________________
AliMUONVTrackerData*
AliMUONVPainter::InteractiveReadOutConfig() const
{
  /// get the interactive readout config object
  return AliMUONPainterDataRegistry::Instance()->InteractiveReadOutConfig();
}

//_____________________________________________________________________________
AliMUONVPainter* 
AliMUONVPainter::CreatePainter(const char* className, 
                               const AliMUONAttPainter& att, 
                               Int_t id1, Int_t id2)
{
  /// Create a painter (factory method)
  
  TClass* c = TClass::GetClass(className);
  
  if (!c)
  {
    AliErrorClass(Form("Cannot get class %s",className));
    return 0x0;
  }
  
  Int_t n(0);
  
  TMethodCall call;
  
  call.InitWithPrototype(c,className,"AliMUONAttPainter&,Int_t");
  
  if (call.IsValid()) n = 1;
  else
  {
    call.InitWithPrototype(c,className,"AliMUONAttPainter&,Int_t,Int_t");
    
    if ( call.IsValid() ) n = 2;
  }
  
  Long_t returnLong(0x0);
  
  if ( n ==1 ) 
  {
    Long_t params[] = { (Long_t)(&att), (Long_t)(id1) };
    call.SetParamPtrs((void*)(params));
    call.Execute((void*)(0x0),returnLong);
  }
  else if ( n == 2 ) 
  {
    Long_t params[] = { (Long_t)(&att), (Long_t)(id1), (Long_t)(id2) };
    call.SetParamPtrs((void*)(params));
    call.Execute((void*)(0x0),returnLong);
  }
  
  if (!returnLong)
  {
    AliErrorClass(Form("Cannot create a painter of class %s",className));
  }
  
  AliMUONVPainter* rv = reinterpret_cast<AliMUONVPainter*> (returnLong);
  
  if (!rv->IsValid()) 
  {
    AliErrorClass(Form("Painter of class %s is not valid",className));
    delete rv;
    rv = 0x0;
  }
  return rv;
}

//_____________________________________________________________________________
void
AliMUONVPainter::PaintArea(Int_t fillColor)
{
  /// Draw a filled area
  AliMUONContourPainter::Paint(*(Contour()),-1,-1,fillColor);
}
