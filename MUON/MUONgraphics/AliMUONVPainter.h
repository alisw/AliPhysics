#ifndef ALIMUONVPAINTER_H
#define ALIMUONVPAINTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONVPainter
/// \brief Base class for a graphical object representing some part of the
/// MUON tracking system
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ALIMUONATTPAINTER_H
#  include "AliMUONAttPainter.h"
#endif
#ifndef ROOT_TQObject
#  include <TQObject.h>
#endif
#ifndef ROOT_TString
#  include "TString.h"
#endif
#ifndef ROOT_TObject
#  include "TObject.h"
#endif
#ifndef AL_MP_AREA_H
#  include "AliMpArea.h"
#endif
#include <float.h>

class AliMUONContour;
class AliMUONPainterGroup;
class AliMUONVTrackerData;
class AliMpArea;
class TCollection;
class TH1;
class TList;
class TMap;
class TObjArray;
class TVirtualPad;

class AliMUONVPainter : public TObject, public TQObject
{
public:  

  AliMUONVPainter(TRootIOCtor*);
  AliMUONVPainter(const char* type="");
  AliMUONVPainter(const AliMUONVPainter& rhs);
  AliMUONVPainter& operator=(const AliMUONVPainter& rhs);
  virtual ~AliMUONVPainter();

  /// Add a painter to our list of children. We adopt this painter (i.e. we become owner).
  void Add(AliMUONVPainter* painter);
  
  /// Return the area containing this painter
  AliMpArea Area() const;
  
  virtual void SetAttributes(const AliMUONAttPainter& attributes);
  
  /// Convert attributes so they are valid ones for us.
  virtual AliMUONAttPainter Validate(const AliMUONAttPainter& attributes) const { return attributes; }
  
  /// Get our attributes
  const AliMUONAttPainter& Attributes() const { return fAttributes; }
  
  virtual void ComputeDataRange(const AliMUONVTrackerData& data, Int_t dataIndex, 
                                Double_t& dataMin, Double_t& dataMax) const;
  
  virtual void Copy(TObject& object) const;
  
  AliMUONVPainter* Detach() const;
  
  /// Whether this painter can be detached from the current view.
  virtual Bool_t CanBeDetached() const { return kTRUE; }
  
  /// Whether we are valid or not
  virtual Bool_t IsValid() const { return fIsValid; }
  
  /// Mark us as not valid
  void Invalidate() { fIsValid = kFALSE; }
  
  Int_t Depth() const;
  
  virtual Int_t	DistancetoPrimitive(Int_t px, Int_t py);
  
  virtual void Draw(Option_t* opt="");
  
  virtual void ExecuteEvent(Int_t event, Int_t px, Int_t py);
  
  /// Return the contour representing the outline of this object
  AliMUONContour* Contour() const { return fContour; }

  /// Get our name
  virtual const char* GetName() const { return Name().Data(); }
  
  /// Get our name
  virtual TString Name() const { return fName; }
  
  /// Get our path name (aka fullname)
  virtual TString PathName() const { return fPathName; }
  
  virtual TString ContourName() const;

  virtual char*	GetObjectInfo(Int_t px, Int_t py) const;
  
  void GetTypes(TObjArray& types) const;

  /// Return our mother group
  AliMUONPainterGroup* MotherGroup() const { return fGroup; }
  
  /// Return specific name at a given position, if needed.
  virtual TString NameAtPosition(Double_t /*x*/, Double_t /*y*/) const { return GetName(); }
  
  AliMUONPainterGroup* Group(const char* type) const;
  
  AliMUONPainterGroup* Group(Int_t depth) const;
  
  /// Whether we handle mouse motion or not
  virtual Bool_t HandleMouseMotion() const { return kFALSE; }
  
  Bool_t IsResponder() const;
  
  Bool_t IsInside(Double_t x, Double_t y) const;

  /// Return our mother (0 if we're the top node)
  AliMUONVPainter* Mother() const { return fMother; }

  virtual void Paint(Option_t* opt="");

  virtual void Print(Option_t* opt="") const;

  /// Return the plotter group
  AliMUONPainterGroup* PlotterGroup() const { return fPlotterGroup; }
  
  /// Return the responder group
  AliMUONPainterGroup* ResponderGroup() const { return fResponderGroup; }

  /// Set out contour
  void SetContour(AliMUONContour* contour);
  
  void SetData(const char* pattern, AliMUONVTrackerData* data, Int_t dataIndex);

  void SetLine(Int_t depth, Int_t lineColor, Int_t lineWidth);
  
  /// Set our mother group
  void SetMotherGroup(AliMUONPainterGroup* group) { fGroup = group; }
  
  void SetMother(AliMUONVPainter* painter);
  
  void SetOutlined(const char* pattern, Bool_t flag);
  
  virtual void SetResponder(const char* pattern);
  
  virtual void SetResponder(Int_t depth);
  
  void SetVisible(const char* pattern, Bool_t flag);
  
  /// Return our type (e.g. PCB, Chamber, DE, MANU, etc...)
  const char* Type() const { return fType.Data(); }

  // SIGNALS
  
  void Clicked(AliMUONVPainter* painter, Double_t* pos); // *SIGNAL*

  void ShiftClicked(AliMUONVPainter* painter, Double_t* pos); // *SIGNAL*

  void DoubleClicked(AliMUONVPainter* painter, Double_t* pos); // *SIGNAL*
    
  static void PixelToPad(Int_t px, Int_t py, Double_t& x, Double_t& y);

  virtual void PaintOutline(Int_t color=-1, Int_t width=-1, Double_t x=FLT_MAX, Double_t y=FLT_MAX);

  virtual void PaintArea(Int_t fillColor);
  
  virtual void PaintArea(const AliMUONVTrackerData& data, Int_t dataIndex,
                         Double_t min, Double_t max);
    
  /// Get the pad in which we are plotted
  TVirtualPad* Pad() const { return fPad; }
  
  /// Get our line color
  Int_t GetLineColor() const { return fLineColor; }
  
  /// Get our line width
  Int_t GetLineWidth() const { return fLineWidth; }
  
  /// Set our line color
  void SetLineColor(Int_t lineColor) { fLineColor = lineColor; }
  
  /// Set our line width
  void SetLineWidth(Int_t lineWidth) { fLineWidth = lineWidth; }
  
  /// Set our name
  void SetName(const char* name) { fName = name; }
  
  /// Set our path name (aka fullname)
  void SetPathName(const char* pathName) { fPathName = pathName; }
  
  static AliMUONVPainter* CreatePainter(const char* className, 
                                        const AliMUONAttPainter& att,
                                        Int_t id1, Int_t id2);
    
  /// Get our first ID
  Int_t ID0() const { return fID[0]; }
  /// Get our second ID
  Int_t ID1() const { return fID[1]; }
  
  /// Set our IDs
  void SetID(Int_t id0, Int_t id1) { fID[0] = id0; fID[1] = id1; }
  
  virtual TString Describe(const AliMUONVTrackerData& data, Int_t dataIndex,
                           Double_t x=FLT_MAX, Double_t y=FLT_MAX);

  void UpdateGroupsFrom(const AliMUONVPainter& painter);

  AliMUONVPainter* Master() const;
  
  virtual void DrawHistogram(Double_t* values=0x0) const;

  /// To avoid getting a popup asking for the parameter in the GUI...
  void DrawHistogram0() const { DrawHistogram(0x0); }
  /// To avoid getting a popup asking for the parameter in the GUI...
  void DrawHistogramClone0() const { DrawHistogramClone(0x0); }

  virtual void DrawHistogramClone(Double_t* values=0x0) const;
  
  virtual void DrawInternalHistogram(Int_t dim) const;

  virtual void DrawInternalHistogramClone(Int_t dim) const;

  /// Append (i.e. don't have the right to clear the array !) our list of manus to manuList
  virtual void FillManuList(TObjArray& manuList) const;
  
  /// following kind of stupid lines (SL), because I don't know how to
  /// pass parameters to TClassMenuItem for context menu (don't even
  /// know if that's possible at all)
  /// SL
  void DrawInternalHistogramClone0() { DrawInternalHistogramClone(0); }
  /// SL
  void DrawInternalHistogramClone1() { DrawInternalHistogramClone(1); }
  /// SL
  void DrawInternalHistogramClone2() { DrawInternalHistogramClone(2); }
  /// SL
  void DrawInternalHistogramClone3() { DrawInternalHistogramClone(3); }
  /// SL
  void DrawInternalHistogramClone4() { DrawInternalHistogramClone(4); }
  /// SL
  void DrawInternalHistogramClone5() { DrawInternalHistogramClone(5); }
  /// SL
  void DrawInternalHistogramClone6() { DrawInternalHistogramClone(6); }
  /// SL
  void DrawInternalHistogramClone7() { DrawInternalHistogramClone(7); }
  /// SL
  void DrawInternalHistogramClone8() { DrawInternalHistogramClone(8); }
  /// SL
  void DrawInternalHistogramClone9() { DrawInternalHistogramClone(9); }
  
  /// Whether or not the part of the detector represented by this painter should be included in readout.
  virtual Bool_t IsIncluded() const = 0;
  
  /// Whether or not the part of the detector represented by this painter should be excluded from readout.
  Bool_t IsExcluded() const { return ! IsIncluded(); }

  virtual void Include();
  
  virtual void Exclude();

protected:
    
  virtual TCollection* Children() const;

  void CreateGroups();

  AliMUONVTrackerData* InteractiveReadOutConfig() const;
  
  mutable TH1* fHistogram; //!< histogram
  
  TMap* fPainterGroups; ///< map of groups
  AliMUONPainterGroup* fResponderGroup; ///< the responder group

private:
  
  void FlatList(TList& list);

  AliMUONPainterGroup* CreateGroup(const char* type, Int_t depth);
  
  void GetBoundingBox(Double_t& x1, Double_t& y1, 
                      Double_t& x2, Double_t& y2) const;
  
  AliMUONVPainter* GetPainter(Int_t px, Int_t py, Double_t& x, Double_t& y) const;
  
  void WriteIROC(Double_t value);
  
  void GetIROCManuList(TObjArray& manuList);

  TString fName; ///< our (short) name
  TString fPathName; ///< our long name
  TString fType; ///< our type (DE, Chamber, MANU, etc...)
  AliMUONVPainter* fMother;  ///< our mother
  AliMUONPainterGroup* fGroup; ///< our group
  AliMUONContour* fContour;  ///< our contour
  TObjArray* fChildren; ///< our children
  AliMUONPainterGroup* fPlotterGroup; ///< the plotter group
  Double_t fBorderFactor; ///< border factor for visu 
  TVirtualPad* fPad; ///< the pad we're drawn in
  AliMUONAttPainter fAttributes; ///< our attributes (e.g. view type)
  Int_t fLineColor; ///< our outline color
  Int_t fLineWidth; ///< our outline width
  Int_t fID[2]; ///< our ids
  Bool_t fIsValid; ///< whether we were properly initialized
  
  ClassDef(AliMUONVPainter,3) // Base class of a graphical object for MUON spectrometer
};

#endif
