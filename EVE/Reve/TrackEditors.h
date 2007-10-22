// $Header$

#ifndef REVE_TrackEditors_H
#define REVE_TrackEditors_H

#include <TGedFrame.h>

class TGButton;
class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TGComboBox;

namespace Reve {

class RGValuator;
class RGDoubleValuator;
class TrackRnrStyleSubEditor;
class Track;
class TrackList;

/**************************************************************************/
// TrackEditor
/**************************************************************************/
class TrackEditor : public TGedFrame
{
private:
  TrackEditor(const TrackEditor&);            // Not implemented
  TrackEditor& operator=(const TrackEditor&); // Not implemented

protected: 
  Track                          *fM; 
  TGTextButton                   *fRSEditor;
public:
  TrackEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
		  UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  ~TrackEditor(){}

  virtual void SetModel(TObject* obj);
  void DoEditRnrStyle();
 
  ClassDef(TrackEditor, 1); // Editor for Track
}; // endclass TrackEditor

/**************************************************************************/
// TrackListEditor
/**************************************************************************/

class TrackListEditor : public TGedFrame
{
private:
  TrackListEditor(const TrackListEditor&);            // Not implemented
  TrackListEditor& operator=(const TrackListEditor&); // Not implemented

  void CreateRefTab();
protected: 
  TGVerticalFrame                 *fRefs;

  TrackList                       *fTC; // fModel dynamic-casted to TrackListEditor

  TGCheckButton                   *fRnrLine;
  TGCheckButton                   *fRnrPoints;

  RGDoubleValuator                *fPtRange;
  RGDoubleValuator                *fPRange;

  TrackRnrStyleSubEditor          *fRSSubEditor;

public:
  TrackListEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
		  UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  ~TrackListEditor();

  void CreateRefsTab();
  virtual void SetModel(TObject* obj);

  void DoRnrLine();
  void DoRnrPoints();
 
  void DoPtRange();
  void DoPRange();

  ClassDef(TrackListEditor, 1); // Editor for TrackList
}; // endclass TrackListEditor



/**************************************************************************/
// 
/**************************************************************************/

class TrackCounter;

class TrackCounterEditor : public TGedFrame
{
private:
  TrackCounterEditor(const TrackCounterEditor&);            // Not implemented
  TrackCounterEditor& operator=(const TrackCounterEditor&); // Not implemented

protected:
  TrackCounter* fM; // fModel dynamic-casted to TrackCounter

  // Declare widgets
  TGComboBox*    fClickAction;
  TGLabel*       fInfoLabel;
  TGNumberEntry* fEventId;

public:
  TrackCounterEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
		     UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~TrackCounterEditor();

  virtual void SetModel(TObject* obj);

  void DoOrtoXY();
  void DoOrtoZY();
  void DoPersp();

  void DoPrev();
  void DoNext();
  void DoSetEvent();

  void DoPrintReport();
  void DoFileReport();
  void DoShowHistos();

  void DoClickAction(Int_t);

  ClassDef(TrackCounterEditor, 1); // Editor for TrackCounter
}; // endclass TrackCounterEditor

}

#endif
