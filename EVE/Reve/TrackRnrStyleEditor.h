// $Header$

#ifndef REVE_TrackRnrStyleEditor_H
#define REVE_TrackRnrStyleEditor_H

#include <TGedFrame.h>

class TGButton;
class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TGComboBox;
class TGLineWidthComboBox;
class TGLineStyleComboBox;

class TAttMarkerEditor;

namespace Reve {

class TrackRnrStyle;

class RGValuator;
class RGDoubleValuator;
class TrackRnrStyleSubEditor;

class TrackRnrStyleSubEditor : public TGVerticalFrame
{
  friend class TrackRnrStyleEditor;
  friend class TrackListEditor;

private:
  TrackRnrStyleSubEditor(const TrackRnrStyleSubEditor&);            // Not implemented
  TrackRnrStyleSubEditor& operator=(const TrackRnrStyleSubEditor&); // Not implemented

protected:
  TrackRnrStyle       *fM;  
 
  Reve::RGValuator    *fMaxR;
  Reve::RGValuator    *fMaxZ;
  Reve::RGValuator    *fMaxOrbits;
  Reve::RGValuator    *fMinAng;
  Reve::RGValuator    *fDelta;

  TGCheckButton       *fRnrFV;

  TGCompositeFrame    *fPMFrame;
  TGButton            *fFitDaughters;
  TGButton            *fFitReferences;
  TGButton            *fFitDecay;
  TGButton            *fRnrDaughters;
  TGButton            *fRnrReferences;
  TGButton            *fRnrDecay;

  TGCompositeFrame    *fRefsCont;

  TAttMarkerEditor    *fPMAtt;
  TAttMarkerEditor    *fFVAtt;

public:
  TrackRnrStyleSubEditor(const TGWindow* p);
  virtual ~TrackRnrStyleSubEditor() {}

  void SetModel(TrackRnrStyle* m);

  void Changed(); //*SIGNAL*

  void DoMaxR();
  void DoMaxZ();
  void DoMaxOrbits();
  void DoMinAng();
  void DoDelta();

  void DoFitPM();
  void DoRnrPM();
  
  void DoRnrFV();

  void CreateRefsContainer(TGVerticalFrame* p);

  ClassDef(TrackRnrStyleSubEditor, 0) // Sub-editor for TrackRnrStyle.
};

/**************************************************************************/
// TrackRnrStyleEditor
/**************************************************************************/

class TrackRnrStyleEditor : public TGedFrame
{
private:
  TrackRnrStyleEditor(const TrackRnrStyleEditor&);            // Not implemented
  TrackRnrStyleEditor& operator=(const TrackRnrStyleEditor&); // Not implemented

  void CreateRefTab();
protected: 
  TrackRnrStyle           *fM;           // Model object.
  TrackRnrStyleSubEditor  *fRSSubEditor; // Render-style sub-editor.

public:
  TrackRnrStyleEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
		  UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  ~TrackRnrStyleEditor();

  virtual void SetModel(TObject* obj);

  ClassDef(TrackRnrStyleEditor, 1); // Editor for TrackRnrStyle.
}; // endclass TrackRnrStyleEditor

}

#endif
