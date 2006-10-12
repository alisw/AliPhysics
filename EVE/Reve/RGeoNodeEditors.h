// $Header$

#ifndef REVE_RGeoNodeEditors_H
#define REVE_RGeoNodeEditors_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

namespace Reve {

class GeoNodeRnrEl;
class GeoTopNodeRnrEl;


class GeoNodeRnrElEditor : public TGedFrame
{
  GeoNodeRnrElEditor(const GeoNodeRnrElEditor&);            // Not implemented
  GeoNodeRnrElEditor& operator=(const GeoNodeRnrElEditor&); // Not implemented

protected:
  GeoNodeRnrEl*   fNodeRE;

  TGCheckButton*  fVizNode;
  TGCheckButton*  fVizNodeDaughters;
  TGCheckButton*  fVizVolume;
  TGCheckButton*  fVizVolumeDaughters;

  TGNumberEntry*  fTransparency;

public:
  GeoNodeRnrElEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
		     UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~GeoNodeRnrElEditor() {}

  virtual void SetModel(TObject* obj);

  void DoVizNode();
  void DoVizNodeDaughters();
  void DoVizVolume();
  void DoVizVolumeDaughters();

  void DoTransparency();

  ClassDef(GeoNodeRnrElEditor, 1);
};

/**************************************************************************/

class GeoTopNodeRnrElEditor : public TGedFrame
{
  GeoTopNodeRnrElEditor(const GeoTopNodeRnrElEditor&);            // Not implemented
  GeoTopNodeRnrElEditor& operator=(const GeoTopNodeRnrElEditor&); // Not implemented

protected:
  GeoTopNodeRnrEl*   fTopNodeRE;

  TGNumberEntry*     fVisOption;
  TGNumberEntry*     fVisLevel;

public:
  GeoTopNodeRnrElEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
			UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~GeoTopNodeRnrElEditor() {}

  virtual void SetModel(TObject* obj);

  void DoVisOption();
  void DoVisLevel();

  ClassDef(GeoTopNodeRnrElEditor, 1);
};

}

#endif
