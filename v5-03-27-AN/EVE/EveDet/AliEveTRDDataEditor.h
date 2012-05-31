#ifndef AliEveTRDDataEditor_H
#define AliEveTRDDataEditor_H

#include <TGedFrame.h>

class AliEveTRDHits;
class AliEveTRDHitsEditor : public TGedFrame
{
public:
  AliEveTRDHitsEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
		      UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveTRDHitsEditor() {}

  virtual void SetModel(TObject* obj);

protected:
  AliEveTRDHits* fM; // Model object.

private:
  AliEveTRDHitsEditor(const AliEveTRDHitsEditor&);            // Not implemented
  AliEveTRDHitsEditor& operator=(const AliEveTRDHitsEditor&); // Not implemented

  ClassDef(AliEveTRDHitsEditor, 0); // Editor for AliEveTRDHits.
};


// class AliEveTRDDigits;
// class AliEveTRDDigitsEditor : public TGedFrame
// {
// public:
//   AliEveTRDDigitsEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
// 			UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
//   virtual ~AliEveTRDDigitsEditor() {}
// 
//   virtual void SetModel(TObject* obj);
// 
// protected:
//   AliEveTRDDigits* fM; // Model object.
// 
// private:
//   AliEveTRDDigitsEditor(const AliEveTRDDigitsEditor&);            // Not implemented
//   AliEveTRDDigitsEditor& operator=(const AliEveTRDDigitsEditor&); // Not implemented
// 
//   ClassDef(AliEveTRDDigitsEditor, 0); // Editor for AliEveTRDDigits
// };
#endif

