// $Header$

#ifndef REVE_CLASS_H
#define REVE_CLASS_H

#include <TGedFrame.h>

class TGButton;
class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

namespace Alieve {

class STEM;

class XXCLASS : public TGVerticalFrame
{
private:
   XXCLASS(const XXCLASS&);            // Not implemented
   XXCLASS& operator=(const XXCLASS&); // Not implemented

protected:
   STEM       *fM;

public:
   XXCLASS(const TGWindow* p);
   virtual ~XXCLASS() {}

   void SetModel(STEM* m);

   void Changed(); //*SIGNAL*

   // void DoABCD();

   ClassDef(XXCLASS, 0) // Sub-editor for STEM
};


class CLASS : public TGedFrame
{
private:
   CLASS(const CLASS&);            // Not implemented
   CLASS& operator=(const CLASS&); // Not implemented

protected:
   STEM      *fM;  // fModel dynamic-casted to STEM
   XXCLASS   *fSE;

public:
   CLASS(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
   virtual ~CLASS() {}

   virtual void SetModel(TObject* obj);

   void DoXYZZ();

   ClassDef(CLASS, 0) // Editor for STEM
};

}

#endif
