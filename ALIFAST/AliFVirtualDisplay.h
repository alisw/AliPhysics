#ifndef AliFVirtualDisplay_H
#define AliFVirtualDisplay_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliFVirtualDisplay                                                   //
//                                                                      //
// Virtual base class for AliFast event display                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TObject
#include <TObject.h>
#endif

class AliFTrigger;

class AliFVirtualDisplay : public TObject {

public:
                     AliFVirtualDisplay();
   virtual          ~AliFVirtualDisplay();
   virtual void      Clear(Option_t *option="") = 0;
   virtual void      DisplayButtons() = 0;
   virtual Int_t     DistancetoPrimitive(Int_t px, Int_t py) = 0;
   virtual void      Draw(Option_t *option="") = 0;
   virtual void      DrawAllViews()  = 0;
   virtual Bool_t    DrawParticles() = 0;
   virtual void      DrawTitle(Option_t *option="") = 0;
   virtual void      DrawView(Float_t theta, Float_t phi) = 0;
   virtual void      DrawViewGL() = 0;
   virtual void      DrawViewX3D() = 0;
   virtual void      ExecuteEvent(Int_t event, Int_t px, Int_t py) = 0;
   virtual void      Paint(Option_t *option="") = 0;
   virtual void      PaintFruit(TObject *obj, Float_t eta, Float_t phi, Float_t pt, Int_t type, Option_t *option="") = 0;
   virtual void      PaintParticles(Option_t *option="") = 0;
   virtual Float_t   PTcut() = 0;
   virtual Float_t   PTcutEGMUNU() = 0;
   virtual void      SetView(Float_t theta, Float_t phi) = 0;
   virtual void      ShowNextEvent(Int_t delta=1) = 0;
   virtual void      SizeFruit() const;
   virtual void      SizeParticles() const;

   ClassDef(AliFVirtualDisplay, 0)   //Virtual base class for AliFast event display
};

#endif
