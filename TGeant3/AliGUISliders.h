#ifndef ALIGUISLIDERS_H
#define ALIGUISLIDERS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TGFrame.h"
#include "TGSlider.h"
#include "TGTextEntry.h"
#include "TGTextBuffer.h"
#include "TGLabel.h"

class AliDrawVolume;

class AliGUISliders : public  TGCompositeFrame {
public:
   AliGUISliders(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h);
   virtual ~AliGUISliders();
   virtual void CloseWindow();
   virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);
   virtual void Update();
private:
//
    TGHorizontalFrame *fHframe[8];       // 8 Horizontal frames for sliders
    TGLayoutHints     *fBly, *fBfly1;    // Lay-out hints
    TGHSlider         *fHslider[8];      // 8 Sliders
    TGTextEntry       *fTeh[8];          // Text entries for slider position
    TGTextBuffer      *fTbh[8];          // Text buffer  
    TGLabel           *fLabel[8];        // Slider labels
    Text_t            fLabelText[8];     // Label text 

  AliGUISliders(const AliGUISliders &gs) :
    TGCompositeFrame((const TGCompositeFrame&) gs) {}
  AliGUISliders & operator=(const AliGUISliders &) {return *this;}
  
      
   //   ClassDef(AliGUISliders,1)  // Window containing sliders 
};

R__EXTERN AliDrawVolume  *gCurrentVolume;

#endif
