#ifndef ALIGUIGEOMDIALOG_H
#define ALIGUIGEOMDIALOG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TGFrame.h"

class AliGUISliders;
class TGButton;
class TGComboBox;
class TGLabel;
class TGTab;
class TGTextBuffer;
class TGTextEntry;
class TGDoubleHSlider;
class AliDrawVolume;

class AliGuiGeomDialog : public TGTransientFrame {
public:
   AliGuiGeomDialog(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h,
               UInt_t options = kMainFrame | kVerticalFrame);
   virtual ~AliGuiGeomDialog();
// Destroy this window
   virtual void CloseWindow();
// Process messages from this window    
   virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);
// Update widgets   
   virtual void Update();
private:
    AliGUISliders       *fF1;                                          // Slider for Draw Control
    TGCompositeFrame    *fFrame1, *fF2, *fF3, *fF4;                            // Outer frames
    TGButton            *fOkButton, *fCancelButton;                            // Buttons
    TGButton            *fChk1, *fChk2, *fChk3;                                // Buttons
    TGComboBox          *fCombo, *fCombo2;                                     // Combo Boxes
    TGLabel             *fLabel1, *fLabel2;                                    // Labels
    TGTab               *fTab;                                                 // Tab Entries
    TGLayoutHints       *fL1, *fL2, *fL3, *fL4, *fBly, *fBfly1;                // Layout hints
    TGHorizontalFrame   *fHSframe1, *fHSframe2, *fHSframe3;                    // Horizontal frames
    TGTextBuffer        *fTbh11, *fTbh12, *fTbh21, *fTbh22, *fTbh31, *fTbh32;  // Text buffers
    TGTextEntry         *fTeh11, *fTeh12, *fTeh21, *fTeh22, *fTeh31, *fTeh32;  // Text Entries
    TGDoubleHSlider     *fDslider1, *fDslider2, *fDslider3;                    // Sliders for clip box
    TGLabel             *fSLabel1,  *fSLabel2,  *fSLabel3;                     // Labels

private:
  AliGuiGeomDialog(const AliGuiGeomDialog& gd):
    TGTransientFrame((const TGTransientFrame&)gd) {}
  virtual AliGuiGeomDialog & operator=(const AliGuiGeomDialog &) 
  {return *this;}
};

R__EXTERN AliDrawVolume  *gCurrentVolume;

#endif
