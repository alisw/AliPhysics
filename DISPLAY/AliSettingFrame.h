#ifndef ALISETTINGFRAME_H
#define ALISETTINGFRAME_H
/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////////
// ALICE SETTING FRAME CLASS                                           //
// Author: Mayeul   ROUSSELET                                          //
// e-mail: Mayeul.Rousselet@cern.ch                                    //
// Last update:26/08/2003                                              //
/////////////////////////////////////////////////////////////////////////

#include <Rtypes.h>
#include <RQ_OBJECT.h>

class TGWindow;
class TGCompositeFrame;
class TGLayoutHints;
class TGNumberEntryField;
class TGLabel;
class TGCheckButton;

class AliSettingFrame:public TGTransientFrame{
  //This classe implement the setting frame where the different otption can be set

public:

 AliSettingFrame(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h);
 virtual ~AliSettingFrame();

 //Slots
 void					DoSettings(Int_t id=0) const;

private:

 TGCompositeFrame		*fMainFrame; // Main frame
 TGCompositeFrame		*fZoomStepFrame; // Zoom step frame
 TGLayoutHints			*fZoomStepLayout; // Zoom step layout
 TGNumberEntryField		*fZoomStepEntry; // Zoom step entry
 TGLabel       			*fZoomStepLabel; // zoom step label
 TGCompositeFrame		*fSliderStepFrame; // Slider step frame
 TGLayoutHints			*fSliderStepLayout; // Slider step layout
 TGNumberEntryField		*fSliderStepEntry; // Slider step entry
 TGLabel       			*fSliderStepLabel; // Slider step label
 TGCompositeFrame		*fSliderUpdateFrame; // Slider update frame
 TGLayoutHints			*fSliderUpdateLayout;// Slider update layout
 TGCheckButton                  *fSliderUpdateButton; // Slider update button
 Bool_t                         fIsLoading;//Used when retrieving the state of the check button

 RQ_OBJECT("AliSettingFrame")

 ClassDef(AliSettingFrame,0);
};

#endif
