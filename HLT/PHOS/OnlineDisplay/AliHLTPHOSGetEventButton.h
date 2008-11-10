//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSGETEVENTBUTTON_H
#define ALIHLTPHOSGETEVENTBUTTON_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

#include <TGButton.h>
#include <TGFrame.h>

class AliHLTPHOSOnlineDisplay;

class AliHLTPHOSGetEventButton : public TGTextButton
{
 public:
  AliHLTPHOSGetEventButton();
  AliHLTPHOSGetEventButton(TGGroupFrame *gfPtr, char *name, char opt ='e');
  AliHLTPHOSGetEventButton(TGCompositeFrame *gfPtr, char *name, char opt='e');
  //  AliHLTPHOSGetEventButton(AliHLTPHOSOnlineDisplay *gfPtr, char *name); 
  virtual Bool_t HandleButton(Event_t* event);
 private:
  AliHLTPHOSOnlineDisplay* onlineDisplayPtr;
  char fOption;
};

#endif
