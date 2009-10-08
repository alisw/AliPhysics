//-*- Mode: C++ -*-
// $Id: AliHLTEMCALGetEventButton.h 29824 2008-11-10 13:43:55Z richterm $

#ifndef ALIHLTEMCALGETEVENTBUTTON_H
#define ALIHLTEMCALGETEVENTBUTTON_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

#include <TGButton.h>
#include <TGFrame.h>

class AliHLTEMCALOnlineDisplay;

class AliHLTEMCALGetEventButton : public TGTextButton
{
 public:
  AliHLTEMCALGetEventButton();
  AliHLTEMCALGetEventButton(TGGroupFrame *gfPtr, char *name, char opt ='e');
  AliHLTEMCALGetEventButton(TGCompositeFrame *gfPtr, char *name, char opt='e');
  //  AliHLTEMCALGetEventButton(AliHLTEMCALOnlineDisplay *gfPtr, char *name); 
  virtual Bool_t HandleButton(Event_t* event);
 private:
  AliHLTEMCALOnlineDisplay* onlineDisplayPtr;
  char fOption;
};

#endif
