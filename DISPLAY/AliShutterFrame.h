#ifndef ALISHUTTERFRAME_H
#define ALISHUTTERFRAME_H
/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////////
// ALICE SHUTTER FRAME CLASS                                           //
// Author: Mayeul   ROUSSELET                                          //
// e-mail: Mayeul.Rousselet@cern.ch                                    //
// Last update:26/08/2003                                              //
/////////////////////////////////////////////////////////////////////////

#include <Rtypes.h>
#include <RQ_OBJECT.h>

class TGCompositeFrame;
class TGLayoutHints;
class TGShutter;
class AliDetectorFrame;

class AliShutterFrame{
  //This class implements the shutter frame
public:

 AliShutterFrame(TGCompositeFrame *p, UInt_t w, UInt_t h);
 virtual ~AliShutterFrame();

 TGCompositeFrame*	GetShutterFrame(){return fMainFrame;};

private:

 TGCompositeFrame	*fMainFrame; // Main frame
 TGLayoutHints		*fLayout; // Layout
 TGShutter	       	*fShutter; // Shutter
 AliDetectorFrame	*fDetectorFrame; // Detector frame
 TGLayoutHints		*fDetectorFrameLayout; // Detector frame layout

 RQ_OBJECT("AliShutterFrame")

 ClassDef(AliShutterFrame,0);
};

#endif
