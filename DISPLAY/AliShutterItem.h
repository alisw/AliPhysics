#ifndef ALISHUTTERITEM_H
#define ALISHUTTERITEM_H

/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////////
// ALICE SHUTTER ITEM CLASS                                            //
// Author: Mayeul   ROUSSELET                                          //
// e-mail: Mayeul.Rousselet@cern.ch                                    //
// Last update:26/08/2003                                              //
/////////////////////////////////////////////////////////////////////////

#include <Rtypes.h>
#include <RQ_OBJECT.h>

class TGShutter;
class TGShutterItem;
class TGCompositeFrame;
class TGButton;

class AliShutterItem{
  //This class implements the shutter item, ie the base element of a shutter and provides functions to add button... in the shutter
public:
	
 AliShutterItem(TGShutter *s, const char *text,UInt_t id);
 virtual ~AliShutterItem();

 //Getters
 TGShutterItem*		GetShutterItem()const {return fShutterItem;};
 TGCompositeFrame*	GetShutterItemFrame()const {return fMainFrame;};

 //Fill functions
 void	      		AddTextButton(const char *text, const char *tiptext,  UInt_t idb);
 void	       		AddPictureButton(const char *file, const char *tiptext,UInt_t idb);
 void	       		AddCheckButton(const char *txt,Int_t idb);

 //Slot
 void	       		DoButton(Int_t pos=0) const;

private:

 TGCompositeFrame	*fMainFrame; // Main frame
 TGShutterItem		*fShutterItem; // Shutter item
 TGButton		*fButton; // Button

 RQ_OBJECT("AliShutterItem")

 ClassDef(AliShutterItem,0);
};

#endif
