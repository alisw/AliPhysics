/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/////////////////////////////////////////////////////////////////////////
// ALICE SHUTTER FRAME CLASS                                           //
// Author: Mayeul   ROUSSELET                                          //
// e-mail: Mayeul.Rousselet@cern.ch                                    //
// Last update:26/08/2003                                              //
/////////////////////////////////////////////////////////////////////////

#include <TGFrame.h>
#include <TGLayout.h>
#include <TGShutter.h>

#include "AliDetectorFrame.h"
#include "AliDisplay2.h"
#include "AliShutterItem.h"

#include "AliShutterFrame.h"


ClassImp(AliShutterFrame)

//_____________________________________________________________
AliShutterFrame::AliShutterFrame(TGCompositeFrame *p, UInt_t /*w*/, UInt_t h)
{
  // Constructor
  fShutter = new TGShutter(p,kSunkenFrame);
  fLayout = new TGLayoutHints(kLHintsExpandY | kLHintsTop | kLHintsLeft);
  fMainFrame = (TGCompositeFrame *) fShutter;
  
  //Event Shutter
  AliShutterItem *item = new AliShutterItem(fShutter,"Event",kIdsEVENT);
  
  item->AddPictureButton("next.xpm","Show next event",kIdbNextEVENT);
  item->AddPictureButton("prev.xpm","Show previous event",kIdbPrevEVENT);
  
  //View Shutter
  item = new AliShutterItem(fShutter,"View",kIdsVIEW);
  item->AddPictureButton("top.xpm","Top view",kIdbTOPVIEW);
  item->AddPictureButton("side.xpm","Side view",kIdbSIDEVIEW);
  item->AddPictureButton("front.xpm","Front view",kIdbFRONTVIEW);
  item->AddPictureButton("four.xpm","Four views",kIdbALLVIEW);
  
  //Detector Shutter
  item = new AliShutterItem(fShutter,"Detectors",kIdsDETECTORS);
  TGCompositeFrame *frame = item->GetShutterItemFrame();
  fDetectorFrameLayout = new TGLayoutHints( kLHintsTop | kLHintsLeft| kLHintsExpandX | kLHintsCenterX,5,5,5,5);
  fDetectorFrame = new AliDetectorFrame(frame,200,200,item->GetShutterItem()->GetDefaultFrameBackground());
  frame->AddFrame(fDetectorFrame->GetDetectorFrame(),fDetectorFrameLayout);
  
  //Options Shutter
  item = new AliShutterItem(fShutter,"Options",kIdsOPTIONS);
  item->AddCheckButton("Display Hits",kIdbCheckHITS);
  item->AddCheckButton("Display Clusters",kIdbCheckCLUSTERS);
  item->AddCheckButton("Display HLT Clusters",kIdbCheckHLT);
  //	item->AddCheckButton("Display Tracks",kIdbCheckTRACKS);
  
  fMainFrame->Resize(150,h);
}

//_____________________________________________________________
AliShutterFrame::~AliShutterFrame(void)
{
  // Destructor
  delete fLayout;
  delete fShutter;
  delete fMainFrame;
  delete fDetectorFrame;
  delete fDetectorFrameLayout;
}

