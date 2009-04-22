// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALIEVECASCADELISTEDITOR_H
#define ALIEVECASCADELISTEDITOR_H

//------------------------------------------------------------------------------
// This code defines the List Editor coming with the visualisation of cascades,
// within AliEVE
//
// Origin :   Boris Hippolyte, IPHC (hippolyt@in2p3.fr)
// Modified : Antonin Maire, April 2009, IPHC (antonin.maire@cern.ch)
//------------------------------------------------------------------------------


//class TGButton;
//class TGCheckButton;
//class TGNumberEntry;
//class TGColorSelect;
class TEveGDoubleValuator;
class TGComboBox;

class AliEveCascadeList;

#include "TGedFrame.h"

//______________________________________________________________________________
// Short description of AliEveCascadeListEditor
//

class AliEveCascadeListEditor : public TGedFrame
{
public:
  AliEveCascadeListEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
                     UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveCascadeListEditor() {}

  virtual void SetModel(TObject* obj);

  // Declare callback/slot methods
  void DoMinMaxRCut();
  void DoMinMaxDaughterDCA();
  void DoMinMaxPt();
  
  void DoSelectInvMassHyp(Int_t rInvMassHyp);
  void DoMinMaxInvariantMass();

protected:
  AliEveCascadeList    *fM; //! Model object.

  // Declare widgets
  // TGSomeWidget*   fXYZZ;
  TGComboBox*          fCascadeSpecies;		//! Box meant to choose the mass hyp. to be applied : Xi or Omega ?
  
  TEveGDoubleValuator* fMinMaxRCut;		//! Transverse radius range targeted by the user
  TEveGDoubleValuator* fMinMaxDaughterDCA;	//! DCA (between Xi daughters) range targeted by the user
  TEveGDoubleValuator* fMinMaxPt;		//! Pt range targeted by the user
  TEveGDoubleValuator* fMinMaxInvariantMass;	//! Inv Mass range targeted by the user

private:
  AliEveCascadeListEditor(const AliEveCascadeListEditor&);            // Not implemented
  AliEveCascadeListEditor& operator=(const AliEveCascadeListEditor&); // Not implemented

  ClassDef(AliEveCascadeListEditor, 1); // GUI editor for AliEveCascadeList.
};

#endif
