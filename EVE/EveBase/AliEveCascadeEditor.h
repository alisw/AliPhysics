// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
 
//-----------------------------------------------------------------------------
// This code defines the Editor coming with the visualisation of cascades,
// within AliEVE
//
// Origin :   Boris Hippolyte, IPHC (hippolyt@in2p3.fr)
// Modified : Antonin Maire, April 2009, IPHC (antonin.maire@cern.ch)
//-----------------------------------------------------------------------------
 

#ifndef ALIEVECASCADEEDITOR_H
#define ALIEVECASCADEEDITOR_H

#include "TGedFrame.h"

class TGButton;
class TGLabel;
//class TGCheckButton;
//class TGNumberEntry;
//class TGColorSelect;

class AliEveCascade;

//______________________________________________________________________________
// Short description of AliEveCascadeEditor
//

class AliEveCascadeEditor : public TGedFrame
{
public:
  AliEveCascadeEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
                 UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveCascadeEditor() {}

  virtual void SetModel(TObject* obj);

  // Declare callback/slot methods
  // void DoXYZZ();
  void DisplayDetailed();
  void DisplayMassHyp();

protected:
  AliEveCascade  *fM; //! Model object.

  TGLabel   *fInfoLabelRadius; //! label about transverse radius decay for the cascade
  TGLabel   *fInfoLabelDCA;    //! label about the DCA between Xi daughters
  TGLabel   *fInfoLabelCharge; //! label about the charge of the cascade
  TGLabel   *fInfoLabelPhi;    //! label about phi
  TGLabel   *fInfoLabelTheta;  //! label about theta
  TGLabel   *fInfoLabelPtot;   //! label about total momentum, Ptot
  TGLabel   *fInfoLabelPt;     //! label about transverse momentum, Pt
  TGLabel   *fInfoLabelEta;    //! label about pseudo-rapidity
  
  TGButton  *fXButtonDetailedView;  //! button to get the detailed view 
  TGButton  *fXButtonMassHyp;       //! button to printf the calculation of eff inv mass, under Xi and Omega hypotheses

private:
  AliEveCascadeEditor(const AliEveCascadeEditor&);            // Not implemented
  AliEveCascadeEditor& operator=(const AliEveCascadeEditor&); // Not implemented

  ClassDef(AliEveCascadeEditor, 0); // GUI editor for AliEveCascade.
};

#endif
