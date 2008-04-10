// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveITSModuleStepper_H
#define AliEveITSModuleStepper_H

#include <TEveElement.h>
#include <TGLOverlay.h>
#include <TEveGridStepper.h>
#include <TGLFontManager.h>

#include <vector>


class TEveRGBAPalette;
class AliEveITSDigitsInfo;
class AliEveDigitScaleInfo;

class TGLAxis;


class AliEveITSModuleStepper : public TEveElementList,
			       public TGLOverlayElement
{
public:
  typedef std::vector<UInt_t>           vpInt_t;
  typedef std::vector<UInt_t>::iterator vpInt_i;

  AliEveITSModuleStepper(AliEveITSDigitsInfo* di);
  virtual ~AliEveITSModuleStepper();

  TEveGridStepper* GetStepper()                   { return fStepper; }
  void             SetStepper(TEveGridStepper* s) { fStepper = s; Apply(); }

  // overlay
  virtual  Bool_t MouseEnter(TGLOvlSelectRecord& selRec);
  virtual  Bool_t Handle(TGLRnrCtx& rnrCtx, TGLOvlSelectRecord& selRec, Event_t* event);
  virtual  void   MouseLeave();

 // menu callbacks
  void     DisplayDet(Int_t det, Int_t layer = -1);
  Int_t    GetCurrentPage() const;
  Int_t    GetPages();
  void     Start();
  void     Next();
  void     Previous();
  void     End();
  void     Apply();
  void     Capacity();

  virtual void    Render(TGLRnrCtx& rnrCtx);

protected:
  AliEveITSDigitsInfo    *fDigitsInfo; // Source of data and geometry.
  AliEveDigitScaleInfo   *fScaleInfo;  // Parameters for digit-scaling.
  TEveGridStepper        *fStepper;    // Module placement.

  vpInt_t                 fModuleIDs;  // Vector of module IDs to be displayed.
  UInt_t                  fPosition;  // Position of top corner ITS module in vector fIDs.
  Int_t                   fSubDet;    // Sub-det, 0~SPD, 1~SDD, 2~SSD.

  mutable TGLFont         fModuleFont; // Pixmap font for module ids.
  mutable TGLFont         fTextFont;   // Texture font for text tool bar.
  mutable TGLFont         fSymbolFont; // Webdings font for pager and scale actions. 
  TGLAxis*                fAxis;

  Float_t                 fMenuHeight; // Height of a tool bar.
  Int_t                   fTextSize;   // Size of texture for menu font.
  Color_t                 fTextCol;    // Default text color in menu.
  Color_t                 fActiveCol;  // Color of selected menu item.

  Int_t                   fActiveID;   // Id of active menu item.

  // ITS module ID navigation
  Int_t  Nxy()            const { return fStepper->GetNx()*fStepper->GetNy(); }
  void   AddToList(Int_t modID) { fModuleIDs.push_back(modID);}
  void   ResetList()            { fModuleIDs.clear();}
  void   SetFirst(Int_t first);


private:
  AliEveITSModuleStepper(const AliEveITSModuleStepper&);            // Not implemented
  AliEveITSModuleStepper& operator=(const AliEveITSModuleStepper&); // Not implemented

  // GUI
  void   RenderModuleIDs();
  void   RenderText(const char* tex ,Int_t id, const TGLFont &font, Float_t step=-1);
  void   RenderPalette(TEveRGBAPalette* p);
  void   RenderMenu(Int_t currP, Int_t MaxP, Int_t scaleX, Int_t scaleZ);

  ClassDef(AliEveITSModuleStepper, 0); // Display scaled ITS modules in a paged layout, also providing GL-overaly control GUI.
};

#endif
