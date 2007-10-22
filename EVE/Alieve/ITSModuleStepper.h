// $Header$

#ifndef ALIEVE_ITSModuleStepper_H
#define ALIEVE_ITSModuleStepper_H

#include <TNamed.h>
#include <TGLOverlay.h>

#include <Reve/RenderElement.h>
#include <Reve/GridStepper.h>

#include <vector>

class TGLText;
class TGLAxis;

namespace Alieve {

class ITSDigitsInfo;
class DigitScaleInfo;

class ITSModuleStepper : public Reve::RenderElementList,
                         public TGLOverlayElement
{
  friend class ITSModuleStepperGL;

public:
  typedef std::vector<UInt_t>           vpInt_t;
  typedef std::vector<UInt_t>::iterator vpInt_i;

private:
  vpInt_t                 fIDs;
  UInt_t                  fPosition;  // position of top corner ITS module in vector fIDs

  ITSModuleStepper(const ITSModuleStepper&);            // Not implemented
  ITSModuleStepper& operator=(const ITSModuleStepper&); // Not implemented

protected:
  ITSDigitsInfo*          fDigitsInfo;
  DigitScaleInfo*         fScaleInfo;
  Int_t                   fSubDet;

  Reve::GridStepper*      fStepper;
  TGLAxis*                fAxis;
  TGLText*                fText;
  Float_t                 fTextSize;
  Float_t                 fPagerGap;
  Bool_t                  fRnrFrame;

  // module configuration
  Float_t                 fExpandCell;
  Color_t                 fModuleFrameCol;

  // palette configuratiom
  Float_t                 fPaletteOffset;
  Float_t                 fPaletteLength;  

  // symbol configuration 
  Int_t                   fWActive; 
  Float_t                 fWWidth;
  Float_t                 fWHeight;
  Float_t                 fWOff; ///offset relative to widget size
  Color_t                 fWCol;  
  Int_t                   fWActiveCol;
  Color_t                 fFontCol;

  // wrappers
  Float_t TextLength(const char* txt);
  void    RenderString(TString tex ,Int_t id = -1);
  void    RenderFrame(Float_t dx, Float_t dy, Int_t id);
  void    RenderSymbol(Float_t dx, Float_t dy, Int_t id);
  void    RenderPalette(Float_t dx, Float_t x, Float_t y);
  void    RenderMenu();
  void    RenderCellIDs();

  // module ID navigation
  Int_t  Nxy(){ return fStepper->Nx*fStepper->Ny; }
  void   AddToList( Int_t modID ){ fIDs.push_back(modID);}
  void   ResetList(){ fIDs.clear();}
  void   SetFirst(Int_t first);

public:
  ITSModuleStepper(ITSDigitsInfo* di);
  virtual ~ITSModuleStepper();

  // external functions
  void     DisplayDet(Int_t det, Int_t layer = -1);
  void     DisplayTheta(Float_t min, Float_t max);

  // overlay functions
  virtual  Bool_t MouseEnter(TGLOvlSelectRecord& selRec);
  virtual  Bool_t Handle(TGLRnrCtx& rnrCtx, TGLOvlSelectRecord& selRec,
                        Event_t* event);
  virtual void   MouseLeave();
  virtual void   Render(TGLRnrCtx& rnrCtx);

  // stepper
  Reve::GridStepper*  GetStepper(){return fStepper;}
  void                SetStepper(Reve::GridStepper* s){ fStepper = s; Apply();}

  Int_t    GetCurrentPage();
  Int_t    GetPages();
  void     Start();
  void     Next();
  void     Previous();
  void     End();
  void     Apply();
  void     Capacity();


  // getters/setters
  Color_t  GetWColor(){ return fWCol; };
  void     SetWColor(Color_t c){ fWCol = c; }
  TGLText* GetFont(){ return fText; }
  void     SetGLText(TGLText* t) {fText = t;}

  ClassDef(ITSModuleStepper, 0);
};

} // Alieve namespace

#endif
