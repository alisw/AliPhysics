// $Header$

#ifndef ALIEVE_ITSModuleStepper_H
#define ALIEVE_ITSModuleStepper_H

#include <TNamed.h>
#include <TAtt3D.h>
#include <TAttBBox.h>

#include <Reve/RenderElement.h>
#include <Reve/GridStepper.h>
#include <Reve/ZTrans.h>

#include <vector>

namespace Alieve {

class ITSDigitsInfo;

class ITSModuleStepper : public Reve::RenderElement, 
                         public TNamed,
		         public TAtt3D,
		         public TAttBBox
{
  friend class ITSModuleStepperGL;

public:
  typedef std::vector<UInt_t>           vpInt_t;
  typedef std::vector<UInt_t>::iterator vpInt_i;

  enum PositionType_e { PT_BottomLeft, PT_BottomRight, PT_TopLeft, PT_TopRight };

private:
  ITSModuleStepper(const ITSModuleStepper&);            // Not implemented
  ITSModuleStepper& operator=(const ITSModuleStepper&); // Not implemented

protected:
  ITSDigitsInfo*          fDigitsInfo;
  Reve::GridStepper*      fStepper;  
 
  Float_t                 fExpand;

  vpInt_t                 fIDs;
  UInt_t                  fPosition;

  Reve::ZTrans            fHMTrans; 

  Bool_t                  fRnrFrame;    
  PositionType_e          fWCorner; 
  Color_t                 fWColor;   
  Float_t                 fWWidth;
  Float_t                 fWHeight;

  void                    Apply();
  Int_t                   Nxy(){ return fStepper->Nx*fStepper->Ny; }

  void                    SetFirst(Int_t first);

public:
  ITSModuleStepper(ITSDigitsInfo* di);
  virtual ~ITSModuleStepper();

  void   Start();
  void   Next();
  void   Previous();
  void   End();

  void   SetStepper(Int_t nx, Int_t ny, Float_t dx = -1, Float_t dy = -1);
  Reve::GridStepper*  GetStepper(){ return fStepper; }
  
  void   AddToList( Int_t modID ){ fIDs.push_back(modID);}
  void   ResetList(){ fIDs.clear();}

  void   DisplayDet(Int_t det, Int_t layer = -1);
  void   DisplayTheta(Float_t min, Float_t max);

  Int_t  GetCurrentPage();
  Int_t  GetPages();

  virtual Bool_t CanEditMainColor() { return kTRUE; }

  Bool_t  GetRnrFrame(){ return fRnrFrame; }
  void    SetRnrFrame(Bool_t rnr){ fRnrFrame = rnr; }
  Color_t GetWColor(){ return fWColor; };
  void    SetWColor(Color_t c){ fWColor=c; }

  virtual Reve::ZTrans* PtrMainHMTrans()     { return &fHMTrans; }

  virtual void ComputeBBox();

  virtual void Paint(Option_t* option = "");

  ClassDef(ITSModuleStepper, 0);
};

} // Alieve namespace

#endif
