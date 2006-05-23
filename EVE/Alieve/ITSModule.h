#ifndef ALIEVE_ITSModule_H
#define ALIEVE_ITSModule_H

#include <Reve/QuadSet.h>
#include <Reve/RenderElement.h>

#include <Alieve/ITSDigitsInfo.h>

namespace Alieve {

class ITSModule : public Reve::QuadSet, public Reve::RenderElement
{
private:
  void Init();
  void LoadQuads();

protected:
  virtual void InitModule();
  virtual void SetTrans();

  ITSDigitsInfo* fInfo; 

  Int_t       fID;
  Int_t       fDetID;

  Int_t       fLayer;
  Int_t       fLadder;
  Int_t       fDet;
  
  Float_t     fDx;
  Float_t     fDz;
  Float_t     fDy;

  Color_t     fFrameColor;

public:
  ITSModule(const Text_t* n="ITSModule", const Text_t* t=0, Color_t col=2) :
    QuadSet(n, t), Reve::RenderElement(fFrameColor), fFrameColor(col)
  { Init(); }
  ITSModule(Int_t id, ITSDigitsInfo* info, Color_t col=2);
  virtual ~ITSModule();

  virtual Bool_t CanEditMainColor()  { return true; }
  virtual void SetMainColor(Color_t col);

  virtual void SetID(Int_t id);
  virtual void Print(Option_t* opt="") const;

  static Short_t   fgSDDThreshold;
  static Short_t   fgSDDMaxVal;

  static Short_t   fgSSDThreshold;
  static Short_t   fgSSDMaxVal;

  ClassDef(ITSModule, 1);
}; 
}
#endif
