// $Header$

#ifndef ALIEVE_ITSModuleStepper_H
#define ALIEVE_ITSModuleStepper_H

#include <Reve/RenderElement.h>
#include <Reve/GridStepper.h>

#include <vector>

namespace Alieve {

class ITSDigitsInfo;

class ITSModuleStepper : public Reve::RenderElementList
{
public:
  typedef std::vector<Int_t>           vpInt_t;
  typedef std::vector<Int_t>::iterator vpInt_i;

private:
  ITSModuleStepper(const ITSModuleStepper&);            // Not implemented
  ITSModuleStepper& operator=(const ITSModuleStepper&); // Not implemented

protected:
  ITSDigitsInfo*          fDigitsInfo;
  Reve::GridStepper*      fStepper;  
 
  Float_t                 fExpand;

  vpInt_t                 fIDs;
  vpInt_i                 fPosition;

  void                    Apply();

public:
  ITSModuleStepper(ITSDigitsInfo* di);
  virtual ~ITSModuleStepper();

  void   Start();
  void   Next();
  void   SetStepper(Int_t nx, Int_t ny, Float_t dx = -1, Float_t dy = -1);
  Reve::GridStepper*  GetStepper(){return fStepper;}
  
  void   AddToList( Int_t modID ){ fIDs.push_back(modID);}
  void   ResetList(){ fIDs.clear();}

  ClassDef(ITSModuleStepper, 0);
};

} // Alieve namespace

#endif
