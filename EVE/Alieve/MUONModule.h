#ifndef ALIEVE_MUONModule_H
#define ALIEVE_MUONModule_H

#include <Reve/QuadSet.h>
#include <Reve/RenderElement.h>

#include <Alieve/MUONDigitsInfo.h>

namespace Alieve {

class MUONModule : public Reve::RenderElement,
                   public Reve::QuadSet
{
  MUONModule(const MUONModule&);            // Not implemented
  MUONModule& operator=(const MUONModule&); // Not implemented

public:

  MUONModule(const Text_t* n="MUONModule", const Text_t* t=0, Color_t col=2);
  MUONModule(Int_t id, Int_t cath, MUONDigitsInfo* info, Bool_t dig, Bool_t clus, Color_t col=2);
  virtual ~MUONModule() {}

  virtual void SetID(Int_t id, Int_t cath);

private:
    
  void LoadQuadsChambers(Int_t chamber1, Int_t chamber2, Int_t delElemId = -1, Int_t cat = -1);
  void LoadQuadsDigits();
  void LoadQuadsClusters();
  void LoadQuadsTracks(Int_t id);

protected:

  virtual void InitModule();
  virtual void SetTrans();
 
  MUONDigitsInfo* fInfo; 
  
  Int_t       fID;
  Int_t       fCath;
 
  Bool_t      fShowDigits;
  Bool_t      fShowClusters;
  Bool_t      fShowTracks;

  Color_t     fFrameCol;

  Int_t       fDetElemId;

  ClassDef(MUONModule,1);

};
}

#endif
