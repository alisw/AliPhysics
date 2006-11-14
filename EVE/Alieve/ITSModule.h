#ifndef ALIEVE_ITSModule_H
#define ALIEVE_ITSModule_H

#include <Reve/QuadSet.h>

#include <Alieve/ITSDigitsInfo.h>

namespace Alieve {

class ITSModule : public Reve::QuadSet
{
  ITSModule(const ITSModule&);            // Not implemented
  ITSModule& operator=(const ITSModule&); // Not implemented

private:
  void LoadQuads();

protected:
  virtual void InitModule();
  virtual void SetTrans();

  ITSDigitsInfo* fInfo; 

  Int_t       fID;    // Module   id
  Int_t       fDetID; // Detector id (0~SPD, 1~SDD, 2~SSD)

  Int_t       fLayer;
  Int_t       fLadder;
  Int_t       fDet;
  
  Float_t     fDx;
  Float_t     fDz;
  Float_t     fDy;

public:
  ITSModule(const Text_t* n="ITSModule", const Text_t* t=0);
  ITSModule(Int_t id, ITSDigitsInfo* info);
  virtual ~ITSModule();

  ITSDigitsInfo* GetDigitsInfo() const { return fInfo; }
  void SetDigitsInfo(ITSDigitsInfo* info);

  Int_t GetID() const { return fID; }
  void  SetID(Int_t id);

  virtual void Print(Option_t* opt="") const;

  static Bool_t    fgStaticInitDone;
  static void      InitStatics(ITSDigitsInfo* info);

  static Reve::FrameBox* fgSPDFrameBox;
  static Reve::FrameBox* fgSDDFrameBox;
  static Reve::FrameBox* fgSSDFrameBox;

  static Reve::RGBAPalette* fgSPDPalette;
  static Reve::RGBAPalette* fgSDDPalette;
  static Reve::RGBAPalette* fgSSDPalette;

  ClassDef(ITSModule, 1);
};

}

#endif
