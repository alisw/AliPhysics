#ifndef ALIEVE_ITSModule_H
#define ALIEVE_ITSModule_H

#include <Reve/QuadSet.h>

#include <Alieve/ITSDigitsInfo.h>

namespace Alieve {

class ITSModule : public Reve::QuadSet
{
  ITSModule(const ITSModule&);            // Not implemented
  ITSModule& operator=(const ITSModule&); // Not implemented

protected:
  void LoadQuads();
  void SetTrans();

  ITSDigitsInfo* fInfo; 

  Int_t       fID;    // Module   id
  Int_t       fDetID; // Detector id (0~SPD, 1~SDD, 2~SSD)

  Int_t       fLayer;
  Int_t       fLadder;
  Int_t       fDet;
  
  Float_t     fDx;
  Float_t     fDz;
  Float_t     fDy;

  static Bool_t fgStaticInitDone;

public:
  ITSModule(const Text_t* n="ITSModule", const Text_t* t=0);
  ITSModule(Int_t gid, ITSDigitsInfo* info);
  virtual ~ITSModule();

  static void InitStatics(ITSDigitsInfo* info);

  ITSDigitsInfo* GetDigitsInfo() const { return fInfo; }
  void SetDigitsInfo(ITSDigitsInfo* info);

  Int_t GetID() const { return fID; }
  void  SetID(Int_t gid);

  virtual void Print(Option_t* opt="") const;

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
