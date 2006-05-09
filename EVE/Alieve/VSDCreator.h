// $Header$

#ifndef ALIEVE_VSDCreator_H
#define ALIEVE_VSDCreator_H

#include <Reve/VSD.h>

class AliTPCParam;
class AliRunLoader;

#include <map>

namespace Alieve {

class VSDCreator : public Reve::VSD
{
public:
  enum KineType_e { KT_Standard, KT_ProtonProton };

protected:
  void          MakeItsDigitsInfo();
  Reve::GenInfo* GetGeninfo(Int_t label);
  AliTPCParam*   GetTpcParam(const Reve::Exc_t& eh);

  KineType_e    mKineType;  // X{GS} 7 PhonyEnum()
  TString       mDataDir;   // X{G}
  Int_t         mEvent;     // X{G}

  Float_t       mTRDHitRes;  // X{gs} 
  Float_t       mTPCHitRes;  // X{gs} 

  Int_t         mDebugLevel;

  std::map<Int_t, Reve::GenInfo*> mGenInfoMap; //!

  void init();

public:
  VSDCreator();
  VSDCreator(const Text_t* name, const Text_t* title="");
  virtual ~VSDCreator() {}

  void CreateVSD(const Text_t* data_dir, Int_t event,
                 const Text_t* vsd_file);  // X{Ed} 

  void CreateTrees();

  // --------------------------------------------------------------
  // Conversion functions.

  void ConvertKinematics();
  void ConvertHits();
  void ConvertClusters();
  void ConvertTPCClusters();
  void ConvertITSClusters();
  void ConvertV0();
  void ConvertKinks();
  void ConvertRecTracks();
  void ConvertGenInfo();

  // --------------------------------------------------------------
  // Get/Set crap
  Int_t GetDebugLevel() const   { return mDebugLevel; }
  void  SetDebugLevel(Int_t dl) { mDebugLevel = dl; }

  // --------------------------------------------------------------
  // Globals.

  AliRunLoader* pRunLoader;

  ClassDef(VSDCreator, 1);
}; // endclass VSDCreator

}

#endif
