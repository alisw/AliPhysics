// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALIEVE_VSDCreator_H
#define ALIEVE_VSDCreator_H

#include <TEveVSD.h>

class AliTPCParam;
class AliRunLoader;

#include <map>


class AliEveVSDCreator : public TEveVSD
{
  AliEveVSDCreator(const AliEveVSDCreator&);            // Not implemented
  AliEveVSDCreator& operator=(const AliEveVSDCreator&); // Not implemented

public:
  enum KineType_e { KT_Standard, KT_ProtonProton };

protected:
  void          MakeItsDigitsInfo();
  TEveMCRecCrossRef* GetGeninfo(Int_t label);
  AliTPCParam*   GetTpcParam(const TEveException& eh);

  KineType_e    mKineType;  // X{GS} 7 PhonyEnum()
  TString       mDataDir;   // X{G}
  Int_t         mEvent;     // X{G}

  Float_t       mTPCHitRes;  // X{gs}
  Float_t       mTRDHitRes;  // X{gs}

  Int_t         mDebugLevel;

  std::map<Int_t, TEveMCRecCrossRef*> mGenInfoMap; //!

public:
  AliEveVSDCreator(const Text_t* name="AliEveVSDCreator", const Text_t* title="");
  virtual ~AliEveVSDCreator() {}

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

  ClassDef(AliEveVSDCreator, 1);
}; // endclass AliEveVSDCreator

#endif
