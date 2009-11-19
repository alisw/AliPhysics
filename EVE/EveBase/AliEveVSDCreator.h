// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveVSDCreator_H
#define AliEveVSDCreator_H

#include <TEveVSD.h>

#include <map>

class AliTPCParam;
class AliRunLoader;

class AliEveVSDCreator : public TEveVSD
{
public:
  AliEveVSDCreator(const Text_t* name="AliEveVSDCreator", const Text_t* title="");
  virtual ~AliEveVSDCreator() {}

  void CreateVSD(const Text_t* vsdFile="AliVSD.root");

  void CreateTrees();

  // --------------------------------------------------------------
  // Conversion functions.

  void ConvertKinematics();

  void ConvertHits();
  void ConvertAnyHits(const TString& detector, const TString& clones_class,
                      Int_t det_id, Float_t minDistSqr=0);
  void ConvertTPCHits(Int_t det_id, Float_t minDistSqr=0);

  void ConvertClusters();
  void ConvertAnyClusters(const TString& detector, const TString& branch_name,
                          Int_t det_id);
  void ConvertTPCClusters(Int_t det_id);

  void ConvertRecTracks();

  void ConvertV0();
  void ConvertKinks();

  void ConvertGenInfo();

  // --------------------------------------------------------------

  Int_t GetDebugLevel() const   { return fDebugLevel; }
  void  SetDebugLevel(Int_t dl) { fDebugLevel = dl; }

protected:
  void               MakeItsDigitsInfo();
  TEveMCRecCrossRef* GetGeninfo(Int_t label);

  Float_t       fTPCHitRes;  // Resolution for storing TPC hits.
  Float_t       fTRDHitRes;  // Resolution for storing TRD hits.

  Int_t         fDebugLevel; // Internal debug level.

  AliRunLoader *fRunLoader;  // Internal run-loader.

  std::map<Int_t, TEveMCRecCrossRef*> fGenInfoMap; // Map label to MC-Rec cross-ref data structure.

private:
  AliEveVSDCreator(const AliEveVSDCreator&);            // Not implemented
  AliEveVSDCreator& operator=(const AliEveVSDCreator&); // Not implemented

  ClassDef(AliEveVSDCreator, 0); // Create VSD file from ALICE data.
}; // endclass AliEveVSDCreator

#endif
