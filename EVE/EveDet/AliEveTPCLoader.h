// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALIEVE_TPCLoader_H
#define ALIEVE_TPCLoader_H

#include <TEveElement.h>
#include <vector>

class AliRawReaderRoot;


class AliEveTPCData;
class AliEveTPCSector2D;
class AliEveTPCSector3D;

class AliEveTPCLoader : public TEveElementList
{
  friend class AliEveTPCLoaderEditor;

  AliEveTPCLoader(const AliEveTPCLoader&);            // Not implemented
  AliEveTPCLoader& operator=(const AliEveTPCLoader&); // Not implemented

protected:
  TString           fFile;      // File holding raw-data.
  Int_t             fEvent;     // Current event.
  Bool_t            fDoubleSR;  // Set to true for double sampling-rate.

  TString           fTPCEquipementMap; // Equipement-map file-name, if set passed to raw-reader.
  AliRawReaderRoot *fReader;           // Raw-data reader.
  AliEveTPCData    *fData;             // TPC data container.

  std::vector<AliEveTPCSector2D*> fSec2Ds; // 2D sector representations.
  std::vector<AliEveTPCSector3D*> fSec3Ds; // 3D sector representations.

  Bool_t   fSetInitSectorParams; // If true, initial parameters of 2D and 3D sectors are set from values below.
  Int_t    fInitMinTime;         // Min time for display.
  Int_t    fInitMaxTime;         // Max time for display.
  Int_t    fInitThreshold;       // Threshold.
  Int_t    fInitMaxVal;          // Maximum-signal value (all signals above mapped to saturation color).

public:
  AliEveTPCLoader(const Text_t* n="AliEveTPCLoader", const Text_t* t=0);
  virtual ~AliEveTPCLoader();

  virtual void RemoveElementLocal(TEveElement* el);
  virtual void RemoveElementsLocal();

  void SetFile(const Text_t* f) { fFile = f; }
  void SetDoubleSR(Bool_t d)    { fDoubleSR = d; }

  const Text_t* GetTPCEquipementMap() const  { return fTPCEquipementMap; }
  void SetTPCEquipementMap(const Text_t* em) { fTPCEquipementMap = em; }
  AliRawReaderRoot* GetReader()        const { return fReader; }
  void SetReader(AliRawReaderRoot* reader)   { fReader = reader; }
  AliEveTPCData* GetData() const { return fData; }
  void     SetData(AliEveTPCData* d);

  void OpenFile();
  void LoadEvent();
  void NextEvent(Bool_t rewindOnEnd=kTRUE);
  void GotoEvent(Int_t event);
  static void* LoopEvent(AliEveTPCLoader* loader);

  void UpdateSectors(Bool_t dropNonPresent=kFALSE);
  void ReloadSectors();
  void CreateSectors3D();
  void DeleteSectors3D();

  void SetInitParams(Int_t mint, Int_t maxt, Int_t thr, Int_t maxval=128);

  ClassDef(AliEveTPCLoader, 1); // Front-end for stand-alone inspection of TPC raw-data.
}; // endclass AliEveTPCLoader

#endif
