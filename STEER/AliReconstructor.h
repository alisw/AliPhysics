#ifndef ALIRECONSTRUCTOR_H
#define ALIRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// base class for reconstruction algorithm
// Derived classes should implement a default constructor and
// the virtual methods
//

#include <TObject.h>
#include <TString.h>

class TTree;
class AliRawReader;
class AliVertexer;
class AliTracker;
class AliTrackleter;
class AliESDEvent;
class AliDetectorRecoParam;
class AliRunInfo;
class AliEventInfo;
class AliESDpid;

#include "AliReconstruction.h"

class AliReconstructor: public TObject {
public:
  AliReconstructor(): TObject(), fOption(), fRunInfo(0x0), fEventInfo(0x0) {};
  virtual ~AliReconstructor() {};

  virtual void         Init() {};

  virtual Bool_t       HasDigitConversion() const {return kFALSE;};
  virtual void         ConvertDigits(AliRawReader* rawReader, TTree* digitsTree) const;

  virtual void         Reconstruct(TTree* digitsTree, TTree* clustersTree) const;
  virtual void         Reconstruct(AliRawReader* rawReader, TTree* clustersTree) const;

  virtual AliVertexer* CreateVertexer() const 
    {return NULL;}
  virtual AliTracker*  CreateTracker() const 
    {return NULL;}
  virtual AliTracker*  CreateTrackleter() const 
    {return NULL;}
  virtual AliTrackleter* CreateMultFinder() const 
    {return NULL;}

  virtual void         FillESD(TTree* digitsTree, TTree* clustersTree, 
			       AliESDEvent* esd) const;
  virtual void         FillESD(AliRawReader* rawReader, TTree* clustersTree, 
			       AliESDEvent* esd) const;

  virtual const char*  GetDetectorName() const;

  void                 SetOption(Option_t* option) {fOption = option;};
  virtual Option_t*    GetOption() const {return fOption.Data();};

  void                 SetRunInfo(AliRunInfo *runInfo) {fRunInfo = runInfo;}
  const AliRunInfo*    GetRunInfo() const {return fRunInfo;}
  void                 SetEventInfo(AliEventInfo *evInfo) {fEventInfo = evInfo;}
  const AliEventInfo*  GetEventInfo() const {return fEventInfo;}

  void                               SetRecoParam(const AliDetectorRecoParam *par);
  static const AliDetectorRecoParam* GetRecoParam(Int_t iDet);
  virtual void                 GetPidSettings(AliESDpid *esdPID);

private:

  AliReconstructor(const AliReconstructor &); // Not implemented
  AliReconstructor& operator=(const AliReconstructor &); // Not implemented
  
  TString                            fOption;                                       //! option for reconstruction
  static const AliDetectorRecoParam* fgRecoParam[AliReconstruction::kNDetectors]; //! event reconstruction parameters for all detectors
  AliRunInfo*                        fRunInfo;                                    //! pointer to the run info object
  AliEventInfo*                      fEventInfo;                                  //! pointer to the event info object

  ClassDef(AliReconstructor, 0)   // base class for reconstruction algorithms
};

#endif
