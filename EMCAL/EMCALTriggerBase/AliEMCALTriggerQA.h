/**
 * @file AliEMCALTriggerQA.h
 * @date Nov. 23, 2015
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 */
#ifndef ALIEMCALTRIGGERQA_H
#define ALIEMCALTRIGGERQA_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TNamed.h>
#include <TArrayI.h>
#include <cstring>

#include <AliEMCALTriggerConstants.h>

class AliEMCALTriggerPatchInfo;
class TCollection;
class TObjArray;
class AliEMCALTriggerFastOR;
class AliVCaloCells;
class AliEMCALGeometry;

/**
 * @class AliEMCALTriggerQA
 * @brief Class to generate EMCal trigger QA plots
 */

class AliEMCALTriggerQA : public TNamed {
public:

  typedef EMCALTrigger::EMCalTriggerType_t EMCalTriggerType_t;

  enum PatchTypes_t {
    kOnlinePatch  = 0,
    kRecalcPatch  = 1,
    kOfflinePatch = 2
  };

  struct AliEMCALCellInfo {
    AliEMCALCellInfo() : fAbsId(-1), fEnergy(0.) {}
    void Set(Int_t absId, Double_t e) { fAbsId = absId; fEnergy = e; }

    Short_t  fAbsId;
    Double_t fEnergy;
  };

  AliEMCALTriggerQA();
  AliEMCALTriggerQA(const char* name);
  AliEMCALTriggerQA(const AliEMCALTriggerQA& ref);
  virtual ~AliEMCALTriggerQA();

  void   SetDebugLevel(Int_t l)       { fDebugLevel        = l; }
  void   SetADCperBin(Int_t i)        { fADCperBin         = i; }

  Int_t  GetDebugLevel()        const { return fDebugLevel    ; }

  void   EnablePatchType(PatchTypes_t patchtype, EMCalTriggerType_t triggertype, Bool_t e);
  Bool_t IsPatchTypeEnabled(Int_t patchtype, Int_t triggertype) const;
  void   EnableHistogramsByTimeStamp(UInt_t binWidth = 600){ fTimeStampBinWidth  = binWidth   ; }
  void   SetEMCALGeometry(const AliEMCALGeometry* geom) { fGeom = geom; }

  // This is the minimum set of methods that must be implemented by derived classes
  virtual void           Init() = 0;
  virtual void           ProcessPatch(const AliEMCALTriggerPatchInfo* patch) = 0;
  virtual void           ProcessFastor(const AliEMCALTriggerFastOR* fastor, AliVCaloCells* cells = 0) = 0;
  virtual void           EventCompleted() = 0;
  virtual TCollection*   GetListOfHistograms() = 0;

  // Additional virtual methods that can optionally be overloaded
  virtual void   ExecOnce();
  virtual void   ProcessCell(const AliEMCALCellInfo& /*cell*/) {;}
  virtual void   EventTimeStamp(UInt_t timeStamp);

  // These virtual methods are implemented only for PbPb
  virtual void   ProcessBkgPatch(const AliEMCALTriggerPatchInfo* /*patch*/) {;}
  virtual void   ComputeBackground() {;}
  virtual void   GetEMCalMedian(Double_t /*median*/[3]) const {;}
  virtual void   GetDCalMedian(Double_t /*median*/[3])  const {;}
  virtual void   GetEMCalBkg(Double_t /*bkg*/[3])       const {;}
  virtual void   GetDCalBkg(Double_t /*bkg*/[3])        const {;}

  static Int_t  GetAmplitude(const AliEMCALTriggerPatchInfo* patch, Int_t itype);

  static const Int_t fgkSM = 20;
  static const Int_t fgkNPatchTypes = 3;
  static const Int_t fgkNTriggerTypes = 6;
  static const Int_t fgkNDet = 2;

  static const Int_t      fgkMaxPatchAmp[fgkNTriggerTypes];          ///< Maximum patch amplitude for the histograms
  static const TString    fgkPatchTypes[fgkNPatchTypes];             ///< Patch type names

protected:
  UInt_t                  fEnabledTriggerPatches[fgkNPatchTypes];    ///< Patch types to be plotted
  Int_t                   fFastorL0Th;                               ///< FastOR L0 threshold
  Int_t                   fFastorL1Th;                               ///< FastOR L1 threshold
  Int_t                   fADCperBin;                                ///< ADC counts per bin
  Int_t                   fDebugLevel;                               ///< Debug level
  UInt_t                  fTimeStampBinWidth;                        ///< Time stamp bin width

  const AliEMCALGeometry *fGeom;                                     //!<! Pointer to the EMCal geometry
  UInt_t                  fEventTimeStamp;                           //!<! Time stamp of the current event
  UInt_t                  fEventTimeStampBin;                        //!<! Time stamp bin

private:
  AliEMCALTriggerQA &operator=(const AliEMCALTriggerQA &); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliEMCALTriggerQA, 3);
  /// \endcond
};

#endif
