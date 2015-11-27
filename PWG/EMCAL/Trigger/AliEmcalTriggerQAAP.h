/**
 * @file AliEmcalTriggerQAAP.h
 * @date Nov. 23, 2015
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 */
#ifndef ALIEMCALTRIGGERQAAP_H
#define ALIEMCALTRIGGERQAAP_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TNamed.h>
#include <TArrayI.h>

class AliEmcalTriggerPatchInfoAPV1;
class THashList;
class TObjArray;
class AliEmcalTriggerFastORAP;

/**
 * @class AliEMCALTriggerQA
 * @brief Class to generate EMCal trigger QA plots
 */

class AliEmcalTriggerQAAP : public TNamed {
public:

  typedef EmcalTriggerAP::EMCalTriggerType_t EMCalTriggerType_t;

  enum PatchTypes_t {
    kOnlinePatch,
    kRecalcPatch,
    kOfflinePatch
  };

  AliEmcalTriggerQAAP();
  AliEmcalTriggerQAAP(const char* name);
  virtual ~AliEmcalTriggerQAAP();

  void   SetDebugLevel(Int_t l)       { fDebugLevel        = l; }
  void   SetBkgPatchType(Int_t t)     { fBkgPatchType      = t; }

  Int_t  GetDebugLevel()        const { return fDebugLevel    ; }
  Int_t  GetBkgPatchType()      const { return fBkgPatchType  ; }

  void   EnablePatchType(PatchTypes_t type, Bool_t e = kTRUE);

  void   Init();
  void   ProcessPatch(AliEmcalTriggerPatchInfoAPV1* patch);
  void   ProcessFastor(AliEmcalTriggerFastORAP* fastor);
  void   EventCompleted();

  THashList* GetListOfHistograms()  { return fHistos; }

protected:

  static const Int_t      fgkMaxPatchAmp[6];            ///< Maximum patch amplitude for the histograms
  static const TString    fgkPatchTypes[3];             ///< Patch type names

  Bool_t                  fEnabledPatchTypes[3];        ///< Patch types to be plotted
  Int_t                   fFastorL0Th;                  ///< FastOR L0 threshold
  Int_t                   fFastorL1Th;                  ///< FastOR L1 threshold
  Int_t                   fBkgPatchType;                ///< Background patch type

  Int_t                   fDebugLevel;                  ///< Debug level

  TArrayI                 fBkgADCAmpEMCal[3];           //!<! ADC EMCal amplitudes (0=online, 1=offline) of background patches (will be reset each event)
  Int_t                   fNBkgPatchesEMCal[3];         //!<! Number of processed background patches (will be reset each event)
  Int_t                   fMaxPatchEMCal[6][3];         //!<! EMCal max ADC amplitude (0=online, 1=offline) (will be reset each event)
  TArrayI                 fBkgADCAmpDCal[3];            //!<! ADC DCal amplitudes (0=online, 1=offline) of background patches (will be reset each event)
  Int_t                   fNBkgPatchesDCal[3];          //!<! Number of processed background patches (will be reset each event)
  Int_t                   fMaxPatchDCal[6][3];          //!<! DCal max ADC amplitude (0=online, 1=offline) (will be reset each event)

  THashList              *fHistos;                      //!<! Histograms for QA

private:
  void CreateTProfile(const char *name, const char *title, int nbins, double xmin, double xmax);
  void CreateTH1(const char *name, const char *title, int nbins, double xmin, double xmax);
  void CreateTH2(const char *name, const char *title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax);
  void CreateTH3(const char *name, const char *title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax, int nbinsz, double zmin, double zmax);
  void FillTProfile(const char *name, double x, double y, double weight = 1.);
  void FillTH1(const char *hname, double x, double weight = 1.);
  void FillTH2(const char *hname, double x, double y, double weight = 1.);
  void FillTH3(const char *hname, double x, double y, double z, double weight = 1.);

  TObject *FindObject(const char *name) const;
  virtual TObject *FindObject(const TObject *obj) const;

  AliEmcalTriggerQAAP(const AliEmcalTriggerQAAP &);
  AliEmcalTriggerQAAP &operator=(const AliEmcalTriggerQAAP &);

  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerQAAP, 1);
  /// \endcond
};

#endif
