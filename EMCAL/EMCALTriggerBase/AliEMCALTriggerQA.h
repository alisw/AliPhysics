/**
 * @file AliEMCALTriggerQA.h
 * @date Nov. 1, 2015
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 */
#ifndef ALIEMCALTRIGGERQA_H
#define ALIEMCALTRIGGERQA_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TNamed.h>
#include <TArrayI.h>

class AliEMCALTriggerPatchInfo;
class THashList;
class TObjArray;

/**
 * @class AliEMCALTriggerQA
 * @brief Class to generate EMCal trigger QA plots
 */

class AliEMCALTriggerQA : public TNamed {
public:
  /**
   * \enum TriggerMakerTriggerType_t
   * \brief Definition of different trigger patch types
   *
   * This enumeration defines the different trigger patch types
   * processed by the trigger maker. Each trigger patch type has
   * a certain patch size and therefore a certain length and
   * geometric center
   */
  enum TriggerQATriggerType_t {
    kTMEMCalJet = 0,        ///< EMCAL Jet trigger (16x16)
    kTMEMCalGamma = 1,      ///< EMCAL Gamma trigger (2x2)
    kTMEMCalLevel0 = 2,     ///< EMCAL Level0 patches (2x2)
    kTMEMCalRecalcJet = 3,  ///< EMCAL Jet patches, recalculated (16x16)
    kTMEMCalRecalcGamma = 4,///< EMCAL Gamma patches, recalculated (2x2)
    kTMEMCalBackground = 5  ///< EMCAL Background patches (8x8)
  };
  
  AliEMCALTriggerQA();
  AliEMCALTriggerQA(const char* name);
  virtual ~AliEMCALTriggerQA();

  void   SetDebugLevel(Int_t l)       { fDebugLevel        = l; }
  void   SetBkgPatchType(Int_t t)     { fBkgPatchType      = t; }
  void   SetBkgOfflineAmp(Bool_t b)   { fBkgOfflineAmp     = b; }
  
  Int_t  GetDebugLevel()        const { return fDebugLevel    ; }
  Int_t  GetBkgPatchType()      const { return fBkgPatchType  ; }
  Bool_t GetBkgOfflineAmp()     const { return fBkgOfflineAmp ; }
    
  void Init();
  void ProcessPatch(AliEMCALTriggerPatchInfo* patch);
  void EventCompleted();

  THashList* GetListOfHistograms()  { return fHistos     ; }

protected:

  Int_t                   GetPatchType(AliEMCALTriggerPatchInfo* patch);

  static const Int_t      fgkNTriggerTypes = 6;                       ///< Number of trigger types
  static const TString    fgkTriggerTypeNames[fgkNTriggerTypes];      ///< Histogram name tags

  Bool_t                  fBkgOfflineAmp;                             ///< Use offline amplitudes to calculate the background
  Int_t                   fBkgPatchType;                              ///< Background patch type
  
  Int_t                   fDebugLevel;                                ///< Debug level

  TArrayI                 fADCAmpEMCal[fgkNTriggerTypes];             //!<! ADC EMCal amplitudes (will be reset each event)
  Int_t                   fNPatchesEMCal[fgkNTriggerTypes];           //!<! Number of processed patches (will be reset each event)
  Int_t                   fMaxPatchEMCal[fgkNTriggerTypes];           //!<! Number of processed patches (will be reset each event)
  TArrayI                 fADCAmpDCal[fgkNTriggerTypes];              //!<! ADC DCal amplitudes (will be reset each event)
  Int_t                   fNPatchesDCal[fgkNTriggerTypes];            //!<! Number of processed patches (will be reset each event)
  Int_t                   fMaxPatchDCal[fgkNTriggerTypes];            //!<! Number of processed patches (will be reset each event)

  THashList              *fHistos;                                    //!<! Histograms for QA

 private:
  void CreateTH1(const char *name, const char *title, int nbins, double xmin, double xmax);
  void CreateTH2(const char *name, const char *title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax);
  void CreateTH3(const char *name, const char *title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax, int nbinsz, double zmin, double zmax);
  void SetObject(TObject * const o, const char *group = "/");
  void FillTH1(const char *hname, double x, double weight = 1.);
  void FillTH2(const char *hname, double x, double y, double weight = 1.);
  void FillTH3(const char *hname, double x, double y, double z, double weight = 1.);

  TObject *FindObject(const char *name) const;
  virtual TObject *FindObject(const TObject *obj) const;

  AliEMCALTriggerQA(const AliEMCALTriggerQA &);
  AliEMCALTriggerQA &operator=(const AliEMCALTriggerQA &);
  
  /// \cond CLASSIMP
  ClassDef(AliEMCALTriggerQA, 1);
  /// \endcond
};

#endif
