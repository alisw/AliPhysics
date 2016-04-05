/**
 * @file AliEMCALTriggerOnlineQAPbPb.h
 * @date Nov. 23, 2015
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 */
#ifndef ALIEMCALONLINETRIGGERQAPBPB_H
#define ALIEMCALONLINETRIGGERQAPBPB_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TArrayI.h>
#include <cstring>
#include <THashList.h>

#include "AliEMCALTriggerQA.h"

class AliEMCALTriggerPatchInfo;
class THashList;
class TObjArray;
class AliEMCALTriggerFastOR;
class AliVCaloCells;

/**
 * @class AliEMCALTriggerQAPbPb
 * @brief Class to generate EMCal trigger QA plots in PbPb collision mode
 */

class AliEMCALTriggerOnlineQAPbPb : public AliEMCALTriggerQA {
public:

  AliEMCALTriggerOnlineQAPbPb();
  AliEMCALTriggerOnlineQAPbPb(const char* name);
  AliEMCALTriggerOnlineQAPbPb(const AliEMCALTriggerOnlineQAPbPb& triggerQA);
  virtual ~AliEMCALTriggerOnlineQAPbPb();

  void   SetBkgPatchType(Int_t t)     { fBkgPatchType      = t; }
  Int_t  GetBkgPatchType()      const { return fBkgPatchType  ; }

  void   Init();
  void   ProcessPatch(const AliEMCALTriggerPatchInfo* patch);
  void   ProcessBkgPatch(const AliEMCALTriggerPatchInfo* patch);
  void   ProcessFastor(const AliEMCALTriggerFastOR* fastor, AliVCaloCells* /*cells*/ = 0);
  void   EventCompleted();
  void   ComputeBackground();

  void   GetEMCalMedian(Double_t median[3]) const { memcpy(median, fMedianEMCal, sizeof(Double_t)*3); }
  void   GetDCalMedian(Double_t median[3])  const { memcpy(median, fMedianDCal, sizeof(Double_t)*3) ; }
  void   GetEMCalBkg(Double_t bkg[3])       const { memcpy(bkg, fBkgEMCal, sizeof(Double_t)*3)      ; }
  void   GetDCalBkg(Double_t bkg[3])        const { memcpy(bkg, fBkgDCal, sizeof(Double_t)*3)       ; }

  TCollection* GetListOfHistograms()  { return fHistos; }

protected:
  Int_t                   fBkgPatchType;                ///< Background patch type

  TArrayI                 fBkgADCAmpEMCal[3];           //!<! ADC EMCal amplitudes (0=online, 1=offline) of background patches (will be reset each event)
  Int_t                   fNBkgPatchesEMCal[3];         //!<! Number of processed background patches (will be reset each event)
  Int_t                   fMaxPatchEMCal[6][3];         //!<! EMCal max ADC amplitude (0=online, 1=offline) (will be reset each event)
  TArrayI                 fBkgADCAmpDCal[3];            //!<! ADC DCal amplitudes (0=online, 1=offline) of background patches (will be reset each event)
  Int_t                   fNBkgPatchesDCal[3];          //!<! Number of processed background patches (will be reset each event)
  Int_t                   fMaxPatchDCal[6][3];          //!<! DCal max ADC amplitude (0=online, 1=offline) (will be reset each event)
  Int_t                   fPatchAreas[6];               //!<! Patch sizes retrieved directly during the patch processing
  Double_t                fMedianEMCal[3];              //!<! Median of background patches in the EMCal
  Double_t                fMedianDCal[3];               //!<! Median of background patches in the DCal
  Double_t                fBkgEMCal[3];                 //!<! Background in the EMCal
  Double_t                fBkgDCal[3];                  //!<! Background in the DCal
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

  AliEMCALTriggerOnlineQAPbPb &operator=(const AliEMCALTriggerOnlineQAPbPb &);

  /// \cond CLASSIMP
  ClassDef(AliEMCALTriggerOnlineQAPbPb, 1);
  /// \endcond
};

#endif
