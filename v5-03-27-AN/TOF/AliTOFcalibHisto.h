#ifndef ALITOFCALIBHISTO_H
#define ALITOFCALIBHISTO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// *
// *
// * AliTOFcalibHisto - class to handle TOF calibration histograms,
// *                    map histograms and more
// *
// *
// * autors:   Roberto Preghenella (R+)
// * concacts: preghenella@bo.infn.it
// *
// *

#include "TObject.h"
#include "TString.h"

class TH1D;
class TFile;
class AliESDtrack;

class AliTOFcalibHisto :
public TObject
{

 public:

  enum ECorrection_t {
    kDDLBCcorr,
    kAmphenolCableCorr,
    kFlatCableCorr,
    kInterfaceCardCorr,
    kDDLdelayCorr,
    kHPTDCdelayCorr,
    kFEAchDelayCorr,
    kFEAdelayCorr,
    kTRMdelayCorr,
    kICdelayCorr,
    kStripDelayCorr,
    kIndexDelayCorr,
    kTimeSlewingCorr,
    kNcorrections
  };
  
  enum ECalibConst_t {
    kLHCperiod,
    kAmphenolCableDelay,
    kFlatCableDelay,
    kInterfaceCardDelay,
    kNcalibConsts
  };

  enum ECalibMap_t {
    /* main index */
    kIndex,
    /* EO indices */
    kDDL, 
    kTRM, 
    kChain, 
    kTDC, 
    kChannel, 
    /* DO indices */
    kSector, 
    kPlate, 
    kStrip, 
    kSectorStrip, 
    kPadZ, 
    kPadX, 
    kPad,
    kInterfaceCardIndex,
    /* calib constants */
    kDDLBCshift,
    kFlatCableLength,
    kInterfaceCardLength,
    kAmphenolCableLength,
    /* number of histos */
    kNcalibMaps
  };

  enum ECalibPar_t {
    kDDLdelayPar,
    kHPTDCdelayPar,
    kLeftFEAchDelayPar,
    kRightFEAchDelayPar,
    kFEAdelayPar,
    kICdelayPar,
    kTRMdelayPar,
    kStripDelayPar,
    kIndexDelayPar,
    kTimeSlewingPar,
    kNcalibPars
  };

  enum ECalibStat_t {
    kStripStat,
    kNcalibStats
  };

  AliTOFcalibHisto(); /* default constructor */
  virtual ~AliTOFcalibHisto(); /* default destructor */

  /* getters */
  static const Char_t *GetCalibHistoFileName() {return fgCalibHistoFileName.Data();}; /* get calib histo file name */
  static const Char_t *GetCalibParFileName() {return fgCalibParFileName.Data();}; /* get calib par file name */
  static const Char_t *GetCalibStatFileName() {return fgCalibStatFileName.Data();}; /* get calib stat file name */

  /* setters */
  static void SetCalibHistoFileName(const Char_t *value) {fgCalibHistoFileName = value;}; /* set calib histo file name */
  static void SetCalibParFileName(const Char_t *value) {fgCalibParFileName = value;}; /* set calib par file name */
  static void SetCalibStatFileName(const Char_t *value) {fgCalibStatFileName = value;}; /* set calib stat file name */
  static void SetCableCorrectionFlag(Int_t i, Bool_t flag) {if (i < kNcorrections) fgCableCorrectionFlag[i] = flag;}; // set cable correction flag
  static void SetFullCorrectionFlag(Int_t i, Bool_t flag) {if (i < kNcorrections) fgFullCorrectionFlag[i] = flag;}; // set full correction flag

  /* methods */
  static Int_t GetIndexEO(Int_t ddl, Int_t trm, Int_t chain, Int_t tdc, Int_t channel) {return (channel + 8 * tdc + 120 * chain + 240 * trm + 2400 * ddl);}; /* get index EO */
  void LoadCalibHisto(); /* load calib histo */
  void LoadCalibPar(); /* load calib par */
  void LoadCalibStat(); /* load calib stat */
  void WriteCalibHisto(); /* write calib histo */
  void WriteCalibStat(); /* write calib stat */
  Double_t GetCalibConst(Int_t histo) {return GetHisto(fCalibConst[histo], 0);}; /* get calib const */
  Double_t GetCalibMap(Int_t histo, Int_t index) {return GetHisto(fCalibMap[histo], index);}; /* get calib map */
  Double_t GetCalibPar(Int_t histo, Int_t index) {return GetHisto(fCalibPar[histo], index);}; /* get calib par */
  Bool_t GetCalibStat(Int_t histo, Int_t index) {return GetHisto(fCalibStat[histo], index) != 0. ? kTRUE : kFALSE;}; /* get calib stat */
  Double_t GetCorrection(Int_t corr, Int_t index, Double_t tot = 0.); /* get correction */
  Double_t GetNominalCorrection(Int_t index); /* get nominal correction */
  void ApplyNominalCorrection(AliESDtrack *track); /* apply nominal corrections */

  Double_t GetCableCorrection(Int_t index); /* get cable correction */
  Double_t GetFullCorrection(Int_t index, Double_t tot = 0.); /* get full correction */

  Bool_t GetStatus(Int_t index); /* get status */
  
 private:

  static TFile *fgCalibHistoFile; /* calib histo file */
  static TFile *fgCalibParFile; /* calib par file */
  static TFile *fgCalibStatFile; /* calib stat file */

  static TString fgCalibHistoFileName; /* calib histo file name */
  static TString fgCalibParFileName; /* calib par file name */
  static TString fgCalibStatFileName; /* calib par file name */
  static const TString fgkCalibConstName[kNcalibConsts]; // calib const name array */
  static const TString fgkCalibMapName[kNcalibMaps]; // calib map name array */
  static const TString fgkCalibParName[kNcalibPars]; // calib par name array */
  static const TString fgkCalibStatName[kNcalibStats]; // calib stat name array */

  static const Double_t fgkLHCperiod; /* LHC clock period */
  static const Double_t fgkAmphenolCableDelay; /* Amphenol cable delay */
  static const Double_t fgkFlatCableDelay; /* flat cable delay */
  static const Double_t fgkInterfaceCardDelay; /* interface card delay */

  static const Int_t fgkNchannels; /* number of readout channels (DO) */
  static const Int_t fgkNchannelsEO; /* number of readout channels (EO) */
  static const Int_t fgkDDLBCshift[72]; /* DDL BC shifts due to TTC fibers */
  static const Double_t fgkFlatCableLength[91]; /* strip flat-cable length */
  static const Double_t fgkInterfaceCardLength[48]; /* interface card length */

  static Bool_t fgCableCorrectionFlag[kNcorrections]; // cable correction flag
  static Bool_t fgFullCorrectionFlag[kNcorrections]; // full correction flag

  static const Int_t fgkStripStat[18][91]; // strip status

  TH1D *fCalibConst[kNcalibConsts]; // calib const array
  TH1D *fCalibMap[kNcalibMaps]; // calib map array
  TH1D *fCalibPar[kNcalibPars]; // calib par array
  TH1D *fCalibStat[kNcalibStats]; // calib stat array

  /* methods */
  void LoadHisto(TFile * const file, TH1D **histo, const Char_t *name); /* create histo */
  void CreateHisto(TH1D **histo, const Char_t *name, Int_t size); /* create histo */
  void WriteHisto(TFile *file, TH1D *histo); /* write histo */
  void SetHisto(TH1D *histo, Int_t index, Double_t value); /* set histo */
  Double_t GetHisto(TH1D *histo, Int_t index); /* get histo */
  
  AliTOFcalibHisto(const AliTOFcalibHisto &source) : TObject(source) {}; /* copy constructor */
  AliTOFcalibHisto &operator=(const AliTOFcalibHisto &) {return *this;}; /* operator= */

  ClassDef(AliTOFcalibHisto, 1);
  
};

#endif /* ALITOFCALIBHISTO_H */
