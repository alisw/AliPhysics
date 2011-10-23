#ifndef ALITRDDIGITIZER_H
#define ALITRDDIGITIZER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Produces digits from the hits information                             //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliDigitizer.h"

class TFile;
class TF1;

class AliDigitizationInput;
class AliRunLoader;

class AliTRD;
class AliTRDdigitsManager;
class AliTRDgeometry;
class AliTRDarraySignal;
class AliTRDarrayADC;
class AliTRDmcmSim;

class AliTRDdigitizer : public AliDigitizer {

 public:

  AliTRDdigitizer();
  AliTRDdigitizer(const Text_t *name, const Text_t *title);
  AliTRDdigitizer(AliDigitizationInput* digInput, const Text_t *name, const Text_t *title);
  AliTRDdigitizer(AliDigitizationInput* digInput);
  AliTRDdigitizer(const AliTRDdigitizer &d);
  virtual             ~AliTRDdigitizer();
  AliTRDdigitizer     &operator=(const AliTRDdigitizer &d);

  virtual void         Copy(TObject &d) const;
          Bool_t       InitDetector();
          void         InitOutput(Int_t iEvent);
  virtual void         Digitize(const Option_t * option = 0);  

  virtual Bool_t       Open(const Char_t *file, Int_t nEvent = 0);
  virtual Bool_t       Open(AliRunLoader * const runLoader, Int_t nEvent = 0);
  virtual Bool_t       MakeBranch(TTree *tree) const;
  virtual Bool_t       WriteDigits() const;

  virtual void         AddSDigitsManager(AliTRDdigitsManager *manager);
  virtual void         DeleteSDigitsManager();

  virtual Bool_t       MakeDigits();

          Bool_t       SortHits(Float_t **hits, Int_t *nhit);
          Bool_t       ConvertHits(Int_t det, const Float_t * const hits, Int_t nhit, AliTRDarraySignal *signals);
          Bool_t       ConvertSignals(Int_t det, AliTRDarraySignal *signals);

          Bool_t       Digits2SDigits(AliTRDdigitsManager * const manDig, AliTRDdigitsManager * const manSDig);
          Bool_t       SDigits2Digits();
          Bool_t       MergeSDigits();
          Bool_t       ConvertSDigits();

          Bool_t       Signal2ADC(Int_t det, AliTRDarraySignal *signals);
          Bool_t       Signal2SDigits(Int_t det, AliTRDarraySignal *signals);
          Bool_t       CopyDictionary(Int_t det);
	  void         CompressOutputArrays(Int_t det);

          void         SetCompress(Int_t c = 1)                    { fCompress        = c;   }
          void         SetSDigits(Int_t v = 1)                     { fSDigits         = v;   }
          void         SetEvent(Int_t v = 0)                       { fEvent           = v;   }
          void         SetManager(AliTRDdigitsManager * const man) { fDigitsManager   = man; }
          void         SetGeometry(AliTRDgeometry * const geo)     { fGeo             = geo; }
          void         SetMergeSignalOnly(Bool_t m = kTRUE)        { fMergeSignalOnly = m;   }

  AliTRDdigitsManager *Digits() const                              { return fDigitsManager;  }

          Bool_t       GetCompress() const                         { return fCompress;       }
          Bool_t       GetSDigits() const                          { return fSDigits;        }
          Float_t      GetLorentzFactor(Float_t vdrift);

          Int_t        Diffusion(Float_t vdrift, Double_t absdriftlength
                               , Double_t &lRow, Double_t &lCol, Double_t &lTime);
          Int_t        ExB(Float_t vdrift, Double_t driftlength, Double_t &lRow);
	  void         RunDigitalProcessing(Int_t det = 0);

 protected:

  AliRunLoader        *fRunLoader;          //! Local pointer
  AliTRDdigitsManager *fDigitsManager;      //! Manager for the output digits
  AliTRDdigitsManager *fSDigitsManager;     //! Manager for the summed input s-digits
  TList               *fSDigitsManagerList; //! List of managers of input s-digits
  AliTRD              *fTRD;                //! TRD detector class
  AliTRDgeometry      *fGeo;                //! TRD geometry

  AliTRDmcmSim        *fMcmSim;             //! MCM simulation for digital processing

          Int_t        fEvent;              //! Event number
          Int_t       *fMasks;              //! Masks for the merging
          Bool_t       fCompress;           //  Switch to keep only compressed data in memory
          Bool_t       fSDigits;            //  Switch for the summable digits
          Bool_t       fMergeSignalOnly;    //  Merge only detectors that contain a signal

  ClassDef(AliTRDdigitizer,20)              //  Produces TRD-Digits

};
#endif
