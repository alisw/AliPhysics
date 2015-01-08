#ifndef ALIEMCALTRIGGERBINNINGCOMPONENT_H
#define ALIEMCALTRIGGERBINNINGCOMPONENT_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel
#include <TArrayD.h>
#include <TNamed.h>

class TObjArray;

namespace EMCalTriggerPtAnalysis {

class AliEMCalTriggerBinningDimension : public TNamed{
public:
  AliEMCalTriggerBinningDimension():
    TNamed(),
    fBinning()
  {}
  AliEMCalTriggerBinningDimension(const char *name):
    TNamed(name, ""),
    fBinning()
  {}
  AliEMCalTriggerBinningDimension(const char *name, int nbins, double *binning):
    TNamed(name, ""),
    fBinning(nbins+1, binning)
  {}
  AliEMCalTriggerBinningDimension(const char *name, const TArrayD &binning):
    TNamed(name, ""),
    fBinning(binning.GetSize(), binning.GetArray())
  {}
  ~AliEMCalTriggerBinningDimension() {}

  void Set(int nbins, double *binning) { fBinning.Set(nbins+1, binning); }
  void Set(const TArrayD &binning) { fBinning = binning; }
  const double *GetBinLimits() const { return fBinning.GetArray(); }
  int GetNumberOfBins() const { return fBinning.GetSize() - 1; }
  virtual void Print(Option_t *option="") const;

private:
  TArrayD fBinning;             // Bin limits

  ClassDef(AliEMCalTriggerBinningDimension, 1);
};

class AliEMCalTriggerBinningComponent: public TObject {
public:
  AliEMCalTriggerBinningComponent();
  AliEMCalTriggerBinningComponent(const AliEMCalTriggerBinningComponent &ref);
  AliEMCalTriggerBinningComponent &operator=(const AliEMCalTriggerBinningComponent &ref);
  virtual ~AliEMCalTriggerBinningComponent();

  AliEMCalTriggerBinningDimension *GetBinning(const char *name) const;
  void SetBinning(const char *dimname, int nbins, double *binning);
  void SetBinning(const char *dimname, const TArrayD &binning);

private:
  TObjArray       *fDimensions;           // List of binnings (dimensions)

  ClassDef(AliEMCalTriggerBinningComponent, 1);
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGERBINNINGCOMPONENT_H */
