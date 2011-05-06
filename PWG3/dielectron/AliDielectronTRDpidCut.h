/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//
// TRD PID Cut Class
// Does PID either on a x% electron efficiency basis or on dE/dx
// For more information please check the implementation file
//
#ifndef ALIDIELECTRONTRDPIDCUT_H
#define ALIDIELECTRONTRDPIDCUT_H

#include <AliAnalysisCuts.h>

class AliESDtrack;
class AliVParticle;

class AliDielectronTRDpidCut : public AliAnalysisCuts{
public:
  typedef enum{
    kLQ = 0,
    kNN = 1
  } PIDMethodTRD_t;
  enum{
    kThreshParams = 24
  };
  enum{
    kHistTRDlikeBefore = 0,
    kHistTRDlikeAfter = 1,
    kHistTRDthresholds = 2,
    kHistTRDSigV1 = 3,
    kHistTRDSigV2 = 4,
    kHistOverallSpecies = 5
  };
  enum PIDbitType {kIgnore=0, kRequire, kIfAvailable};
  AliDielectronTRDpidCut();
  AliDielectronTRDpidCut(const Char_t *name);
  AliDielectronTRDpidCut(const AliDielectronTRDpidCut &ref);
  AliDielectronTRDpidCut& operator=(const AliDielectronTRDpidCut &ref);
  virtual ~AliDielectronTRDpidCut();

  virtual Bool_t InitializePID();

  Double_t GetTRDSignalV1(const AliESDtrack *track, Float_t truncation = 0.7) const;
  Double_t GetTRDSignalV2(const AliESDtrack *track, Float_t trucation = 0.7) const;

  Bool_t IsCalculateTRDSignals() const { return TestBit(kTRDsignals); }
  Bool_t IsRenormalizeElPi() const { return TestBit(kTRDrenormalize); }
  void SetPIDBitType(UInt_t pidBitType=AliDielectronTRDpidCut::kRequire) { fRequirePIDbit = pidBitType; };
  void SetPIDMethod(PIDMethodTRD_t method) { fPIDMethod = method; };
  void SetRenormalizeElPi(Bool_t doRenorm = kTRUE) { if(doRenorm) SetBit(kTRDrenormalize, kTRUE); else SetBit(kTRDrenormalize, kFALSE);}
  void SetElectronEfficiency(Double_t electronEfficiency) { fElectronEfficiency = electronEfficiency; }
  void SetThresholdParameters(Double_t electronEff, Double_t *params);
  void SetMinP(Double_t p) { fMinP = p; }
  void CalculateTRDSignals(Bool_t docalc) { SetBit(kTRDsignals, docalc); }

  Double_t GetElectronLikelihood(const AliVParticle *track) const;
  void     GetTRDmomenta(const AliVParticle *track, Double_t *mom) const;
  Double_t GetP(const AliVParticle *track) const;
  Double_t GetTRDthresholds(Double_t electronEff, Double_t p) const;
  Double_t GetChargeLayer(const AliVParticle *track, UInt_t layer) const;

  //
  //Analysis cuts interface
  //
  virtual Bool_t IsSelected(TObject* track);
  virtual Bool_t IsSelected(TList*   /* list */ ) {return kFALSE;}
  
protected:
  enum{
    kTRDsignals = BIT(16),
    kTRDdefaultThresholds = BIT(17),
    kTRDrenormalize = BIT(18)
  };
  void Copy(TObject &ref) const;
  void InitParameters();
  void InitParameters1DLQ();
  void GetParameters(Double_t electronEff, Double_t *parameters) const;
  void SetUseDefaultParameters(Bool_t useDefault = kTRUE) { SetBit(kTRDdefaultThresholds, useDefault); }
  Bool_t UseDefaultParameters() const { return TestBit(kTRDdefaultThresholds); }
  void RenormalizeElPi(const Double_t *likein, Double_t *likeout) const;

private:
  static const Double_t fgkVerySmall;                     // Check for 0
  Double_t fMinP;                                         // Minimum momentum above which TRD PID is applied
  Double_t fElectronEfficiency;                           // Cut on electron efficiency
  PIDMethodTRD_t fPIDMethod;                              // PID Method: 2D Likelihood or Neural Network
  UChar_t  fRequirePIDbit;                                //How to make use of the pid bit (see)
  Double_t fThreshParams[kThreshParams];                  // Threshold parametrisation
  ClassDef(AliDielectronTRDpidCut, 1)     // TRD electron ID class
};

#endif
