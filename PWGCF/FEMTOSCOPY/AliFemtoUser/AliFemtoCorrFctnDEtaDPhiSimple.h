////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnDEtaDPhiSimple - A correlation function that analyzes      //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and pseudorapidity (eta) difference                                        //
//                                                                            //
// Authors: Adam Kisiel Adam.Kisiel@cern.ch                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCORRFCTNDETADPHISIMPLE_H
#define ALIFEMTOCORRFCTNDETADPHISIMPLE_H

#include "AliFemtoCorrFctn.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>

class AliFemtoCorrFctnDEtaDPhiSimple : public AliFemtoCorrFctn {
public:
  enum CorrectionType
  {
    kNone = 0,
    kPt = 1,
    kEta = 2
  };
  typedef enum CorrectionType ReadCorrectionType;

  enum ParticleType
  {
    kNoCorrection = 0,
    kPion = 1,
    kKaon = 2,
    kProton = 3,
    kAll = 4,
    kPionMinus = 5,
    kKaonMinus = 6,
    kProtonMinus = 7,
    kLambda = 8
  };

  AliFemtoCorrFctnDEtaDPhiSimple(const char *title, const int aPhiBins=20, const int aEtaBins=20);
  AliFemtoCorrFctnDEtaDPhiSimple(const AliFemtoCorrFctnDEtaDPhiSimple &aCorrFctn);
  virtual ~AliFemtoCorrFctnDEtaDPhiSimple();

  AliFemtoCorrFctnDEtaDPhiSimple& operator=(const AliFemtoCorrFctnDEtaDPhiSimple &aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair *aPair);
  virtual void AddMixedPair(AliFemtoPair *aPair);

  virtual void Finish();

  void SetParticleTypes(ParticleType partType1, ParticleType partType2);
  void SetParticle1Type(ParticleType partType);
  void SetParticle2Type(ParticleType partType);

  void SetReadHiddenInfo(bool read);

  void WriteHistos();
  virtual TList *GetOutputList();

  virtual AliFemtoCorrFctn *Clone() const
    { return new AliFemtoCorrFctnDEtaDPhiSimple(*this); }

protected:
  TH2D *fDPhiDEtaNumerator;   // Numerator of dEta dPhi function
  TH2D *fDPhiDEtaDenominator; // Denominator of dEta dPhi function

  TH2D *fDPhiDEtaHiddenNumerator;   // Numerator of dEta dPhi function from MC
  TH2D *fDPhiDEtaHiddenDenominator; // Denominator of dEta dPhi function from MC

  TH2D *fDPhiDEtaHiddenPrimaryNumerator;   // Numerator of dEta dPhi function from MC, physical primaries
  TH2D *fDPhiDEtaHiddenPrimaryDenominator; // Denominator of dEta dPhi function from MC, physical primaries

  TH2D *fDPhiDEtaHiddenSecWeakNumerator;   // Numerator of dEta dPhi function from MC, secondaries from weak decays
  TH2D *fDPhiDEtaHiddenSecWeakDenominator; // Denominator of dEta dPhi function from MC, secondaries from weak decays

  TH2D *fDPhiDEtaHiddenSecMatNumerator;   // Numerator of dEta dPhi function from MC, secondaries from material
  TH2D *fDPhiDEtaHiddenSecMatDenominator; // Denominator of dEta dPhi function from MC, secondaries from material

  TH2D *fDPhiDEtaHiddenPrimaryNumeratorData;   // Numerator of dEta dPhi function from MC, physical primaries, filled with data values
  TH2D *fDPhiDEtaHiddenPrimaryDenominatorData; // Denominator of dEta dPhi function from MC, physical primaries, filled with data values

  TH2D *fDPhiDEtaHiddenSecWeakNumeratorData;   // Numerator of dEta dPhi function from MC, secondaries from weak decays, filled with data values
  TH2D *fDPhiDEtaHiddenSecWeakDenominatorData; // Denominator of dEta dPhi function from MC, secondaries from weak decays, filled with data values

  TH2D *fDPhiDEtaHiddenSecMatNumeratorData;   // Numerator of dEta dPhi function from MC, secondaries from material, filled with data values
  TH2D *fDPhiDEtaHiddenSecMatDenominatorData; // Denominator of dEta dPhi function from MC, secondaries from material, filled with data values

  double fphiL;
  double fphiT;

  int fEtaBins;
  int fPhiBins;

  TString fTitle;

  ParticleType part1;
  ParticleType part2;

  bool fReadHiddenInfo;

#ifdef __ROOT__
  ClassDef(AliFemtoCorrFctnDEtaDPhiSimple, 1)
#endif
};

#endif
