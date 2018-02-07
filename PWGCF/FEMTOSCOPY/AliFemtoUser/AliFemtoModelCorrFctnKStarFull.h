///
/// \file AliFemtoModelCorrFctnKStarFull.h
/// \authors Jesse Buxton and Andrew Kubera
///


#ifndef ALIFEMTOMODELCORRFCTNKSTARFULL_H
#define ALIFEMTOMODELCORRFCTNKSTARFULL_H

class TH1F;
class TH2F;
class TH3F;

#include "AliFemtoCorrFctn.h"
#include "AliFemtoModelCorrFctn.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoPair.h"
#include "AliFemtoThreeVector.h"
#include "AliFemtoLorentzVector.h"
#include "AliFemtoModelManager.h"

#include <limits>

/// \class AliFemtoModelCorrFctnKStarFull
/// \brief The correlation function which plots numerators and denominator
///        plots from the real monte caro
///
/// \authors: Jesse Buxton, jesse.thomas.buxton@cern.ch
///           Andrew Kubera, andrew.kubera@cern.ch
///
///

class AliFemtoModelCorrFctnKStarFull : public AliFemtoModelCorrFctn 
{
public:

  AliFemtoModelCorrFctnKStarFull();
  AliFemtoModelCorrFctnKStarFull(const char *title, int aNbins, double aKStarLo, double aKStarHi);
  virtual ~AliFemtoModelCorrFctnKStarFull();

  AliFemtoModelCorrFctnKStarFull(const AliFemtoModelCorrFctnKStarFull& aCorrFctn);
  AliFemtoModelCorrFctnKStarFull& operator=(const AliFemtoModelCorrFctnKStarFull& aCorrFctn);

  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual TList* GetOutputList();
  virtual AliFemtoModelCorrFctn* Clone() const { return new AliFemtoModelCorrFctnKStarFull(*this); }
  virtual void Write();

  /// Sets the PDG particle codes that this correlation function will check
  /// particle true.
  virtual void SetExpectedPDGCodes(const Int_t track1_pdg_code,
                                   const Int_t track2_pdg_code);
  /// Returns true if the pair's particles match the ExpectedPDGCodes provided
  virtual bool PairContainsExpectedTypes(const AliFemtoPair*);
  virtual bool PairContainsExpectedTypes(const AliFemtoModelHiddenInfo*,
                                         const AliFemtoModelHiddenInfo*);

  double CalcKStar(const AliFemtoLorentzVector&, const AliFemtoLorentzVector&);
  double GetKStarTrue(AliFemtoPair* aPair);

  int GetMotherBin(const AliFemtoModelHiddenInfo *info);

  void SetBuildBaseClassHistograms(bool aBuild);
  void SetBuildUnitWeights(bool aBuild);
  void SetBuildParentInfo(bool aBuild);
  void SetBuildTrueVsRec(bool aBuild);
  void SetBuildRotated(bool aBuild);

  //inline
  void SetRemoveMisidentified(bool aSet);


protected:
  const char* fTitle;
  int fNbinsKStar;
  double fKStarLow, fKStarHigh;

  bool fRemoveMisidentified;

  int fExpectedTrack1Code, fExpectedTrack2Code;
  bool fBuildBaseClassHistograms;
  //----- Defined in AliFemtoModelCorrFctn.h -----//
  /*
   *  TH1F *fNumeratorTrue;      // Numerator made with RECONSTRUCTED k* of pairs from SAME event
                                 //   Weight from fManager, if exists (otherwise weight = 1)
   *  TH1F *fNumeratorFake;      // Numerator made with RECONSTRUCTED k* of pairs from MIXED events
                                 //   Weight from fManager, if exists (otherwise weight = 1)
   *  TH1F *fDenominator;        // Denominator made with RECONSTRUCTED k* of pairs from MIXED events
                                 //   Weight = 1 always 

   *  TH1F *fNumeratorTrueIdeal; // Numerator made with TRUE k* of pairs from SAME event
                                 //   Weight from fManager, if exists (otherwise weight = 1)
   *  TH1F *fNumeratorFakeIdeal; // Numerator made with TRUE k* of pairs from MIXED events
                                 //   Weight from fManager, if exists (otherwise weight = 1)
   *  TH1F *fDenominatorIdeal;   // Denominator made with TRUE k* of pairs from MIXED events
                                 //   Weight = 1 always
  */
  //----- END: Defined in AliFemtoModelCorrFctn.h -----//

  bool fBuildUnitWeights;
  TH1F *fNumTrueUnitWeights;      // Numerator made with RECONSTRUCTED k* of pairs from SAME event with unit weights
  TH1F *fNumTrueIdealUnitWeights; // Numerator made with TRUE k* of pairs from SAME event with unit weights

  bool fBuildParentInfo;
  TH3F *fNumTrueIdealwParentInfo;  //Same as fNumeratorTrueIdeal but additionally binned according to parent PIDs
  TH3F *fDenIdealwParentInfo;      //Same as fDenominatorIdeal but additionally binned according to parent PIDs

  bool fBuildTrueVsRec;
  TH2F *fKTrueVsKRecSame;           // 2D histogram of k*_{True} vs k*_{Reconstructed} of pairs from SAME event
  TH2F *fKTrueVsKRecMixed;          // 2D histogram of k*_{True} vs k*_{Reconstructed} of pairs from MIXED events

  bool fBuildRotated;
  TH2F *fKTrueVsKRecRotSame;        // fKTrueVsKRecSame rotated by 45 degrees
  TH2F *fKTrueVsKRecRotMixed;       // fKTrueVsKRecMixed rotated by 45 degrees


#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoModelCorrFctnKStarFull, 1);
  /// \endcond
#endif

};

inline void AliFemtoModelCorrFctnKStarFull::SetRemoveMisidentified(bool aSet) {fRemoveMisidentified=aSet;}


#endif
