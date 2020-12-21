////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelWeightGenerator - abstract base class for femtoscopic       ///
/// weight generator                                                         ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOMODELWEIGHTGENERATOR_H
#define ALIFEMTOMODELWEIGHTGENERATOR_H

#include "TRandom2.h"
class AliFemtoPair;

class AliFemtoModelWeightGenerator
{
 public:
  AliFemtoModelWeightGenerator();
  AliFemtoModelWeightGenerator(const AliFemtoModelWeightGenerator &aModel);
  virtual ~AliFemtoModelWeightGenerator();
  AliFemtoModelWeightGenerator& operator=(const AliFemtoModelWeightGenerator &aModel);
  virtual Double_t GenerateWeight(AliFemtoPair *aPair) = 0;

  virtual void     SetPairType(Int_t aPairType);
  virtual void     SetPairTypeFromPair(AliFemtoPair *aPair);
  virtual Int_t    GetPairType() const;
  virtual Int_t    GetPairTypeFromPair(AliFemtoPair *aPair);

  virtual Double_t GetKStar() const;
  virtual Double_t GetKStarOut() const;
  virtual Double_t GetKStarSide() const;
  virtual Double_t GetKStarLong() const;
  virtual Double_t GetRStar() const;
  virtual Double_t GetRStarOut() const;
  virtual Double_t GetRStarSide() const;
  virtual Double_t GetRStarLong() const;

  virtual AliFemtoModelWeightGenerator* Clone() const;

  static Int_t PionPlusPionPlus();
  static Int_t PionPlusPionMinus();
  static Int_t KaonPlusKaonPlus();
  static Int_t KaonPlusKaonMinus();
  static Int_t ProtonProton();
  static Int_t ProtonAntiproton();
  static Int_t PionPlusKaonPlus();
  static Int_t PionPlusKaonMinus();
  static Int_t PionPlusProton();
  static Int_t PionPlusAntiproton();
  static Int_t KaonPlusProton();
  static Int_t KaonPlusAntiproton();
  static Int_t PairTypeNone();
  static Int_t LambdaLambda();
  static Int_t AntilambdaAntilambda();
  static Int_t LambdaAntilambda();

 protected:
  static const Int_t fgkPairTypeNone;        ///< no pair type set - read from model
  static const Int_t fgkPionPlusPionPlus;    ///< identical pion pair
  static const Int_t fgkPionPlusPionMinus;   ///< non-identical pion pair
  static const Int_t fgkKaonPlusKaonPlus;    ///< identical kaon pair
  static const Int_t fgkKaonPlusKaonMinus;   ///< non-identical kaon pair
  static const Int_t fgkProtonProton;        ///< identical proton pair
  static const Int_t fgkProtonAntiproton;    ///< non-identical proton pair
  static const Int_t fgkPionPlusKaonPlus;    ///< same-charge pion kaon pair
  static const Int_t fgkPionPlusKaonMinus;   ///< opposite-charge pion kaon pair
  static const Int_t fgkPionPlusProton;      ///< same-charge pion proton pair
  static const Int_t fgkPionPlusAntiproton;  ///< opposite-chare pion proton pair
  static const Int_t fgkKaonPlusProton;      ///< same-charge kaon proton pair
  static const Int_t fgkKaonPlusAntiproton;  ///< opposite-charge kaon proton pair
  static const Int_t fgkLambdaLambda;        ///< same-type lambdas
  static const Int_t fgkAntilambdaAntilambda; ///< same-type antilambdas
  static const Int_t fgkLambdaAntilambda;    ///< non-same-type lambdas

  Int_t fPairType;     ///< Type of the pair for which the calculation is done

  Double_t fKStarOut;  ///< relative momentum out component in PRF
  Double_t fKStarSide; ///< relative momentum side component in PRF
  Double_t fKStarLong; ///< relative momentum long component in PRF
  Double_t fKStar;     ///< relative momentum magnitude

  Double_t fRStarOut;  ///< relative separation out component in PRF
  Double_t fRStarSide; ///< relative separation side component in PRF
  Double_t fRStarLong; ///< relative separation long component in PRF
  Double_t fRStar;     ///< relative separation magnitude
 private:

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoModelWeightGenerator, 1);
  /// \endcond
#endif

};

inline Double_t AliFemtoModelWeightGenerator::GetKStar() const { return fKStar; }
inline Double_t AliFemtoModelWeightGenerator::GetKStarOut() const { return fKStarOut; }
inline Double_t AliFemtoModelWeightGenerator::GetKStarSide() const { return fKStarSide; }
inline Double_t AliFemtoModelWeightGenerator::GetKStarLong() const { return fKStarLong; }
inline Double_t AliFemtoModelWeightGenerator::GetRStar() const { return fRStar; }
inline Double_t AliFemtoModelWeightGenerator::GetRStarOut() const { return fRStarOut; }
inline Double_t AliFemtoModelWeightGenerator::GetRStarSide() const { return fRStarSide; }
inline Double_t AliFemtoModelWeightGenerator::GetRStarLong() const { return fRStarLong; }

inline Int_t AliFemtoModelWeightGenerator::PairTypeNone() { return fgkPairTypeNone; }
inline Int_t AliFemtoModelWeightGenerator::PionPlusPionPlus() { return fgkPionPlusPionPlus; }
inline Int_t AliFemtoModelWeightGenerator::PionPlusPionMinus() { return fgkPionPlusPionMinus; }
inline Int_t AliFemtoModelWeightGenerator::KaonPlusKaonPlus() { return fgkKaonPlusKaonPlus; }
inline Int_t AliFemtoModelWeightGenerator::KaonPlusKaonMinus() { return fgkKaonPlusKaonMinus; }
inline Int_t AliFemtoModelWeightGenerator::ProtonProton() { return fgkProtonProton; }
inline Int_t AliFemtoModelWeightGenerator::ProtonAntiproton() { return fgkProtonAntiproton; }
inline Int_t AliFemtoModelWeightGenerator::PionPlusKaonPlus() { return fgkPionPlusKaonPlus; }
inline Int_t AliFemtoModelWeightGenerator::PionPlusKaonMinus() { return fgkPionPlusKaonMinus; }
inline Int_t AliFemtoModelWeightGenerator::PionPlusProton() { return fgkPionPlusProton; }
inline Int_t AliFemtoModelWeightGenerator::PionPlusAntiproton() { return fgkPionPlusAntiproton; }
inline Int_t AliFemtoModelWeightGenerator::KaonPlusProton() { return fgkKaonPlusProton; }
inline Int_t AliFemtoModelWeightGenerator::KaonPlusAntiproton() { return fgkKaonPlusAntiproton; }
inline Int_t AliFemtoModelWeightGenerator::LambdaLambda() { return fgkLambdaLambda; }
inline Int_t AliFemtoModelWeightGenerator::AntilambdaAntilambda() { return fgkAntilambdaAntilambda; }
inline Int_t AliFemtoModelWeightGenerator::LambdaAntilambda() { return fgkLambdaAntilambda; }

#endif


