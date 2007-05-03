////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelWeightGenerator - abstract base class for femtoscopic       ///
/// weight generator                                                         ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoModelWeightGenerator_hh
#define AliFemtoModelWeightGenerator_hh

#include "TRandom2.h"
#include "AliFemtoPair.h"

class AliFemtoModelWeightGenerator 
{
 public:
  AliFemtoModelWeightGenerator();
  AliFemtoModelWeightGenerator(const AliFemtoModelWeightGenerator &aModel);
  virtual ~AliFemtoModelWeightGenerator();
  virtual Double_t GenerateWeight(AliFemtoPair *aPair) = 0;

  virtual void     SetPairType(Int_t aPairType);
  virtual void     SetPairTypeFromPair(AliFemtoPair *aPair);
  virtual Int_t    GetPairType();

  virtual Double_t GetKStar();
  virtual Double_t GetKStarOut();
  virtual Double_t GetKStarSide();
  virtual Double_t GetKStarLong();
  virtual Double_t GetRStar();
  virtual Double_t GetRStarOut();
  virtual Double_t GetRStarSide();
  virtual Double_t GetRStarLong();

  virtual AliFemtoModelWeightGenerator* Clone() const;

  static const Int_t kPionPlusPionPlus;
  static const Int_t kPionPlusPionMinus;
  static const Int_t kKaonPlusKaonPlus;
  static const Int_t kKaonPlusKaonMinus;
  static const Int_t kProtonProton;
  static const Int_t kProtonAntiproton;
  static const Int_t kPionPlusKaonPlus;
  static const Int_t kPionPlusKaonMinus;
  static const Int_t kPionPlusProton;
  static const Int_t kPionPlusAntiproton;
  static const Int_t kKaonPlusProton;
  static const Int_t kKaonPlusAntiproton;

 protected:
  Int_t fPairType;

  Double_t fKStarOut;  // relative momentum out component in PRF
  Double_t fKStarSide; // relative momentum side component in PRF
  Double_t fKStarLong; // relative momentum long component in PRF
  Double_t fKStar;     // relative momentum magnitude

  Double_t fRStarOut;  // relative separation out component in PRF
  Double_t fRStarSide; // relative separation side component in PRF
  Double_t fRStarLong; // relative separation long component in PRF
  Double_t fRStar;     // relative separation magnitude
 private:
  
#ifdef __ROOT__
  ClassDef(AliFemtoModelWeightGenerator, 1)
#endif

    };
  
inline Double_t AliFemtoModelWeightGenerator::GetKStar() { return fKStar; }
inline Double_t AliFemtoModelWeightGenerator::GetKStarOut() { return fKStarOut; }
inline Double_t AliFemtoModelWeightGenerator::GetKStarSide() { return fKStarSide; }
inline Double_t AliFemtoModelWeightGenerator::GetKStarLong() { return fKStarLong; }
inline Double_t AliFemtoModelWeightGenerator::GetRStar() { return fRStar; }
inline Double_t AliFemtoModelWeightGenerator::GetRStarOut() { return fRStarOut; }
inline Double_t AliFemtoModelWeightGenerator::GetRStarSide() { return fRStarSide; }
inline Double_t AliFemtoModelWeightGenerator::GetRStarLong() { return fRStarLong; }


#endif


