//////////////////////////////////////////////////////////////////////////////////////
//                                                                                  //
//  AliFemtoModelWeightGeneratorBasicLednicky - basic femtoscopic weight            //
//  generator returing a simple Lednicky weight.  Does not nearly include all       //
//  AliFemtoModelWeightGeneratorLednicky, but will be more user friendly and        //
//  will not use FORTRAN                                                            //
//                                                                                  //
//  To use:                                                                         //
//  1.  Create generator object, and set the 8 variables as arguments in the        //
//      normal constructor or via the various setters (SetIdenticalParticles(bool ) //
//      SetParamLambda(double ), ... , SetParamNorm(doube ) )                       //
//                                                                                  //
//  2.  Load generator into AliFemtoModelManager via                                //
//      AliFemtoModelManager::AcceptWeightGenerator()                               //
//                                                                                  //
//  3.  Connect manager to AliFemtoModelCorrFctn via                                //
//      AliFemtoModelCorrFctn::ConnectToManager                                     //
//                                                                                  //
//  4.  Obtain weight in Cf via AliFemtoModelManager::GetWeight(AliFemtoPair* )     //
//                                                                                  //
//  Author: Jesse Buxton jesse.thomas.buxton@cern.ch                                //
//////////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOMODELWEIGHTGENERATORBASICLEDNICKY_H
#define ALIFEMTOMODELWEIGHTGENERATORBASICLEDNICKY_H

#include "TRandom2.h"
#include "AliFemtoPair.h"
#include "AliFemtoModelWeightGenerator.h"

#include <complex>
#include <math.h>
#include "TMath.h"

#include "Faddeeva.h"

#include <limits>

const double hbarc = 0.197327;
const std::complex<double> ImI (0.,1.);

class AliFemtoModelWeightGeneratorBasicLednicky : public AliFemtoModelWeightGenerator
{
public:


  AliFemtoModelWeightGeneratorBasicLednicky();
  AliFemtoModelWeightGeneratorBasicLednicky(bool aAreIdentical, double aLambda, double aAlpha, double aRadius, double aRef0, double aImf0, double ad0, double aNorm);
  virtual ~AliFemtoModelWeightGeneratorBasicLednicky();

  AliFemtoModelWeightGeneratorBasicLednicky(const AliFemtoModelWeightGeneratorBasicLednicky &aGenerator);
  AliFemtoModelWeightGeneratorBasicLednicky& operator=(const AliFemtoModelWeightGeneratorBasicLednicky &aGenerator);

  double GetLednickyF1(double z);
  double GetLednickyF2(double z);
  double GetLednickyWeight(double aKStar);
  double GetKStarTrue(AliFemtoPair *aPair);

  virtual Double_t GenerateWeight(AliFemtoPair *aPair);

  void SetIdenticalParticles(bool aAreIdentical);
  void SetParamLambda(double aParam);
  void SetParamAlpha(double aParam);
  void SetParamRadius(double aParam);
  void SetParamRef0(double aParam);
  void SetParamImf0(double aParam);
  void SetParamd0(double aParam);
  void SetParamNorm(double aParam);


protected:
  bool fIdenticalParticles;
  double fParamLambda, fParamAlpha, fParamRadius, fParamRef0, fParamImf0, fParamd0, fParamNorm;



#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoModelWeightGeneratorBasicLednicky, 1);
  /// \endcond
#endif

};

inline void AliFemtoModelWeightGeneratorBasicLednicky::SetIdenticalParticles(bool aAreIdentical) {fIdenticalParticles = aAreIdentical;}
inline void AliFemtoModelWeightGeneratorBasicLednicky::SetParamLambda(double aParam) {fParamLambda = aParam;}
inline void AliFemtoModelWeightGeneratorBasicLednicky::SetParamAlpha(double aParam) {fParamAlpha = aParam;}
inline void AliFemtoModelWeightGeneratorBasicLednicky::SetParamRadius(double aParam) {fParamRadius = aParam;}
inline void AliFemtoModelWeightGeneratorBasicLednicky::SetParamRef0(double aParam) {fParamRef0 = aParam;}
inline void AliFemtoModelWeightGeneratorBasicLednicky::SetParamImf0(double aParam) {fParamImf0 = aParam;}
inline void AliFemtoModelWeightGeneratorBasicLednicky::SetParamd0(double aParam) {fParamd0 = aParam;}
inline void AliFemtoModelWeightGeneratorBasicLednicky::SetParamNorm(double aParam) {fParamNorm = aParam;}

#endif
