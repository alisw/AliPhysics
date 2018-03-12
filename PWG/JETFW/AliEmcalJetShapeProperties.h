#ifndef ALIEMCALJETSHAPEPROPERTIES_H
#define ALIEMCALJETSHAPEPROPERTIES_H

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TArrayF.h>
#include <Rtypes.h>
#include <TString.h>

/// \class AliEmcalJetShapeProperties
/// \brief This class contains the derivative subtraction operators for jet shapes
///
/// This class was designed to contain the derivative subtraction operators
/// for jet shapes, as an extension of the AliEmcalJet class.
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>, University of Birmingham
/// \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
/// \date Apr 30, 2016
class AliEmcalJetShapeProperties {

public:

  AliEmcalJetShapeProperties();
  AliEmcalJetShapeProperties(const AliEmcalJetShapeProperties &jetshape);
  AliEmcalJetShapeProperties& operator=(const AliEmcalJetShapeProperties &jetshape);

  //jet shape derivatives
  //jet mass
  void              SetFirstDerivative(Double_t d)                             { fJetShapeMassFirstDer = d                       ; }
  void              SetSecondDerivative(Double_t d)                            { fJetShapeMassSecondDer = d                      ; }
  void              SetFirstOrderSubtracted(Double_t d)                        { fJetShapeMassFirstSub = d                       ; }
  void              SetSecondOrderSubtracted(Double_t d)                       { fJetShapeMassSecondSub = d                      ; }
  Double_t          GetFirstDerivative()                                 const { return fJetShapeMassFirstDer                    ; }
  Double_t          GetSecondDerivative()                                const { return fJetShapeMassSecondDer                   ; }
  Double_t          GetFirstOrderSubtracted()                            const { return fJetShapeMassFirstSub                    ; }
  Double_t          GetSecondOrderSubtracted()                           const { return fJetShapeMassSecondSub                   ; }

  //jet structure function
  TArrayF           GetGRNumerator()                                     const { return fGRNumerator                             ; }
  TArrayF           GetGRDenominator()                                   const { return fGRDenominator                           ; }
  TArrayF           GetGRNumeratorSub()                                  const { return fGRNumeratorSub                          ; }
  TArrayF           GetGRDenominatorSub()                                const { return fGRDenominatorSub                        ; }
  void              AddGRNumAt(Float_t num, Int_t idx)                         { fGRNumerator.AddAt(num, idx)                    ; }
  void              AddGRDenAt(Float_t den, Int_t idx)                         { fGRDenominator.AddAt(den, idx)                  ; }
  void              SetGRNumSize(UInt_t s)                                     { fGRNumerator.Set(s)                             ; }
  void              SetGRDenSize(UInt_t s)                                     { fGRDenominator.Set(s)                           ; }

  void              AddGRNumSubAt(Float_t num, Int_t idx)                      { fGRNumeratorSub.AddAt(num, idx)                 ; }
  void              AddGRDenSubAt(Float_t den, Int_t idx)                      { fGRDenominatorSub.AddAt(den, idx)               ; }
  void              SetGRNumSubSize(UInt_t s)                                  { fGRNumeratorSub.Set(s)                          ; }
  void              SetGRDenSubSize(UInt_t s)                                  { fGRDenominatorSub.Set(s)                        ; }
  void              PrintGR();

  //Angularity
  void              SetFirstDerivativeAngularity(Double_t d)                   { fJetShapeAngularityFirstDer = d                 ; }
  void              SetSecondDerivativeAngularity(Double_t d)                  { fJetShapeAngularitySecondDer = d                ; }
  void              SetFirstOrderSubtractedAngularity(Double_t d)              { fJetShapeAngularityFirstSub = d                 ; }
  void              SetSecondOrderSubtractedAngularity(Double_t d)             { fJetShapeAngularitySecondSub = d                ; }
  Double_t          GetFirstDerivativeAngularity()                       const { return fJetShapeAngularityFirstDer              ; }
  Double_t          GetSecondDerivativeAngularity()                      const { return fJetShapeAngularitySecondDer             ; }
  Double_t          GetFirstOrderSubtractedAngularity()                  const { return fJetShapeAngularityFirstSub              ; }
  Double_t          GetSecondOrderSubtractedAngularity()                 const { return fJetShapeAngularitySecondSub             ; }

  //pTD
  void              SetFirstDerivativepTD(Double_t d)                          { fJetShapepTDFirstDer = d                        ; }
  void              SetSecondDerivativepTD(Double_t d)                         { fJetShapepTDSecondDer = d                       ; }
  void              SetFirstOrderSubtractedpTD(Double_t d)                     { fJetShapepTDFirstSub = d                        ; }
  void              SetSecondOrderSubtractedpTD(Double_t d)                    { fJetShapepTDSecondSub = d                       ; }
  Double_t          GetFirstDerivativepTD()                              const { return fJetShapepTDFirstDer                     ; }
  Double_t          GetSecondDerivativepTD()                             const { return fJetShapepTDSecondDer                    ; }
  Double_t          GetFirstOrderSubtractedpTD()                         const { return fJetShapepTDFirstSub                     ; }
  Double_t          GetSecondOrderSubtractedpTD()                        const { return fJetShapepTDSecondSub                    ; }

  //Circularity
  void              SetFirstDerivativeCircularity(Double_t d)                  { fJetShapeCircularityFirstDer = d                ; }
  void              SetSecondDerivativeCircularity(Double_t d)                 { fJetShapeCircularitySecondDer = d               ; }
  void              SetFirstOrderSubtractedCircularity(Double_t d)             { fJetShapeCircularityFirstSub = d                ; }
  void              SetSecondOrderSubtractedCircularity(Double_t d)            { fJetShapeCircularitySecondSub = d               ; }
  Double_t          GetFirstDerivativeCircularity()                      const { return fJetShapeCircularityFirstDer             ; }
  Double_t          GetSecondDerivativeCircularity()                     const { return fJetShapeCircularitySecondDer            ; }
  Double_t          GetFirstOrderSubtractedCircularity()                 const { return fJetShapeCircularityFirstSub             ; }
  Double_t          GetSecondOrderSubtractedCircularity()                const { return fJetShapeCircularitySecondSub            ; }

  //Sigma2
  void              SetFirstDerivativeSigma2(Double_t d)                       { fJetShapeSigma2FirstDer = d                     ; }
  void              SetSecondDerivativeSigma2(Double_t d)                      { fJetShapeSigma2SecondDer = d                    ; }
  void              SetFirstOrderSubtractedSigma2(Double_t d)                  { fJetShapeSigma2FirstSub = d                     ; }
  void              SetSecondOrderSubtractedSigma2(Double_t d)                 { fJetShapeSigma2SecondSub = d                    ; }
  Double_t          GetFirstDerivativeSigma2()                           const { return fJetShapeSigma2FirstDer                  ; }
  Double_t          GetSecondDerivativeSigma2()                          const { return fJetShapeSigma2SecondDer                 ; }
  Double_t          GetFirstOrderSubtractedSigma2()                      const { return fJetShapeSigma2FirstSub                  ; }
  Double_t          GetSecondOrderSubtractedSigma2()                     const { return fJetShapeSigma2SecondSub                 ; }


  //number of contituents
  void              SetFirstDerivativeConstituent(Double_t d)                  { fJetShapeConstituentFirstDer = d                ; }
  void              SetSecondDerivativeConstituent(Double_t d)                 { fJetShapeConstituentSecondDer = d               ; }
  void              SetFirstOrderSubtractedConstituent(Double_t d)             { fJetShapeConstituentFirstSub = d                ; }
  void              SetSecondOrderSubtractedConstituent(Double_t d)            { fJetShapeConstituentSecondSub = d               ; }
  Double_t          GetFirstDerivativeConstituent()                      const { return fJetShapeConstituentFirstDer             ; }
  Double_t          GetSecondDerivativeConstituent()                     const { return fJetShapeConstituentSecondDer            ; }
  Double_t          GetFirstOrderSubtractedConstituent()                 const { return fJetShapeConstituentFirstSub             ; }
  Double_t          GetSecondOrderSubtractedConstituent()                const { return fJetShapeConstituentSecondSub            ; }

  //leading minus subleading constituent
  void              SetFirstDerivativeLeSub(Double_t d)                        { fJetShapeLeSubFirstDer = d                      ; }
  void              SetSecondDerivativeLeSub(Double_t d)                       { fJetShapeLeSubSecondDer = d                     ; }
  void              SetFirstOrderSubtractedLeSub(Double_t d)                   { fJetShapeLeSubFirstSub = d                      ; }
  void              SetSecondOrderSubtractedLeSub(Double_t d)                  { fJetShapeLeSubSecondSub = d                     ; }
  Double_t          GetFirstDerivativeLeSub()                            const { return fJetShapeLeSubFirstDer                   ; }
  Double_t          GetSecondDerivativeLeSub()                           const { return fJetShapeLeSubSecondDer                  ; }
  Double_t          GetFirstOrderSubtractedLeSub()                       const { return fJetShapeLeSubFirstSub                   ; }
  Double_t          GetSecondOrderSubtractedLeSub()                      const { return fJetShapeLeSubSecondSub                  ; }

  //1subjettiness_kt
  void              SetFirstDerivative1subjettiness_kt(Double_t d)             { fJetShape1subjettinessktFirstDer = d                      ; }
  void              SetSecondDerivative1subjettiness_kt(Double_t d)            { fJetShape1subjettinessktSecondDer = d                     ; }
  void              SetFirstOrderSubtracted1subjettiness_kt(Double_t d)        { fJetShape1subjettinessktFirstSub = d                      ; }
  void              SetSecondOrderSubtracted1subjettiness_kt(Double_t d)       { fJetShape1subjettinessktSecondSub = d                     ; }
  Double_t          GetFirstDerivative1subjettiness_kt()                 const { return fJetShape1subjettinessktFirstDer                   ; }
  Double_t          GetSecondDerivative1subjettiness_kt()                const { return fJetShape1subjettinessktSecondDer                  ; }
  Double_t          GetFirstOrderSubtracted1subjettiness_kt()            const { return fJetShape1subjettinessktFirstSub                   ; }
  Double_t          GetSecondOrderSubtracted1subjettiness_kt()           const { return fJetShape1subjettinessktSecondSub                  ; }

  //2subjettiness_kt
  void              SetFirstDerivative2subjettiness_kt(Double_t d)             { fJetShape2subjettinessktFirstDer = d                      ; }
  void              SetSecondDerivative2subjettiness_kt(Double_t d)            { fJetShape2subjettinessktSecondDer = d                     ; }
  void              SetFirstOrderSubtracted2subjettiness_kt(Double_t d)        { fJetShape2subjettinessktFirstSub = d                      ; }
  void              SetSecondOrderSubtracted2subjettiness_kt(Double_t d)       { fJetShape2subjettinessktSecondSub = d                     ; }
  Double_t          GetFirstDerivative2subjettiness_kt()                 const { return fJetShape2subjettinessktFirstDer                   ; }
  Double_t          GetSecondDerivative2subjettiness_kt()                const { return fJetShape2subjettinessktSecondDer                  ; }
  Double_t          GetFirstOrderSubtracted2subjettiness_kt()            const { return fJetShape2subjettinessktFirstSub                   ; }
  Double_t          GetSecondOrderSubtracted2subjettiness_kt()           const { return fJetShape2subjettinessktSecondSub                  ; }

  //3subjettiness_kt
  void              SetFirstDerivative3subjettiness_kt(Double_t d)             { fJetShape3subjettinessktFirstDer = d                      ; }
  void              SetSecondDerivative3subjettiness_kt(Double_t d)            { fJetShape3subjettinessktSecondDer = d                     ; }
  void              SetFirstOrderSubtracted3subjettiness_kt(Double_t d)        { fJetShape3subjettinessktFirstSub = d                      ; }
  void              SetSecondOrderSubtracted3subjettiness_kt(Double_t d)       { fJetShape3subjettinessktSecondSub = d                     ; }
  Double_t          GetFirstDerivative3subjettiness_kt()                 const { return fJetShape3subjettinessktFirstDer                   ; }
  Double_t          GetSecondDerivative3subjettiness_kt()                const { return fJetShape3subjettinessktSecondDer                  ; }
  Double_t          GetFirstOrderSubtracted3subjettiness_kt()            const { return fJetShape3subjettinessktFirstSub                   ; }
  Double_t          GetSecondOrderSubtracted3subjettiness_kt()           const { return fJetShape3subjettinessktSecondSub                  ; }

  //OpeningAngle_kt
  void              SetFirstDerivativeOpeningAngle_kt(Double_t d)             { fJetShapeOpeningAnglektFirstDer = d                      ; }
  void              SetSecondDerivativeOpeningAngle_kt(Double_t d)            { fJetShapeOpeningAnglektSecondDer = d                     ; }
  void              SetFirstOrderSubtractedOpeningAngle_kt(Double_t d)        { fJetShapeOpeningAnglektFirstSub = d                      ; }
  void              SetSecondOrderSubtractedOpeningAngle_kt(Double_t d)       { fJetShapeOpeningAnglektSecondSub = d                     ; }
  Double_t          GetFirstDerivativeOpeningAngle_kt()                 const { return fJetShapeOpeningAnglektFirstDer                   ; }
  Double_t          GetSecondDerivativeOpeningAngle_kt()                const { return fJetShapeOpeningAnglektSecondDer                  ; }
  Double_t          GetFirstOrderSubtractedOpeningAngle_kt()            const { return fJetShapeOpeningAnglektFirstSub                   ; }
  Double_t          GetSecondOrderSubtractedOpeningAngle_kt()           const { return fJetShapeOpeningAnglektSecondSub                  ; }

  //1subjettiness_ca
  void              SetFirstDerivative1subjettiness_ca(Double_t d)             { fJetShape1subjettinesscaFirstDer = d                      ; }
  void              SetSecondDerivative1subjettiness_ca(Double_t d)            { fJetShape1subjettinesscaSecondDer = d                     ; }
  void              SetFirstOrderSubtracted1subjettiness_ca(Double_t d)        { fJetShape1subjettinesscaFirstSub = d                      ; }
  void              SetSecondOrderSubtracted1subjettiness_ca(Double_t d)       { fJetShape1subjettinesscaSecondSub = d                     ; }
  Double_t          GetFirstDerivative1subjettiness_ca()                 const { return fJetShape1subjettinesscaFirstDer                   ; }
  Double_t          GetSecondDerivative1subjettiness_ca()                const { return fJetShape1subjettinesscaSecondDer                  ; }
  Double_t          GetFirstOrderSubtracted1subjettiness_ca()            const { return fJetShape1subjettinesscaFirstSub                   ; }
  Double_t          GetSecondOrderSubtracted1subjettiness_ca()           const { return fJetShape1subjettinesscaSecondSub                  ; }

  //2subjettiness_ca
  void              SetFirstDerivative2subjettiness_ca(Double_t d)             { fJetShape2subjettinesscaFirstDer = d                      ; }
  void              SetSecondDerivative2subjettiness_ca(Double_t d)            { fJetShape2subjettinesscaSecondDer = d                     ; }
  void              SetFirstOrderSubtracted2subjettiness_ca(Double_t d)        { fJetShape2subjettinesscaFirstSub = d                      ; }
  void              SetSecondOrderSubtracted2subjettiness_ca(Double_t d)       { fJetShape2subjettinesscaSecondSub = d                     ; }
  Double_t          GetFirstDerivative2subjettiness_ca()                 const { return fJetShape2subjettinesscaFirstDer                   ; }
  Double_t          GetSecondDerivative2subjettiness_ca()                const { return fJetShape2subjettinesscaSecondDer                  ; }
  Double_t          GetFirstOrderSubtracted2subjettiness_ca()            const { return fJetShape2subjettinesscaFirstSub                   ; }
  Double_t          GetSecondOrderSubtracted2subjettiness_ca()           const { return fJetShape2subjettinesscaSecondSub                  ; }

  //OpeningAngle_ca
  void              SetFirstDerivativeOpeningAngle_ca(Double_t d)             { fJetShapeOpeningAnglecaFirstDer = d                      ; }
  void              SetSecondDerivativeOpeningAngle_ca(Double_t d)            { fJetShapeOpeningAnglecaSecondDer = d                     ; }
  void              SetFirstOrderSubtractedOpeningAngle_ca(Double_t d)        { fJetShapeOpeningAnglecaFirstSub = d                      ; }
  void              SetSecondOrderSubtractedOpeningAngle_ca(Double_t d)       { fJetShapeOpeningAnglecaSecondSub = d                     ; }
  Double_t          GetFirstDerivativeOpeningAngle_ca()                 const { return fJetShapeOpeningAnglecaFirstDer                   ; }
  Double_t          GetSecondDerivativeOpeningAngle_ca()                const { return fJetShapeOpeningAnglecaSecondDer                  ; }
  Double_t          GetFirstOrderSubtractedOpeningAngle_ca()            const { return fJetShapeOpeningAnglecaFirstSub                   ; }
  Double_t          GetSecondOrderSubtractedOpeningAngle_ca()           const { return fJetShapeOpeningAnglecaSecondSub                  ; }

  //1subjettiness_akt02
  void              SetFirstDerivative1subjettiness_akt02(Double_t d)             { fJetShape1subjettinessakt02FirstDer = d                      ; }
  void              SetSecondDerivative1subjettiness_akt02(Double_t d)            { fJetShape1subjettinessakt02SecondDer = d                     ; }
  void              SetFirstOrderSubtracted1subjettiness_akt02(Double_t d)        { fJetShape1subjettinessakt02FirstSub = d                      ; }
  void              SetSecondOrderSubtracted1subjettiness_akt02(Double_t d)       { fJetShape1subjettinessakt02SecondSub = d                     ; }
  Double_t          GetFirstDerivative1subjettiness_akt02()                 const { return fJetShape1subjettinessakt02FirstDer                   ; }
  Double_t          GetSecondDerivative1subjettiness_akt02()                const { return fJetShape1subjettinessakt02SecondDer                  ; }
  Double_t          GetFirstOrderSubtracted1subjettiness_akt02()            const { return fJetShape1subjettinessakt02FirstSub                   ; }
  Double_t          GetSecondOrderSubtracted1subjettiness_akt02()           const { return fJetShape1subjettinessakt02SecondSub                  ; }

  //2subjettiness_akt02
  void              SetFirstDerivative2subjettiness_akt02(Double_t d)             { fJetShape2subjettinessakt02FirstDer = d                      ; }
  void              SetSecondDerivative2subjettiness_akt02(Double_t d)            { fJetShape2subjettinessakt02SecondDer = d                     ; }
  void              SetFirstOrderSubtracted2subjettiness_akt02(Double_t d)        { fJetShape2subjettinessakt02FirstSub = d                      ; }
  void              SetSecondOrderSubtracted2subjettiness_akt02(Double_t d)       { fJetShape2subjettinessakt02SecondSub = d                     ; }
  Double_t          GetFirstDerivative2subjettiness_akt02()                 const { return fJetShape2subjettinessakt02FirstDer                   ; }
  Double_t          GetSecondDerivative2subjettiness_akt02()                const { return fJetShape2subjettinessakt02SecondDer                  ; }
  Double_t          GetFirstOrderSubtracted2subjettiness_akt02()            const { return fJetShape2subjettinessakt02FirstSub                   ; }
  Double_t          GetSecondOrderSubtracted2subjettiness_akt02()           const { return fJetShape2subjettinessakt02SecondSub                  ; }

  //OpeningAngle_akt02
  void              SetFirstDerivativeOpeningAngle_akt02(Double_t d)             { fJetShapeOpeningAngleakt02FirstDer = d                      ; }
  void              SetSecondDerivativeOpeningAngle_akt02(Double_t d)            { fJetShapeOpeningAngleakt02SecondDer = d                     ; }
  void              SetFirstOrderSubtractedOpeningAngle_akt02(Double_t d)        { fJetShapeOpeningAngleakt02FirstSub = d                      ; }
  void              SetSecondOrderSubtractedOpeningAngle_akt02(Double_t d)       { fJetShapeOpeningAngleakt02SecondSub = d                     ; }
  Double_t          GetFirstDerivativeOpeningAngle_akt02()                 const { return fJetShapeOpeningAngleakt02FirstDer                   ; }
  Double_t          GetSecondDerivativeOpeningAngle_akt02()                const { return fJetShapeOpeningAngleakt02SecondDer                  ; }
  Double_t          GetFirstOrderSubtractedOpeningAngle_akt02()            const { return fJetShapeOpeningAngleakt02FirstSub                   ; }
  Double_t          GetSecondOrderSubtractedOpeningAngle_akt02()           const { return fJetShapeOpeningAngleakt02SecondSub                  ; }

  //1subjettiness_onepassca
  void              SetFirstDerivative1subjettiness_onepassca(Double_t d)             { fJetShape1subjettinessonepasscaFirstDer = d                      ; }
  void              SetSecondDerivative1subjettiness_onepassca(Double_t d)            { fJetShape1subjettinessonepasscaSecondDer = d                     ; }
  void              SetFirstOrderSubtracted1subjettiness_onepassca(Double_t d)        { fJetShape1subjettinessonepasscaFirstSub = d                      ; }
  void              SetSecondOrderSubtracted1subjettiness_onepassca(Double_t d)       { fJetShape1subjettinessonepasscaSecondSub = d                     ; }
  Double_t          GetFirstDerivative1subjettiness_onepassca()                 const { return fJetShape1subjettinessonepasscaFirstDer                   ; }
  Double_t          GetSecondDerivative1subjettiness_onepassca()                const { return fJetShape1subjettinessonepasscaSecondDer                  ; }
  Double_t          GetFirstOrderSubtracted1subjettiness_onepassca()            const { return fJetShape1subjettinessonepasscaFirstSub                   ; }
  Double_t          GetSecondOrderSubtracted1subjettiness_onepassca()           const { return fJetShape1subjettinessonepasscaSecondSub                  ; }

  //2subjettiness_onepassca
  void              SetFirstDerivative2subjettiness_onepassca(Double_t d)             { fJetShape2subjettinessonepasscaFirstDer = d                      ; }
  void              SetSecondDerivative2subjettiness_onepassca(Double_t d)            { fJetShape2subjettinessonepasscaSecondDer = d                     ; }
  void              SetFirstOrderSubtracted2subjettiness_onepassca(Double_t d)        { fJetShape2subjettinessonepasscaFirstSub = d                      ; }
  void              SetSecondOrderSubtracted2subjettiness_onepassca(Double_t d)       { fJetShape2subjettinessonepasscaSecondSub = d                     ; }
  Double_t          GetFirstDerivative2subjettiness_onepassca()                 const { return fJetShape2subjettinessonepasscaFirstDer                   ; }
  Double_t          GetSecondDerivative2subjettiness_onepassca()                const { return fJetShape2subjettinessonepasscaSecondDer                  ; }
  Double_t          GetFirstOrderSubtracted2subjettiness_onepassca()            const { return fJetShape2subjettinessonepasscaFirstSub                   ; }
  Double_t          GetSecondOrderSubtracted2subjettiness_onepassca()           const { return fJetShape2subjettinessonepasscaSecondSub                  ; }

  //OpeningAngle_onepassca
  void              SetFirstDerivativeOpeningAngle_onepassca(Double_t d)             { fJetShapeOpeningAngleonepasscaFirstDer = d                      ; }
  void              SetSecondDerivativeOpeningAngle_onepassca(Double_t d)            { fJetShapeOpeningAngleonepasscaSecondDer = d                     ; }
  void              SetFirstOrderSubtractedOpeningAngle_onepassca(Double_t d)        { fJetShapeOpeningAngleonepasscaFirstSub = d                      ; }
  void              SetSecondOrderSubtractedOpeningAngle_onepassca(Double_t d)       { fJetShapeOpeningAngleonepasscaSecondSub = d                     ; }
  Double_t          GetFirstDerivativeOpeningAngle_onepassca()                 const { return fJetShapeOpeningAngleonepasscaFirstDer                   ; }
  Double_t          GetSecondDerivativeOpeningAngle_onepassca()                const { return fJetShapeOpeningAngleonepasscaSecondDer                  ; }
  Double_t          GetFirstOrderSubtractedOpeningAngle_onepassca()            const { return fJetShapeOpeningAngleonepasscaFirstSub                   ; }
  Double_t          GetSecondOrderSubtractedOpeningAngle_onepassca()           const { return fJetShapeOpeningAngleonepasscaSecondSub                  ; }

  //SoftDrop
  void              SetSoftDropZg(Double_t d)                                 { fSoftDropZg = d                      ; }
  void              SetSoftDropdR(Double_t d)                                 { fSoftDropdR = d                      ; }
  void              SetSoftDropPtfrac(Double_t d)                             { fSoftDropPtfrac = d                  ; }
  void              SetSoftDropDropCount(Int_t d)                             { fSoftDropDropCount = d               ; }
  Double_t          GetSoftDropZg()                                     const { return fSoftDropZg                   ; }
  Double_t          GetSoftDropdR()                                     const { return fSoftDropdR                   ; }
  Double_t          GetSoftDropPtfrac()                                 const { return fSoftDropPtfrac               ; }
  Int_t             GetSoftDropDropCount()                              const { return fSoftDropDropCount;           ; }

  void              PrintGR() const;

protected:

  Double_t          fJetShapeMassFirstDer;                   //!<!   result from shape derivatives for jet mass: 1st derivative
  Double_t          fJetShapeMassSecondDer;                  //!<!   result from shape derivatives for jet mass: 2nd derivative
  Double_t          fJetShapeMassFirstSub;                   //!<!   result from shape derivatives for jet mass: 1st order subtracted
  Double_t          fJetShapeMassSecondSub;                  //!<!   result from shape derivatives for jet mass: 2nd order subtracted

  TArrayF           fGRNumerator;                            //!<!   array with angular structure function numerator
  TArrayF           fGRDenominator;                          //!<!   array with angular structure function denominator
  TArrayF           fGRNumeratorSub;                         //!<!   array with angular structure function numerator
  TArrayF           fGRDenominatorSub;                       //!<!   array with angular structure function denominator

  Double_t          fJetShapeAngularityFirstDer;             //!<!   result from shape derivatives for jet Angularity: 1st derivative
  Double_t          fJetShapeAngularitySecondDer;            //!<!   result from shape derivatives for jet Angularity: 2nd derivative
  Double_t          fJetShapeAngularityFirstSub;             //!<!   result from shape derivatives for jet Angularity: 1st order subtracted
  Double_t          fJetShapeAngularitySecondSub;            //!<!   result from shape derivatives for jet Angularity: 2nd order subtracted

  Double_t          fJetShapepTDFirstDer;                    //!<!   result from shape derivatives for jet pTD: 1st derivative
  Double_t          fJetShapepTDSecondDer;                   //!<!   result from shape derivatives for jet pTD: 2nd derivative
  Double_t          fJetShapepTDFirstSub;                    //!<!   result from shape derivatives for jet pTD: 1st order subtracted
  Double_t          fJetShapepTDSecondSub;                   //!<!   result from shape derivatives for jet pTD: 2nd order subtracted

  Double_t          fJetShapeCircularityFirstDer;            //!<!   result from shape derivatives for jet circularity: 1st derivative
  Double_t          fJetShapeCircularitySecondDer;           //!<!   result from shape derivatives for jet circularity: 2nd derivative
  Double_t          fJetShapeCircularityFirstSub;            //!<!   result from shape derivatives for jet circularity: 1st order subtracted
  Double_t          fJetShapeCircularitySecondSub;           //!<!   result from shape derivatives for jetcircularity: 2nd order subtracted

  Double_t          fJetShapeSigma2FirstDer;                 //!<!   result from shape derivatives for jet sigma2: 1st derivative
  Double_t          fJetShapeSigma2SecondDer;                //!<!   result from shape derivatives for jet sigma2: 2nd derivative
  Double_t          fJetShapeSigma2FirstSub;                 //!<!   result from shape derivatives for jet sigma2: 1st order subtracted
  Double_t          fJetShapeSigma2SecondSub;                //!<!   result from shape derivatives for jetsigma2: 2nd order subtracted

  Double_t          fJetShapeConstituentFirstDer;            //!<!   result from shape derivatives for jet const: 1st derivative
  Double_t          fJetShapeConstituentSecondDer;           //!<!   result from shape derivatives for jet const: 2nd derivative
  Double_t          fJetShapeConstituentFirstSub;            //!<!   result from shape derivatives for jet const: 1st order subtracted
  Double_t          fJetShapeConstituentSecondSub;           //!<!   result from shape derivatives for jet const: 2nd order subtracted

  Double_t          fJetShapeLeSubFirstDer;                  //!<!   result from shape derivatives for jet LeSub: 1st derivative
  Double_t          fJetShapeLeSubSecondDer;                 //!<!   result from shape derivatives for jet LeSub: 2nd derivative
  Double_t          fJetShapeLeSubFirstSub;                  //!<!   result from shape derivatives for jet LeSub: 1st order subtracted
  Double_t          fJetShapeLeSubSecondSub;                 //!<!   result from shape derivatives for jet LeSub: 2nd order subtracted

  Double_t          fJetShape1subjettinessktFirstDer;        //!<!   result from shape derivatives for jet 1subjettiness_kt: 1st derivative
  Double_t          fJetShape1subjettinessktSecondDer;       //!<!   result from shape derivatives for jet 1subjettiness_kt: 2nd derivative
  Double_t          fJetShape1subjettinessktFirstSub;        //!<!   result from shape derivatives for jet 1subjettiness_kt: 1st order subtracted
  Double_t          fJetShape1subjettinessktSecondSub;       //!<!   result from shape derivatives for jet 1subjettiness_kt: 2nd order subtracted
  
  Double_t          fJetShape2subjettinessktFirstDer;        //!<!   result from shape derivatives for jet 2subjettiness_kt: 1st derivative
  Double_t          fJetShape2subjettinessktSecondDer;       //!<!   result from shape derivatives for jet 2subjettiness_kt: 2nd derivative
  Double_t          fJetShape2subjettinessktFirstSub;        //!<!   result from shape derivatives for jet 2subjettiness_kt: 1st order subtracted
  Double_t          fJetShape2subjettinessktSecondSub;       //!<!   result from shape derivatives for jet 2subjettiness_kt: 2nd order subtracted

  Double_t          fJetShape3subjettinessktFirstDer;        //!<!   result from shape derivatives for jet 3subjettiness_kt: 1st derivative
  Double_t          fJetShape3subjettinessktSecondDer;       //!<!   result from shape derivatives for jet 3subjettiness_kt: 2nd derivative
  Double_t          fJetShape3subjettinessktFirstSub;        //!<!   result from shape derivatives for jet 3subjettiness_kt: 1st order subtracted
  Double_t          fJetShape3subjettinessktSecondSub;       //!<!   result from shape derivatives for jet 3subjettiness_kt: 2nd order subtracted

  Double_t          fJetShapeOpeningAnglektFirstDer;         //!<!   result from shape derivatives for jet OpeningAngle_kt: 1st derivative
  Double_t          fJetShapeOpeningAnglektSecondDer;        //!<!   result from shape derivatives for jet OpeningAngle_kt: 2nd derivative
  Double_t          fJetShapeOpeningAnglektFirstSub;         //!<!   result from shape derivatives for jet OpeningAngle_kt: 1st order subtracted
  Double_t          fJetShapeOpeningAnglektSecondSub;        //!<!   result from shape derivatives for jet OpeningAngle_kt: 2nd order subtracted

  Double_t          fJetShape1subjettinesscaFirstDer;        //!<!   result from shape derivatives for jet 1subjettiness_ca: 1st derivative
  Double_t          fJetShape1subjettinesscaSecondDer;       //!<!   result from shape derivatives for jet 1subjettiness_ca: 2nd derivative
  Double_t          fJetShape1subjettinesscaFirstSub;        //!<!   result from shape derivatives for jet 1subjettiness_ca: 1st order subtracted
  Double_t          fJetShape1subjettinesscaSecondSub;       //!<!   result from shape derivatives for jet 1subjettiness_ca: 2nd order subtracted
  
  Double_t          fJetShape2subjettinesscaFirstDer;        //!<!   result from shape derivatives for jet 2subjettiness_ca: 1st derivative
  Double_t          fJetShape2subjettinesscaSecondDer;       //!<!   result from shape derivatives for jet 2subjettiness_ca: 2nd derivative
  Double_t          fJetShape2subjettinesscaFirstSub;        //!<!   result from shape derivatives for jet 2subjettiness_ca: 1st order subtracted
  Double_t          fJetShape2subjettinesscaSecondSub;       //!<!   result from shape derivatives for jet 2subjettiness_ca: 2nd order subtracted

  Double_t          fJetShapeOpeningAnglecaFirstDer;         //!<!   result from shape derivatives for jet OpeningAngle_ca: 1st derivative
  Double_t          fJetShapeOpeningAnglecaSecondDer;        //!<!   result from shape derivatives for jet OpeningAngle_ca: 2nd derivative
  Double_t          fJetShapeOpeningAnglecaFirstSub;         //!<!   result from shape derivatives for jet OpeningAngle_ca: 1st order subtracted
  Double_t          fJetShapeOpeningAnglecaSecondSub;        //!<!   result from shape derivatives for jet OpeningAngle_ca: 2nd order subtracted

  Double_t          fJetShape1subjettinessakt02FirstDer;        //!<!   result from shape derivatives for jet 1subjettiness_akt02: 1st derivative
  Double_t          fJetShape1subjettinessakt02SecondDer;       //!<!   result from shape derivatives for jet 1subjettiness_akt02: 2nd derivative
  Double_t          fJetShape1subjettinessakt02FirstSub;        //!<!   result from shape derivatives for jet 1subjettiness_akt02: 1st order subtracted
  Double_t          fJetShape1subjettinessakt02SecondSub;       //!<!   result from shape derivatives for jet 1subjettiness_akt02: 2nd order subtracted
  
  Double_t          fJetShape2subjettinessakt02FirstDer;        //!<!   result from shape derivatives for jet 2subjettiness_akt02: 1st derivative
  Double_t          fJetShape2subjettinessakt02SecondDer;       //!<!   result from shape derivatives for jet 2subjettiness_akt02: 2nd derivative
  Double_t          fJetShape2subjettinessakt02FirstSub;        //!<!   result from shape derivatives for jet 2subjettiness_akt02: 1st order subtracted
  Double_t          fJetShape2subjettinessakt02SecondSub;       //!<!   result from shape derivatives for jet 2subjettiness_akt02: 2nd order subtracted

  Double_t          fJetShapeOpeningAngleakt02FirstDer;         //!<!   result from shape derivatives for jet OpeningAngle_akt02: 1st derivative
  Double_t          fJetShapeOpeningAngleakt02SecondDer;        //!<!   result from shape derivatives for jet OpeningAngle_akt02: 2nd derivative
  Double_t          fJetShapeOpeningAngleakt02FirstSub;         //!<!   result from shape derivatives for jet OpeningAngle_akt02: 1st order subtracted
  Double_t          fJetShapeOpeningAngleakt02SecondSub;        //!<!   result from shape derivatives for jet OpeningAngle_akt02: 2nd order subtracted

  Double_t          fJetShape1subjettinessonepasscaFirstDer;        //!<!   result from shape derivatives for jet 1subjettiness_onepassca: 1st derivative
  Double_t          fJetShape1subjettinessonepasscaSecondDer;       //!<!   result from shape derivatives for jet 1subjettiness_onepassca: 2nd derivative
  Double_t          fJetShape1subjettinessonepasscaFirstSub;        //!<!   result from shape derivatives for jet 1subjettiness_onepassca: 1st order subtracted
  Double_t          fJetShape1subjettinessonepasscaSecondSub;       //!<!   result from shape derivatives for jet 1subjettiness_onepassca: 2nd order subtracted
  
  Double_t          fJetShape2subjettinessonepasscaFirstDer;        //!<!   result from shape derivatives for jet 2subjettiness_onepassca: 1st derivative
  Double_t          fJetShape2subjettinessonepasscaSecondDer;       //!<!   result from shape derivatives for jet 2subjettiness_onepassca: 2nd derivative
  Double_t          fJetShape2subjettinessonepasscaFirstSub;        //!<!   result from shape derivatives for jet 2subjettiness_onepassca: 1st order subtracted
  Double_t          fJetShape2subjettinessonepasscaSecondSub;       //!<!   result from shape derivatives for jet 2subjettiness_onepassca: 2nd order subtracted

  Double_t          fJetShapeOpeningAngleonepasscaFirstDer;         //!<!   result from shape derivatives for jet OpeningAngle_onepassca: 1st derivative
  Double_t          fJetShapeOpeningAngleonepasscaSecondDer;        //!<!   result from shape derivatives for jet OpeningAngle_onepassca: 2nd derivative
  Double_t          fJetShapeOpeningAngleonepasscaFirstSub;         //!<!   result from shape derivatives for jet OpeningAngle_onepassca: 1st order subtracted
  Double_t          fJetShapeOpeningAngleonepasscaSecondSub;        //!<!   result from shape derivatives for jet OpeningAngle_onepassca: 2nd order subtracted

  Double_t          fSoftDropZg;                             //!<!   SoftDrop groomed momentum fraction
  Double_t          fSoftDropdR;                             //!<!   SoftDrop deltaR
  Double_t          fSoftDropPtfrac;                         //!<!   SoftDrop pt fraction after grooming
  Int_t             fSoftDropDropCount;                      //!<!   SoftDrop number of dropped branches [requires set_verbose_structure(bool enable=true)]
};

#endif
