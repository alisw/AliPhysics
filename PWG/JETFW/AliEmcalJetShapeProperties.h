#ifndef ALIEMCALJETSHAPEPROPERTIES_H
#define ALIEMCALJETSHAPEPROPERTIES_H

#include <TArrayF.h>
#include <Rtypes.h>
#include <TString.h>

class AliEmcalJetShapeProperties {

public:

  AliEmcalJetShapeProperties();
  AliEmcalJetShapeProperties(const AliEmcalJetShapeProperties &jetshape);
  AliEmcalJetShapeProperties& operator=(const AliEmcalJetShapeProperties &jetshape);

  //jet shape derivatives
  //jet mass
  void              SetFirstDerivative(Double_t d)                  { fJetShapeMassFirstDer = d                       ; }
  void              SetSecondDerivative(Double_t d)                 { fJetShapeMassSecondDer = d                      ; }
  void              SetFirstOrderSubtracted(Double_t d)             { fJetShapeMassFirstSub = d                       ; }
  void              SetSecondOrderSubtracted(Double_t d)            { fJetShapeMassSecondSub = d                      ; }
  Double_t          GetFirstDerivative()                      const { return fJetShapeMassFirstDer                    ; }
  Double_t          GetSecondDerivative()                     const { return fJetShapeMassSecondDer                   ; }
  Double_t          GetFirstOrderSubtracted()                 const { return fJetShapeMassFirstSub                    ; }
  Double_t          GetSecondOrderSubtracted()                const { return fJetShapeMassSecondSub                   ; }

  //jet structure function
  TArrayF           GetGRNumerator()                          const { return fGRNumerator                             ; }
  TArrayF           GetGRDenominator()                        const { return fGRDenominator                           ; }
  TArrayF           GetGRNumeratorSub()                       const { return fGRNumeratorSub                          ; }
  TArrayF           GetGRDenominatorSub()                     const { return fGRDenominatorSub                        ; }
  void              AddGRNumAt(Float_t num, Int_t idx)              { fGRNumerator.AddAt(num, idx)                    ; }
  void              AddGRDenAt(Float_t den, Int_t idx)              { fGRDenominator.AddAt(den, idx)                  ; }
  void              SetGRNumSize(UInt_t s)                          { fGRNumerator.Set(s)                             ; }
  void              SetGRDenSize(UInt_t s)                          { fGRDenominator.Set(s)                           ; }

  void              AddGRNumSubAt(Float_t num, Int_t idx)           { fGRNumeratorSub.AddAt(num, idx)                 ; }
  void              AddGRDenSubAt(Float_t den, Int_t idx)           { fGRDenominatorSub.AddAt(den, idx)               ; }
  void              SetGRNumSubSize(UInt_t s)                       { fGRNumeratorSub.Set(s)                          ; }
  void              SetGRDenSubSize(UInt_t s)                       { fGRDenominatorSub.Set(s)                        ; }
  void              PrintGR();

  //Angularity
  void              SetFirstDerivativeAngularity(Double_t d)        { fJetShapeAngularityFirstDer = d                 ; }
  void              SetSecondDerivativeAngularity(Double_t d)       { fJetShapeAngularitySecondDer = d                ; }
  void              SetFirstOrderSubtractedAngularity(Double_t d)   { fJetShapeAngularityFirstSub = d                 ; }
  void              SetSecondOrderSubtractedAngularity(Double_t d)  { fJetShapeAngularitySecondSub = d                ; }
  Double_t          GetFirstDerivativeAngularity()            const { return fJetShapeAngularityFirstDer              ; }
  Double_t          GetSecondDerivativeAngularity()           const { return fJetShapeAngularitySecondDer             ; }
  Double_t          GetFirstOrderSubtractedAngularity()       const { return fJetShapeAngularityFirstSub              ; }
  Double_t          GetSecondOrderSubtractedAngularity()      const { return fJetShapeAngularitySecondSub             ; }

  //pTD
  void              SetFirstDerivativepTD(Double_t d)               { fJetShapepTDFirstDer = d                        ; }
  void              SetSecondDerivativepTD(Double_t d)              { fJetShapepTDSecondDer = d                       ; }
  void              SetFirstOrderSubtractedpTD(Double_t d)          { fJetShapepTDFirstSub = d                        ; }
  void              SetSecondOrderSubtractedpTD(Double_t d)         { fJetShapepTDSecondSub = d                       ; }
  Double_t          GetFirstDerivativepTD()                   const { return fJetShapepTDFirstDer                     ; }
  Double_t          GetSecondDerivativepTD()                  const { return fJetShapepTDSecondDer                    ; }
  Double_t          GetFirstOrderSubtractedpTD()              const { return fJetShapepTDFirstSub                     ; }
  Double_t          GetSecondOrderSubtractedpTD()             const { return fJetShapepTDSecondSub                    ; }

  //Circularity
  void              SetFirstDerivativeCircularity(Double_t d)       { fJetShapeCircularityFirstDer = d                ; }
  void              SetSecondDerivativeCircularity(Double_t d)      { fJetShapeCircularitySecondDer = d               ; }
  void              SetFirstOrderSubtractedCircularity(Double_t d)  { fJetShapeCircularityFirstSub = d                ; }
  void              SetSecondOrderSubtractedCircularity(Double_t d) { fJetShapeCircularitySecondSub = d               ; }
  Double_t          GetFirstDerivativeCircularity()           const { return fJetShapeCircularityFirstDer             ; }
  Double_t          GetSecondDerivativeCircularity()          const { return fJetShapeCircularitySecondDer            ; }
  Double_t          GetFirstOrderSubtractedCircularity()      const { return fJetShapeCircularityFirstSub             ; }
  Double_t          GetSecondOrderSubtractedCircularity()     const { return fJetShapeCircularitySecondSub            ; }

  //Sigma2
  void              SetFirstDerivativeSigma2(Double_t d)            { fJetShapeSigma2FirstDer = d                     ; }
  void              SetSecondDerivativeSigma2(Double_t d)           { fJetShapeSigma2SecondDer = d                    ; }
  void              SetFirstOrderSubtractedSigma2(Double_t d)       { fJetShapeSigma2FirstSub = d                     ; }
  void              SetSecondOrderSubtractedSigma2(Double_t d)      { fJetShapeSigma2SecondSub = d                    ; }
  Double_t          GetFirstDerivativeSigma2()           const      { return fJetShapeSigma2FirstDer                  ; }
  Double_t          GetSecondDerivativeSigma2()          const      { return fJetShapeSigma2SecondDer                 ; }
  Double_t          GetFirstOrderSubtractedSigma2()      const      { return fJetShapeSigma2FirstSub                  ; }
  Double_t          GetSecondOrderSubtractedSigma2()     const      { return fJetShapeSigma2SecondSub                 ; }


  //number of contituents
  void              SetFirstDerivativeConstituent(Double_t d)       { fJetShapeConstituentFirstDer = d                ; }
  void              SetSecondDerivativeConstituent(Double_t d)      { fJetShapeConstituentSecondDer = d               ; }
  void              SetFirstOrderSubtractedConstituent(Double_t d)  { fJetShapeConstituentFirstSub = d                ; }
  void              SetSecondOrderSubtractedConstituent(Double_t d) { fJetShapeConstituentSecondSub = d               ; }
  Double_t          GetFirstDerivativeConstituent()           const { return fJetShapeConstituentFirstDer             ; }
  Double_t          GetSecondDerivativeConstituent()          const { return fJetShapeConstituentSecondDer            ; }
  Double_t          GetFirstOrderSubtractedConstituent()      const { return fJetShapeConstituentFirstSub             ; }
  Double_t          GetSecondOrderSubtractedConstituent()     const { return fJetShapeConstituentSecondSub            ; }

  //leading minus subleading constituent
  void              SetFirstDerivativeLeSub(Double_t d)             { fJetShapeLeSubFirstDer = d                      ; }
  void              SetSecondDerivativeLeSub(Double_t d)            { fJetShapeLeSubSecondDer = d                     ; }
  void              SetFirstOrderSubtractedLeSub(Double_t d)        { fJetShapeLeSubFirstSub = d                      ; }
  void              SetSecondOrderSubtractedLeSub(Double_t d)       { fJetShapeLeSubSecondSub = d                     ; }
  Double_t          GetFirstDerivativeLeSub()                 const { return fJetShapeLeSubFirstDer                   ; }
  Double_t          GetSecondDerivativeLeSub()                const { return fJetShapeLeSubSecondDer                  ; }
  Double_t          GetFirstOrderSubtractedLeSub()            const { return fJetShapeLeSubFirstSub                   ; }
  Double_t          GetSecondOrderSubtractedLeSub()           const { return fJetShapeLeSubSecondSub                  ; }

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

  void              PrintGR() const;

protected:

  Double_t          fJetShapeMassFirstDer;         //!   result from shape derivatives for jet mass: 1st derivative
  Double_t          fJetShapeMassSecondDer;        //!   result from shape derivatives for jet mass: 2nd derivative
  Double_t          fJetShapeMassFirstSub;         //!   result from shape derivatives for jet mass: 1st order subtracted
  Double_t          fJetShapeMassSecondSub;        //!   result from shape derivatives for jet mass: 2nd order subtracted

  TArrayF           fGRNumerator;                  //!   array with angular structure function numerator
  TArrayF           fGRDenominator;                //!   array with angular structure function denominator
  TArrayF           fGRNumeratorSub;               //!   array with angular structure function numerator
  TArrayF           fGRDenominatorSub;             //!   array with angular structure function denominator

  Double_t          fJetShapeAngularityFirstDer;   //!   result from shape derivatives for jet Angularity: 1st derivative
  Double_t          fJetShapeAngularitySecondDer;  //!   result from shape derivatives for jet Angularity: 2nd derivative
  Double_t          fJetShapeAngularityFirstSub;   //!   result from shape derivatives for jet Angularity: 1st order subtracted
  Double_t          fJetShapeAngularitySecondSub;  //!   result from shape derivatives for jet Angularity: 2nd order subtracted

  Double_t          fJetShapepTDFirstDer;          //!   result from shape derivatives for jet pTD: 1st derivative
  Double_t          fJetShapepTDSecondDer;         //!   result from shape derivatives for jet pTD: 2nd derivative
  Double_t          fJetShapepTDFirstSub;          //!   result from shape derivatives for jet pTD: 1st order subtracted
  Double_t          fJetShapepTDSecondSub;         //!   result from shape derivatives for jet pTD: 2nd order subtracted

  Double_t          fJetShapeCircularityFirstDer;  //!   result from shape derivatives for jet circularity: 1st derivative
  Double_t          fJetShapeCircularitySecondDer; //!   result from shape derivatives for jet circularity: 2nd derivative
  Double_t          fJetShapeCircularityFirstSub;  //!   result from shape derivatives for jet circularity: 1st order subtracted
  Double_t          fJetShapeCircularitySecondSub; //!   result from shape derivatives for jetcircularity: 2nd order subtracted

  Double_t          fJetShapeSigma2FirstDer;       //!   result from shape derivatives for jet sigma2: 1st derivative
  Double_t          fJetShapeSigma2SecondDer;      //!   result from shape derivatives for jet sigma2: 2nd derivative
  Double_t          fJetShapeSigma2FirstSub;       //!   result from shape derivatives for jet sigma2: 1st order subtracted
  Double_t          fJetShapeSigma2SecondSub;      //!   result from shape derivatives for jetsigma2: 2nd order subtracted

  Double_t          fJetShapeConstituentFirstDer;  //!   result from shape derivatives for jet const: 1st derivative
  Double_t          fJetShapeConstituentSecondDer; //!   result from shape derivatives for jet const: 2nd derivative
  Double_t          fJetShapeConstituentFirstSub;  //!   result from shape derivatives for jet const: 1st order subtracted
  Double_t          fJetShapeConstituentSecondSub; //!   result from shape derivatives for jet const: 2nd order subtracted

  Double_t          fJetShapeLeSubFirstDer;        //!   result from shape derivatives for jet LeSub: 1st derivative
  Double_t          fJetShapeLeSubSecondDer;       //!   result from shape derivatives for jet LeSub: 2nd derivative
  Double_t          fJetShapeLeSubFirstSub;        //!   result from shape derivatives for jet LeSub: 1st order subtracted
  Double_t          fJetShapeLeSubSecondSub;       //!   result from shape derivatives for jet LeSub: 2nd order subtracted

  Double_t          fJetShape1subjettinessktFirstDer;        //!   result from shape derivatives for jet 1subjettiness_kt: 1st derivative
  Double_t          fJetShape1subjettinessktSecondDer;       //!   result from shape derivatives for jet 1subjettiness_kt: 2nd derivative
  Double_t          fJetShape1subjettinessktFirstSub;        //!   result from shape derivatives for jet 1subjettiness_kt: 1st order subtracted
  Double_t          fJetShape1subjettinessktSecondSub;       //!   result from shape derivatives for jet 1subjettiness_kt: 2nd order subtracted
  
  Double_t          fJetShape2subjettinessktFirstDer;        //!   result from shape derivatives for jet 2subjettiness_kt: 1st derivative
  Double_t          fJetShape2subjettinessktSecondDer;       //!   result from shape derivatives for jet 2subjettiness_kt: 2nd derivative
  Double_t          fJetShape2subjettinessktFirstSub;        //!   result from shape derivatives for jet 2subjettiness_kt: 1st order subtracted
  Double_t          fJetShape2subjettinessktSecondSub;       //!   result from shape derivatives for jet 2subjettiness_kt: 2nd order subtracted

  Double_t          fJetShape3subjettinessktFirstDer;        //!   result from shape derivatives for jet 3subjettiness_kt: 1st derivative
  Double_t          fJetShape3subjettinessktSecondDer;       //!   result from shape derivatives for jet 3subjettiness_kt: 2nd derivative
  Double_t          fJetShape3subjettinessktFirstSub;        //!   result from shape derivatives for jet 3subjettiness_kt: 1st order subtracted
  Double_t          fJetShape3subjettinessktSecondSub;       //!   result from shape derivatives for jet 3subjettiness_kt: 2nd order subtracted

  Double_t          fJetShapeOpeningAnglektFirstDer;         //!   result from shape derivatives for jet OpeningAngle_kt: 1st derivative
  Double_t          fJetShapeOpeningAnglektSecondDer;        //!   result from shape derivatives for jet OpeningAngle_kt: 2nd derivative
  Double_t          fJetShapeOpeningAnglektFirstSub;         //!   result from shape derivatives for jet OpeningAngle_kt: 1st order subtracted
  Double_t          fJetShapeOpeningAnglektSecondSub;        //!   result from shape derivatives for jet OpeningAngle_kt: 2nd order subtracted
};

#endif
