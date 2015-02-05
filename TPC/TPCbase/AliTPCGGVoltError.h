#ifndef ALI_TPCGG_VOLT_ERROR_H
#define ALI_TPCGG_VOLT_ERROR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \class AliTPCGGVoltError
/// \brief AliTPCGGVoltError class
///
/// The class calculates the electric field and the resulting space point distortions
/// due a Gating Grid (GG) voltage error. It uses the analytical solution for such a problem.
///
/// The input is the effective GG voltage residual in respect to the ideal setting. The effective
/// residual voltage can be set via the functions SetDeltaVGGx. Note that this effective
/// voltage-residuals are approx. a factor 0.9 lower than the actual difference in the setting
/// of the GG due to the fact that the voltage on the GG is partially screened by the wire
/// structure. The calculation has to be performed with the observable effective voltage difference.
///
/// Unfortunately, the class is not capable of calculation the $dz$ offset due to possible changes
/// of the drift velocity in dependence of the electric field. The special case of the numerical
/// approximation (AliTPCBoundaryVoltError), which is capable of calculating the same effect, should
/// be used for this purpose.
/// ![Picture from ROOT macro](AliTPCGGVoltError_h_359f7ec.png)
///
/// \author Jim Thomas, Stefan Rossegger, Magnus Mager
/// \date 27/04/2010


#include "AliTPCCorrection.h"

class AliTPCGGVoltError : public AliTPCCorrection {
public:
  AliTPCGGVoltError();
  virtual ~AliTPCGGVoltError();

  // initialization and update functions
  virtual void Init();
  virtual void Update(const TTimeStamp &timeStamp);

  // common setters and getters for ExB
  virtual void SetOmegaTauT1T2(Float_t omegaTau,Float_t t1,Float_t t2) {
    fT1=t1; fT2=t2;
    const Double_t wt0=t2*omegaTau;     fC0=1./(1.+wt0*wt0);
    const Double_t wt1=t1*omegaTau;     fC1=wt1/(1.+wt1*wt1);
  };

  void SetC0C1(Double_t c0,Double_t c1) {fC0=c0;fC1=c1;} // CAUTION: USE WITH CARE
  Float_t GetC0() const {return fC0;}
  Float_t GetC1() const {return fC1;}

  // setters and getters for GG
  void SetDeltaVGGA(Double_t deltaVGGA) {fDeltaVGGA=deltaVGGA;}
  void SetDeltaVGGC(Double_t deltaVGGC) {fDeltaVGGC=deltaVGGC;}
  Double_t GetDeltaVGGA() const {return fDeltaVGGA;}
  Double_t GetDeltaVGGC() const {return fDeltaVGGC;}

  void InitGGVoltErrorDistortion();

  Float_t GetIntErOverEz(const Float_t x[],const Short_t roc);

  virtual void Print(const Option_t* option="") const;

protected:
  virtual void GetCorrection(const Float_t x[],const Short_t roc, Float_t dx[]);
private:

  Float_t fC0; ///< coefficient C0                 (compare Jim Thomas's notes for definitions)
  Float_t fC1; ///< coefficient C1                 (compare Jim Thomas's notes for definitions)

  Double_t fDeltaVGGA;            ///< Missmatch of gating grid voltage on A-side [V]
  Double_t fDeltaVGGC;            ///< Missmatch of gating grid voltage on C-side [V]
  Double_t fGGVoltErrorER[kNZ][kNR]; ///< Array to store electric field for GGVoltError calculation

  Bool_t fInitLookUp;             ///< flag to check it the Look Up table was created

  /// \cond CLASSIMP
  ClassDef(AliTPCGGVoltError,1);
  /// \endcond
};

#endif
