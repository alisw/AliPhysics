#ifndef ALI_TPC_EX_B_TWIST_H
#define ALI_TPC_EX_B_TWIST_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \class AliTPCExBTwist
/// \brief AliTPCExBTwist class
///
/// The class calculates the space point distortions due to a mismatch
/// of the E and B field axis (original code from STAR)
///
/// It assumes a homogeneous magnetic and electric field. The calculations performed
/// within this class can be written in a compact form as well (compare notes of Jim Thomas)
///
/// The inputs to this class are the twist angles in the x and y plane ($\theta_x,\theta_y$)
/// in cartesian coordinates. They can be set via the member functions SetXTwist and SetYTwist.
///
/// ![Picture from ROOT macro](AliTPCExBTwist_h_76e03fd.png)
///
/// Support since 2010: mrain.ivanov@cern.ch
///
/// \author Jim Thomas, Magnus Mager, Stefan Rossegger
/// \date 27/04/2010

#include "AliTPCCorrection.h"

class AliTPCExBTwist : public AliTPCCorrection {
public:
  AliTPCExBTwist();
  virtual ~AliTPCExBTwist();
  virtual Bool_t AddCorrectionCompact(AliTPCCorrection* corr, Double_t weight);
  // initialization and update functions
  virtual void Init();
  virtual void Update(const TTimeStamp &timeStamp);


  // common setters and getters for ExB
  virtual void SetOmegaTauT1T2(Float_t omegaTau,Float_t t1,Float_t t2) {
    fT1=t1; fT2=t2;
    const Float_t wt1=t1*omegaTau;     fC1=wt1/(1.+wt1*wt1);
    const Float_t wt2=t2*omegaTau;     fC2=wt2*wt2/(1.+wt2*wt2);
  };
  void SetC1C2(Float_t c1,Float_t c2) {fC1=c1;fC2=c2;} // CAUTION: USE WITH CARE
  Float_t GetC1() const {return fC1;}
  Float_t GetC2() const {return fC2;}

  // setters and getters for twist
  void SetXTwist(Float_t xTwist) {fXTwist=xTwist;}
  void SetYTwist(Float_t yTwist) {fYTwist=yTwist;}
  Float_t GetXTwist() const {return fXTwist;}
  Float_t GetYTwist() const {return fYTwist;}

  virtual void Print(const Option_t* option="") const;

protected:
  virtual void GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]);

private:
  Float_t fC1; ///< coefficient C1                 (compare Jim Thomas's notes for definitions)
  Float_t fC2; ///< coefficient C2                 (compare Jim Thomas's notes for definitions)

  Float_t fXTwist;               ///< Twist of E to B field in X-Z [rad]
  Float_t fYTwist;               ///< Twist of E to B field in Y-Z [rad]

  /// \cond CLASSIMP
  ClassDef(AliTPCExBTwist,1);
  /// \endcond
};

#endif
