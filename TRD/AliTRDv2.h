#ifndef TRDv2_H
#define TRDv2_H
////////////////////////////////////////////////////////
//  Manager and hits classes for set:TRD version 2    //
////////////////////////////////////////////////////////
 
#include "AliTRD.h"

class AliTRDv2 : public AliTRD {

public:
  AliTRDv2() {}
  AliTRDv2(const char *name, const char *title);
  virtual      ~AliTRDv2() {}
  virtual void  CreateGeometry();
  virtual void  CreateMaterials();
  virtual Int_t IsVersion() const {return 2;}
  virtual void  StepManager();
  virtual void  Init();
  virtual void  DrawDetector();

protected:
  Int_t        fIdSenO1;    // Sensitive volume identifier for outer chambers
  Int_t        fIdSenO2;    // Sensitive volume identifier for outer chambers
  Int_t        fIdSenO3;    // Sensitive volume identifier for outer chambers
  Int_t        fIdSenO4;    // Sensitive volume identifier for outer chambers
  Int_t        fIdSenO5;    // Sensitive volume identifier for outer chambers
  Int_t        fIdSenO6;    // Sensitive volume identifier for outer chambers
  Int_t        fIdSenI1;    // Sensitive volume identifier for inner chambers
  Int_t        fIdSenI2;    // Sensitive volume identifier for inner chambers
  Int_t        fIdSenI3;    // Sensitive volume identifier for inner chambers
  Int_t        fIdSenI4;    // Sensitive volume identifier for inner chambers
  Int_t        fIdSenI5;    // Sensitive volume identifier for inner chambers
  Int_t        fIdSenI6;    // Sensitive volume identifier for inner chambers
            
private:
  // Inline functions for AliTRDv2
  
  inline Float_t Eloss(Float_t rndm)
    {
      //
      // Calculates the energy loss 
      // 1/E^2. distribution for the fluctuations
      //
      // Exponent of the distribution for the energy loss 
      // 2.0 is the apropriate value for Argon, 2.2 would be for Neon,
      // and how about Xenon? We take the Argon value for the time being.
      const Float_t kEexp = 2.0;
      // First ionization potential for the gas mixture (90% Xe + 10% CO2)
      // taken from: Ionization Measurements in High Energy Physics, Springer
      const Float_t kPoti = 12.3E-9;
      // Maximum energy (10 keV);
      const Float_t kEend = 10.0E-6;
      
      Float_t ex   = 1. - kEexp;
      Float_t xpot = TMath::Power(kPoti,ex);
      Float_t xend = TMath::Power(kEend,ex);
      Float_t elos = (1. - rndm)*xpot + rndm*xend;
      
      return(TMath::Power(elos,(1./ex)) - kPoti);
      
    }
  
  inline Float_t BetheBloch(Float_t xx) 
    {
      //
      // Parametrization of the Bethe-Bloch-curve
      // The parametrization is the same as for the TPC and is taken from
      // Lehrhaus.
      // The parameters have been adjusted to Xe-data found in:
      // Allison & Cobb, Ann. Rev. Nucl. Sci. (1980), 30, 253
      //
      const Float_t kP1 = 0.76176E-1;
      const Float_t kP2 = 10.632;
      const Float_t kP3 = 3.17983E-6;
      const Float_t kP4 = 1.8631;
      const Float_t kP5 = 1.9479;
      
      Float_t yy = xx / TMath::Sqrt(1. + xx*xx);
      Float_t aa = TMath::Power(yy,kP4);
      Float_t bb = TMath::Power((1./xx),kP5);
      bb = TMath::Log(kP3 + bb);
      
      return((kP2 - aa - bb)*kP1 / aa);
      
    }
  
  ClassDef(AliTRDv2,1)     // Transition Radiation Detector version 2
};

#endif
