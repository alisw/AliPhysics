#ifndef ALITRDPARTID_H
#define ALITRDPARTID_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>

class AliESDtrack;
class TProfile;
class TF1;


class AliTRDPartID: public TObject {
  public: 
    AliTRDPartID();
    AliTRDPartID(TF1* betheBloch, Double_t res, Double_t range);
    virtual ~AliTRDPartID();

    Bool_t          MakePID(AliESDtrack* track);

    void            FitBetheBloch(TProfile* dEdxVsBetaGamma);
    inline TF1*     GetBetheBloch() {return fBetheBloch;};
    TF1*            CreateBetheBloch(Double_t mass);

  private:
    static Double_t fcnBetheBloch(Double_t* xx, Double_t* par);
    static Double_t fcnBetheBlochMass(Double_t* xx, Double_t* par);

    TF1*            fBetheBloch;   // parametrized bethe bloch function
    Double_t        fRes;          // relative dE/dx resolution
    Double_t        fRange;        // cut off in standard deviations

    ClassDef(AliTRDPartID,1)   // TRD PID class
};

#endif


