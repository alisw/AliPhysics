#ifndef ALIMUONCLUSTERINPUT_H
#define ALIMUONCLUSTERINPUT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id $*/

#include <TObject.h>
#include <TClonesArray.h>

class TMinuit;
class AliMUONDigit;
class AliMUONRawCluster;
class AliMUONSegmentation;
class AliMUONResponse;



class AliMUONClusterInput : public TObject {
 public:
    static AliMUONClusterInput* Instance();
//  Setters
    void SetDigits(Int_t chamber, TClonesArray* dig1, TClonesArray* dig2);
    void SetDigits(Int_t chamber, TClonesArray* dig);
    void SetCluster(AliMUONRawCluster* cluster);
// Access functions
    AliMUONDigit* Digit(Int_t cath, Int_t i) {return (AliMUONDigit*) (fDigits[cath]->UncheckedAt(i));}
    TMinuit*      Fitter() {return fgMinuit;}
    Float_t       TotalCharge(Int_t cath) {return fChargeTot[cath];}
    Float_t       Charge(Int_t dig, Int_t cath) {return fCharge[dig][cath];}
    Int_t         Ix(Int_t dig, Int_t cath) {return fix[dig][cath];}
    Int_t         Iy(Int_t dig, Int_t cath) {return fiy[dig][cath];}
    Int_t         Nmul(Int_t cath)  {return fNmul[cath];}
//  Helpers for Fit     
    Float_t DiscrChargeS1(Int_t i,Double_t *par);
    Float_t DiscrChargeCombiS1(Int_t i,Double_t *par, Int_t cath);
    Float_t DiscrChargeS2(Int_t i,Double_t *par);
    Float_t DiscrChargeCombiS2(Int_t i,Double_t *par, Int_t cath);
// 
 protected:
    AliMUONClusterInput(){;}
 private:
    static AliMUONClusterInput* fgClusterInput;
    // Digits
    TClonesArray*        fDigits[2];       // ! Array of pointers to digits
    AliMUONSegmentation* fSegmentation[2]; // ! Segmentation per cathode
    AliMUONResponse*     fResponse;        // ! Response
    Int_t                fNseg;            // ! number of cathode planes
    // Current cluster
    AliMUONRawCluster*   fCluster;         // ! current cluster
    Int_t                fNmul[2];         // ! current cluster multiplicity
    // Digits contribuing to current cluster
    Int_t                fix[500][2];      // ! List of x-positions for current cluster
    Int_t                fiy[500][2];      // ! List of y-positions for current cluster
    Float_t              fCharge[500][2];  // ! List of charges for current cluster
    Int_t                fChargeTot[2];    // ! Total charge
    Float_t              fQtot[2];         // ! Total charge
    // Fitter
    static TMinuit*      fgMinuit;          // ! Fitter
    ClassDef(AliMUONClusterInput, 1) // Class definition in ROOT context
};
#endif

