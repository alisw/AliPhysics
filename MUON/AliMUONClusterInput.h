#ifndef ALIMUONCLUSTERINPUT_H
#define ALIMUONCLUSTERINPUT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id $*/

class TMinuit;
class AliMUONDigit;
class AliMUONRawCluster;
class AliSegmentation;
class AliMUONResponse;

#include <TClonesArray.h> // needed for inline function Digit


class AliMUONClusterInput : public TObject {
 public:
    static AliMUONClusterInput* Instance();
//  Configuration
    void SetDigits(Int_t chamber, TClonesArray* dig1, TClonesArray* dig2);
    void SetDigits(Int_t chamber, TClonesArray* dig);
    void SetCluster(AliMUONRawCluster* cluster);
// Access functions
    Int_t Chamber()  {return fChamber;}
    AliMUONDigit* Digit(Int_t cath, Int_t i) {return (AliMUONDigit*) (fDigits[cath]->UncheckedAt(i));}
    TClonesArray* Digits(Int_t cath) {return fDigits[cath];}
    Int_t NDigits(Int_t cath) {return fNDigits[cath];}
    AliSegmentation* Segmentation(Int_t cath)  {return fSegmentation[cath];}
    AliMUONResponse* Response()  {return fResponse;}    
// Fitting    
    TMinuit*      Fitter() {return fgMinuit;}
// Current cluster information    
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
    AliMUONClusterInput();
    AliMUONClusterInput(const AliMUONClusterInput& clusterInput){;}
    AliMUONClusterInput & operator = (const AliMUONClusterInput& rhs);
 private:
    static AliMUONClusterInput* fgClusterInput; // singleton instance
    // Digits
    TClonesArray*        fDigits[2];       // ! Array of pointers to digits
    Int_t                fNDigits[2];      // ! Number of digits
    AliSegmentation*     fSegmentation[2]; // ! Segmentation per cathode
    AliMUONResponse*     fResponse;        // ! Response
    Int_t                fNseg;            // ! number of cathode planes
    Int_t                fChamber;         // ! Current chamber number
    
    // Current cluster
    AliMUONRawCluster*   fCluster;         // ! current cluster
    Int_t                fNmul[2];         // ! current cluster multiplicity
    // Digits contribuing to current cluster
    Int_t                fix[500][2];      // ! List of x-positions for current cluster
    Int_t                fiy[500][2];      // ! List of y-positions for current cluster
    Float_t              fCharge[500][2];  // ! List of charges for current cluster
    Int_t                fChargeTot[2];    // ! Total charge
    Float_t              fQtot[2];         // ! Total charge
    Float_t              fZ;                // ! Current z-position
    // Fitter
    static TMinuit*      fgMinuit;          // ! Fitter
    ClassDef(AliMUONClusterInput, 0)        // Global data service for hit reconstruction
};
#endif

