#ifndef ALIMUONCONSTANTS_H
#define ALIMUONCONSTANTS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$*/
// Revision of includes 07/05/2004

#include <TObject.h>

class AliMUONConstants : public TObject {
 public:
    // return number of chambers
    static Int_t    NCh() {return fgNCh;}
    // return number of tracking chambers
    static Int_t    NTrackingCh() {return fgNTrackingCh;}
    // return number of trigger chambers
    static Int_t    NTriggerCh() {return fgNTriggerCh;}
    // return number of trigger circuits
    static Int_t    NTriggerCircuit() {return fgNTriggerCircuit;}
    // return position of chamber i
    static Float_t  DefaultChamberZ(Int_t i) {return fgDefaultChamberZ[i];}
    // return pointer to array of positions
    static Float_t* DefaultChamberZ() {return fgDefaultChamberZ;}
    // return chamber i inner diameter
    static Float_t  Dmin(Int_t i) {return fgDmin[i];}
    // return chamber i outer diameter
    static Float_t  Dmax(Int_t i) {return fgDmax[i];}
    // return maximum zoom for event display
    static Int_t    MaxZoom() {return fgMaxZoom;}
    // return half-distance between two half-chambers
    static Float_t    DzCh() {return fgDzCh;}
    // return half-distance between two slats
    static Float_t    DzSlat() {return fgDzSlat;}
    static  Int_t ChamberNumber(Float_t z); 
    // return SqrtKx3 and SqrtKy3 for Slat
    static Float_t SqrtKx3Slat() {return fgSqrtKx3Slat;}
    static Float_t SqrtKy3Slat() {return fgSqrtKy3Slat;}
    // return SqrtKx3 and SqrtKy3 for Station 1 & 2
    static Float_t SqrtKx3St12() {return fgSqrtKx3St12;}
    static Float_t SqrtKy3St12() {return fgSqrtKy3St12;}
       // return charge correlation (needed for response and for cluster finder !?)
    static Float_t ChargeCorrelSlat() {return fgChargeCorrelSlat;}
    static Float_t ChargeCorrelSt12() {return fgChargeCorrelSt12;}
     // return wire pitch
    static Float_t PitchSlat() {return fgPitchSlat;}
    static Float_t PitchSt12() {return fgPitchSt12;}

 protected:
    AliMUONConstants() : TObject() {}
    virtual ~AliMUONConstants(){}

 private:
    static Int_t  fgNCh;                //  Number of Chambers    
    static Int_t  fgNTrackingCh;        //  Number of Tracking Chambers
    static Int_t  fgNTriggerCh;         //  Number of Trigger Chambers
    static Int_t  fgNTriggerCircuit;    //  Number of Trigger Circuits
//
    static Float_t  fgDefaultChamberZ[14];    // ! Z-positions of chambers
    static Float_t  fgDmin[7];                // ! inner diameter
    static Float_t  fgDmax[7];                // ! outer diameter

    static Float_t  fgDzCh;             // half-distance between two half-chambers 
    static Float_t  fgDzSlat;           // half-distance between two slat on the same chamber
    static Float_t  fgSqrtKx3Slat;      // SqrtKx3 for Slat
    static Float_t  fgSqrtKy3Slat;      // SqrtKy3 for Slat
    static Float_t  fgSqrtKx3St12;      // SqrtKx3 for Station 1 & 2
    static Float_t  fgSqrtKy3St12;      // SqrtKy3 for Station 1 & 2
 
    static Float_t  fgChargeCorrelSlat;      // charge correlation for Slats
    static Float_t  fgChargeCorrelSt12;      // charge correlation for Station 1 & 2

    static Float_t  fgPitchSlat;      // wire pitch for Slats
    static Float_t  fgPitchSt12;      // wire pitch for Station 1 & 2

//
    static Int_t    fgMaxZoom;                // Maximum Zoom for event display
    ClassDef(AliMUONConstants, 0)             // MUON global constants 
};
	
#endif








