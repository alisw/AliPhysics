#ifndef ALIMUONCONSTANTS_H
#define ALIMUONCONSTANTS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id $*/

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
    // return poistion of chamber i
    static Float_t  DefaultChamberZ(Int_t i) {return fgDefaultChamberZ[i];}
    // return pointer to array of positions
    static Float_t* DefaultChamberZ() {return fgDefaultChamberZ;}
    // return chamber i inner diameter
    static Float_t  Dmin(Int_t i) {return fgDmin[i];}
    // return chamber i outer diameter
    static Float_t  Dmax(Int_t i) {return fgDmax[i];}
    // return maximum zoom for event display
    static Int_t    MaxZoom() {return fgMaxZoom;}
 private:
    AliMUONConstants(){}
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
//
    static Int_t    fgMaxZoom;                // Maximum Zoom for event display
    ClassDef(AliMUONConstants, 1)             // MUON global constants 
};
	
#endif








