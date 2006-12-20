#ifndef ALIMUONCONSTANTS_H
#define ALIMUONCONSTANTS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$*/
// Revision of includes 07/05/2004
//
/// \ingroup base
/// \class AliMUONConstants
/// \brief MUON global constants

#include <TObject.h>

class AliMUONConstants : public TObject {
 public:
    /// Return number of chambers
    static Int_t    NCh() {return fgNCh;}
    /// Return number of tracking chambers
    static Int_t    NTrackingCh() {return fgNTrackingCh;}
    /// Return number of tracking stations
    static Int_t    NTrackingSt() {return fgNTrackingSt;}
    /// Return number of trigger chambers
    static Int_t    NTriggerCh() {return fgNTriggerCh;}
    /// Return number of trigger circuits
    static Int_t    NTriggerCircuit() {return fgNTriggerCircuit;}
    /// Return number of detection element
    static Int_t    NDetElem() {return fgNDetElem;}
    /// Return number of geometry modules
    static Int_t    NGeomModules() {return fgNGeomModules;}
    /// Return position of chamber i
    static Float_t  DefaultChamberZ(Int_t i) {return fgDefaultChamberZ[i];}
    /// Return ratio between trigger chambers
    static Float_t  DefaultRatioTriggerChamber(Int_t i) {return fgDefaultRatioTriggerChamber[i];}
    /// Return Inclination with respect the vertical axis of stations 345
    static Float_t  St345Inclination()       {return fgSt345inclination;}
    /// Return pointer to array of positions
    static Float_t* DefaultChamberZ()        {return fgDefaultChamberZ;}
    /// Return chamber i inner diameter
    static Float_t  Dmin(Int_t i)            {return fgDmin[i];}
    /// Return chamber i inner radius
    static Float_t  Rmin(Int_t i)            {return Dmin(i)/2.0;}
    /// Return chamber i outer diameter
    static Float_t  Dmax(Int_t i)            {return fgDmax[i];}
    /// Return chamber i outer radius
    static Float_t  Rmax(Int_t i)            {return Dmax(i)/2.0;}
    /// Return maximum zoom for event display
    static Int_t    MaxZoom()                {return fgMaxZoom;}
    /// Return half-distance between two half-chambers
    static Float_t  DzCh()                   {return fgDzCh;}
    /// Return half-distance between two slats
    static Float_t  DzSlat()                 {return fgDzSlat;}
    /// Return chamber number according z position of hit.
    static  Int_t ChamberNumber(Float_t z); 
    /// Return SqrtKx3 for Slat
    static Float_t SqrtKx3()                 {return fgSqrtKx3;}
    /// Return SqrtKy3 for Slat
    static Float_t SqrtKy3()                 {return fgSqrtKy3;}
    /// Return SqrtKx3 for Station 1 & 2
    static Float_t SqrtKx3St1()              {return fgSqrtKx3St1;}
    /// Return SqrtKy3 for Station 1 & 2
    static Float_t SqrtKy3St1()              {return fgSqrtKy3St1;}
    /// Return charge correlation (needed for response and for cluster finder !?)
    static Float_t ChargeCorrel()            {return fgChargeCorrel;}
    /// Return charge correlation for Station 1 & 2 (needed for response and for cluster finder !?)
    static Float_t ChargeCorrelSt1()         {return fgChargeCorrelSt1;}
    /// Return wire pitch
    static Float_t Pitch()    {return fgPitch;}
    /// Return wire pitch for Station 1 & 2
    static Float_t PitchSt1() {return fgPitchSt1;}
    /// Return coil z-position
    static Double_t CoilZ() {return fgCoilZ;}
    /// Return coil lenght
    static Double_t CoilL() {return fgCoilL;}
    /// Return yoke z-position
    static Double_t YokeZ() {return fgYokeZ;}
    /// Return yoke lenght
    static Double_t YokeL() {return fgYokeL;}
    /// Return the z position of the end of the absorber
    static Double_t ZAbsorberEnd() {return fgZAbsorberEnd;}
    /// Return the number of elements making the absorber
    static Int_t    NAbsorberElements() {return fgNAbsorberElements;}
    /// Return the z position of the element "iElement" making the absorber
    static Double_t ZAbsorberElement(Int_t iElement) {return ((iElement >= 0) && (iElement < fgNAbsorberElements)) ? fgZAbsorberElement[iElement] : 0.;}
    /// Return the radiation lenght of the element "iElement" making the inner part of the absorber
    static Double_t X0AbsorberIn(Int_t iElement) {return ((iElement >= 0) && (iElement < fgNAbsorberElements)) ? fgX0AbsorberIn[iElement] : -1.;}
    /// Return the radiation lenght of the element "iElement" making the outer part of the absorber
    static Double_t X0AbsorberOut(Int_t iElement) {return ((iElement >= 0) && (iElement < fgNAbsorberElements)) ? fgX0AbsorberOut[iElement] : -1.;}
    /// Return chamber thickness in X0
    static Double_t ChamberThicknessInX0() {return fgChamberThicknessInX0;}
    /// Return Trigger ToF Limit (75 ns)
    static Float_t TriggerTofLimit() {return fgkTriggerTofLimit;}
 
 protected:
    /// Default constructor
    AliMUONConstants() : TObject() {}
    /// Destructor
    virtual ~AliMUONConstants(){}

 private:
    static Int_t  fgNCh;                ///<  Number of Chambers    
    static Int_t  fgNTrackingCh;        ///<  Number of Tracking Chambers
    static Int_t  fgNTrackingSt;        ///<  Number of Tracking Stations
    static Int_t  fgNTriggerCh;         ///<  Number of Trigger Chambers
    static Int_t  fgNTriggerCircuit;    ///<  Number of Trigger Circuits
    static Int_t  fgNDetElem;           ///<  Number of Detection Elements.
    static Int_t  fgNGeomModules;       ///<  Number of Geometry modules   

    static Float_t  fgDefaultChamberZ[14];		//!< Z-positions of chambers
    static Float_t  fgDefaultRatioTriggerChamber[4];	///< Ratio between trigger chambers
    static Float_t  fgSt345inclination;			//!< Inclination with respect the vertical axis of stations 345
    static Float_t  fgDmin[7];				//!< Inner diameter
    static Float_t  fgDmax[7];				//!< Outer diameter

    static Float_t  fgDzCh;             ///< Half-distance between two half-chambers 
    static Float_t  fgDzSlat;           ///< Half-distance between two slat on the same chamber
    static Float_t  fgSqrtKx3;          ///< SqrtKx3 for St2 & Slat
    static Float_t  fgSqrtKy3;          ///< SqrtKy3 for St2 & Slat
    static Float_t  fgSqrtKx3St1;       ///< SqrtKx3 for Station 1 
    static Float_t  fgSqrtKy3St1;       ///< SqrtKy3 for Station 1 
 
    static Float_t  fgChargeCorrel;      ///< Charge correlation for St2 & Slats
    static Float_t  fgChargeCorrelSt1;   ///< Charge correlation for Station 1

    static Float_t  fgPitch;             ///< Wire pitch for St2 & Slats
    static Float_t  fgPitchSt1;          ///< Wire pitch for Station 1

    static Double_t  fgChamberThicknessInX0; ///< default chamber thickness in X0 for reconstruction
    
    static Double_t fgCoilZ; ///< Coil z-position
    static Double_t fgCoilL; ///< Coil lenght
    static Double_t fgYokeZ; ///< Yoke z-position
    static Double_t fgYokeL; ///< Yoke lenght

    static Double_t fgZAbsorberEnd;		///< Z position of the end of the absorber
    static Int_t    fgNAbsorberElements;	///< Number of elements making the absorber
    static Double_t fgZAbsorberElement[5];	///< Z positions of elements making the absorber
    static Double_t fgX0AbsorberIn[5];		///< Radiation lenght of materials in the inner part of the absorber (theta > 3 degrees)
    static Double_t fgX0AbsorberOut[5];		///< Radiation lenght of materials in the outer part of the absorber (theta > 3 degrees)
    
    static Int_t    fgMaxZoom;          ///< Maximum Zoom for event display
    static Float_t  fgkTriggerTofLimit; ///< Particle above this threshold are discarded in trigger algorithm
    
    ClassDef(AliMUONConstants, 0) // MUON global constants 
};
	
#endif








