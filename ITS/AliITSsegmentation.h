#ifndef ALIITSSEGMENTATION_H
#define ALIITSSEGMENTATION_H

#include <TObject.h>

class TF1;
class AliITSgeom;
//----------------------------------------------
//
// ITS  segmentation virtual base class
//
class AliITSsegmentation :
public TObject {
 public:
    virtual ~AliITSsegmentation() {}
    // Set Detector Segmentation Parameters
    //
    // Detector size  
    virtual void    SetDetSize(Float_t,Float_t,Float_t) {}
    // Cell size   
    virtual void    SetPadSize(Float_t,Float_t) {}
    // Maximum number of cells along the two coordinates  
    virtual void    SetNPads(Int_t,Int_t) {}
    // Returns the maximum number of cells (digits) posible
    virtual Int_t   GetNPads(){return 0;}
    // Set angles - find a generic name fit for other detectors as well
    // might be useful for beam test setups (3 angles ?)
    virtual void    SetAngles(Float_t, Float_t) {}
    // Set layer
    virtual void SetLayer(Int_t) {}
    // Transform from real to cell coordinates
    virtual void    GetPadIxz(Float_t,Float_t,Int_t &,Int_t &) {}
    // Transform from cell to real coordinates
    virtual void    GetPadCxz(Int_t,Int_t,Float_t &,Float_t &) {}
    // Transform from real global to local coordinates
    virtual void    GetLocal(Int_t,Float_t *,Float_t *) {}
    // Transform from real local to global coordinates
    virtual void    GetGlobal(Int_t,Float_t *,Float_t *) {}
    // Local transformation of real local coordinates -
    virtual void    GetPadTxz(Float_t &,Float_t &) {}
    // Transformation from Geant cm detector center local coordinates
    // to detector segmentation/cell coordiantes starting from (0,0).
    virtual void    LocalToDet(Float_t,Float_t,Int_t &,Int_t &){}
    // Transformation from detector segmentation/cell coordiantes starting
    // from (0,0) to Geant cm detector center local coordinates.
    virtual void    DetToLocal(Int_t,Int_t,Float_t &,Float_t &){}
    // Initialisation
    virtual void Init() {}
    //
    // Get member data
    //
    // Detector type geometry
    virtual AliITSgeom* Geometry() {return 0;}
    // Detector length
    virtual Float_t Dx() {return 0.;}
    // Detector width
    virtual Float_t Dz() {return 0.;}
    // Detector thickness
    virtual Float_t Dy() {return 0.;}
    // Cell size in x
    virtual Float_t Dpx(Int_t) {return 0.;}
    // Cell size in z 
    virtual Float_t Dpz(Int_t) {return 0.;}
    // Maximum number of Cells in x
    virtual Int_t    Npx() {return 0;}
    // Maximum number of Cells in z
    virtual Int_t    Npz() {return 0;}
    // Layer
    virtual Int_t GetLayer() const {return 0;}
    // Angles 
    virtual void Angles(Float_t &, Float_t&) {}
    // Set cell position
    virtual void     SetPad(Int_t, Int_t) {}
    // Set hit position
    virtual void     SetHit(Float_t, Float_t) {}
    
    //
    // Iterate over cells 
    // Initialiser
    virtual void  FirstPad(Float_t,Float_t,Float_t,Float_t) {}
    // Stepper
    virtual void  NextPad() {}
    // Condition
    virtual Int_t MorePads() {return 0;}
    //
    // Get next neighbours 
    virtual void Neighbours(Int_t,Int_t,Int_t*,Int_t[10],Int_t[10]) {}
    //
    // Current cell cursor during disintegration
    // x-coordinate
    virtual Int_t  Ix() {return 0;}
    // z-coordinate
    virtual Int_t  Iz() {return 0;}
    //
    // Signal Generation Condition during Stepping
    virtual Int_t SigGenCond(Float_t,Float_t,Float_t) {return 0;}
    // Initialise signal generation at coord (x,y,z)
    virtual void  SigGenInit(Float_t,Float_t,Float_t) {}
    // Current integration limits 
    virtual void  IntegrationLimits(Float_t&,Float_t&,Float_t&,Float_t&) {}
    // Test points for auto calibration
    virtual void GiveTestPoints(Int_t &,Float_t *,Float_t *) {}
    // Function for systematic corrections
    // Set the correction function
    virtual void SetCorrFunc(Int_t, TF1*) {}
    // Get the correction Function
    virtual TF1* CorrFunc(Int_t) {return 0;}
	    
    ClassDef(AliITSsegmentation,1) //Segmentation virtual base class 
};
#endif
