#ifndef ALIITSSEGMENTATIONSDD_H
#define ALIITSSEGMENTATIONSDD_H


#include "AliITSsegmentation.h"
#include "AliITSresponse.h"

// segmentation for SDD

class AliITSsegmentationSDD :
public AliITSsegmentation {
 public:


    AliITSsegmentationSDD();
    AliITSsegmentationSDD(AliITSgeom *gm, AliITSresponse *resp);
    AliITSsegmentationSDD(AliITSsegmentationSDD &source);
    virtual ~AliITSsegmentationSDD(){}
    AliITSsegmentationSDD& operator=(AliITSsegmentationSDD &source);

    // Set Detector Segmentation Parameters
    //
    // Detector size : x,z,y
  virtual  void   SetDetSize
          (Float_t p1=35000., Float_t p2=75264., Float_t p3= 300.) 
          {fDx=p1; fDz=p2; fDy=p3;}

    // Cell size dz*dx  
    virtual void    SetPadSize(Float_t pitch=294., Float_t clock=40.) 
                         {fPitch=pitch;fTimeStep=1000./clock;}

    // Maximum number of cells along the two coordinates z,x (anodes,samples) 
    virtual void    SetNPads(Int_t p1=256, Int_t p2=256) 
                         {fNanodes=2*p1;fNsamples=p2;}

    // Transform from real local to cell coordinates
    virtual void    GetPadIxz(Float_t x ,Float_t z ,Int_t   &ix,Int_t   &iz);
    // Transform from cell to real local coordinates
    virtual void    GetPadCxz(Int_t   ix,Int_t   iz,Float_t &x ,Float_t &z );
    // Transform from real global to local coordinates
    virtual void    GetLocal(Int_t module,Float_t *g ,Float_t *l);
    // Transform from real local to global coordinates
    virtual void    GetGlobal(Int_t module,Float_t *l ,Float_t *g);
    // Get anode and time bucket as floats - numbering from 0
    virtual void    GetPadTxz(Float_t &x ,Float_t &z);
    //
    // Initialisation
    virtual void Init();
    //
    // Get member data
    //
    // Detector type geometry
    virtual AliITSgeom* Geometry() {return fGeom;}
    // Detector length
    virtual Float_t Dx() {return fDx;}
    // Detector width
    virtual Float_t Dz()  {return fDz;}  
    // Detector thickness
    virtual Float_t Dy() {return fDy;}
    // Cell size in x
    virtual Float_t Dpx(Int_t dummy) {return fTimeStep;}
    // Cell size in z 
    virtual Float_t Dpz(Int_t dummy) {return fPitch;} 

    // Maximum number of samples in x
    virtual Int_t    Npx() {return fNsamples;}
    // Maximum number of anodes in z
    virtual Int_t    Npz() {return fNanodes;}

    //
    // Get next neighbours 
    virtual void Neighbours
      (Int_t iX, Int_t iZ, Int_t* Nlist, Int_t Xlist[10], Int_t Zlist[10]);

    // Set cell position
    virtual void     SetPad(Int_t, Int_t) {}
    // Set hit position
    virtual void     SetHit(Float_t, Float_t) {}
    
    //
    // Iterate over cells 
    // Initialiser
    virtual void  FirstPad
          (Float_t xhit, Float_t zhit, Float_t dx, Float_t dz) {}
    // Stepper
    virtual void  NextPad() {}
    // Condition
    virtual Int_t MorePads() {return 0;}
    //
    // Current cell cursor during disintegration
    // x-coordinate
    virtual Int_t  Ix() {return 0;}
    // z-coordinate
    virtual Int_t  Iz() {return 0;}
    //
    // Signal Generation Condition during Stepping
    virtual Int_t SigGenCond(Float_t x, Float_t y, Float_t z) {return 0;}
    // Initialise signal generation at coord (x,y,z)
    virtual void  SigGenInit(Float_t x, Float_t y, Float_t z) {}
    // Current integration limits 
    virtual void  IntegrationLimits
    (Float_t& x1, Float_t& x2, Float_t& z1, Float_t& z2) {}
    // Test points for auto calibration
    virtual void GiveTestPoints(Int_t &n, Float_t *x, Float_t *z) {}
    // Function for systematic corrections
    // Set the correction function
    virtual void SetCorrFunc(Int_t, TF1*) {}
    // Get the correction Function
    virtual TF1* CorrFunc(Int_t) {return 0;}
	    
  protected:

    Int_t      fNsamples;      // Number of time samples in x
    Int_t      fNanodes;       // Summed # of anodes in the two det halves (z)
    Float_t    fPitch;         // Anode pitch - microns
    Float_t    fTimeStep;      // Sampling time - ns
    Float_t    fDx   ;         // Full width of the detector (x axis) - microns
    Float_t    fDz    ;        // Length of half-detector (z axis) - microns
    Float_t    fDy;            // Full thickness of the detector (y axis)

    AliITSgeom *fGeom;         //! pointer to the geometry class
    AliITSresponse *fResponse; // pointer to the response class
   
    TF1*       fCorr;          // correction function

    ClassDef(AliITSsegmentationSDD,1) // SDD segmentation
};

#endif
