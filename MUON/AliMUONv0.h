#ifndef MUONv0_H
#define MUONv0_H
/////////////////////////////////////////////////////////
//  Manager and hits classes for set:MUON version 0    //
/////////////////////////////////////////////////////////
 
#include "AliMUON.h"

class AliMUONv0 : public AliMUON {
 
public:
   AliMUONv0();
   AliMUONv0(const char *name, const char *title);
   virtual       ~AliMUONv0() {}
   virtual void   CreateGeometry();
   virtual void   CreateMaterials();
   virtual void   Init();
   virtual Int_t  IsVersion() const {return 0;}
   virtual void   StepManager();
   virtual void   Trigger(Float_t (*)[4], Float_t (*)[4], Int_t& iflag);
private:
   ClassDef(AliMUONv0,1)  //Hits manager for set:MUON version 0

};
class AliMUONsegmentationV0 :
public AliMUONsegmentation {
 public:
    AliMUONsegmentationV0(){}
    virtual ~AliMUONsegmentationV0(){}
    //    
    // Set Chamber Segmentation Parameters
    virtual  void    SetPADSIZ(Float_t p1, Float_t p2);
    virtual  void    SetDAnod(Float_t D) {fWireD = D;};
    //
    // Transform from pad (wire) to real coordinates and vice versa  
    virtual Float_t GetAnod(Float_t xhit);
    virtual void    GetPadIxy(Float_t x ,Float_t y ,Int_t   &ix,Int_t   &iy);
    virtual void    GetPadCxy(Int_t   ix,Int_t   iy,Float_t &x ,Float_t &y );
    //
    // Initialisation
    virtual void Init(AliMUONchamber*);
    //
    // Get member data
    virtual Float_t Dpx(){return fDpx;}
    virtual Float_t Dpy(){return fDpy;}
    virtual Int_t   Npx(){return fNpx;}
    virtual Int_t   Npy(){return fNpy;}
    //
    // Iterate over pads
    virtual void  FirstPad(Float_t xhit, Float_t yhit, Float_t dx, Float_t dy);
    virtual void  NextPad();
    virtual Int_t MorePads();
    // Get next neighbours 
    virtual void Neighbours
	(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10]);
    // Provisory RecCluster coordinates reconstructor
    virtual void FitXY(AliMUONRecCluster* Cluster,TClonesArray* MUONdigits);
    //
    // Current Pad during Integration
    virtual Int_t  Ix(){return fix;}
    virtual Int_t  Iy(){return fiy;}    
    virtual Int_t  ISector(){return 1;}
    //
    // Signal Generation Condition during Stepping
    virtual Int_t SigGenCond(Float_t x, Float_t y, Float_t z);
    virtual void  SigGenInit(Float_t x, Float_t y, Float_t z);
    virtual void IntegrationLimits
	(Float_t& x1, Float_t& x2, Float_t& y1, Float_t& y2);
    //
    // Identification
    ClassDef(AliMUONsegmentationV0,1)
	 protected:
    //
    // Implementation of the segmentation data
    // Version 0 models rectangular pads with the same dimensions all
    // over the cathode plane
    //
    //  geometry
    //
    Float_t    fDpx;           // x pad width per sector  
    Float_t    fDpy;           // y pad base width
    Int_t      fNpx;
    Int_t      fNpy;           // Number of pads in y
    Float_t    fWireD;         // wire pitch
    
    // Chamber region consideres during disintegration (lower left and upper right corner)
    //
    Int_t fixmin;
    Int_t fixmax;
    Int_t fiymin;
    Int_t fiymax;
    //
    // Current pad during integration (cursor for disintegration)
    Int_t fix;
    Int_t fiy;
    Float_t fx;
    Float_t fy;
    //
    // Current pad and wire during tracking (cursor at hit centre)
    Int_t fixt;
    Int_t fiyt;
    Int_t fiwt;
    Float_t fxt;
    Float_t fyt;
};

class AliMUONresponseV0 : //Mathieson response
public AliMUONresponse {
 public:
    AliMUONresponseV0(){}
    virtual ~AliMUONresponseV0(){}
    //
    // Configuration methods
    // 
    virtual void   SetRSIGM(Float_t p1) {fNsigma=p1;} 
    virtual void   SetMUCHSP(Float_t p1) {fChslope=p1;}
    virtual void   SetMUSIGM(Float_t p1, Float_t p2) {fChwX=p1; fChwY=p2;}
    virtual void   SetMAXADC(Float_t p1) {fadc_satm=p1;}
    // Mathieson parameters
    virtual void   SetSqrtKx3(Float_t p1) {fSqrtKx3=p1;};
    virtual void   SetKx2(Float_t p1) {fKx2=p1;};
    virtual void   SetKx4(Float_t p1) {fKx4=p1;};
    virtual void   SetSqrtKy3(Float_t p1) {fSqrtKy3=p1;};
    virtual void   SetKy2(Float_t p1) {fKy2=p1;};
    virtual void   SetKy4(Float_t p1) {fKy4=p1;};
    virtual void   SetPitch(Float_t p1) {fPitch=p1;};
    //
    // Get member data
    virtual Float_t Chslope() {return fChslope;}
    virtual Float_t ChwX() {return fChwX;}    
    virtual Float_t ChwY() {return fChwY;}        
    virtual Float_t Nsigma() {return fNsigma;}    
    virtual Float_t adc_satm() {return fadc_satm;}
    //  
    // Chamber response methods
    // Pulse height from scored quantity (eloss)
    virtual Float_t IntPH(Float_t eloss);
    // Charge disintegration
    virtual Float_t IntXY(AliMUONsegmentation * segmentation);
    // Identification
    //
    ClassDef(AliMUONresponseV0,1)
 protected:
    Float_t fChslope;         // Slope of the charge distribution
    Float_t fChwX;            // Width of the charge distribution in x
    Float_t fChwY;            // Width of the charge distribution in y
    Float_t fNsigma;          // Number of sigma's used for charge distribution
    Float_t fadc_satm;        // Maximum ADC channel
    Float_t fSqrtKx3;         // Mathieson parameters for x
    Float_t fKx2;
    Float_t fKx4;
    Float_t fSqrtKy3;         // Mathieson parameters for y
    Float_t fKy2;
    Float_t fKy4;
    Float_t fPitch;           //anode-cathode pitch
};

#endif







