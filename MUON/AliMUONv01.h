#ifndef MUONv01_H
#define MUONv01_H
/////////////////////////////////////////////////////
//  Segmentation and Response classes version 01   //
/////////////////////////////////////////////////////
 
#include "AliMUON.h"
#include "TArrayF.h"
#include "TArrayI.h"
class AliMUONsegmentationV01 :
public AliMUONsegmentation {
 public:
    AliMUONsegmentationV01();
    virtual ~AliMUONsegmentationV01(){}
    //    
    // Set Chamber Segmentation Parameters
    virtual  void    SetPADSIZ(Float_t p1, Float_t p2);
    virtual  void    SetDAnod(Float_t D) {fWireD = D;};
    virtual  void    SetPadDivision(Int_t ndiv[4]);
    virtual  void    SetSegRadii(Float_t  r[4]);
    //
    // Initialisation
    virtual void Init(AliMUONchamber*);
    //
    // Get member data
    virtual Float_t Dpx(){return fDpx;}
    virtual Float_t Dpy(){return fDpy;}
    virtual Int_t   Npx(){return fNpx[fNsec-1][1]+1;}
    virtual Int_t   Npy(){return fNpy;}
    //
    // Transform from pad (wire) to real coordinates and vice versa  
    virtual Float_t GetAnod(Float_t xhit);
    virtual void    GetPadIxy(Float_t x ,Float_t y ,Int_t   &ix,Int_t   &iy);
    virtual void    GetPadCxy(Int_t   ix,Int_t   iy,Float_t &x ,Float_t &y );
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
    virtual Int_t Ix()       {return fix;}
    virtual Int_t Iy()       {return fiy;}    
    virtual Int_t ISector()  {return fSector;}
    //
    // Integration
    virtual void IntegrationLimits
	(Float_t& x1, Float_t& x2, Float_t& y1, Float_t& y2);
    Int_t SigGenCond(Float_t x, Float_t y, Float_t z);
    void  SigGenInit(Float_t x, Float_t y, Float_t z);
    //
    // Identification
    ClassDef(AliMUONsegmentationV01,1)
    
 private:
    // Signal Generation Condition during Stepping
    void  GetSuperPadIxy(Float_t x, Float_t y, Int_t  &ix, Int_t  &iy);
    void  GetSuperPadCxy(Int_t  ix, Int_t  iy, Float_t &x, Float_t &y);
    Int_t Sector(Int_t ix, Int_t iy);
 protected:
    //
    // Implementation of the segmentation data
    // Version 0 models rectangular pads with the same dimensions all
    // over the cathode plane
    //
    //  geometry
    //
    Int_t      fNsec;          // Number of sectors
    TArrayF    fRSec;          // sector outer radia
    TArrayI    fNDiv;          // pad size division
    Float_t    fDpx;           // x pad width per sector  
    Float_t    fDpy;           // y pad base width
    TArrayF    fDpxD;          // y pad width per sector
    Int_t      fNpx[10][1000]; // Number of pads per sector in x
    Float_t    fCx[10][1000];  // pad-sector contour x vs y  
    Int_t      fNpy;           // Number of pads in y
    Float_t    fWireD;         // wire pitch
    
    // Chamber region consideres during disintegration
    // (lower left and upper right corner)
    //
    Int_t  fixmin;
    Int_t  fixmax;
    Int_t  fiymin;
    Int_t  fiymax;
    Float_t fxmin;
    Float_t fxmax;
    //
    // Current pad during integration (cursor for disintegration)
    Int_t fix;
    Int_t fiy;
    Float_t fx;
    Float_t fy;
    Int_t   fSector;
    //
    // Current pad and wire during tracking (cursor at hit centre)
    Int_t fixt;
    Int_t fiyt;
    Int_t fiwt;
    Float_t fxt;
    Float_t fyt;

};
#endif






