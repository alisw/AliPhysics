#ifndef MUONsegmentv1_H
#define MUONsegmentv1_H
/////////////////////////////////////////////////////////
//  Manager and hits classes for set:MUON version 0    //
/////////////////////////////////////////////////////////
 
#include "AliMUON.h"

const Int_t NZONE = 3; // Specific for chamber with equal pads
const Int_t NZONEm1 = 2; // NZONE - 1
const Int_t NZONECUT = 30;

class AliMUONsegmentationV1 :
public AliMUONsegmentation {
 public:
    AliMUONsegmentationV1();
    virtual ~AliMUONsegmentationV1(){}
    //    
    // Set Chamber Segmentation Parameters
    void SetNzone(Int_t N) {fNzone = N;};
    virtual  void    SetPADSIZ(Float_t p1, Float_t p2);
    void SetSensOffset(Float_t Offset) {fSensOffset = Offset;};
    void SetDAnod(Float_t D) {fDAnod = D;};
      // max x and y for the zone in number of pads units 
      //(WARNING : first pad is labelled 0 !!) 
    virtual void AddCut(Int_t Zone, Int_t nX, Int_t nY); 
    virtual void DefaultCut(void);

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
    // know the zone of segmentation
    virtual Int_t GetZone(Float_t X, Float_t Y);
    virtual Int_t GetZone(Int_t X, Int_t Y);
    //
    // Transform from pad (wire) to real coordinates and vice versa  
    virtual Int_t GetiAnod(Float_t xhit);
    virtual Float_t GetAnod(Float_t xhit);
    virtual void    GetPadIxy(Float_t x ,Float_t y ,Int_t   &ix,Int_t   &iy);
    virtual void    GetPadCxy(Int_t   ix,Int_t   iy,Float_t &x ,Float_t &y );
    //
    // Iterate over pads
    virtual void SetPadCoord(Int_t iX, Int_t iY);
    virtual void  FirstPad(Float_t xhit, Float_t yhit, Float_t dx, Float_t dy);
    virtual void  NextPad();
    virtual Int_t MorePads();
    // Get next neighbours 
    virtual void Neighbours
	(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10]);
    // Provisory RecCluster coordinates reconstructor
    virtual void FitXY(AliMUONRecCluster* Cluster,TClonesArray* MUONdigits);
    //
    // Channel number expressed in pad coordinates (stored in Cluster)
    virtual Int_t Ix();
    virtual Int_t Iy(){return fiy;}
    // Actual number of pad in the chain
    virtual Int_t ISector();
    //
    // Signal Generation Condition during Stepping
    Int_t SigGenCond(Float_t x, Float_t y, Float_t z);
    void  SigGenInit(Float_t x, Float_t y, Float_t z);
    virtual void IntegrationLimits
	(Float_t& x1, Float_t& x2, Float_t& y1, Float_t& y2);
    //
    // Identification
    virtual char* YourName() {return fName;}

    ClassDef(AliMUONsegmentationV1,1)
 protected:
    //
    // Implementation of the segmentation data
    // Version This models rectangular pads with the same dimensions all
    // over the cathode plane but let the possibilit for different design
    //
    //  geometry
    Int_t fNzone; // Number of differents sensitive zones
    Float_t fDpx;         // X pad width
    Float_t fDpy;         // Y pad width
    Int_t   fNZoneCut[NZONEm1];    // Number of cuts for given zone 
    Int_t fZoneX[NZONEm1][NZONECUT]; // X descriptor of zone segmentations
    Int_t fZoneY[NZONEm1][NZONECUT]; // Y descriptor of zone segmentations
    Float_t frSensMax2; // square of maximum sensitive radius
    Float_t frSensMin2; // square of minimum sensitive radius
    Int_t   fNpx;         // Maximum number of pads along x
    Int_t   fNpy;         // Maximum number of pads along y
//    Int_t   fNwire;       // Number of wires per pad
    Float_t fDAnod;       // Anod gap
    Float_t fSensOffset;  // Offset of sensitive zone with respect to quadrant (positive)
    
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
    //
    // Version Identifier
    char    *fName;       
};

#endif


