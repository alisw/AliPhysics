#ifndef MUONSegResV1_H
#define MUONSegResV1_H
/////////////////////////////////////////////////////////
//  Manager and hits classes for set:MUON version 0    //
/////////////////////////////////////////////////////////
 
#include "AliMUONSegRes.h"

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
    virtual Float_t Dpx(Int_t ){return fDpx;}
    virtual Float_t Dpy(Int_t ){return fDpy;}
    virtual Int_t   Npx(){return fNpx;}
    virtual Int_t   Npy(){return fNpy;}
    virtual Float_t GetRealDpx(Int_t ) {return fDpx;}
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
    // set pad position
    virtual void     SetPad(Int_t, Int_t);
    // set hit position
    virtual void     SetHit(Float_t, Float_t);
    //
    // Iterate over pads
    virtual void SetPadCoord(Int_t iX, Int_t iY);
    virtual void  FirstPad(Float_t xhit, Float_t yhit, Float_t dx, Float_t dy);
    virtual void  NextPad();
    virtual Int_t MorePads();
    // Get next neighbours 
    virtual void Neighbours // implementation Neighbours function
	(Int_t iX, Int_t iY, Int_t* Nlist, Int_t *Xlist, Int_t *Ylist);
    virtual void NeighboursDiag // with diagonal elements
	(Int_t iX, Int_t iY, Int_t* Nlist, Int_t *Xlist, Int_t *Ylist);
    virtual void NeighboursNonDiag // without diagonal elements
	(Int_t iX, Int_t iY, Int_t* Nlist, Int_t *Xlist, Int_t *Ylist);
    void CleanNeighbours(Int_t* Nlist, Int_t *Xlist, Int_t *Ylist);
    // Channel number expressed in pad coordinates (stored in Cluster)
    virtual Int_t Ix(Int_t trueX, Int_t trueY);
    virtual Int_t Ix();
    virtual Int_t Iy(){return fiy;}
    // Actual number of pad in the chain
    virtual Int_t ISector();
    virtual Int_t Sector(Int_t , Int_t ) {return 1;}
    // Position of pad in perellel read-out
    virtual Int_t IsParallel2(Int_t iX, Int_t iY);
    virtual Int_t IsParallel3(Int_t iX, Int_t iY);
    // Number of pads read in parallel
    virtual Int_t NParallel2(Int_t iX, Int_t iY);
    virtual Int_t NParallel3(Int_t iX, Int_t iY);
    //
    // Number of pads read in parallel and offset to add to x
    virtual void GetNParallelAndOffset(Int_t iX, Int_t iY,
    	Int_t *Nparallel, Int_t *Offset);
    // Minimum distance between 1 pad and a position
    virtual Float_t Distance2AndOffset(Int_t iX, Int_t iY, Float_t X, Float_t Y, Int_t *Offset);
    //
    // Signal Generation Condition during Stepping
    Int_t SigGenCond(Float_t x, Float_t y, Float_t z);
    void  SigGenInit(Float_t x, Float_t y, Float_t z);
    void  GiveTestPoints(Int_t &n, Float_t *x, Float_t *y);
    virtual void IntegrationLimits
	(Float_t& x1, Float_t& x2, Float_t& y1, Float_t& y2);
    //
    virtual void Draw(Option_t *){;}
    // Function for systematic corrections
    virtual void SetCorrFunc(Int_t , TF1* func) {fCorr=func;}
    virtual TF1* CorrFunc(Int_t) {return fCorr;} 

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
    Float_t fxhit;
    Float_t fyhit;
    
    TF1* fCorr;
    
    //
    // Version Identifier
    char    *fName;       



};

#endif


