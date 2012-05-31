#ifndef ALIITSSEGMENTATIONSDD_H
#define ALIITSSEGMENTATIONSDD_H


#include "AliITSsegmentation.h"

//////////////////////////////////////////////////////
// Segmentation class for SDD                       //
//                                                  //
//////////////////////////////////////////////////////
class AliITSsegmentationSDD :
public AliITSsegmentation {
 public:


    AliITSsegmentationSDD(Option_t *opt="");
    AliITSsegmentationSDD(const AliITSsegmentationSDD &source);
    AliITSsegmentationSDD& operator=(const AliITSsegmentationSDD &source);
    virtual ~AliITSsegmentationSDD(){}

    // Set Detector Segmentation Parameters
    //
    // get the drift speed used for this segmentation
    virtual Float_t GetDriftSpeed()const{return fDriftSpeed;}
    // Set the drift speed needed for this segmentation
    virtual void SetDriftSpeed(Float_t ds){
      fDriftSpeed = ds;
      fSetDriftSpeed=kTRUE;
    }
    // Cell size dz*dx  
    virtual void    SetPadSize(Float_t pitch, Float_t clock) 
                         {fPitch=pitch;fTimeStep=1000./clock;}

    // Maximum number of cells along the two coordinates z,x (anodes,samples) 
    virtual void    SetNPads(Int_t p1, Int_t p2) 
                         {fNanodes=2*p1;fNsamples=p2;}
    // Returns the maximum number of cells (digits) posible
    virtual Int_t   GetNPads() const {return fNanodes*fNsamples;}

    // Transform from real local to cell coordinates
    virtual void    GetPadIxz(Float_t x ,Float_t z ,Int_t   &ix,Int_t   &iz) const;
    // Transform from cell to real local coordinates
    virtual void    GetPadCxz(Int_t   ix,Int_t   iz,Float_t &x ,Float_t &z ) const;    

    // Get anode and time bucket as floats - numbering from 0
    virtual void    GetPadTxz(Float_t &x ,Float_t &z) const;
    // Transformation from Geant cm detector center local coordinates
    // to detector segmentation/cell coordiantes starting from (0,0).
    virtual Bool_t  LocalToDet(Float_t x,Float_t z,Int_t &ix,Int_t &iz) const;
    // Transformation from detector segmentation/cell coordiantes starting
    // from (0,0) to Geant cm detector center local coordinates.
    virtual void    DetToLocal(Int_t ix,Int_t iz,Float_t &x,Float_t &z) const;
    //
    virtual Float_t GetAnodeFromLocal(Float_t x,Float_t z) const;
    virtual Int_t   GetSideFromLocalX(Float_t xloc) const {
      if(xloc>0) return 0;   // left side (channel 0) positive xloc
      else return 1;         // right side (channel 1) negative xloc
    }
    virtual Float_t GetLocalZFromAnode(Int_t nAnode) const;
    virtual Float_t GetLocalZFromAnode(Float_t zAnode) const;
    virtual Float_t GetDriftTimeFromTb(Int_t tb) const {
      Float_t xtb=(Float_t)tb+0.5;
      return GetDriftTimeFromTb(xtb);
    }
    virtual Float_t GetDriftTimeFromTb(Float_t xtb) const {
      return xtb*fTimeStep;
    }
    virtual Int_t    GetNumberOfChips() const {
      return fgkNchipsPerHybrid;
    }
    virtual Int_t    GetMaximumChipIndex() const{
      return fgkNchipsPerHybrid*2-1;
    }
    virtual Int_t    GetChipFromLocal(Float_t xloc, Float_t zloc) const;
    virtual Int_t    GetChipsInLocalWindow(Int_t* array, Float_t zmin, Float_t zmax, Float_t xmin, Float_t xmax) const;
      

    virtual Int_t    GetChipFromChannel(Int_t ix, Int_t iz) const; 

    // Initialisation
    virtual void Init();

    //
    // Get member data
    //
    // Cell size in x
    virtual Float_t Dpx(Int_t) const {return fTimeStep;}
    // Cell size in z 
    virtual Float_t Dpz(Int_t) const {return fPitch;} 

    // Maximum number of samples in x
    virtual Int_t    Npx() const {return fNsamples;}
    // Maximum number of anodes in z
    virtual Int_t    Npz() const {return fNanodes;}
    // Maximum number of anodes in z
    virtual Int_t    NpzHalf() const {return fNanodes/2;}

    //
    // Get next neighbours 
    virtual void Neighbours(Int_t iX,Int_t iZ,Int_t* Nlist,Int_t Xlist[10],
			    Int_t Zlist[10]) const;

    // Current integration limits 
    virtual void  IntegrationLimits(Float_t&,Float_t&,Float_t&,Float_t&) {}
    // Print Parameters
    virtual void    Print(Option_t *opt="") const;
// Print default parameters (static const data members, if any)
    virtual void PrintDefaultParameters() const;	    

    static Int_t GetNAnodesPerHybrid()  { return fgkHalfNanodesDefault;}
  protected:

    virtual void Copy(TObject &obj) const;

    Int_t      fNsamples; // Number of time samples in x
    Int_t      fNanodes;  // Summed # of anodes in the two det halves (z)
    Float_t    fPitch;    // Anode pitch - microns
    Float_t    fTimeStep; // Sampling time - ns
    Float_t    fDriftSpeed;  // Drift speed 
    Bool_t     fSetDriftSpeed; // flag for usage of drift speed
    static const Float_t fgkDxDefault;  // Default value for fDx
    static const Float_t fgkDzDefault;  // Default value for fDz
    static const Float_t fgkDyDefault;  // Default value for fDy
    static const Float_t fgkPitchDefault; //Default value for fPitch
    static const Float_t fgkClockDefault; //Default value for the clock freq.
    static const Int_t fgkHalfNanodesDefault; //Default value for fNanodes/2
    static const Int_t fgkNsamplesDefault; //Default value for fNsamples
    static const Int_t fgkNchipsPerHybrid;    //number of chips per hybrid
    static const Int_t fgkNanodesPerChip;    //number of chips per hybrid
    static const Float_t fgkCm2Micron;  // conversion from cm to micron
    static const Float_t fgkMicron2Cm;  // conversion from micron to cm
    ClassDef(AliITSsegmentationSDD,6) // SDD segmentation
};

#endif
