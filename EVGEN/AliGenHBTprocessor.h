#ifndef ALIGENHBTPROCESSOR_H
#define ALIGENHBTPROCESSOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Implementation of the interface for THBTprocessor
// Author: Piotr Krzysztof Skowronski <Piotr.Skowronski@cern.ch>

#include "AliGenerator.h"
#include <TFile.h>
#include <TTree.h>
#include <AliPDG.h>
#include "THBTprocessor.h"

enum {kHBTPMaxParticleTypes = 50};

class AliGenHBTprocessor : public AliGenerator { 

  public:
    AliGenHBTprocessor();
    virtual ~AliGenHBTprocessor();

    virtual void Init();
    virtual void Generate();
    virtual void GetParticles(TClonesArray * particles);
    Int_t        IdFromPDG(Int_t) const;
    Int_t        PDGFromId(Int_t) const;

    Int_t   GetHbtPStatusCode(Int_t part) const; 
    void    SetHbtPStatusCode(Int_t hbtstatcode, Int_t part);
    static const Int_t fgkHBTPMAXPART;
/************* S E T T E R S ******************/  

    virtual void SetTrackRejectionFactor(Float_t trf = 1.0);

    virtual void SetRefControl(Int_t rc =2);
    virtual void SetPIDs(Int_t pid1 = kPiPlus,Int_t pid2 = kPiMinus); //PDG Codes of particles to be processed, default \\Pi^{+} and \\Pi^{-}
    virtual void SetNPIDtypes(Int_t npidt = 2); //Number ofparticle types to be processed
    virtual void SetDeltap(Float_t deltp = 0.1); //maximum range for random momentum shifts in GeV/c;
                                                 //px,py,pz independent; Default = 0.1 GeV/c.
    virtual void SetMaxIterations(Int_t maxiter = 50);
    virtual void SetDelChi(Float_t dc = 0.1);
    virtual void SetIRand(Int_t irnd = 76564) ;
     
    virtual void SetLambda(Float_t lam = 0.6);
    virtual void SetR1d(Float_t r = 7.0) ;
    virtual void SetRSide(Float_t rs = 6.0);
    virtual void SetROut(Float_t ro = 7.0) ;
    virtual void SetRLong(Float_t rl = 4.0) ;
    virtual void SetRPerp(Float_t rp = 6.0);
    virtual void SetRParallel(Float_t rprl = 4.0);
    virtual void SetR0(Float_t r0 = 4.0) ;
    virtual void SetQ0(Float_t q0 = 9.0) ;
    virtual void SetSwitch1D(Int_t s1d = 3);
    virtual void SetSwitch3D(Int_t s3d = 0) ;
    virtual void SetSwitchType(Int_t st = 3);
    virtual void SetSwitchCoherence(Int_t sc = 0);
    virtual void SetSwitchCoulomb(Int_t scol = 2);
    virtual void SetSwitchFermiBose(Int_t sfb = 1);
    
    virtual void SetMomentumRange(Float_t pmin=0, Float_t pmax=0); //Dummy method
    virtual void SetPtRange(Float_t ptmin = 0.1, Float_t ptmax = 0.98);
    virtual void SetPxRange(Float_t pxmin = -1.0, Float_t pxmax = 1.0);
    virtual void SetPyRange(Float_t pymin = -1.0, Float_t pymax = 1.0);  
    virtual void SetPzRange(Float_t pzmin = -3.6, Float_t pzmax = 3.6);
    
    virtual void SetPhiRange(Float_t phimin = 0.0, Float_t phimax = 360.0);//Angle in degrees
                                                                                   //coherent with AliGenCocktail
                                                                                   //incohernet with AliGenerator
    virtual void SetEtaRange(Float_t etamin = -1.5, Float_t etamax = 1.5);//Pseudorapidity
    
    virtual void SetNPtBins(Int_t nptbin = 50);
    virtual void SetNPhiBins(Int_t nphibin = 50);
    virtual void SetNEtaBins(Int_t netabin = 50);
    virtual void SetNPxBins(Int_t npxbin = 20);
    virtual void SetNPyBins(Int_t npybin = 20);
    virtual void SetNPzBins(Int_t npzbin = 70);
   
    
    virtual void SetNBins1DFineMesh(Int_t n = 10);
    virtual void SetBinSize1DFineMesh(Float_t x=0.01);
      
    virtual void SetNBins1DCoarseMesh(Int_t n =2 );
    virtual void SetBinSize1DCoarseMesh(Float_t x=0.05);
      
    virtual void SetNBins3DFineMesh(Int_t n = 8);
    virtual void SetBinSize3DFineMesh(Float_t x=0.01);
      
    virtual void SetNBins3DCoarseMesh(Int_t n = 2);
    virtual void SetBinSize3DCoarseMesh(Float_t x=0.08);
      
    virtual void SetNBins3DFineProjectMesh(Int_t n =3 );
/***********************************************************************/
/* * * * * * *    P R O T E C T E D   A R E A    * * * * * * * * * * * */ 
/***********************************************************************/
  protected:
    
    THBTprocessor * fHBTprocessor;       //pointer to generator (TGenerator)
    Int_t         **fHbtPStatCodes;      //! hbtp status codes of particles
    Int_t           fNPDGCodes;          //! Number of defined particles   
    Int_t           fPDGCode[kHBTPMaxParticleTypes]; //! PDG codes (for conversion PDG<->Geant)
    void            DefineParticles();   //initiates array with PDG codes
    void            InitStatusCodes();   //Initiates status codes (allocates memory and sets everything to zero) 
    void            CleanStatusCodes();   //deletes array with status codes
    /**********   P A R A M E T E R S  OF THE GENERATOR****************/
	   
    Float_t fTrackRejectionFactor; //variates in range 0.0 <-> 1.0
                                   //Describes the factor of particles rejected from the output.
                                   //Used only in case of low muliplicity particles e.g. lambdas.
                                   //Processor generates addisional particles and builds the 
                                   //correletions on such a statistics.
                                   //At the end these particels are left in the event according 
                                   //to this factor: 1==all particles are left
                                   //                0==all are removed
      Int_t fReferenceControl;     //switch wether read reference histograms from file =1
                                   //              compute from input events =2 - default
      Int_t fPrintFull;             // Full print out option - each event
      Int_t fPrintSectorData;       // Print sector overflow diagnostics
      Int_t fNPidTypes;             // # particle ID types to correlate
      Int_t fPid[2];                // Geant particle ID #s, max of 2 types
      Int_t fNevents ;              // # events in input event text file
      Int_t fSwitch_1d;              // Include 1D correlations
      Int_t fSwitch_3d;              // Include 3D correlations
      Int_t fSwitch_type ;           // For like, unlike or both PID pairs
      Int_t fSwitch_coherence;       // To include incoh/coher mixed source
      Int_t fSwitch_coulomb;         // Coulomb correction selection options
      Int_t fSwitch_fermi_bose;      // For fermions or bosons

//   Numbers of particles and pairs:

      Int_t fN_part_1_trk;           // Total # PID #1 in 'trk', all flags
      Int_t fN_part_2_trk;           // Total # PID #2 in 'trk', all flags
      Int_t fN_part_tot_trk;         // Total # all part. in 'trk', all flgs
      Int_t fN_part_used_1_trk;      // # PID#1, used (flag=0) in 'trk'
      Int_t fN_part_used_2_trk;      // # PID#2, used (flag=0) in 'trk'

      Int_t fN_part_1_trk2;          // Total # PID #1 in 'trk2', all flags
      Int_t fN_part_2_trk2;          // Total # PID #2 in 'trk2', all flags
      Int_t fN_part_tot_trk2;        // Total # all part. in 'trk2', all flgs
      Int_t fN_part_used_1_trk2;     // # PID#1, used (flag=0) in 'trk2'
      Int_t fN_part_used_2_trk2;     // # PID#2, used (flag=0) in 'trk2'

      Int_t fN_part_used_1_ref;      // # PID#1, used (flag=0) in Reference
      Int_t fN_part_used_2_ref;      // # PID#2, used (flag=0) in Reference
      Int_t fN_part_used_1_inc;      // # PID#1, used (flag=0) in Inclusive 
      Int_t fN_part_used_2_inc;      // # PID#2, used (flag=0) in Inclusive

      Int_t fNum_pairs_like;         // # like pairs used (flag=0) in fit 
      Int_t fNum_pairs_unlike;       // # unlike pairs used (flag=0) in fit
      Int_t fNum_pairs_like_ref;     // # like pairs used (flag=0) in Ref. 
      Int_t fNum_pairs_unlike_ref;   // # unlike pairs used (flag=0) in Ref. 
      Int_t fNum_pairs_like_inc;     // # like pairs used (flag=0) in Incl. 
      Int_t fNum_pairs_unlike_inc;   // # unlike pairs used (flag=0) in Incl. 

//   Counters:

      Int_t fEvent_line_counter;     // Input event text file line counter
      Int_t fMaxit;                  // Max # iterations in track adjustment
      Int_t fIrand;                  // Random # starting seed (Def=12345)      
      Int_t fFile10_line_counter;    // Output, correlated event text file
//                                    //    line counter

//   Correlation Model Parameters:

      Float_t    fLambda;               // Chaoticity parameter
      Float_t    fR_1d;                   // Spherical source radius (fm)
      Float_t    fRside;                  // 3D Bertsch-Pratt source 'side' R (fm)
      Float_t    fRout;                   // 3D Bertsch-Pratt source 'out'  R (fm)
      Float_t    fRlong;                  // 3D Bertsch-Pratt source 'long' R (fm)
      Float_t    fRperp;                  // 3D YKP source transverse radius  (fm)
      Float_t    fRparallel;              // 3D YKP source longitudinal radius(fm)
      Float_t    fR0;                     // 3D YKP source emission time durat(fm)
      Float_t    fQ0;                     // NA35 Coulomb parameter (GeV/c) or
//                                    // Coul radius for Pratt finite src (fm)

//   Search Control Parameters:


      Float_t    fDeltap;                 // Max limit for x,y,z momt shifts(GeV/c)
      Float_t    fDelchi;                 // Min% change in Chi-Sq to stop iterat.


//   Chi-Square Values:

      Float_t    fChisq_wt_like_1d;          // 1D, Like pairs
      Float_t    fChisq_wt_unlike_1d;        // 1D, Unlike pairs
      Float_t    fChisq_wt_like_3d_fine;     // 3D, Like pairs, Fine Mesh
      Float_t    fChisq_wt_unlike_3d_fine;   // 3D, Unlike pairs, Fine Mesh
      Float_t    fChisq_wt_like_3d_coarse;   // 3D, Like pairs, Coarse Mesh
      Float_t    fChisq_wt_unlike_3d_coarse; // 3D, Unlike pairs, Coarse Mesh
      Float_t    fChisq_wt_hist1_1;          // One-body, particle ID type #1
      Float_t    fChisq_wt_hist1_2;          // One-body, particle ID type #2

//   Particle Masses:

      Float_t    fMass1, fMass2;           // Particle ID# 1 and 2 masses (GeV)


  /**********   M E S H  ****************/      


      Int_t fN_pt_bins;                  // # one-body pt bins
      Int_t fN_phi_bins;                 // # one-body phi bins
      Int_t fN_eta_bins;                 // # one-body eta bins
     
      Int_t fN_1d_fine;                  // # bins for 1D, Fine Mesh
      Int_t fN_1d_coarse;                // # bins for 1D, Coarse Mesh
      Int_t fN_1d_total;                 // Total # bins for 1D
      Int_t fN_3d_fine ;                 // # bins for 3D, Fine Mesh
      Int_t fN_3d_coarse;                // # bins for 3D, Coarse Mesh
      Int_t fN_3d_total;                 // Total # bins for 3D
      Int_t fN_3d_fine_project;          // # 3D fine mesh bins to sum over for

//   Momentum Space Sectors for Track Sorting:

      Int_t fN_px_bins;                  // # sector bins in px
      Int_t fN_py_bins;                  // # sector bins in py
      Int_t fN_pz_bins;                  // # sector bins in pz
      Int_t fN_sectors;                  // Total # sectors in 3D momentum space

//   Temporary Momentum Space Sector information storage during trk adjust.

      Int_t fOld_sec_ntrk;               // Old sector # tracks
      Int_t fOld_sec_flag;               // Old sector flag value
      Int_t fOld_sec_trkid[MAX_TRK_SAVE];         // Old sector track id array

      Int_t fNew_sec_ntrk;               // New sector # tracks
      Int_t fNew_sec_flag;               // New sector flag value
      Int_t fNew_sec_trkid[MAX_TRK_SAVE];// New sector track id array
      Int_t fNew_sec_save;               // New sector ID value
      Int_t fNld_sec_save;               // Old sector ID value
     
      Float_t    fPt_bin_size ;          // One-body pt bin size in (GeV/c)

      
      Float_t    fPhi_bin_size;          // One-body phi bin size in (degrees)
      
      Float_t    fEta_bin_size ;         // One-body eta bin size
      Float_t    fEta_min;               // One-body eta min/max
      Float_t    fEta_max;
//   Two-Body Histograms and Correlation Mesh for 1D and 3D distributions:
//                                       // projections onto single axis.

      Float_t    fBinsize_1d_fine;       // Bin Size - 1D, Fine Mesh in (GeV/c)
      Float_t    fBinsize_1d_coarse;     // Bin Size - 1D, Coarse Mesh in (GeV/c)
      Float_t    fQmid_1d;               // q (GeV/c) at fine-coarse mesh boundary
      Float_t    fQmax_1d;               // Max q (GeV/c) for 1D distributions
      Float_t    fBinsize_3d_fine;       // Bin Size - 3D, Fine Mesh in (GeV/c)
      Float_t    fBinsize_3d_coarse;     // Bin Size - 3D, Coarse Mesh in (GeV/c)
      Float_t    fQmid_3d;               // q (GeV/c) at fine-coarse mesh boundary
      Float_t    fQmax_3d;               // Max q (GeV/c) for 3D distributions

      Float_t    fPx_min;                // Sector range in px in GeV/c
      Float_t    fPx_max;                //--//--
      Float_t    fDelpx;                 // Mom. space sector cell size - px(GeV/c)     
      
      Float_t    fPy_min;                // Sector range in py in GeV/c 
      Float_t    fPy_max;                // --//--
      Float_t    fDelpy;                 // Mom. space sector cell size - py(GeV/c)     

      Float_t    fPz_min;                // Sector range in pz in GeV/c min
      Float_t    fPz_max;                // Sector range in pz in GeV/c max
      Float_t    fDelpz;                 // Mom. space sector cell size - pz(GeV/c)

  public:  
    ClassDef(AliGenHBTprocessor,1) // Interface class for AliMevsim
    
};
#endif
