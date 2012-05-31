// Implementation of the interface for THBTprocessor
// which is a wrapper itself to Fortran 
// program "HBT processor" written by Lanny Ray
// Author: Piotr Krzysztof Skowronski <Piotr.Skowronski@cern.ch>
// 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#ifndef ALIGENHBTPROCESSOR_H
#define ALIGENHBTPROCESSOR_H

#include "AliGenerator.h"
#include "TPDGCode.h"

class THBTprocessor;
class TClonesArray;
class TParticle;

enum {kHBTPMaxParticleTypes = 50};

class AliGenHBTprocessor : public AliGenerator 
{ 
//Wrapper class for THBTProcessor 
//which is a wrapper itself to Fortran 
//program "HBT processor" written by Lanny Ray
//
//Piotr.Skowronski@cern.ch

  public:
    AliGenHBTprocessor();
    virtual ~AliGenHBTprocessor();

    virtual void Init();
    virtual void Generate();
    virtual void GetParticles(TClonesArray * particles) const;
    Int_t        IdFromPDG(Int_t pdg) const;
    Int_t        PDGFromId(Int_t id) const;

    Int_t        GetHbtPStatusCode(Int_t part) const; 
    void         SetHbtPStatusCode(Int_t hbtstatcode, Int_t part);
/************* S E T T E R S ******************/  

    virtual void SetTrackRejectionFactor(Float_t trf = 1.0);

    virtual void SetRefControl(Int_t rc =2);
    virtual void SetPIDs(Int_t pid1 = kPiPlus,Int_t pid2 = kPiMinus); //PDG Codes of particles to be processed, default \\Pi^{+} and \\Pi^{-}
    virtual void SetNPIDtypes(Int_t npidt = 2); //Number ofparticle types to be processed
    virtual void SetDeltap(Float_t deltp = 0.1); //maximum range for random momentum shifts in GeV/c;
                                                 //px,py,pz independent; Default = 0.1 GeV/c.
    virtual void SetMaxIterations(Int_t maxiter = 50);//
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
    
    virtual void SetPhiRange(Float_t phimin = -180.0, Float_t phimax = 180.0);//Phi angle
    virtual void SetEtaRange(Float_t etamin = -1.5, Float_t etamax = 1.5);//Pseudorapidity
    void SetThetaRange(Float_t thetamin = 0, Float_t thetamax = 180); //Azimuthal angle, override AliGenerator method
                                                                      //which uses this, core fortran HBTProcessor uses Eta (pseudorapidity)
                                                                      //so these methods has to be synchronized         

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

    virtual void SetPrintFull(Int_t flag = 1);
    
/************* E V E N T   M E R G E ******************/  
    
    Int_t        GetNumberOfEvents();
    Int_t        GetNumberOfTracks();
    void         SetActiveEventNumber(Int_t n);
    TParticle*   GetTrack(Int_t n);
    void         SetNEventsToMerge(Int_t nev);
    
    
    //conveerts Eta (pseudorapidity) to etha(azimuthal angle). Returns radians 
    static Double_t EtaToTheta(Double_t arg){return 2.*TMath::ATan(TMath::Exp(-arg));}
    //converts tetha(azimuthal angle) to Eta (pseudorapidity). Argument in radians
    static Double_t ThetaToEta(Double_t arg);
    //converts Degrees To Radians
    static Double_t DegreesToRadians(Double_t arg){return arg*TMath::Pi()/180.;}
    //converts Radians To Degrees 
    static Double_t RadiansToDegrees(Double_t arg){return arg*180./TMath::Pi();}
    
    static Int_t  GetDebug() {return fgDebug;}
//    static Int_t  GetDebug() {return fgDebug;}
    
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
      Int_t fSwitch1d;              // Include 1D correlations
      Int_t fSwitch3d;              // Include 3D correlations
      Int_t fSwitchType ;           // For like, unlike or both PID pairs
      Int_t fSwitchCoherence;       // To include incoh/coher mixed source
      Int_t fSwitchCoulomb;         // Coulomb correction selection options
      Int_t fSwitchFermiBose;      // For fermions or bosons

//   Counters:

      Int_t fEventLineCounter;     // Input event text file line counter
      Int_t fMaxit;                  // Max # iterations in track adjustment
      Int_t fIrand;                  // Random # starting seed (Def=12345)      
//                                    //    line counter

//   Correlation Model Parameters:

      Float_t    fLambda;               // Chaoticity parameter
      Float_t    fR1d;                   // Spherical source radius (fm)
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


//   Particle Masses:


  /**********   M E S H  ****************/      


      Int_t fNPtBins;                  // # one-body pt bins
      Int_t fNPhiBins;                 // # one-body phi bins
      Int_t fNEtaBins;                 // # one-body eta bins
     
      Int_t fN1dFine;                  // # bins for 1D, Fine Mesh
      Int_t fN1dCoarse;                // # bins for 1D, Coarse Mesh
      Int_t fN1dTotal;                 // Total # bins for 1D
      Int_t fN3dFine ;                 // # bins for 3D, Fine Mesh
      Int_t fN3dCoarse;                // # bins for 3D, Coarse Mesh
      Int_t fN3dTotal;                 // Total # bins for 3D
      Int_t fN3dFineProject;          // # 3D fine mesh bins to sum over for

//   Momentum Space Sectors for Track Sorting:

      Int_t fNPxBins;                  // # sector bins in px
      Int_t fNPyBins;                  // # sector bins in py
      Int_t fNPzBins;                  // # sector bins in pz
      Int_t fNSectors;                  // Total # sectors in 3D momentum space

     
      Float_t    fPtBinSize ;          // One-body pt bin size in (GeV/c)

      
      Float_t    fPhiBinSize;          // One-body phi bin size in (degrees)
      
      Float_t    fEtaBinSize ;         // One-body eta bin size
      Float_t    fEtaMin;              // One-body eta min
      Float_t    fEtaMax;              // One-body eta max
//   Two-Body Histograms and Correlation Mesh for 1D and 3D distributions:
//                                       // projections onto single axis.

      Float_t    fBinsize1dFine;       // Bin Size - 1D, Fine Mesh in (GeV/c)
      Float_t    fBinsize1dCoarse;     // Bin Size - 1D, Coarse Mesh in (GeV/c)
      Float_t    fQmid1d;               // q (GeV/c) at fine-coarse mesh boundary
      Float_t    fQmax1d;               // Max q (GeV/c) for 1D distributions
      Float_t    fBinsize3dFine;       // Bin Size - 3D, Fine Mesh in (GeV/c)
      Float_t    fBinsize3dCoarse;     // Bin Size - 3D, Coarse Mesh in (GeV/c)
      Float_t    fQmid3d;               // q (GeV/c) at fine-coarse mesh boundary
      Float_t    fQmax3d;               // Max q (GeV/c) for 3D distributions

      Float_t    fPxMin;                // Sector range in px in GeV/c
      Float_t    fPxMax;                //--//--
      Float_t    fDelpx;                 // Mom. space sector cell size - px(GeV/c)     
      
      Float_t    fPyMin;                // Sector range in py in GeV/c 
      Float_t    fPyMax;                // --//--
      Float_t    fDelpy;                 // Mom. space sector cell size - py(GeV/c)     

      Float_t    fPzMin;                // Sector range in pz in GeV/c min
      Float_t    fPzMax;                // Sector range in pz in GeV/c max
      Float_t    fDelpz;                 // Mom. space sector cell size - pz(GeV/c)


      Int_t fEventMerge;                //number of events that are masked as an one event
      Int_t fActiveStack;               //current active stack

      static Int_t    fgDebug;          //debug level

      /******* P R O T E C T E D   M E T H O D S  *****/
      void GetTrackEventIndex(Int_t n, Int_t &evno, Int_t &index) const; //returns event(stack) number and 

 private:
      AliGenHBTprocessor(const AliGenHBTprocessor& in);
      AliGenHBTprocessor & operator=(const AliGenHBTprocessor& in);

    ClassDef(AliGenHBTprocessor,1) // Interface class for AliMevsim
    
};
#endif
