#ifndef ALILHC_H
#define ALILHC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Class for a simple description of the LHC.
// The LHC is described by two beams,
// the interaction regions and the
// beam loss processes.
// Run paramters can be set in order to simulate the time evolution
// of emittance, number of particles per bunch and luminosity.
// Author: Andreas Morsch
// andreas.morsch@cern.ch

#include <TObject.h>
#include <TList.h>
#include <TObjArray.h>

class AliLhcIRegion;
class AliLhcProcess;
class AliLhcBeam;


class AliLHC : public TObject
{
 public:
    AliLHC();
    AliLHC(const AliLHC &lhc);
    virtual ~AliLHC();

    virtual void AddIRegion(AliLhcIRegion *region);
    virtual void AddProcess(AliLhcProcess *process);
    virtual void SetBeams(AliLhcBeam* beam1, AliLhcBeam* beam2);
    virtual void SetTime(Float_t dt, Float_t tmax) 
	{fTimeStep = dt; fTimeMax = tmax;}
    virtual void SetFillingTime(Float_t t) {fFillingTime = t;}
    virtual void SetSetUpTime(Float_t t)   {fSetUpTime = t;}

    virtual void SetRadius(Float_t r) {fRadius = r;}
    virtual void SetAverageBeta(Float_t b) {fAverageBeta = b;}
    virtual void SetAverageDisp(Float_t b) {fAverageDisp = b;}
    
    virtual Float_t Radius()      const {return fRadius;}
    virtual Float_t AverageBeta() const {return fAverageBeta;}
    virtual Float_t AverageDisp() const {return fAverageDisp;}
    virtual Float_t SetUpTime()   const {return fSetUpTime;}
    virtual Float_t FillingTime() const {return fFillingTime;}


    virtual AliLhcBeam* Beam(Int_t i) 
      {return (AliLhcBeam*) (*fBeams)[i];}
    virtual TList* IRegions()   const {return fIRegions;}
    virtual void Init();
    virtual void EvolveTime();
    virtual void Evaluate();
    virtual Float_t  Time()     const {return fTime;}
    virtual Float_t  TimeStep() const {return fTimeStep;}
    virtual Float_t* TimeA()    const {return fTimeA;}
    virtual Int_t    Nt()       const {return fNt;}
    
    AliLHC & operator=(const AliLHC & rhs);
    
 protected:
    Int_t           fNRegions;              // Number of IR
    Int_t           fNProcesses;            // Number of processes
    TList*          fIRegions;              // List of intercation regions
    TList*          fProcesses;             // Beam processes
    TObjArray*      fBeams;                 // Lhc beams
    Float_t         fRadius;                // Radius (cm)
    Float_t         fAverageBeta;           // Average beta (cm)
    Float_t         fAverageDisp;           // Average dispersion (cm)

    Int_t    fNt;                  // Number of time steps
    Int_t  fNmax;                  // Max. Number of time steps

    Float_t  fTime;                // Current time
    Float_t* fTimeA;               // [fNmax] Current time
    Float_t  fTimeStep;            // Current time step
    Float_t  fTimeMax;             // Maximal time 
    //
    Float_t  fFillingTime;         // Filling Time
    Float_t  fSetUpTime;           // Set-up time
//
    ClassDef(AliLHC,1) // LHC manager class
};

#endif









