#ifndef ALIPYTHIA_H
#define ALIPYTHIA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TPythia6.h>
#include <AliRndm.h>
#include <AliStructFuncType.h>
typedef enum
{kPyCharm, kPyBeauty, kPyCharmUnforced, kPyBeautyUnforced, kPyJpsi, kPyJpsiChi, kPyMb, kPyMbNonDiffr, kPyJets, kPyDirectGamma, kPyCharmPbMNR, kPyD0PbMNR, kPyBeautyPbMNR}
Process_t;
/*
typedef enum
{
    kDOSet1     = 1006,
    kGRVLO      = 5005,
    kGRVHO      = 5006,
    kMRSDminus  = 3031,
    kMRSD0      = 3030,
    kMRSG       = 3041,
    kCTEQ2pM    = 4024,
    kCTEQ4L     = 4032,
    kCTEQ4M     = 4034,
    kMRSTcgLO   = 3072,
    kCTEQ5L     = 4046,
    kGRVLO98    = 5012
}
StrucFunc_t;
*/
class AliPythia : public TPythia6, public AliRndm
{

 public:
    virtual ~AliPythia(){;}
    // convert to compressed code and print result (for debugging only)
    virtual Int_t CheckedLuComp(Int_t kf);
    // Pythia initialisation for selected processes
    virtual void ProcInit
	(Process_t process, Float_t energy, StrucFunc_t strucfunc);
    // treat protons as inside nuclei
    virtual void    SetNuclei(Int_t a1, Int_t a2);
    // Print particle properties
    virtual void PrintParticles();
    virtual void ResetDecayTable();
    virtual void SetDecayTable();
    // return instance of the singleton
    static  AliPythia* Instance();

 protected:
    Process_t     fProcess;           // Process type
    Float_t       fEcms;              // Centre of mass energy
    StrucFunc_t   fStrucFunc;         // Structure function
    Int_t         fDefMDCY[501];      //  ! Default decay switches per particle
    Int_t         fDefMDME[2000];     //  ! Default decay switches per mode
    static AliPythia*    fgAliPythia; // Pointer to single instance
 private: 
    AliPythia();

    ClassDef(AliPythia,1) //ALICE UI to PYTHIA
};

#endif



