#ifndef ALIPYTHIA_H
#define ALIPYTHIA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TPythia6.h>
#include <AliRndm.h>

typedef enum
{kPyCharm, kPyBeauty, kPyCharmUnforced, kPyBeautyUnforced, kPyJpsi, kPyJpsiChi, kPyMb, kPyJets, kPyDirectGamma}
Process_t;

typedef enum
{
    kDO_Set_1=1006,
    kGRV_LO=5005,
    kGRV_HO=5006,
    kMRS_D_minus=3031,
    kMRS_D0=3030,
    kMRS_G=3041,
    kCTEQ_2pM=4024,
    kCTEQ_4M=4034
}
StrucFunc_t;

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



