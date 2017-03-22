#ifndef ALIPYTHIA8_H
#define ALIPYTHIA8_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliPythia.h,v 1.22 2007/10/09 08:43:24 morsch Exp $ */

#include "Analysis.h"
#include "AliPythiaBase.h"
#include "AliTPythia8.h"

class AliStack;
class AliPythia8 :public AliTPythia8, public AliPythiaBase
{

 public:
    AliPythia8();
    AliPythia8(const AliPythia8& pythia);
    virtual ~AliPythia8() {;}
    virtual Int_t Version() {return (8);}
    // convert to compressed code and print result (for debugging only)
    virtual Int_t CheckedLuComp(Int_t /*kf*/) {return -1;}
    // Pythia initialisation for selected processes
    virtual void  ProcInit (Process_t process, Float_t energy, StrucFunc_t strucfunc, Int_t tune);
    virtual void  SetSeed(UInt_t seed);
    virtual void  GenerateEvent();
    virtual void  GenerateMIEvent();
    virtual void  HadronizeEvent();
    virtual Int_t GetNumberOfParticles() {return GetN();}
    virtual void  SetNumberOfParticles(Int_t i);
    virtual void  EditEventList(Int_t i);
    virtual void  PrintStatistics();
    virtual void  EventListing();    
    virtual Int_t GetParticles(TClonesArray *particles) {return ImportParticles(particles, "All");}
    // Treat protons as inside nuclei
    virtual void  SetNuclei(Int_t a1, Int_t a2);
    // Print particle properties
    virtual void PrintParticles();
    // Reset the decay table
    virtual void ResetDecayTable();
    //
    // Common Physics Configuration
    virtual void SetPtHardRange(Float_t ptmin, Float_t ptmax);
    virtual void SetYHardRange(Float_t ymin, Float_t ymax);
    virtual void SetFragmentation(Int_t flag);
    virtual void SetInitialAndFinalStateRadiation(Int_t flag1, Int_t flag2);
    virtual void SetIntrinsicKt(Float_t kt);
    virtual void SwitchHFOff();
    virtual void SetPycellParameters(Float_t etamax, Int_t neta, Int_t nphi,
				     Float_t thresh, Float_t etseed, Float_t minet, Float_t r);
    virtual void GetJet(Int_t i, Float_t& px, Float_t& py, Float_t& pz, Float_t& e);
    virtual void ModifiedSplitting();
    virtual void InitQuenching(Float_t bmin, Float_t bmax, Float_t k,
			       Int_t iECMethod, Float_t zmax = 0.97,
			       Int_t ngmax = 30);
    virtual void SwitchHadronisationOff();
    virtual void SwitchHadronisationOn();
    //
    // Common Getters
    virtual void    GetXandQ(Float_t& x1, Float_t& x2, Float_t& q);
    virtual Float_t GetXSection();
    virtual Float_t GetPtHard();
    virtual Int_t GetNMPI() { return fLastNMPI; }
    virtual Int_t GetNSuperpositions() { return fLastNSuperposition; }

    //
    //
    virtual void  SetDecayTable();
    virtual void  Pystat(Int_t /*i*/){;}
    virtual void  Pylist(Int_t /*i*/){;}
    virtual Int_t ProcessCode();
    virtual void  Pycell(Int_t& nclus);
    virtual void  Pyclus(Int_t& nclus);
    virtual void  Pyshow(Int_t /*ip1*/, Int_t /*ip2*/, Double_t /*qmax*/) {;}
    virtual void  Pyrobo(Int_t /*imi*/, Int_t /*ima*/, Double_t /*the*/,
			 Double_t /*phi*/, Double_t /*bex*/, Double_t /*bey*/, Double_t /*bez*/) {;}
    static  AliPythia8* Instance();
    virtual void Quench() {;}
    // Assignment Operator
    AliPythia8 & operator=(const AliPythia8 & rhs);
    void Copy(TObject&) const;
    //
    // Not yet implemented
    //
    virtual void Pyquen(Double_t, Int_t, Double_t);
    virtual void Pyexec() {;}
    virtual void GetQuenchingParameters(Double_t& xp, Double_t& yp, Double_t z[4]);
    virtual void LoadEvent(AliStack* stack, Int_t flag, Int_t reHadr);
    virtual void ConfigHeavyFlavor();
    virtual void AtlasTuning();
 protected:
    Process_t             fProcess;           // Process type
    Float_t               fEcms;              // Centre of mass energy
    StrucFunc_t           fStrucFunc;         // Structure function
    Int_t                 fDefMDCY[501];      //  ! Default decay switches per particle
    Int_t                 fDefMDME[2001];     //  ! Default decay switches per mode
    Double_t              fZQuench[4];        //  ! Quenching fractions for this even
    Pythia8::CellJet      fCellJet;           //  ! Cell jet object
    Float_t               fEtSeed;            //  ! Et seed for cell jets 
    Float_t               fMinEtJet;          //  ! Min jet et 
    Float_t               fRJet;              //  ! Radius for cell jets
    Pythia8::ClusterJet   fClusterJet;        //  ! Cluster jet object
    Float_t               fYScale;            //  ! cut-off joining scale
    Float_t               fPtScale;           //  ! cut-off joining scale
    Int_t                 fNJetMin;           //  ! min. number of jets
    Int_t                 fNJetMax;           //  ! max. number of jets
    static AliPythia8*    fgAliPythia8;       //    Pointer to single instance

    ClassDef(AliPythia8, 1) //ALICE UI to PYTHIA8
};

#endif





