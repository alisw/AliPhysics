#ifndef ALIPYTHIABASE_H
#define ALIPYTHIABASE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
#include "AliRndm.h"
#include <TObject.h>
#include "AliStructFuncType.h"
#include "PythiaProcesses.h"

class AliFastGlauber;
class AliQuenchingWeights;
class AliStack;
class TClonesArray;

class AliPythiaBase : public AliRndm 
{

 public:
    AliPythiaBase();
    virtual ~AliPythiaBase(){;}
    void Dummy(){;}
    virtual Int_t Version() {return -1;}
    // Convert to compressed code and print result (for debugging only)
    virtual Int_t CheckedLuComp(Int_t /*kf*/) {return -1;}   
    // Pythia initialisation for selected processes
    virtual void  ProcInit (Process_t /*process*/, Float_t /*energy*/, StrucFunc_t /*strucfunc*/, Int_t /* tune */) {;}
    virtual void  SetSeed(UInt_t seed);
    virtual void  GenerateEvent() {;}
    virtual void  GenerateMIEvent() {;}
    virtual Int_t GetNumberOfParticles() {return -1;};
    virtual void  SetNumberOfParticles(Int_t /*i*/){;}
    virtual void  EditEventList(Int_t /*i*/) {;}
    virtual void  HadronizeEvent() {;}
    virtual Int_t GetParticles(TClonesArray */*particles*/){return -1;}
    virtual void  PrintStatistics() {;}
    virtual void  EventListing() {;}
    // Treat protons as inside nuclei
    virtual void  SetNuclei(Int_t /*a1*/, Int_t /*a2*/) {;}
    // Print particle properties
    virtual void PrintParticles() {;}
    // Reset the decay table
    virtual void ResetDecayTable() {;}
    //
    // Common Physics Configuration
    virtual void SetPtHardRange(Float_t /*ptmin*/, Float_t /*ptmax*/) {;}
    virtual void SetYHardRange(Float_t /*ymin*/, Float_t /*ymax*/) {;}
    virtual void SetFragmentation(Int_t /*flag*/) {;}
    virtual void SetInitialAndFinalStateRadiation(Int_t /*flag1*/, Int_t /*flag2*/) {;}
    virtual void SetIntrinsicKt(Float_t /*kt*/) {;}
    virtual void SwitchHFOff() {;}
    virtual void SetPycellParameters(Float_t /*etamax*/, Int_t /*neta*/, Int_t /*nphi*/,
				     Float_t /*thresh*/, Float_t /*etseed*/, Float_t /*minet*/, Float_t /*r*/) {;}
    virtual void ModifiedSplitting() {;}    
    virtual void SwitchHadronisationOff() {;}
    virtual void SwitchHadronisationOn() {;} 
    //
    // Common Getters
    virtual void    GetXandQ(Float_t& /*x1*/, Float_t& /*x2*/, Float_t& /*q*/) {;}
    virtual Float_t GetXSection() {return -1.;}
    virtual Int_t   ProcessCode() {return -1;}
    virtual Float_t GetPtHard() {return -1.;}
    virtual Int_t   GetNMPI() {return -1;}
    virtual Int_t   GetNSuperpositions() { return 1; }
    //
    //
    virtual void SetDecayTable() {;}
    virtual void Pycell(Int_t& /*nclus*/) {;}
    virtual void Pyclus(Int_t& /*nclus*/) {;}
    virtual void GetJet(Int_t /*i*/, Float_t& /*px*/, Float_t& /*py*/, Float_t& /*pz*/, Float_t& /*e*/){;}
    virtual void LoadEvent(AliStack* /*stack*/, Int_t /*flag*/, Int_t /*reHadr*/){;}
    virtual void Pyshow(Int_t /*ip1*/, Int_t /*ip2*/, Double_t /*qmax*/) {;}
    virtual void Pyrobo(Int_t /*imi*/, Int_t /*ima*/, Double_t /*the*/, Double_t /*phi*/, Double_t /*bex*/, 
			Double_t /*bey*/, Double_t /*bez*/){;}
    virtual void InitQuenching(Float_t /*bmin*/, Float_t /*bmax*/, Float_t /*k*/, Int_t /*iECMethod*/, 
			       Float_t /*zmax*/, Int_t /*ngmax*/) {;}
    virtual void Pyquen(Double_t /*a*/, Int_t /*ibf*/, Double_t /*b*/) {;}
    virtual void GetQuenchingParameters(Double_t& /*xp*/, Double_t& /*yp*/, Double_t z[4]) {;}
    // return instance of the singleton
    virtual void Quench() {;}
    virtual void ConfigHeavyFlavor() {;}
    virtual void AtlasTuning() {;}
    ClassDef(AliPythiaBase, 1) //ALICE UI to PYTHIA
};

#endif





