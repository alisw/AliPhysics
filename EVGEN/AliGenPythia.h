#ifndef AliGenPythia_H
#define AliGenPythia_H
/////////////////////////////////////////////////////////
//  Manager and hits classes for set:MUON version 0    //
/////////////////////////////////////////////////////////
#include "AliGenerator.h"
#include "TNamed.h"
#include "TF1.h"
#include "TArrayF.h"
#include "TTree.h"
#include "AliPythia.h"
#include "TArrayI.h"
#include "TParticle.h"

class AliGenPythia : public AliGenerator
{
 protected:
    Process_t   fProcess;
    StrucFunc_t fStrucFunc;
    Decay_t     fForceDecay;
    Float_t     fEnergyCMS;
    Float_t     fKineBias;
    Int_t       fTrials;
    TArrayI     fParentSelect;
    TArrayI     fChildSelect;
    Float_t     fXsection;
    AliPythia   *fPythia;
    Float_t     fPtHardMin;
    Float_t     fPtHardMax;    

 private:
    // check if particle is selected as parent particle
    Bool_t ParentSelected(Int_t ip);
    // check if particle is selected as child particle
    Bool_t ChildSelected(Int_t ip);
    // all kinematic selection cuts go here 
    Bool_t KinematicSelection(TParticle *particle);
    // adjust the weight from kinematic cuts
    void   AdjustWeights();
 public:
    AliGenPythia();
    AliGenPythia(Int_t npart);
    virtual ~AliGenPythia();
    virtual void    Generate();
    virtual void    Init();
    // select process type
    virtual void    SetProcess(Process_t proc=charm) {fProcess=proc;}
    // select structure function
    virtual void    SetStrucFunc(StrucFunc_t func=GRV_HO) {fStrucFunc=func;}
    // select pt of hard scattering 
    virtual void    SetPtHard(Float_t ptmin=0, Float_t ptmax=1.e10)
	{fPtHardMin=ptmin; fPtHardMax=ptmax; }
    // set centre of mass energy
    virtual void    SetEnergyCMS(Float_t energy=5500) {fEnergyCMS=energy;}
    // force decay type
    virtual void    ForceDecay(Decay_t decay=semimuonic) {fForceDecay=decay;}
    // get cross section of process
    virtual Float_t GetXsection() {return fXsection;}
    // Check PDG code
    virtual Int_t CheckPDGCode(Int_t pdgcode);

    ClassDef(AliGenPythia,1)
};
#endif





