#ifndef ALIAODPARTICLECUT_H
#define ALIAODPARTICLECUT_H
//__________________________________________________________________________
////////////////////////////////////////////////////////////////////////////
//                                                                        //
// class AliAODParticleCut                                                //
//                                                                        //
// Classes for single particle cuts.                                      //
// User should use mainly AliAODParticleCut interface methods,            //
// eventually EmptyCut which passes all particles.                        //
//                                                                        //
// There is all interface for setting cuts on all particle properties     //
// The main method is Rejected - which returns                            //
//         True to reject particle                                        //
//         False in case it meets all the criteria of the given cut       //
//                                                                        //
// This class has the list of base particle  cuts that perform check on   //
// single property. Particle  is rejected if any of cuts rejects it.      //
// There are implemented logical base cuts that perform logical           //
// operations on results of two other base cuts. Using them user can      //
// create a tree structure of a base cuts that performs sophisticated     //
// cut.                                                                   //
//                                                                        //
// User can also implement a base cut that performs complicated           //
// calculations, if it is only more convenient and/or efficint.           //
//                                                                        //
// User should delete created cuts  himself                               //
// because when setting a cut, other objects (functions,analyses,         //
// readers, other cuts) make their own copy of a cut.                     //
//                                                                        //
// more info: http://aliweb.cern.ch/people/skowron/analyzer/index.html    //
// responsible: Piotr Skowronski@cern.ch                                  //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include "AliVAODParticle.h"
#include "AliAODParticleBaseCut.h"


class AliAODParticleEmptyCut;
class AliAODParticleCut;
class AliAODParticleBaseCut;


/******************************************************************/
/******************************************************************/
/******************************************************************/

/******************************************************************/
/******************************************************************/
/******************************************************************/

class AliAODParticleCut: public TObject
{
//Class describing cut on particle
  public:

    AliAODParticleCut();
    AliAODParticleCut(const AliAODParticleCut& in);
    virtual ~AliAODParticleCut();
    AliAODParticleCut& operator = (const AliAODParticleCut& in);
    
    virtual Bool_t Rejected(AliVAODParticle* p) const;
    Bool_t IsEmpty() const {return kFALSE;}
    
    void AddBasePartCut(AliAODParticleBaseCut* basecut);
    
    Int_t GetPID() const { return fPID;}
    void SetPID(Int_t pid){fPID=pid;}
    void SetMomentumRange(Double_t min, Double_t max);
    void SetPRange(Double_t min, Double_t max){SetMomentumRange(min,max);}
    void SetPtRange(Double_t min, Double_t max);
    void SetEnergyRange(Double_t min, Double_t max);
    void SetRapidityRange(Double_t min, Double_t max);
    void SetYRange(Double_t min, Double_t max){SetRapidityRange(min,max);}
    void SetPseudoRapidityRange(Double_t min, Double_t max);
    void SetPxRange(Double_t min, Double_t max);
    void SetPyRange(Double_t min, Double_t max);
    void SetPzRange(Double_t min, Double_t max);
    void SetPhiRange(Double_t min, Double_t max);
    void SetThetaRange(Double_t min, Double_t max);
    void SetVxRange(Double_t min, Double_t max);
    void SetVyRange(Double_t min, Double_t max);
    void SetVzRange(Double_t min, Double_t max);
    
    void Print(void) const;
  protected:
     
    AliAODParticleBaseCut* FindCut(AliAODParticleBaseCut::EAODCutProperty property);

    AliAODParticleBaseCut ** fCuts;//! Array with cuts
    Int_t fNCuts; //number of base cuts stored in fCuts

    Int_t fPID; //particle PID  - if=0 (rootino) all pids are accepted
          
  private:
    static const Int_t fgkMaxCuts; //Size of the fCuts array

    ClassDef(AliAODParticleCut,1)
};
/******************************************************************/
/******************************************************************/
/******************************************************************/

class AliAODParticleEmptyCut:  public AliAODParticleCut
{
//Empty - it passes possitively all particles - it means returns always False
//Class describing cut on particles
  public:
    AliAODParticleEmptyCut(){};
    virtual ~AliAODParticleEmptyCut(){};
    
    Bool_t Rejected(AliVAODParticle*) const {return kFALSE;} //accept everything
    Bool_t IsEmpty() const {return kTRUE;}

    ClassDef(AliAODParticleEmptyCut,1)
 
};

/******************************************************************/
/******************************************************************/
/******************************************************************/


#endif
