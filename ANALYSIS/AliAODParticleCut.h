#ifndef ALIAODPARTICLECUT_H
#define ALIAODPARTICLECUT_H
//__________________________________________________________________________
////////////////////////////////////////////////////////////////////////////
//                                                                        //
// class AliAODParticleCut                                                //
//                                                                        //
// Classes for single particle cuts                                       //
// User should use only AliAODParticleCut, eventually                     //
// EmptyCut which passes all particles                                    //
// There is all interface for setting cuts on all particle properties     //
// The main method is Pass - which returns                                //
//         True to reject particle                                        //
//         False in case it meets all the criteria of the given cut       //
//                                                                        //
// User should create (and also destroy) cuts himself                     // 
// and then pass them to the Analysis And Function by a proper method     //
//                                                                        //
//                                                                        //
// more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html   //
// responsible: Piotr Skowronski@cern.ch                                   //
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
