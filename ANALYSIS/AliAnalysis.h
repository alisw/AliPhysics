#ifndef ALIANALYSIS_H
#define ALIANALYSIS_H
//________________________________
///////////////////////////////////////////////////////////
//
// class AliAnalysis
//
// Base class for analysis
//
//
// Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////

#include <TTask.h>
#include <AliAODParticleCut.h>
#include <AliAODPairCut.h>

class AliAOD;
class AliStack;
class AliEventCut;
class AliVAODParticle;
class AliAODPair;
 
class AliAnalysis: public TTask
{
  public: 
    AliAnalysis();
    AliAnalysis(const char* name,const char* title);
    virtual ~AliAnalysis();
    
    virtual Int_t Init() = 0;
    virtual Int_t ProcessEvent(AliAOD* aodrec, AliAOD* aodsim = 0x0) = 0;
    virtual Int_t Finish() = 0;

    void          SetCutsOnRec();
    void          SetCutsOnSim();
    void          SetCutsOnRecAndSim();
    
    void          SetEventCut(AliEventCut* evcut);
    void          SetPairCut(AliAODPairCut* cut);
    
  protected:
    Bool_t        Rejected(AliAOD* recevent, AliAOD* simevent);
    AliEventCut*  fEventCut;//event cut

    Bool_t        fCutOnSim;//flag indicating that event cut is performed on simulated particles 
    Bool_t        fCutOnRec;//flag indicating that event cut is performed on reconstructed tracks

    AliAODPairCut*   fPairCut;// Pair cut applied for all mixed particles

    /**********************************************/
    /*                C U T S                     */
    /**********************************************/

    Bool_t (AliAnalysis::*fkPass)(AliAODPair* partpair, AliAODPair* trackpair) const;//Pointer to function that performes pair cut
    Bool_t (AliAnalysis::*fkPass1)(AliVAODParticle* partpair, AliVAODParticle* trackpair) const;//Pointer to function that performes cut on first particle
    Bool_t (AliAnalysis::*fkPass2)(AliVAODParticle* partpair, AliVAODParticle* trackpair) const;//Pointer to function that performes cut on second particle
    Bool_t (AliAnalysis::*fkPassPairProp)(AliAODPair* partpair, AliAODPair* trackpair) const;//Pointer to function that performes pair cut

    Bool_t PassPartAndTrack (AliAODPair* partpair, AliAODPair* trackpair) const {return (fPairCut->Rejected((AliAODPair*)partpair))?kTRUE:fPairCut->Rejected((AliAODPair*)trackpair);}
    Bool_t PassPartAndTrack1(AliVAODParticle* part, AliVAODParticle* track) const;
    Bool_t PassPartAndTrack2(AliVAODParticle* part, AliVAODParticle* track) const;
    Bool_t PassPairPropPartAndTrack (AliAODPair* partpair, AliAODPair* trackpair) const {return (fPairCut->PassPairProp((AliAODPair*)partpair))?kTRUE:fPairCut->PassPairProp((AliAODPair*)trackpair);}

    Bool_t PassPart (AliAODPair* partpair, AliAODPair* /*trackpair*/) const {return fPairCut->Rejected((AliAODPair*)partpair);}
    Bool_t PassPart1(AliVAODParticle* part, AliVAODParticle* /*track*/) const {return fPairCut->GetFirstPartCut()->Rejected(part);}
    Bool_t PassPart2(AliVAODParticle* part, AliVAODParticle* /*track*/) const {return fPairCut->GetSecondPartCut()->Rejected(part);}
    Bool_t PassPairPropPart (AliAODPair* partpair, AliAODPair* /*trackpair*/) const {return fPairCut->PassPairProp((AliAODPair*)partpair);}

    Bool_t PassTrack (AliAODPair* /*partpair*/, AliAODPair* trackpair) const {return fPairCut->Rejected((AliAODPair*)trackpair);}
    Bool_t PassTrack1(AliVAODParticle* /*part*/, AliVAODParticle* track) const {return fPairCut->GetFirstPartCut()->Rejected(track);}
    Bool_t PassTrack2(AliVAODParticle* /*part*/, AliVAODParticle* track) const {return fPairCut->GetSecondPartCut()->Rejected(track);}
    Bool_t PassPairPropTrack (AliAODPair* /*partpair*/, AliAODPair* trackpair) const {return fPairCut->PassPairProp((AliAODPair*)trackpair);}

  private:
    ClassDef(AliAnalysis,1)
};


inline Bool_t AliAnalysis::PassPartAndTrack1(AliVAODParticle* part,AliVAODParticle* track) const
{
//Checks first particle from both, particle and track pairs
  AliAODParticleCut* pc = fPairCut->GetFirstPartCut();
  return (pc->Rejected(part))?kTRUE:pc->Rejected(track);
}
/*************************************************************************************/

inline Bool_t AliAnalysis::PassPartAndTrack2(AliVAODParticle* part,AliVAODParticle* track) const
{
//Checks second particle from both, particle and track pairs
  AliAODParticleCut* pc = fPairCut->GetSecondPartCut();
  return (pc->Rejected(part))?kTRUE:pc->Rejected(track);
}
/*************************************************************************************/ 

#endif
