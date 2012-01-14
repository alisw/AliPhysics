#ifndef ALIANALYSISTASKMCPARTICLEFILTER_H
#define ALIANALYSISTASKMCPARTICLEFILTER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//
//  Analysis task for Kinematic filtering
//  Fill AOD tracks from Kinematic stack
//

#include "AliAnalysisTaskSE.h"

class AliAnalysisFilter;
class TString;
class TList;
class AliAODMCHeader;
class AliMCParticle;
class TClonesArray;

class AliAnalysisTaskMCParticleFilter : public AliAnalysisTaskSE
{
 public:
                                  AliAnalysisTaskMCParticleFilter();
                                  AliAnalysisTaskMCParticleFilter( const char* name );
    virtual                      ~AliAnalysisTaskMCParticleFilter();
    
    // Implementation of interface methods
    virtual                void   UserCreateOutputObjects();
    virtual                void   UserExec( Option_t *option );
    virtual                Bool_t Notify();
    virtual                void   Terminate( Option_t *option );
    // Setters
    virtual                void   SetTrackFilterMother(AliAnalysisFilter* trackF) { fTrackFilterMother = trackF; }
    
 private:
    Bool_t Select(TParticle* part, Float_t rv, Float_t zv);                
    
    // pivate c'tors to prevent misuse
    AliAnalysisTaskMCParticleFilter&   operator=(const AliAnalysisTaskMCParticleFilter& other);
    AliAnalysisTaskMCParticleFilter(const AliAnalysisTaskMCParticleFilter& obj);
    void PrintMCParticle(const AliMCParticle *mcp,Int_t np);


    AliAnalysisFilter*  fTrackFilterMother;   //  Track Filter
    AliAODMCHeader *fAODMcHeader;                 // ! MC Header
    TClonesArray *fAODMcParticles;            // ! tclones array of MC Particles
    TList *fHistList;                         // list to store e histograms, only as exchange

    ClassDef( AliAnalysisTaskMCParticleFilter, 4 ); // Analysis task for Kinematic filtering
};
 
#endif
