#ifndef ALIANALYSISTASKMCPARTICLEFILTER_H
#define ALIANALYSISTASKMCPARTICLEFILTER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//
//  Analysis task for Kinematic filtering
//  Fill AOD tracks from Kinematic stack
//

#include <TList.h> 
#include "AliAnalysisTaskSE.h"

class AliAnalysisFilter;
class TString;

class AliAnalysisTaskMCParticleFilter : public AliAnalysisTaskSE
{
 public:
                                  AliAnalysisTaskMCParticleFilter();
                                  AliAnalysisTaskMCParticleFilter( const char* name );
                                  AliAnalysisTaskMCParticleFilter(const AliAnalysisTaskMCParticleFilter& obj);
    virtual                      ~AliAnalysisTaskMCParticleFilter();
     AliAnalysisTaskMCParticleFilter&   operator=(const AliAnalysisTaskMCParticleFilter& other);
    
    // Implementation of interface methods
    virtual                void   UserCreateOutputObjects();
    virtual                void   UserExec( Option_t *option );
    
    // Setters
    virtual                void   SetTrackFilterMother(AliAnalysisFilter* trackF) { fTrackFilterMother = trackF; }
    
 private:
    Bool_t Select(TParticle* part, Float_t rv, Float_t zv);                

    AliAnalysisFilter*  fTrackFilterMother;   //  Track Filter
                
    ClassDef( AliAnalysisTaskMCParticleFilter, 1 ); // Analysis task for Kinematic filtering
};
 
#endif
