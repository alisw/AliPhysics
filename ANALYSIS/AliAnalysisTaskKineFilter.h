#ifndef ALIANALYSISTASKKINEFILTER_H
#define ALIANALYSISTASKKINEFILTER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//
//  Analysis task for Kinematic filtering
//  Fill AOD tracks from Kinematic stack
//

#include <TList.h> 
#include "AliAnalysisTaskSE.h"

class AliMCEvent;
class TChain;
class AliAODEvent;
class AliAnalysisFilter;
class TParticle;
class AliAODTrack;  
class AliAODVertex;

class AliAnalysisTaskKineFilter : public AliAnalysisTaskSE
{
 public:
                                  AliAnalysisTaskKineFilter();
                                  AliAnalysisTaskKineFilter( const char* name );
                                  AliAnalysisTaskKineFilter(const AliAnalysisTaskKineFilter& obj);
    virtual                      ~AliAnalysisTaskKineFilter();
     AliAnalysisTaskKineFilter&   operator=(const AliAnalysisTaskKineFilter& other);
    
    // Implementation of interface methods
    virtual                void   UserCreateOutputObjects();
    virtual                void   Exec( Option_t *option );
    
    // Setters
    virtual                void   SetTrackFilter(AliAnalysisFilter* trackF) { fTrackFilter = trackF; }
    
 private:
                          Int_t   LoopOverSecondaries(TParticle *mother, Int_t& jTracks, Int_t& jVertices, Int_t& nPos, Int_t& nNeg );
                           void   SetChargeAndPID(Int_t pdgCode, AliAODTrack *track);
                           void   SetVertexType(TParticle *part, AliAODVertex *vertex);
                         
              AliAnalysisFilter*  fTrackFilter;          //  Track Filter
                    
    ClassDef( AliAnalysisTaskKineFilter, 1 ); // Analysis task for Kinematic filtering
};
 
#endif
