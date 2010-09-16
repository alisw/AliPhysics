#ifndef ALIANALYSISETRECONSTRUCTED_H
#define ALIANALYSISETRECONSTRUCTED_H
//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for ESD analysis
//  - reconstruction output
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________

#include "AliAnalysisEt.h"

class AliVParticle;

class AliAnalysisEtReconstructed : public AliAnalysisEt
{

public:
   
    AliAnalysisEtReconstructed();
    virtual ~AliAnalysisEtReconstructed();

    virtual Int_t AnalyseEvent(AliVEvent* event);

    virtual void Init();
    
protected:

    bool CheckGoodVertex(AliVParticle *track);
    virtual bool TrackHitsCalorimeter(AliVParticle *track, Double_t magField);

    Double_t fTrackDistanceCut; // cut on track distance    
    Char_t fClusterType; // selection on cluster type
    
};

#endif // ALIANALYSISETRECONSTRUCTED_H
