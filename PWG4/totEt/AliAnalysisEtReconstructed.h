#ifndef ALIANALYSISETRECONSTRUCTED_H
#define ALIANALYSISETRECONSTRUCTED_H

#include "AliAnalysisEt.h"

class AliVParticle;

class AliAnalysisEtReconstructed : public AliAnalysisEt
{

public:
   
    AliAnalysisEtReconstructed();
   
    virtual Int_t AnalyseEvent(AliVEvent* event);

    virtual void Init();
    
protected:

    bool CheckGoodVertex(AliVParticle *track);
    virtual bool TrackHitsCalorimeter(AliVParticle *track, Double_t magField);

    Int_t fNTpcClustersCut;
    Int_t fNItsClustersCut;
   
    Double_t fTrackDistanceCut;
    
    Char_t fClusterType;
    
};

#endif // ALIANALYSISETRECONSTRUCTED_H
