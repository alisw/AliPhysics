#ifndef ALIANALYSISETRECONSTRUCTEDPHOS_H
#define ALIANALYSISETRECONSTRUCTEDPHOS_H

#include "AliAnalysisEtReconstructed.h"


class AliAnalysisEtReconstructedPhos : public AliAnalysisEtReconstructed
{

public:
   
    AliAnalysisEtReconstructedPhos();

    virtual void Init();
    
   protected:
      
      virtual bool TrackHitsCalorimeter(AliVParticle *track, Double_t magField);
};

#endif // ALIANALYSISETRECONSTRUCTEDPHOS_H
