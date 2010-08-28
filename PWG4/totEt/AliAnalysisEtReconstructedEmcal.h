#ifndef ALIANALYSISETRECONSTRUCTEDEMCAL_H
#define ALIANALYSISETRECONSTRUCTEDEMCAL_H

#include "AliAnalysisEtReconstructed.h"


class AliAnalysisEtReconstructedEmcal : public AliAnalysisEtReconstructed
{

public:
   
    AliAnalysisEtReconstructedEmcal();

    virtual void Init();
    
   protected:
      
      virtual bool TrackHitsCalorimeter(AliVParticle *track, Double_t magField);
};

#endif // ALIANALYSISETRECONSTRUCTEDEMCAL_H
