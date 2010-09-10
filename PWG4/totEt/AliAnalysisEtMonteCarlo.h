#ifndef ALIANALYSISETMONTECARLO_H
#define ALIANALYSISETMONTECARLO_H

#include "AliAnalysisEt.h"
#include "TParticle.h"

class AliAnalysisEtMonteCarlo : public AliAnalysisEt
{

public:
   
    AliAnalysisEtMonteCarlo() {}
   
    virtual Int_t AnalyseEvent(AliVEvent* event);

    virtual void Init();

protected:

    virtual bool TrackHitsCalorimeter(TParticle *part, Double_t magField=0.5);

};

#endif // ALIANALYSISETMONTECARLO_H
