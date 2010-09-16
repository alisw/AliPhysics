#ifndef ALIANALYSISETMONTECARLO_H
#define ALIANALYSISETMONTECARLO_H
//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for MC analysis
//  - MC output
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________

#include "AliAnalysisEt.h"
class TParticle;

class AliAnalysisEtMonteCarlo : public AliAnalysisEt
{

public:
   
    AliAnalysisEtMonteCarlo() {}
    virtual ~AliAnalysisEtMonteCarlo() {}

    virtual Int_t AnalyseEvent(AliVEvent* event);

    virtual void Init();

protected:

    virtual bool TrackHitsCalorimeter(TParticle *part, Double_t magField=0.5);

};

#endif // ALIANALYSISETMONTECARLO_H
