#ifndef ALIANALYSISETMONTECARLOEMCAL_H
#define ALIANALYSISETMONTECARLOEMCAL_H
//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for MC analysis, for EMCAL
//  - MC output
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________
#include "AliAnalysisEtMonteCarlo.h"


class AliAnalysisEtMonteCarloEmcal : public AliAnalysisEtMonteCarlo
{

public:
   
    AliAnalysisEtMonteCarloEmcal();
    virtual ~AliAnalysisEtMonteCarloEmcal();

    virtual void Init();
    void CreateHistograms();

 private:

    ClassDef(AliAnalysisEtMonteCarloEmcal, 1);
};

#endif // ALIANALYSISETRECONSTRUCTEDEMCAL_H
