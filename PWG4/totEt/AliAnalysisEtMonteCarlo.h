#ifndef ALIANALYSISETMONTECARLO_H
#define ALIANALYSISETMONTECARLO_H

#include "AliAnalysisEt.h"


class AliAnalysisEtMonteCarlo : public AliAnalysisEt
{

public:
   
    AliAnalysisEtMonteCarlo() {}
   
    virtual Int_t AnalyseEvent(AliVEvent* event);

    virtual void Init();
    
};

#endif // ALIANALYSISETMONTECARLO_H
