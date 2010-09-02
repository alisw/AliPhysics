//Create by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville
#ifndef ALIANALYSISHADETMONTECARLO_H
#define ALIANALYSISHADETMONTECARLO_H

#include "AliAnalysisHadEt.h"
#include "AliESDtrackCuts.h"
#include <iostream>

class AliAnalysisHadEtMonteCarlo : public AliAnalysisHadEt
{

public:
   
    AliAnalysisHadEtMonteCarlo() {}
   
    virtual Int_t AnalyseEvent(AliVEvent* event);
    virtual Int_t AnalyseEvent(AliVEvent* event,AliVEvent* event2);

    //void FillHistograms();
    void CreateHistograms();
    virtual void Init();
    
 private:
    //AliESDtrackCuts* esdtrackCutsITSTPC;
};

#endif // ALIANALYSISHADETMONTECARLO_H
