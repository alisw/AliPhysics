//_________________________________________________________________________
//  Utility Class for transverse energy studies, charged hadrons
//  Base class for MC analysis
//  - MC output
// 
//Created by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville
//_________________________________________________________________________
#ifndef ALIANALYSISHADETMONTECARLO_H
#define ALIANALYSISHADETMONTECARLO_H

#include "AliAnalysisHadEt.h"
class AliVEvent;

class AliAnalysisHadEtMonteCarlo : public AliAnalysisHadEt
{

public:
   
    AliAnalysisHadEtMonteCarlo() {}
    virtual ~AliAnalysisHadEtMonteCarlo() {}
   
    virtual Int_t AnalyseEvent(AliVEvent* event);
    virtual Int_t AnalyseEvent(AliVEvent* event,AliVEvent* event2);

    //void FillHistograms();
    void CreateHistograms();
    virtual void Init();
    
 private:
    //Declare it private to avoid compilation warning
    AliAnalysisHadEtMonteCarlo & operator = (const AliAnalysisHadEtMonteCarlo & g) ;//cpy assignment
    AliAnalysisHadEtMonteCarlo(const AliAnalysisHadEtMonteCarlo & g) ; // cpy ctor

    ClassDef(AliAnalysisHadEtMonteCarlo, 1);
};

#endif // ALIANALYSISHADETMONTECARLO_H
