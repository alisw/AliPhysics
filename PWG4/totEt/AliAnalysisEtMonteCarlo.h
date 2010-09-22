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
   
  AliAnalysisEtMonteCarlo();
  virtual ~AliAnalysisEtMonteCarlo();

    virtual Int_t AnalyseEvent(AliVEvent* event);

    virtual void Init();
    virtual void ResetEventValues();
    virtual void CreateHistograms();

protected:

    virtual bool TrackHitsCalorimeter(TParticle *part, Double_t magField=0.5);

protected:

    Double_t fImpactParameter; // b(fm), for Hijing; 0 otherwise
    Int_t fNcoll; // Ncoll, for Hijing; 1 otherwise
    Int_t fNpart; // Ncoll, for Hijing; 2 otherwise

 private:

    ClassDef(AliAnalysisEtMonteCarlo, 1);
};

#endif // ALIANALYSISETMONTECARLO_H
