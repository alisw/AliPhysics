#ifndef ALIANALYSISETRECONSTRUCTEDEMCAL_H
#define ALIANALYSISETRECONSTRUCTEDEMCAL_H
//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for ESD analysis, for EMCAL
//  - reconstruction output
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________
#include "AliAnalysisEtReconstructed.h"


class AliAnalysisEtReconstructedEmcal : public AliAnalysisEtReconstructed
{

public:
   
    AliAnalysisEtReconstructedEmcal();
    virtual ~AliAnalysisEtReconstructedEmcal();

    virtual void Init();
    
   protected:
      
      virtual bool TrackHitsCalorimeter(AliVParticle *track, Double_t magField);
 private:

      ClassDef(AliAnalysisEtReconstructedEmcal, 1);
};

#endif // ALIANALYSISETRECONSTRUCTEDEMCAL_H
