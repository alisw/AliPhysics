#ifndef ALIANALYSISETRECONSTRUCTEDPHOS_H
#define ALIANALYSISETRECONSTRUCTEDPHOS_H
//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for ESD analysis, for PHOS
//  - reconstruction output
//  implementation file
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________

#include "AliAnalysisEtReconstructed.h"


class AliAnalysisEtReconstructedPhos : public AliAnalysisEtReconstructed
{

public:
   
    AliAnalysisEtReconstructedPhos();
    virtual ~AliAnalysisEtReconstructedPhos();

    virtual void Init();
    
   protected:
      
      virtual bool TrackHitsCalorimeter(AliVParticle *track, Double_t magField);
 private:

      ClassDef(AliAnalysisEtReconstructedPhos, 1);
};

#endif // ALIANALYSISETRECONSTRUCTEDPHOS_H
