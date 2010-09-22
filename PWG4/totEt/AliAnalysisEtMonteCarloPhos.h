#ifndef ALIANALYSISETMONTECARLOPHOS_H
#define ALIANALYSISETMONTECARLOPHOS_H
//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for MC analysis, for PHOS
//  - MC output
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________

#include "AliAnalysisEtMonteCarlo.h"


class AliAnalysisEtMonteCarloPhos : public AliAnalysisEtMonteCarlo
{

public:
   
    AliAnalysisEtMonteCarloPhos();
    virtual ~AliAnalysisEtMonteCarloPhos();

    virtual void Init();

 private:

   ClassDef(AliAnalysisEtMonteCarloPhos, 1); 
};

#endif // ALIANALYSISETRECONSTRUCTEDPHOS_H
