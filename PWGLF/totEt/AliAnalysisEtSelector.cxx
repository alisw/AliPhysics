#include "AliAnalysisEtSelector.h"
#include "AliAnalysisEtCuts.h"
#include <AliESDCaloCluster.h>

ClassImp(AliAnalysisEtSelector)

AliAnalysisEtSelector::AliAnalysisEtSelector(AliAnalysisEtCuts *cuts) :
fEvent(0)
,fCuts(cuts)
,fRunNumber(0)
{

}

AliAnalysisEtSelector::~AliAnalysisEtSelector()
{

}

