#include "AliEmcalContainerUtils.h"

// TEMP
// Make template explicit for now
#include <AliVCluster.h>
#include <AliVParticle.h>

#include "AliClusterContainer.h"
#include "AliParticleContainer.h"
template class AliEmcalContainerUtils<TClonesArray, AliVCluster>;
template class AliEmcalContainerUtils<TClonesArray, AliVParticle>;
template class AliEmcalContainerUtils<AliClusterContainer, AliVCluster>;
template class AliEmcalContainerUtils<AliParticleContainer, AliVParticle>;
// ENDTEMP

