#include "AliFMDReconstParticles.h"

ClassImp(AliFMDReconstParticles)

AliFMDReconstParticles::AliFMDReconstParticles()
{
     fNumOfDet=0;
     fNumOfSector=0;
     fNumOfRing=0;
     fNumOfReconstParticles=0;
}
AliFMDReconstParticles::AliFMDReconstParticles(Int_t *RecParticles)
{
     fNumOfDet=RecParticles[0];
     fNumOfSector=RecParticles[1];
     fNumOfRing=RecParticles[2];
     fNumOfReconstParticles=RecParticles[3];
}
Int_t AliFMDReconstParticles::GetVolume(){return fNumOfDet;}
Int_t AliFMDReconstParticles::GetNumberOfSector() {return fNumOfSector;}
Int_t AliFMDReconstParticles::GetNumberOfRing() {return fNumOfRing;}
Int_t AliFMDReconstParticles::GetNumberOfReconstParticles() {return fNumOfReconstParticles;}

