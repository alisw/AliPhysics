#include "AliFMDReconstParticles.h"

ClassImp(AliFMDReconstParticles)

AliFMDReconstParticles::AliFMDReconstParticles()
{
     fNumOfDet=0;
     fNumOfMinSector=0;
     fNumOfMaxSector=0;
     fNumOfMinRing=0;
     fNumOfMaxRing=0;
     fNumOfReconstParticles=0;
}
AliFMDReconstParticles::AliFMDReconstParticles(Int_t *RecParticles)
{
     fNumOfDet=RecParticles[0];
     fNumOfMinSector=RecParticles[1];
     fNumOfMaxSector=RecParticles[2];
     fNumOfMinRing=RecParticles[3];
     fNumOfMaxRing=RecParticles[4];
     fNumOfReconstParticles=RecParticles[5];
}
Int_t AliFMDReconstParticles::GetVolume(){return fNumOfDet;}
Int_t AliFMDReconstParticles::GetNumberOfMinSector() {return fNumOfMinSector;}
Int_t AliFMDReconstParticles::GetNumberOfMaxSector() {return fNumOfMaxSector;}
Int_t AliFMDReconstParticles::GetNumberOfMinRing() {return fNumOfMinRing;}
Int_t AliFMDReconstParticles::GetNumberOfMaxRing() {return fNumOfMaxRing;}
Int_t AliFMDReconstParticles::GetNumberOfReconstParticles() {return fNumOfReconstParticles;}

