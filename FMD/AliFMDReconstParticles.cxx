
 //////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Forward Multiplicity Detector have to be reconstructed                   //
//  number of particles in fixed pseudorapidity interval from                //
// fNumOfMinRing  to fNumOfMaxRing                                           //
// and phi interval                                                          //
//  from fNumOfMinSector to fNumOfMaxSector                                  //
 //////////////////////////////////////////////////////////////////////////////
 

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
