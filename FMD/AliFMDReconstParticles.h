#ifndef AliFMDReconstParticles_h
#define AliFMDReconstParticles_h

#include <TObject.h>
class TClonesArray; 

class AliFMDReconstParticles: public TObject
{
  /*
 Reconstracted Particles Class:
 has number of reconstructed particles in 
 sectors from NumOfMinSector to NumberOfMaxSector()
 rings from NumOfMinRing to NumOfMaxRing for each FMDvolume
  */


 public:
  AliFMDReconstParticles();
  AliFMDReconstParticles (Int_t *RecParticles);
  virtual ~AliFMDReconstParticles(){};
  Int_t GetVolume() const  {return fNumOfDet;}  //number of volume 
  Int_t GetNumberOfMinSector() const  {return fNumOfMinSector;} //number of min sector
  Int_t GetNumberOfMaxSector() const {return fNumOfMaxSector;}  // number of max sector
  Int_t GetNumberOfMinRing() const {return fNumOfMinRing;} //number of min ring
  Int_t GetNumberOfMaxRing() const {return fNumOfMaxRing;}  //number of max ring 
  Int_t GetNumberOfReconstParticles()const {return fNumOfReconstParticles;} //number of reconstructed particles
 private:
  Int_t fNumOfDet;                       //Number of FMD disk;
  Int_t fNumOfMinSector;                    //Number of min sector 
  Int_t fNumOfMaxSector;                    //Number of max sector
  Int_t fNumOfMinRing;                   //Number of min ring
  Int_t fNumOfMaxRing;                      //Number of max ring
  Int_t fNumOfReconstParticles;          //Number of reconstructed particles

  ClassDef(AliFMDReconstParticles,2)
};
#endif
