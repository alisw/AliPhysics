#ifndef AliFMDReconstParticles_h
#define AliFMDReconstParticles_h

#include <TClonesArray.h> 
#include <TObject.h>
class AliFMDReconstParticles: public TObject
{
 //Reconstracted Particles Class
 private:
  Int_t fNumOfDet;                       //Number of FMD disk;
  Int_t fNumOfSector;                    //Number of sector
  Int_t fNumOfRing;                      //Number of ring
  Int_t fNumOfReconstParticles;          //Number of reconstructed particles

 public:
  AliFMDReconstParticles();
  AliFMDReconstParticles (Int_t *RecParticles);
  virtual ~AliFMDReconstParticles(){};
  Int_t GetVolume();                    //Return the number of volume 
  Int_t GetNumberOfSector();            //Return the number of sector
  Int_t GetNumberOfRing();              //Return the number of ring 
  Int_t GetNumberOfReconstParticles();  //Returnthe number of reconstructed particles
  ClassDef(AliFMDReconstParticles,1)
};
#endif
