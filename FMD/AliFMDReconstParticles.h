#ifndef AliFMDReconstParticles_h
#define AliFMDReconstParticles_h

#include <TClonesArray.h> 
#include <TObject.h>
class AliFMDReconstParticles: public TObject
{
 //Reconstracted Particles Class
 private:
  Int_t fNumOfDet;                       //Number of FMD disk;
  Int_t fNumOfMinSector;                    //Number of min sector 
  Int_t fNumOfMaxSector;                    //Number of max sector
  Int_t fNumOfMinRing;                   //Number of min ring
  Int_t fNumOfMaxRing;                      //Number of max ring
  Int_t fNumOfReconstParticles;          //Number of reconstructed particles

 public:
  AliFMDReconstParticles();
  AliFMDReconstParticles (Int_t *RecParticles);
  virtual ~AliFMDReconstParticles(){};
  Int_t GetVolume();                    //Return the number of volume 
  Int_t GetNumberOfMinSector();           //Return the number of min sector
  Int_t GetNumberOfMaxSector();           //Return the number of max sector
  Int_t GetNumberOfMinRing();             //Return the number of min ring
  Int_t GetNumberOfMaxRing();              //Return the number of max ring 
  Int_t GetNumberOfReconstParticles();  //Returnthe the number of reconstructed particles
  ClassDef(AliFMDReconstParticles,2)
};
#endif
