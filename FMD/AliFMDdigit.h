#ifndef ALIFMDDIGIT_H
#define ALIFMDDIGIT_H

////////////////////////////////////////////////
//  Digits classes for set:FMD                //
////////////////////////////////////////////////
#include <TClonesArray.h>
#include <TObject.h>


//___________________________________________
class AliFMDdigit: public TObject  {
  
private:
  Int_t fNumOfDet ;       //Number of FMD disk
  Int_t fNumOfSector ;    //Number of sector
  Int_t fNumOfRing;       //Number of ring
  Int_t fNelectrons;      // real charge
  Int_t fADC;             // ADC signal for this charge
 
public:
  AliFMDdigit() {
    // constructor
    fNumOfDet=fNumOfSector=fNumOfRing=fNelectrons=fADC=0;
  }
  AliFMDdigit(Int_t *digits);
  virtual ~AliFMDdigit() {
    // destructor
  }
  Int_t Volume();
  Int_t NumberOfSector();
  Int_t NumberOfRing();
  Int_t Charge();
  Int_t ADCsignal();
  
  ClassDef(AliFMDdigit,2)     // Real data digit object for set:ITS
};
inline Int_t AliFMDdigit::Volume(){return fNumOfDet;}
inline Int_t AliFMDdigit::NumberOfSector() {return fNumOfSector;} 
inline Int_t AliFMDdigit::NumberOfRing() {return fNumOfRing;}
inline Int_t AliFMDdigit::Charge() {return fNelectrons;}
inline Int_t AliFMDdigit::ADCsignal() {return fADC;}


#endif
