#ifndef ALIFMDDIGITIZERALLA_H
#define ALIFMDDIGITIZERALLA_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//
// Resurection of Alla's old digitizer
// This is to investigate the changes we've seen in the reconstructed 
// particle multiplicity 
//
#include <AliDigitizer.h>
#include <AliRunDigitizer.h>
class TClonesArray;

class AliFMDDigitizerAlla : public AliDigitizer 
{
public:
  AliFMDDigitizerAlla();
  AliFMDDigitizerAlla(AliRunDigitizer * manager);
  virtual ~AliFMDDigitizerAlla();
  virtual Bool_t Init();
  // Do the main work
  void  Exec(Option_t* option=0) ;
  Int_t PutNoise(Int_t charge) {return (Int_t)(gRandom->Gaus(charge,500));}
  TClonesArray *Digits() const {return fDigits;}
  TClonesArray *Hits()   const {return fHits;}
  enum {kBgTag = -1};
private:
  TClonesArray *fDigits;               // ! array with digits
  TClonesArray *fHits;                 // List of hits
  AliRunDigitizer* GetManager(){return fManager;}
  ClassDef(AliFMDDigitizerAlla,0)
};    
#endif
//
// Local Variables:
//  mode: C++
// End:
//




