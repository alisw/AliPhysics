#ifndef ALIMUONTRIGGERCIRCUIT_H
#define ALIMUONTRIGGERCIRCUIT_H

#include "TObjArray.h"
#include <iostream.h>
#include "AliMUONSegmentationTrigger.h"

class AliMUONSegmentationTrigger;
//----------------------------------------------
class AliMUONTriggerCircuit : 
public TObject {
 public: 
  AliMUONTriggerCircuit();  
  virtual ~AliMUONTriggerCircuit(){;} 
  // copy constructor
  AliMUONTriggerCircuit(const AliMUONTriggerCircuit& AliMUONTriggerCircuit); 
  // assignment operator
  AliMUONTriggerCircuit& operator=(const AliMUONTriggerCircuit& AliMUONTriggerCircuit); 

  // initializations
  void Init(Int_t iCircuit);    
  
  // get calculated pt
  Float_t PtCal(Int_t istripX, Int_t idev, Int_t istripY);

  //--- methods which return member data related info
  Int_t   GetIdCircuit();
  Int_t   GetIdModule();
  Int_t   GetNstripX();
  Int_t   GetNstripY();
  Int_t   GetPosCircuit();
  Int_t   GetIdCircuitD();
  Int_t   GetICircuitD();
  Int_t   GetIdCircuitU();
  Int_t   GetICircuitU();
  Int_t   GetX2m();
  Int_t   GetX2ud();
  void    GetOrMud(Int_t orMud[2]);
  Int_t   GetXcode(Int_t chamber, Int_t istrip);
  Int_t   GetYcode(Int_t chamber, Int_t istrip);
  Float_t GetY11Pos(Int_t istrip);
  Float_t GetY21Pos(Int_t istrip);
  Float_t GetX11Pos(Int_t istrip);
 
  //  Get reference to segmentation model
  virtual AliSegmentation*  SegmentationModel(Int_t isec) {
      return (AliSegmentation *) (*fSegmentation)[isec-1];
  }

 protected:
  TObjArray            *fSegmentation;    // pointer to segmentation

 private:
  Int_t CircuitNumber(Int_t idCircuit);
  Int_t ModuleNumber(Int_t idModule); 
  Int_t Module(Int_t idCircuit);
  Int_t Position(Int_t idCircuit);
  void LoadX2();
  void LoadXCode();
  void LoadYCode();
  void LoadYPos();
  void LoadXPos();
  
  ClassDef(AliMUONTriggerCircuit,1) // Trigger Circuit class
    
 private:    
  Int_t fIdCircuit;            // circuit Id number
  Int_t fX2m;                  // internal info needed by TriggerDecision
  Int_t fX2ud;                 // internal info needed by TriggerDecision
  Int_t fOrMud[2];             // internal info needed by TriggerDecision
  Int_t fXcode[4][32];         // code of X strips
  Int_t fYcode[4][32];         // code of Y strips 
  Float_t fXpos11[16];         // X position of Y strips in MC11
  Float_t fYpos11[31];         // Y position of X strips in MC11
  Float_t fYpos21[63];         // Y position of X strips in MC21

};
#endif



