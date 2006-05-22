#ifndef ALIMUONTRIGGERCIRCUIT_H
#define ALIMUONTRIGGERCIRCUIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004
//
/// \ingroup base
/// \class AliMUONTriggerCircuit
/// \brief MUON Trigger circuit
///
/// \author Philippe Crochet (LPCCFd)

#include <TObject.h>
#include <TObjArray.h>

class TObjArray;

class AliMUONTriggerCircuit : public TObject 
{
 public: 
  AliMUONTriggerCircuit();  
  virtual ~AliMUONTriggerCircuit(){;} 

  // initializations
  void Init(Int_t iCircuit);    
  
  // get calculated pt
  Float_t PtCal(Int_t istripX, Int_t idev, Int_t istripY);

  //--- methods which return member data related info
  Int_t   GetIdCircuit() const;
  Int_t   GetIdModule() const;
  Int_t   GetNstripX() const;
  Int_t   GetNstripY() const;
  Int_t   GetPosCircuit() const;
  Int_t   GetIdCircuitD() const;
  Int_t   GetICircuitD() const;
  Int_t   GetIdCircuitU() const;
  Int_t   GetICircuitU() const;
  Int_t   GetX2m() const;
  Int_t   GetX2ud() const;
  void    GetOrMud(Int_t orMud[2]) const;
  Int_t   GetXcode(Int_t chamber, Int_t istrip) const;
  Int_t   GetYcode(Int_t chamber, Int_t istrip) const;
  Float_t GetY11Pos(Int_t istrip) const;
  Float_t GetY21Pos(Int_t istrip) const;
  Float_t GetX11Pos(Int_t istrip) const;
  Int_t   DetElemId(Int_t ichamber, Int_t idModule);

  void Print(Option_t* opt="") const;
  
 protected:
  // copy constructor
  AliMUONTriggerCircuit(const AliMUONTriggerCircuit& AliMUONTriggerCircuit); 
  // assignment operator
  AliMUONTriggerCircuit& operator=(const AliMUONTriggerCircuit& AliMUONTriggerCircuit); 

 private:
  Int_t CircuitNumber(Int_t idCircuit) const;
  Int_t ModuleNumber(Int_t idModule) const; 
  Int_t Module(Int_t idCircuit) const;
  Int_t Position(Int_t idCircuit) const;
  void LoadX2();
  void LoadXCode();
  void LoadYCode();
  void LoadYPos2();
  void LoadXPos2();
   
  Int_t fIdCircuit;            ///< circuit Id number
  Int_t fX2m;                  ///< internal info needed by TriggerDecision
  Int_t fX2ud;                 ///< internal info needed by TriggerDecision
  Int_t fOrMud[2];             ///< internal info needed by TriggerDecision
  Int_t fXcode[4][32];         ///< code of X strips
  Int_t fYcode[4][32];         ///< code of Y strips 
  Float_t fXpos11[16];         ///< X position of Y strips in MC11
  Float_t fYpos11[31];         ///< Y position of X strips in MC11
  Float_t fYpos21[63];         ///< Y position of X strips in MC21

  ClassDef(AliMUONTriggerCircuit,1) // Trigger Circuit class
};
#endif



