#ifndef ALIMUONTRIGGERELECTRONICS_H
#define ALIMUONTRIGGERELECTRONICS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*-- Author: Rachid Guernane (LPCCFd)

#ifndef ROOT_TTask
#  include "TTask.h"
#endif

#ifndef ROOT_TArrayI
#  include "TArrayI.h"
#endif

#ifndef ROOT_TString
#  include "TString.h"
#endif

class AliMUONTriggerCrate;
class AliMUONData;
class AliMUONGlobalTriggerBoard;
class TClonesArray;

class AliMUONTriggerElectronics : public TTask
{
   public:
      AliMUONTriggerElectronics(AliMUONData* Data = 0);
      virtual ~AliMUONTriggerElectronics();

      virtual void Exec(Option_t*);
      
//    CRATE CONFIG FROM ASCII FILE
      virtual void SetDataSource(TString SourceFile = "$ALICE_ROOT/MUON/data/CRATE.TXT") 
      {fSourceFileName = SourceFile;}

      virtual void Factory();
      virtual void AddCrate(char *name);

      virtual AliMUONTriggerCrate* Crate(char *name);

      virtual void Feed();
      virtual void Feed(UShort_t pattern[2][4]);

      virtual void BoardName(Int_t ix, Int_t iy, char *name);

      virtual void Reset();

      virtual void Scan(Option_t *option);

      virtual void LocalResponse();
      virtual void RegionalResponse();
      virtual void GlobalResponse();

      virtual void BuildName(Int_t icirc, char name[20]);

      virtual void DumpOS();

      virtual void Digits2Trigger();
      virtual void Trigger();
      virtual void ClearDigitNumbers();
      virtual void DigitFiredCircuit(Int_t circuit, Int_t cathode, Int_t chamber, Int_t digit);

//      virtual AliMUONData* GetMUONData() {return fMUONData;}

   private:
      TString                    fSourceFileName;
      TClonesArray              *fCrates;
      AliMUONGlobalTriggerBoard *fGlobalTriggerBoard;
      Int_t                      fNCrates;
      static const Int_t         fgkNCrates;
      UShort_t                   fLocal[16][16];     // 16 crates of 16 LB
      UShort_t                   fRegional[16];
      UShort_t                   fGlobal;
      AliMUONData               *fMUONData;          //! Data container for MUON subsystem 
      TArrayI                    fDigitNumbers[234]; //! The digit number that fired a circuit.

   ClassDef(AliMUONTriggerElectronics,1)
};
#endif
