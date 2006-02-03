#ifndef ALIMUONTRIGGERCRATE_H
#define ALIMUONTRIGGERCRATE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*-- Author: Rachid Guernane (LPCCFd)

#include <TNamed.h>
#include <TObjArray.h>

class AliMUONTriggerBoard;

class AliMUONTriggerCrate : public TNamed
{
   public:

      AliMUONTriggerCrate();
      AliMUONTriggerCrate(const AliMUONTriggerCrate &entry);
      AliMUONTriggerCrate(const char *name, Int_t n = 17); // 16 + 1
      virtual ~AliMUONTriggerCrate();

//    CRATE CONFIG FROM ASCII FILE
      virtual void SetDataSource(TString SourceFile) {fSourceFileName = SourceFile;}

      virtual void AddBoard(AliMUONTriggerBoard *board, Int_t i);

      virtual TObjArray* Boards() {return fBoards;}

      AliMUONTriggerCrate& operator=(const AliMUONTriggerCrate &rhs);

   protected:

      void Copy(TObject&) const;

   private:

      Int_t     fNslots;
      Int_t     fNboards;

      TObjArray *fBoards;
      TString   fSourceFileName;

   ClassDef(AliMUONTriggerCrate,1)
};
#endif
