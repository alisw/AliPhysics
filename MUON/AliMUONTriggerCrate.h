#ifndef ALIMUONTRIGGERCRATE_H
#define ALIMUONTRIGGERCRATE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup sim
/// \class AliMUONTriggerCrate
/// \brief Trigger Crate
///
/// \author Rachid Guernane (LPCCFd)

#include <TNamed.h>

class AliMUONTriggerBoard;
class TObjArray;

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

      Int_t     fNslots;          ///< NUMBER OF SLOTS
      Int_t     fNboards;         ///< NUMBER OF BOARDS

      TObjArray *fBoards;         ///< POINTER TO BOARD OBJECTS
      TString   fSourceFileName;  ///< SOURCE FILE

   ClassDef(AliMUONTriggerCrate,1)
};
#endif
