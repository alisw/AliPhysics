#ifndef ALIMUONTRIGGERCRATE_H
#define ALIMUONTRIGGERCRATE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup trigger
/// \class AliMUONTriggerCrate
/// \brief Trigger Crate
///
//  Author Rachid Guernane (LPCCFd)

#include <TNamed.h>

class AliMUONTriggerBoard;
class TObjArray;

class AliMUONTriggerCrate : public TNamed
{
   public:
      AliMUONTriggerCrate();
      AliMUONTriggerCrate(const char *name, Int_t n = 17); // 16 + 1
      virtual ~AliMUONTriggerCrate();

      /// Crate config from ascii file
      virtual void SetDataSource(TString SourceFile) {fSourceFileName = SourceFile;}

      virtual void AddBoard(AliMUONTriggerBoard *board, Int_t i);

      /// Return pointer to board objects
      virtual TObjArray* Boards() {return fBoards;}



   private:
      /// Not implemented
      AliMUONTriggerCrate(const AliMUONTriggerCrate &entry);
      /// Not implemented
      AliMUONTriggerCrate& operator=(const AliMUONTriggerCrate &rhs);

      Int_t     fNslots;          ///< Number of slots
      Int_t     fNboards;         ///< Number of boards

      TObjArray *fBoards;         ///< Pointer to board objects
      TString   fSourceFileName;  ///< Source file

   ClassDef(AliMUONTriggerCrate,1) //Trigger Crate
};
#endif
