#ifndef ALIMUONDDLTRIGGER_H
#define ALIMUONDDLTRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup raw
/// \class AliMUONDDLTrigger
/// \brief MUON DDL Trigger
///
/// \author Christian Finck

#include <TObject.h>
#include <TClonesArray.h>

class AliMUONLocalStruct;
class AliMUONRegHeader;
class AliMUONDarcHeader;

class AliMUONDDLTrigger : public TObject {
 
public:
   AliMUONDDLTrigger();
   virtual ~AliMUONDDLTrigger();


   void    AddLocStruct(const AliMUONLocalStruct& loc, Int_t iReg);
   void    AddRegHeader(const AliMUONRegHeader& regHeader);

   // get AliMUONDarcHeader
   AliMUONDarcHeader*  GetDarcHeader() const {return fDarcHeader;}

 private:

   AliMUONDarcHeader* fDarcHeader;  ///< pointer of darc header
   

   ClassDef(AliMUONDDLTrigger,1)  // MUON DDL Trigger
};
#endif
