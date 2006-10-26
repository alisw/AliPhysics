#ifndef ALIMUONDDLTRACKER_H
#define ALIMUONDDLTRACKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup raw
/// \class AliMUONDDLTracker
/// \brief MUON DDL tracker
///
/// \author Christian Finck

#include <TObject.h>
#include <TClonesArray.h>

class AliMUONBusStruct;
class AliMUONDspHeader;
class AliMUONBlockHeader;

class AliMUONDDLTracker : public TObject {

public:
   AliMUONDDLTracker();
   virtual ~AliMUONDDLTracker();
 
   void    AddBusPatch(const AliMUONBusStruct& busPatch, Int_t iBlock, Int_t iDsp);
   void    AddDspHeader(const AliMUONDspHeader& dspHeader, Int_t iBlock);
   void    AddBlkHeader(const AliMUONBlockHeader& blkHeader);

   // get TClonesArray
   TClonesArray*  GetBlkHeaderArray() const {return fBlkHeaderArray;}

   // get entries
   Int_t GetBlkHeaderEntries() const {return fBlkHeaderArray->GetEntriesFast();}
 
   // get entry
   AliMUONBlockHeader* GetBlkHeaderEntry(Int_t i) const {
     return (AliMUONBlockHeader*)fBlkHeaderArray->At(i);}

   // clear
   void Clear(Option_t* opt);


 private:

   TClonesArray* fBlkHeaderArray;  ///< array of block header
 
   AliMUONDDLTracker(const AliMUONDDLTracker& event);
   AliMUONDDLTracker& operator=(const AliMUONDDLTracker& event);
 
   ClassDef(AliMUONDDLTracker,1)  // MUON DDL Tracker
};
#endif
