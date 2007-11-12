#ifndef ALI_MP_CDB_H
#define ALI_MP_CDB_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$ 

/// \ingroup management
/// \class AliMpCDB
/// \brief Manager class for mapping CDB IO
///
/// \author Ivana Hrivnacova, IPN Orsay

#ifndef ROOT_TObject
#  include <TObject.h>
#endif

class AliMpSegmentation;
class AliMpDDLStore;

class AliMpCDB : public  TObject {

  public:
    // static methods
    static Bool_t LoadMpSegmentation(Bool_t warn = false);
    static Bool_t LoadDDLStore(Bool_t warn = false);

    static Bool_t WriteMpSegmentation(Bool_t readData = true);
    static Bool_t WriteDDLStore(Bool_t readData= true);
     
  private:
    /// Not implemented
    AliMpCDB();
    /// Not implemented
    AliMpCDB(const AliMpCDB& rhs);
    /// Not implemented
    AliMpCDB& operator=(const AliMpCDB& rhs);
    
  ClassDef(AliMpCDB,0)  // The factory for building mapping segmentations
};

#endif //ALI_MP_CDB_H















