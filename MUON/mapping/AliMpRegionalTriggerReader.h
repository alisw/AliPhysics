/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $MpId: $ 

/// \ingroup mptrigger
/// \class AliMpRegionalTriggerReader
/// \brief The class to read regional trigger crate file
///
/// \author Ch. Finck, Subatech Nantes

#ifndef ALI_MP_REGIONAL_TRIGGER_READER_H
#define ALI_MP_REGIONAL_TRIGGER_READER_H

#include <TObject.h>
#include <TList.h>

class AliMpRegionalTriggerReader : public  TObject{

  public:
    AliMpRegionalTriggerReader();
    virtual ~AliMpRegionalTriggerReader();
    
    // methods
    static Int_t ReadData(TList& list, istream& in);
    

  private:
    /// Not implemented
    AliMpRegionalTriggerReader(const AliMpRegionalTriggerReader& rhs);
    /// Not implemented
    AliMpRegionalTriggerReader& operator=(const AliMpRegionalTriggerReader& rhs);

 
  ClassDef(AliMpRegionalTriggerReader,0) // Regional trigger crate reader
};


#endif //ALI_MP_REGIONAL_CRATE_READER_H














