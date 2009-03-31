/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $MpId: $ 

/// \ingroup management
/// \class AliMpFrtCrocusConstants
/// \brief The class defines the properties of CROCUS FRT
///
/// \author Ch. Finck, Subatech Nantes

#ifndef ALI_MP_FRT_CROCUS_CONSTANTS_H
#define ALI_MP_FRT_CROCUS_CONSTANTS_H

#include <TObject.h>
#include <TString.h>

#include "AliMpArrayI.h"
#include "AliMpEncodePair.h"

class AliMpFrtCrocusConstants : public  TObject {

  public:
    AliMpFrtCrocusConstants();
    virtual ~AliMpFrtCrocusConstants();

   // static methods
    static Int_t GetGlobalFrtID(Int_t localID, Int_t ddlID);
    static Int_t GetLocalFrtID(Int_t globalID, Int_t ddlID);
    
    // get methods
    static Int_t  GetNofDsps();
    static Int_t  GetNofBusPatches();    
    static MpPair_t GetLinkPortId(Int_t index);
    
    // return VME top address
    static UInt_t GetTopAddress(Int_t id);
    static Int_t  GetIdFromTopAddress(UInt_t add);
    
    // return VME bottom address
    static UInt_t GetBotAddress(Int_t id) ;
    static Int_t  GetIdFromBotAddress(UInt_t add);


  private:
    /// Not implemented
    AliMpFrtCrocusConstants(const AliMpFrtCrocusConstants& rhs);
    /// Not implemented
    AliMpFrtCrocusConstants& operator=(const AliMpFrtCrocusConstants& rhs);

   // static data members	
    static const Int_t  fgkOffset;         ///< Offset for conversion global/local ID  
    static const Int_t  fgkLinkPorts[10];  ///< Link port Ids connected to this crocus
    static const Int_t  fgkNofDsps;        ///< Number of Dsps  connected to this crocus
    static const Int_t  fgkNofBusPatches;  ///< Number of Dsps  connected to this crocus
    static const UInt_t fgkBaseAddress;    ///< VME base address for FRT crocus
    static const UInt_t fgkAddressOffset;  ///< VME address offset for FRT crocus

  ClassDef(AliMpFrtCrocusConstants,1)  // The class collectiong electronics properties of CROCUS FRT
};

#endif //ALI_FRT_CROCUS_CONSTANTS_H














