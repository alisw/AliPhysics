#ifndef CEPRawEventBuffer_H
#define CEPRawEventBuffer_H

#include "TObject.h"
#include "TObjArray.h"
#include "CEPRawTrackBuffer.h"
#include "CEPRawADCellBuffer.h"
#include "CEPRawV0Buffer.h"
#include "CEPRawCaloBuffer.h"
#include "CEPRawCaloClusterTrack.h"
#include "CEPRawFMDBuffer.h"
#include "AliESDEvent.h"

class CEPRawEventBuffer : public TObject 
{

  private:
    // important for cross checks with CEPEventBuffer object
    /// Event number
    Int_t                       fEventNumber;

    /// summary track information
    Int_t                       fnTracks;       // number of tracks in fCEPTracks
    Int_t                       fnCaloTracks;   // # calo cluster in fCEPRawCaloClusterTracks
    
    /// List of raw tracks
    TObjArray*                  fCEPRawTracks;
    TObjArray*                  fCEPRawCaloClusterTracks;

    /// raw detector information stored in objects
    CEPRawADCellBuffer*         fADCellBuffer;
    CEPRawV0Buffer*             fV0Buffer;
    CEPRawCaloBuffer*           fEMCalBuffer;
    CEPRawCaloBuffer*           fPHOSBuffer;
    CEPRawFMDBuffer*            fFMDBuffer;
    
  public:
                                CEPRawEventBuffer();
                                ~CEPRawEventBuffer();
    /// Modifiers
    void                        Reset();
    void                        SetEventNumber(Int_t evnum)    { fEventNumber = evnum; }

    void                        AddTrack(CEPRawTrackBuffer* trk);
    void                        AddTrack(CEPRawCaloClusterTrack* caloTrk);
    
    /// Setter
    void                        SetEventVariables(AliESDEvent* ESDobj);

    // Accessors
    Int_t                       GetEventNumber()        const { return fEventNumber; }
    Int_t                       GetnTracksTotal()       const { return fnTracks;     }
    Int_t                       GetnCaloClusterTotal()  const { return fnCaloTracks; }

    /// Track accessors
    CEPRawTrackBuffer*          GetTrack(UInt_t ind);
    CEPRawCaloClusterTrack*     GetCaloClusterTrack(UInt_t ind);
    /// Detector object accessors
    CEPRawADCellBuffer*         GetRawADBuffer()        const { return fADCellBuffer; }
    CEPRawV0Buffer*             GetRawV0Buffer()        const { return fV0Buffer;     }
    CEPRawCaloBuffer*           GetRawEMCalBuffer()     const { return fEMCalBuffer;  }
    CEPRawCaloBuffer*           GetRawPHOSBuffer()      const { return fPHOSBuffer;   }
    CEPRawFMDBuffer*            GetRawFMDBuffer()       const { return fFMDBuffer;    }
 
    /// Track removers
    Bool_t                      RemoveTrack(UInt_t ind);
    Bool_t                      RemoveCaloCluster(UInt_t ind);


    ClassDef(CEPRawEventBuffer, 1)     // CEP event buffer
};

#endif
