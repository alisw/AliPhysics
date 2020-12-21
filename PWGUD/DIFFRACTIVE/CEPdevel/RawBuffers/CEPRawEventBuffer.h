#ifndef CEPRawEventBuffer_H
#define CEPRawEventBuffer_H

#include "TObject.h"
#include "TObjArray.h"
#include "CEPRawTrackBuffer.h"
#include "CEPRawADBuffer.h"
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
    UInt_t                      fEventNumber;

    /// summary track information
    UInt_t                      fnTracks;       // number of tracks in fCEPTracks
    UInt_t                      fnCaloTracks;   // # calo cluster in fCEPRawCaloClusterTracks
    
    /// Detector cell sum of amplitude, times, etc.
    // AD
    Float_t                     fADTotalMult;
    Float_t                     fADTotalTime;
    Float_t                     fADTotalCharge;
    // FMD
    Float_t                     fFMDTotalMult;
    // V0 
    Float_t                     fV0TotalMult;
    Float_t                     fV0TotalTime;
    Float_t                     fV0TotalCharge;
    Float_t                     fV0TotalSigWidth;
    // Calo buffers
    Float_t                     fEMCTotalAmplitude;
    Float_t                     fEMCTotalTime;
    Float_t                     fPHOSTotalAmplitude;
    Float_t                     fPHOSTotalTime;

    /// List of raw tracks
    TObjArray*                  fCEPRawTracks;
    TObjArray*                  fCEPRawCaloClusterTracks;

    /// raw detector information stored in objects
    CEPRawADBuffer*             fADCellBuffer;
    CEPRawV0Buffer*             fV0Buffer;
    CEPRawCaloBuffer*           fEMCalBuffer;
    CEPRawCaloBuffer*           fPHOSBuffer;
    CEPRawFMDBuffer*            fFMDBuffer;
    
  public:
                                CEPRawEventBuffer();
    virtual                     ~CEPRawEventBuffer();
    
    void                        AddTrack(CEPRawTrackBuffer* trk);
    void                        AddCaloTrack(CEPRawCaloClusterTrack* caloTrk);
    
    /// Setter
    void                        SetEventNumber(Int_t evnum)         { fEventNumber = evnum; }
    // AD total setters
    void                        SetTotalADMult(Float_t mult)        { fADTotalMult   = mult; }
    void                        SetTotalADTime(Float_t tme)         { fADTotalTime   = tme;  }
    void                        SetTotalADCharge(Float_t chrg)      { fADTotalCharge = chrg; }
    // FMD total setters
    void                        SetTotalFMDMult(Float_t mult)       { fFMDTotalMult = mult; }
    // V0 total setters
    void                        SetTotalV0Mult(Float_t mult)        { fV0TotalMult   = mult; }
    void                        SetTotalV0Time(Float_t tme)         { fV0TotalTime   = tme;  }
    void                        SetTotalV0Charge(Float_t chrg)      { fV0TotalCharge = chrg; }
    void                        SetTotalV0SigWidth(Float_t sigW)    { fV0TotalSigWidth = sigW; }
    // Calo total setters
    void                        SetTotalEMCAmplitude(Float_t ampl)  { fEMCTotalAmplitude  = ampl; }
    void                        SetTotalEMCTime(Float_t tme)        { fEMCTotalTime       = tme;  }
    void                        SetTotalPHOSAmplitude(Float_t ampl) { fPHOSTotalAmplitude = ampl; }
    void                        SetTotalPHOSTime(Float_t tme)       { fPHOSTotalTime      = tme;  }
        
    /// Global variable setter
    void                        SetEventVariables(AliESDEvent* ESDobj, TArrayI* TTindices);

    /// Modifiers
    void                        Reset();


    // Accessors
    UInt_t                      GetEventNumber()        const { return fEventNumber; }
    UInt_t                      GetnTracksTotal()       const { return fnTracks;     }
    UInt_t                      GetnCaloClusterTotal()  const { return fnCaloTracks; }
    // AD total getters
    Float_t                     GetTotalADMult()        const { return fADTotalMult; }
    Float_t                     GetTotalADTime()        const { return fADTotalTime; }
    Float_t                     GetTotalADCharge()      const { return fADTotalCharge; }
    // FMD total getters
    Float_t                     GetTotalFMDMult()       const { return fFMDTotalMult; }
    // V0 total getters
    Float_t                     GetTotalV0Mult()        const { return fV0TotalMult;     }
    Float_t                     GetTotalV0Time()        const { return fV0TotalTime;     }
    Float_t                     GetTotalV0Charge()      const { return fV0TotalCharge;   }
    Float_t                     GetTotalV0SigWidth()    const { return fV0TotalSigWidth; }
    // Calo total getters
    Float_t                     GetTotalEMCAmplitude()  const { return fEMCTotalAmplitude;  }
    Float_t                     GetTotalEMCTime()       const { return fEMCTotalTime;       }
    Float_t                     GetTotalPHOSAmplitude() const { return fPHOSTotalAmplitude; }
    Float_t                     GetTotalPHOSTime()      const { return fPHOSTotalTime;      }
        
    /// Track accessors
    CEPRawTrackBuffer*          GetTrack(UInt_t ind);
    CEPRawCaloClusterTrack*     GetCaloClusterTrack(UInt_t ind);
    /// Detector object accessors
    CEPRawADBuffer*             GetRawADBuffer()        const { return fADCellBuffer; }
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
