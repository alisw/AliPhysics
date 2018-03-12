#ifndef CEPRawTrackBuffer_H
#define CEPRawTrackBuffer_H

#include "TObject.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"

class CEPRawTrackBuffer : public TObject 
{

  private:
    /// Get track lenght of TOF clusters
    Double_t        fTrackLength;
    /// Get Chi2 measures
    Double_t        fGlobalChi2;
        
    /// PID information 
    Double_t        fPIDITSsignal;
    Double_t        fPIDHMPIDsignal;
    Double_t        fPIDTRDsignal;
    Double_t        fPIDTOFsignalRaw;
    Double_t        fPIDITSsignalTuned;
    Double_t        fPIDTOFsignalTuned;
    Double_t        fPIDTPCsignalTuned;
    
    Double_t        fITSChi2;  
    Double_t        fTPCChi2;  
    /// Impact parameters
    Double_t        fXY;                // in the xy plane
    Double_t        fZ;                 // in the z  plane

    Double_t        fDCAx;
    Double_t        fDCAy;
    Double_t        fDCAz;

    Double_t        fDx;                // local x of tracks impact on the TOF pad [cm]
    Double_t        fDz;                // local z of tracks impact on the TOF pad [cm]

    /// eta, phi, pt and p after EMCal propagation
    Double_t        fTrkPhiOnEMCal;
    Double_t        fTrkEtaOnEMCal;
    Double_t        fTrkPtOnEMCal;
    Double_t        fTrkPOnEMCal;    

  public:
                    CEPRawTrackBuffer();
    virtual         ~CEPRawTrackBuffer() {};
    /// Modifiers
    void            Reset();

    /// Setter functions
    void            SetTrackLenght (Double_t trackLen)         { fTrackLength      = trackLen;  }
    void            SetGlobalChi2  (Double_t globChi2)         { fGlobalChi2       = globChi2;  }

    void            SetPIDITSsignal(Double_t pidITSsig)        { fPIDITSsignal     = pidITSsig; }
    void            SetPIDHMPsignal(Double_t pidHMPsig)        { fPIDHMPIDsignal   = pidHMPsig; }
    void            SetPIDTRDsignal(Double_t pidTRDsig)        { fPIDTRDsignal     = pidTRDsig; }
    void            SetPIDTOFsignalRaw(Double_t pidTOFrawsig)  { fPIDTOFsignalRaw  = pidTOFrawsig;}

    void            SetPIDITSsigTunedOnData(Double_t itsTuned) { itsTuned = fPIDITSsignalTuned; }
    void            SetPIDTOFsigTunedOnData(Double_t tofTuned) { tofTuned = fPIDTOFsignalTuned; }
    void            SetPIDTPCsigTunedOnData(Double_t tpcTuned) { tpcTuned = fPIDTPCsignalTuned; }

    void            SetITSChi2(Double_t itschi2) { fITSChi2  = itschi2; }
    void            SetTPCChi2(Double_t tpcchi2) { fTPCChi2  = tpcchi2; }

    void            SetImpactXY(Double_t impXY) { fXY = impXY; }
    void            SetImpactZ (Double_t impZ)  { fZ  = impZ;  }

    void            SetDCAx (Double_t dca_x) { fDCAx = dca_x; }
    void            SetDCAy (Double_t dca_y) { fDCAy = dca_y; }
    void            SetDCAz (Double_t dca_z) { fDCAz = dca_z; }

    void            SetLocalTOFImpX (Double_t tofImpX)  { fDx = tofImpX; }
    void            SetLocalTOFImpZ (Double_t tofImpZ)  { fDz = tofImpZ; }

    void            SetTrkPhiOnEMC (Double_t phiOnEMC)  { fTrkPhiOnEMCal = phiOnEMC; }
    void            SetTrkEtaOnEMC (Double_t etaOnEMC)  { fTrkEtaOnEMCal = etaOnEMC; }
    void            SetTrkPtOnEMC  (Double_t ptOnEMC)   { fTrkPtOnEMCal  = ptOnEMC;  }
    void            SetTrkPOnEMC   (Double_t POnEMC)    { fTrkPOnEMCal   = POnEMC;   }

    /// Setting the setters
    void            SetTrackVariables(AliESDtrack* track, AliESDVertex* vertex);
    
    /// Accessors
    Double_t        GetTrackLength()  const { return fTrackLength; }
    Double_t        GetGlobalChi2()   const { return fGlobalChi2;  }

    Double_t        GetPIDITSsig()    const { return fPIDITSsignal;    }
    Double_t        GetPIDHMPsig()    const { return fPIDHMPIDsignal;  }
    Double_t        GetPIDTRDsig()    const { return fPIDTRDsignal;    }
    Double_t        GetPIDTOFsigRaw() const { return fPIDTOFsignalRaw; }

    Double_t        GetPIDITSsigTunedOnData()  const { return fPIDITSsignalTuned; }
    Double_t        GetPIDTOFsigTunedOnData()  const { return fPIDTOFsignalTuned; }
    Double_t        GetPIDTPCsigTunedOnData()  const { return fPIDTPCsignalTuned; }

    Double_t        GetITSChi2()  const { return fITSChi2; }
    Double_t        GetTPCChi2()  const { return fTPCChi2; }
    
    Double_t        GetImpactXY()  const { return fXY; }
    Double_t        GetImpactZ()   const { return fZ;  }

    Double_t        GetDCAx()   const { return fDCAx; }
    Double_t        GetDCAy()   const { return fDCAy; }
    Double_t        GetDCAz()   const { return fDCAz; }

    Double_t        GetTOFImpactDx()  const { return fDx; }
    Double_t        GetTOFImpactDz()  const { return fDz; }

    Double_t        GetTrkPhiOnEMC()  const { return fTrkPhiOnEMCal; }
    Double_t        GetTrkEtaOnEMC()  const { return fTrkEtaOnEMCal; }
    Double_t        GetTrkPtOnEMC()   const { return fTrkPtOnEMCal;  }
    Double_t        GetTrkPOnEMC()    const { return fTrkPOnEMCal;   }
   

    ClassDef(CEPRawTrackBuffer,1)     // CEP raw track buffer

};

#endif
