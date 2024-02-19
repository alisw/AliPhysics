#ifndef AliResoNanoTrack_H
#define AliResoNanoTrack_H

#include <Rtypes.h>
#include <TLorentzVector.h>

// Basic class for track information
// Ref: AliJBaseTrack.h

class AliResoNanoTrack : public TLorentzVector
{
public:
    AliResoNanoTrack();
    AliResoNanoTrack(float px, float py, float pz, float e, Int_t id, Short_t ptype, Char_t charge);
    AliResoNanoTrack(const AliResoNanoTrack &a);
    AliResoNanoTrack(const TLorentzVector &a);
    virtual ~AliResoNanoTrack();

    double EtaAbs() { return TMath::Abs(Eta()); }
    TLorentzVector GetLorentzVector() { return TLorentzVector(Px(), Py(), Pz(), E()); }

    Int_t GetID() const { return fID; }
    Int_t GetLabel() const { return fLabel; }
    Short_t GetParticleType() const { return fParticleType; }
    ULong64_t GetStatus() const { return fStatus; }
    Char_t GetCharge() const { return fCharge; }
    UInt_t GetFlags() const { return fFlags; }
    Double_t GetDCAxy() const { return fDCAxy; }
    Double_t GetDCAz() const { return fDCAz; }
    Double_t GetTOFNSigmaPi() const { return fTOFNSigmaPi; }
    Double_t GetTPCNSigmaPi() const { return fTPCNSigmaPi; }
    Double_t GetTOFNSigmaKa() const { return fTOFNSigmaKa; }
    Double_t GetTPCNSigmaKa() const { return fTPCNSigmaKa; }
    Double_t GetTOFNSigmaPr() const { return fTOFNSigmaPr; }
    Double_t GetTPCNSigmaPr() const { return fTPCNSigmaPr; }
    Int_t GetMCPDGCode() const { return fMCPDGCode; }
    Int_t GetMCMotherPDGCode() const { return fMCMotherPDGCode; }
    Int_t GetMCMotherID() const { return fMCMotherID; }

    void SetID(const int id) { fID = id; }
    void SetLabel(const Int_t label) { fLabel = label; }
    void SetParticleType(const Short_t ptype) { fParticleType = ptype; }
    void SetStatus(const ULong64_t status) { fStatus = status; }
    void SetCharge(const Char_t charge) { fCharge = charge; }
    void SetFlags(const UInt_t bits) { fFlags = bits; } // MC, is primary flag

    void SetDCAxy(const Double_t dca) { fDCAxy = dca; }
    void SetDCAz(const Double_t dca) { fDCAz = dca; }
    void SetTOFNSigmaPi(const Double_t nsigma) { fTOFNSigmaPi = nsigma; }
    void SetTPCNSigmaPi(const Double_t nsigma) { fTPCNSigmaPi = nsigma; }
    void SetTOFNSigmaKa(const Double_t nsigma) { fTOFNSigmaKa = nsigma; }
    void SetTPCNSigmaKa(const Double_t nsigma) { fTPCNSigmaKa = nsigma; }
    void SetTOFNSigmaPr(const Double_t nsigma) { fTOFNSigmaPr = nsigma; }
    void SetTPCNSigmaPr(const Double_t nsigma) { fTPCNSigmaPr = nsigma; }
    void SetMCPDGCode(const Int_t pdg) { fMCPDGCode = pdg; }
    void SetMCMotherPDGCode(const Int_t pdg) { fMCMotherPDGCode = pdg; }
    void SetMCMotherID(const Int_t id) { fMCMotherID = id; }

private:
    Int_t fID;              // track id
    Short_t fParticleType;  // particle type
    Char_t fCharge;         // charge
    Int_t fLabel;           // label
    ULong64_t fStatus;      // status
    UInt_t fFlags;          // flags
    Double_t fDCAxy;        // DCA xy
    Double_t fDCAz;         // DCA z
    Double_t fTOFNSigmaPi;  // TOF n sigma pion
    Double_t fTPCNSigmaPi;  // TPC n sigma pion
    Double_t fTOFNSigmaKa;  // TOF n sigma kaon
    Double_t fTPCNSigmaKa;  // TPC n sigma kaon
    Double_t fTOFNSigmaPr;  // TOF n sigma proton
    Double_t fTPCNSigmaPr;  // TPC n sigma proton
    Int_t fMCPDGCode;       // MC PDG code
    Int_t fMCMotherPDGCode; // MC mother PDG code
    Int_t fMCMotherID;      // MC mother id

    ClassDef(AliResoNanoTrack, 1);
};

#endif
