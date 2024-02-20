#ifndef AliResoNanoMCParticle_H
#define AliResoNanoMCParticle_H

#include <Rtypes.h>
#include <TNamed.h>

// Basic class for mc particle information

class AliResoNanoMCParticle : public TNamed
{
public:
    AliResoNanoMCParticle();
    virtual ~AliResoNanoMCParticle();

    void SetID(Int_t id) { fID = id; }
    void SetEventID(Int_t id) { fEventID = id; }
    void SetMotherID(Int_t id) { fMotherID = id; }
    void SetDaughter(Int_t i, Int_t id) { fDaughter[i] = id; }
    void SetPDGCode(Int_t pdg) { fPDGCode = pdg; }
    void SetPt(Double32_t pt) { fPt = pt; }
    void SetEta(Double32_t eta) { fEta = eta; }
    void SetPhi(Double32_t phi) { fPhi = phi; }
    void SetRap(Double32_t rap) { fRap = rap; }
    void SetCharge(Char_t charge) { fCharge = charge; }
    void SetStatus(Int_t status) { fStatus = status; }

    Int_t GetID() const { return fID; }
    Int_t GetEventID() const { return fEventID; }
    Int_t GetMotherID() const { return fMotherID; }
    Int_t GetDaughter(Int_t i) const { return fDaughter[i]; }
    Int_t GetPDGCode() const { return fPDGCode; }
    Double32_t GetPt() const { return fPt; }
    Double32_t GetEta() const { return fEta; }
    Double32_t GetPhi() const { return fPhi; }
    Double32_t GetRap() const { return fRap; }
    Char_t GetCharge() const { return fCharge; }
    Int_t GetStatus() const { return fStatus; }

private:
    Int_t fID;          // id
    Int_t fEventID;     // event id
    Int_t fMotherID;    // mother id
    Int_t fDaughter[2]; // daughter id
    Int_t fPDGCode;     // pdg code
    Double32_t fPt;     // pt
    Double32_t fEta;    // eta
    Double32_t fPhi;    // phi
    Double32_t fRap;    // rapidity
    Char_t fCharge;     // charge
    Int_t fStatus;      // status

    ClassDef(AliResoNanoMCParticle, 1);
};

#endif
