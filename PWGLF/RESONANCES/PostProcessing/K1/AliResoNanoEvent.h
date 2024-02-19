#ifndef AliResoNanoEvent_H
#define AliResoNanoEvent_H

#include <Rtypes.h>
#include <TNamed.h>

// Basic class for event information
// Ref: AliJBaseEventHeader.h

class AliResoNanoEvent : public TNamed
{
public:
    AliResoNanoEvent();
    virtual ~AliResoNanoEvent();

    void SetEventID(Int_t id) { fEventID = id; }
    void SetCentrality(Double32_t c) { fCentrality = c; }
    void SetVertex(Double32_t x, Double32_t y, Double32_t z, Double32_t zerr = 0.0)
    {
        fVtxX = x;
        fVtxY = y;
        fVtxZ = z;
        fVtxZErr = zerr;
    }
    void SetVertexMC(Double32_t x, Double32_t y, Double32_t z)
    {
        fVtxMCX = x;
        fVtxMCY = y;
        fVtxMCZ = z;
    }

    Int_t GetEventID() const { return fEventID; }
    Double32_t GetCentrality() const { return fCentrality; }
    Double32_t GetVertexX() const { return fVtxX; }
    Double32_t GetVertexY() const { return fVtxY; }
    Double32_t GetVertexZ() const { return fVtxZ; }
    Double32_t GetVertexZErr() const { return fVtxZErr; }
    Double32_t GetVertexMCX() const { return fVtxMCX; }
    Double32_t GetVertexMCY() const { return fVtxMCY; }
    Double32_t GetVertexMCZ() const { return fVtxMCZ; }

private:
    Int_t fEventID;         // event id
    Double32_t fCentrality; // centrality
    Double32_t fVtxX;       // vertex X
    Double32_t fVtxY;       // vertex Y
    Double32_t fVtxZ;       // vertex Z
    Double32_t fVtxZErr;    // vertex error

    Double32_t fVtxMCX; // vertex X MC
    Double32_t fVtxMCY; // vertex Y MC
    Double32_t fVtxMCZ; // vertex Z MC

    ClassDef(AliResoNanoEvent, 1);
};

#endif