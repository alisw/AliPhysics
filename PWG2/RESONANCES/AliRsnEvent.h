/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

/* $Id$ */

//
// *** Class AliRsnEvent ***
//
// A container for a collection of AliRsnDaughter objects from an event.
// Contains also the primary vertex, useful for some cuts.
// In order to retrieve easily the tracks which have been identified
// as a specific type and charge, there is an array of indexes which
// allows to avoid to loop on all tracks and have only the neede ones.
//
// authors: A. Pulvirenti (email: alberto.pulvirenti@ct.infn.it)
//          M. Vala (email: martin.vala@cern.ch)
//

#ifndef ALIRSNEVENT_H
#define ALIRSNEVENT_H

#include <TNamed.h>
#include <TArrayI.h>
#include <TClonesArray.h>

#include "AliRsnPID.h"
#include "AliRsnPIDIndex.h"
#include "AliRsnDaughter.h"

class AliRsnEvent : public TNamed
{
  public:

    AliRsnEvent();
    AliRsnEvent(const AliRsnEvent& copy);
    AliRsnEvent& operator= (const AliRsnEvent& copy);
    virtual ~AliRsnEvent();

    // Array management
    void            Init();
    void            Clear(Option_t *option = "");
    AliRsnDaughter* AddTrack(AliRsnDaughter track);
    AliRsnDaughter* GetTrack(Int_t index);
    AliRsnDaughter* GetLeadingParticle(Double_t ptMin = 0.0, AliRsnPID::EType type = AliRsnPID::kUnknown, Bool_t realistic = kTRUE);
    Int_t           GetLastFastTrack(Double_t ptMin, AliRsnPID::EType type = AliRsnPID::kUnknown, Bool_t realistic = kTRUE);
    TClonesArray*   GetTracks() {return fTracks;}
    TArrayI*        GetCharged(Char_t sign);
    TArrayI*        GetTracksArray(AliRsnDaughter::EPIDMethod method, Char_t sign, AliRsnPID::EType type);
    void            FillPIDArrays();
    void            SortTracks() {fTracks->Sort();}
    void            Print(Option_t *option = "") const;

    // Primary vertex
    Double_t GetPrimaryVertexX() const {return fPVx;}
    Double_t GetPrimaryVertexY() const {return fPVy;}
    Double_t GetPrimaryVertexZ() const {return fPVz;}
    void     GetPrimaryVertex(Double_t &x, Double_t &y, Double_t &z) const {x=fPVx;y=fPVy;z=fPVz;}
    void     SetPrimaryVertexX(Double_t value) {fPVx = value;}
    void     SetPrimaryVertexY(Double_t value) {fPVy = value;}
    void     SetPrimaryVertexZ(Double_t value) {fPVz = value;}
    void     SetPrimaryVertex(Double_t x, Double_t y, Double_t z) {fPVx=x;fPVy=y;fPVz=z;}
    void     CorrectByPrimaryVertex();

    // Multiplicity
    Int_t GetMultiplicity() const;
    Int_t GetNCharged(Char_t sign);

  private:

    Int_t ChargeIndex(Char_t sign) const;
    Int_t Fill(TObjArray *array);

    Double_t        fPVx;                 // position of
    Double_t        fPVy;                 // primary
    Double_t        fPVz;                 // vertex

    TClonesArray   *fTracks;              // collection of particles

    AliRsnPIDIndex *fNoPID;               // array index only for charged tracks
    AliRsnPIDIndex *fPerfectPID;          // array index for perfect PID
    AliRsnPIDIndex *fRealisticPID;        // array index for realistic PID (largest prob)

    ClassDef(AliRsnEvent, 2);
};

#endif
