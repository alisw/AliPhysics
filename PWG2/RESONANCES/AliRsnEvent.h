/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

/* $Id$ */

//-------------------------------------------------------------------------
//                      Class AliRsnEvent
//  Simple collection of reconstructed tracks, selected from an ESD event
// 
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//-------------------------------------------------------------------------

#ifndef ALIRSNEVENT_H
#define ALIRSNEVENT_H

class TTree;
class TParticle;
class TRefArray;
class TClonesArray;

class AliESDtrack;
class AliAODTrack;
class AliRsnDaughter;

#include "AliRsnPID.h"

class AliRsnEvent : public TObject
{
public:

    enum ESource {
        kUnknown = 0,
        kESD,
        kAOD,
        kMC
    };
	
    AliRsnEvent();
    AliRsnEvent(const AliRsnEvent& copy);
    AliRsnEvent& operator=(const AliRsnEvent& copy);
    virtual ~AliRsnEvent();

    // Array management
    void            Init();
    AliRsnDaughter* AddTrack(AliRsnDaughter track);
    void            Clear(Option_t *option = "");
    void            Print(Option_t *option = "") const;
    TClonesArray*   GetAllTracks() {return fTracks;}
    TRefArray*      GetTracks(Char_t sign, AliRsnPID::EType type);
    TRefArray*      GetCharged(Char_t sign) {if (sign=='+') return fPos; else if (sign=='-') return fNeg; else return 0x0;}
    void            FillPIDArrays();

    // Primary vertex
    Double_t GetPrimaryVertexX() const {return fPVx;}
	Double_t GetPrimaryVertexY() const {return fPVy;}
	Double_t GetPrimaryVertexZ() const {return fPVz;}
    void     GetPrimaryVertex(Double_t &x, Double_t &y, Double_t &z) const {x=fPVx;y=fPVy;z=fPVz;}
    void     GetPrimaryVertex(Double_t* &v) {GetPrimaryVertex(v[0], v[1], v[2]);}
    void     SetPrimaryVertexX(Double_t value) {fPVx = value;}
	void     SetPrimaryVertexY(Double_t value) {fPVy = value;}
	void     SetPrimaryVertexZ(Double_t value) {fPVz = value;}
    void     SetPrimaryVertex(Double_t x, Double_t y, Double_t z) {fPVx=x;fPVy=y;fPVz=z;}
    void     SetPrimaryVertex(Double_t* &v) {SetPrimaryVertex(v[0], v[1], v[2]);}

    // Multiplicity
    Int_t GetMultiplicity() const;
    Int_t GetNPos() const;
    Int_t GetNNeg() const;
    
    // Source type
    void    SetSource(ESource source) {fSource = source;}
    ESource GetSource() {return fSource;}
    
private:

    ESource        fSource;                          // type of source event

    Double_t       fPVx;                             // position of
    Double_t       fPVy;                             // primary
    Double_t       fPVz;                             // vertex

    TClonesArray  *fTracks;                          // collection of particles
    TRefArray     *fPos;                             // ref array to all positive particles
    TRefArray     *fNeg;                             // ref array to all negative particles
    TRefArray     *fPosID[AliRsnPID::kSpecies + 1];  // ref array to pos particles of each PID
    TRefArray     *fNegID[AliRsnPID::kSpecies + 1];  // ref array to pos particles of each PID

    ClassDef(AliRsnEvent,1);
};

#endif
