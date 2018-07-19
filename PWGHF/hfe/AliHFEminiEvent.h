/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
// Originator:  M.Fasel <M.Fasel@gsi.de>
// Base class for AliHFEminiEventCreator.cxx/h
// Contributors: Nirbhay K. Behera, Jiyeon Kwon, Jonghan Park

#ifndef ALIHFEMINIEVENT_H
#define ALIHFEMINIEVENT_H

#include <TObject.h>

class TObjArray;
class AliHFEreducedMCParticle;
class AliHFEminiTrack;


class AliHFEminiEvent : public TObject{
 public:
  AliHFEminiEvent();
  AliHFEminiEvent(const AliHFEminiEvent &ref);
  AliHFEminiEvent &operator=(const AliHFEminiEvent &ref);
  ~AliHFEminiEvent();
  
  //getters-------
  void AddTrack(const AliHFEminiTrack *track);
  const AliHFEminiTrack *GetTrack(int itrk) const;
  Int_t GetNumberOfTracks() const { return fNtracks; }
  void AddMCParticle(const AliHFEreducedMCParticle *mctrack);
  const AliHFEreducedMCParticle *GetMCParticle(int itrk) const;
  Int_t GetNumberOfMCParticles() const { return fNmcparticles; }
  
  Float_t GetVZ() const { return fVZ; }
  Float_t GetV0AMultiplicity() const { return fV0Multiplicity[0]; }
  Float_t GetV0CMultiplicity() const { return fV0Multiplicity[1]; }
  Float_t GetV0MMultiplicity() const { return fV0Multiplicity[0] + fV0Multiplicity[1]; }
  Int_t   GetSPDMultiplicity() const { return fSPDMultiplicity; }
  Int_t   GetCentralityBin() const   { return fCentralityBin; }
  Int_t   GetPrimaryNchMC()    const { return fNprimaryNchMC; }   //Get the Phy. Primary track MC
  
  //setters---
  void SetVZ(Float_t vz) { fVZ = vz; }
  void SetV0Multiplicity(Float_t v0A, Float_t v0C) {
    fV0Multiplicity[0] = v0A;
    fV0Multiplicity[1] = v0C;
  }
  void SetSPDMultiplicity(Int_t mult) { fSPDMultiplicity = mult; }
  void SetCentralityBin(Int_t centralitybin){ fCentralityBin = centralitybin;}
  void SetPrimaryNchMC( Int_t primaryTrack) { fNprimaryNchMC = primaryTrack; }
  
 private:
  TObjArray      *fTracks;           // Array with reconstructed tracks
  TObjArray      *fMCparticles;      // Array with MC particles
  Int_t          fNtracks;               // Number of tracks
  Int_t          fNmcparticles;          // Number of MC Particles
  
  Float_t        fVZ;               // Vertex Z
  Float_t        fV0Multiplicity[2];  // V0 multiplicity
  Int_t          fSPDMultiplicity;    // SPD tracklet multiplicity
  Int_t          fCentralityBin;
  Int_t          fNprimaryNchMC;  //Phy. Primary from MC
  
  ClassDef(AliHFEminiEvent, 6)
    };

#endif

