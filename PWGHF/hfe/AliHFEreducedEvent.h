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
//
// Debug event to look at the distribution of the variable we are cutting on
//
//
#ifndef ALIHFEREDUCEDEVENT_H
#define ALIHFEREDUCEDEVENT_H

#include <TObject.h>

class TObjArray;
class AliHFEreducedTrack;
class AliHFEreducedMCParticle;

class AliHFEreducedEvent : public TObject{
 public:
  AliHFEreducedEvent();
  AliHFEreducedEvent(const AliHFEreducedEvent &ref);
  AliHFEreducedEvent &operator=(const AliHFEreducedEvent &ref);
  ~AliHFEreducedEvent();
  
  void AddTrack(const AliHFEreducedTrack *track);
  const AliHFEreducedTrack *GetTrack(int itrk) const;
  Int_t GetNumberOfTracks() const { return fNtracks; }
  void AddMCParticle(const AliHFEreducedMCParticle *mctrack);
  const AliHFEreducedMCParticle *GetMCParticle(int itrk) const;
  Int_t GetNumberOfMCParticles() const { return fNmcparticles; }
  
  Float_t GetVX() const { return fVX[0]; }
  Float_t GetVY() const { return fVY[0]; }
  Float_t GetVZ() const { return fVZ[0]; }
  Float_t GetVXSPD() const { return fVX[1]; }
  Float_t GetVYSPD() const { return fVY[1]; }
  Float_t GetVZSPD() const { return fVZ[1]; }
  Int_t GetNContribVertex() const { return fNContrib[0]; }
  Int_t GetNContribVertexSPD() const { return fNContrib[1]; }
  Bool_t HasPrimaryVertex() const { return fNContrib[0] > 0; }
  Bool_t HasPrimaryVertexSPD() const { return fNContrib[1] > 0; }
  Float_t GetVertexZResolution() const { return fVertexResolution[0]; };
  Float_t GetVertexZResolutionSPD() const { return fVertexResolution[1]; };
  Int_t GetRunNumber() const { return fRunNumber; }
  Double_t GetCentrality() const { return fCentrality[0]; }
  Double_t GetCentralityV0M() const { return fCentrality[0]; }
  Double_t GetCentralityV0A() const { return fCentrality[1]; }
  Double_t GetCentralityV0C() const { return fCentrality[2]; }
  Double_t GetCentralityTracklets() const { return fCentrality[3]; }
  Double_t GetCentralityTracks() const { return fCentrality[4]; }
  Double_t GetCentralityZNA() const { return fCentrality[5]; }
  Double_t GetCentralityZNC() const { return fCentrality[6]; }
  Double_t GetCentralityCL0() const { return fCentrality[7]; }
  Double_t GetCentralityCL1() const { return fCentrality[8]; }
  Double_t GetCentralityCND() const { return fCentrality[9]; }
  Float_t GetV0AMultiplicity() const { return fV0Multiplicity[0]; }
  Float_t GetV0CMultiplicity() const { return fV0Multiplicity[1]; }
  Float_t GetV0MMultiplicity() const { return fV0Multiplicity[0] + fV0Multiplicity[1]; }
  Float_t GetZNAEnergy() const { return fZDCEnergy[0]; }
  Float_t GetZNCEnergy() const { return fZDCEnergy[1]; }
  Float_t GetZPAEnergy() const { return fZDCEnergy[2]; }
  Float_t GetZPCEnergy() const { return fZDCEnergy[3]; }
  Float_t GetZDCNEnergySum() const { return fZDCEnergy[0] + fZDCEnergy[1]; }
  Float_t GetZDCNEnergyDifference() const { return fZDCEnergy[0] - fZDCEnergy[1]; }
  Float_t GetZDCNEnergyAsymmetry() const { return fZDCEnergy[0] + fZDCEnergy[1] != 0 ? (fZDCEnergy[0] - fZDCEnergy[1])/(fZDCEnergy[0] + fZDCEnergy[1]) : 1.; }
  Float_t GetZDCPEnergySum() const { return fZDCEnergy[2] + fZDCEnergy[3]; }
  Float_t GetZDCPEnergyDifference() const { return fZDCEnergy[2] - fZDCEnergy[3]; }
  Float_t GetZDCPEnergyAsymmetry() const { return fZDCEnergy[2] + fZDCEnergy[3] != 0 ? (fZDCEnergy[2] - fZDCEnergy[3])/(fZDCEnergy[2] + fZDCEnergy[3]) : 1.; }
  Int_t   GetSPDMultiplicity() const { return fSPDMultiplicity; }
  
  void SetVX(Float_t vx) { fVX[0] = vx; }
  void SetVY(Float_t vy) { fVY[0] = vy; }
  void SetVZ(Float_t vz) { fVZ[0] = vz; }
  void SetVXSPD(Float_t vx) { fVX[1] = vx; }
  void SetVYSPD(Float_t vy) { fVY[1] = vy; }
  void SetVZSPD(Float_t vz) { fVZ[1] = vz; }
  void SetRunNumber(Int_t runnumber) { fRunNumber = runnumber; }
  inline void SetCentrality(
        Float_t centV0M, 
        Float_t centV0A, 
        Float_t centV0C, 
        Float_t centTLS, 
        Float_t centTrks, 
        Float_t centZNA, 
        Float_t centZNC,
        Float_t centCL0,
        Float_t centCL1,
        Float_t centCND
  );
  void SetNContribVertex(Int_t ncontrib) { fNContrib[0] = ncontrib; }
  void SetNContribVertexSPD(Int_t ncontrib) { fNContrib[1] = ncontrib; }
  void SetVertexResolution(Float_t res) { fVertexResolution[0] = res; }
  void SetVertexResolutionSPD(Float_t res) { fVertexResolution[1] = res; }
  
  Bool_t IsMBTrigger() const { return TESTBIT(fTrigger, kMB); }
  Bool_t IsSemiCentralTrigger() const { return TESTBIT(fTrigger, kSemiCentral); }
  Bool_t IsCentralTrigger() const { return TESTBIT(fTrigger, kCentral); }
  Bool_t IsEMCalTrigger() const { return TESTBIT(fTrigger, kEMCAL); }
  Bool_t IsTRDSETrigger() const { return TESTBIT(fTrigger, kTRDSE); }
  Bool_t IsTRDDQTrigger() const { return TESTBIT(fTrigger, kTRDDQ); }
  Bool_t IsINTTrigger() const { return TESTBIT(fTrigger, kINTTRG); }

  void SetMBTrigger() { SETBIT(fTrigger, kMB); }
  void SetSemiCentralTrigger() { SETBIT(fTrigger, kSemiCentral); }
  void SetCentralTrigger() { SETBIT(fTrigger, kCentral); }
  void SetEMCALTrigger() { SETBIT(fTrigger, kEMCAL); }
  void SetTRDSETrigger() { SETBIT(fTrigger, kTRDSE); }
  void SetTRDDQTrigger() { SETBIT(fTrigger, kTRDDQ); }
  void SetINTTrigger() { SETBIT(fTrigger, kINTTRG); }

  void SetV0Multiplicity(Float_t v0A, Float_t v0C) {
    fV0Multiplicity[0] = v0A;
    fV0Multiplicity[1] = v0C;
  }
  inline void SetZDCEnergy(Float_t zna, Float_t znc, Float_t zpa, Float_t zpc);
  void SetSPDMultiplicity(Int_t mult) { fSPDMultiplicity = mult; }
  
 private:
  typedef enum{
    kMB = 0,
    kSemiCentral = 1,
    kCentral = 2,
    kEMCAL = 3,
    kTRDSE = 4,
    kTRDDQ = 5,
    kINTTRG = 6 
  } Trigger_t;
  enum{
    kCentBuff = 15
  };
  TObjArray *fTracks;           // Array with reconstructed tracks
  TObjArray *fMCparticles;      // Array with MC particles
  Int_t fNtracks;               // Number of tracks
  Int_t fNmcparticles;          // Number of MC Particles
  Int_t fRunNumber;             // Run Number
  Float_t  fCentrality[kCentBuff];     // Centrality (V0M, V0A, V0C, TLS, TRK, ZNA, ZNC, CL0, CL1, CND)
  Int_t fTrigger;               // Trigger bits
  Float_t fVX[2];               // Vertex X
  Float_t fVY[2];               // Vertex Y
  Float_t fVZ[2];               // Vertex Z
  Int_t    fNContrib[2];        // Number of vertex contributors
  Float_t  fVertexResolution[2];// z-Vertex resolution 
  Float_t  fV0Multiplicity[2];  // V0 multiplicity
  Float_t  fZDCEnergy[4];       // ZDC Energy (n,p)
  Int_t    fSPDMultiplicity;    // SPD tracklet multiplicity
  
  ClassDef(AliHFEreducedEvent, 4)
};

//____________________________________________________________
void AliHFEreducedEvent::SetCentrality(
        Float_t centV0M, 
        Float_t centV0A, 
        Float_t centV0C, 
        Float_t centTLS, 
        Float_t centTrks, 
        Float_t centZNA, 
        Float_t centZNC,
        Float_t centCL0,
        Float_t centCL1,
        Float_t centCND
  ) 
{ 
    fCentrality[0] = centV0M; 
    fCentrality[1] = centV0A; 
    fCentrality[2] = centV0C; 
    fCentrality[3] = centTLS; 
    fCentrality[4] = centTrks; 
    fCentrality[5] = centZNA; 
    fCentrality[6] = centZNC; 
    fCentrality[7] = centCL0; 
    fCentrality[8] = centCL1; 
    fCentrality[9] = centCND; 
}

//_________________________________________________________
void AliHFEreducedEvent::SetZDCEnergy(Float_t zna, Float_t znc, Float_t zpa, Float_t zpc){
  fZDCEnergy[0] = zna;
  fZDCEnergy[1] = znc;
  fZDCEnergy[2] = zpa;
  fZDCEnergy[3] = zpc;
}
#endif
