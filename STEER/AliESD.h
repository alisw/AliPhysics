#ifndef ALIESD_H
#define ALIESD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//   Alice ESD object                                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TObjArray.h"
#include "TArrayI.h"
#include "TArrayF.h"
#include "TArrayD.h"

class AliESDTrack : public TObject
{
public:
  AliESDTrack();
  virtual ~AliESDTrack() {}
  
protected:
  Int_t     fTrackID;        // Track number

  // Track parameters at Vertex
  TArrayD   fPVertex;        // (5) Track parameters
  TArrayD   fPEVertex;       // (15) Track parameter errors
  
  // Track parameters at first measured point
  TArrayD   fPFMeasPoint;        // (6) Track parameters
  TArrayD   fPFMeasPointErr;     // (15) Track parameter error
  
  // Track parameters at last measured point
  TArrayD   fPLMeasPoint;        // (6) Track parameters
  TArrayD   fPLMeasPointErr;     // (15) Track parameter error

  Float_t   fTrackLength;        // Track length
  Float_t   fTrackLengthErr;     // Track length error
  Int_t     fStopVertex;         // Index of stop vertex
  
  Int_t     fNPointsITS;         // Number of points in ITS
  Int_t     fNPointsTPC;         // Number of points in TPC
  Int_t     fNPointsTRD;         // Number of points in TRD
  Float_t   fMeanResITS;         // Mean residual in ITS
  Float_t   fMeanResTPC;         // Mean residual in TPC
  Float_t   fMeanResTRD;         // Mean residual in TRD
  Float_t   fGlobalChi2;         // Global chi square
  Int_t     fParticleType;       // PDG code

  Float_t   fPIDprobPi;          // PID probability for pi
  Float_t   fPIDprobK;           // PID probability for K
  Float_t   fPIDprobP;           // PID probability for p
  Float_t   fPIDprobE;           // PID probability for e

private:
  AliESDTrack(const AliESDTrack &);
  AliESDTrack & operator=(const AliESDTrack &) {return (*this);}
  
  ClassDef(AliESDTrack,1)  //ESDTrack 
};


class AliESDVertex : public TObject
{
public:
  AliESDVertex();
  virtual ~AliESDVertex() {}
  
protected:
  Int_t        fNPrimary;               // Number of primary tracks
  TArrayF      fCoordinates;            // (3) Vertex coordinates
  TArrayF      fErrorMatrix;            // (6) Error Matrix
  TObjArray    fPrimaryTracks;          // List of primary tracks
  Float_t      fEffectiveMass;          // Effective Mass
  Float_t      fEffectiveMassError;     // Effective Mass Error
private:
  AliESDVertex(const AliESDVertex &);
  AliESDVertex & operator=(const AliESDVertex &) {return (*this);}
  
  ClassDef(AliESDVertex,1)  //ESDVertex 
};

class AliESD : public TObject
{
public:
  AliESD();
  virtual ~AliESD() {}

  Int_t EventNumber() const {return fEventNumber;}
  Int_t RunNumber() const {return fRunNumber;}
  Long_t Trigger() const {return fTrigger;}
  
  Int_t BitDDL() const {return fBitDDL;}
  Int_t NSecVertex() const {return fNSecVertex;}
  Float_t NParticipants() const {return fNParticipants;}
  
  
protected:

  // Event Identification
  Int_t        fEventNumber;            // Event Number
  Int_t        fRunNumber;              // Run Number
  Long_t       fTrigger;                // Trigger Type (cfg Transverse Energy&Max trans ch mom)
  Int_t        fRecoVersion;            // Version of reconstruction 

  // Summary Information
  Int_t        fBitDDL;                 // Bitmap of active DDL
  Int_t        fNSecVertex;             // Number of Secondary Vertexes
  Float_t      fNParticipants;          // Estimated Number of participants
  Float_t      fNPartError;             // N of participant error
  Int_t        fNElectron;              // N of electrons
  Int_t        fNMuons;                 // N of muons
  Int_t        fNPions;                 // N of pions
  Int_t        fNKaons;                 // N of kaons
  Int_t        fNProtons;               // N of protons
  Int_t        fNPHOSPhotons;           // N of photons in PHOS
  Int_t        fNPHOSNeutrons;          // N of neutrons in PHOS
  Int_t        fNPHOSCCluster;          // N of charged clusters in PHOS
  Int_t        fNEMCALCluster;          // N of clusters in EMCAL
  Int_t        fNPMDCluster;            // N of clusters in PMD
  Float_t      fTMaxClusterEnergy;      // Transverse energy of biggest cluster
  Float_t      fTMaxPCharged;           // Biggest transverse momentum of charged particles
  TArrayI      fNCharged;               // Charged Multiplicity
  Float_t      fTotTranEnergy;          // Total transverse energy

  // Primary Vertex Object
  AliESDVertex fESDVertex;              // Primary Vertex Object
  TObjArray    fSecVertex;              // List secondary vertexes
  TObjArray    fNonAssTrack;            // List of non assigned tracks
  TObjArray    fPhoton;                 // List of photons
  TObjArray    fNeutron;                // List of neutrons
  TObjArray    fEMCALCluster;           // List of EMCAL clusters
  TObjArray    fPMDCluster;             // List of PMD clusters

private:
  AliESD(const AliESD &);
  AliESD & operator=(const AliESD &) {return (*this);}
  
  ClassDef(AliESD,1)  //ESD 
};

#endif 

