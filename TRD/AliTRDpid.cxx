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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//   The TRD particle identification base class                              //
//                                                                           //
//   Its main purposes are:                                                  //
//      - Provide I/O framework for all neccessary files                     //
//      - Assignment of a e/pi propability to a given track                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <math.h>

#include <TROOT.h>
#include <TH1.h>
#include <TObjArray.h>
#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>

#include "AliRun.h"
#include "AliTRDpid.h"
#include "AliTRDcluster.h"
#include "AliTRDtrack.h"
#include "AliTRDtracker.h"
#include "AliTRDgeometry.h"
#include "AliMC.h"

ClassImp(AliTRDpid)

//_____________________________________________________________________________
AliTRDpid::AliTRDpid():TNamed()
{
  //
  // AliTRDpid default constructor
  // 

  fTrackArray    = NULL;
  fClusterArray  = NULL;
  fGeometry      = NULL;
  fFileKine      = NULL;

  fPIDratioMin   = 0.0;
  fPIDpurePoints = kFALSE;
  fPIDindexMin   = 0;
  fPIDindexMax   = 0;

  fEvent         = 0;

  fThreePadOnly  = kFALSE;

}

//_____________________________________________________________________________
AliTRDpid::AliTRDpid(const char* name, const char* title):TNamed(name,title)
{
  //
  // AliTRDpid constructor
  //

  fTrackArray   = NULL;
  fClusterArray = NULL;
  fGeometry     = NULL;
  fFileKine     = NULL;

  fEvent        = 0;

  Init();

}

//_____________________________________________________________________________
AliTRDpid::AliTRDpid(const AliTRDpid &p):TNamed(p)
{
  //
  // AliTRDpid copy constructor
  //

  ((AliTRDpid &) p).Copy(*this);

}

//_____________________________________________________________________________
AliTRDpid::~AliTRDpid()
{
  //
  // AliTRDpid destructor
  //

  if (fClusterArray) {
    fClusterArray->Delete();
    delete fClusterArray;
    fClusterArray = NULL;
  }

  if (fTrackArray) {
    fTrackArray->Delete();
    delete fTrackArray;
    fTrackArray = NULL;
  }

  if (fFileKine) {
    delete fFileKine;
    fFileKine = NULL;
  }

}

//_____________________________________________________________________________
AliTRDpid &AliTRDpid::operator=(const AliTRDpid &p)
{
  //
  // Assignment operator
  //

  if (this != &p) ((AliTRDpid &) p).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDpid::Copy(TObject &p) const
{
  //
  // Copy function
  //

  ((AliTRDpid &) p).fTrackArray    = NULL;    
  ((AliTRDpid &) p).fClusterArray  = NULL;    
  ((AliTRDpid &) p).fGeometry      = NULL;    
  ((AliTRDpid &) p).fFileKine      = NULL;
  ((AliTRDpid &) p).fPIDratioMin   = fPIDratioMin;
  ((AliTRDpid &) p).fPIDpurePoints = fPIDpurePoints;
  ((AliTRDpid &) p).fPIDindexMin   = fPIDindexMin;
  ((AliTRDpid &) p).fPIDindexMax   = fPIDindexMax;
  ((AliTRDpid &) p).fThreePadOnly  = fThreePadOnly;
  ((AliTRDpid &) p).fEvent         = fEvent;

}

//_____________________________________________________________________________
Bool_t AliTRDpid::Init()
{
  //
  // Initializes the PID object 
  //

  fClusterArray  = new TObjArray();
  fTrackArray    = new TObjArray();

  fPIDratioMin   = 0.7;
  fPIDpurePoints = kTRUE;
  fPIDindexMin   = -1;
  fPIDindexMax   = -1;

  fThreePadOnly  = kFALSE;

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDpid::AssignLikelihood()
{
  //
  // Assigns the e / pi likelihood to all tracks.
  //

  return AssignLikelihood(fTrackArray);

}

//_____________________________________________________________________________
Bool_t AliTRDpid::AssignLikelihood(TObjArray *tarray)
{
  //
  // Assigns the e / pi likelihood to all tracks in the array <tarray>.
  //

  Bool_t status = kTRUE;

  AliTRDtrack *track;

  TIter nextTrack(tarray);  
  while ((track = (AliTRDtrack *) nextTrack())) {
    if (!AssignLikelihood(track)) {
      status = kFALSE;
      continue;
    }
  }

  return status;

}

//_____________________________________________________________________________
Bool_t AliTRDpid::FillSpectra()
{
  //
  // Fills the energy loss distributions with all tracks.
  //

  return FillSpectra(fTrackArray);

}
  
//_____________________________________________________________________________
Bool_t AliTRDpid::FillSpectra(TObjArray *tarray)
{
  //
  // Fills the energy loss distributions with all tracks in <tarray>.
  //
  
  Bool_t status = kTRUE;

  AliTRDtrack *track;

  TIter nextTrack(tarray);
  while ((track = (AliTRDtrack *) nextTrack())) {
    if (!FillSpectra(track)) {
      status = kFALSE;
      continue;
    }
  }  

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDpid::Open(const Char_t *name, Int_t event)
{
  //
  // Opens and reads the kine tree, the geometry, the cluster array
  // and the track array from the file <name>.
  //

  Bool_t status = kTRUE;

  status = ReadCluster(name);
  status = ReadTracks(name);
  status = ReadKine(name,event);

  return status;

}

//_____________________________________________________________________________
Bool_t AliTRDpid::Open(const Char_t *namekine
                     , const Char_t *namecluster
                     , const Char_t *nametracks
                     , Int_t event)
{
  //
  // Opens and reads the kine tree and the geometry from file <namekine>,
  // the cluster array from file <namecluster>,
  // and the track array from the file <nametracks>.
  //

  Bool_t status = kTRUE;

  status = ReadCluster(namecluster);
  status = ReadTracks(nametracks);
  status = ReadKine(namekine,event);

  return status;

}

//_____________________________________________________________________________
Bool_t AliTRDpid::ReadKine(const Char_t *name, Int_t event)
{
  //
  // Opens and reads the kine tree and the geometry from the file <name>.
  //

  TFile *fFileKine = (TFile *) gROOT->GetListOfFiles()->FindObject(name);
  if (!fFileKine) {
    printf("AliTRDpid::ReadKine -- ");
    printf("Open file %s\n",name);
    fFileKine = new TFile(name);
    if (!fFileKine) {
      printf("AliTRDpid::ReadKine -- ");
      printf("Cannot open file %s\n",name);
      return kFALSE;
    }
  }

  gAlice = (AliRun *) fFileKine->Get("gAlice");
  if (!gAlice) {
    printf("AliTRDpid::ReadKine -- ");
    printf("No AliRun object found\n");    
    return kFALSE;
  }
  gAlice->GetEvent(event);

  fGeometry = AliTRDgeometry::GetGeometry(gAlice->GetRunLoader());
  if (!fGeometry) {
    printf("AliTRDpid::ReadKine -- ");
    printf("No TRD geometry found\n");
    return kFALSE;
  }

  fEvent = event;

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDpid::ReadCluster(const Char_t *name)
{
  //
  // Opens and reads the cluster array from the file <name>.
  //
 
  TDirectory *savedir = gDirectory;                                                   

  if (fClusterArray) {
    fClusterArray->Delete();
  }
  else {
    fClusterArray = new TObjArray();
  }

  printf("AliTRDpid::ReadCluster -- ");
  printf("Open file %s\n",name);

  AliTRDtracker *tracker = new AliTRDtracker();
  TFile* file = TFile::Open(name);
  file->cd("Event0");
  TTree* tree = (TTree*) file->Get("TreeD");
  tracker->ReadClusters(fClusterArray,tree);
  file->Close();
  delete file;

  if (!fClusterArray) {
    printf("AliTRDpid::ReadCluster -- ");
    printf("Error reading the cluster array from file %s\n",name);
    return kFALSE;
  }

  delete tracker;

  savedir->cd();

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDpid::ReadTracks(const Char_t *name)
{
  //
  // Opens and reads the track array from the file <name>.
  //
 
  TDirectory *savedir = gDirectory;                                                   

  if (fTrackArray) {
    fTrackArray->Delete();
  }
  else {
    fTrackArray = new TObjArray();
  }

  TFile *file = (TFile *) gROOT->GetListOfFiles()->FindObject(name);
  if (!file) {
    printf("AliTRDpid::ReadTracks -- ");
    printf("Open file %s\n",name);
    file = new TFile(name);
    if (!file) {
      printf("AliTRDpid::ReadTracks -- ");
      printf("Cannot open file %s\n",name);
      return kFALSE;
    }
  }

  Char_t treeName[12];
  sprintf(treeName,"TreeT%d_TRD",fEvent);
  TTree   *trackTree   = (TTree *) file->Get(treeName);
  TBranch *trackBranch = trackTree->GetBranch("tracks");

  Int_t nEntry = ((Int_t) trackTree->GetEntries());
  for (Int_t iEntry = 0; iEntry < nEntry; iEntry++) {
    AliTRDtrack *track = new AliTRDtrack();
    trackBranch->SetAddress(&track);
    trackTree->GetEvent(iEntry);
    fTrackArray->AddLast(track);
  }

  file->Close();

  savedir->cd();

  return kTRUE;

}

//_____________________________________________________________________________
Int_t AliTRDpid::MCpid(const AliTRDtrack *t)
{
  //
  // Determines the pid of the MC track <t>.
  //

  // PDG codes
  const Int_t kPdgEl =  11;
  const Int_t kPdgPi = 211;

  AliTRDcluster *cluster;
  TParticle     *particle;

  Int_t   nClusterEl   = 0;
  Int_t   nClusterPi   = 0;
  Int_t   nClusterPure = 0;
  Int_t   nClusterMix  = 0;

  Float_t ratioEl      = 0;
  Float_t ratioPi      = 0;

  Int_t   ipid         = -1;

  if (!fClusterArray) {
    printf("AliTRDpid::MCpid -- ");
    printf("ClusterArray not defined. Initialize the PID object first\n");
    return -1;  
  }
  
  // Loop through all clusters associated to this track
  Int_t nCluster = t->GetNumberOfClusters();
  for (Int_t iCluster = 0; iCluster < nCluster; iCluster++) {

    // Get a cluster
    Int_t index = t->GetClusterIndex(iCluster);
    if (!(cluster = (AliTRDcluster *) fClusterArray->UncheckedAt(index))) {
      break;
    } 

    // Get the first two MC track indices
    Int_t track0 = cluster->GetLabel(0);
    Int_t track1 = cluster->GetLabel(1);

    // Check on the track index to find the right primaries
    if ((track0 >  fPIDindexMin) && 
        (track0 <= fPIDindexMax)) {

      // Check whether it is a pure cluster, i.e. not overlapping
      Bool_t accept = kTRUE;
      if ((fPIDpurePoints) && (track1 != -1)) {
        accept = kFALSE;
      }
      if (accept) {

        particle = gAlice->GetMCApp()->Particle(track0);
        if (particle->GetFirstMother() == -1) {
          switch (TMath::Abs(particle->GetPdgCode())) {
          case kPdgEl:
            nClusterEl++;
            break;
          case kPdgPi:
            nClusterPi++;
            break;
          };
          nClusterPure++;
        }

      }
      else {
 
        nClusterMix++;

      }

    }

  }
  
  if (nCluster) {
    ratioEl = ((Float_t) nClusterEl) / ((Float_t) nCluster);
    ratioPi = ((Float_t) nClusterPi) / ((Float_t) nCluster);
    if      (ratioEl > fPIDratioMin) {
      ipid = kElectron;
    }
    else if (ratioPi > fPIDratioMin) {
      ipid = kPion;
    }
  }
//   printf("AliTRDpid::MCpid -- ");
//   printf("nCluster = %d, nClusterEl = %d, nClusterPi = %d\n"
//   	 ,nCluster,nClusterEl,nClusterPi);
//   printf("AliTRDpid::MCpid -- ");
//   printf("nClusterPure = %d, nClusterMix = %d\n"
//   	 ,nClusterPure,nClusterMix);
//   printf("AliTRDpid::MCpid -- ");
//   printf("ratioEl = %f, ratioPi = %f, ipid = %d\n",ratioEl,ratioPi,ipid);

  return ipid;

}

//_____________________________________________________________________________
Int_t AliTRDpid::MCpid(const AliTRDtrack *t, Int_t *pdg
                                           , Int_t *nFound
                                           , Int_t *indices)
{
  //
  // Determines the pid of the MC track <t>.
  // Returns the number of different MC particles that contribute to this
  // track. <pdg> contains their PDG code, ordered by their frequency.
  // <nFound> contains the corresponding number of cluster that contain
  // contributions to this cluster.
  //

  const Int_t kNtrack = 3;
  const Int_t kNpart  = 200;

  AliTRDcluster *cluster;
  TParticle     *particle;

  Int_t nPart    = 0;
  Int_t iPart    = 0;
  Int_t iTrack   = 0;
  Int_t iCluster = 0;

  if (!fClusterArray) {
    printf("AliTRDpid::MCpid -- ");
    printf("ClusterArray not defined. Initialize the PID object first\n");
    return -1;  
  }
  
  Int_t sort[kNpart];
  Int_t pdgTmp[kNpart];
  Int_t nFoundTmp[kNpart];
  Int_t indicesTmp[kNpart];
  for (iPart = 0; iPart < kNpart; iPart++) {
    pdg[iPart]        = 0;
    nFound[iPart]     = 0;
    indices[iPart]    = 0;
    pdgTmp[iPart]     = 0;
    nFoundTmp[iPart]  = 0;
    indicesTmp[iPart] = 0;
    sort[iPart]       = 0;
  }

  // Loop through all clusters associated to this track
  Int_t nCluster = t->GetNumberOfClusters();
  for (iCluster = 0; iCluster < nCluster; iCluster++) {

    // Get a cluster
    Int_t index = t->GetClusterIndex(iCluster);
    if (!(cluster = (AliTRDcluster *) fClusterArray->UncheckedAt(index))) {
      break;
    } 

    // Get the MC track indices
    for (iTrack = 0; iTrack < kNtrack; iTrack++) {

      Int_t trackIndex = cluster->GetLabel(iTrack);
      if (trackIndex >= 0) {
        particle = gAlice->GetMCApp()->Particle(trackIndex);
        Int_t  pdgCode = particle->GetPdgCode(); 
        Bool_t newPart = kTRUE;
        for (iPart = 0; iPart < nPart; iPart++) {
         if (trackIndex == indicesTmp[iPart]) {
            nFoundTmp[iPart]++;
            newPart = kFALSE;
            break;
  	  }
        }
        if ((newPart) && (nPart < kNpart)) {
          indicesTmp[nPart] = trackIndex;
          pdgTmp[nPart]     = pdgCode;
          nFoundTmp[nPart]++;
          nPart++;
        }
      }

    }

  }

  // Sort the arrays by the frequency of the appearance of a MC track
  TMath::Sort(nPart,nFoundTmp,sort,kTRUE);
  for (iPart = 0; iPart < nPart; iPart++) {
    pdg[iPart]     = pdgTmp[sort[iPart]];
    nFound[iPart]  = nFoundTmp[sort[iPart]];
    indices[iPart] = indicesTmp[sort[iPart]];
//     printf(" iPart = %d, pdg = %d, nFound = %d, indices = %d\n"
// 	  ,iPart,pdg[iPart],nFound[iPart],indices[iPart]);
  }

  return nPart;

}

//_____________________________________________________________________________
Bool_t AliTRDpid::SumCharge(const AliTRDtrack *t
                           , Float_t *charge
                           , Int_t   *nCluster)
{
  //
  // Sums up the charge in each plane for track <t>.
  //

  Bool_t status = kTRUE;

  AliTRDcluster *cluster;

  const Int_t kNplane = AliTRDgeometry::Nplan();
  for (Int_t iPlane = 0; iPlane < kNplane; iPlane++) {
    charge[iPlane]   = 0.0;
    nCluster[iPlane] = 0;
  }

  if (!fClusterArray) {
    printf("AliTRDpid::SumCharge -- ");
    printf("ClusterArray not defined. Initialize the PID object first\n");
    return kFALSE;  
  }
  if (!fGeometry) {
    printf("AliTRDpid::SumCharge -- ");
    printf("TRD geometry not defined. Initialize the PID object first\n");
    return kFALSE;  
  }
  
  // Loop through all clusters associated to this track
  Int_t nClus = t->GetNumberOfClusters();
  for (Int_t iClus = 0; iClus < nClus; iClus++) {

    // Get a cluster
    Int_t index = t->GetClusterIndex(iClus);
    if (!(cluster = (AliTRDcluster *) fClusterArray->UncheckedAt(index))) {
      status = kFALSE;
      break;
    } 

    if (fThreePadOnly) {
      if ((!cluster->From2pad()) &&
          (!cluster->From3pad())) continue;
    }

    // Sum up the charge
    Int_t plane    = fGeometry->GetPlane(cluster->GetDetector());
    //charge[plane] += cluster->GetQ();
    // Corrected charge
    charge[plane] += t->GetClusterdQdl(iClus);
    nCluster[plane]++;

  }
//   printf("AliTRDpid::SumCharge -- ");
//   printf("charge = %d, %d, %d, %d, %d, %d\n"
//   	 ,(Int_t) charge[0]
//          ,(Int_t) charge[1]
//          ,(Int_t) charge[2]
//          ,(Int_t) charge[3]
//          ,(Int_t) charge[4]
//          ,(Int_t) charge[5]);
//   printf("AliTRDpid::SumCharge -- ");
//   printf("nCluster = %d, %d, %d, %d, %d, %d\n"
//          ,nCluster[0]
//          ,nCluster[1]
//          ,nCluster[2]
//          ,nCluster[3]
//          ,nCluster[4]
//          ,nCluster[5]);

  return status;

}
