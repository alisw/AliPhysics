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

/*
$Log$
*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//   The TRD particle identification class                                   //
//                                                                           //
//   Its main purposes are:                                                  //
//      - Creation and bookkeeping of the propability distributions          //
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
#include "AliTRD.h"
#include "AliTRDpid.h"
#include "AliTRDcluster.h"
#include "AliTRDtrack.h"
#include "AliTRDtracker.h"
#include "AliTRDgeometry.h"

ClassImp(AliTRDpid)

//_____________________________________________________________________________
AliTRDpid::AliTRDpid():TNamed()
{
  //
  // AliTRDpid default constructor
  // 

  fNMom   = 0;
  fMinMom = 0;
  fMaxMom = 0;
  fWidMom = 0;

  fQHist        = NULL;
  fLQHist       = NULL;
  fTrackArray   = NULL;
  fClusterArray = NULL;
  fGeometry     = NULL;

}

//_____________________________________________________________________________
AliTRDpid::AliTRDpid(const char* name, const char* title):TNamed(name,title)
{
  //
  // AliTRDpid constructor
  //

  fNMom   = 0;
  fMinMom = 0;
  fMaxMom = 0;
  fWidMom = 0;

  fQHist        = NULL;
  fLQHist       = NULL;
  fTrackArray   = NULL;
  fClusterArray = NULL;
  fGeometry     = NULL;

  Init();

}

//_____________________________________________________________________________
AliTRDpid::AliTRDpid(const AliTRDpid &p)
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
  }

  if (fTrackArray) {
    fTrackArray->Delete();
    delete fTrackArray;
  }

  if (fQHist) {
    fQHist->Delete();
    delete fQHist;
  }

  if (fLQHist) {
    fLQHist->Delete();
    delete fLQHist;
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
void AliTRDpid::Copy(TObject &p)
{
  //
  // Copy function
  //

  fQHist->Copy(*((AliTRDpid &) p).fQHist);
  fLQHist->Copy(*((AliTRDpid &) p).fLQHist);

  ((AliTRDpid &) p).fTrackArray   = NULL;    
  ((AliTRDpid &) p).fClusterArray = NULL;    
  ((AliTRDpid &) p).fGeometry     = NULL;    

  ((AliTRDpid &) p).fNMom   = fNMom;
  ((AliTRDpid &) p).fMinMom = fMinMom;
  ((AliTRDpid &) p).fMaxMom = fMaxMom;
  ((AliTRDpid &) p).fWidMom = fWidMom;

}

//_____________________________________________________________________________
Bool_t AliTRDpid::Init()
{
  //
  // Initializes the PID object 
  //

  fClusterArray = new TObjArray();
  fTrackArray   = new TObjArray();

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDpid::AssignLQ(TObjArray *tarray)
{
  //
  // Assigns the e / pi Q-likelihood to all tracks in the array
  //

  Bool_t status = kTRUE;

  AliTRDtrack *track;

  TIter nextTrack(tarray);  
  while ((track = (AliTRDtrack *) nextTrack())) {
    if (!AssignLQ(track)) {
      status = kFALSE;
      break;
    }
  }

  return status;

}

//_____________________________________________________________________________
Bool_t AliTRDpid::AssignLQ(AliTRDtrack *t)
{
  //
  // Assigns the e / pi Q-likelihood to a given track
  //

  const Int_t kNplane = AliTRDgeometry::Nplan();
  Float_t charge[kNplane];

  // Calculate the total charge in each plane
  if (!SumCharge(t,charge)) return kFALSE;

  // Assign the likelihoods 
  t->SetLikelihoodPion(LQPion(charge));
  t->SetLikelihoodElectron(LQElectron(charge));

  return kTRUE;  

}

//_____________________________________________________________________________
Bool_t AliTRDpid::CreateHistograms(const Int_t nmom
                                 , const Float_t minmom
                                 , const Float_t maxmom)
{
  //
  // Creates the likelihood histograms
  //

  Int_t imom;
  Int_t ipid;
  Int_t ipla;

  const Int_t kNpla = AliTRDgeometry::Nplan();

  fNMom   = nmom;
  fMinMom = minmom;
  fMaxMom = maxmom;
  fWidMom = (maxmom - minmom) / ((Float_t) nmom);

  // The L-Q distributions
  fLQHist = new TObjArray(kNpid * nmom);
  for (imom = 0; imom < nmom;  imom++) {
    for (ipid = 0; ipid < kNpid; ipid++) {
      Int_t index = GetIndexLQ(imom,ipid);
      Char_t name[10];
      Char_t title[80];
      sprintf(name ,"hLQ%03d",index);
      if (ipid == kElectron) {
        sprintf(title,"L-Q electron p-bin %03d",imom);
      }
      else {
        sprintf(title,"L-Q pion p-bin %03d",imom);
      }
      TH1F *hTmp  = new TH1F(name,title,416,-0.02,1.02);
      fLQHist->AddAt(hTmp,index);
    }
  }

  // The Q-distributions
  fQHist = new TObjArray(kNpla * kNpid * nmom);
  for (imom = 0; imom < nmom;  imom++) {
    for (ipid = 0; ipid < kNpid; ipid++) {
      for (ipla = 0; ipla < kNpla; ipla++) {
        Int_t index = GetIndexQ(imom,ipla,ipid);
        Char_t name[10];
        Char_t title[80];
        sprintf(name ,"hQ%03d",index);
        if (ipid == kElectron) {
          sprintf(title,"Q electron p-bin %03d plane %d",imom,ipla);
        }
        else {
          sprintf(title,"Q pion p-bin %03d plane %d",imom,ipla);
        }
        TH1F *hTmp  = new TH1F(name,title,100,0.0,5000.0);
        fQHist->AddAt(hTmp,index);
      }
    }
  }

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDpid::FillQspectra()
{
  //
  // Fills the Q-spectra histograms.
  //
  
  Bool_t status = kTRUE;

  AliTRDtrack *track;

  TIter nextTrack(fTrackArray);
  while ((track = (AliTRDtrack *) nextTrack())) {
    if (!FillQspectra(track)) {
      status = kFALSE;
      break;
    }
  }

  return status;

}

//_____________________________________________________________________________
Bool_t AliTRDpid::FillQspectra(const AliTRDtrack *t)
{
  //
  // Fills the Q-spectra histograms with track <t>.
  //

  const Int_t kNpla = AliTRDgeometry::Nplan();

  if (isnan(t->GetP())) return kTRUE;

  printf("----------------------------------------------------------\n");
  Float_t charge[kNpla];
  Float_t mom  = t->GetP();
  Int_t   ipid = Pid(t);

  if (!SumCharge(t,charge)) return kFALSE;

  printf(" ipid   = %d\n",ipid);
  printf(" mom    = %f\n",mom);
  printf(" charge = %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f\n"
        ,charge[0],charge[1],charge[2],charge[3],charge[4],charge[5]);

  for (Int_t ipla = 0; ipla < kNpla; ipla++) {
    Int_t index = GetIndexQ(mom,ipla,ipid);    
    if (index > -1) {
      TH1F *hTmp = (TH1F *) fQHist->UncheckedAt(index);
      hTmp->Fill(charge[ipla]);
    }
  }  

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDpid::Open(const Char_t *name, Int_t event)
{
  //
  // Opens and reads the kine tree, the geometry, the cluster array.
  // and the track array from the file <name>
  //

  Bool_t status = kTRUE;

  status = ReadKine(name,event);
  status = ReadCluster(name);
  status = ReadTracks(name);

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
  // and the track array from the file <nametracks>
  //

  Bool_t status = kTRUE;

  status = ReadKine(namekine,event);
  status = ReadCluster(namecluster);
  status = ReadTracks(nametracks);

  return status;

}

//_____________________________________________________________________________
Bool_t AliTRDpid::ReadKine(const Char_t *name, Int_t event)
{
  //
  // Opens and reads the kine tree and the geometry from the file <name>.
  //

  TFile *file = (TFile *) gROOT->GetListOfFiles()->FindObject(name);
  if (!file) {
    printf("AliTRDpid::ReadKine -- ");
    printf("Open file %s\n",name);
    file = new TFile(name);
    if (!file) {
      printf("AliTRDpid::ReadKine -- ");
      printf("Cannot open file %s\n",name);
      return kFALSE;
    }
  }

  gAlice = (AliRun *) file->Get("gAlice");
  if (!gAlice) {
    printf("AliTRDpid::ReadKine -- ");
    printf("No AliRun object found\n");    
    return kFALSE;
  }
  gAlice->GetEvent(event);

  AliTRD *trd = (AliTRD *) gAlice->GetDetector("TRD");
  if (!trd) {
    printf("AliTRDpid::ReadKine -- ");
    printf("No TRD object found\n");    
    return kFALSE;
  }

  fGeometry = trd->GetGeometry();
  if (!fGeometry) {
    printf("AliTRDpid::ReadKine -- ");
    printf("No TRD geometry found\n");
    return kFALSE;
  }

  file->Close(); 

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDpid::ReadCluster(const Char_t *name)
{
  //
  // Opens and reads the cluster array from the file <name>.
  //

  fClusterArray->Delete();

  printf("AliTRDpid::ReadCluster -- ");
  printf("Open file %s\n",name);

  AliTRDtracker *tracker = new AliTRDtracker("dummy","dummy");
  tracker->ReadClusters(fClusterArray,name);

  if (!fClusterArray) {
    printf("AliTRDpid::ReadCluster -- ");
    printf("Error reading the cluster array from file %s\n",name);
    return kFALSE;
  }

  delete tracker;

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDpid::ReadTracks(const Char_t *name)
{
  //
  // Opens and reads the track array from the file <name>.
  //

  fTrackArray->Delete();

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

  TTree   *trackTree   = (TTree *) file->Get("TreeT");
  TBranch *trackBranch = trackTree->GetBranch("tracks");

  Int_t nEntry = ((Int_t) trackTree->GetEntries());
  for (Int_t iEntry = 0; iEntry < nEntry; iEntry++) {
    AliTRDtrack *track = new AliTRDtrack();
    trackBranch->SetAddress(&track);
    trackTree->GetEvent(iEntry);
    fTrackArray->AddLast(track);
  }

  file->Close();

  return kTRUE;

}

//_____________________________________________________________________________
Int_t AliTRDpid::GetIndexQ(const Float_t mom, const Int_t ipla, const Int_t ipid)
{
  //
  // Returns the Q-spectrum histogram index
  //

  if ((mom < fMinMom) || (mom > fMaxMom))  return -1;
  Int_t imom = ((Int_t) ((mom - fMinMom) / fWidMom));
  return GetIndexQ(imom,ipla,ipid);

}

//_____________________________________________________________________________
Int_t AliTRDpid::GetIndexQ(const Int_t imom, const Int_t ipla, const Int_t ipid)
{
  //
  // Returns the Q-spectrum histogram index
  //

  const Int_t kNpla = AliTRDgeometry::Nplan();
  if ((ipid < 0) || (ipid >= kNpid)) return -1;
  return imom * (kNpid * kNpla) + ipla * kNpid + ipid; 

}

//_____________________________________________________________________________
Int_t AliTRDpid::GetIndexLQ(const Float_t mom, const Int_t ipid)
{
  //
  // Returns the Q-likelihood histogram index
  //

  if ((mom < fMinMom) || (mom > fMaxMom))  return -1;
  Int_t imom = ((Int_t) ((mom - fMinMom) / fWidMom));
  return GetIndexLQ(imom,ipid);

}

//_____________________________________________________________________________
Int_t AliTRDpid::GetIndexLQ(const Int_t imom, const Int_t ipid)
{
  //
  // Returns the Q-likelihood histogram index
  //

  if ((ipid < 0) || (ipid >= kNpid)) return -1;
  return imom * kNpid + ipid; 

}

//_____________________________________________________________________________
Int_t AliTRDpid::Pid(const AliTRDtrack *t)
{
  //
  // Determines the pid of the track <t>
  // For a given track to be assigned an electron or pion pid it requires
  // that more than a fraction of 0.7 of all associated cluster are 'pure'
  // clusters -- meaning to have only contribution from one particle --
  // of the specific particle type.
  //

  // PDG codes
  const Int_t kPdgEl =  11;
  const Int_t kPdgPi = 211;

  // Minimum fraction of cluster from one particle
  const Float_t kRatioMin = 0.7;

  AliTRDcluster *cluster;
  TParticle     *particle;

  Int_t   nClusterEl = 0;
  Int_t   nClusterPi = 0;
  Int_t   nClusterNo = 0;

  Float_t ratioEl    = 0;
  Float_t ratioPi    = 0;

  Int_t   ipid       = -1;

  if (!fClusterArray) {
    printf("AliTRDpid::Pid -- ");
    printf("ClusterArray not defined. Initialize the PID object first\n");
    return -1;  
  }
  
  // Loop through all clusters associated to this track
  Int_t nCluster = t->GetNclusters();
  for (Int_t iCluster = 0; iCluster < nCluster; iCluster++) {

    // Get a cluster
    Int_t index = t->GetClusterIndex(iCluster);
    if (!(cluster = (AliTRDcluster *) fClusterArray->UncheckedAt(index))) {
      break;
    } 

    // Get the first two MC track indices
    Int_t track0 = cluster->GetTrackIndex(0);
    Int_t track1 = cluster->GetTrackIndex(1);
    printf(" track0 = %d, track1 = %d\n",track0,track1);

    // Take only 'pure' cluster 
    if ((track0 > -1) && (track1 == -1)) {

      particle = gAlice->Particle(track0);
      switch (TMath::Abs(particle->GetPdgCode())) {
      case kPdgEl:
        nClusterEl++;
        break;
      case kPdgPi:
        nClusterPi++;
        break;
      default:
        nClusterNo++;
        break;
      };

    }

  }
  
  if (nCluster) {
    ratioEl = ((Float_t) nClusterEl) / ((Float_t) nCluster);
    ratioPi = ((Float_t) nClusterPi) / ((Float_t) nCluster);
    if      (ratioEl > kRatioMin) {
      ipid = kElectron;
    }
    else if (ratioPi > kRatioMin) {
      ipid = kPion;
    }
  }
  printf(" nCluster = %d, nClusterEl = %d, nClusterPi = %d, nClusterNo = %d\n"
	 ,nCluster,nClusterEl,nClusterPi,nClusterNo);
  printf(" ratioEl = %f, ratioPi = %f\n",ratioEl,ratioPi);

  return ipid;

}

//_____________________________________________________________________________
Bool_t AliTRDpid::SumCharge(const AliTRDtrack *t, Float_t *charge)
{
  //
  // Sums up the charge in each plane for track <t>.
  //

  Bool_t status = kTRUE;

  AliTRDcluster *cluster;

  const Int_t kNplane = AliTRDgeometry::Nplan();
  for (Int_t iPlane = 0; iPlane < kNplane; iPlane++) {
    charge[iPlane] = 0.0;
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
  Int_t nCluster = t->GetNclusters();
  for (Int_t iCluster = 0; iCluster < nCluster; iCluster++) {

    // Get a cluster
    Int_t index = t->GetClusterIndex(iCluster);
    if (!(cluster = (AliTRDcluster *) fClusterArray->UncheckedAt(index))) {
      status = kFALSE;
      break;
    } 

    // Sum up the charge
    Int_t plane    = fGeometry->GetPlane(cluster->GetDetector());
    charge[plane] += cluster->GetQ();

  }

  return status;

}

//_____________________________________________________________________________
Float_t AliTRDpid::LQElectron(const Float_t *charge)
{
  //
  // Returns the Q-likelihood to be a electron
  //

  return 1;

}

//_____________________________________________________________________________
Float_t AliTRDpid::LQPion(const Float_t *charge)
{
  //
  // Returns the Q-likelihood to be a pion
  //

  return 0;

}

