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
Revision 1.12  1999/11/02 16:35:56  fca
New version of TRD introduced

Revision 1.11  1999/11/01 20:41:51  fca
Added protections against using the wrong version of FRAME

Revision 1.10  1999/09/29 09:24:35  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Transition Radiation Detector version 0 -- fast simulator                //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliTRDfullClass.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TRandom.h>
#include <TVector.h>

#include "AliTRDv0.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliConst.h"
  
ClassImp(AliTRDv0)

//_____________________________________________________________________________
AliTRDv0::AliTRDv0(const char *name, const char *title) 
         :AliTRD(name, title) 
{
  //
  // Standard constructor for Transition Radiation Detector version 0
  //

  fIdSens     = 0;
  fHitsOn     = 0;

  fIdChamber1 = 0;
  fIdChamber2 = 0;
  fIdChamber3 = 0;

  fRphiSigma  = 0;
  fRphiDist   = 0;

}
 
//_____________________________________________________________________________
void AliTRDv0::CreateGeometry()
{
  //
  // Create the GEANT geometry for the Transition Radiation Detector - Version 0
  // This version covers the full azimuth. 
  //

  // Check that FRAME is there otherwise we have no place where to put the TRD
  AliModule* FRAME = gAlice->GetModule("FRAME");
  if (!FRAME) return;

  // Define the chambers
  AliTRD::CreateGeometry();

}

//_____________________________________________________________________________
void AliTRDv0::CreateMaterials()
{
  //
  // Create materials for the Transition Radiation Detector
  //

  AliTRD::CreateMaterials();

}

//_____________________________________________________________________________
void AliTRDv0::Hits2Clusters()
{
  // A simple cluster generator. It takes the hits from the
  // fast simulator (one hit per plane) and transforms them
  // into cluster, by applying position smearing and merging
  // of nearby cluster. The searing is done uniformly in z-direction
  // over the length of a readout pad. In rphi-direction a Gaussian
  // smearing is applied with a sigma given by fRphiSigma.
  // Clusters are considered as overlapping when they are closer in
  // rphi-direction than the value defined in fRphiDist.
  // Use the macro fastClusterCreate.C to create the cluster.

  printf("AliTRDv0::Hits2Clusters -- Start creating cluster\n");

  Int_t nBytes = 0;

  AliTRDhit *TRDhit;
  
  // Get the pointer to the hit tree
  TTree *HitTree     = gAlice->TreeH();
  // Get the pointer to the reconstruction tree
  TTree *ClusterTree = gAlice->TreeD();

  TObjArray *Chamber = new TObjArray();

  // Get the number of entries in the hit tree
  // (Number of primary particles creating a hit somewhere)
  Int_t nTrack = (Int_t) HitTree->GetEntries();

  // Loop through all the chambers
  for (Int_t icham = 0; icham < kNcham; icham++) {
    for (Int_t iplan = 0; iplan < kNplan; iplan++) {
      for (Int_t isect = 0; isect < kNsect; isect++) {

        // Loop through all entries in the tree
        for (Int_t iTrack = 0; iTrack < nTrack; iTrack++) {

          gAlice->ResetHits();
          nBytes += HitTree->GetEvent(iTrack);

          // Get the number of hits in the TRD created by this particle
          Int_t nHit = fHits->GetEntriesFast();

          // Loop through the TRD hits  
          for (Int_t iHit = 0; iHit < nHit; iHit++) {

            if (!(TRDhit = (AliTRDhit *) fHits->UncheckedAt(iHit))) 
              continue;

            Float_t x       = TRDhit->fX;
            Float_t y       = TRDhit->fY;
            Float_t z       = TRDhit->fZ;
            Int_t   track   = TRDhit->fTrack;
            Int_t   plane   = TRDhit->fPlane;
            Int_t   sector  = TRDhit->fSector;
            Int_t   chamber = TRDhit->fChamber;        

            if ((sector  != isect+1) ||
                (plane   != iplan+1) ||
                (chamber != icham+1)) 
              continue;

            // Rotate the sectors on top of each other
            Float_t phi  = 2.0 * kPI /  (Float_t) kNsect 
                               * ((Float_t) sector - 0.5);
            Float_t xRot = -x * TMath::Cos(phi) + y * TMath::Sin(phi);
            Float_t yRot =  x * TMath::Sin(phi) + y * TMath::Cos(phi);
            Float_t zRot =  z;

            // Add this cluster to the temporary cluster-array for this chamber
            Int_t   tracks[3];
            tracks[0] = track;
            Int_t   clusters[5];
            clusters[0] = sector;
            clusters[1] = chamber;
            clusters[2] = plane;
            clusters[3] = 0;
            clusters[4] = 0;
            Float_t position[3];
            position[0] = zRot;
            position[1] = yRot;
            position[2] = xRot;
	    AliTRDcluster *Cluster = new AliTRDcluster(tracks,clusters,position);
            Chamber->Add(Cluster);

	  }

	}
  
        // Loop through the temporary cluster-array
        for (Int_t iClus1 = 0; iClus1 < Chamber->GetEntries(); iClus1++) {

          AliTRDcluster *Cluster1 = (AliTRDcluster *) Chamber->UncheckedAt(iClus1);
          Float_t x1 = Cluster1->fX;
          Float_t y1 = Cluster1->fY;
          Float_t z1 = Cluster1->fZ;

          if (!(z1)) continue;             // Skip marked cluster  

          const Int_t nSave = 2;
          Int_t idxSave[nSave];
          Int_t iSave = 0;

          Int_t tracks[3];
          tracks[0] = Cluster1->fTracks[0];

          // Check the other cluster to see, whether there are close ones
          for (Int_t iClus2 = iClus1 + 1; iClus2 < Chamber->GetEntries(); iClus2++) {
            AliTRDcluster *Cluster2 = (AliTRDcluster *) Chamber->UncheckedAt(iClus2);
            Float_t x2 = Cluster2->fX;
            Float_t y2 = Cluster2->fY;
            if ((TMath::Abs(x1 - x2) < fRowPadSize) ||
                (TMath::Abs(y1 - y2) <   fRphiDist)) {
              if (iSave == nSave) { 
                printf("AliTRDv0::Hits2Clusters -- Boundary error: iSave = %d, nSave = %d\n"
                      ,iSave,nSave);
	      }
              else {                
                idxSave[iSave]  = iClus2;
                tracks[iSave+1] = Cluster2->fTracks[0];
	      }
              iSave++;
	    }
	  }
     
          // Merge close cluster
          Float_t yMerge = y1;
          Float_t xMerge = x1;
          if (iSave) {
            for (Int_t iMerge = 0; iMerge < iSave; iMerge++) {
              AliTRDcluster *Cluster2 = (AliTRDcluster *) Chamber->UncheckedAt(idxSave[iMerge]);
              xMerge += Cluster2->fX;
              yMerge += Cluster2->fY;
              Cluster2->fZ = 0;            // Mark merged cluster
	    }
            xMerge /= (iSave + 1);
            yMerge /= (iSave + 1);
          }

          // The position smearing in z-direction (uniform over pad width)
          Int_t row = (Int_t) ((xMerge - fRow0[iplan][icham][isect]) / fRowPadSize);
          Float_t xSmear = (row + gRandom->Rndm()) * fRowPadSize 
                         + fRow0[iplan][icham][isect];

          // The position smearing in rphi-direction (Gaussian)
          Float_t ySmear = 0;
          do
            ySmear = gRandom->Gaus(yMerge,fRphiSigma);
          while ((ySmear < fCol0[iplan])                               ||
                 (ySmear > fCol0[iplan] + fColMax[iplan] * fColPadSize));

          // Time direction stays unchanged
          Float_t zSmear = z1;

          Int_t   clusters[5];
          clusters[0] = Cluster1->fSector;
          clusters[1] = Cluster1->fChamber;
          clusters[2] = Cluster1->fPlane;
          clusters[3] = 0;
          clusters[4] = 0;
          Float_t position[3];
          // Rotate the sectors back into their real position
          Float_t phi = 2*kPI / kNsect * ((Float_t) Cluster1->fSector - 0.5);
          position[0] = -zSmear * TMath::Cos(phi) + ySmear * TMath::Sin(phi);
          position[1] =  zSmear * TMath::Sin(phi) + ySmear * TMath::Cos(phi);
          position[2] =  xSmear;

          // Add the smeared cluster to the output array 
          AddCluster(tracks,clusters,position);

	}

        // Clear the temporary cluster-array and delete the cluster
        Chamber->Delete();

      }
    }
  }

  printf("AliTRDv0::Hits2Clusters -- Found %d cluster\n",fClusters->GetEntries());
  printf("AliTRDv0::Hits2Clusters -- Fill the cluster tree\n");
  ClusterTree->Fill();

}

//_____________________________________________________________________________
void AliTRDv0::Init() 
{
  //
  // Initialise Transition Radiation Detector after geometry is built
  //

  AliTRD::Init();

  // Identifier of the sensitive volume (amplification region)
  fIdSens     = gMC->VolId("UL06");

  // Identifier of the TRD-driftchambers
  fIdChamber1 = gMC->VolId("UCIO");
  fIdChamber2 = gMC->VolId("UCIM");
  fIdChamber3 = gMC->VolId("UCII");

  // Parameter for Hits2Cluster

  // Position resolution in rphi-direction
  fRphiSigma  = 0.02;
  // Minimum distance of non-overlapping cluster
  fRphiDist   = 1.0;

  printf("          Fast simulator\n");
  for (Int_t i = 0; i < 80; i++) printf("*");
  printf("\n");
  
}

//_____________________________________________________________________________
void AliTRDv0::StepManager()
{
  //
  // Procedure called at every step in the TRD
  // Fast simulator. If switched on, a hit is produced when a track
  // crosses the border between amplification region and pad plane.
  //

  Int_t   vol[3]; 
  Int_t   iIdSens, icSens; 
  Int_t   iIdChamber, icChamber;

  Float_t hits[4];

  TLorentzVector p;
  TClonesArray  &lhits = *fHits;

  // Writing out hits enabled?
  if (!(fHitsOn)) return;

  // Use only charged tracks and count them only once per volume
  if (gMC->TrackCharge()    && 
      gMC->IsTrackExiting()) {
    
    // Check on sensitive volume
    iIdSens = gMC->CurrentVolID(icSens);
    if (iIdSens == fIdSens) { 

      gMC->TrackPosition(p);
      for (Int_t i = 0; i < 3; i++) hits[i] = p[i];
      // No charge created
      hits[3] = 0;

      // The sector number
      Float_t phi = hits[1] != 0 ? kRaddeg*TMath::ATan2(hits[0],hits[1]) : (hits[0] > 0 ? 180. : 0.);
      vol[0] = ((Int_t) (phi / 20)) + 1;

      // The chamber number 
      //   1: outer left
      //   2: middle left
      //   3: inner
      //   4: middle right
      //   5: outer right
      iIdChamber = gMC->CurrentVolOffID(1,icChamber);
      if      (iIdChamber == fIdChamber1)
        vol[1] = (hits[2] < 0 ? 1 : 5);
      else if (iIdChamber == fIdChamber2)       
        vol[1] = (hits[2] < 0 ? 2 : 4);
      else if (iIdChamber == fIdChamber3)       
        vol[1] = 3;

      // The plane number
      vol[2] = icChamber - TMath::Nint((Float_t) (icChamber / 7)) * 6;

      new(lhits[fNhits++]) AliTRDhit(fIshunt,gAlice->CurrentTrack(),vol,hits);

    }

  }  

}
