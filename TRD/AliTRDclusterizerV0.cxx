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
// TRD cluster finder for the fast simulator. It takes the hits from the     //
// fast simulator (one hit per plane) and transforms them                    //
// into cluster, by applying position smearing and merging                   //
// of nearby cluster. The searing is done uniformly in z-direction           //
// over the length of a readout pad. In rphi-direction a Gaussian            //
// smearing is applied with a sigma given by fRphiSigma.                     //
// Clusters are considered as overlapping when they are closer in            //
// rphi-direction than the value defined in fRphiDist.                       //
// Use the macro fastClusterCreate.C to create the cluster.                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TRandom.h>

#include "AliTRDclusterizerV0.h"
#include "AliTRDconst.h"
#include "AliTRDgeometry.h"
#include "AliTRDrecPoint.h"

ClassImp(AliTRDclusterizerV0)

//_____________________________________________________________________________
AliTRDclusterizerV0::AliTRDclusterizerV0():AliTRDclusterizer()
{
  //
  // AliTRDclusterizerV0 default constructor
  //

}

//_____________________________________________________________________________
AliTRDclusterizerV0::AliTRDclusterizerV0(const Text_t* name, const Text_t* title)
                    :AliTRDclusterizer(name,title)
{
  //
  // AliTRDclusterizerV0 default constructor
  //

  Init();

}

//_____________________________________________________________________________
AliTRDclusterizerV0::~AliTRDclusterizerV0()
{

}

//_____________________________________________________________________________
void AliTRDclusterizerV0::Init()
{
  //
  // Initializes the cluster finder
  //

  // Position resolution in rphi-direction
  fRphiSigma  = 0.02;
  // Minimum distance of non-overlapping cluster
  fRphiDist   = 1.0;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizerV0::MakeCluster()
{
  //
  // Generates the cluster
  //

  // Get the pointer to the detector class and check for version 1
  AliTRD *TRD = (AliTRD*) gAlice->GetDetector("TRD");
  if (TRD->IsVersion() != 0) {
    printf("AliTRDclusterizerV0::MakeCluster -- ");
    printf("TRD must be version 0 (fast simulator).\n");
    return kFALSE; 
  }

  // Get the geometry
  AliTRDgeometry *Geo = TRD->GetGeometry();
  
  printf("AliTRDclusterizerV0::MakeCluster -- ");
  printf("Start creating cluster.\n");

  Int_t nBytes = 0;

  AliTRDhit      *Hit;
  
  // Get the pointer to the hit tree
  TTree *HitTree     = gAlice->TreeH();
  // Get the pointer to the reconstruction tree
  TTree *ClusterTree = gAlice->TreeR();

  TObjArray *Chamber = new TObjArray();

  // Get the number of entries in the hit tree
  // (Number of primary particles creating a hit somewhere)
  Int_t nTrack = (Int_t) HitTree->GetEntries();

  // Loop through all the chambers
  for (Int_t icham = 0; icham < kNcham; icham++) {
    for (Int_t iplan = 0; iplan < kNplan; iplan++) {
      for (Int_t isect = 0; isect < kNsect; isect++) {

        Int_t   nColMax    = Geo->GetColMax(iplan);
        Float_t row0       = Geo->GetRow0(iplan,icham,isect);
        Float_t col0       = Geo->GetCol0(iplan);

        Float_t rowPadSize = Geo->GetRowPadSize();
        Float_t colPadSize = Geo->GetColPadSize();

        // Loop through all entries in the tree
        for (Int_t iTrack = 0; iTrack < nTrack; iTrack++) {

          gAlice->ResetHits();
          nBytes += HitTree->GetEvent(iTrack);

          // Get the number of hits in the TRD created by this particle
          Int_t nHit = TRD->Hits()->GetEntriesFast();

          // Loop through the TRD hits  
          for (Int_t iHit = 0; iHit < nHit; iHit++) {

            if (!(Hit = (AliTRDhit *) TRD->Hits()->UncheckedAt(iHit))) 
              continue;

            Float_t pos[3];
                    pos[0]   = Hit->fX;
                    pos[1]   = Hit->fY;
                    pos[2]   = Hit->fZ;
            Int_t   track    = Hit->fTrack;
            Int_t   detector = Hit->fDetector;
            Int_t   plane    = Geo->GetPlane(detector);
            Int_t   sector   = Geo->GetSector(detector);
            Int_t   chamber  = Geo->GetChamber(detector);        

            if ((sector  != isect) ||
                (plane   != iplan) ||
                (chamber != icham)) 
              continue;

            // Rotate the sectors on top of each other
            Float_t rot[3];
            Geo->Rotate(detector,pos,rot);

            // Add this cluster to the temporary cluster-array for this chamber
            Int_t   tracks[3];
            tracks[0] = track;
            Int_t   clusters[2];
            clusters[0] = detector;
            clusters[1] = 0;
            Float_t position[3];
            position[0] = rot[2];
            position[1] = rot[1];
            position[2] = rot[0];
	    AliTRDcluster *Cluster = new AliTRDcluster(tracks,clusters,0,position);
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
            if ((TMath::Abs(x1 - x2) < rowPadSize) ||
                (TMath::Abs(y1 - y2) <  fRphiDist)) {
              if (iSave == nSave) {
                printf("AliTRDclusterizerV0::MakeCluster -- ");
                printf("Boundary error: iSave = %d, nSave = %d.\n"
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
              AliTRDcluster *Cluster2 = 
                (AliTRDcluster *) Chamber->UncheckedAt(idxSave[iMerge]);
              xMerge += Cluster2->fX;
              yMerge += Cluster2->fY;
              Cluster2->fZ = 0;            // Mark merged cluster
	    }
            xMerge /= (iSave + 1);
            yMerge /= (iSave + 1);
          }

          Float_t smear[3];

          // The position smearing in z-direction (uniform over pad width)            
          Int_t row = (Int_t) ((xMerge - row0) / rowPadSize);
          smear[0]  = (row + gRandom->Rndm()) * rowPadSize + row0;

          // The position smearing in rphi-direction (Gaussian)
          smear[1] = 0;
          do
            smear[1] = gRandom->Gaus(yMerge,fRphiSigma);
          while ((smear[1] < col0                        ) ||
                 (smear[1] > col0 + nColMax * colPadSize));

          // Time direction stays unchanged
          smear[2] = z1;

          // Add the smeared cluster to the output array 
          Int_t detector  = Cluster1->fDetector;
          Int_t digits[3] = {0};
          TRD->AddRecPoint(smear,digits,detector,0.0);

	}

        // Clear the temporary cluster-array and delete the cluster
        Chamber->Delete();

      }
    }
  }

  printf("AliTRDclusterizerV0::MakeCluster -- ");
  printf("Found %d points.\n",TRD->RecPoints()->GetEntries());
  printf("AliTRDclusterizerV0::MakeCluster -- ");
  printf("Fill the cluster tree.\n");
  ClusterTree->Fill();

  return kTRUE;

}
