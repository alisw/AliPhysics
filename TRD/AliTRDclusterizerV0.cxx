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
Revision 1.3  2000/06/07 16:27:01  cblume
Try to remove compiler warnings on Sun and HP

Revision 1.2  2000/05/08 16:17:27  cblume
Merge TRD-develop

Revision 1.1.4.1  2000/05/08 15:08:41  cblume
Replace AliTRDcluster by AliTRDrecPoint

Revision 1.1  2000/02/28 18:58:33  cblume
Add new TRD classes

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
  //
  // AliTRDclusterizerV0 destructor
  //

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
  AliTRD *trd = (AliTRD*) gAlice->GetDetector("TRD");
  if (trd->IsVersion() != 0) {
    printf("AliTRDclusterizerV0::MakeCluster -- ");
    printf("TRD must be version 0 (fast simulator).\n");
    return kFALSE; 
  }

  // Get the geometry
  AliTRDgeometry *geo = trd->GetGeometry();
  
  printf("AliTRDclusterizerV0::MakeCluster -- ");
  printf("Start creating cluster.\n");

  Int_t nBytes = 0;

  AliTRDhit *hit;
  
  // Get the pointer to the hit tree
  TTree     *hitTree      = gAlice->TreeH();
  // Get the pointer to the reconstruction tree
  TTree     *clusterTree  = gAlice->TreeR();

  TObjArray *chamberArray = new TObjArray();

  // Get the number of entries in the hit tree
  // (Number of primary particles creating a hit somewhere)
  Int_t nTrack = (Int_t) hitTree->GetEntries();

  // Loop through all the chambers
  for (Int_t icham = 0; icham < kNcham; icham++) {
    for (Int_t iplan = 0; iplan < kNplan; iplan++) {
      for (Int_t isect = 0; isect < kNsect; isect++) {

        Int_t   nColMax     = geo->GetColMax(iplan);
        Float_t row0        = geo->GetRow0(iplan,icham,isect);
        Float_t col0        = geo->GetCol0(iplan);
        Float_t time0       = geo->GetTime0(iplan);

        Float_t rowPadSize  = geo->GetRowPadSize();
        Float_t colPadSize  = geo->GetColPadSize();
        Float_t timeBinSize = geo->GetTimeBinSize();

        // Loop through all entries in the tree
        for (Int_t iTrack = 0; iTrack < nTrack; iTrack++) {

          gAlice->ResetHits();
          nBytes += hitTree->GetEvent(iTrack);

          // Get the number of hits in the TRD created by this particle
          Int_t nHit = trd->Hits()->GetEntriesFast();

          // Loop through the TRD hits  
          for (Int_t iHit = 0; iHit < nHit; iHit++) {

            if (!(hit = (AliTRDhit *) trd->Hits()->UncheckedAt(iHit))) 
              continue;

            Float_t pos[3];
                    pos[0]   = hit->fX;
                    pos[1]   = hit->fY;
                    pos[2]   = hit->fZ;
            Int_t   track    = hit->fTrack;
            Int_t   detector = hit->GetDetector();
            Int_t   plane    = geo->GetPlane(detector);
            Int_t   sector   = geo->GetSector(detector);
            Int_t   chamber  = geo->GetChamber(detector);        

            if ((sector  != isect) ||
                (plane   != iplan) ||
                (chamber != icham)) 
              continue;

            // Rotate the sectors on top of each other
            Float_t rot[3];
            geo->Rotate(detector,pos,rot);

            // Add this recPoint to the temporary array for this chamber
            AliTRDrecPoint *recPoint = new AliTRDrecPoint();
            recPoint->SetLocalRow(rot[2]);
            recPoint->SetLocalCol(rot[1]);
            recPoint->SetLocalTime(rot[0]);
            recPoint->SetEnergy(0);
            recPoint->SetDetector(detector);
            recPoint->AddDigit(track);
            chamberArray->Add(recPoint);

	  }

	}
  
        // Loop through the temporary cluster-array
        for (Int_t iClus1 = 0; iClus1 < chamberArray->GetEntries(); iClus1++) {

          AliTRDrecPoint *recPoint1 = (AliTRDrecPoint *) 
                                      chamberArray->UncheckedAt(iClus1);
          Float_t row1  = recPoint1->GetLocalRow();
          Float_t col1  = recPoint1->GetLocalCol();
          Float_t time1 = recPoint1->GetLocalTime();

          if (recPoint1->GetEnergy() < 0) continue;        // Skip marked cluster  

          const Int_t kNsave  = 5;
          Int_t idxSave[kNsave];
          Int_t iSave = 0;

          const Int_t kNsaveTrack = 3;
          Int_t tracks[kNsaveTrack];
          tracks[0] = recPoint1->GetDigit(0);

          // Check the other cluster to see, whether there are close ones
          for (Int_t iClus2 = iClus1 + 1; iClus2 < chamberArray->GetEntries(); iClus2++) {

            AliTRDrecPoint *recPoint2 = (AliTRDrecPoint *) 
                                        chamberArray->UncheckedAt(iClus2);
            Float_t row2 = recPoint2->GetLocalRow();
            Float_t col2 = recPoint2->GetLocalCol();

            if ((TMath::Abs(row1 - row2) < rowPadSize) ||
                (TMath::Abs(col1 - col2) <  fRphiDist)) {
              if (iSave == kNsave) {
                printf("AliTRDclusterizerV0::MakeCluster -- ");
                printf("Boundary error: iSave = %d, kNsave = %d.\n"
                      ,iSave,kNsave);
	      }
              else {                              
                idxSave[iSave]  = iClus2;
                iSave++;
                if (iSave < kNsaveTrack) tracks[iSave] = recPoint2->GetDigit(0);
	      }
	    }
	  }
     
          // Merge close cluster
          Float_t rowMerge = row1;
          Float_t colMerge = col1;
          if (iSave) {
            for (Int_t iMerge = 0; iMerge < iSave; iMerge++) {
              AliTRDrecPoint *recPoint2 =
                (AliTRDrecPoint *) chamberArray->UncheckedAt(idxSave[iMerge]);
              rowMerge += recPoint2->GetLocalRow();
              colMerge += recPoint2->GetLocalCol();
              recPoint2->SetEnergy(-1);     // Mark merged cluster
	    }
            rowMerge /= (iSave + 1);
            colMerge /= (iSave + 1);
          }

          Float_t smear[3];

          // The position smearing in row-direction (uniform over pad width)            
          Int_t row = (Int_t) ((rowMerge - row0) / rowPadSize);
          smear[0]  = (row + gRandom->Rndm()) * rowPadSize + row0;

          // The position smearing in rphi-direction (Gaussian)
          smear[1] = 0;
          do
            smear[1] = gRandom->Gaus(colMerge,fRphiSigma);
          while ((smear[1] < col0                        ) ||
                 (smear[1] > col0 + nColMax * colPadSize));

          // Time direction stays unchanged
          smear[2] = time1;
         
	  // Transform into local coordinates
          smear[0] = (Int_t) ((smear[0] -  row0) /  rowPadSize);
          smear[1] = (Int_t) ((smear[1] -  col0) /  colPadSize);
          smear[2] = (Int_t) ((smear[2] - time0) / timeBinSize);

          // Add the smeared cluster to the output array 
          Int_t detector  = recPoint1->GetDetector();
          Int_t digits[3] = {0};
          trd->AddRecPoint(smear,digits,detector,0.0);

	}

        // Clear the temporary cluster-array and delete the cluster
        chamberArray->Delete();

      }
    }
  }

  printf("AliTRDclusterizerV0::MakeCluster -- ");
  printf("Found %d points.\n",trd->RecPoints()->GetEntries());
  printf("AliTRDclusterizerV0::MakeCluster -- ");
  printf("Fill the cluster tree.\n");
  clusterTree->Fill();

  return kTRUE;

}
