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
Revision 1.13  1999/11/02 17:14:51  fca
Correct ansi scoping not accepted by HP compilers

Revision 1.12  1999/11/02 16:35:56  fca
New version of TRD introduced

Revision 1.11  1999/11/01 20:41:51  fca
Added protections against using the wrong version of FRAME

Revision 1.10  1999/09/29 09:24:35  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Transition Radiation Detector version 2 -- slow simulator                //
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
#include <TVector.h>
#include <TRandom.h>

#include "AliTRDv1.h"
#include "AliTRDmatrix.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliConst.h"

ClassImp(AliTRDv1)

//_____________________________________________________________________________
AliTRDv1::AliTRDv1(const char *name, const char *title) 
         :AliTRD(name, title) 
{
  //
  // Standard constructor for Transition Radiation Detector version 2
  //

  fIdSens        = 0;

  fIdChamber1    = 0;
  fIdChamber2    = 0;
  fIdChamber3    = 0;

  fSensSelect    = 0;
  fSensPlane     = 0;
  fSensChamber   = 0;
  fSensSector    = 0;

  fGasGain       = 0;
  fNoise         = 0;
  fChipGain      = 0;
  fADCoutRange   = 0;
  fADCinRange    = 0;
  fADCthreshold  = 0;

  fDiffusionT    = 0;
  fDiffusionL    = 0;

  fClusMaxThresh = 0;
  fClusSigThresh = 0;
  fClusMethod    = 0;

  fDeltaE        = NULL;

  SetBufferSize(128000);

}

//_____________________________________________________________________________
AliTRDv1::~AliTRDv1()
{

  if (fDeltaE) delete fDeltaE;

}
 
//_____________________________________________________________________________
void AliTRDv1::CreateGeometry()
{
  //
  // Create the GEANT geometry for the Transition Radiation Detector - Version 2
  // This version covers the full azimuth. 
  //

  // Check that FRAME is there otherwise we have no place where to put the TRD
  AliModule* FRAME = gAlice->GetModule("FRAME");
  if (!FRAME) return;

  // Define the chambers
  AliTRD::CreateGeometry();

}

//_____________________________________________________________________________
void AliTRDv1::CreateMaterials()
{
  //
  // Create materials for the Transition Radiation Detector version 2
  //

  AliTRD::CreateMaterials();

}

//_____________________________________________________________________________
void AliTRDv1::Diffusion(Float_t driftlength, Float_t *xyz)
{
  //
  // Applies the diffusion smearing to the position of a single electron
  //

  if ((driftlength >        0) && 
      (driftlength < kDrThick)) {
    Float_t driftSqrt = TMath::Sqrt(driftlength);
    Float_t sigmaT = driftSqrt * fDiffusionT;
    Float_t sigmaL = driftSqrt * fDiffusionL;
    xyz[0] = gRandom->Gaus(xyz[0], sigmaL);
    xyz[1] = gRandom->Gaus(xyz[1], sigmaT);
    xyz[2] = gRandom->Gaus(xyz[2], sigmaT);
  }
  else {
    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
  }

}

//_____________________________________________________________________________
void AliTRDv1::Hits2Digits()
{
  //
  // Creates TRD digits from hits. This procedure includes the following:
  //      - Diffusion
  //      - Gas gain including fluctuations
  //      - Pad-response (simple Gaussian approximation)
  //      - Electronics noise
  //      - Electronics gain
  //      - Digitization
  //      - ADC threshold
  // The corresponding parameter can be adjusted via the various Set-functions.
  // If these parameters are not explicitly set, default values are used (see
  // Init-function).
  // To produce digits from a root-file with TRD-hits use the
  // slowDigitsCreate.C macro.
  //

  printf("AliTRDv1::Hits2Digits -- Start creating digits\n");

  ///////////////////////////////////////////////////////////////
  // Parameter 
  ///////////////////////////////////////////////////////////////

  // Converts number of electrons to fC
  const Float_t el2fC = 1.602E-19 * 1.0E15; 

  ///////////////////////////////////////////////////////////////

  Int_t nBytes = 0;

  Int_t iRow;

  AliTRDhit *TRDhit;

  // Get the pointer to the hit tree
  TTree *HitTree    = gAlice->TreeH();
  // Get the pointer to the digits tree
  TTree *DigitsTree = gAlice->TreeD();

  // Get the number of entries in the hit tree
  // (Number of primary particles creating a hit somewhere)
  Int_t nTrack = (Int_t) HitTree->GetEntries();

  Int_t chamBeg = 0;
  Int_t chamEnd = kNcham;
  if (fSensChamber) chamEnd = chamBeg = fSensChamber;
  Int_t planBeg = 0;
  Int_t planEnd = kNplan;
  if (fSensPlane)   planEnd = planBeg = fSensPlane;
  Int_t sectBeg = 0;
  Int_t sectEnd = kNsect;
  if (fSensSector)  sectEnd = sectBeg = fSensSector;

  // Loop through all the chambers
  for (Int_t icham = chamBeg; icham < chamEnd; icham++) {
    for (Int_t iplan = planBeg; iplan < planEnd; iplan++) {
      for (Int_t isect = sectBeg; isect < sectEnd; isect++) {

        Int_t nDigits = 0;

        printf("AliTRDv1::Hits2Digits -- Digitizing chamber %d, plane %d, sector %d\n"
              ,icham+1,iplan+1,isect+1);

        // Create a detector matrix to keep the signal and track numbers
        AliTRDmatrix *matrix = new AliTRDmatrix(fRowMax[iplan][icham][isect]
                                               ,fColMax[iplan]
                                               ,fTimeMax
                                               ,isect+1,icham+1,iplan+1);

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
            Float_t q       = TRDhit->fQ;
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

            // The hit position in pad coordinates (center pad)
            // The pad row (z-direction)
            Int_t rowH  = (Int_t) ((zRot -  fRow0[iplan][icham][isect]) / fRowPadSize);
            // The pad column (rphi-direction)  
            Int_t colH  = (Int_t) ((yRot -  fCol0[iplan]              ) / fColPadSize);
            // The time bucket
            Int_t timeH = (Int_t) ((xRot - fTime0[iplan]              ) / fTimeBinSize);

            // Array to sum up the signal in a box surrounding the
            // hit postition
            const Int_t timeBox = 5;
            const Int_t  colBox = 7;
            const Int_t  rowBox = 5;
            Float_t signalSum[rowBox][colBox][timeBox];
            for (iRow  = 0;  iRow <  rowBox; iRow++ ) {
              for (Int_t iCol  = 0;  iCol <  colBox; iCol++ ) {
                for (Int_t iTime = 0; iTime < timeBox; iTime++) {
                  signalSum[iRow][iCol][iTime] = 0;
		}
	      }
	    }

            // Loop over all electrons of this hit
            Int_t nEl = (Int_t) q;
            for (Int_t iEl = 0; iEl < nEl; iEl++) {

              // Apply the diffusion smearing
              Float_t driftlength = xRot - fTime0[iplan];
              Float_t xyz[3];
              xyz[0] = xRot;
              xyz[1] = yRot;
              xyz[2] = zRot;
              Diffusion(driftlength,xyz);

              // At this point absorption effects that depend on the 
	      // driftlength could be taken into account.              

              // The electron position and the distance to the hit position
	      // in pad units
              // The pad row (z-direction)
              Int_t  rowE = (Int_t) ((xyz[2] -  fRow0[iplan][icham][isect]) / fRowPadSize);
              Int_t  rowD =  rowH -  rowE;
              // The pad column (rphi-direction)
              Int_t  colE = (Int_t) ((xyz[1] -  fCol0[iplan]              ) / fColPadSize);
              Int_t  colD =  colH -  colE;
              // The time bucket
              Int_t timeE = (Int_t) ((xyz[0] - fTime0[iplan]              ) / fTimeBinSize);
              Int_t timeD = timeH - timeE;

              // Apply the gas gain including fluctuations
              Int_t signal = (Int_t) (-fGasGain * TMath::Log(gRandom->Rndm()));

	      // The distance of the electron to the center of the pad 
	      // in units of pad width
              Float_t dist = (xyz[1] - fCol0[iplan] - (colE + 0.5) * fColPadSize) 
                           / fColPadSize;

              // Sum up the signal in the different pixels
              // and apply the pad response
              Int_t  rowIdx =  rowD + (Int_t) ( rowBox / 2);
              Int_t  colIdx =  colD + (Int_t) ( colBox / 2);
              Int_t timeIdx = timeD + (Int_t) (timeBox / 2);
              signalSum[rowIdx][colIdx-1][timeIdx] += PadResponse(dist-1.) * signal;
              signalSum[rowIdx][colIdx  ][timeIdx] += PadResponse(dist   ) * signal;
              signalSum[rowIdx][colIdx+1][timeIdx] += PadResponse(dist+1.) * signal;

            }

            // Add the padcluster to the detector matrix
            for (iRow  = 0;  iRow <  rowBox; iRow++ ) {
              for (Int_t iCol  = 0;  iCol <  colBox; iCol++ ) {
                for (Int_t iTime = 0; iTime < timeBox; iTime++) {

                  Int_t  rowB =  rowH + iRow  - (Int_t) ( rowBox / 2); 
                  Int_t  colB =  colH + iCol  - (Int_t) ( colBox / 2);
                  Int_t timeB = timeH + iTime - (Int_t) (timeBox / 2);

                  Float_t signalB = signalSum[iRow][iCol][iTime];
                  if (signalB > 0.0) {
                    matrix->AddSignal(rowB,colB,timeB,signalB);
                    if (!(matrix->AddTrack(rowB,colB,timeB,track))) 
                      printf(" More than three tracks in a pixel!\n");
		  }

		}
	      }
	    }

          }

	}

        // Create the hits for this chamber
        for (Int_t iRow  = 0; iRow  <  fRowMax[iplan][icham][isect]; iRow++ ) {
          for (Int_t iCol  = 0; iCol  <  fColMax[iplan]              ; iCol++ ) {
            for (Int_t iTime = 0; iTime < fTimeMax                     ; iTime++) {         

              Float_t signalAmp = matrix->GetSignal(iRow,iCol,iTime);

              // Add the noise
              signalAmp  = TMath::Max(gRandom->Gaus(signalAmp,fNoise),(Float_t) 0.0);
	      // Convert to fC
              signalAmp *= el2fC;
              // Convert to mV
              signalAmp *= fChipGain;
	      // Convert to ADC counts
              Int_t adc  = (Int_t) (signalAmp * (fADCoutRange / fADCinRange));

	      // Apply threshold on ADC value
              if (adc > fADCthreshold) {

                Int_t trackSave[3];
                for (Int_t ii = 0; ii < 3; ii++) {
                  trackSave[ii] = matrix->GetTrack(iRow,iCol,iTime,ii);
	        }

                Int_t digits[7];
                digits[0] = matrix->GetSector();
                digits[1] = matrix->GetChamber();
                digits[2] = matrix->GetPlane();
                digits[3] = iRow;
                digits[4] = iCol;
                digits[5] = iTime;
                digits[6] = adc;

		// Add this digit to the TClonesArray
                AddDigit(trackSave,digits);
                nDigits++;

	      }

	    }
	  }
	}

        printf("AliTRDv1::Hits2Digits -- Number of digits found: %d\n",nDigits);

	// Clean up
        delete matrix;

      }
    }
  }

  // Fill the digits-tree
  printf("AliTRDv1::Hits2Digits -- Fill the digits tree\n");
  DigitsTree->Fill();

}

//_____________________________________________________________________________
void AliTRDv1::Digits2Clusters()
{

  //
  // Method to convert AliTRDdigits created by AliTRDv1::Hits2Digits()
  // into AliTRDclusters
  // To produce cluster from a root-file with TRD-digits use the
  // slowClusterCreate.C macro.
  //
  
  Int_t row;

  printf("AliTRDv1::Digits2Clusters -- Start creating clusters\n");

  AliTRDdigit  *TRDdigit;
  TClonesArray *TRDDigits;

  // Parameters
  Float_t maxThresh        = fClusMaxThresh;   // threshold value for maximum
  Float_t signalThresh     = fClusSigThresh;   // threshold value for digit signal
  Int_t   clusteringMethod = fClusMethod;      // clustering method option (for testing)

  const Float_t epsilon    = 0.01;             // iteration limit for unfolding procedure

  // Get the pointer to the digits tree
  TTree *DigitTree   = gAlice->TreeD();
  // Get the pointer to the cluster tree
  TTree *ClusterTree = gAlice->TreeD();

  // Get the pointer to the digits container
  TRDDigits = Digits();

  Int_t chamBeg = 0;
  Int_t chamEnd = kNcham;
  if (fSensChamber) chamEnd = chamBeg = fSensChamber;
  Int_t planBeg = 0;
  Int_t planEnd = kNplan;
  if (fSensPlane)   planEnd = planBeg = fSensPlane;
  Int_t sectBeg = 0;
  Int_t sectEnd = kNsect;
  if (fSensSector)  sectEnd = sectBeg = fSensSector;

  // Import the digit tree
  gAlice->ResetDigits();
  Int_t nbytes;
  nbytes += DigitTree->GetEvent(1);

  // Get the number of digits in the detector
  Int_t nTRDDigits = TRDDigits->GetEntriesFast();

  // *** Start clustering *** in every chamber
  for (Int_t icham = chamBeg; icham < chamEnd; icham++) {
    for (Int_t iplan = planBeg; iplan < planEnd; iplan++) {
      for (Int_t isect = sectBeg; isect < sectEnd; isect++) {

        Int_t nClusters = 0;
        printf("AliTRDv1::Digits2Clusters -- Finding clusters in chamber %d, plane %d, sector %d\n"
               ,icham+1,iplan+1,isect+1);

        // Create a detector matrix to keep maxima
        AliTRDmatrix *digitMatrix = new AliTRDmatrix(fRowMax[iplan][icham][isect]
                                                     ,fColMax[iplan]
                                                     ,fTimeMax,isect+1
                                                     ,icham+1,iplan+1);
        // Create a matrix to contain maximum flags
        AliTRDmatrix *maximaMatrix = new AliTRDmatrix(fRowMax[iplan][icham][isect]
                                                      ,fColMax[iplan]
                                                      ,fTimeMax
                                                      ,isect+1,icham+1,iplan+1);

        // Loop through all TRD digits
        for (Int_t iTRDDigits = 0; iTRDDigits < nTRDDigits; iTRDDigits++) {

          // Get the information for this digit
          TRDdigit = (AliTRDdigit*) TRDDigits->UncheckedAt(iTRDDigits);
          Int_t   signal  = TRDdigit->fSignal;
          Int_t   sector  = TRDdigit->fSector;
          Int_t   chamber = TRDdigit->fChamber;
          Int_t   plane   = TRDdigit->fPlane;
          Int_t   row     = TRDdigit->fRow;
          Int_t   col     = TRDdigit->fCol;
          Int_t   time    = TRDdigit->fTime;

          Int_t   track[3];
          for (Int_t iTrack = 0; iTrack < 3; iTrack++) {
            track[iTrack]  = TRDdigit->AliDigit::fTracks[iTrack];
          }

          if ((sector  != isect+1) ||
              (plane   != iplan+1) ||
              (chamber != icham+1))
            continue;

          // Fill the detector matrix
          if (signal > signalThresh) {
            digitMatrix->SetSignal(row,col,time,signal);
            for (Int_t iTrack = 0; iTrack < 3; iTrack++) {
              if (track[iTrack] > 0) {
                digitMatrix->AddTrack(row,col,time,track[iTrack]);
              }
            }
          }

        }

        // Loop chamber and find maxima in digitMatrix
        for (row = 0;  row < fRowMax[iplan][icham][isect];  row++) {
          for (Int_t  col = 1;  col < fColMax[iplan]              ;  col++) {
            for (Int_t time = 0; time < fTimeMax                    ; time++) {

              if (digitMatrix->GetSignal(row,col,time) 
                  < digitMatrix->GetSignal(row,col - 1,time)) {
                // really maximum?
                if (col > 1) {
                  if (digitMatrix->GetSignal(row,col - 2,time)
                      < digitMatrix->GetSignal(row,col - 1,time)) {
                    // yes, so set maximum flag
                    maximaMatrix->SetSignal(row,col - 1,time,1);
                  }
                  else maximaMatrix->SetSignal(row,col - 1,time,0);
                }
              }

            }   // time
          }     // col
        }       // row

        // now check maxima and calculate cluster position
        for (row = 0;  row < fRowMax[iplan][icham][isect];  row++) {
          for (Int_t  col = 1;  col < fColMax[iplan]              ;  col++) {
            for (Int_t time = 0; time < fTimeMax                    ; time++) {

              if ((maximaMatrix->GetSignal(row,col,time) > 0)
                  && (digitMatrix->GetSignal(row,col,time) > maxThresh)) {

                Int_t   clusters[5]        = {0};   // cluster-object data

                Float_t ratio              =  0;    // ratio resulting from unfolding
                Float_t padSignal[5]       = {0};   // signals on max and neighbouring pads
                Float_t clusterSignal[3]   = {0};   // signals from cluster
                Float_t clusterPos[3]      = {0};   // cluster in ALICE refFrame coords
                Float_t clusterPads[6]     = {0};   // cluster pad info

                // setting values
                clusters[0] = isect+1;    // = isect ????
                clusters[1] = icham+1;    // = ichamber ????
                clusters[2] = iplan+1;    // = iplane ????
                clusters[3] = time;

                clusterPads[0] = icham+1;
                clusterPads[1] = isect+1;
                clusterPads[2] = iplan+1;

                for (Int_t iPad = 0; iPad < 3; iPad++) {
                  clusterSignal[iPad] = digitMatrix->GetSignal(row,col-1+iPad,time);
                }

                // neighbouring maximum on right side?
                if (col < fColMax[iplan] - 2) {
                  if (maximaMatrix->GetSignal(row,col + 2,time) > 0) {
                    for (Int_t iPad = 0; iPad < 5; iPad++) {
                      padSignal[iPad] = digitMatrix->GetSignal(row,col-1+iPad,time);
                    }

                    // unfold:
                    ratio = Unfold(epsilon, padSignal);

                    // set signal on overlapping pad to ratio
                    clusterSignal[2] *= ratio;
                  }
                }

                switch (clusteringMethod) {
                case 1:
                  // method 1: simply center of mass
                  clusterPads[3] = row + 0.5;
                  clusterPads[4] = col - 0.5 + (clusterSignal[2] - clusterSignal[0]) /
                                   (clusterSignal[1] + clusterSignal[2] + clusterSignal[3]);
                  clusterPads[5] = time + 0.5;

                  nClusters++;
                  break;
                case 2:
                  // method 2: integral gauss fit on 3 pads
                  TH1F *hPadCharges = new TH1F("hPadCharges", "Charges on center 3 pads"
	                                                    , 5, -1.5, 3.5);
                  for (Int_t iCol = -1; iCol <= 3; iCol++) {
                    if (clusterSignal[iCol] < 1) clusterSignal[iCol] = 1;
                    hPadCharges->Fill(iCol, clusterSignal[iCol]);
                  }
                  hPadCharges->Fit("gaus", "IQ", "SAME", -0.5, 2.5);
                  TF1     *fPadChargeFit = hPadCharges->GetFunction("gaus");
                  Double_t  colMean = fPadChargeFit->GetParameter(1);

                  clusterPads[3] = row + 0.5;
                  clusterPads[4] = col - 1.5 + colMean;
                  clusterPads[5] = time + 0.5;

                  delete hPadCharges;

                  nClusters++;
                  break;
                }

                Float_t clusterCharge =   clusterSignal[0]
                                        + clusterSignal[1]
                                        + clusterSignal[2];
                clusters[4] = (Int_t)clusterCharge;

                Int_t trackSave[3];
                for (Int_t iTrack = 0; iTrack < 3; iTrack++) {
                  trackSave[iTrack] = digitMatrix->GetTrack(row,col,time,iTrack);
                }

                // Calculate cluster position in ALICE refFrame coords
                // and set array clusterPos to calculated values
                Pads2XYZ(clusterPads, clusterPos);

                // Add cluster to reconstruction tree
                AddCluster(trackSave,clusters,clusterPos);

              }

            }  // time
          }    // col
        }      // row

        printf("AliTRDv1::Digits2Clusters -- Number of clusters found: %d\n",nClusters);

        delete digitMatrix;
        delete maximaMatrix;

      }          // isect
    }            // iplan
  }              // icham

  // Fill the cluster-tree
  printf("AliTRDv1::Digits2Clusters -- Total number of clusters found: %d\n"
        ,fClusters->GetEntries());
  printf("AliTRDv1::Digits2Clusters -- Fill the cluster tree\n");
  ClusterTree->Fill();

}

//_____________________________________________________________________________
void AliTRDv1::Init() 
{
  //
  // Initialise Transition Radiation Detector after geometry has been built.
  // Includes the default settings of all parameter for the digitization.
  //

  AliTRD::Init();

  printf("          Slow simulator\n");
  if (fSensPlane)
    printf("          Only plane %d is sensitive\n",fSensPlane);
  if (fSensChamber)   
    printf("          Only chamber %d is sensitive\n",fSensChamber);
  if (fSensSector)
    printf("          Only sector %d is sensitive\n",fSensSector);

  // First ionization potential (eV) for the gas mixture (90% Xe + 10% CO2)
  const Float_t kPoti = 12.1;
  // Maximum energy (50 keV);
  const Float_t kEend = 50000.0;
  // Ermilova distribution for the delta-ray spectrum
  Float_t Poti = TMath::Log(kPoti);
  Float_t Eend = TMath::Log(kEend);
  fDeltaE  = new TF1("deltae",Ermilova,Poti,Eend,0);

  // Identifier of the sensitive volume (drift region)
  fIdSens     = gMC->VolId("UL05");

  // Identifier of the TRD-driftchambers
  fIdChamber1 = gMC->VolId("UCIO");
  fIdChamber2 = gMC->VolId("UCIM");
  fIdChamber3 = gMC->VolId("UCII");

  // The default parameter for the digitization
  if (!(fGasGain))       fGasGain       = 2.0E3;
  if (!(fNoise))         fNoise         = 3000.;
  if (!(fChipGain))      fChipGain      = 10.;
  if (!(fADCoutRange))   fADCoutRange   = 255.;
  if (!(fADCinRange))    fADCinRange    = 2000.;
  if (!(fADCthreshold))  fADCthreshold  = 1;

  // Transverse and longitudinal diffusion coefficients (Xe/Isobutane)
  if (!(fDiffusionT))    fDiffusionT    = 0.060;
  if (!(fDiffusionL))    fDiffusionL    = 0.017;

  // The default parameter for the clustering
  if (!(fClusMaxThresh)) fClusMaxThresh = 5.0;
  if (!(fClusSigThresh)) fClusSigThresh = 2.0;
  if (!(fClusMethod))    fClusMethod    = 1;

  for (Int_t i = 0; i < 80; i++) printf("*");
  printf("\n");

}

//_____________________________________________________________________________
Float_t AliTRDv1::PadResponse(Float_t x)
{
  //
  // The pad response for the chevron pads. 
  // We use a simple Gaussian approximation which should be good
  // enough for our purpose.
  //

  // The parameters for the response function
  const Float_t aa  =  0.8872;
  const Float_t bb  = -0.00573;
  const Float_t cc  =  0.454;
  const Float_t cc2 =  cc*cc;

  Float_t pr = aa * (bb + TMath::Exp(-x*x / (2. * cc2)));

  return (pr);

}

//_____________________________________________________________________________
void AliTRDv1::SetSensPlane(Int_t iplane)
{
  //
  // Defines the hit-sensitive plane (1-6)
  //

  if ((iplane < 0) || (iplane > 6)) {
    printf("Wrong input value: %d\n",iplane);
    printf("Use standard setting\n");
    fSensPlane  = 0;
    fSensSelect = 0;
    return;
  }

  fSensSelect = 1;
  fSensPlane  = iplane;

}

//_____________________________________________________________________________
void AliTRDv1::SetSensChamber(Int_t ichamber)
{
  //
  // Defines the hit-sensitive chamber (1-5)
  //

  if ((ichamber < 0) || (ichamber > 5)) {
    printf("Wrong input value: %d\n",ichamber);
    printf("Use standard setting\n");
    fSensChamber = 0;
    fSensSelect  = 0;
    return;
  }

  fSensSelect  = 1;
  fSensChamber = ichamber;

}

//_____________________________________________________________________________
void AliTRDv1::SetSensSector(Int_t isector)
{
  //
  // Defines the hit-sensitive sector (1-18)
  //

  if ((isector < 0) || (isector > 18)) {
    printf("Wrong input value: %d\n",isector);
    printf("Use standard setting\n");
    fSensSector = 0;
    fSensSelect = 0;
    return;
  }

  fSensSelect = 1;
  fSensSector = isector;

}

//_____________________________________________________________________________
void AliTRDv1::StepManager()
{
  //
  // Called at every step in the Transition Radiation Detector version 2.
  // Slow simulator. Every charged track produces electron cluster as hits 
  // along its path across the drift volume. The step size is set acording
  // to Bethe-Bloch. The energy distribution of the delta electrons follows
  // a spectrum taken from Ermilova et al.
  //

  Int_t    iIdSens, icSens;
  Int_t    iIdSpace, icSpace;
  Int_t    iIdChamber, icChamber;
  Int_t    vol[3]; 
  Int_t    iPid;

  Float_t  hits[4];
  Float_t  random[1];
  Float_t  charge;
  Float_t  aMass;

  Double_t pTot;
  Double_t qTot;
  Double_t eDelta;
  Double_t betaGamma, pp;

  TLorentzVector pos, mom;
  TClonesArray  &lhits = *fHits;

  const Double_t kBig = 1.0E+12;

  // Ionization energy
  const Float_t kWion    = 22.04;
  // Maximum energy for e+ e- g for the step-size calculation
  const Float_t kPTotMax = 0.002;
  // Plateau value of the energy-loss for electron in xenon
  // taken from: Allison + Comb, Ann. Rev. Nucl. Sci. (1980), 30, 253
  //const Double_t kPlateau = 1.70;
  // the averaged value (26/3/99)
  const Float_t kPlateau = 1.55;
  // dN1/dx|min for the gas mixture (90% Xe + 10% CO2)
  const Float_t kPrim    = 48.0;
  // First ionization potential (eV) for the gas mixture (90% Xe + 10% CO2)
  const Float_t kPoti    = 12.1;

  // Set the maximum step size to a very large number for all 
  // neutral particles and those outside the driftvolume
  gMC->SetMaxStep(kBig); 

  // Use only charged tracks 
  if (( gMC->TrackCharge()       ) &&
      (!gMC->IsTrackStop()       ) && 
      (!gMC->IsTrackDisappeared())) {

    // Inside a sensitive volume?
    iIdSens = gMC->CurrentVolID(icSens);
    if (iIdSens == fIdSens) { 

      iIdSpace   = gMC->CurrentVolOffID(4,icSpace  );
      iIdChamber = gMC->CurrentVolOffID(1,icChamber);

      // Calculate the energy of the delta-electrons
      eDelta = TMath::Exp(fDeltaE->GetRandom()) - kPoti;
      eDelta = TMath::Max(eDelta,0.0);

      // The number of secondary electrons created
      qTot = (Double_t) ((Int_t) (eDelta / kWion) + 1);

      // The hit coordinates and charge
      gMC->TrackPosition(pos);
      hits[0] = pos[0];
      hits[1] = pos[1];
      hits[2] = pos[2];
      hits[3] = qTot;

      // The sector number
      Float_t phi = pos[1] != 0 ? TMath::Atan2(pos[0],pos[1]) : (pos[0] > 0 ? 180. : 0.);
      vol[0] = ((Int_t) (phi / 20)) + 1;

      // The chamber number 
      //   1: outer left
      //   2: middle left
      //   3: inner
      //   4: middle right
      //   5: outer right
      if      (iIdChamber == fIdChamber1)
        vol[1] = (hits[2] < 0 ? 1 : 5);
      else if (iIdChamber == fIdChamber2)       
        vol[1] = (hits[2] < 0 ? 2 : 4);
      else if (iIdChamber == fIdChamber3)       
        vol[1] = 3;

      // The plane number
      vol[2] = icChamber - TMath::Nint((Float_t) (icChamber / 7)) * 6;

      // Check on selected volumes
      Int_t addthishit = 1;
      if (fSensSelect) {
        if ((fSensPlane)   && (vol[2] != fSensPlane  )) addthishit = 0;
        if ((fSensChamber) && (vol[1] != fSensChamber)) addthishit = 0;
        if ((fSensSector)  && (vol[0] != fSensSector )) addthishit = 0;
      }

      // Add this hit
      if (addthishit) {

        new(lhits[fNhits++]) AliTRDhit(fIshunt,gAlice->CurrentTrack(),vol,hits);

        // The energy loss according to Bethe Bloch
        gMC->TrackMomentum(mom);
        pTot = mom.Rho();
        iPid = gMC->TrackPid();
        if ( (iPid >  3) ||
	    ((iPid <= 3) && (pTot < kPTotMax))) {
          aMass     = gMC->TrackMass();
          betaGamma = pTot / aMass;
          pp        = kPrim * BetheBloch(betaGamma);
	  // Take charge > 1 into account
          charge = gMC->TrackCharge();
          if (TMath::Abs(charge) > 1) pp = pp * charge*charge;
        }
        // Electrons above 20 Mev/c are at the plateau
        else {
          pp = kPrim * kPlateau;
        }
      
        // Calculate the maximum step size for the next tracking step
        if (pp > 0) {
          do 
            gMC->Rndm(random,1);
          while ((random[0] == 1.) || (random[0] == 0.));
          gMC->SetMaxStep( - TMath::Log(random[0]) / pp);
	}

      }
      else {
        // set step size to maximal value
        gMC->SetMaxStep(kBig); 
      }

    }

  }

}

//_____________________________________________________________________________
Double_t AliTRDv1::BetheBloch(Double_t bg) 
{
  //
  // Parametrization of the Bethe-Bloch-curve
  // The parametrization is the same as for the TPC and is taken from Lehrhaus.
  //

  // This parameters have been adjusted to averaged values from GEANT
  const Double_t kP1 = 7.17960e-02;
  const Double_t kP2 = 8.54196;
  const Double_t kP3 = 1.38065e-06;
  const Double_t kP4 = 5.30972;
  const Double_t kP5 = 2.83798;

  // This parameters have been adjusted to Xe-data found in:
  // Allison & Cobb, Ann. Rev. Nucl. Sci. (1980), 30, 253
  //const Double_t kP1 = 0.76176E-1;
  //const Double_t kP2 = 10.632;
  //const Double_t kP3 = 3.17983E-6;
  //const Double_t kP4 = 1.8631;
  //const Double_t kP5 = 1.9479;

  if (bg > 0) {
    Double_t yy = bg / TMath::Sqrt(1. + bg*bg);
    Double_t aa = TMath::Power(yy,kP4);
    Double_t bb = TMath::Power((1./bg),kP5);
             bb = TMath::Log(kP3 + bb);
    return ((kP2 - aa - bb)*kP1 / aa);
  }
  else
    return 0;

}

//_____________________________________________________________________________
Double_t Ermilova(Double_t *x, Double_t *)
{
  //
  // Calculates the delta-ray energy distribution according to Ermilova.
  // Logarithmic scale !
  //

  Double_t energy;
  Double_t dpos;
  Double_t dnde;

  Int_t    pos1, pos2;

  const Int_t nV = 31;

  Float_t vxe[nV] = { 2.3026, 2.9957, 3.4012, 3.6889, 3.9120  
                    , 4.0943, 4.2485, 4.3820, 4.4998, 4.6052
                    , 4.7005, 5.0752, 5.2983, 5.7038, 5.9915
                    , 6.2146, 6.5221, 6.9078, 7.3132, 7.6009
                    , 8.0064, 8.5172, 8.6995, 8.9872, 9.2103
                    , 9.4727, 9.9035,10.3735,10.5966,10.8198
                    ,11.5129 };

  Float_t vye[nV] = { 80.0  , 31.0  , 23.3  , 21.1  , 21.0
                    , 20.9  , 20.8  , 20.0  , 16.0  , 11.0
                    ,  8.0  ,  6.0  ,  5.2  ,  4.6  ,  4.0
                    ,  3.5  ,  3.0  ,  1.4  ,  0.67 ,  0.44
                    ,  0.3  ,  0.18 ,  0.12 ,  0.08 ,  0.056
                    ,  0.04 ,  0.023,  0.015,  0.011,  0.01
		    ,  0.004 };

  energy = x[0];

  // Find the position 
  pos1 = pos2 = 0;
  dpos = 0;
  do {
    dpos = energy - vxe[pos2++];
  } 
  while (dpos > 0);
  pos2--; 
  if (pos2 > nV) pos2 = nV;
  pos1 = pos2 - 1;

  // Differentiate between the sampling points
  dnde = (vye[pos1] - vye[pos2]) / (vxe[pos2] - vxe[pos1]);

  return dnde;

}

//_____________________________________________________________________________
void AliTRDv1::Pads2XYZ(Float_t *pads, Float_t *pos)
{
  // Method to convert pad coordinates (row,col,time)
  // into ALICE reference frame coordinates (x,y,z)

  Int_t        chamber   = (Int_t) pads[0];     // chamber info (1-5)
  Int_t        sector    = (Int_t) pads[1];     // sector info  (1-18)
  Int_t        plane     = (Int_t) pads[2];     // plane info   (1-6)

  Int_t        icham     = chamber - 1;         // chamber info (0-4)
  Int_t        isect     = sector  - 1;         // sector info  (0-17)
  Int_t        iplan     = plane   - 1;         // plane info   (0-5)

  Float_t      padRow    = pads[3];             // Pad Row position
  Float_t      padCol    = pads[4];             // Pad Column position
  Float_t      timeSlice = pads[5];             // Time "position"

  // calculate (x,y) position in rotated chamber
  Float_t yRot = fCol0[iplan]               + padCol    * fColPadSize;
  Float_t xRot = fTime0[iplan]              + timeSlice * fTimeBinSize;
  // calculate z-position:
  Float_t z    = fRow0[iplan][icham][isect] + padRow    * fRowPadSize;

  /**
    rotate chamber back to original position
    1. mirror at y-axis, 2. rotate back to position (-phi)
    / cos(phi) -sin(phi) \   / -1 0 \     / -cos(phi) -sin(phi) \
    \ sin(phi)  cos(phi) / * \  0 1 /  =  \ -sin(phi)  cos(phi) /
  **/
  //Float_t phi = 2*kPI / kNsect * ((Float_t) sector - 0.5);
  //Float_t x = -xRot * TMath::Cos(phi) - yRot * TMath::Sin(phi);
  //Float_t y = -xRot * TMath::Sin(phi) + yRot * TMath::Cos(phi);
  Float_t phi = 2*kPI / kNsect * ((Float_t) sector - 0.5);
  Float_t x   = -xRot * TMath::Cos(phi) + yRot * TMath::Sin(phi);
  Float_t y   =  xRot * TMath::Sin(phi) + yRot * TMath::Cos(phi);

  // Setting values
  pos[0] = x;
  pos[1] = y;
  pos[2] = z;

}

//_____________________________________________________________________________
Float_t AliTRDv1::Unfold(Float_t eps, Float_t* padSignal)
{
  // Method to unfold neighbouring maxima.
  // The charge ratio on the overlapping pad is calculated
  // until there is no more change within the range given by eps.
  // The resulting ratio is then returned to the calling method.

  Int_t   itStep            = 0;      // count iteration steps

  Float_t ratio             = 0.5;    // start value for ratio
  Float_t prevRatio         = 0;      // store previous ratio

  Float_t newLeftSignal[3]  = {0};    // array to store left cluster signal
  Float_t newRightSignal[3] = {0};    // array to store right cluster signal

  // start iteration:
  while ((TMath::Abs(prevRatio - ratio) > eps) && (itStep < 10)) {

    itStep++;
    prevRatio = ratio;

    // cluster position according to charge ratio
    Float_t maxLeft  = (ratio*padSignal[2] - padSignal[0]) /
                       (padSignal[0] + padSignal[1] + ratio*padSignal[2]);
    Float_t maxRight = (padSignal[4] - (1-ratio)*padSignal[2]) /
                       ((1-ratio)*padSignal[2] + padSignal[3] + padSignal[4]);

    // set cluster charge ratio
    Float_t ampLeft  = padSignal[1];
    Float_t ampRight = padSignal[3];

    // apply pad response to parameters
    newLeftSignal[0] = ampLeft*PadResponse(-1 - maxLeft);
    newLeftSignal[1] = ampLeft*PadResponse( 0 - maxLeft);
    newLeftSignal[2] = ampLeft*PadResponse( 1 - maxLeft);

    newRightSignal[0] = ampRight*PadResponse(-1 - maxRight);
    newRightSignal[1] = ampRight*PadResponse( 0 - maxRight);
    newRightSignal[2] = ampRight*PadResponse( 1 - maxRight);

    // calculate new overlapping ratio
    ratio = newLeftSignal[2]/(newLeftSignal[2] + newRightSignal[0]);

  }

  return ratio;

}

