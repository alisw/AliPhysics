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
Revision 1.11  2001/05/21 16:45:47  hristov
Last minute changes (C.Blume)

Revision 1.10  2001/05/07 08:06:44  cblume
Speedup of the code. Create only AliTRDcluster

Revision 1.9  2000/11/01 14:53:20  cblume
Merge with TRD-develop

Revision 1.1.4.5  2000/10/15 23:40:01  cblume
Remove AliTRDconst

Revision 1.1.4.4  2000/10/06 16:49:46  cblume
Made Getters const

Revision 1.1.4.3  2000/10/04 16:34:58  cblume
Replace include files by forward declarations

Revision 1.1.4.2  2000/09/22 14:49:49  cblume
Adapted to tracking code

Revision 1.8  2000/10/02 21:28:19  fca
Removal of useless dependecies via forward declarations

Revision 1.7  2000/06/27 13:08:50  cblume
Changed to Copy(TObject &A) to appease the HP-compiler

Revision 1.6  2000/06/09 11:10:07  cblume
Compiler warnings and coding conventions, next round

Revision 1.5  2000/06/08 18:32:58  cblume
Make code compliant to coding conventions

Revision 1.4  2000/06/07 16:27:01  cblume
Try to remove compiler warnings on Sun and HP

Revision 1.3  2000/05/08 16:17:27  cblume
Merge TRD-develop

Revision 1.1.4.1  2000/05/08 15:09:01  cblume
Introduce AliTRDdigitsManager

Revision 1.1  2000/02/28 18:58:54  cblume
Add new TRD classes

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TRD cluster finder for the slow simulator. 
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TF1.h>
#include <TTree.h>
#include <TH1.h>

#include "AliRun.h"

#include "AliTRD.h"
#include "AliTRDclusterizerV1.h"
#include "AliTRDmatrix.h"
#include "AliTRDgeometry.h"
#include "AliTRDdigitizer.h"
#include "AliTRDdataArrayF.h"
#include "AliTRDdataArrayI.h"
#include "AliTRDdigitsManager.h"

ClassImp(AliTRDclusterizerV1)

//_____________________________________________________________________________
AliTRDclusterizerV1::AliTRDclusterizerV1():AliTRDclusterizer()
{
  //
  // AliTRDclusterizerV1 default constructor
  //

  fDigitsManager = NULL;

  fClusMaxThresh = 0;
  fClusSigThresh = 0;

  fUseLUT        = kFALSE;

}

//_____________________________________________________________________________
AliTRDclusterizerV1::AliTRDclusterizerV1(const Text_t* name, const Text_t* title)
                    :AliTRDclusterizer(name,title)
{
  //
  // AliTRDclusterizerV1 default constructor
  //

  fDigitsManager = new AliTRDdigitsManager();

  Init();

}

//_____________________________________________________________________________
AliTRDclusterizerV1::AliTRDclusterizerV1(const AliTRDclusterizerV1 &c)
{
  //
  // AliTRDclusterizerV1 copy constructor
  //

  ((AliTRDclusterizerV1 &) c).Copy(*this);

}

//_____________________________________________________________________________
AliTRDclusterizerV1::~AliTRDclusterizerV1()
{
  //
  // AliTRDclusterizerV1 destructor
  //

  if (fDigitsManager) {
    delete fDigitsManager;
  }

}

//_____________________________________________________________________________
AliTRDclusterizerV1 &AliTRDclusterizerV1::operator=(const AliTRDclusterizerV1 &c)
{
  //
  // Assignment operator
  //

  if (this != &c) ((AliTRDclusterizerV1 &) c).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDclusterizerV1::Copy(TObject &c)
{
  //
  // Copy function
  //

  ((AliTRDclusterizerV1 &) c).fUseLUT        = fUseLUT;
  ((AliTRDclusterizerV1 &) c).fClusMaxThresh = fClusMaxThresh;
  ((AliTRDclusterizerV1 &) c).fClusSigThresh = fClusSigThresh;
  ((AliTRDclusterizerV1 &) c).fDigitsManager = NULL;
  for (Int_t ilut = 0; ilut < kNlut; ilut++) {
    ((AliTRDclusterizerV1 &) c).fLUT[ilut] = fLUT[ilut];
  }

  AliTRDclusterizer::Copy(c);

}

//_____________________________________________________________________________
void AliTRDclusterizerV1::Init()
{
  //
  // Initializes the cluster finder
  //

  // The default parameter for the clustering
  fClusMaxThresh = 3;
  fClusSigThresh = 1;

  // Use the lookup table for the position determination
  fUseLUT        = kTRUE;

  // The lookup table from Bogdan
  Float_t lut[128] = {  
    0.0068,  0.0198,  0.0318,  0.0432,  0.0538,  0.0642,  0.0742,  0.0838,
    0.0932,  0.1023,  0.1107,  0.1187,  0.1268,  0.1347,  0.1423,  0.1493,  
    0.1562,  0.1632,  0.1698,  0.1762,  0.1828,  0.1887,  0.1947,  0.2002,  
    0.2062,  0.2118,  0.2173,  0.2222,  0.2278,  0.2327,  0.2377,  0.2428,  
    0.2473,  0.2522,  0.2567,  0.2612,  0.2657,  0.2697,  0.2743,  0.2783,  
    0.2822,  0.2862,  0.2903,  0.2943,  0.2982,  0.3018,  0.3058,  0.3092,  
    0.3128,  0.3167,  0.3203,  0.3237,  0.3268,  0.3302,  0.3338,  0.3368,  
    0.3402,  0.3433,  0.3462,  0.3492,  0.3528,  0.3557,  0.3587,  0.3613,  
    0.3643,  0.3672,  0.3702,  0.3728,  0.3758,  0.3783,  0.3812,  0.3837,  
    0.3862,  0.3887,  0.3918,  0.3943,  0.3968,  0.3993,  0.4017,  0.4042,  
    0.4067,  0.4087,  0.4112,  0.4137,  0.4157,  0.4182,  0.4207,  0.4227,  
    0.4252,  0.4272,  0.4293,  0.4317,  0.4338,  0.4358,  0.4383,  0.4403,  
    0.4423,  0.4442,  0.4462,  0.4482,  0.4502,  0.4523,  0.4543,  0.4563,  
    0.4582,  0.4602,  0.4622,  0.4638,  0.4658,  0.4678,  0.4697,  0.4712,  
    0.4733,  0.4753,  0.4767,  0.4787,  0.4803,  0.4823,  0.4837,  0.4857,  
    0.4873,  0.4888,  0.4908,  0.4922,  0.4942,  0.4958,  0.4972,  0.4988  
  }; 
  for (Int_t ilut = 0; ilut < kNlut; ilut++) {
    fLUT[ilut] = lut[ilut];
  }

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizerV1::ReadDigits()
{
  //
  // Reads the digits arrays from the input aliroot file
  //

  if (!fInputFile) {
    printf("AliTRDclusterizerV1::ReadDigits -- ");
    printf("No input file open\n");
    return kFALSE;
  }

  // Read in the digit arrays
  return (fDigitsManager->ReadDigits());  

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizerV1::MakeClusters()
{
  //
  // Generates the cluster.
  //

  Int_t row, col, time;

  if (fTRD->IsVersion() != 1) {
    printf("AliTRDclusterizerV1::MakeCluster -- ");
    printf("TRD must be version 1 (slow simulator).\n");
    return kFALSE; 
  }

  // Get the geometry
  AliTRDgeometry *geo = fTRD->GetGeometry();

  printf("AliTRDclusterizerV1::MakeCluster -- ");
  printf("Start creating clusters.\n");

  AliTRDdataArrayI *digits;
  AliTRDdataArrayI *track0;
  AliTRDdataArrayI *track1;
  AliTRDdataArrayI *track2; 

  // Threshold value for the maximum
  Int_t maxThresh = fClusMaxThresh;   
  // Threshold value for the digit signal
  Int_t sigThresh = fClusSigThresh;   

  // Iteration limit for unfolding procedure
  const Float_t kEpsilon = 0.01;             

  const Int_t   kNclus   = 3;  
  const Int_t   kNsig    = 5;
  const Int_t   kNtrack  = 3 * kNclus;

  // For the LUT
  const Float_t kLUTmin  = 0.106113;
  const Float_t kLUTmax  = 0.995415;

  Int_t   iType          = 0;
  Int_t   iUnfold        = 0;

  Float_t ratioLeft      = 1.0;
  Float_t ratioRight     = 1.0;

  Float_t padSignal[kNsig];   
  Float_t clusterSignal[kNclus];
  Float_t clusterPads[kNclus];   
  Int_t   clusterDigit[kNclus];
  Int_t   clusterTracks[kNtrack];   

  Int_t chamBeg = 0;
  Int_t chamEnd = AliTRDgeometry::Ncham();
  if (fTRD->GetSensChamber()  >= 0) {
    chamBeg = fTRD->GetSensChamber();
    chamEnd = chamBeg + 1;
  }
  Int_t planBeg = 0;
  Int_t planEnd = AliTRDgeometry::Nplan();
  if (fTRD->GetSensPlane()    >= 0) {
    planBeg = fTRD->GetSensPlane();
    planEnd = planBeg + 1;
  }
  Int_t sectBeg = 0;
  Int_t sectEnd = AliTRDgeometry::Nsect();

  // Start clustering in every chamber
  for (Int_t icham = chamBeg; icham < chamEnd; icham++) {
    for (Int_t iplan = planBeg; iplan < planEnd; iplan++) {
      for (Int_t isect = sectBeg; isect < sectEnd; isect++) {

        if (fTRD->GetSensSector() >= 0) {
          Int_t sens1 = fTRD->GetSensSector();
          Int_t sens2 = sens1 + fTRD->GetSensSectorRange();
          sens2 -= ((Int_t) (sens2 / AliTRDgeometry::Nsect())) 
                 * AliTRDgeometry::Nsect();
          if (sens1 < sens2) {
            if ((isect < sens1) || (isect >= sens2)) continue;
	  }
          else {
            if ((isect < sens1) && (isect >= sens2)) continue;
	  }
	}

        Int_t idet = geo->GetDetector(iplan,icham,isect);

        Int_t nClusters      = 0;
        Int_t nClusters2pad  = 0;
        Int_t nClusters3pad  = 0;
        Int_t nClusters4pad  = 0;
        Int_t nClusters5pad  = 0;
        Int_t nClustersLarge = 0;

        printf("AliTRDclusterizerV1::MakeCluster -- ");
        printf("Analyzing chamber %d, plane %d, sector %d.\n"
              ,icham,iplan,isect);

        Int_t nRowMax     = geo->GetRowMax(iplan,icham,isect);
        Int_t nColMax     = geo->GetColMax(iplan);
        Int_t nTimeBefore = geo->GetTimeBefore();
        Int_t nTimeTotal  = geo->GetTimeTotal();  

        // Get the digits
        digits = fDigitsManager->GetDigits(idet);
        digits->Expand();
        track0 = fDigitsManager->GetDictionary(idet,0);
        track0->Expand();
        track1 = fDigitsManager->GetDictionary(idet,1);
        track1->Expand();
        track2 = fDigitsManager->GetDictionary(idet,2); 
        track2->Expand();

        // Loop through the chamber and find the maxima 
        for ( row = 0;  row <  nRowMax;    row++) {
          for ( col = 2;  col <  nColMax;    col++) {
            for (time = 0; time < nTimeTotal; time++) {

              Int_t signalL = digits->GetDataUnchecked(row,col  ,time);
              Int_t signalM = digits->GetDataUnchecked(row,col-1,time);
              Int_t signalR = digits->GetDataUnchecked(row,col-2,time);
 
	      // Look for the maximum
              if (signalM >= maxThresh) {
                if (((signalL >= sigThresh) &&
                     (signalL <  signalM))  ||
                    ((signalR >= sigThresh) &&
                     (signalR <  signalM))) {
                  // Maximum found, mark the position by a negative signal
                  digits->SetDataUnchecked(row,col-1,time,-signalM);
		}
	      }

            }  
          }    
        }      

        // Now check the maxima and calculate the cluster position
        for ( row = 0;  row <  nRowMax  ;  row++) {
          for (time = 0; time < nTimeTotal; time++) {
            for ( col = 1;  col <  nColMax-1;  col++) {

              // Maximum found ?             
              if (digits->GetDataUnchecked(row,col,time) < 0) {

                Int_t iPad;
                for (iPad = 0; iPad < kNclus; iPad++) {
                  Int_t iPadCol = col - 1 + iPad;
                  clusterSignal[iPad]     = TMath::Abs(digits->GetDataUnchecked(row
                                                                               ,iPadCol
                                                                               ,time));
                  clusterDigit[iPad]      = digits->GetIndexUnchecked(row,iPadCol,time);
                  clusterTracks[3*iPad  ] = track0->GetDataUnchecked(row,iPadCol,time) - 1;
		  clusterTracks[3*iPad+1] = track1->GetDataUnchecked(row,iPadCol,time) - 1;
		  clusterTracks[3*iPad+2] = track2->GetDataUnchecked(row,iPadCol,time) - 1;
                }

		// Count the number of pads in the cluster
                Int_t nPadCount = 0;
                Int_t ii        = 0;
                while (TMath::Abs(digits->GetDataUnchecked(row,col-ii  ,time))
                                                                  >= sigThresh) {
                  nPadCount++;
                  ii++;
                  if (col-ii   <        0) break;
		}
                ii = 0;
                while (TMath::Abs(digits->GetDataUnchecked(row,col+ii+1,time))
                                                                  >= sigThresh) {
                  nPadCount++;
                  ii++;
                  if (col+ii+1 >= nColMax) break;
		}

                nClusters++;
                switch (nPadCount) {
                case 2:
                  iType = 0;
                  nClusters2pad++;
                  break;
                case 3:
                  iType = 1;
                  nClusters3pad++;
                  break;
                case 4:
                  iType = 2;
                  nClusters4pad++;
                  break;
                case 5:
                  iType = 3;
                  nClusters5pad++;
                  break;
                default:
                  iType = 4;
                  nClustersLarge++;
                  break;
		};

		// Don't analyze large clusters
                //if (iType == 4) continue;

                // Look for 5 pad cluster with minimum in the middle
                Bool_t fivePadCluster = kFALSE;
                if (col < nColMax-3) {
                  if (digits->GetDataUnchecked(row,col+2,time) < 0) {
                    fivePadCluster = kTRUE;
	  	  }
                  if ((fivePadCluster) && (col < nColMax-5)) {
                    if (digits->GetDataUnchecked(row,col+4,time) >= sigThresh) {
                      fivePadCluster = kFALSE;
		    }
		  }
                  if ((fivePadCluster) && (col >         1)) {
                    if (digits->GetDataUnchecked(row,col-2,time) >= sigThresh) {
                      fivePadCluster = kFALSE;
		    }
		  }
		}

		// 5 pad cluster
                // Modify the signal of the overlapping pad for the left part 
		// of the cluster which remains from a previous unfolding
                if (iUnfold) {
                  clusterSignal[0] *= ratioLeft;
                  iType   = 3;
                  iUnfold = 0;
		}

		// Unfold the 5 pad cluster
                if (fivePadCluster) {
                  for (iPad = 0; iPad < kNsig; iPad++) {
                    padSignal[iPad] = TMath::Abs(digits->GetDataUnchecked(row
                                                                         ,col-1+iPad
                                                                         ,time));
                  }
                  // Unfold the two maxima and set the signal on 
                  // the overlapping pad to the ratio
                  ratioRight        = Unfold(kEpsilon,padSignal);
                  ratioLeft         = 1.0 - ratioRight; 
                  clusterSignal[2] *= ratioRight;
                  iType   = 3;
                  iUnfold = 1;
                }

                Float_t clusterCharge = clusterSignal[0]
                                      + clusterSignal[1]
                                      + clusterSignal[2];
                
		// The position of the cluster
                clusterPads[0] = row + 0.5;
		// Take the shift of the additional time bins into account
                clusterPads[2] = time - nTimeBefore + 0.5;

                if (fUseLUT) {

  		  // Calculate the position of the cluster by using the
		  // lookup table method
                  Float_t ratioLUT;
                  Float_t signLUT;
                  Float_t lut = 0.0;
                  if (clusterSignal[0] > clusterSignal[2]) {
                    ratioLUT = clusterSignal[0] / clusterSignal[1];
                    signLUT  = -1.0;
		  }
                  else {
                    ratioLUT = clusterSignal[2] / clusterSignal[1];
                    signLUT  =  1.0;
		  }
                  if      (ratioLUT < kLUTmin) {
                    lut = 0.0;
		  }
                  else if (ratioLUT > kLUTmax) {
                    lut = 0.5;
		  }
                  else {
                    Int_t indexLUT = TMath::Nint ((kNlut-1) * (ratioLUT - kLUTmin)  
						            / (kLUTmax  - kLUTmin)); 
                    lut = fLUT[indexLUT];
		  }
                  clusterPads[1] = col + 0.5 + signLUT * lut;

		}
		else {

  		  // Calculate the position of the cluster by using the
		  // center of gravity method
                  clusterPads[1] = col + 0.5 
                                 + (clusterSignal[2] - clusterSignal[0]) 
		                 / clusterCharge;

		}

                if (fVerbose) {
                  printf("-----------------------------------------------------------\n");
                  printf("Create cluster no. %d\n",nClusters);
                  printf("Position: row = %f, col = %f, time = %f\n",clusterPads[0]
			                                            ,clusterPads[1]
                                                                    ,clusterPads[2]);
                  printf("Indices: %d, %d, %d\n",clusterDigit[0]
			                        ,clusterDigit[1]
                                                ,clusterDigit[2]);
                  printf("Total charge = %f\n",clusterCharge);
                  printf("Tracks: pad0 %d, %d, %d\n",clusterTracks[0]
			                            ,clusterTracks[1]
                                                    ,clusterTracks[2]);
                  printf("        pad1 %d, %d, %d\n",clusterTracks[3]
			                            ,clusterTracks[4]
                                                    ,clusterTracks[5]);
                  printf("        pad2 %d, %d, %d\n",clusterTracks[6]
			                            ,clusterTracks[7]
                                                    ,clusterTracks[8]);
                  printf("Type = %d, Number of pads = %d\n",iType,nPadCount);
                }

                // Add the cluster to the output array 
                fTRD->AddCluster(clusterPads
                                ,clusterDigit
                                ,idet
                                ,clusterCharge
                                ,clusterTracks
                                ,iType);

              }
            } 
          }   
        }     

	// Compress the arrays
        digits->Compress(1,0);
        track0->Compress(1,0);
        track1->Compress(1,0);
        track2->Compress(1,0);

        // Write the cluster and reset the array
	WriteClusters(idet);
	fTRD->ResetRecPoints();

        printf("AliTRDclusterizerV1::MakeCluster -- ");
        printf("Found %d clusters in total.\n"
              ,nClusters);
        printf("                                    2pad:  %d\n",nClusters2pad);
        printf("                                    3pad:  %d\n",nClusters3pad);
        printf("                                    4pad:  %d\n",nClusters4pad);
        printf("                                    5pad:  %d\n",nClusters5pad);
        printf("                                    Large: %d\n",nClustersLarge);

      }    
    }      
  }        

  printf("AliTRDclusterizerV1::MakeCluster -- ");
  printf("Done.\n");

  return kTRUE;

}

//_____________________________________________________________________________
Float_t AliTRDclusterizerV1::Unfold(Float_t eps, Float_t* padSignal)
{
  //
  // Method to unfold neighbouring maxima.
  // The charge ratio on the overlapping pad is calculated
  // until there is no more change within the range given by eps.
  // The resulting ratio is then returned to the calling method.
  //

  Int_t   itStep            = 0;      // Count iteration steps

  Float_t ratio             = 0.5;    // Start value for ratio
  Float_t prevRatio         = 0;      // Store previous ratio

  Float_t newLeftSignal[3]  = {0};    // Array to store left cluster signal
  Float_t newRightSignal[3] = {0};    // Array to store right cluster signal

  // Start the iteration
  while ((TMath::Abs(prevRatio - ratio) > eps) && (itStep < 10)) {

    itStep++;
    prevRatio = ratio;

    // Cluster position according to charge ratio
    Float_t maxLeft  = (ratio*padSignal[2] - padSignal[0]) 
                     / (padSignal[0] + padSignal[1] + ratio*padSignal[2]);
    Float_t maxRight = (padSignal[4] - (1-ratio)*padSignal[2]) 
                     / ((1-ratio)*padSignal[2] + padSignal[3] + padSignal[4]);

    // Set cluster charge ratio
    Float_t ampLeft  = padSignal[1];
    Float_t ampRight = padSignal[3];

    // Apply pad response to parameters
    newLeftSignal[0]  = ampLeft  * PadResponse(-1 - maxLeft);
    newLeftSignal[1]  = ampLeft  * PadResponse( 0 - maxLeft);
    newLeftSignal[2]  = ampLeft  * PadResponse( 1 - maxLeft);

    newRightSignal[0] = ampRight * PadResponse(-1 - maxRight);
    newRightSignal[1] = ampRight * PadResponse( 0 - maxRight);
    newRightSignal[2] = ampRight * PadResponse( 1 - maxRight);

    // Calculate new overlapping ratio
    ratio = TMath::Min((Float_t)1.0,newLeftSignal[2] / 
                          (newLeftSignal[2] + newRightSignal[0]));

  }

  return ratio;

}

//_____________________________________________________________________________
Float_t AliTRDclusterizerV1::PadResponse(Float_t x)
{
  //
  // The pad response for the chevron pads. 
  // We use a simple Gaussian approximation which should be good
  // enough for our purpose.
  // Updated for new PRF 1/5/01.
  //

  // The parameters for the response function
  const Float_t kA  =  0.8303; 
  const Float_t kB  = -0.00392;
  const Float_t kC  =  0.472 * 0.472;
  const Float_t kD  =  2.19;

  Float_t pr = kA * (kB + TMath::Exp(-TMath::Power(x*x,kD) / (2.*kC)));

  return (pr);

}
