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
// TRD cluster finder for the slow simulator. 
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TF1.h>
#include <TTree.h>
#include <TH1.h>
#include <TFile.h>

#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliLoader.h"

#include "AliTRDclusterizerMI.h"
#include "AliTRDmatrix.h"
#include "AliTRDgeometry.h"
#include "AliTRDdataArrayF.h"
#include "AliTRDdataArrayI.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDparameter.h"
#include "AliTRDclusterMI.h"
#include "AliTRDpadPlane.h"
#include "AliTRDcalibDB.h"
#include "AliTRDRecParam.h"
#include "AliTRDCommonParam.h"

ClassImp(AliTRDclusterizerMI)

//_____________________________________________________________________________
AliTRDclusterizerMI::AliTRDclusterizerMI():AliTRDclusterizerV1()
{
  //
  // AliTRDclusterizerMI default constructor
  //
}

//_____________________________________________________________________________
AliTRDclusterizerMI::AliTRDclusterizerMI(const Text_t* name, const Text_t* title)
                    :AliTRDclusterizerV1(name,title)
{
  //
  // AliTRDclusterizerMI default constructor
  //
}

//_____________________________________________________________________________
AliTRDclusterizerMI::~AliTRDclusterizerMI()
{
  //
  // AliTRDclusterizerMI destructor
  //
}

//_____________________________________________________________________________
AliTRDclusterMI *  AliTRDclusterizerMI::AddCluster()
{
  //
  // Adds cluster
  //

  AliTRDclusterMI *c = new AliTRDclusterMI();
  fClusterContainer->Add(c);
  return c;

}

//_____________________________________________________________________________
void AliTRDclusterizerMI::SetCluster(AliTRDclusterMI * cl, Double_t *pos, Int_t timebin
                                   , Int_t det, Double_t amp
		                   , Int_t *tracks, Double_t *sig, Int_t iType
                                   , Double_t sigmay, Double_t relpad)
{
  //
  // Sets cluster
  //

  cl->SetDetector(det);
  cl->AddTrackIndex(tracks);
  cl->SetQ(amp);
  cl->SetX(pos[2]);
  cl->SetY(pos[0]);
  cl->SetZ(pos[1]);
  cl->SetSigmaY2(sig[0]);   
  cl->SetSigmaZ2(sig[1]);
  cl->SetLocalTimeBin(timebin);
  cl->SetNPads(iType);
  cl->SetRelPos(relpad);
  cl->fRmsY = sigmay;
}

//_____________________________________________________________________________
void AliTRDclusterizerMI::MakeCluster(Double_t * padSignal, Double_t * pos
                                    , Double_t &sigma, Double_t & relpad)
{
  //
  // Does something with the cluster  
  //

  Double_t sum   = 0;
  Double_t sumx  = 0;
  Double_t sumx2 = 0;
  Double_t signal[3]={padSignal[0],padSignal[1],padSignal[2]};
  if ( signal[0]<2){
    signal[0] = 0.015*(signal[1]*signal[1])/(signal[2]+0.5);
    if (signal[0]>2) signal[0]=2;
  }
  if ( signal[2]<2){
    signal[2] = 0.015*(signal[1]*signal[1])/(signal[0]+0.5);
    if (signal[2]>2) signal[2]=2;
  }

  for (Int_t i=-1;i<=1;i++){
    sum   +=signal[i+1];
    sumx  +=signal[i+1]*float(i);
    sumx2 +=signal[i+1]*float(i)*float(i);
  }
  
  pos[0] = sumx/sum;
  sigma  = sumx2/sum-pos[0]*pos[0];
  relpad = pos[0];
}

//_____________________________________________________________________________
Bool_t AliTRDclusterizerMI::MakeClusters()
{
  //
  // Generates the cluster.
  //

  //////////////////////
  //STUPIDITY to be fixed later
  fClusterContainer = RecPoints();

  Int_t row, col, time;

  /*
  if (fTRD->IsVersion() != 1) {
    printf("<AliTRDclusterizerMI::MakeCluster> ");
    printf("TRD must be version 1 (slow simulator).\n");
    return kFALSE; 
  }
  */

  // Get the geometry
  AliTRDgeometry *geo = AliTRDgeometry::GetGeometry(fRunLoader);

  // Create a default parameter class if none is defined
  if (!fPar) {
    fPar = new AliTRDparameter("TRDparameter","Standard TRD parameter");
    printf("<AliTRDclusterizerMI::MakeCluster> ");
    printf("Create the default parameter object.\n");
  }

  AliTRDcalibDB* calibration = AliTRDcalibDB::Instance();
  if (!calibration)
  {
    printf("<AliTRDclusterizerMI::MakeCluster> ");
    printf("ERROR getting instance of AliTRDcalibDB");
    return kFALSE;  
  }
  
  AliTRDRecParam* recParam = AliTRDRecParam::Instance();
  if (!recParam)
  {
    printf("<AliTRDclusterizerMI::MakeCluster> ");
    printf("ERROR getting instance of AliTRDRecParam");
    return kFALSE;  
  }
  
  AliTRDCommonParam* commonParam = AliTRDCommonParam::Instance();
  if (!commonParam)
  {
    printf("<AliTRDdigitizer::MakeDigits> ");
    printf("Could not get common params\n");
    return kFALSE;
  }
      
  // Half of ampl.region
  const Float_t kAmWidth = AliTRDgeometry::AmThick()/2.; 

  if (fVerbose > 0) {
    //printf("<AliTRDclusterizerMI::MakeCluster> ");
    //printf("OmegaTau = %f \n",omegaTau);
    printf("<AliTRDclusterizerMI::MakeCluster> ");
    printf("Start creating clusters.\n");
  } 

  AliTRDdataArrayI *digits;
  AliTRDdataArrayI *track0;
  AliTRDdataArrayI *track1;
  AliTRDdataArrayI *track2; 

  // Threshold value for the maximum
  Int_t maxThresh = recParam->GetClusMaxThresh();   
  // Threshold value for the digit signal
  Int_t sigThresh = recParam->GetClusSigThresh();   

  // Iteration limit for unfolding procedure
  const Double_t kEpsilon = 0.01;             

  const Int_t    kNclus   = 3;  
  const Int_t    kNsig    = 5;
  const Int_t    kNtrack  = 3 * kNclus;

  Int_t    iType          = 0;
  Int_t    iUnfold        = 0;

  Double_t ratioLeft      = 1.0;
  Double_t ratioRight     = 1.0;

  Double_t padSignal[kNsig];   
  Double_t clusterSignal[kNclus];
  Double_t clusterPads[kNclus];   
  Int_t    clusterDigit[kNclus];
  Int_t    clusterTracks[kNtrack];   

  Int_t chamBeg = 0;
  Int_t chamEnd = AliTRDgeometry::Ncham();
  Int_t planBeg = 0;
  Int_t planEnd = AliTRDgeometry::Nplan();
  Int_t sectBeg = 0;
  Int_t sectEnd = AliTRDgeometry::Nsect();

  // Start clustering in every chamber
  for (Int_t icham = chamBeg; icham < chamEnd; icham++) {
    for (Int_t iplan = planBeg; iplan < planEnd; iplan++) {
      for (Int_t isect = sectBeg; isect < sectEnd; isect++) {

        Int_t idet = geo->GetDetector(iplan,icham,isect);

        Int_t nClusters      = 0;
        Int_t nClusters2pad  = 0;
        Int_t nClusters3pad  = 0;
        Int_t nClusters4pad  = 0;
        Int_t nClusters5pad  = 0;
        Int_t nClustersLarge = 0;

        if (fVerbose > 0) {
          printf("<AliTRDclusterizerMI::MakeCluster> ");
          printf("Analyzing chamber %d, plane %d, sector %d.\n"
                ,icham,iplan,isect);
	}

        Int_t   nRowMax     = commonParam->GetRowMax(iplan,icham,isect);
        Int_t   nColMax     = commonParam->GetColMax(iplan);
        Int_t   nTimeTotal  = calibration->GetNumberOfTimeBins();

        AliTRDpadPlane *padPlane = commonParam->GetPadPlane(iplan,icham);

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

              Int_t signalL = TMath::Abs(digits->GetDataUnchecked(row,col  ,time));
              Int_t signalM = TMath::Abs(digits->GetDataUnchecked(row,col-1,time));
              Int_t signalR = TMath::Abs(digits->GetDataUnchecked(row,col-2,time));
 
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
          for ( col = 1;  col <  nColMax-1;  col++) {
            for (time = 0; time < nTimeTotal; time++) {

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
                  ratioRight        = Unfold(kEpsilon,iplan,padSignal);
                  ratioLeft         = 1.0 - ratioRight; 
                  clusterSignal[2] *= ratioRight;
                  iType   = 3;
                  iUnfold = 1;
                }

                Double_t clusterCharge = clusterSignal[0]
                                       + clusterSignal[1]
                                       + clusterSignal[2];
                
		// The position of the cluster
                clusterPads[0] = row + 0.5;
		// Take the shift of the additional time bins into account
                clusterPads[2] = time + 0.5;
                
                // correct for t0
                clusterPads[2] -= calibration->GetT0(idet, col, row);

                if (recParam->LUTOn()) {

  		  // Calculate the position of the cluster by using the
		  // lookup table method
//                   clusterPads[1] = col + 0.5
//                                  + fPar->LUTposition(iplan,clusterSignal[0]
//                                                           ,clusterSignal[1]
// 					                  ,clusterSignal[2]);
                  clusterPads[1] = 0.5
                      + recParam->LUTposition(iplan,clusterSignal[0]
                                                          ,clusterSignal[1]
					                  ,clusterSignal[2]);

		}
		else {

  		  // Calculate the position of the cluster by using the
		  // center of gravity method
//                   clusterPads[1] = col + 0.5 
//                                  + (clusterSignal[2] - clusterSignal[0]) 
// 		                 / clusterCharge;
                  clusterPads[1] = 0.5 
                                 + (clusterSignal[2] - clusterSignal[0]) 
		                 / clusterCharge;

		}

                Double_t q0 = clusterSignal[0];
		Double_t q1 = clusterSignal[1];
                Double_t q2 = clusterSignal[2];
                Double_t clusterSigmaY2 = (q1*(q0+q2)+4*q0*q2) /
                                         (clusterCharge*clusterCharge);

		// Calculate the position and the error
                Float_t vdrift = calibration->GetVdrift(idet, col, row);
                Double_t clusterPos[3];
//                 clusterPos[0] = clusterPads[1] * colSize + col0;
//                 clusterPos[1] = clusterPads[0] * rowSize + row0;
                clusterPos[0] = padPlane->GetColPos(col) - clusterPads[1];
                clusterPos[1] = padPlane->GetRowPos(row) - clusterPads[0];
                clusterPos[2] = CalcXposFromTimebin(clusterPads[2], vdrift);
                Double_t clusterSig[2];
                Double_t colSize = padPlane->GetColSize(col);
                Double_t rowSize = padPlane->GetRowSize(row);
                clusterSig[0] = (clusterSigmaY2 + 1./12.) * colSize*colSize;
                clusterSig[1] = rowSize * rowSize / 12.;

                Int_t    local_time_bin = (Int_t) clusterPads[2];
                
                // Correct for ExB displacement
                if (commonParam->ExBOn()) { 
                  Float_t omegaTau = calibration->GetOmegaTau(vdrift);
                  Double_t timeBinSize = vdrift / calibration->GetSamplingFrequency();
                  
                  Double_t driftLength    = local_time_bin * timeBinSize + kAmWidth;
                  Double_t deltaY         = omegaTau * driftLength;
                  clusterPos[1]           = clusterPos[1] - deltaY;
                }

		//
		//
		AliTRDclusterMI * cluster = AddCluster();
		Double_t sigma, relpos;
		MakeCluster(clusterSignal, clusterPos, sigma,relpos);

// 		   clusterPos[0] = clusterPads[1] * colSize + col0;
//                 clusterPos[1] = clusterPads[0] * rowSize + row0;
                clusterPos[2] = CalcXposFromTimebin(clusterPads[2], vdrift);
                clusterPos[0] = padPlane->GetColPos(col) - clusterPads[1];
                clusterPos[1] = padPlane->GetRowPos(row) - clusterPads[0];
                SetCluster(cluster, clusterPos, (Int_t) clusterPads[2],idet,clusterCharge,clusterTracks,clusterSig,iType,sigma,relpos);
                // Add the cluster to the output array 
		//                fTRD->AddCluster(clusterPos
                //                ,idet
                //                ,clusterCharge
                //                ,clusterTracks
		//		,clusterSig
                //                ,iType);

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
	ResetRecPoints();

        if (fVerbose > 0) {
          printf("<AliTRDclusterizerMI::MakeCluster> ");
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
  }        

  if (fVerbose > 0) {
    printf("<AliTRDclusterizerMI::MakeCluster> ");
    printf("Done.\n");
  }

  return kTRUE;

}

//_____________________________________________________________________________
Double_t AliTRDclusterizerMI::Unfold(Double_t eps, Int_t plane, Double_t* padSignal)
{
  //
  // Method to unfold neighbouring maxima.
  // The charge ratio on the overlapping pad is calculated
  // until there is no more change within the range given by eps.
  // The resulting ratio is then returned to the calling method.
  //

  AliTRDCommonParam* commonParam = AliTRDCommonParam::Instance();
  if (!commonParam)
  {
    printf("<AliTRDdigitizer::MakeDigits> ");
    printf("Could not get common params\n");
    return kFALSE;
  }
  
  Int_t   irc                = 0;
  Int_t   itStep             = 0;      // Count iteration steps

  Double_t ratio             = 0.5;    // Start value for ratio
  Double_t prevRatio         = 0;      // Store previous ratio

  Double_t newLeftSignal[3]  = {0};    // Array to store left cluster signal
  Double_t newRightSignal[3] = {0};    // Array to store right cluster signal
  Double_t newSignal[3]      = {0};

  // Start the iteration
  while ((TMath::Abs(prevRatio - ratio) > eps) && (itStep < 10)) {

    itStep++;
    prevRatio = ratio;

    // Cluster position according to charge ratio
    Double_t maxLeft  = (ratio*padSignal[2] - padSignal[0]) 
                      / (padSignal[0] + padSignal[1] + ratio*padSignal[2]);
    Double_t maxRight = (padSignal[4] - (1-ratio)*padSignal[2]) 
                      / ((1-ratio)*padSignal[2] + padSignal[3] + padSignal[4]);

    // Set cluster charge ratio
    irc = commonParam->PadResponse(1.0,maxLeft ,plane,newSignal);
    Double_t ampLeft  = padSignal[1] / newSignal[1];
    irc = commonParam->PadResponse(1.0,maxRight,plane,newSignal);
    Double_t ampRight = padSignal[3] / newSignal[1];

    // Apply pad response to parameters
    irc = commonParam->PadResponse(ampLeft ,maxLeft ,plane,newLeftSignal );
    irc = commonParam->PadResponse(ampRight,maxRight,plane,newRightSignal);

    // Calculate new overlapping ratio
    ratio = TMath::Min((Double_t)1.0,newLeftSignal[2] / 
                          (newLeftSignal[2] + newRightSignal[0]));

  }

  return ratio;

}



