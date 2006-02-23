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
#include "AliRawReader.h"

#include "AliTRDclusterizerV1.h"
#include "AliTRDmatrix.h"
#include "AliTRDgeometry.h"
#include "AliTRDdataArrayF.h"
#include "AliTRDdataArrayI.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDparameter.h"
#include "AliTRDpadPlane.h"
#include "AliTRDrawData.h"
#include "AliTRDcalibDB.h"
#include "AliTRDRecParam.h"
#include "AliTRDCommonParam.h"
#include "AliTRDcluster.h"

ClassImp(AliTRDclusterizerV1)

//_____________________________________________________________________________
AliTRDclusterizerV1::AliTRDclusterizerV1():AliTRDclusterizer()
{
  //
  // AliTRDclusterizerV1 default constructor
  //

  fDigitsManager = 0;

}

//_____________________________________________________________________________
AliTRDclusterizerV1::AliTRDclusterizerV1(const Text_t* name, const Text_t* title)
                    :AliTRDclusterizer(name,title)
{
  //
  // AliTRDclusterizerV1 default constructor
  //

  fDigitsManager = new AliTRDdigitsManager();
  fDigitsManager->CreateArrays();

}

//_____________________________________________________________________________
AliTRDclusterizerV1::AliTRDclusterizerV1(const AliTRDclusterizerV1 &c)
:AliTRDclusterizer(c)
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
    fDigitsManager = NULL;
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
void AliTRDclusterizerV1::Copy(TObject &c) const
{
  //
  // Copy function
  //

  ((AliTRDclusterizerV1 &) c).fDigitsManager = 0;

  AliTRDclusterizer::Copy(c);

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizerV1::ReadDigits()
{
  //
  // Reads the digits arrays from the input aliroot file
  //

  if (!fRunLoader) {
    printf("<AliTRDclusterizerV1::ReadDigits> ");
    printf("No input file open\n");
    return kFALSE;
  }
  AliLoader* loader = fRunLoader->GetLoader("TRDLoader");
  if (!loader->TreeD()) loader->LoadDigits();

  // Read in the digit arrays
  return (fDigitsManager->ReadDigits(loader->TreeD()));

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizerV1::ReadDigits(AliRawReader* rawReader)
{
  //
  // Reads the digits arrays from the ddl file
  //

  AliTRDrawData *raw = new AliTRDrawData();
  raw->SetDebug(1);

  fDigitsManager = raw->Raw2Digits(rawReader);

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDclusterizerV1::MakeClusters()
{
  //
  // Generates the cluster.
  //

  Int_t row, col, time;

  /*
  if (fTRD->IsVersion() != 1) {
    printf("<AliTRDclusterizerV1::MakeCluster> ");
    printf("TRD must be version 1 (slow simulator).\n");
    return kFALSE; 
  }
  */

  // Get the geometry
  AliTRDgeometry *geo = AliTRDgeometry::GetGeometry(fRunLoader);

  // Create a default parameter class if none is defined
  if (!fPar) {
    fPar = new AliTRDparameter("TRDparameter","Standard TRD parameter");
    printf("<AliTRDclusterizerV1::MakeCluster> ");
    printf("Create the default parameter object.\n");
  }
  fPar->Init();
  
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
    
  //Float_t timeBinSize = fPar->GetDriftVelocity()
  //                    / fPar->GetSamplingFrequency();
  // Half of ampl.region
  //  const Float_t kAmWidth = AliTRDgeometry::AmThick()/2.; 

  //Float_t omegaTau = fPar->GetOmegaTau();
  if (fVerbose > 0) {
    //printf("<AliTRDclusterizerV1::MakeCluster> ");
    //printf("OmegaTau = %f \n",omegaTau);
    printf("<AliTRDclusterizerV1::MakeCluster> ");
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
  const Float_t kEpsilon = 0.01;             

  const Int_t   kNclus   = 3;  
  const Int_t   kNsig    = 5;
  const Int_t   kNtrack  = 3 * kNclus;

  Int_t    iType         = 0;
  Int_t    iUnfold       = 0;  
  Double_t ratioLeft     = 1.0;
  Double_t ratioRight    = 1.0;

  //
  Double_t padSignal[kNsig];   
  Double_t clusterSignal[kNclus];
  Double_t clusterPads[kNclus];   
  Int_t    clusterDigit[kNclus];
  Int_t    clusterTracks[kNtrack];   

  Int_t    chamBeg = 0;
  Int_t    chamEnd = AliTRDgeometry::Ncham();
  Int_t    planBeg = 0;
  Int_t    planEnd = AliTRDgeometry::Nplan();
  Int_t    sectBeg = 0;
  Int_t    sectEnd = AliTRDgeometry::Nsect();

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
          printf("<AliTRDclusterizerV1::MakeCluster> ");
          printf("Analyzing chamber %d, plane %d, sector %d.\n"
                ,icham,iplan,isect);
	}

        Int_t    nRowMax     = commonParam->GetRowMax(iplan,icham,isect);
        Int_t    nColMax     = commonParam->GetColMax(iplan);
        Int_t    nTimeTotal  = calibration->GetNumberOfTimeBins();

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
            //for ( col = 4;  col <  nColMax-2;    col++) {
            for (time = 0; time < nTimeTotal; time++) {

              Int_t signalL = TMath::Abs(digits->GetDataUnchecked(row,col  ,time));
              Int_t signalM = TMath::Abs(digits->GetDataUnchecked(row,col-1,time));
              Int_t signalR = TMath::Abs(digits->GetDataUnchecked(row,col-2,time));
 
// 	      // Look for the maximum
//               if (signalM >= maxThresh) {
//                 if (((signalL >= sigThresh) &&
//                      (signalL <  signalM))  ||
//                     ((signalR >= sigThresh) &&
//                      (signalR <  signalM))) {
//                   // Maximum found, mark the position by a negative signal
//                   digits->SetDataUnchecked(row,col-1,time,-signalM);
// 		}
// 	      }
	      // Look for the maximum
              if (signalM >= maxThresh) {
                if ( (TMath::Abs(signalL)<=signalM) && (TMath::Abs(signalR)<=signalM) && 
		     (TMath::Abs(signalL)+TMath::Abs(signalR))>sigThresh ) {
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
                  iType   = 5;
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
                  iType   = 5;
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
                  clusterPads[1] =
                      recParam->LUTposition(iplan,clusterSignal[0]
                                                          ,clusterSignal[1]
					                  ,clusterSignal[2]);
		}
		else {
  		  // Calculate the position of the cluster by using the
		  // center of gravity method
		  for (Int_t i=0;i<5;i++) padSignal[i]=0;
		  padSignal[2] = TMath::Abs(digits->GetDataUnchecked(row,col,time));   // central  pad
		  padSignal[1] = TMath::Abs(digits->GetDataUnchecked(row,col-1,time)); // left     pad
		  padSignal[3] = TMath::Abs(digits->GetDataUnchecked(row,col+1,time)); // right    pad
		  if (col>2 &&TMath::Abs(digits->GetDataUnchecked(row,col-2,time)<padSignal[1])){
		    padSignal[0] = TMath::Abs(digits->GetDataUnchecked(row,col-2,time));
		  }
		  if (col<nColMax-3 &&TMath::Abs(digits->GetDataUnchecked(row,col+2,time)<padSignal[3])){
		    padSignal[4] = TMath::Abs(digits->GetDataUnchecked(row,col+2,time));
		  }		  
		  clusterPads[1] =  GetCOG(padSignal);
		  //Double_t check = fPar->LUTposition(iplan,clusterSignal[0]
                  //                                        ,clusterSignal[1]
		  //	  		                    ,clusterSignal[2]);
		  //		  Float_t diff = clusterPads[1] -  check;

		}

                Double_t q0 = clusterSignal[0];
		Double_t q1 = clusterSignal[1];
                Double_t q2 = clusterSignal[2];
                Double_t clusterSigmaY2 = (q1*(q0+q2)+4*q0*q2) /
                                          (clusterCharge*clusterCharge);

                Float_t vdrift = calibration->GetVdrift(idet, col, row);  
		
                // Calculate the position and the error
                Double_t colSize = padPlane->GetColSize(col);
                Double_t rowSize = padPlane->GetRowSize(row);
                Double_t clusterPos[3];
		clusterPos[0] = padPlane->GetColPos(col) - (clusterPads[1]+0.5)*colSize;  // MI change
		clusterPos[1] = padPlane->GetRowPos(row) -0.5*rowSize; //MI change
                clusterPos[2] = CalcXposFromTimebin(clusterPads[2], vdrift);
                Double_t clusterSig[2];
                clusterSig[0] = (clusterSigmaY2 + 1./12.) * colSize*colSize;
                clusterSig[1] = rowSize * rowSize / 12.;                                       
                
                
                // Add the cluster to the output array 
                AliTRDcluster * cluster = AddCluster(clusterPos
                          ,(Int_t) clusterPads[2]
                          ,idet
			  ,clusterCharge
			  ,clusterTracks
			  ,clusterSig
			   ,iType,clusterPads[1]);
		//
		//
		Short_t signals[7]={0,0,0,0,0,0,0};
		for (Int_t jPad = col-3;jPad<=col+3;jPad++){
		  if (jPad<0 ||jPad>=nColMax-1) continue;
		  signals[jPad-col+3] =  TMath::Abs(digits->GetDataUnchecked(row,jPad,time));
		}
		cluster->SetSignals(signals);
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
      }    
    }      
  }        

  if (fVerbose > 0) {
    printf("<AliTRDclusterizerV1::MakeCluster> ");
    printf("Done.\n");
  }

  return kTRUE;

}

Double_t AliTRDclusterizerV1::GetCOG(Double_t signal[5])
{
  //
  // get COG position
  // used for clusters with more than 3 pads - where LUT not applicable
  Double_t sum = signal[0]+signal[1]+signal[2]+signal[3]+signal[4];
  Double_t res = (0.0*(-signal[0]+signal[4])+(-signal[1]+signal[3]))/sum;
  return res;		  
}



//_____________________________________________________________________________
Double_t AliTRDclusterizerV1::Unfold(Double_t eps, Int_t plane, Double_t* padSignal)
{
  //
  // Method to unfold neighbouring maxima.
  // The charge ratio on the overlapping pad is calculated
  // until there is no more change within the range given by eps.
  // The resulting ratio is then returned to the calling method.
  //

  AliTRDcalibDB* calibration = AliTRDcalibDB::Instance();
  if (!calibration)
  {
    printf("<AliTRDclusterizerMI::Unfold> ");
    printf("ERROR getting instance of AliTRDcalibDB");
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
    irc = calibration->PadResponse(1.0,maxLeft ,plane,newSignal);
    Double_t ampLeft  = padSignal[1] / newSignal[1];
    irc = calibration->PadResponse(1.0,maxRight,plane,newSignal);
    Double_t ampRight = padSignal[3] / newSignal[1];

    // Apply pad response to parameters
    irc = calibration->PadResponse(ampLeft ,maxLeft ,plane,newLeftSignal );
    irc = calibration->PadResponse(ampRight,maxRight,plane,newRightSignal);

    // Calculate new overlapping ratio
    ratio = TMath::Min((Double_t)1.0,newLeftSignal[2] / 
                          (newLeftSignal[2] + newRightSignal[0]));

  }

  return ratio;

}

