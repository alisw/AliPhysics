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

#include "AliTRDclusterizerV1.h"
#include "AliTRDmatrix.h"
#include "AliTRDgeometry.h"
#include "AliTRDdigitizer.h"
#include "AliTRDrecPoint.h"
#include "AliTRDdataArrayF.h"

ClassImp(AliTRDclusterizerV1)

//_____________________________________________________________________________
AliTRDclusterizerV1::AliTRDclusterizerV1():AliTRDclusterizer()
{
  //
  // AliTRDclusterizerV1 default constructor
  //

  fDigitsManager = NULL;

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

  ((AliTRDclusterizerV1 &) c).fClusMaxThresh = fClusMaxThresh;
  ((AliTRDclusterizerV1 &) c).fClusSigThresh = fClusSigThresh;
  ((AliTRDclusterizerV1 &) c).fClusMethod    = fClusMethod;
  ((AliTRDclusterizerV1 &) c).fDigitsManager = NULL;

  AliTRDclusterizer::Copy(c);

}

//_____________________________________________________________________________
void AliTRDclusterizerV1::Init()
{
  //
  // Initializes the cluster finder
  //

  // The default parameter for the clustering
  fClusMaxThresh = 5.0;
  fClusSigThresh = 2.0;
  fClusMethod    = 1;

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
Bool_t AliTRDclusterizerV1::MakeCluster()
{
  //
  // Generates the cluster.
  //

  Int_t row, col, time;

  // Get the pointer to the detector class and check for version 1
  AliTRD *trd = (AliTRD*) gAlice->GetDetector("TRD");
  if (trd->IsVersion() != 1) {
    printf("AliTRDclusterizerV1::MakeCluster -- ");
    printf("TRD must be version 1 (slow simulator).\n");
    return kFALSE; 
  }

  // Get the geometry
  AliTRDgeometry *geo = trd->GetGeometry();

  printf("AliTRDclusterizerV1::MakeCluster -- ");
  printf("Start creating clusters.\n");

  AliTRDdataArrayI *digits;

  // Parameters
  Float_t maxThresh        = fClusMaxThresh;   // threshold value for maximum
  Float_t signalThresh     = fClusSigThresh;   // threshold value for digit signal
  Int_t   clusteringMethod = fClusMethod;      // clustering method option (for testing)

  // Iteration limit for unfolding procedure
  const Float_t kEpsilon = 0.01;             

  const Int_t   kNclus   = 3;  
  const Int_t   kNsig    = 5;

  Int_t chamBeg = 0;
  Int_t chamEnd = kNcham;
  if (trd->GetSensChamber()  >= 0) {
    chamBeg = trd->GetSensChamber();
    chamEnd = chamBeg + 1;
  }
  Int_t planBeg = 0;
  Int_t planEnd = kNplan;
  if (trd->GetSensPlane()    >= 0) {
    planBeg = trd->GetSensPlane();
    planEnd = planBeg + 1;
  }
  Int_t sectBeg = 0;
  Int_t sectEnd = kNsect;

  // *** Start clustering *** in every chamber
  for (Int_t icham = chamBeg; icham < chamEnd; icham++) {
    for (Int_t iplan = planBeg; iplan < planEnd; iplan++) {
      for (Int_t isect = sectBeg; isect < sectEnd; isect++) {

        if (trd->GetSensSector() >= 0) {
          Int_t sens1 = trd->GetSensSector();
          Int_t sens2 = sens1 + trd->GetSensSectorRange();
          sens2 -= ((Int_t) (sens2 / kNsect)) * kNsect;
          if (sens1 < sens2) {
            if ((isect < sens1) || (isect >= sens2)) continue;
	  }
          else {
            if ((isect < sens1) && (isect >= sens2)) continue;
	  }
	}

        Int_t idet = geo->GetDetector(iplan,icham,isect);

        Int_t nClusters = 0;
        printf("AliTRDclusterizerV1::MakeCluster -- ");
        printf("Analyzing chamber %d, plane %d, sector %d.\n"
               ,icham,iplan,isect);

        Int_t   nRowMax  = geo->GetRowMax(iplan,icham,isect);
        Int_t   nColMax  = geo->GetColMax(iplan);
        Int_t   nTimeMax = geo->GetTimeMax();

        // Create a detector matrix to keep maxima
        AliTRDmatrix *digitMatrix  = new AliTRDmatrix(nRowMax,nColMax,nTimeMax
                                                     ,isect,icham,iplan);
        // Create a matrix to contain maximum flags
        AliTRDmatrix *maximaMatrix = new AliTRDmatrix(nRowMax,nColMax,nTimeMax
                                                     ,isect,icham,iplan);

        // Read in the digits
        digits = fDigitsManager->GetDigits(idet);

        // Loop through the detector pixel
        for (time = 0; time < nTimeMax; time++) {
          for ( col = 0;  col <  nColMax;  col++) {
            for ( row = 0;  row <  nRowMax;  row++) {

              Int_t signal = digits->GetData(row,col,time);
              Int_t index  = digits->GetIndex(row,col,time);

              // Fill the detector matrix
              if (signal > signalThresh) {
	        // Store the signal amplitude
                digitMatrix->SetSignal(row,col,time,signal);
	        // Store the digits number
                digitMatrix->AddTrack(row,col,time,index);
              }

	    }
	  }
	}

        // Loop chamber and find maxima in digitMatrix
        for ( row = 0;  row <  nRowMax;  row++) {
          for ( col = 1;  col <  nColMax;  col++) {
            for (time = 0; time < nTimeMax; time++) {

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
        for ( row = 0;  row <  nRowMax;  row++) {
          for ( col = 1;  col <  nColMax;  col++) {
            for (time = 0; time < nTimeMax; time++) {

              if ((maximaMatrix->GetSignal(row,col,time) > 0)
                  && (digitMatrix->GetSignal(row,col,time) > maxThresh)) {

                // Ratio resulting from unfolding
                Float_t ratio                 =  0;    
                // Signals on max and neighbouring pads
                Float_t padSignal[kNsig]      = {0};   
                // Signals from cluster
                Float_t clusterSignal[kNclus] = {0};
                // Cluster pad info
                Float_t clusterPads[kNclus]   = {0};   
                // Cluster digit info
                Int_t   clusterDigit[kNclus]  = {0};

                Int_t iPad;
                for (iPad = 0; iPad < kNclus; iPad++) {
                  clusterSignal[iPad] = digitMatrix->GetSignal(row,col-1+iPad,time);
                  clusterDigit[iPad]  = digitMatrix->GetTrack(row,col-1+iPad,time,0);
                }

                // neighbouring maximum on right side?
                if (col < nColMax - 2) {
                  if (maximaMatrix->GetSignal(row,col + 2,time) > 0) {

                    for (iPad = 0; iPad < 5; iPad++) {
                      padSignal[iPad] = digitMatrix->GetSignal(row,col-1+iPad,time);
                    }

                    // unfold:
                    ratio = Unfold(kEpsilon, padSignal);

                    // set signal on overlapping pad to ratio
                    clusterSignal[2] *= ratio;

                  }
                }
                
		// Calculate the position of the cluster
                switch (clusteringMethod) {
                case 1:
                  // method 1: simply center of mass
                  clusterPads[0] = row + 0.5;
                  clusterPads[1] = col - 0.5 + (clusterSignal[2] - clusterSignal[0]) /
                                   (clusterSignal[0] + clusterSignal[1] + clusterSignal[2]);
                  clusterPads[2] = time + 0.5;

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

                  clusterPads[0] = row + 0.5;
                  clusterPads[1] = col - 1.5 + colMean;
                  clusterPads[2] = time + 0.5;

                  delete hPadCharges;

                  nClusters++;
                  break;
                }

                Float_t clusterCharge =   clusterSignal[0]
                                        + clusterSignal[1]
                                        + clusterSignal[2];

                // Add the cluster to the output array 
                trd->AddRecPoint(clusterPads,clusterDigit,idet,clusterCharge);

              }
            }  // time
          }    // col
        }      // row

        printf("AliTRDclusterizerV1::MakeCluster -- ");
        printf("Number of clusters found: %d\n",nClusters);

        delete digitMatrix;
        delete maximaMatrix;

      }          // isect
    }            // iplan
  }              // icham

  printf("AliTRDclusterizerV1::MakeCluster -- ");
  printf("Total number of points found: %d\n"
        ,trd->RecPoints()->GetEntries());

  // Get the pointer to the cluster branch
  TTree *clusterTree = gAlice->TreeR(); 

  // Fill the cluster-branch
  printf("AliTRDclusterizerV1::MakeCluster -- ");
  printf("Fill the cluster tree.\n");
  clusterTree->Fill();
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

//_____________________________________________________________________________
Float_t AliTRDclusterizerV1::PadResponse(Float_t x)
{
  //
  // The pad response for the chevron pads. 
  // We use a simple Gaussian approximation which should be good
  // enough for our purpose.
  //

  // The parameters for the response function
  const Float_t kA  =  0.8872;
  const Float_t kB  = -0.00573;
  const Float_t kC  =  0.454;
  const Float_t kC2 =  kC*kC;

  Float_t pr = kA * (kB + TMath::Exp(-x*x / (2. * kC2)));

  return (pr);

}
