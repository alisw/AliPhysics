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
Revision 1.1  2000/02/28 19:00:13  cblume
Add new TRD classes

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Creates and handles digits from TRD hits                                 //
//                                                                           //
//  The following effects are included:                                      //
//      - Diffusion                                                          //
//      - ExB effects                                                        //
//      - Gas gain including fluctuations                                    //
//      - Pad-response (simple Gaussian approximation)                       //
//      - Electronics noise                                                  //
//      - Electronics gain                                                   //
//      - Digitization                                                       //
//      - ADC threshold                                                      //
//  The corresponding parameter can be adjusted via the various              //
//  Set-functions. If these parameters are not explicitly set, default       //
//  values are used (see Init-function).                                     //
//  To produce digits from a root-file with TRD-hits use the                 //
//  slowDigitsCreate.C macro.                                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TVector.h>
#include <TRandom.h>

#include "AliTRD.h"
#include "AliTRDdigitizer.h"
#include "AliTRDmatrix.h"

ClassImp(AliTRDdigitizer)

//_____________________________________________________________________________
AliTRDdigitizer::AliTRDdigitizer():TNamed()
{
  //
  // AliTRDdigitizer default constructor
  //

  fInputFile = NULL;
  fEvent     = 0;

}

//_____________________________________________________________________________
AliTRDdigitizer::AliTRDdigitizer(const Text_t *name, const Text_t *title)
                :TNamed(name,title)
{
  //
  // AliTRDdigitizer default constructor
  //

  fInputFile   = NULL;
  fEvent       = 0;

  fDigitsArray   = new AliTRDsegmentArray(kNsect*kNplan*kNcham);
  for (Int_t iDict = 0; iDict < kNDict; iDict++) {
    fDictionary[iDict] = new AliTRDsegmentArray(kNsect*kNplan*kNcham);
  }

  Init();

}

//_____________________________________________________________________________
AliTRDdigitizer::~AliTRDdigitizer()
{

  if (fInputFile) {
    fInputFile->Close();
    delete fInputFile;
  }

  if (fDigitsArray) {
    fDigitsArray->Delete();
    delete fDigitsArray;
  }

  for (Int_t iDict = 0; iDict < kNDict; iDict++) {
    fDictionary[iDict]->Delete();
    delete fDictionary[iDict];
  }

}

//_____________________________________________________________________________
Int_t AliTRDdigitizer::Diffusion(Float_t driftlength, Float_t *xyz)
{
  //
  // Applies the diffusion smearing to the position of a single electron
  //

  Float_t driftSqrt = TMath::Sqrt(driftlength);
  Float_t sigmaT = driftSqrt * fDiffusionT;
  Float_t sigmaL = driftSqrt * fDiffusionL;
  xyz[0] = gRandom->Gaus(xyz[0], sigmaL * fLorentzFactor);
  xyz[1] = gRandom->Gaus(xyz[1], sigmaT * fLorentzFactor);
  xyz[2] = gRandom->Gaus(xyz[2], sigmaT);
  return 1;

}

//_____________________________________________________________________________
Int_t AliTRDdigitizer::ExB(Float_t driftlength, Float_t *xyz)
{
  //
  // Applies E x B effects to the position of a single electron
  //

  xyz[0] = xyz[0];
  xyz[1] = xyz[1] + fLorentzAngle * driftlength;
  xyz[2] = xyz[2];

  return 1;

}

//_____________________________________________________________________________
void AliTRDdigitizer::Init()
{
  //
  // Initializes the digitization procedure with standard values
  //

  // The default parameter for the digitization
  fGasGain       = 2.0E3;
  fNoise         = 3000.;
  fChipGain      = 10.;
  fADCoutRange   = 255.;
  fADCinRange    = 2000.;
  fADCthreshold  = 1;

  // Transverse and longitudinal diffusion coefficients (Xe/Isobutane)
  fDiffusionOn   = 1;
  fDiffusionT    = 0.060;
  fDiffusionL    = 0.017;

  // Propability for electron attachment
  fElAttachOn    = 0;
  fElAttachProp  = 0.0;

  // E x B effects
  fExBOn         = 0;
  // omega * tau. (tau ~ 12 * 10^-12, B = 0.2T)
  fLorentzAngle  = 17.6 * 12.0 * 0.2 * 0.01;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::Open(const Char_t *name, Int_t nEvent)
{
  //
  // Opens a ROOT-file with TRD-hits and reads in the hit-tree
  //

  // Connect the AliRoot file containing Geometry, Kine, and Hits
  fInputFile = (TFile*) gROOT->GetListOfFiles()->FindObject(name);
  if (!fInputFile) {
    printf("AliTRDdigitizer::Open -- ");
    printf("Open the ALIROOT-file %s.\n",name);
    fInputFile = new TFile(name,"UPDATE");
  }
  else {
    printf("AliTRDdigitizer::Open -- ");
    printf("%s is already open.\n",name);
  }

  // Get AliRun object from file or create it if not on file
  //if (!gAlice) {
    gAlice = (AliRun*) fInputFile->Get("gAlice");
    if (gAlice) {
      printf("AliTRDdigitizer::Open -- ");
      printf("AliRun object found on file.\n");
    }
    else {
      printf("AliTRDdigitizer::Open -- ");
      printf("Could not find AliRun object.\n");
      return kFALSE;
    }
  //}

  fEvent = nEvent;

  // Import the Trees for the event nEvent in the file
  Int_t nparticles = gAlice->GetEvent(fEvent);
  if (nparticles <= 0) {
    printf("AliTRDdigitizer::Open -- ");
    printf("No entries in the trees for event %d.\n",fEvent);
    return kFALSE;
  }

  return kTRUE;

}

//_____________________________________________________________________________
Float_t AliTRDdigitizer::PadResponse(Float_t x)
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
Bool_t AliTRDdigitizer::MakeDigits()
{
  //
  // Loops through the TRD-hits and creates the digits.
  //

  // Get the pointer to the detector class and check for version 1
  AliTRD *TRD = (AliTRD*) gAlice->GetDetector("TRD");
  if (TRD->IsVersion() != 1) {
    printf("AliTRDdigitizer::MakeDigits -- ");
    printf("TRD must be version 1 (slow simulator).\n");
    return kFALSE; 
  }

  // Get the geometry
  AliTRDgeometry *Geo = TRD->GetGeometry();
  printf("AliTRDdigitizer::MakeDigits -- ");
  printf("Geometry version %d\n",Geo->IsVersion());

  printf("AliTRDdigitizer::MakeDigits -- ");
  printf("Start creating digits.\n");

  ///////////////////////////////////////////////////////////////
  // Parameter 
  ///////////////////////////////////////////////////////////////

  // Converts number of electrons to fC
  const Float_t el2fC  = 1.602E-19 * 1.0E15; 

  ///////////////////////////////////////////////////////////////

  Int_t   iRow, iCol, iTime;
  Int_t   nBytes = 0;

  Int_t   totalSizeDigits = 0;
  Int_t   totalSizeDict0  = 0;
  Int_t   totalSizeDict1  = 0;
  Int_t   totalSizeDict2  = 0;

  AliTRDhit       *Hit;
  AliTRDdataArray *Digits;
  AliTRDdataArray *Dictionary[kNDict];

  // Get the pointer to the hit tree
  TTree *HitTree    = gAlice->TreeH();

  // The Lorentz factor
  if (fExBOn) {
    fLorentzFactor = 1.0 / (1.0 + fLorentzAngle*fLorentzAngle);
  }
  else {
    fLorentzFactor = 1.0;
  }

  // Get the number of entries in the hit tree
  // (Number of primary particles creating a hit somewhere)
  Int_t nTrack = (Int_t) HitTree->GetEntries();

  Int_t chamBeg = 0;
  Int_t chamEnd = kNcham;
  if (TRD->GetSensChamber() >= 0) {
    chamBeg = TRD->GetSensChamber();
    chamEnd = chamEnd + 1;
  }
  Int_t planBeg = 0;
  Int_t planEnd = kNplan;
  if (TRD->GetSensPlane()   >= 0) {
    planBeg = TRD->GetSensPlane();
    planEnd = planBeg + 1;
  }
  Int_t sectBeg = 0;
  Int_t sectEnd = kNsect;
  if (TRD->GetSensSector()  >= 0) {
    sectBeg = TRD->GetSensSector();
    sectEnd = sectBeg + 1;
  }

  // Loop through all the chambers
  for (Int_t iCham = chamBeg; iCham < chamEnd; iCham++) {
    for (Int_t iPlan = planBeg; iPlan < planEnd; iPlan++) {
      for (Int_t iSect = sectBeg; iSect < sectEnd; iSect++) {

        Int_t nDigits = 0;

        Int_t iDet = Geo->GetDetector(iPlan,iCham,iSect);

        printf("AliTRDdigitizer::MakeDigits -- ");
        printf("Digitizing chamber %d, plane %d, sector %d.\n"
              ,iCham,iPlan,iSect);

        Int_t   nRowMax     = Geo->GetRowMax(iPlan,iCham,iSect);
        Int_t   nColMax     = Geo->GetColMax(iPlan);
        Int_t   nTimeMax    = Geo->GetTimeMax();
        Float_t row0        = Geo->GetRow0(iPlan,iCham,iSect);
        Float_t col0        = Geo->GetCol0(iPlan);
        Float_t time0       = Geo->GetTime0(iPlan);
        Float_t rowPadSize  = Geo->GetRowPadSize();
        Float_t colPadSize  = Geo->GetColPadSize();
        Float_t timeBinSize = Geo->GetTimeBinSize();

        // Create a detector matrix to keep the signal and track numbers
        AliTRDmatrix *Matrix = new AliTRDmatrix(nRowMax,nColMax,nTimeMax
                                               ,iSect,iCham,iPlan);

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
            Float_t q        = Hit->fQ;
            Int_t   track    = Hit->fTrack;
            Int_t   detector = Hit->fDetector;
            Int_t   plane    = Geo->GetPlane(detector);
            Int_t   sector   = Geo->GetSector(detector);
            Int_t   chamber  = Geo->GetChamber(detector);

            if ((sector  != iSect) ||
                (plane   != iPlan) ||
                (chamber != iCham)) 
              continue;

            // Rotate the sectors on top of each other
            Float_t rot[3];
            Geo->Rotate(detector,pos,rot);

            // The hit position in pad coordinates (center pad)
            // The pad row (z-direction)
            Int_t  rowH = (Int_t) ((rot[2] -  row0) /  rowPadSize);
            // The pad column (rphi-direction)  
            Int_t  colH = (Int_t) ((rot[1] -  col0) /  colPadSize);
            // The time bucket
            Int_t timeH = (Int_t) ((rot[0] - time0) / timeBinSize);

            // Array to sum up the signal in a box surrounding the
            // hit postition
            const Int_t timeBox = 7;
            const Int_t  colBox = 9;
            const Int_t  rowBox = 7;
            Float_t signalSum[rowBox][colBox][timeBox];
            for (iRow  = 0;  iRow <  rowBox; iRow++ ) {
              for (iCol  = 0;  iCol <  colBox; iCol++ ) {
                for (iTime = 0; iTime < timeBox; iTime++) {
                  signalSum[iRow][iCol][iTime] = 0;
		}
	      }
	    }

            // Loop over all electrons of this hit
            Int_t nEl = (Int_t) q;
            for (Int_t iEl = 0; iEl < nEl; iEl++) {

	      // The driftlength
              Float_t driftlength = rot[0] - time0;
              if ((driftlength <        0) || 
                  (driftlength > kDrThick)) break;
              Float_t driftlengthL = driftlength;
              if (fExBOn) driftlengthL /= TMath::Sqrt(fLorentzFactor);
              Float_t xyz[3];
              xyz[0] = rot[0];
              xyz[1] = rot[1];
              xyz[2] = rot[2];

              // Electron attachment
              if (fElAttachOn) {
                if (gRandom->Rndm() < (driftlengthL * fElAttachProp / 100.)) continue;
	      }

              // Apply the diffusion smearing
              if (fDiffusionOn) {
                if (!(Diffusion(driftlengthL,xyz))) continue;
	      }

              // Apply E x B effects
              if (fExBOn) { 
                if (!(ExB(driftlength,xyz))) continue;   
	      }

              // The electron position and the distance to the hit position
	      // in pad units
              // The pad row (z-direction)
              Int_t  rowE = (Int_t) ((xyz[2] -  row0) /  rowPadSize);
              Int_t  rowD =  rowH -  rowE;
              // The pad column (rphi-direction)
              Int_t  colE = (Int_t) ((xyz[1] -  col0) /  colPadSize);
              Int_t  colD =  colH -  colE;
              // The time bucket
              Int_t timeE = (Int_t) ((xyz[0] - time0) / timeBinSize);
              Int_t timeD = timeH - timeE;

              // Apply the gas gain including fluctuations
              Int_t signal = (Int_t) (-fGasGain * TMath::Log(gRandom->Rndm()));

	      // The distance of the electron to the center of the pad 
	      // in units of pad width
              Float_t dist = (xyz[1] - col0 - (colE + 0.5) * colPadSize) 
                           / colPadSize;

              // Sum up the signal in the different pixels
              // and apply the pad response
              Int_t  rowIdx =  rowD + (Int_t) ( rowBox / 2);
              Int_t  colIdx =  colD + (Int_t) ( colBox / 2);
              Int_t timeIdx = timeD + (Int_t) (timeBox / 2);

              if (( rowIdx < 0) || ( rowIdx >  rowBox)) {
                printf("AliTRDdigitizer::MakeDigits -- ");
                printf("Boundary error. rowIdx = %d (%d)\n", rowIdx, rowBox);
                continue;
	      }
              if (( colIdx < 0) || ( colIdx >  colBox)) {
                printf("AliTRDdigitizer::MakeDigits -- ");
                printf("Boundary error. colIdx = %d (%d)\n", colIdx, colBox);
                continue;
	      }
              if ((timeIdx < 0) || (timeIdx > timeBox)) {
                printf("AliTRDdigitizer::MakeDigits -- ");
                printf("Boundary error. timeIdx = %d (%d)\n",timeIdx,timeBox);
                continue;
	      }
              signalSum[rowIdx][colIdx-1][timeIdx] += PadResponse(dist-1.) * signal;
              signalSum[rowIdx][colIdx  ][timeIdx] += PadResponse(dist   ) * signal;
              signalSum[rowIdx][colIdx+1][timeIdx] += PadResponse(dist+1.) * signal;

            }

            // Add the padcluster to the detector matrix
            for (iRow  = 0;  iRow <  rowBox; iRow++ ) {
              for (iCol  = 0;  iCol <  colBox; iCol++ ) {
                for (iTime = 0; iTime < timeBox; iTime++) {

                  Int_t  rowB =  rowH + iRow  - (Int_t) ( rowBox / 2); 
                  Int_t  colB =  colH + iCol  - (Int_t) ( colBox / 2);
                  Int_t timeB = timeH + iTime - (Int_t) (timeBox / 2);

                  Float_t signalB = signalSum[iRow][iCol][iTime];
                  if (signalB > 0.0) {
                    Matrix->AddSignal(rowB,colB,timeB,signalB);
                    if (!(Matrix->AddTrack(rowB,colB,timeB,track))) { 
                      printf("AliTRDdigitizer::MakeDigits -- ");
                      printf("More than three tracks in a pixel!\n");
	            }
		  }

		}
	      }
	    }

          }

 	}

        // Add a container for the digits of this detector
        Digits = (AliTRDdataArray *) fDigitsArray->At(iDet);        
        // Allocate memory space for the digits buffer
        Digits->Allocate(nRowMax,nColMax,nTimeMax);

        for (Int_t iDict = 0; iDict < kNDict; iDict++) {
          Dictionary[iDict] = (AliTRDdataArray *) fDictionary[iDict]->At(iDet);
          Dictionary[iDict]->Allocate(nRowMax,nColMax,nTimeMax);
	}

        // Create the hits for this chamber
        for (iRow  = 0; iRow  <  nRowMax; iRow++ ) {
          for (iCol  = 0; iCol  <  nColMax; iCol++ ) {
            for (iTime = 0; iTime < nTimeMax; iTime++) {         

              Float_t signalAmp = Matrix->GetSignal(iRow,iCol,iTime);

              // Add the noise
              signalAmp  = TMath::Max((Float_t) gRandom->Gaus(signalAmp,fNoise)
                                     ,(Float_t) 0.0);
	      // Convert to fC
              signalAmp *= el2fC;
              // Convert to mV
              signalAmp *= fChipGain;
	      // Convert to ADC counts
              Int_t adc  = (Int_t) (signalAmp * (fADCoutRange / fADCinRange));

              // Store the amplitude of the digit
              Digits->SetData(iRow,iCol,iTime,adc);

              // Store the track index in the dictionary
              // Note: We store index+1 in order to allow the array to be compressed
              for (Int_t iDict = 0; iDict < kNDict; iDict++) {
                Dictionary[iDict]->SetData(iRow,iCol,iTime
                                          ,Matrix->GetTrack(iRow,iCol,iTime,iDict)+1);
	      }

              if (adc > fADCthreshold) nDigits++;

	    }
	  }
	}

        // Compress the arrays
        Digits->Compress(1,fADCthreshold);
        for (Int_t iDict = 0; iDict < kNDict; iDict++) {
          Dictionary[iDict]->Compress(1,0);
	}

        totalSizeDigits += Digits->GetSize();
        totalSizeDict0  += Dictionary[0]->GetSize();
        totalSizeDict1  += Dictionary[1]->GetSize();
        totalSizeDict2  += Dictionary[2]->GetSize();

        printf("AliTRDdigitizer::MakeDigits -- ");
        printf("Number of digits found: %d.\n",nDigits);

	// Clean up
        if (Matrix) delete Matrix;

      }
    }
  }

  printf("AliTRDdigitizer::MakeDigits -- ");
  printf("Total digits data size = %d, %d, %d, %d\n",totalSizeDigits
                                                    ,totalSizeDict0
                                                    ,totalSizeDict1
                                                    ,totalSizeDict2);        

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::MakeBranch()
{
  //
  // Creates the branches for the digits and the dictionary
  //

  Int_t buffersize = 64000;

  Bool_t status = kTRUE;

  if (gAlice->TreeD()) {

    // Make the branch for the digits
    if (fDigitsArray) {
      const AliTRDdataArray *Digits = 
           (AliTRDdataArray *) fDigitsArray->At(0);
      if (Digits) {
        gAlice->TreeD()->Branch("TRDdigits",Digits->IsA()->GetName()
                                           ,&Digits,buffersize,1);
        printf("AliTRDdigitizer::MakeBranch -- ");
        printf("Making branch TRDdigits\n");
      }
      else {
        status = kFALSE;
      }
    }
    else {
      status = kFALSE;
    }

    // Make the branches for the dictionaries
    for (Int_t iDict = 0; iDict < kNDict; iDict++) {

      Char_t branchname[15];
      sprintf(branchname,"TRDdictionary%d",iDict);
      if (fDictionary[iDict]) {
        const AliTRDdataArray *Dictionary = 
             (AliTRDdataArray *) fDictionary[iDict]->At(0);
        if (Dictionary) {
          gAlice->TreeD()->Branch(branchname,Dictionary->IsA()->GetName()
                                            ,&Dictionary,buffersize,1);
          printf("AliTRDdigitizer::MakeBranch -- ");
          printf("Making branch %s\n",branchname);
	}
        else {
          status = kFALSE;
	}
      }
      else {
        status = kFALSE;
      }
    }

  }
  else {
    status = kFALSE;
  }

  return status;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::WriteDigits()
{
  //
  // Writes out the TRD-digits and the dictionaries
  //

  // Create the branches
  if (!(gAlice->TreeD()->GetBranch("TRDdigits"))) { 
    if (!MakeBranch()) return kFALSE;
  }

  // Store the contents of the segment array in the tree
  if (!fDigitsArray->StoreArray("TRDdigits")) {
    printf("AliTRDdigitizer::WriteDigits -- ");
    printf("Error while storing digits in branch TRDdigits\n");
    return kFALSE;
  }
  for (Int_t iDict = 0; iDict < kNDict; iDict++) {
    Char_t branchname[15];
    sprintf(branchname,"TRDdictionary%d",iDict);
    if (!fDictionary[iDict]->StoreArray(branchname)) {
      printf("AliTRDdigitizer::WriteDigits -- ");
      printf("Error while storing dictionary in branch %s\n",branchname);
      return kFALSE;
    }
  }

  // Write the new tree into the input file (use overwrite option)
  Char_t treeName[7];
  sprintf(treeName,"TreeD%d",fEvent);
  printf("AliTRDdigitizer::WriteDigits -- ");
  printf("Write the digits tree %s for event %d.\n"
        ,treeName,fEvent);
  gAlice->TreeD()->Write(treeName,2);
 
  return kTRUE;

}

ClassImp(AliTRDdigit)

//_____________________________________________________________________________
AliTRDdigit::AliTRDdigit(Int_t *digits):AliDigitNew()
{
  //
  // Create a TRD digit
  //

  // Store the volume hierarchy
  fDetector  = digits[0];

  // Store the row, pad, and time bucket number
  fRow       = digits[1];
  fCol       = digits[2];
  fTime      = digits[3];

  // Store the signal amplitude
  fAmplitude = digits[4];

}
