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
Revision 1.31  2002/02/11 14:27:11  cblume
New pad plane design, new TRF+PRF, tail cancelation, cross talk

Revision 1.30  2001/11/19 08:44:08  cblume
Fix bugs reported by Rene

Revision 1.29  2001/11/14 19:44:25  hristov
Numeric const casted (Alpha)

Revision 1.28  2001/11/14 16:35:58  cblume
Inherits now from AliDetector

Revision 1.27  2001/11/14 10:50:45  cblume
Changes in digits IO. Add merging of summable digits

Revision 1.26  2001/11/06 17:19:41  cblume
Add detailed geometry and simple simulator

Revision 1.25  2001/06/27 09:54:44  cblume
Moved fField initialization to InitDetector()

Revision 1.24  2001/05/21 16:45:47  hristov
Last minute changes (C.Blume)

Revision 1.23  2001/05/07 08:04:48  cblume
New TRF and PRF. Speedup of the code. Digits from amplification region included

Revision 1.22  2001/03/30 14:40:14  cblume
Update of the digitization parameter

Revision 1.21  2001/03/13 09:30:35  cblume
Update of digitization. Moved digit branch definition to AliTRD

Revision 1.20  2001/02/25 20:19:00  hristov
Minor correction: loop variable declared only once for HP, Sun

Revision 1.19  2001/02/14 18:22:26  cblume
Change in the geometry of the padplane

Revision 1.18  2001/01/26 19:56:57  hristov
Major upgrade of AliRoot code

Revision 1.17  2000/12/08 12:53:27  cblume
Change in Copy() function for HP-compiler

Revision 1.16  2000/12/07 12:20:46  cblume
Go back to array compression. Use sampled PRF to speed up digitization

Revision 1.15  2000/11/23 14:34:08  cblume
Fixed bug in expansion routine of arrays (initialize buffers properly)

Revision 1.14  2000/11/20 08:54:44  cblume
Switch off compression as default

Revision 1.13  2000/11/10 14:57:52  cblume
Changes in the geometry constants for the DEC compiler

Revision 1.12  2000/11/01 14:53:20  cblume
Merge with TRD-develop

Revision 1.1.4.9  2000/10/26 17:00:22  cblume
Fixed bug in CheckDetector()

Revision 1.1.4.8  2000/10/23 13:41:35  cblume
Added protection against Log(0) in the gas gain calulation

Revision 1.1.4.7  2000/10/17 02:27:34  cblume
Get rid of global constants

Revision 1.1.4.6  2000/10/16 01:16:53  cblume
Changed timebin 0 to be the one closest to the readout

Revision 1.1.4.5  2000/10/15 23:34:29  cblume
Faster version of the digitizer

Revision 1.1.4.4  2000/10/06 16:49:46  cblume
Made Getters const

Revision 1.1.4.3  2000/10/04 16:34:58  cblume
Replace include files by forward declarations

Revision 1.1.4.2  2000/09/22 14:41:10  cblume
Bug fix in PRF. Included time response. New structure

Revision 1.10  2000/10/05 07:27:53  cblume
Changes in the header-files by FCA

Revision 1.9  2000/10/02 21:28:19  fca
Removal of useless dependecies via forward declarations

Revision 1.8  2000/06/09 11:10:07  cblume
Compiler warnings and coding conventions, next round

Revision 1.7  2000/06/08 18:32:58  cblume
Make code compliant to coding conventions

Revision 1.6  2000/06/07 16:27:32  cblume
Try to remove compiler warnings on Sun and HP

Revision 1.5  2000/05/09 16:38:57  cblume
Removed PadResponse(). Merge problem

Revision 1.4  2000/05/08 15:53:45  cblume
Resolved merge conflict

Revision 1.3  2000/04/28 14:49:27  cblume
Only one declaration of iDict in MakeDigits()

Revision 1.1.4.1  2000/05/08 14:42:04  cblume
Introduced AliTRDdigitsManager

Revision 1.1  2000/02/28 19:00:13  cblume
Add new TRD classes

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Creates and handles digits from TRD hits                                 //
//  Author: C. Blume (C.Blume@gsi.de)                                        //
//                                                                           //
//  The following effects are included:                                      //
//      - Diffusion                                                          //
//      - ExB effects                                                        //
//      - Gas gain including fluctuations                                    //
//      - Pad-response (simple Gaussian approximation)                       //
//      - Time-response                                                      //
//      - Electronics noise                                                  //
//      - Electronics gain                                                   //
//      - Digitization                                                       //
//      - ADC threshold                                                      //
//  The corresponding parameter can be adjusted via the various              //
//  Set-functions. If these parameters are not explicitly set, default       //
//  values are used (see Init-function).                                     //
//  As an example on how to use this class to produce digits from hits       //
//  have a look at the macro hits2digits.C                                   //
//  The production of summable digits is demonstrated in hits2sdigits.C      //
//  and the subsequent conversion of the s-digits into normal digits is      //
//  explained in sdigits2digits.C.                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>

#include <TMath.h>
#include <TVector.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TList.h>
#include <TTask.h>

#include "AliRun.h"
#include "AliMagF.h"
#include "AliRunDigitizer.h"

#include "AliTRD.h"
#include "AliTRDhit.h"
#include "AliTRDdigitizer.h"
#include "AliTRDdataArrayI.h"
#include "AliTRDdataArrayF.h"
#include "AliTRDsegmentArray.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDgeometry.h"

ClassImp(AliTRDdigitizer)

//_____________________________________________________________________________
AliTRDdigitizer::AliTRDdigitizer()
{
  //
  // AliTRDdigitizer default constructor
  //

  fInputFile          = NULL;
  fDigitsManager      = NULL;
  fSDigitsManagerList = NULL;
  fSDigitsManager     = NULL;
  fTRD                = NULL;
  fGeo                = NULL;
  fPRFsmp             = NULL;
  fTRFsmp             = NULL;
  fCTsmp              = NULL;

  fMasks              = 0;

  fEvent              = 0;
  fGasGain            = 0.0;
  fNoise              = 0.0;
  fChipGain           = 0.0;
  fADCoutRange        = 0.0;
  fADCinRange         = 0.0;
  fADCthreshold       = 0;
  fDiffusionOn        = 0;
  fDiffusionT         = 0.0;
  fDiffusionL         = 0.0;
  fElAttachOn         = 0;
  fElAttachProp       = 0.0;
  fExBOn              = 0;
  fOmegaTau           = 0.0;
  fPRFOn              = 0;
  fTRFOn              = 0;
  fCTOn               = 0;
  fTCOn               = 0;
  fDriftVelocity      = 0.0;
  fPadCoupling        = 0.0;
  fTimeCoupling       = 0.0;
  fTimeBinWidth       = 0.0;
  fField              = 0.0;
  fTiltingAngle       = 0.0;

  fPRFbin             = 0;
  fPRFlo              = 0.0;
  fPRFhi              = 0.0;
  fPRFwid             = 0.0;
  fPRFpad             = 0;
  fTRFbin             = 0;
  fTRFlo              = 0.0;
  fTRFhi              = 0.0;
  fTRFwid             = 0.0;
  fTCnexp             = 0;

  fCompress           = kTRUE;
  fDebug              = 0;
  fSDigits            = kFALSE;
  fSDigitsScale       = 0.0;

}

//_____________________________________________________________________________
AliTRDdigitizer::AliTRDdigitizer(const Text_t *name, const Text_t *title)
                :AliDigitizer(name,title)
{
  //
  // AliTRDdigitizer constructor
  //

  fInputFile          = NULL;

  fDigitsManager      = NULL;
  fSDigitsManager     = NULL;
  fSDigitsManagerList = NULL;

  fTRD                = NULL;
  fGeo                = NULL;
  fPRFsmp             = NULL;
  fTRFsmp             = NULL;
  fCTsmp              = NULL;

  fMasks              = 0;

  fEvent              = 0;

  fCompress           = kTRUE;
  fDebug              = 0;
  fSDigits            = kFALSE;

  Init();

}

//_____________________________________________________________________________
AliTRDdigitizer::AliTRDdigitizer(AliRunDigitizer *manager
                                , const Text_t *name, const Text_t *title)
                :AliDigitizer(manager,name,title)
{
  //
  // AliTRDdigitizer constructor
  //

  fInputFile          = NULL;

  fDigitsManager      = NULL;
  fSDigitsManager     = NULL;
  fSDigitsManagerList = NULL;

  fTRD                = NULL;
  fGeo                = NULL;
  fPRFsmp             = NULL;
  fTRFsmp             = NULL;
  fCTsmp              = NULL;

  fMasks              = 0;

  fEvent              = 0;

  fCompress           = kTRUE;
  fDebug              = 0;
  fSDigits            = kFALSE;

  Init();

}

//_____________________________________________________________________________
AliTRDdigitizer::AliTRDdigitizer(AliRunDigitizer *manager)
                :AliDigitizer(manager,"AliTRDdigitizer","TRD digitizer")
{
  //
  // AliTRDdigitizer constructor
  //

  fInputFile          = NULL;

  fDigitsManager      = NULL;
  fSDigitsManager     = NULL;
  fSDigitsManagerList = NULL;

  fTRD                = NULL;
  fGeo                = NULL;
  fPRFsmp             = NULL;
  fTRFsmp             = NULL;
  fCTsmp              = NULL;

  fMasks              = 0;

  fEvent              = 0;

  fCompress           = kTRUE;
  fDebug              = 0;
  fSDigits            = kFALSE;

  Init();

}

//_____________________________________________________________________________
AliTRDdigitizer::AliTRDdigitizer(const AliTRDdigitizer &d)
{
  //
  // AliTRDdigitizer copy constructor
  //

  ((AliTRDdigitizer &) d).Copy(*this);

}

//_____________________________________________________________________________
AliTRDdigitizer::~AliTRDdigitizer()
{
  //
  // AliTRDdigitizer destructor
  //

  if (fInputFile) {
    fInputFile->Close();
    delete fInputFile;
    fInputFile = NULL;
  }

  if (fDigitsManager) {
    delete fDigitsManager;
    fDigitsManager = NULL;
  }

  if (fSDigitsManager) {
    delete fSDigitsManager;
    fSDigitsManager = NULL;
  }

  if (fSDigitsManagerList) {
    fSDigitsManagerList->Delete();
    delete fSDigitsManagerList;
    fSDigitsManagerList = NULL;
  }

  if (fTRFsmp) {
    delete [] fTRFsmp;
    fTRFsmp = NULL;
  }

  if (fPRFsmp) {
    delete [] fPRFsmp;
    fPRFsmp = NULL;
  }

  if (fCTsmp) {
    delete [] fCTsmp;
    fCTsmp  = NULL;
  }

  if (fMasks) {
    delete [] fMasks;
    fMasks = 0;
  }

}

//_____________________________________________________________________________
AliTRDdigitizer &AliTRDdigitizer::operator=(const AliTRDdigitizer &d)
{
  //
  // Assignment operator
  //

  if (this != &d) ((AliTRDdigitizer &) d).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDdigitizer::Copy(TObject &d)
{
  //
  // Copy function
  //

  Int_t iBin;

  ((AliTRDdigitizer &) d).fInputFile          = NULL;
  ((AliTRDdigitizer &) d).fSDigitsManagerList = NULL;
  ((AliTRDdigitizer &) d).fSDigitsManager     = NULL;
  ((AliTRDdigitizer &) d).fDigitsManager      = NULL;
  ((AliTRDdigitizer &) d).fTRD                = NULL;
  ((AliTRDdigitizer &) d).fGeo                = NULL;

  ((AliTRDdigitizer &) d).fMasks              = 0;

  ((AliTRDdigitizer &) d).fEvent              = 0;

  ((AliTRDdigitizer &) d).fGasGain            = fGasGain;
  ((AliTRDdigitizer &) d).fNoise              = fNoise;
  ((AliTRDdigitizer &) d).fChipGain           = fChipGain;
  ((AliTRDdigitizer &) d).fADCoutRange        = fADCoutRange;
  ((AliTRDdigitizer &) d).fADCinRange         = fADCinRange;
  ((AliTRDdigitizer &) d).fADCthreshold       = fADCthreshold;
  ((AliTRDdigitizer &) d).fDiffusionOn        = fDiffusionOn; 
  ((AliTRDdigitizer &) d).fDiffusionT         = fDiffusionT;
  ((AliTRDdigitizer &) d).fDiffusionL         = fDiffusionL;
  ((AliTRDdigitizer &) d).fElAttachOn         = fElAttachOn;
  ((AliTRDdigitizer &) d).fElAttachProp       = fElAttachProp;
  ((AliTRDdigitizer &) d).fExBOn              = fExBOn;
  ((AliTRDdigitizer &) d).fOmegaTau           = fOmegaTau;
  ((AliTRDdigitizer &) d).fLorentzFactor      = fLorentzFactor;
  ((AliTRDdigitizer &) d).fDriftVelocity      = fDriftVelocity;
  ((AliTRDdigitizer &) d).fPadCoupling        = fPadCoupling;
  ((AliTRDdigitizer &) d).fTimeCoupling       = fTimeCoupling;
  ((AliTRDdigitizer &) d).fTimeBinWidth       = fTimeBinWidth;
  ((AliTRDdigitizer &) d).fField              = fField;
  ((AliTRDdigitizer &) d).fPRFOn              = fPRFOn;
  ((AliTRDdigitizer &) d).fTRFOn              = fTRFOn;
  ((AliTRDdigitizer &) d).fCTOn               = fCTOn;
  ((AliTRDdigitizer &) d).fTCOn               = fTCOn;
  ((AliTRDdigitizer &) d).fTiltingAngle       = fTiltingAngle;

  ((AliTRDdigitizer &) d).fCompress           = fCompress;
  ((AliTRDdigitizer &) d).fDebug              = fDebug  ;
  ((AliTRDdigitizer &) d).fSDigits            = fSDigits;
  ((AliTRDdigitizer &) d).fSDigitsScale       = fSDigitsScale;

  ((AliTRDdigitizer &) d).fPRFbin             = fPRFbin;
  ((AliTRDdigitizer &) d).fPRFlo              = fPRFlo;
  ((AliTRDdigitizer &) d).fPRFhi              = fPRFhi;
  ((AliTRDdigitizer &) d).fPRFwid             = fPRFwid;
  ((AliTRDdigitizer &) d).fPRFpad             = fPRFpad;
  if (((AliTRDdigitizer &) d).fPRFsmp) delete [] ((AliTRDdigitizer &) d).fPRFsmp;
  ((AliTRDdigitizer &) d).fPRFsmp = new Float_t[fPRFbin];
  for (iBin = 0; iBin < fPRFbin; iBin++) {
    ((AliTRDdigitizer &) d).fPRFsmp[iBin] = fPRFsmp[iBin];
  }                                                                             
  ((AliTRDdigitizer &) d).fTRFbin             = fTRFbin;
  ((AliTRDdigitizer &) d).fTRFlo              = fTRFlo;
  ((AliTRDdigitizer &) d).fTRFhi              = fTRFhi;
  ((AliTRDdigitizer &) d).fTRFwid             = fTRFwid;
  if (((AliTRDdigitizer &) d).fTRFsmp) delete [] ((AliTRDdigitizer &) d).fTRFsmp;
  ((AliTRDdigitizer &) d).fTRFsmp = new Float_t[fTRFbin];
  for (iBin = 0; iBin < fTRFbin; iBin++) {
    ((AliTRDdigitizer &) d).fTRFsmp[iBin] = fTRFsmp[iBin];
  }                                      

  if (((AliTRDdigitizer &) d).fCTsmp)  delete [] ((AliTRDdigitizer &) d).fCTsmp;
  ((AliTRDdigitizer &) d).fCTsmp  = new Float_t[fTRFbin];
  for (iBin = 0; iBin < fTRFbin; iBin++) {
    ((AliTRDdigitizer &) d).fCTsmp[iBin]  = fCTsmp[iBin];
  }                                      

  ((AliTRDdigitizer &) d).fTCnexp             = fTCnexp;
                                      
}

//_____________________________________________________________________________
Float_t AliTRDdigitizer::CrossTalk(Float_t time)
{
  //
  // Applies the pad-pad capacitive cross talk
  //

  Int_t iBin = ((Int_t) ((time - fTRFlo) / fTRFwid)); 
  if ((iBin >= 0) && (iBin < fTRFbin)) {
    return fCTsmp[iBin];
  }
  else {
    return 0.0;
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
  xyz[1] = xyz[1] + fOmegaTau * driftlength;
  xyz[2] = xyz[2];

  return 1;

}

//_____________________________________________________________________________
Int_t AliTRDdigitizer::PadResponse(Float_t signal, Float_t dist
                                 , Int_t plane, Float_t *pad)
{
  //
  // Applies the pad response
  //

  const Int_t kNplan = AliTRDgeometry::kNplan;

  Int_t iBin  = ((Int_t) (( - dist - fPRFlo) / fPRFwid));
  Int_t iOff  = plane * fPRFbin;

  Int_t iBin0 = iBin - fPRFpad + iOff;
  Int_t iBin1 = iBin           + iOff;
  Int_t iBin2 = iBin + fPRFpad + iOff;

  if ((iBin0 >= 0) && (iBin2 < (fPRFbin*kNplan))) {

    pad[0] = signal * fPRFsmp[iBin0];
    pad[1] = signal * fPRFsmp[iBin1];
    pad[2] = signal * fPRFsmp[iBin2];

    return 1;

  }
  else {

    return 0;

  }

}

//_____________________________________________________________________________
Float_t AliTRDdigitizer::TimeResponse(Float_t time)
{
  //
  // Applies the preamp shaper time response
  //

  Int_t iBin = ((Int_t) ((time - fTRFlo) / fTRFwid)); 
  if ((iBin >= 0) && (iBin < fTRFbin)) {
    return fTRFsmp[iBin];
  }
  else {
    return 0.0;
  }    

}

//_____________________________________________________________________________
Float_t AliTRDdigitizer::Col0Tilted(Float_t col0, Float_t rowOffset
                                  , Int_t plane)
{
  //
  // Calculates col0 for tilted pads
  //

  Float_t diff = fTiltingAngle * rowOffset;
  return (col0 + TMath::Power(-1.0,plane) * diff);

}

//_____________________________________________________________________________
void AliTRDdigitizer::Exec(Option_t* option)
{
  //
  // Executes the merging
  //

  Int_t iInput;

  AliTRDdigitsManager *sdigitsManager;

  TString optionString = option;
  if (optionString.Contains("deb")) {
    fDebug = 1;
    if (optionString.Contains("2")) {
      fDebug = 2;
    }
    printf("<AliTRDdigitizer::Exec> ");
    printf("Called with debug option %d\n",fDebug);
  }

  Int_t nInput = fManager->GetNinputs();
  fMasks = new Int_t[nInput];
  for (iInput = 0; iInput < nInput; iInput++) {
    fMasks[iInput] = fManager->GetMask(iInput);
  }
  
  // Set the event number
  fEvent = gAlice->GetEvNumber();

  // Initialization
  InitDetector();

  for (iInput = 0; iInput < nInput; iInput++) {

    if (fDebug > 0) {
      printf("<AliTRDdigitizer::Exec> ");
      printf("Add input stream %d\n",iInput);
    }

    // Read the s-digits via digits manager
    sdigitsManager = new AliTRDdigitsManager();
    sdigitsManager->SetDebug(fDebug);
    sdigitsManager->SetSDigits(kTRUE);
    sdigitsManager->ReadDigits(fManager->GetInputTreeTRDS(iInput));

    // Add the s-digits to the input list 
    AddSDigitsManager(sdigitsManager);

  }

  // Convert the s-digits to normal digits
  if (fDebug > 0) {
    printf("<AliTRDdigitizer::Exec> ");
    printf("Do the conversion\n");
  }
  SDigits2Digits();

  // Store the digits
  if (fDebug > 0) {
    printf("<AliTRDdigitizer::Exec> ");
    printf("Write the digits\n");
  }
  fDigitsManager->MakeBranch(fManager->GetTreeDTRD());
  fDigitsManager->WriteDigits();
  if (fDebug > 0) {
    printf("<AliTRDdigitizer::Exec> ");
    printf("Done\n");
  }

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::Init()
{
  //
  // Initializes the digitization procedure with standard values
  //

  // The default parameter for the digitization
  fGasGain        = 2800.;
  fChipGain       = 6.1;
  fNoise          = 1000.;
  fADCoutRange    = 1023.;          // 10-bit ADC
  fADCinRange     = 1000.;          // 1V input range
  fADCthreshold   = 1;

  // For the summable digits
  fSDigitsScale   = 100.;

  // The drift velocity (cm / mus)
  fDriftVelocity  = 1.5;

  // Diffusion on
  fDiffusionOn    = 1;

  // E x B effects
  fExBOn          = 0;

  // Propability for electron attachment
  fElAttachOn     = 0;
  fElAttachProp   = 0.0;

  // The pad response function
  fPRFOn          = 1;

  // The time response function
  fTRFOn          = 1;

  // The cross talk
  fCTOn           = 0;

  // The tail cancelation
  fTCOn           = 1;
  
  // The number of exponentials
  fTCnexp         = 2;

  // The pad coupling factor (same number as for the TPC)
  fPadCoupling    = 0.5;

  // The time coupling factor (same number as for the TPC)
  fTimeCoupling   = 0.4;

  // The tilting angle for the readout pads
  SetTiltingAngle(5.0);

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::ReInit()
{
  //
  // Reinitializes the digitization procedure after a change in the parameter
  //

  if (!fGeo) {
    printf("AliTRDdigitizer::ReInit -- ");
    printf("No geometry defined. Run InitDetector() first\n");
    return kFALSE;
  }

  // Calculate the time bin width in ns
  fTimeBinWidth   = fGeo->GetTimeBinSize() / fDriftVelocity * 1000.0;

  // The range and the binwidth for the sampled TRF 
  fTRFbin = 100;
  // Start 0.2 mus before the signal
  fTRFlo  = -0.2 * fDriftVelocity;
  // End the maximum driftlength after the signal 
  fTRFhi  = AliTRDgeometry::DrThick() 
          + fGeo->GetTimeAfter() * fGeo->GetTimeBinSize();
  fTRFwid = (fTRFhi - fTRFlo) / ((Float_t) fTRFbin);

  // Transverse and longitudinal diffusion coefficients (Xe/CO2)
  fDiffusionT     = GetDiffusionT(fDriftVelocity,fField);
  fDiffusionL     = GetDiffusionL(fDriftVelocity,fField);

  // omega * tau.= tan(Lorentz-angle)
  fOmegaTau       = GetOmegaTau(fDriftVelocity,fField);

  // The Lorentz factor
  if (fExBOn) {
    fLorentzFactor = 1.0 / (1.0 + fOmegaTau*fOmegaTau);
  }
  else {
    fLorentzFactor = 1.0;
  }

  return kTRUE;

}

//_____________________________________________________________________________
void AliTRDdigitizer::SampleTRF()
{
  //
  // Samples the time response function
  //
  // New TRF from Venelin Angelov, simulated with CADENCE
  // Pad-ground capacitance = 25 pF
  // Pad-pad cross talk capacitance = 6 pF   
  //

  Int_t   ipos1;
  Int_t   ipos2;
  Float_t diff;

  const Int_t kNpasa     = 252;

  Float_t time[kNpasa]   = { -0.220000, -0.210000, -0.200000, -0.190000 
                           , -0.180000, -0.170000, -0.160000, -0.150000 
                           , -0.140000, -0.130000, -0.120000, -0.110000 
                           , -0.100000, -0.090000, -0.080000, -0.070000 
                           , -0.060000, -0.050000, -0.040000, -0.030000 
                           , -0.020000, -0.010000, -0.000000,  0.010000 
                           ,  0.020000,  0.030000,  0.040000,  0.050000 
                           ,  0.060000,  0.070000,  0.080000,  0.090000 
                           ,  0.100000,  0.110000,  0.120000,  0.130000 
                           ,  0.140000,  0.150000,  0.160000,  0.170000 
                           ,  0.180000,  0.190000,  0.200000,  0.210000 
                           ,  0.220000,  0.230000,  0.240000,  0.250000 
                           ,  0.260000,  0.270000,  0.280000,  0.290000 
                           ,  0.300000,  0.310000,  0.320000,  0.330000 
                           ,  0.340000,  0.350000,  0.360000,  0.370000 
                           ,  0.380000,  0.390000,  0.400000,  0.410000 
                           ,  0.420000,  0.430000,  0.440000,  0.450000 
                           ,  0.460000,  0.470000,  0.480000,  0.490000 
                           ,  0.500000,  0.510000,  0.520000,  0.530000 
                           ,  0.540000,  0.550000,  0.560000,  0.570000 
                           ,  0.580000,  0.590000,  0.600000,  0.610000 
                           ,  0.620000,  0.630000,  0.640000,  0.650000 
                           ,  0.660000,  0.670000,  0.680000,  0.690000 
                           ,  0.700000,  0.710000,  0.720000,  0.730000 
                           ,  0.740000,  0.750000,  0.760000,  0.770000 
                           ,  0.780000,  0.790000,  0.800000,  0.810000 
                           ,  0.820000,  0.830000,  0.840000,  0.850000 
                           ,  0.860000,  0.870000,  0.880000,  0.890000 
                           ,  0.900000,  0.910000,  0.920000,  0.930000 
                           ,  0.940000,  0.950000,  0.960000,  0.970000 
                           ,  0.980000,  0.990000,  1.000000,  1.010000 
                           ,  1.020000,  1.030000,  1.040000,  1.050000 
                           ,  1.060000,  1.070000,  1.080000,  1.090000 
                           ,  1.100000,  1.110000,  1.120000,  1.130000 
                           ,  1.140000,  1.150000,  1.160000,  1.170000 
                           ,  1.180000,  1.190000,  1.200000,  1.210000 
                           ,  1.220000,  1.230000,  1.240000,  1.250000 
                           ,  1.260000,  1.270000,  1.280000,  1.290000 
                           ,  1.300000,  1.310000,  1.320000,  1.330000 
                           ,  1.340000,  1.350000,  1.360000,  1.370000 
                           ,  1.380000,  1.390000,  1.400000,  1.410000 
                           ,  1.420000,  1.430000,  1.440000,  1.450000 
                           ,  1.460000,  1.470000,  1.480000,  1.490000 
                           ,  1.500000,  1.510000,  1.520000,  1.530000 
                           ,  1.540000,  1.550000,  1.560000,  1.570000 
                           ,  1.580000,  1.590000,  1.600000,  1.610000 
                           ,  1.620000,  1.630000,  1.640000,  1.650000 
                           ,  1.660000,  1.670000,  1.680000,  1.690000 
                           ,  1.700000,  1.710000,  1.720000,  1.730000 
                           ,  1.740000,  1.750000,  1.760000,  1.770000 
                           ,  1.780000,  1.790000,  1.800000,  1.810000 
                           ,  1.820000,  1.830000,  1.840000,  1.850000 
                           ,  1.860000,  1.870000,  1.880000,  1.890000 
                           ,  1.900000,  1.910000,  1.920000,  1.930000 
                           ,  1.940000,  1.950000,  1.960000,  1.970000 
                           ,  1.980000,  1.990000,  2.000000,  2.010000 
                           ,  2.020000,  2.030000,  2.040000,  2.050000 
                           ,  2.060000,  2.070000,  2.080000,  2.090000 
                           ,  2.100000,  2.110000,  2.120000,  2.130000 
                           ,  2.140000,  2.150000,  2.160000,  2.170000 
                           ,  2.180000,  2.190000,  2.200000,  2.210000 
                           ,  2.220000,  2.230000,  2.240000,  2.250000 
                           ,  2.260000,  2.270000,  2.280000,  2.290000 };

  Float_t signal[kNpasa] = {  0.000000,  0.000000,  0.000000,  0.000000 
                           ,  0.000000,  0.000000,  0.000000,  0.000396 
                           ,  0.005096,  0.022877,  0.061891,  0.126614 
                           ,  0.215798,  0.324406,  0.444507,  0.566817 
                           ,  0.683465,  0.787089,  0.873159,  0.937146 
                           ,  0.979049,  0.999434,  1.000000,  0.983579 
                           ,  0.954134,  0.913364,  0.866365,  0.813703 
                           ,  0.759910,  0.706116,  0.653454,  0.603624 
                           ,  0.556625,  0.514156,  0.475085,  0.439977 
                           ,  0.408834,  0.380578,  0.355549,  0.333352 
                           ,  0.313647,  0.296093,  0.280351,  0.266195 
                           ,  0.253397,  0.241789,  0.231257,  0.221574 
                           ,  0.212627,  0.204417,  0.196772,  0.189581 
                           ,  0.182956,  0.176784,  0.171008,  0.165515 
                           ,  0.160419,  0.155606,  0.151076,  0.146716 
                           ,  0.142639,  0.138845,  0.135221,  0.131767 
                           ,  0.128482,  0.125368,  0.122424,  0.119592 
                           ,  0.116931,  0.114326,  0.111891,  0.109513 
                           ,  0.107248,  0.105096,  0.103058,  0.101019 
                           ,  0.099151,  0.097282,  0.095527,  0.093715 
                           ,  0.092129,  0.090544,  0.088958,  0.087429 
                           ,  0.086014,  0.084598,  0.083239,  0.081880 
                           ,  0.080634,  0.079388,  0.078143,  0.077010 
                           ,  0.075878,  0.074745,  0.073669,  0.072593 
                           ,  0.071574,  0.070612,  0.069649,  0.068686 
                           ,  0.067780,  0.066874,  0.066025,  0.065176 
                           ,  0.064326,  0.063533,  0.062684,  0.061948 
                           ,  0.061212,  0.060419,  0.059740,  0.059003 
                           ,  0.058324,  0.057644,  0.057022,  0.056342 
                           ,  0.055663,  0.055096,  0.054473,  0.053851 
                           ,  0.053284,  0.052718,  0.052152,  0.051585 
                           ,  0.051019,  0.050566,  0.050000,  0.049490 
                           ,  0.048981,  0.048528,  0.048018,  0.047508 
                           ,  0.047055,  0.046602,  0.046149,  0.045696 
                           ,  0.045300,  0.044904,  0.044451,  0.044054 
                           ,  0.043658,  0.043205,  0.042865,  0.042469 
                           ,  0.042072,  0.041733,  0.041336,  0.040997 
                           ,  0.040657,  0.040260,  0.039921,  0.039581 
                           ,  0.039241,  0.038958,  0.038618,  0.038335 
                           ,  0.037995,  0.037656,  0.037373,  0.037089 
                           ,  0.036806,  0.036467,  0.036183,  0.035900 
                           ,  0.035617,  0.035334,  0.035108,  0.034824 
                           ,  0.034541,  0.034315,  0.034032,  0.033805 
                           ,  0.033522,  0.033296,  0.033069,  0.032786 
                           ,  0.032559,  0.032333,  0.032106,  0.031880 
                           ,  0.031653,  0.031427,  0.031200,  0.030974 
                           ,  0.030804,  0.030578,  0.030351,  0.030125 
                           ,  0.029955,  0.029785,  0.029558,  0.029332 
                           ,  0.029162,  0.028992,  0.028766,  0.028596 
                           ,  0.028426,  0.028199,  0.028086,  0.027860 
                           ,  0.027746,  0.027633,  0.027463,  0.027293 
                           ,  0.027180,  0.027067,  0.026954,  0.026954 
                           ,  0.026840,  0.026727,  0.026727,  0.026614 
                           ,  0.026614,  0.026614,  0.026557,  0.026501 
                           ,  0.026501,  0.026501,  0.026501,  0.026501 
                           ,  0.026501,  0.026501,  0.026501,  0.026387 
                           ,  0.026387,  0.026387,  0.026387,  0.026387 
                           ,  0.026387,  0.026387,  0.026387,  0.026387 
                           ,  0.026387,  0.026387,  0.026387,  0.026387 
                           ,  0.026387,  0.026274,  0.026274,  0.026274 
                           ,  0.026274,  0.026274,  0.026274,  0.026274 
                           ,  0.026274,  0.026274,  0.026274,  0.026274 
                           ,  0.026274,  0.026274,  0.026274,  0.026161 };

  Float_t xtalk[kNpasa]  = {  0.000000,  0.000000,  0.000000,  0.000000 
                           ,  0.000000,  0.000000,  0.000000,  0.000113 
                           ,  0.000793,  0.003058,  0.007305,  0.013194 
                           ,  0.019706,  0.025821,  0.030634,  0.033465 
                           ,  0.034145,  0.032729,  0.029615,  0.025198 
                           ,  0.019989,  0.014496,  0.009003,  0.003964 
                           , -0.000510, -0.004190, -0.007191, -0.009400 
                           , -0.010872, -0.011835, -0.012288, -0.012288 
                           , -0.012005, -0.011495, -0.010872, -0.010136 
                           , -0.009343, -0.008607, -0.007871, -0.007191 
                           , -0.006512, -0.005946, -0.005379, -0.004926 
                           , -0.004473, -0.004077, -0.003737, -0.003398 
                           , -0.003114, -0.002831, -0.002605, -0.002378 
                           , -0.002208, -0.002039, -0.001869, -0.001699 
                           , -0.001585, -0.001472, -0.001359, -0.001246 
                           , -0.001132, -0.001019, -0.001019, -0.000906 
                           , -0.000906, -0.000793, -0.000793, -0.000680 
                           , -0.000680, -0.000680, -0.000566, -0.000566 
                           , -0.000566, -0.000566, -0.000453, -0.000453 
                           , -0.000453, -0.000453, -0.000453, -0.000453 
                           , -0.000340, -0.000340, -0.000340, -0.000340 
                           , -0.000340, -0.000340, -0.000340, -0.000340 
                           , -0.000340, -0.000340, -0.000340, -0.000340 
                           , -0.000340, -0.000227, -0.000227, -0.000227 
                           , -0.000227, -0.000227, -0.000227, -0.000227 
                           , -0.000227, -0.000227, -0.000227, -0.000227 
                           , -0.000227, -0.000227, -0.000227, -0.000227 
                           , -0.000227, -0.000227, -0.000227, -0.000227 
                           , -0.000227, -0.000227, -0.000227, -0.000227 
                           , -0.000227, -0.000227, -0.000227, -0.000227 
                           , -0.000227, -0.000227, -0.000227, -0.000227 
                           , -0.000227, -0.000227, -0.000227, -0.000227 
                           , -0.000227, -0.000227, -0.000227, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113,  0.000000 
                           ,  0.000000,  0.000000,  0.000000,  0.000000 
                           ,  0.000000,  0.000000,  0.000000,  0.000000 
                           ,  0.000000,  0.000000,  0.000000,  0.000000 
                           ,  0.000000,  0.000000,  0.000000,  0.000000 
                           ,  0.000000,  0.000000,  0.000000,  0.000000 
                           ,  0.000000,  0.000000,  0.000000,  0.000000 
                           ,  0.000000,  0.000000,  0.000000,  0.000000 
                           ,  0.000000,  0.000000,  0.000000,  0.000000 
                           ,  0.000000,  0.000000,  0.000000,  0.000000 
                           ,  0.000000,  0.000000,  0.000000,  0.000000 };

  // increase CrossTalk to measurements
  for (Int_t ipasa = 0; ipasa < kNpasa; ipasa++) {
    xtalk[ipasa] *= 1.75;
  }

  if (fTRFsmp) delete [] fTRFsmp;
  fTRFsmp = new Float_t[fTRFbin];
  if (fCTsmp)  delete [] fCTsmp;
  fCTsmp  = new Float_t[fTRFbin];

  Float_t loTRF    = TMath::Max(fTRFlo / fDriftVelocity,time[0]);
  Float_t hiTRF    = TMath::Min(fTRFhi / fDriftVelocity,time[kNpasa-1]);
  Float_t binWidth = (hiTRF - loTRF) / ((Float_t) fTRFbin);

  // Take the linear interpolation
  for (Int_t iBin = 0; iBin < fTRFbin; iBin++) {

    Float_t bin = (((Float_t) iBin) + 0.5) * binWidth + loTRF;
    ipos1 = ipos2 = 0;
    diff  = 0;
    do {
      diff = bin - time[ipos2++];
    } while (diff > 0);
    ipos2--;
    if (ipos2 >= kNpasa) ipos2 = kNpasa - 1;
    ipos1 = ipos2 - 1;

    fTRFsmp[iBin] = signal[ipos2] 
                  + diff * (signal[ipos2] - signal[ipos1]) 
                         / (  time[ipos2] -   time[ipos1]);

    fCTsmp[iBin]  = xtalk[ipos2] 
                  + diff * (xtalk[ipos2]  -  xtalk[ipos1]) 
                         / (  time[ipos2] -   time[ipos1]);

  }

}

//_____________________________________________________________________________
void AliTRDdigitizer::SamplePRF()
{
  //
  // Samples the pad response function
  //

  const Int_t kNplan  = AliTRDgeometry::kNplan;
  const Int_t kPRFbin = 61;

  Float_t prf[kNplan][kPRFbin] = { { 0.018570, 0.022270, 0.026710, 0.032010
                                   , 0.038350, 0.045920, 0.054930, 0.065650
                                   , 0.078370, 0.093420, 0.111150, 0.131940
                                   , 0.156160, 0.184160, 0.216220, 0.252470
                                   , 0.292860, 0.337030, 0.384330, 0.433750
                                   , 0.484010, 0.533630, 0.581150, 0.625200
                                   , 0.664710, 0.698860, 0.727130, 0.749230
                                   , 0.765050, 0.774540, 0.777700, 0.774540
                                   , 0.765050, 0.749230, 0.727130, 0.698860
                                   , 0.664710, 0.625200, 0.581150, 0.533630
                                   , 0.484010, 0.433750, 0.384330, 0.337030
                                   , 0.292860, 0.252470, 0.216220, 0.184160
                                   , 0.156160, 0.131940, 0.111150, 0.093420
                                   , 0.078370, 0.065650, 0.054930, 0.045920
                                   , 0.038350, 0.032010, 0.026710, 0.022270
				   , 0.018570                               }
                                 , { 0.015730, 0.019040, 0.023030, 0.027840
                                   , 0.033650, 0.040650, 0.049060, 0.059160
                                   , 0.071260, 0.085710, 0.102910, 0.123270
                                   , 0.147240, 0.175220, 0.207590, 0.244540
                                   , 0.286090, 0.331920, 0.381350, 0.433290
                                   , 0.486290, 0.538710, 0.588870, 0.635280
                                   , 0.676760, 0.712460, 0.741890, 0.764810
                                   , 0.781150, 0.790930, 0.794180, 0.790930
                                   , 0.781150, 0.764810, 0.741890, 0.712460
                                   , 0.676760, 0.635280, 0.588870, 0.538710
                                   , 0.486290, 0.433290, 0.381350, 0.331920
                                   , 0.286090, 0.244540, 0.207590, 0.175220
                                   , 0.147240, 0.123270, 0.102910, 0.085710
                                   , 0.071260, 0.059160, 0.049060, 0.040650
                                   , 0.033650, 0.027840, 0.023030, 0.019040
				   , 0.015730                               }
				 , { 0.013330, 0.016270, 0.019850, 0.024210
                                   , 0.029510, 0.035960, 0.043790, 0.053280
                                   , 0.064740, 0.078580, 0.095190, 0.115070
                                   , 0.138700, 0.166570, 0.199120, 0.236660
                                   , 0.279260, 0.326660, 0.378140, 0.432540
                                   , 0.488260, 0.543440, 0.596200, 0.644900
                                   , 0.688240, 0.725380, 0.755840, 0.779470
                                   , 0.796260, 0.806280, 0.809610, 0.806280
                                   , 0.796260, 0.779470, 0.755840, 0.725380
                                   , 0.688240, 0.644900, 0.596200, 0.543440
                                   , 0.488260, 0.432540, 0.378140, 0.326660
                                   , 0.279260, 0.236660, 0.199120, 0.166570
                                   , 0.138700, 0.115070, 0.095190, 0.078580
                                   , 0.064740, 0.053280, 0.043790, 0.035960
                                   , 0.029510, 0.024210, 0.019850, 0.016270
				   , 0.013330                               }
                                 , { 0.011280, 0.013890, 0.017090, 0.021030
                                   , 0.025870, 0.031800, 0.039060, 0.047940
                                   , 0.058790, 0.071980, 0.087990, 0.107330
                                   , 0.130550, 0.158220, 0.190850, 0.228870
                                   , 0.272410, 0.321270, 0.374740, 0.431560
                                   , 0.489960, 0.547870, 0.603180, 0.654080
                                   , 0.699190, 0.737640, 0.769030, 0.793260
                                   , 0.810410, 0.820620, 0.824010, 0.820620
                                   , 0.810410, 0.793260, 0.769030, 0.737640
                                   , 0.699190, 0.654080, 0.603180, 0.547870
                                   , 0.489960, 0.431560, 0.374740, 0.321270
                                   , 0.272410, 0.228870, 0.190850, 0.158220
                                   , 0.130550, 0.107330, 0.087990, 0.071980
                                   , 0.058790, 0.047940, 0.039060, 0.031800
                                   , 0.025870, 0.021030, 0.017090, 0.013890
				   , 0.011280                               }
                                 , { 0.009550, 0.011860, 0.014720, 0.018270
                                   , 0.022660, 0.028100, 0.034820, 0.043120
                                   , 0.053340, 0.065900, 0.081280, 0.100040
                                   , 0.122800, 0.150180, 0.182800, 0.221170
                                   , 0.265550, 0.315790, 0.371180, 0.430370
                                   , 0.491430, 0.552030, 0.609840, 0.662860
                                   , 0.709630, 0.749290, 0.781490, 0.806220
                                   , 0.823650, 0.834000, 0.837430, 0.834000
                                   , 0.823650, 0.806220, 0.781490, 0.749290
                                   , 0.709630, 0.662860, 0.609840, 0.552030
                                   , 0.491430, 0.430370, 0.371180, 0.315790
                                   , 0.265550, 0.221170, 0.182800, 0.150180
                                   , 0.122800, 0.100040, 0.081280, 0.065900
                                   , 0.053340, 0.043120, 0.034820, 0.028100
                                   , 0.022660, 0.018270, 0.014720, 0.011860
				   , 0.009550                               }
                                 , { 0.008080, 0.010120, 0.012670, 0.015860
                                   , 0.019840, 0.024820, 0.031030, 0.038760
                                   , 0.048370, 0.060300, 0.075040, 0.093200
                                   , 0.115430, 0.142450, 0.174980, 0.213610
                                   , 0.258720, 0.310250, 0.367480, 0.429010
                                   , 0.492690, 0.555950, 0.616210, 0.671280
                                   , 0.719600, 0.760350, 0.793250, 0.818380
                                   , 0.836020, 0.846460, 0.849920, 0.846460
                                   , 0.836020, 0.818380, 0.793250, 0.760350
                                   , 0.719600, 0.671280, 0.616210, 0.555950
                                   , 0.492690, 0.429010, 0.367480, 0.310250
                                   , 0.258720, 0.213610, 0.174980, 0.142450
                                   , 0.115430, 0.093200, 0.075040, 0.060300
                                   , 0.048370, 0.038760, 0.031030, 0.024820
                                   , 0.019840, 0.015860, 0.012670, 0.010120
				   , 0.008080                               } };

  fPRFlo  = -1.5;
  fPRFhi  =  1.5;
  fPRFbin = kPRFbin;
  fPRFwid = (fPRFhi - fPRFlo) / ((Float_t) fPRFbin);
  fPRFpad = ((Int_t) (1.0 / fPRFwid));

  if (fPRFsmp) delete [] fPRFsmp;
  fPRFsmp = new Float_t[kNplan*fPRFbin];
  for (Int_t iPla = 0; iPla < kNplan; iPla++) {
    for (Int_t iBin = 0; iBin < fPRFbin; iBin++) {
      fPRFsmp[iPla*kPRFbin+iBin] = prf[iPla][iBin];
    }
  } 

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::Open(const Char_t *file, Int_t nEvent)
{
  //
  // Opens a ROOT-file with TRD-hits and reads in the hit-tree
  //

  // Connect the AliRoot file containing Geometry, Kine, and Hits
  fInputFile = (TFile*) gROOT->GetListOfFiles()->FindObject(file);
  if (!fInputFile) {
    if (fDebug > 0) {
      printf("<AliTRDdigitizer::Open> ");
      printf("Open the AliROOT-file %s.\n",file);
    }
    fInputFile = new TFile(file,"UPDATE");
  }
  else {
    if (fDebug > 0) {
      printf("<AliTRDdigitizer::Open> ");
      printf("%s is already open.\n",file);
    }
  }

  gAlice = (AliRun*) fInputFile->Get("gAlice");
  if (gAlice) {
    if (fDebug > 0) {
      printf("<AliTRDdigitizer::Open> ");
      printf("AliRun object found on file.\n");
    }
  }
  else {
    printf("<AliTRDdigitizer::Open> ");
    printf("Could not find AliRun object.\n");
    return kFALSE;
  }

  fEvent = nEvent;

  // Import the Trees for the event nEvent in the file
  Int_t nparticles = gAlice->GetEvent(fEvent);
  if (nparticles <= 0) {
    printf("<AliTRDdigitizer::Open> ");
    printf("No entries in the trees for event %d.\n",fEvent);
    return kFALSE;
  }

  if (InitDetector()) {
    return MakeBranch();
  }
  else {
    return kFALSE;
  }

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::InitDetector()
{
  //
  // Sets the pointer to the TRD detector and the geometry
  //

  // Get the pointer to the detector class and check for version 1
  fTRD = (AliTRD*) gAlice->GetDetector("TRD");
  if (fTRD->IsVersion() != 1) {
    printf("<AliTRDdigitizer::InitDetector> ");
    printf("TRD must be version 1 (slow simulator).\n");
    exit(1);
  }

  // Get the geometry
  fGeo = fTRD->GetGeometry();
  if (fDebug > 0) {
    printf("<AliTRDdigitizer::InitDetector> ");
    printf("Geometry version %d\n",fGeo->IsVersion());
  }

  // The magnetic field strength in Tesla
  fField = 0.2 * gAlice->Field()->Factor();

  // Create a digits manager
  fDigitsManager = new AliTRDdigitsManager();
  fDigitsManager->SetSDigits(fSDigits);
  fDigitsManager->CreateArrays();
  fDigitsManager->SetEvent(fEvent);
  fDigitsManager->SetDebug(fDebug);

  // The list for the input s-digits manager to be merged
  fSDigitsManagerList = new TList();

  return ReInit();

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::MakeBranch(const Char_t *file)
{
  // 
  // Create the branches for the digits array
  //

  return fDigitsManager->MakeBranch(file);

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::MakeDigits()
{
  //
  // Creates digits.
  //

  ///////////////////////////////////////////////////////////////
  // Parameter 
  ///////////////////////////////////////////////////////////////

  // Converts number of electrons to fC
  const Double_t kEl2fC  = 1.602E-19 * 1.0E15; 

  ///////////////////////////////////////////////////////////////

  // Number of pads included in the pad response
  const Int_t kNpad  = 3;

  // Number of track dictionary arrays
  const Int_t kNDict = AliTRDdigitsManager::kNDict;

  // Half the width of the amplification region
  const Float_t kAmWidth = AliTRDgeometry::AmThick() / 2.;

  Int_t   iRow, iCol, iTime, iPad;
  Int_t   iDict  = 0;
  Int_t   nBytes = 0;

  Int_t   totalSizeDigits = 0;
  Int_t   totalSizeDict0  = 0;
  Int_t   totalSizeDict1  = 0;
  Int_t   totalSizeDict2  = 0;

  Int_t   timeTRDbeg = 0;
  Int_t   timeTRDend = 1;

  Float_t pos[3];
  Float_t rot[3];
  Float_t xyz[3];
  Float_t padSignal[kNpad];
  Float_t signalOld[kNpad];

  AliTRDdataArrayF *signals = 0;
  AliTRDdataArrayI *digits  = 0;
  AliTRDdataArrayI *dictionary[kNDict];

  // Create a container for the amplitudes
  AliTRDsegmentArray *signalsArray 
                     = new AliTRDsegmentArray("AliTRDdataArrayF"
                                             ,AliTRDgeometry::Ndet());

  if (fTRFOn) {
    timeTRDbeg = ((Int_t) (-fTRFlo / fGeo->GetTimeBinSize())) - 1;
    timeTRDend = ((Int_t) ( fTRFhi / fGeo->GetTimeBinSize())) - 1;
    if (fDebug > 0) {
      printf("<AliTRDdigitizer::MakeDigits> ");
      printf("Sample the TRF between -%d and %d\n",timeTRDbeg,timeTRDend);
    }
  }

  Float_t elAttachProp = fElAttachProp / 100.; 

  // Create the sampled PRF
  SamplePRF();

  // Create the sampled TRF
  SampleTRF();

  if (!fGeo) {
    printf("<AliTRDdigitizer::MakeDigits> ");
    printf("No geometry defined\n");
    return kFALSE;
  }

  if (fDebug > 0) {
    printf("<AliTRDdigitizer::MakeDigits> ");
    printf("Start creating digits.\n");
  }

  // Get the pointer to the hit tree
  TTree *HitTree = gAlice->TreeH();

  // Get the number of entries in the hit tree
  // (Number of primary particles creating a hit somewhere)
  Int_t nTrack = (Int_t) HitTree->GetEntries();
  if (fDebug > 0) {
    printf("<AliTRDdigitizer::MakeDigits> ");
    printf("Found %d primary particles\n",nTrack);
  } 

  Int_t detectorOld = -1;
  Int_t countHits   =  0; 

  // Loop through all entries in the tree
  for (Int_t iTrack = 0; iTrack < nTrack; iTrack++) {

    gAlice->ResetHits();
    nBytes += HitTree->GetEvent(iTrack);

    // Loop through the TRD hits
    Int_t iHit = 0;
    AliTRDhit *hit = (AliTRDhit *) fTRD->FirstHit(-1);
    while (hit) {
 
      countHits++;
      iHit++;

              pos[0]      = hit->X();
              pos[1]      = hit->Y();
              pos[2]      = hit->Z();
      Float_t q           = hit->GetCharge();
      Int_t   track       = hit->Track();
      Int_t   detector    = hit->GetDetector();
      Int_t   plane       = fGeo->GetPlane(detector);
      Int_t   sector      = fGeo->GetSector(detector);
      Int_t   chamber     = fGeo->GetChamber(detector);
      Int_t   nRowMax     = fGeo->GetRowMax(plane,chamber,sector);
      Int_t   nColMax     = fGeo->GetColMax(plane);
      Int_t   nTimeMax    = fGeo->GetTimeMax();
      Int_t   nTimeBefore = fGeo->GetTimeBefore();
      Int_t   nTimeAfter  = fGeo->GetTimeAfter();
      Int_t   nTimeTotal  = fGeo->GetTimeTotal();
      Float_t row0        = fGeo->GetRow0(plane,chamber,sector);
      Float_t col0        = fGeo->GetCol0(plane);
      Float_t time0       = fGeo->GetTime0(plane);
      Float_t rowPadSize  = fGeo->GetRowPadSize(plane,chamber,sector);
      Float_t colPadSize  = fGeo->GetColPadSize(plane);
      Float_t timeBinSize = fGeo->GetTimeBinSize();
      Float_t divideRow   = 1.0 / rowPadSize;
      Float_t divideCol   = 1.0 / colPadSize;
      Float_t divideTime  = 1.0 / timeBinSize;

      if (fDebug > 1) {
        printf("Analyze hit no. %d ",iHit);
        printf("-----------------------------------------------------------\n");
        hit->Dump();
        printf("plane = %d, sector = %d, chamber = %d\n"
              ,plane,sector,chamber);
        printf("nRowMax = %d, nColMax = %d, nTimeMax = %d\n" 
              ,nRowMax,nColMax,nTimeMax);
        printf("nTimeBefore = %d, nTimeAfter = %d, nTimeTotal = %d\n"
	      ,nTimeBefore,nTimeAfter,nTimeTotal);
        printf("row0 = %f, col0 = %f, time0 = %f\n"
              ,row0,col0,time0);
        printf("rowPadSize = %f, colPadSize = %f, timeBinSize = %f\n"
	       ,rowPadSize,colPadSize,timeBinSize); 
      }
       
      // Don't analyze test hits and switched off detectors
      if ((CheckDetector(plane,chamber,sector)) &&
          (((Int_t) q) != 0)) {

        if (detector != detectorOld) {

          if (fDebug > 1) {
            printf("<AliTRDdigitizer::MakeDigits> ");
            printf("Get new container. New det = %d, Old det = %d\n"
                  ,detector,detectorOld);
	  }
          // Compress the old one if enabled
          if ((fCompress) && (detectorOld > -1)) {
            if (fDebug > 1) {
              printf("<AliTRDdigitizer::MakeDigits> ");
              printf("Compress the old container ...");
	    }
            signals->Compress(1,0);
            for (iDict = 0; iDict < kNDict; iDict++) {
              dictionary[iDict]->Compress(1,0);
	    }
            if (fDebug > 1) printf("done\n");
	  }
	  // Get the new container
          signals = (AliTRDdataArrayF *) signalsArray->At(detector);
          if (signals->GetNtime() == 0) {
            // Allocate a new one if not yet existing
            if (fDebug > 1) {
              printf("<AliTRDdigitizer::MakeDigits> ");
              printf("Allocate a new container ... ");
	    }
            signals->Allocate(nRowMax,nColMax,nTimeTotal);
	  }
          else {
	    // Expand an existing one
            if (fCompress) {
              if (fDebug > 1) {
                printf("<AliTRDdigitizer::MakeDigits> ");
                printf("Expand an existing container ... ");
	      }
              signals->Expand();
	    }
	  }
	  // The same for the dictionary
          for (iDict = 0; iDict < kNDict; iDict++) {       
            dictionary[iDict] = fDigitsManager->GetDictionary(detector,iDict);
            if (dictionary[iDict]->GetNtime() == 0) {
              dictionary[iDict]->Allocate(nRowMax,nColMax,nTimeTotal);
	    }
            else {
              if (fCompress) dictionary[iDict]->Expand();
	    }
          }      
          if (fDebug > 1) printf("done\n");
          detectorOld = detector;
        }

        // Rotate the sectors on top of each other
        fGeo->Rotate(detector,pos,rot);

        // The driftlength. It is negative if the hit is in the 
        // amplification region.
        Float_t driftlength = time0 - rot[0];

        // Take also the drift in the amplification region into account
        // The drift length is at the moment still the same, regardless of
        // the position relativ to the wire. This non-isochronity needs still
        // to be implemented.
        Float_t driftlengthL = TMath::Abs(driftlength + kAmWidth);
        if (fExBOn) driftlengthL /= TMath::Sqrt(fLorentzFactor);

        // Loop over all electrons of this hit
        // TR photons produce hits with negative charge
        Int_t nEl = ((Int_t) TMath::Abs(q));
        for (Int_t iEl = 0; iEl < nEl; iEl++) {

          xyz[0] = rot[0];
          xyz[1] = rot[1];
          xyz[2] = rot[2];

          // Electron attachment
          if (fElAttachOn) {
            if (gRandom->Rndm() < (driftlengthL * elAttachProp)) 
              continue;
          }

          // Apply the diffusion smearing
          if (fDiffusionOn) {
            if (!(Diffusion(driftlengthL,xyz))) continue;
	  }

          // Apply E x B effects (depends on drift direction)
          if (fExBOn) { 
            if (!(ExB(driftlength+kAmWidth,xyz))) continue;   
	  }

          // The electron position after diffusion and ExB in pad coordinates 
          // The pad row (z-direction)
          Float_t rowDist   = xyz[2] - row0;
          Int_t   rowE      = ((Int_t) (rowDist * divideRow));
          if ((rowE < 0) || (rowE >= nRowMax)) continue;   
          Float_t rowOffset = ((((Float_t) rowE) + 0.5) * rowPadSize) - rowDist;

          // The pad column (rphi-direction)
          Float_t col0tilt  =  Col0Tilted(col0,rowOffset,plane);
          Float_t colDist   = xyz[1] - col0tilt;
          Int_t   colE      = ((Int_t) (colDist * divideCol));
          if ((colE < 0) || (colE >= nColMax)) continue;   
          Float_t colOffset = ((((Float_t) colE) + 0.5) * colPadSize) - colDist;    

          // The time bin (negative for hits in the amplification region)
	  // In the amplification region the electrons drift from both sides
	  // to the middle (anode wire plane)
          Float_t timeDist   = time0 - xyz[0];
          Float_t timeOffset = 0;
          Int_t   timeE      = 0;
          if (timeDist > 0) {
	    // The time bin
            timeE      = ((Int_t) (timeDist * divideTime));
            // The distance of the position to the middle of the timebin
            timeOffset = ((((Float_t) timeE) + 0.5) * timeBinSize) - timeDist;
	  }
          else {
	    // Difference between half of the amplification gap width and
	    // the distance to the anode wire
            Float_t anodeDist = kAmWidth - TMath::Abs(timeDist + kAmWidth);
            // The time bin
            timeE      = -1 * (((Int_t ) (anodeDist * divideTime)) + 1);
            // The distance of the position to the middle of the timebin
            timeOffset = ((((Float_t) timeE) + 0.5) * timeBinSize) + anodeDist;
	  }
 
          // Apply the gas gain including fluctuations
          Float_t ggRndm = 0.0;
          do {
            ggRndm = gRandom->Rndm();
	  } while (ggRndm <= 0);
          Int_t signal = (Int_t) (-fGasGain * TMath::Log(ggRndm));

          // Apply the pad response 
          if (fPRFOn) {
  	    // The distance of the electron to the center of the pad 
	    // in units of pad width
            Float_t dist = - colOffset * divideCol;
            if (!(PadResponse(signal,dist,plane,padSignal))) continue;
	  }
	  else {
            padSignal[0] = 0.0;
            padSignal[1] = signal;
            padSignal[2] = 0.0;
	  }

	  // Sample the time response inside the drift region
	  // + additional time bins before and after.
          // The sampling is done always in the middle of the time bin
          for (Int_t iTimeBin = TMath::Max(timeE-timeTRDbeg,        -nTimeBefore) 
  	            ;iTimeBin < TMath::Min(timeE+timeTRDend,nTimeMax+nTimeAfter ) 
        	    ;iTimeBin++) {

     	    // Apply the time response
            Float_t timeResponse = 1.0;
            Float_t crossTalk    = 0.0;
            Float_t time         = (iTimeBin - timeE) * timeBinSize + timeOffset;
            if (fTRFOn) {
              timeResponse = TimeResponse(time);
	    }
            if (fCTOn) {
              crossTalk    = CrossTalk(time);
            }

            signalOld[0] = 0.0;
            signalOld[1] = 0.0;
            signalOld[2] = 0.0;

            for (iPad = 0; iPad < kNpad; iPad++) {

              Int_t colPos = colE + iPad - 1;
              if (colPos <        0) continue;
              if (colPos >= nColMax) break;

              // Add the signals
              // Note: The time bin number is shifted by nTimeBefore to avoid negative
              // time bins. This has to be subtracted later.
              Int_t iCurrentTimeBin = iTimeBin + nTimeBefore;
              signalOld[iPad]  = signals->GetDataUnchecked(rowE,colPos,iCurrentTimeBin);
              if( colPos != colE ) {
                signalOld[iPad] += padSignal[iPad] * (timeResponse + crossTalk);
              } 
              else {
                signalOld[iPad] += padSignal[iPad] * timeResponse;
              }
              signals->SetDataUnchecked(rowE,colPos,iCurrentTimeBin,signalOld[iPad]);

              // Store the track index in the dictionary
              // Note: We store index+1 in order to allow the array to be compressed
              if (signalOld[iPad] > 0) {
                for (iDict = 0; iDict < kNDict; iDict++) {
                  Int_t oldTrack = dictionary[iDict]->GetDataUnchecked(rowE
                                                                      ,colPos
                                                                      ,iCurrentTimeBin);
                  if (oldTrack == track+1) break;
                  if (oldTrack ==       0) {
                    dictionary[iDict]->SetDataUnchecked(rowE,colPos,iCurrentTimeBin,track+1);
                    break;
                  }
                }
              }

	    } // Loop: pads

	  } // Loop: time bins

        } // Loop: electrons of a single hit

      } // If: detector and test hit

      hit = (AliTRDhit *) fTRD->NextHit();   

    } // Loop: hits of one primary track

  } // Loop: primary tracks

  if (fDebug > 0) {
    printf("<AliTRDdigitizer::MakeDigits> ");
    printf("Finished analyzing %d hits\n",countHits);
  }

  // The total conversion factor
  Float_t convert = kEl2fC * fPadCoupling * fTimeCoupling * fChipGain;

  // Loop through all chambers to finalize the digits
  for (Int_t iDet = 0; iDet < AliTRDgeometry::Ndet(); iDet++) {

    Int_t plane       = fGeo->GetPlane(iDet);
    Int_t sector      = fGeo->GetSector(iDet);
    Int_t chamber     = fGeo->GetChamber(iDet);
    Int_t nRowMax     = fGeo->GetRowMax(plane,chamber,sector);
    Int_t nColMax     = fGeo->GetColMax(plane);
    Int_t nTimeMax    = fGeo->GetTimeMax();
    Int_t nTimeTotal  = fGeo->GetTimeTotal();

    Double_t *inADC  = new Double_t[nTimeTotal];
    Double_t *outADC = new Double_t[nTimeTotal];

    if (fDebug > 0) {
      printf("<AliTRDdigitizer::MakeDigits> ");
      printf("Digitization for chamber %d\n",iDet);
    }

    // Add a container for the digits of this detector
    digits = fDigitsManager->GetDigits(iDet);        
    // Allocate memory space for the digits buffer
    digits->Allocate(nRowMax,nColMax,nTimeTotal);

    // Get the signal container
    signals = (AliTRDdataArrayF *) signalsArray->At(iDet);
    if (signals->GetNtime() == 0) {
      // Create missing containers
      signals->Allocate(nRowMax,nColMax,nTimeTotal);      
    }
    else {
      // Expand the container if neccessary
      if (fCompress) signals->Expand();
    }
    // Create the missing dictionary containers
    for (iDict = 0; iDict < kNDict; iDict++) {       
      dictionary[iDict] = fDigitsManager->GetDictionary(iDet,iDict);
      if (dictionary[iDict]->GetNtime() == 0) {
        dictionary[iDict]->Allocate(nRowMax,nColMax,nTimeTotal);
      }
    }

    Int_t nDigits = 0;

    // Don't create noise in detectors that are switched off
    if (CheckDetector(plane,chamber,sector)) {

      // Create the digits for this chamber
      for (iRow  = 0; iRow  <  nRowMax;   iRow++ ) {
        for (iCol  = 0; iCol  <  nColMax;   iCol++ ) {

	  // Create summable digits
          if (fSDigits) {

            for (iTime = 0; iTime < nTimeTotal; iTime++) {         
              Float_t signalAmp = signals->GetDataUnchecked(iRow,iCol,iTime);
              signalAmp *= fSDigitsScale;
              signalAmp  = TMath::Min(signalAmp,(Float_t) 1.0e9);
              Int_t adc  = (Int_t) signalAmp;
              nDigits++;
              digits->SetDataUnchecked(iRow,iCol,iTime,adc);
	    }

	  }
	  // Create normal digits
          else {

            for (iTime = 0; iTime < nTimeTotal; iTime++) {         
              Float_t signalAmp = signals->GetDataUnchecked(iRow,iCol,iTime);
              // Add the noise
              signalAmp  = TMath::Max((Double_t) gRandom->Gaus(signalAmp,fNoise),0.0);
              // Convert to mV
              signalAmp *= convert;
	      // Convert to ADC counts. Set the overflow-bit fADCoutRange if the 
	      // signal is larger than fADCinRange
              Int_t adc  = 0;
              if (signalAmp >= fADCinRange) {
                adc = ((Int_t) fADCoutRange);
	      }
              else {
                adc = ((Int_t) (signalAmp * (fADCoutRange / fADCinRange)));
	      }
              inADC[iTime]  = adc;
              outADC[iTime] = adc;
	    }

	    // Apply the tail cancelation via the digital filter
            if (fTCOn) {
              DeConvExp(inADC,outADC,nTimeTotal,fTCnexp);
	    }

            for (iTime = 0; iTime < nTimeTotal; iTime++) {   
              // Store the amplitude of the digit if above threshold
              if (outADC[iTime] > fADCthreshold) {
                if (fDebug > 2) {
                  printf("  iRow = %d, iCol = %d, iTime = %d, adc = %f\n"
                        ,iRow,iCol,iTime,outADC[iTime]);
	        }
                nDigits++;
                digits->SetDataUnchecked(iRow,iCol,iTime,outADC[iTime]);
  	      }
	    }

	  }

        }
      }

    }

    // Compress the arrays
    digits->Compress(1,0);
    for (iDict = 0; iDict < kNDict; iDict++) {
      dictionary[iDict]->Compress(1,0);
    }

    totalSizeDigits += digits->GetSize();
    totalSizeDict0  += dictionary[0]->GetSize();
    totalSizeDict1  += dictionary[1]->GetSize();
    totalSizeDict2  += dictionary[2]->GetSize();

    Float_t nPixel = nRowMax * nColMax * nTimeMax;
    if (fDebug > 0) {
      printf("<AliTRDdigitizer::MakeDigits> ");
      printf("Found %d digits in detector %d (%3.0f).\n"
            ,nDigits,iDet
            ,100.0 * ((Float_t) nDigits) / nPixel);
    } 

    if (fCompress) signals->Compress(1,0);

    delete [] inADC;
    delete [] outADC;

  }

  if (fDebug > 0) {
    printf("<AliTRDdigitizer::MakeDigits> ");
    printf("Total number of analyzed hits = %d\n",countHits);
    printf("<AliTRDdigitizer::MakeDigits> ");
    printf("Total digits data size = %d, %d, %d, %d\n",totalSizeDigits
                                                      ,totalSizeDict0
                                                      ,totalSizeDict1
                                                      ,totalSizeDict2);        
  }

  return kTRUE;

}

//_____________________________________________________________________________
void AliTRDdigitizer::AddSDigitsManager(AliTRDdigitsManager *man)
{
  //
  // Add a digits manager for s-digits to the input list.
  //

  fSDigitsManagerList->Add(man);

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::ConvertSDigits()
{
  //
  // Converts s-digits to normal digits
  //

  // Number of track dictionary arrays
  const Int_t    kNDict = AliTRDdigitsManager::kNDict;

  // Converts number of electrons to fC
  const Double_t kEl2fC = 1.602E-19 * 1.0E15; 

  Int_t iDict = 0;
  Int_t iRow;
  Int_t iCol;
  Int_t iTime;

  if (fDebug > 0) {
    this->Dump();
  }

  Double_t sDigitsScale = 1.0 / GetSDigitsScale();
  Double_t noise        = GetNoise();
  Double_t padCoupling  = GetPadCoupling();
  Double_t timeCoupling = GetTimeCoupling();
  Double_t chipGain     = GetChipGain();
  Double_t convert      = kEl2fC * padCoupling * timeCoupling * chipGain;;
  Double_t adcInRange   = GetADCinRange();
  Double_t adcOutRange  = GetADCoutRange();
  Int_t    adcThreshold = GetADCthreshold();

  AliTRDdataArrayI *digitsIn;
  AliTRDdataArrayI *digitsOut;
  AliTRDdataArrayI *dictionaryIn[kNDict];
  AliTRDdataArrayI *dictionaryOut[kNDict];

  // Loop through the detectors
  for (Int_t iDet = 0; iDet < AliTRDgeometry::Ndet(); iDet++) {

    if (fDebug > 0) {
      printf("<AliTRDdigitizer::ConvertSDigits> ");
      printf("Convert detector %d to digits.\n",iDet);
    }

    Int_t plane      = fGeo->GetPlane(iDet);
    Int_t sector     = fGeo->GetSector(iDet);
    Int_t chamber    = fGeo->GetChamber(iDet);
    Int_t nRowMax    = fGeo->GetRowMax(plane,chamber,sector);
    Int_t nColMax    = fGeo->GetColMax(plane);
    Int_t nTimeTotal = fGeo->GetTimeTotal();

    Double_t *inADC  = new Double_t[nTimeTotal];
    Double_t *outADC = new Double_t[nTimeTotal];

    digitsIn  = fSDigitsManager->GetDigits(iDet);
    digitsIn->Expand();
    digitsOut = fDigitsManager->GetDigits(iDet);
    digitsOut->Allocate(nRowMax,nColMax,nTimeTotal);
    for (iDict = 0; iDict < kNDict; iDict++) {
      dictionaryIn[iDict]  = fSDigitsManager->GetDictionary(iDet,iDict);
      dictionaryIn[iDict]->Expand();
      dictionaryOut[iDict] = fDigitsManager->GetDictionary(iDet,iDict);
      dictionaryOut[iDict]->Allocate(nRowMax,nColMax,nTimeTotal);
    }

    for (iRow  = 0; iRow  <  nRowMax;   iRow++ ) {
      for (iCol  = 0; iCol  <  nColMax;   iCol++ ) {

        for (iTime = 0; iTime < nTimeTotal; iTime++) {         
          Double_t signal = (Double_t) digitsIn->GetDataUnchecked(iRow,iCol,iTime);
          signal *= sDigitsScale;
          // Add the noise
          signal  = TMath::Max((Double_t) gRandom->Gaus(signal,noise),0.0);
          // Convert to mV
          signal *= convert;
	  // Convert to ADC counts. Set the overflow-bit adcOutRange if the 
	  // signal is larger than adcInRange
          Int_t adc  = 0;
          if (signal >= adcInRange) {
            adc = ((Int_t) adcOutRange);
	  }
          else {
            adc = ((Int_t) (signal * (adcOutRange / adcInRange)));
	  }
          inADC[iTime]  = adc;
          outADC[iTime] = adc;
	}

	// Apply the tail cancelation via the digital filter
        if (fTCOn) {
          DeConvExp(inADC,outADC,nTimeTotal,fTCnexp);
	}

        for (iTime = 0; iTime < nTimeTotal; iTime++) {   
          // Store the amplitude of the digit if above threshold
          if (outADC[iTime] > adcThreshold) {
            digitsOut->SetDataUnchecked(iRow,iCol,iTime,outADC[iTime]);
  	    // Copy the dictionary
            for (iDict = 0; iDict < kNDict; iDict++) { 
              Int_t track = dictionaryIn[iDict]->GetDataUnchecked(iRow,iCol,iTime);
              dictionaryOut[iDict]->SetDataUnchecked(iRow,iCol,iTime,track);
	    }
	  }
	}

      }
    }

    if (fCompress) {
      digitsIn->Compress(1,0);
      digitsOut->Compress(1,0);
      for (iDict = 0; iDict < kNDict; iDict++) {
        dictionaryIn[iDict]->Compress(1,0);
        dictionaryOut[iDict]->Compress(1,0);
      }
    }

    delete [] inADC;
    delete [] outADC;

  }    

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::MergeSDigits()
{
  //
  // Merges the input s-digits:
  //   - The amplitude of the different inputs are summed up.
  //   - Of the track IDs from the input dictionaries only one is
  //     kept for each input. This works for maximal 3 different merged inputs.
  //

  // Number of track dictionary arrays
  const Int_t kNDict = AliTRDdigitsManager::kNDict;

  Int_t iDict = 0;
  Int_t jDict = 0;

  AliTRDdataArrayI *digitsA;
  AliTRDdataArrayI *digitsB;
  AliTRDdataArrayI *dictionaryA[kNDict];
  AliTRDdataArrayI *dictionaryB[kNDict];

  // Get the first s-digits
  fSDigitsManager = (AliTRDdigitsManager *) fSDigitsManagerList->First();
  if (!fSDigitsManager) return kFALSE;

  // Loop through the other sets of s-digits
  AliTRDdigitsManager *mergeSDigitsManager;
  mergeSDigitsManager = (AliTRDdigitsManager *) 
                        fSDigitsManagerList->After(fSDigitsManager);

  if (fDebug > 0) {
    if (mergeSDigitsManager) {
      printf("<AliTRDdigitizer::MergeSDigits> ");
      printf("Merge %d input files.\n",fSDigitsManagerList->GetSize());
    }
    else {
      printf("<AliTRDdigitizer::MergeSDigits> ");
      printf("Only one input file.\n");
    }
  }

  Int_t iMerge = 0;
  while (mergeSDigitsManager) {

    iMerge++;

    // Loop through the detectors
    for (Int_t iDet = 0; iDet < AliTRDgeometry::Ndet(); iDet++) {

      Int_t plane      = fGeo->GetPlane(iDet);
      Int_t sector     = fGeo->GetSector(iDet);
      Int_t chamber    = fGeo->GetChamber(iDet);
      Int_t nRowMax    = fGeo->GetRowMax(plane,chamber,sector);
      Int_t nColMax    = fGeo->GetColMax(plane);
      Int_t nTimeTotal = fGeo->GetTimeTotal();

      // Loop through the pixels of one detector and add the signals
      digitsA = fSDigitsManager->GetDigits(iDet);
      digitsB = mergeSDigitsManager->GetDigits(iDet);
      digitsA->Expand();
      digitsB->Expand();
      for (iDict = 0; iDict < kNDict; iDict++) {
        dictionaryA[iDict] = fSDigitsManager->GetDictionary(iDet,iDict);
        dictionaryB[iDict] = mergeSDigitsManager->GetDictionary(iDet,iDict);
        dictionaryA[iDict]->Expand();
        dictionaryB[iDict]->Expand();
      }

      if (fDebug > 0) {
        printf("<AliTRDdigitizer::MergeSDigits> ");
        printf("Merge detector %d of input no.%d\n",iDet,iMerge+1);
      }

      for (Int_t iRow  = 0; iRow  <  nRowMax;   iRow++ ) {
        for (Int_t iCol  = 0; iCol  <  nColMax;   iCol++ ) {
          for (Int_t iTime = 0; iTime < nTimeTotal; iTime++) {         

	    // Add the amplitudes of the summable digits 
            Int_t ampA = digitsA->GetDataUnchecked(iRow,iCol,iTime);
            Int_t ampB = digitsB->GetDataUnchecked(iRow,iCol,iTime);
            ampA += ampB;
            digitsA->SetDataUnchecked(iRow,iCol,iTime,ampA);

	    // Add the mask to the track id if defined.
            for (iDict = 0; iDict < kNDict; iDict++) {
              Int_t trackB = dictionaryB[iDict]->GetDataUnchecked(iRow,iCol,iTime);
              if ((fMasks) && (trackB > 0)) {
                for (jDict = 0; jDict < kNDict; jDict++) { 
                  Int_t trackA = dictionaryA[iDict]->GetDataUnchecked(iRow,iCol,iTime);
                  if (trackA == 0) {
                    trackA = trackB + fMasks[iMerge];
                    dictionaryA[iDict]->SetDataUnchecked(iRow,iCol,iTime,trackA);
		  }
		}
	      }
	    }

	  }
	}
      }

      if (fCompress) {
        digitsA->Compress(1,0);
        digitsB->Compress(1,0);
        for (iDict = 0; iDict < kNDict; iDict++) {
          dictionaryA[iDict]->Compress(1,0);
          dictionaryB[iDict]->Compress(1,0);
        }
      }

    }    

    // The next set of s-digits
    mergeSDigitsManager = (AliTRDdigitsManager *) 
                          fSDigitsManagerList->After(mergeSDigitsManager);

  }

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::SDigits2Digits()
{
  //
  // Merges the input s-digits and converts them to normal digits
  //

  if (!MergeSDigits()) return kFALSE;

  return ConvertSDigits();

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::CheckDetector(Int_t plane, Int_t chamber, Int_t sector)
{
  //
  // Checks whether a detector is enabled
  //

  if ((fTRD->GetSensChamber() >=       0) &&
      (fTRD->GetSensChamber() != chamber)) return kFALSE;
  if ((fTRD->GetSensPlane()   >=       0) &&
      (fTRD->GetSensPlane()   !=   plane)) return kFALSE;
  if ( fTRD->GetSensSector()  >=       0) {
    Int_t sens1 = fTRD->GetSensSector();
    Int_t sens2 = sens1 + fTRD->GetSensSectorRange();
    sens2 -= ((Int_t) (sens2 / AliTRDgeometry::Nsect())) 
           * AliTRDgeometry::Nsect();
    if (sens1 < sens2) {
      if ((sector < sens1) || (sector >= sens2)) return kFALSE;
    }
    else {
      if ((sector < sens1) && (sector >= sens2)) return kFALSE;
    }
  }

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::WriteDigits()
{
  //
  // Writes out the TRD-digits and the dictionaries
  //

  // Store the digits and the dictionary in the tree
  return fDigitsManager->WriteDigits();

}

//_____________________________________________________________________________
void AliTRDdigitizer::SetTiltingAngle(Float_t v)
{
  //
  // Set the tilting angle for the readout pads
  //

  fTiltingAngle = TMath::Tan(TMath::Pi()/180.0 * v);

}

//_____________________________________________________________________________
Float_t AliTRDdigitizer::GetTiltingAngle() const
{
  //
  // Get the tilting angle for the readout pads
  //

  return TMath::ATan(180.0/TMath::Pi() * fTiltingAngle);

}

//_____________________________________________________________________________
Float_t AliTRDdigitizer::GetDiffusionL(Float_t vd, Float_t b)
{
  //
  // Returns the longitudinal diffusion coefficient for a given drift 
  // velocity <vd> and a B-field <b> for Xe/CO2 (15%).
  // The values are according to a GARFIELD simulation.
  //

  const Int_t kNb = 5;
  Float_t p0[kNb] = {  0.007440,  0.007493,  0.007513,  0.007672,  0.007831 };
  Float_t p1[kNb] = {  0.019252,  0.018912,  0.018636,  0.018012,  0.017343 };
  Float_t p2[kNb] = { -0.005042, -0.004926, -0.004867, -0.004650, -0.004424 };
  Float_t p3[kNb] = {  0.000195,  0.000189,  0.000195,  0.000182,  0.000169 };

  Int_t ib = ((Int_t) (10 * (b - 0.15)));
  ib       = TMath::Max(  0,ib);
  ib       = TMath::Min(kNb,ib);

  Float_t diff = p0[ib] 
               + p1[ib] * vd
               + p2[ib] * vd*vd
               + p3[ib] * vd*vd*vd;

  return diff;

}

//_____________________________________________________________________________
Float_t AliTRDdigitizer::GetDiffusionT(Float_t vd, Float_t b)
{
  //
  // Returns the transverse diffusion coefficient for a given drift 
  // velocity <vd> and a B-field <b> for Xe/CO2 (15%).
  // The values are according to a GARFIELD simulation.
  //

  const Int_t kNb = 5;
  Float_t p0[kNb] = {  0.009550,  0.009599,  0.009674,  0.009757,  0.009850 };
  Float_t p1[kNb] = {  0.006667,  0.006539,  0.006359,  0.006153,  0.005925 };
  Float_t p2[kNb] = { -0.000853, -0.000798, -0.000721, -0.000635, -0.000541 };
  Float_t p3[kNb] = {  0.000131,  0.000122,  0.000111,  0.000098,  0.000085 };

  Int_t ib = ((Int_t) (10 * (b - 0.15)));
  ib       = TMath::Max(  0,ib);
  ib       = TMath::Min(kNb,ib);

  Float_t diff = p0[ib] 
               + p1[ib] * vd
               + p2[ib] * vd*vd
               + p3[ib] * vd*vd*vd;

  return diff;

}

//_____________________________________________________________________________
Float_t AliTRDdigitizer::GetOmegaTau(Float_t vd, Float_t b)
{
  //
  // Returns omega*tau (tan(Lorentz-angle)) for a given drift velocity <vd> 
  // and a B-field <b> for Xe/CO2 (15%).
  // The values are according to a GARFIELD simulation.
  //

  const Int_t kNb = 5;
  Float_t p0[kNb] = {  0.004810,  0.007412,  0.010252,  0.013409,  0.016888 };
  Float_t p1[kNb] = {  0.054875,  0.081534,  0.107333,  0.131983,  0.155455 };
  Float_t p2[kNb] = { -0.008682, -0.012896, -0.016987, -0.020880, -0.024623 };
  Float_t p3[kNb] = {  0.000155,  0.000238,  0.000330,  0.000428,  0.000541 };

  Int_t ib = ((Int_t) (10 * (b - 0.15)));
  ib       = TMath::Max(  0,ib);
  ib       = TMath::Min(kNb,ib);

  Float_t alphaL = p0[ib] 
                 + p1[ib] * vd
                 + p2[ib] * vd*vd
                 + p3[ib] * vd*vd*vd;

  return TMath::Tan(alphaL);

}

//_____________________________________________________________________________
void AliTRDdigitizer::DeConvExp(Double_t *source, Double_t *target
                              , Int_t n, Int_t nexp) 
{
  //
  // Does the deconvolution by the digital filter.
  //
  // Author:        Marcus Gutfleisch, KIP Heidelberg
  // Optimized for: New TRF from Venelin Angelov, simulated with CADENCE
  //                Pad-ground capacitance = 25 pF
  //                Pad-pad cross talk capacitance = 6 pF
  //                For 10 MHz digitization, corresponding to 20 time bins
  //                in the drift region
  //

  Double_t rates[2];
  Double_t coefficients[2];

  /* initialize (coefficient = alpha, rates = lambda) */
  
  if( nexp == 1 ) {
    rates[0] = 0.466998;
    /* no rescaling */
    coefficients[0] = 1.0;
  }
  if( nexp == 2 ) {
    rates[0] = 0.8988162;
    coefficients[0] = 0.11392069;
    rates[1] = 0.3745688;
    coefficients[1] = 0.8860793;
    /* no rescaling */
    Float_t sumc = coefficients[0]+coefficients[1];
    coefficients[0] /= sumc;
    coefficients[1] /= sumc;
  }
      
  Int_t i, k;
  Double_t reminder[2];
  Double_t correction, result;

  /* attention: computation order is important */
  correction=0.0;
  for ( k = 0; k < nexp; k++ ) reminder[k]=0.0;
    
  for ( i = 0; i < n; i++ ) {
    result = ( source[i] - correction );    /* no rescaling */
    target[i] = result;
    
    for ( k = 0; k < nexp; k++ ) reminder[k] = rates[k] 
                             * ( reminder[k] + coefficients[k] * result);
      
    correction=0.0;
    for ( k = 0; k < nexp; k++ ) correction += reminder[k];
  }
  
}

