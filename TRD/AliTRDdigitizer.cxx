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

#include <stdlib.h>

#include <TMath.h>
#include <TVector.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>

#include "AliRun.h"

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
AliTRDdigitizer::AliTRDdigitizer():TNamed()
{
  //
  // AliTRDdigitizer default constructor
  //

  fInputFile      = NULL;
  fDigits         = NULL;
  fTRD            = NULL;
  fGeo            = NULL;
  fPRF            = NULL;
  fPRFsmp         = NULL;
  fTRFsmp         = NULL;

  fEvent          = 0;
  fGasGain        = 0.0;
  fNoise          = 0.0;
  fChipGain       = 0.0;
  fSinRange       = 0.0;
  fSoutRange      = 0.0;
  fADCoutRange    = 0.0;
  fADCinRange     = 0.0;
  fADCthreshold   = 0;
  fDiffusionOn    = 0;
  fDiffusionT     = 0.0;
  fDiffusionL     = 0.0;
  fElAttachOn     = 0;
  fElAttachProp   = 0.0;
  fExBOn          = 0;
  fOmegaTau       = 0.0;
  fPRFOn          = 0;
  fTRFOn          = 0;
  fDriftVelocity  = 0.0;
  fPadCoupling    = 0.0;
  fTimeCoupling   = 0.0;
  fTimeBinWidth   = 0.0;

  fPRFbin         = 0;
  fPRFlo          = 0.0;
  fPRFhi          = 0.0;
  fPRFwid         = 0.0;
  fPRFpad         = 0;
  fTRFbin         = 0;
  fTRFlo          = 0.0;
  fTRFhi          = 0.0;
  fTRFwid         = 0.0;

  fCompress       = kTRUE;
  fVerbose        = 1;
  fSDigits        = kFALSE;

}

//_____________________________________________________________________________
AliTRDdigitizer::AliTRDdigitizer(const Text_t *name, const Text_t *title)
                :TNamed(name,title)
{
  //
  // AliTRDdigitizer default constructor
  //

  fInputFile     = NULL;
  fDigits        = NULL;
  fTRD           = NULL;
  fGeo           = NULL;
  fPRF           = NULL;
  fPRFsmp        = NULL;
  fTRFsmp        = NULL;

  fEvent         = 0;

  fCompress      = kTRUE;
  fVerbose       = 1;
  fSDigits       = kFALSE;

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
  }

  if (fDigits) {
    delete fDigits;
  }

  if (fPRF) delete fPRF;

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

  ((AliTRDdigitizer &) d).fInputFile      = NULL;
  ((AliTRDdigitizer &) d).fDigits         = NULL;
  ((AliTRDdigitizer &) d).fTRD            = NULL;
  ((AliTRDdigitizer &) d).fGeo            = NULL;

  ((AliTRDdigitizer &) d).fEvent          = 0;

  ((AliTRDdigitizer &) d).fGasGain        = fGasGain;
  ((AliTRDdigitizer &) d).fNoise          = fNoise;
  ((AliTRDdigitizer &) d).fChipGain       = fChipGain;
  ((AliTRDdigitizer &) d).fSoutRange      = fSoutRange;
  ((AliTRDdigitizer &) d).fSinRange       = fSinRange;
  ((AliTRDdigitizer &) d).fADCoutRange    = fADCoutRange;
  ((AliTRDdigitizer &) d).fADCinRange     = fADCinRange;
  ((AliTRDdigitizer &) d).fADCthreshold   = fADCthreshold;
  ((AliTRDdigitizer &) d).fDiffusionOn    = fDiffusionOn; 
  ((AliTRDdigitizer &) d).fDiffusionT     = fDiffusionT;
  ((AliTRDdigitizer &) d).fDiffusionL     = fDiffusionL;
  ((AliTRDdigitizer &) d).fElAttachOn     = fElAttachOn;
  ((AliTRDdigitizer &) d).fElAttachProp   = fElAttachProp;
  ((AliTRDdigitizer &) d).fExBOn          = fExBOn;
  ((AliTRDdigitizer &) d).fOmegaTau       = fOmegaTau;
  ((AliTRDdigitizer &) d).fLorentzFactor  = fLorentzFactor;
  ((AliTRDdigitizer &) d).fDriftVelocity  = fDriftVelocity;
  ((AliTRDdigitizer &) d).fPadCoupling    = fPadCoupling;
  ((AliTRDdigitizer &) d).fTimeCoupling   = fTimeCoupling;
  ((AliTRDdigitizer &) d).fTimeBinWidth   = fTimeBinWidth;
  ((AliTRDdigitizer &) d).fPRFOn          = fPRFOn;
  ((AliTRDdigitizer &) d).fTRFOn          = fTRFOn;

  ((AliTRDdigitizer &) d).fCompress       = fCompress;
  ((AliTRDdigitizer &) d).fVerbose        = fVerbose;
  ((AliTRDdigitizer &) d).fSDigits        = fSDigits;

  fPRF->Copy(*((AliTRDdigitizer &) d).fPRF);

  ((AliTRDdigitizer &) d).fPRFbin         = fPRFbin;
  ((AliTRDdigitizer &) d).fPRFlo          = fPRFlo;
  ((AliTRDdigitizer &) d).fPRFhi          = fPRFhi;
  ((AliTRDdigitizer &) d).fPRFwid         = fPRFwid;
  ((AliTRDdigitizer &) d).fPRFpad         = fPRFpad;
  if (((AliTRDdigitizer &) d).fPRFsmp) delete ((AliTRDdigitizer &) d).fPRFsmp;
  ((AliTRDdigitizer &) d).fPRFsmp = new Float_t[fPRFbin];
  for (iBin = 0; iBin < fPRFbin; iBin++) {
    ((AliTRDdigitizer &) d).fPRFsmp[iBin] = fPRFsmp[iBin];
  }                                                                             
  ((AliTRDdigitizer &) d).fTRFbin         = fTRFbin;
  ((AliTRDdigitizer &) d).fTRFlo          = fTRFlo;
  ((AliTRDdigitizer &) d).fTRFhi          = fTRFhi;
  ((AliTRDdigitizer &) d).fTRFwid         = fTRFwid;
  if (((AliTRDdigitizer &) d).fTRFsmp) delete ((AliTRDdigitizer &) d).fTRFsmp;
  ((AliTRDdigitizer &) d).fTRFsmp = new Float_t[fTRFbin];
  for (iBin = 0; iBin < fTRFbin; iBin++) {
    ((AliTRDdigitizer &) d).fTRFsmp[iBin] = fTRFsmp[iBin];
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
Int_t AliTRDdigitizer::PadResponse(Float_t signal, Float_t dist, Float_t *pad)
{
  //
  // Applies the pad response
  //

  Int_t iBin =  ((Int_t) (( - dist - fPRFlo) / fPRFwid));

  Int_t iBin0 = iBin - fPRFpad;
  Int_t iBin1 = iBin;
  Int_t iBin2 = iBin + fPRFpad;

  if ((iBin0 >= 0) && (iBin2 < fPRFbin)) {

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
void AliTRDdigitizer::Init()
{
  //
  // Initializes the digitization procedure with standard values
  //

  // The default parameter for the digitization
  // Updated to new TRF 200 ns   
  fGasGain        = 1600.;
  fChipGain       = 8.0;
  fNoise          = 1000.;
  fADCoutRange    = 1023.;          // 10-bit ADC
  fADCinRange     = 1000.;          // 1V input range
  fADCthreshold   = 1;

  // For the summable digits
  fSinRange       = 1000000.;
  fSoutRange      = 1000000.;

  // Transverse and longitudinal diffusion coefficients (Xe/Isobutane)
  fDiffusionOn    = 1;
  fDiffusionT     = 0.060;
  fDiffusionL     = 0.017;

  // Propability for electron attachment
  fElAttachOn     = 0;
  fElAttachProp   = 0.0;

  // E x B effects
  fExBOn          = 0;
  // omega * tau.= arctan(Lorentz-angle)
  fOmegaTau       = 0.19438031;

  // The pad response function
  fPRFOn          =  1;
  fPRFlo          = -3.0;
  fPRFhi          =  3.0;
  fPRFbin         = 120;
  fPRFwid         = (fPRFhi - fPRFlo) / ((Float_t) fPRFbin);
  fPRFpad         = ((Int_t) (1.0 / fPRFwid));
  // New PRF from Bogdan  25/04/01
  fPRF            = new TF1("PRF"
                           ,"[0]*([1]+exp(-pow(sqrt(x*x),[3])/(2.0*[2])))"
                           ,fPRFlo,fPRFhi);
  fPRF->SetParameter(0, 0.8303); 
  fPRF->SetParameter(1,-0.00392); 
  fPRF->SetParameter(2, 0.472 * 0.472); 
  fPRF->SetParameter(3, 2.19); 

  // The time response function
  fTRFOn          =  1;

  // The drift velocity (cm / mus)
  fDriftVelocity  = 2.0;

  // The pad coupling factor (same number as for the TPC)
  fPadCoupling    = 0.5;

  // The time coupling factor (same number as for the TPC)
  fTimeCoupling   = 0.4;

}

//_____________________________________________________________________________
void AliTRDdigitizer::ReInit()
{
  //
  // Reinitializes the digitization procedure after a change in the parameter
  //

  if (!fGeo) {
    printf("AliTRDdigitizer::ReInit -- ");
    printf("No geometry defined. Run InitDetector() first\n");
    exit(1);
  }

  // Calculate the time bin width in ns
  fTimeBinWidth   = fGeo->GetTimeBinSize() / fDriftVelocity * 1000.0;

  // The range and the binwidth for the sampled TRF 
  fTRFbin = 100;
  // Start 0.8 mus before the signal
  fTRFlo  = -0.8 * fDriftVelocity;
  // End the maximum driftlength after the signal 
  fTRFhi  = AliTRDgeometry::DrThick() 
          + fGeo->GetTimeAfter() * fGeo->GetTimeBinSize();
  fTRFwid = (fTRFhi - fTRFlo) / ((Float_t) fTRFbin);

  // The Lorentz factor
  if (fExBOn) {
    fLorentzFactor = 1.0 / (1.0 + fOmegaTau*fOmegaTau);
  }
  else {
    fLorentzFactor = 1.0;
  }

}

//_____________________________________________________________________________
void AliTRDdigitizer::SampleTRF()
{
  //
  // Samples the time response function
  // It is defined according to Vasiles simulation of the preamp shaper
  // output and includes the effect of the ion tail (based on Tariqs 
  // Garfield simulation) and a shaping time of 125 ns FWHM
  //

  Int_t   ipos1;
  Int_t   ipos2;
  Float_t diff;

  const Float_t kShift = 0.0;
  const Float_t kScale = 0.5;
  //const Float_t kScale = 1.0;

  const Int_t kNpasa = 36;
  Float_t time[kNpasa]   = { -2.80,     -2.40,     -2.00,     -1.60
                           , -1.20,     -0.80,     -0.60,     -0.40
                           , -0.30,     -0.20,     -0.10,      0.00
                           ,  0.10,      0.20,      0.30,      0.40
                           ,  0.60,      0.80,      1.20,      1.60
                           ,  2.00,      2.40,      2.80,      3.20
                           ,  3.60,      4.00,      4.40,      4.80
                           ,  5.20,      5.60,      7.20,      9.20
                           , 11.20,     13.20,     15.20,     17.20     };
  Float_t signal[kNpasa] = {  0.000000,  0.000000,  0.000000,  0.000000
                           ,  0.000000,  0.000000,  0.015385,  0.086154
                           ,  0.236923,  0.452308,  0.726154,  1.003077
                           ,  0.953846,  0.652307,  0.332308,  0.181539
                           ,  0.120000,  0.083077,  0.049231,  0.024615
                           ,  0.015385,  0.009231,  0.003077,  0.000000
                           , -0.003077, -0.006154, -0.009231, -0.012308
                           , -0.015385, -0.018462, -0.018462, -0.018462
			   , -0.015385, -0.012308, -0.009231, -0.006154 };
  for (Int_t ipasa = 0; ipasa < kNpasa; ipasa++) {
    time[ipasa] = kScale * time[ipasa] + kShift; 
  }

  if (fTRFsmp) delete fTRFsmp;
  fTRFsmp = new Float_t[fTRFbin];

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
    if (ipos2 > kNpasa) ipos2 = kNpasa - 1;
    ipos1 = ipos2 - 1;

    fTRFsmp[iBin] = signal[ipos2] 
                  + diff * (signal[ipos2] - signal[ipos1]) 
                         / (  time[ipos2] -   time[ipos1]);

  }

}

//_____________________________________________________________________________
void AliTRDdigitizer::SamplePRF()
{
  //
  // Samples the pad response function
  //

  if (fPRFsmp) delete fPRFsmp;
  fPRFsmp = new Float_t[fPRFbin];
  for (Int_t iBin = 0; iBin < fPRFbin; iBin++) {
    Float_t bin = (((Float_t ) iBin) + 0.5) * fPRFwid + fPRFlo;
    fPRFsmp[iBin] = TMath::Max(fPRF->Eval(bin),0.0);
  }

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

  fEvent = nEvent;

  // Import the Trees for the event nEvent in the file
  Int_t nparticles = gAlice->GetEvent(fEvent);
  if (nparticles <= 0) {
    printf("AliTRDdigitizer::Open -- ");
    printf("No entries in the trees for event %d.\n",fEvent);
    return kFALSE;
  }

  return InitDetector();

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
    printf("AliTRDdigitizer::InitDetector -- ");
    printf("TRD must be version 1 (slow simulator).\n");
    exit(1);
  }

  // Get the geometry
  fGeo = fTRD->GetGeometry();
  printf("AliTRDdigitizer::InitDetector -- ");
  printf("Geometry version %d\n",fGeo->IsVersion());

  ReInit();

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::SumSDigits()
{
  //
  // Sums up the summable digits and creates final digits
  // Not yet implemented
  //

  return kFALSE;

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

  // Create a digits manager
  fDigits = new AliTRDdigitsManager();

  // Create a container for the amplitudes
  AliTRDsegmentArray *signalsArray 
                     = new AliTRDsegmentArray("AliTRDdataArrayF",AliTRDgeometry::Ndet());

  if (fTRFOn) {
    timeTRDbeg = ((Int_t) (-fTRFlo / fGeo->GetTimeBinSize())) - 1;
    timeTRDend = ((Int_t) ( fTRFhi / fGeo->GetTimeBinSize())) - 1;
    printf("AliTRDdigitizer::MakeDigits -- ");
    printf("Sample the TRF between -%d and %d\n",timeTRDbeg,timeTRDend);
  }

  Float_t elAttachProp = fElAttachProp / 100.; 

  // Create the sampled PRF
  SamplePRF();

  // Create the sampled TRF
  SampleTRF();

  if (!fGeo) {
    printf("AliTRDdigitizer::MakeDigits -- ");
    printf("No geometry defined\n");
    return kFALSE;
  }

  printf("AliTRDdigitizer::MakeDigits -- ");
  printf("Start creating digits.\n");
  if (fVerbose > 0) this->Dump();

  // Get the pointer to the hit tree
  TTree *HitTree = gAlice->TreeH();

  // Get the number of entries in the hit tree
  // (Number of primary particles creating a hit somewhere)
  Int_t nTrack = (Int_t) HitTree->GetEntries();
  if (fVerbose > 0) {
    printf("AliTRDdigitizer::MakeDigits -- ");
    printf("Found %d primary particles\n",nTrack);
  } 

  Int_t detectorOld = -1;
  Int_t countHits   =  0; 

  // Loop through all entries in the tree
  for (Int_t iTrack = 0; iTrack < nTrack; iTrack++) {

    gAlice->ResetHits();
    nBytes += HitTree->GetEvent(iTrack);

    // Get the number of hits in the TRD created by this particle
    Int_t nHit = fTRD->Hits()->GetEntriesFast();
    if (fVerbose > 0) {
      printf("AliTRDdigitizer::MakeDigits -- ");
      printf("Found %d hits for primary particle %d\n",nHit,iTrack);
    }

    // Loop through the TRD hits  
    for (Int_t iHit = 0; iHit < nHit; iHit++) {

      countHits++;

      AliTRDhit *hit = (AliTRDhit *) fTRD->Hits()->UncheckedAt(iHit);
              pos[0]   = hit->X();
              pos[1]   = hit->Y();
              pos[2]   = hit->Z();
      Float_t q        = hit->GetCharge();
      Int_t   track    = hit->Track();
      Int_t   detector = hit->GetDetector();
      Int_t   plane    = fGeo->GetPlane(detector);
      Int_t   sector   = fGeo->GetSector(detector);
      Int_t   chamber  = fGeo->GetChamber(detector);

      if (!(CheckDetector(plane,chamber,sector))) continue;

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

      if (fVerbose > 1) {
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
       
      // Don't analyze test hits
      if (hit->FromTest()) continue;

      if (detector != detectorOld) {

        if (fVerbose > 1) {
          printf("AliTRDdigitizer::MakeDigits -- ");
          printf("Get new container. New det = %d, Old det = %d\n"
                ,detector,detectorOld);
	}
        // Compress the old one if enabled
        if ((fCompress) && (detectorOld > -1)) {
          if (fVerbose > 1) {
            printf("AliTRDdigitizer::MakeDigits -- ");
            printf("Compress the old container ...");
	  }
          signals->Compress(1,0);
          for (iDict = 0; iDict < kNDict; iDict++) {
            dictionary[iDict]->Compress(1,0);
	  }
          if (fVerbose > 1) printf("done\n");
	}
	// Get the new container
        signals = (AliTRDdataArrayF *) signalsArray->At(detector);
        if (signals->GetNtime() == 0) {
          // Allocate a new one if not yet existing
          if (fVerbose > 1) {
            printf("AliTRDdigitizer::MakeDigits -- ");
            printf("Allocate a new container ... ");
	  }
          signals->Allocate(nRowMax,nColMax,nTimeTotal);
	}
        else {
	  // Expand an existing one
          if (fCompress) {
            if (fVerbose > 1) {
              printf("AliTRDdigitizer::MakeDigits -- ");
              printf("Expand an existing container ... ");
	    }
            signals->Expand();
	  }
	}
	// The same for the dictionary
        for (iDict = 0; iDict < kNDict; iDict++) {       
          dictionary[iDict] = fDigits->GetDictionary(detector,iDict);
          if (dictionary[iDict]->GetNtime() == 0) {
            dictionary[iDict]->Allocate(nRowMax,nColMax,nTimeTotal);
	  }
          else {
            if (fCompress) dictionary[iDict]->Expand();
	  }
        }      
        if (fVerbose > 1) printf("done\n");
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
        Int_t  rowE = ((Int_t) ((xyz[2] -  row0) * divideRow));
        if ((rowE < 0) || (rowE >= nRowMax)) continue;   

        // The pad column (rphi-direction)
        Int_t  colE = ((Int_t) ((xyz[1] -  col0) * divideCol));
        if ((colE < 0) || (colE >= nColMax)) continue;   
  
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
          Float_t dist = (xyz[1] - col0 - (colE + 0.5) * colPadSize) 
                       * divideCol;
          if (!(PadResponse(signal,dist,padSignal))) continue;
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
          if (fTRFOn) {
            Float_t time = (iTimeBin - timeE) * timeBinSize + timeOffset;
            timeResponse = TimeResponse(time);
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
	    // time bins. This has to be subtracted lateron.
            Int_t iCurrentTimeBin = iTimeBin + nTimeBefore;
            signalOld[iPad]  = signals->GetDataUnchecked(rowE,colPos,iCurrentTimeBin);
            signalOld[iPad] += padSignal[iPad] * timeResponse;
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

	  }

	}

      }

    }

  } // All hits finished

  printf("AliTRDdigitizer::MakeDigits -- ");
  printf("Finished analyzing %d hits\n",countHits);

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

    if (fVerbose > 0) {
      printf("AliTRDdigitizer::MakeDigits -- ");
      printf("Digitization for chamber %d\n",iDet);
    }

    // Add a container for the digits of this detector
    digits = fDigits->GetDigits(iDet);        
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
      dictionary[iDict] = fDigits->GetDictionary(iDet,iDict);
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
          for (iTime = 0; iTime < nTimeTotal; iTime++) {         

	    // Create summable digits
            if (fSDigits) {

              Float_t signalAmp = signals->GetDataUnchecked(iRow,iCol,iTime);
              Int_t adc  = 0;
              if (signalAmp >= fSinRange) {
                adc = ((Int_t) fSoutRange);
	      }
              else {
                adc = ((Int_t) (signalAmp * (fSoutRange / fSinRange)));
	      }
              nDigits++;
              digits->SetDataUnchecked(iRow,iCol,iTime,adc);

	    }
	    // Create normal digits
            else {

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

              // Store the amplitude of the digit if above threshold
              if (adc > fADCthreshold) {
                if (fVerbose > 2) {
                  printf("  iRow = %d, iCol = %d, iTime = %d\n"
                        ,iRow,iCol,iTime);
                  printf("  signal = %f, adc = %d\n",signalAmp,adc);
	        }
                nDigits++;
                digits->SetDataUnchecked(iRow,iCol,iTime,adc);
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
    printf("AliTRDdigitizer::MakeDigits -- ");
    printf("Found %d digits in detector %d (%3.0f).\n"
          ,nDigits,iDet
          ,100.0 * ((Float_t) nDigits) / nPixel);
 
    if (fCompress) signals->Compress(1,0);

  }

  printf("AliTRDdigitizer::MakeDigits -- ");
  printf("Total number of analyzed hits = %d\n",countHits);

  printf("AliTRDdigitizer::MakeDigits -- ");
  printf("Total digits data size = %d, %d, %d, %d\n",totalSizeDigits
                                                    ,totalSizeDict0
                                                    ,totalSizeDict1
                                                    ,totalSizeDict2);        

  return kTRUE;

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

  // Create the branches
  if (!(gAlice->TreeD()->GetBranch("TRDdigits"))) { 
    return kFALSE;
  }

  // Store the digits and the dictionary in the tree
  fDigits->WriteDigits();

  // Write the new tree into the input file (use overwrite option)
  Char_t treeName[15];
  sprintf(treeName,"TreeD%d",fEvent);
  printf("AliTRDdigitizer::WriteDigits -- ");
  printf("Write the digits tree %s for event %d.\n"
        ,treeName,fEvent);
  gAlice->TreeD()->Write(treeName,TObject::kOverwrite);
 
  return kTRUE;

}

//_____________________________________________________________________________
void AliTRDdigitizer::SetPRF(TF1 *prf)
{
  //
  // Defines a new pad response function
  //

  if (fPRF) delete fPRF;
  fPRF = prf;     

}

