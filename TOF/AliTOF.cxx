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
Revision 1.13  2000/05/10 16:52:18  vicinanz
New TOF version with holes for PHOS/RICH

Revision 1.11.2.1  2000/05/10 09:37:15  vicinanz
New version with Holes for PHOS/RICH

Revision 1.11  1999/11/05 22:39:06  fca
New hits structure

Revision 1.10  1999/11/01 20:41:57  fca
Added protections against using the wrong version of FRAME

Revision 1.9  1999/10/15 15:35:19  fca
New version for frame1099 with and without holes

Revision 1.9  1999/09/29 09:24:33  fca
Introduction of the Copyright and cvs Log

*/
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Of Flight                               FCA                         //
//  This class contains the basic functions for the Time Of Flight           //
//  detector. Functions specific to one particular geometry are              //
//  contained in the derived classes                                         //
//
//  VERSIONE WITH 5 SYMMETRIC MODULES ALONG Z AXIS
//  ==============================================
//  
//  VERSION WITH HOLES FOR PHOS AND TRD IN SPACEFRAME WITH HOLES
//
//  Volume sensibile : FPAD
//
//
//
// Begin_Html
/*
<img src="picts/AliTOFClass.gif">
*/
//End_Html
//             
//
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <iostream.h>

#include "AliTOF.h"
#include "AliTOFD.h"
#include "TBRIK.h"
#include "TNode.h"
#include "TObject.h"
#include "TRandom.h"
#include "AliRun.h"
#include "AliConst.h"

 
ClassImp(AliTOF)
 
//_____________________________________________________________________________
AliTOF::AliTOF()
{
  //
  // Default constructor
  //
  fIshunt   = 0;
}
 
//_____________________________________________________________________________
AliTOF::AliTOF(const char *name, const char *title)
       : AliDetector(name,title)
{
  //
  // AliTOF standard constructor
  // 
  // Here are fixed some important parameters
  //

  // Initialization of hits and digits array
  //
  fHits   = new TClonesArray("AliTOFhit",  405);
  gAlice->AddHitList(fHits);
  fIshunt  = 0;
  fDigits = new TClonesArray("AliTOFdigit",405);
  //
  // Digitization parameters
  //
  // (Transfer Functions to be inserted here)
  //
  SetMarkerColor(7);
  SetMarkerStyle(2);
  SetMarkerSize(0.4);

// General Geometrical Parameters
  fNTof    =  18;
  fRmax    = 399.0;//cm
  fRmin    = 370.0;//cm
  fZlenC   = 177.5;//cm
  fZlenB   = 141.0;//cm
  fZlenA   = 106.0;//cm
  fZtof    = 370.5;//cm

// Strip Parameters
  fStripLn = 122.0;//cm  Strip Length
  fSpace   =   5.5;//cm  Space Beetween the strip and the bottom of the plate 
  fDeadBndZ=   1.5;//cm  Dead Boundaries of a Strip along Z direction (width)
  fDeadBndX=   1.0;//cm  Dead Boundaries of a Strip along X direction (length)
  fXpad    =   2.5;//cm  X size of a pad
  fZpad    =   3.5;//cm  Z size of a pad
  fGapA    =   4.; //cm  Gap beetween tilted strip in A-type plate
  fGapB    =   6.; //cm  Gap beetween tilted strip in B-type plate
  fOverSpc =  15.3;//cm Space available for sensitive layers in radial direction
  fNpadX   =  48;  // Number of pads in a strip along the X direction
  fNpadZ   =   2;  // Number of pads in a strip along the Z direction
  fPadXStr = fNpadX*fNpadZ; //Number of pads per strip
  fNStripA = 0;
  fNStripB = 0;
  fNStripC = 0;
 
// Physical performances
  fTimeRes = 100.;//ps
  fChrgRes = 100.;//pC

// DAQ characteristics
  fPadXSector = 1932;
  fNRoc       = 14;
  fNFec       = 32;
  fNTdc       = 32;
  fNPadXRoc   = (Int_t)fPadXSector/fNRoc;
}

//_____________________________________________________________________________
void AliTOF::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add a TOF hit
  //
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliTOFhit(fIshunt, track, vol, hits);
}
 
//_____________________________________________________________________________
void AliTOF::AddDigit(Int_t *tracks, Int_t *vol, Float_t *digits)
{
  //
  // Add a TOF digit
  //
  TClonesArray &ldigits = *fDigits;
  new (ldigits[fNdigits++]) AliTOFdigit(tracks, vol, digits);
}


//_____________________________________________________________________________
void AliTOF::CreateGeometry()
{
  //
  // Common geometry code 
  //
  //Begin_Html
  /*
    <img src="picts/AliTOFv23.gif">
  */
  //End_Html
  //
  const Double_t kPi=TMath::Pi();
  const Double_t kDegrad=kPi/180.;
  //
  Float_t xTof, yTof, Wall;

  // frame inbetween TOF modules
  Wall = 4.;//cm

  // Sizes of TOF module with its support etc..
  xTof = 2.*(fRmin*TMath::Tan(10*kDegrad)-Wall/2-.5);
  yTof = fRmax-fRmin;

//  TOF module internal definitions 
  TOFpc(xTof, yTof, fZlenC, fZlenB, fZlenA, fZtof);
}

//_____________________________________________________________________________
void AliTOF::DrawModule()
{
  //
  // Draw a shaded view of the common part of the TOF geometry
  //

   cout << " Drawing of AliTOF"<< endl; 
  // Set everything unseen
  gMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible
  gMC->Gsatt("FTOA","SEEN",1);
  gMC->Gsatt("FTOB","SEEN",1);
  gMC->Gsatt("FTOC","SEEN",1);
  gMC->Gsatt("FLTA","SEEN",1);
  gMC->Gsatt("FLTB","SEEN",1);
  gMC->Gsatt("FLTC","SEEN",1);
  gMC->Gsatt("FSTR","SEEN",1);
  //
  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 0, 1000, -1000, 1000, -1000, 1000);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 40, 30, 0, 12, 9.5, .02, .02);
  gMC->Gdhead(1111, "Time Of Flight");
  gMC->Gdman(18, 4, "MAN");
  gMC->Gdopt("hide","off");
}

//_____________________________________________________________________________
void AliTOF::CreateMaterials()
{
  //
  // Defines TOF materials for all versions
  // Authors :   Maxim Martemianov, Boris Zagreev (ITEP) 
  //            18/09/98 
  //
  Int_t   ISXFLD = gAlice->Field()->Integ();
  Float_t SXMGMX = gAlice->Field()->Max();
  //
  //--- Quartz (SiO2) 
  Float_t   aq[2] = { 28.0855,15.9994 };
  Float_t   zq[2] = { 14.,8. };
  Float_t   wq[2] = { 1.,2. };
  Float_t   dq = 2.20;
  Int_t nq = -2;
  // --- Freon
  Float_t afre[2]  = {12.011,18.9984032 };
  Float_t zfre[2]  = { 6., 9.};
  Float_t wfre[2]  = { 5.,12.};
  Float_t densfre  = 1.5;
  Int_t nfre = -2;
  // --- CO2 
  Float_t ac[2]   = {12.,16.};
  Float_t zc[2]   = { 6., 8.};
  Float_t wc[2]   = { 1., 2.};
  Float_t dc = .001977;
  Int_t nc = -2;
   // For mylar (C5H4O2) 
  Float_t amy[3] = { 12., 1., 16. };
  Float_t zmy[3] = {  6., 1.,  8. };
  Float_t wmy[3] = {  5., 4.,  2. };
  Float_t dmy    = 1.39;
  Int_t nmy = -3;
 // For polyethilene (CH2) for honeycomb!!!!
  Float_t ape[2] = { 12., 1. };
  Float_t zpe[2] = {  6., 1. };
  Float_t wpe[2] = {  1., 2. };
  Float_t dpe    = 0.935*0.479; //To have 1%X0 for 1cm as for honeycomb
  Int_t npe = -2;
  // --- G10 
  Float_t ag10[4] = { 12.,1.,16.,28. };
  Float_t zg10[4] = {  6.,1., 8.,14. };
  Float_t wmatg10[4] = { .259,.288,.248,.205 };
  Float_t densg10  = 1.7;
  Int_t nlmatg10 = -4;
  // --- DME 
  Float_t adme[5] = { 12.,1.,16.,19.,79. };
  Float_t zdme[5] = {  6.,1., 8., 9.,35. };
  Float_t wmatdme[5] = { .4056,.0961,.2562,.1014,.1407 };
  Float_t densdme  = .00205;
  Int_t nlmatdme = 5;
  // ---- ALUMINA (AL203) 
  Float_t aal[2] = { 27.,16.};
  Float_t zal[2] = { 13., 8.};
  Float_t wmatal[2] = { 2.,3. };
  Float_t densal  = 2.3;
  Int_t nlmatal = -2;
  // -- Water
  Float_t awa[2] = {  1., 16. };
  Float_t zwa[2] = {  1.,  8. };
  Float_t wwa[2] = {  2.,  1. };
  Float_t dwa    = 1.0;
  Int_t nwa = -2;
  //
  //AliMaterial(0, "Vacuum$", 1e-16, 1e-16, 1e-16, 1e16, 1e16);
  AliMaterial( 1, "Air$",14.61,7.3,0.001205,30423.24,67500.);
  AliMaterial( 2, "Cu $",  63.54, 29.0, 8.96, 1.43, 14.8);
  AliMaterial( 3, "C  $",  12.01,  6.0, 2.265,18.8, 74.4);
  AliMixture ( 4, "Polyethilene$", ape, zpe, dpe, npe, wpe);
  AliMixture ( 5, "G10$", ag10, zg10, densg10, nlmatg10, wmatg10);
  AliMixture ( 6, "DME ", adme, zdme, densdme, nlmatdme, wmatdme);
  AliMixture ( 7, "CO2$", ac, zc, dc, nc, wc);
  AliMixture ( 8, "ALUMINA$", aal, zal, densal, nlmatal, wmatal);
  AliMaterial( 9, "Al $", 26.98, 13., 2.7, 8.9, 37.2);
  AliMaterial(10, "C-TRD$", 12.01, 6., 2.265*18.8/69.282*15./100, 18.8, 74.4); // for 15%
  AliMixture (11, "Mylar$",  amy, zmy, dmy, nmy, wmy);
  AliMixture (12, "Freon$",  afre, zfre, densfre, nfre, wfre);
  AliMixture (13, "Quartz$", aq, zq, dq, nq, wq);
  AliMixture (14, "Water$",  awa, zwa, dwa, nwa, wwa);

  Float_t epsil, stmin, deemax, stemax;
 
  //   Previous data
  //       EPSIL  =  0.1   ! Tracking precision, 
  //       STEMAX = 0.1      ! Maximum displacement for multiple scattering
  //       DEEMAX = 0.1    ! Maximum fractional energy loss, DLS 
  //       STMIN  = 0.1 
  //
  //   New data  
  epsil  = .001;  // Tracking precision,
  stemax = -1.;   // Maximum displacement for multiple scattering
  deemax = -.3;   // Maximum fractional energy loss, DLS
  stmin  = -.8;

  AliMedium( 1, "Air$"  ,  1, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium( 2, "Cu $"  ,  2, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium( 3, "C  $"  ,  3, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium( 4, "Pol$"  ,  4, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium( 5, "G10$"  ,  5, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium( 6, "DME$"  ,  6, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium( 7, "CO2$"  ,  7, 0, ISXFLD, SXMGMX, 10., -.01, -.1, .01, -.01);
  AliMedium( 8,"ALUMINA$", 8, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium( 9,"Al Frame$",9, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium(10, "DME-S$",  6, 1, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium(11, "C-TRD$", 10, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium(12, "Myl$"  , 11, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium(13, "Fre$"  , 12, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium(14, "Fre-S$", 12, 1, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium(15, "Glass$", 13, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium(16, "Water$", 14, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
}

//_____________________________________________________________________________
Int_t AliTOF::DistancetoPrimitive(Int_t , Int_t )
{
  //
  // Returns distance from mouse pointer to detector, default version
  //
  return 9999;
}
 
//_____________________________________________________________________________
void AliTOF::Init()
{
  //
  // Initialise TOF detector after it has been built
  //
  // Set id of TOF sensitive volume
  if (IsVersion() !=0) fIdSens=gMC->VolId("FPAD");
  //
}

//____________________________________________________________________________
void AliTOF::MakeBranch(Option_t* option)
//
// Initializes the Branches of the TOF inside the 
// trees written for each event. 
//  AliDetector::MakeBranch initializes just the 
// Branch inside TreeH. Here we add the branch in 
// TreeD.
//
{

  AliDetector::MakeBranch(option);

  Int_t buffersize = 4000;
  Char_t branchname[10];
  sprintf(branchname,"%s",GetName());
  char *D = strstr(option,"D");

  if (fDigits && gAlice->TreeD() && D){
     gAlice->TreeD()->Branch(branchname,&fDigits,buffersize);
     printf("Making Branch %s for digits \n",branchname);
  }
}

//____________________________________________________________________________
void AliTOF::FinishEvent()
{
//  Hits2Digits();
}


//____________________________________________________________________________
void AliTOF::Hits2Digits(Int_t evNumber)
//
//  Starting from the Hits Tree (TreeH), this
// function writes the Digits Tree (TreeD) storing 
// the digits informations.
// Has to be called just at the end of an event or 
// at the end of a whole run.
//  It could  also be called by AliTOF::Finish Event()
// but it can be too heavy.
// Just for MC events. 
//
// Called by the macro H2D.C
//
{
  AliTOFhit* currentHit;
  TTree    *TD, *TH;
  Int_t    tracks[3];
  Int_t    vol[5];
  Float_t  digit[2];
  TClonesArray* TOFhits=this->Hits();

  Int_t nparticles = gAlice->GetEvent(evNumber); 
  if (nparticles <= 0) return;

  TD = gAlice->TreeD();
  TH = gAlice->TreeH();
  Int_t    ntracks =(Int_t) TH->GetEntries();
  Int_t    nbytes, nhits;
  TRandom *rnd = new TRandom();

  for (Int_t ntk=0; ntk<ntracks; ntk++){

     nbytes = TH->GetEvent(ntk);
     nhits  = TOFhits->GetEntriesFast();

     for (Int_t hit=0; hit<nhits; hit++){

        currentHit = (AliTOFhit*)(TOFhits->At(hit));

        vol[0] = currentHit->GetSector();
        vol[1] = currentHit->GetPlate();
        vol[2] = currentHit->GetPad_x();
        vol[3] = currentHit->GetPad_z();
        vol[4] = currentHit->GetStrip();

        Float_t IdealTime = currentHit->GetTof();
        Float_t TDCTime   = rnd->Gaus(IdealTime, fTimeRes);	
        digit[0] = TDCTime;

        Float_t IdealCharge = currentHit->GetEdep();
        Float_t ADCcharge = rnd->Gaus(IdealCharge, fChrgRes);
        digit[1] = ADCcharge;

        Int_t Track = currentHit -> GetTrack();
        tracks[0] = Track;
	tracks[1] = 0;
	tracks[2] = 0;

        Bool_t Overlap = CheckOverlap(vol, digit, Track);
        if(!Overlap) AddDigit(tracks, vol, digit);
     }
  }
  TD->Fill();
  TD->Write();  
}

//___________________________________________________________________________
Bool_t AliTOF::CheckOverlap(Int_t* vol, Float_t* digit,Int_t Track)
//
// Checks if 2 or more hits belong to the same pad.
// In this case the data assigned to the digit object
// are the ones of the first hit in order of Time.
//
// Called only by Hits2Digits.
//
{
        Bool_t Overlap = 0;
        Int_t  vol2[5];

        for (Int_t ndig=0; ndig<fNdigits; ndig++){
	   AliTOFdigit* currentDigit = (AliTOFdigit*)(fDigits->UncheckedAt(ndig));
           currentDigit->GetLocation(vol2);
           Bool_t Idem=1;
           for (Int_t i=0;i<=4;i++){
               if (vol[i]!=vol2[i]) Idem=0;}
           if (Idem){
 	      Float_t TDC2 = digit[0];
              Float_t TDC1 = currentDigit->GetTdc();
              if (TDC1>TDC2){
                  currentDigit->SetTdc(TDC2); 
                  currentDigit->SetAdc(digit[1]);
              }
              currentDigit->AddTrack(Track);
              Overlap = 1;
           }
        }
        return Overlap;
}


//____________________________________________________________________________
void AliTOF::Digits2Raw(Int_t evNumber)
//
// Starting from digits, writes the 
// Raw Data objects, i.e. a 
// TClonesArray of 18 AliTOFRawSector objects
//

{
  TTree* TD;

  Int_t nparticles = gAlice->GetEvent(evNumber); 
  if (nparticles <= 0) return;

  TD = gAlice->TreeD();
  
  TClonesArray* TOFdigits = this->Digits();
  Int_t ndigits = TOFdigits->GetEntriesFast();

  TClonesArray* Raw = new TClonesArray("AliTOFRawSector",fNTof+2); 

  for (Int_t isect=1;isect<=fNTof;isect++){
     AliTOFRawSector* currentSector = (AliTOFRawSector*)Raw->UncheckedAt(isect);
     TClonesArray* RocData = (TClonesArray*)currentSector->GetRocData();

     for (Int_t digit=0; digit<ndigits; digit++){
        AliTOFdigit* currentDigit = (AliTOFdigit*)TOFdigits->UncheckedAt(digit);
        Int_t sector = currentDigit->GetSector();
        if (sector==isect){
	    Int_t   Pad    = currentDigit -> GetTotPad();
	    Int_t   Roc    = (Int_t)(Pad/fNPadXRoc)-1;
	    if (Roc>=fNRoc) printf("Wrong n. of ROC ! Roc = %i",Roc);
            Int_t   PadRoc = (Int_t) Pad%fNPadXRoc;
	    Int_t   Fec    = (Int_t)(PadRoc/fNFec)-1;
            Int_t   Tdc    = (Int_t)(PadRoc%fNFec)-1;
            Float_t Time   = currentDigit->GetTdc();
            Float_t Charge = currentDigit->GetAdc();
	    AliTOFRoc* currentROC = (AliTOFRoc*)RocData->UncheckedAt(Roc);
	    Int_t Error    = 0;
            currentROC->AddItem(Fec, Tdc, Error, Charge, Time);
	}
     }
     
     UInt_t TotSize=16,RocSize=0;
     UInt_t RocHead[14],RocChek[14];
     UInt_t GlobalCheckSum=0;

     for (UInt_t iRoc = 1; iRoc<(UInt_t)fNRoc; iRoc++){
        AliTOFRoc* currentRoc = (AliTOFRoc*)RocData->UncheckedAt(iRoc); 
	RocSize  = currentRoc->Items*2+1;
	TotSize += RocSize*4;
	if (RocSize>=pow(2,16)) RocSize=0;
	RocHead[iRoc]   = iRoc<<28;
	RocHead[iRoc]  += RocSize;
	RocChek[iRoc]   = currentRoc->GetCheckSum();
        Int_t HeadCheck = currentRoc->BitCount(RocHead[iRoc]);
	GlobalCheckSum += HeadCheck;
	GlobalCheckSum += RocChek[iRoc];
     }
     
     AliTOFRoc* DummyRoc = new AliTOFRoc();
     TotSize *= 4;
     if (TotSize>=pow(2,24)) TotSize=0;
     UInt_t Header = TotSize;
     UInt_t SectId = ((UInt_t)isect)<<24;
     Header += SectId;
     GlobalCheckSum += DummyRoc->BitCount(Header);
     currentSector->SetGlobalCS(GlobalCheckSum);
     currentSector->SetHeader(Header);
  }  
}
 
//____________________________________________________________________________
void AliTOF::Raw2Digits(Int_t evNumber)
//
//  Converts Raw Data objects into digits objects.
//  We schematize the raw data with a 
//  TClonesArray of 18 AliTOFRawSector objects
//
{
  TTree    *TD;
  Int_t    vol[5];
  Int_t    tracks[3];
  Float_t  digit[2];
 
  tracks[0]=0;
  tracks[1]=0;
  tracks[2]=0;
 
  Int_t nparticles = gAlice->GetEvent(evNumber); 
  if (nparticles <= 0) return;

  TD = gAlice->TreeD();
  
  TClonesArray* Raw = new TClonesArray("AliTOFRawSector",fNTof+2);
  
  for(Int_t nSec=1; nSec<=fNTof; nSec++){
     AliTOFRawSector* currentSector = (AliTOFRawSector*)Raw->UncheckedAt(nSec);
     TClonesArray* RocData = (TClonesArray*)currentSector->GetRocData();
     for(Int_t nRoc=1; nRoc<=14; nRoc++){
        AliTOFRoc* currentRoc = (AliTOFRoc*)RocData->UncheckedAt(nRoc);
        Int_t currentItems = currentRoc->GetItems();
        for(Int_t item=1; item<currentItems; item++){ 
           Int_t nPad = currentRoc->GetTotPad(item);        
	   vol[0] = nSec;
	   Int_t nStrip = (Int_t)(nPad/fPadXStr)+1;
	   Int_t nPlate = 5;
	   if (nStrip<=fNStripC+2*fNStripB+fNStripA) nPlate = 4;
	   if (nStrip<=fNStripC+fNStripB+fNStripA)   nPlate = 3;
	   if (nStrip<=fNStripC+fNStripB)            nPlate = 2;
	   if (nStrip<=fNStripC)                     nPlate=1;
	   vol[1] = nPlate;
	   switch (nPlate){
	   case 1: break;
	   case 2: nStrip -= (fNStripC);
	           break;
	   case 3: nStrip -= (fNStripC+fNStripB);
	           break;
	   case 4: nStrip -= (fNStripC+fNStripB+fNStripA);
	           break;
	   case 5: nStrip -= (fNStripC+2*fNStripB+fNStripA);
	           break;
	   }
           vol[2] = nStrip;
           Int_t Pad = nPad%fPadXStr;
	   if (Pad==0) Pad=fPadXStr;
	   Int_t nPadX=0, nPadZ=0;
	   (Pad>fNpadX)? nPadX -= fNpadX : nPadX = Pad ;
	   vol[3] = nPadX;
	   (Pad>fNpadX)? nPadZ = 2 : nPadZ = 1 ;
	   vol[4] = nPadZ;
	   UInt_t error=0;
	   Float_t TDC = currentRoc->GetTime(item,error);
	   if (!error) digit[0]=TDC;
	   digit[1] = currentRoc->GetCharge(item);
	   AddDigit(tracks,vol,digit);
        }
     }
  }
  TD->Fill();
  TD->Write();
} 


/******************************************************************************/

ClassImp(AliTOFhit)
 
//______________________________________________________________________________
AliTOFhit::AliTOFhit(Int_t shunt, Int_t track, Int_t *vol,
                     Float_t *hits)
:AliHit(shunt, track)
//
// Constructor of hit object
//
{
  //
  // Store a TOF hit
  // _______________
  //
  // Hit Volume
  // 
  fSector= vol[0];
  fPlate = vol[1];
  fStrip = vol[2];
  fPad_x = vol[3];
  fPad_z = vol[4];
  //
  //Position
  fX = hits[0];
  fY = hits[1];
  fZ = hits[2];
  //
  // Momentum
  fPx  = hits[3];
  fPy  = hits[4];
  fPz  = hits[5];
  fPmom= hits[6];
  //
  // Time Of Flight
  fTof = hits[7];   //TOF[s]
  //
  // Other Data
  fDx  = hits[8];   //Distance from the edge along x axis
  fDy  = hits[9];   //Y cohordinate of the hit
  fDz  = hits[10];  //Distance from the edge along z axis
  fIncA= hits[11];  //Incidence angle
  fEdep= hits[12];  //Energy loss in TOF pad
}

//******************************************************************************

ClassImp(AliTOFdigit)

//______________________________________________________________________________
AliTOFdigit::AliTOFdigit(Int_t *tracks, Int_t *vol,Float_t *digit)
:AliDigit(tracks)
//
// Constructor of digit object
//
{
  fSector = vol[0];
  fPlate  = vol[1];
  fStrip  = vol[2];
  fPad_x  = vol[3];
  fPad_z  = vol[4];
  fTdc    = digit[0];
  fAdc    = digit[1];
}

//______________________________________________________________________________
void AliTOFdigit::GetLocation(Int_t *Loc)
//
// Get the cohordinates of the digit
// in terms of Sector - Plate - Strip - Pad
//
{
   Loc[0]=fSector;
   Loc[1]=fPlate;
   Loc[2]=fStrip;
   Loc[3]=fPad_x;
   Loc[4]=fPad_z;
}

//______________________________________________________________________________
Int_t AliTOFdigit::GetTotPad()
//
// Get the "total" index of the pad inside a Sector
// starting from the digits data.
//
{
  AliTOF* TOF;
  
  if(gAlice){
     TOF =(AliTOF*) gAlice->GetDetector("TOF");
  }else{
     printf("AliTOFdigit::GetTotPad - No AliRun object present, exiting");
     return 0;
  }
  
  Int_t Pad = fPad_x+TOF->fNpadX*(fPad_z-1);
  Int_t Before=0;

  switch(fPlate){ 
  case 1: Before = 0;
          break;
  case 2: Before = TOF->fNStripC;
          break;
  case 3: Before = TOF->fNStripB + TOF->fNStripC;
          break;
  case 4: Before = TOF->fNStripA + TOF->fNStripB + TOF->fNStripC;
          break;
  case 5: Before = TOF->fNStripA + 2*TOF->fNStripB + TOF->fNStripC;
          break;
  }
  
  Int_t Strip = fStrip+Before;
  Int_t PadTot = TOF->fPadXStr*(Strip-1)+Pad;
  return PadTot;
}

//______________________________________________________________________________
void AliTOFdigit::AddTrack(Int_t track)
//
// Add a track to the digit 
//
{
  if (fTracks[1]==0){
     fTracks[1] = track;
  }else if (fTracks[2]==0){
     fTracks[2] = track;
  }else{
     printf("AliTOFdigit::AddTrack ERROR: Too many Tracks (>3) \n");
  }
}

