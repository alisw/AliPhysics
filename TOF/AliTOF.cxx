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
Revision 1.20  2001/05/04 10:09:47  vicinanz
Major upgrades to the strip structure

Revision 1.19  2001/03/12 17:47:25  hristov
Changes needed on Sun with CC 5.0

Revision 1.18  2001/01/26 19:57:42  hristov
Major upgrade of AliRoot code

Revision 1.17  2000/10/19 09:58:14  vicinanz
Updated Hits2Digit procedure

Revision 1.16  2000/10/02 21:28:17  fca
Removal of useless dependecies via forward declarations

Revision 1.15  2000/05/18 14:33:01  vicinanz
Modified to be full HP compliant

Revision 1.14  2000/05/15 19:32:36  fca
Add AddHitList !!

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
#include "TTree.h"
#include "TFile.h"

#include "AliRun.h"
#include "AliMC.h"
#include "AliMagF.h"
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
  fNTof    =  18;  // number of sectors
  fRmax    = 399.0;//cm 
  fRmin    = 370.0;//cm
  fZlenC   = 177.5;//cm length of module C
  fZlenB   = 141.0;//cm length of module B
  fZlenA   = 106.0;//cm length of module A
  fZtof    = 370.5;//cm total semi-length of TOF detector

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
  fNStripA = 15; // number of strips in A type module 
  fNStripB = 19; // number of strips in B type module
  fNStripC = 20; // number of strips in C type module
 
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
  Float_t xTof, yTof, wall;

  // frame inbetween TOF modules
  wall = 4.;//cm

  // Sizes of TOF module with its support etc..
  xTof = 2.*(fRmin*TMath::Tan(10*kDegrad)-wall/2-.5);
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
  // Revision: F. Pierella 5-3-2001
  // Bologna University
  //
  Int_t   isxfld = gAlice->Field()->Integ();
  Float_t sxmgmx = gAlice->Field()->Max();
  //
  //--- Quartz (SiO2) 
  Float_t   aq[2] = { 28.0855,15.9994 };
  Float_t   zq[2] = { 14.,8. };
  Float_t   wq[2] = { 1.,2. };
  Float_t   dq = 2.20;
  Int_t nq = -2;
  // --- Freon C2F4H2 (TOF-TDR pagg.)
  // Geant Manual CONS110-1, pag. 43 (Geant, Detector Description and Simulation Tool)
  Float_t afre[3]  = {12.011,18.998,1.007};
  Float_t zfre[3]  = { 6., 9., 1.}; 
  Float_t wfre[3]  = { 2., 4., 2.};
  Float_t densfre  = 0.00375;   
// http://www.fi.infn.it/sezione/prevprot/gas/freon.html
  Int_t nfre = -3; 
/*
  //-- Isobutane quencher C4H10 (5% in the sensitive mixture)
  Float_t aiso[2]  = {12.011,1.007};
  Float_t ziso[2]  = { 6.,  1.};
  Float_t wiso[2]  = { 4., 10.};
  Float_t densiso  = .......;  // (g/cm3) density
  Int_t nfre = -2; // < 0 i.e. proportion by number of atoms of each kind
  //-- SF6 (5% in the sensitive mixture)
  Float_t asf[3]  = {32.066,18.998};
  Float_t zsf[3]  = { 16., 9.};
  Float_t wsf[3]  = {  1., 6.}; 
  Float_t denssf  = .....;   // (g/cm3) density
  Int_t nfre = -2; // < 0 i.e. proportion by number of atoms of each kind
*/
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
 // For polyethilene (CH2) - honeycomb -
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

  AliMedium( 1, "Air$"  ,  1, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 2, "Cu $"  ,  2, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 3, "C  $"  ,  3, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 4, "Pol$"  ,  4, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 5, "G10$"  ,  5, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 6, "DME$"  ,  6, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 7, "CO2$"  ,  7, 0, isxfld, sxmgmx, 10., -.01, -.1, .01, -.01);
  AliMedium( 8,"ALUMINA$", 8, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 9,"Al Frame$",9, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(10, "DME-S$",  6, 1, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(11, "C-TRD$", 10, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(12, "Myl$"  , 11, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(13, "Fre$"  , 12, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(14, "Fre-S$", 12, 1, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(15, "Glass$", 13, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(16, "Water$", 14, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
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
void AliTOF::MakeBranch(Option_t* option, const char *file)
{
 //
 // Initializes the Branches of the TOF inside the 
 // trees written for each event. 
 //  AliDetector::MakeBranch initializes just the 
 // Branch inside TreeH. Here we add the branch in 
 // TreeD.
 //

  AliDetector::MakeBranch(option,file);

  Int_t buffersize = 4000;
  Char_t branchname[10];
  sprintf(branchname,"%s",GetName());
  
  const char *D = strstr(option,"D");

  if (fDigits && gAlice->TreeD() && D){
     MakeBranchInTree(gAlice->TreeD(), 
                              branchname, &fDigits,buffersize, file) ;
  }
}

//____________________________________________________________________________
void AliTOF::FinishEvent()
{
//  Hits2Digits();
}

//___________________________________________
void AliTOF::SDigits2Digits()
{
//
// Genneratedigits
//
    int nparticles = gAlice->GetNtrack();
    cout << "Particles       :" <<nparticles<<endl;
    if (nparticles > 0 ) {
      Hits2Digits(0);
    }
}

//____________________________________________________________________________
void AliTOF::Hits2Digits(Int_t evNumber)
{
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

  AliTOFhit* currentHit;
  TTree    *tD, *tH;
  Int_t    tracks[3];
  Int_t    vol[5];
  Float_t  digit[2];
  TClonesArray* tofhits=this->Hits();

  Int_t nparticles =  gAlice->GetNtrack();
  if (nparticles <= 0) return;

  tD = gAlice->TreeD();
  tH = gAlice->TreeH();
  Int_t    ntracks =(Int_t) tH->GetEntries();
  Int_t    nbytes, nhits;
  TRandom* rnd = new TRandom();

  for (Int_t ntk=0; ntk<ntracks; ntk++){

     nbytes = tH->GetEvent(ntk);
     nhits  = tofhits->GetEntriesFast();

     for (Int_t hit=0; hit<nhits; hit++){

        currentHit = (AliTOFhit*)(tofhits->At(hit));

        vol[0] = currentHit->GetSector();
        vol[1] = currentHit->GetPlate();
        vol[2] = currentHit->GetPadx();
        vol[3] = currentHit->GetPadz();
        vol[4] = currentHit->GetStrip();

        Float_t idealtime = currentHit->GetTof();
        Float_t tdctime   = rnd->Gaus(idealtime, fTimeRes);	
        digit[0] = tdctime;

        Float_t idealcharge = currentHit->GetEdep();
        Float_t adccharge = rnd->Gaus(idealcharge, fChrgRes);
        digit[1] = adccharge;

        Int_t tracknum = currentHit -> GetTrack();
        tracks[0] = tracknum;
	tracks[1] = 0;
	tracks[2] = 0;

        Bool_t overlap = CheckOverlap(vol, digit, tracknum);
        if(!overlap) AddDigit(tracks, vol, digit);
     }
  }
  delete rnd;
  rnd = 0;
  tD->Fill();
  tD->Write();  
}

//___________________________________________________________________________
Bool_t AliTOF::CheckOverlap(Int_t* vol, Float_t* digit,Int_t Track)
{
//
// Checks if 2 or more hits belong to the same pad.
// In this case the data assigned to the digit object
// are the ones of the first hit in order of Time.
//
// Called only by Hits2Digits.
//

        Bool_t overlap = 0;
        Int_t  vol2[5];

        for (Int_t ndig=0; ndig<fNdigits; ndig++){
	   AliTOFdigit* currentDigit = (AliTOFdigit*)(fDigits->UncheckedAt(ndig));
           currentDigit->GetLocation(vol2);
           Bool_t idem=1;
           for (Int_t i=0;i<=4;i++){
               if (vol[i]!=vol2[i]) idem=0;}
           if (idem){
 	      Float_t tdc2 = digit[0];
              Float_t tdc1 = currentDigit->GetTdc();
              if (tdc1>tdc2){
                  currentDigit->SetTdc(tdc2); 
                  currentDigit->SetAdc(digit[1]);
              }
              currentDigit->AddTrack(Track);
              overlap = 1;
           }
        }
        return overlap;
}


//____________________________________________________________________________
void AliTOF::Digits2Raw(Int_t evNumber)
{
//
// Starting from digits, writes the 
// Raw Data objects, i.e. a 
// TClonesArray of 18 AliTOFRawSector objects
//

  TTree* tD;

  Int_t nparticles = gAlice->GetEvent(evNumber); 
  if (nparticles <= 0) return;

  tD = gAlice->TreeD();
  
  TClonesArray* tofdigits = this->Digits();
  Int_t ndigits = tofdigits->GetEntriesFast();

  TClonesArray* rawsectors = new TClonesArray("AliTOFRawSector",fNTof+2); 

  for (Int_t isect=1;isect<=fNTof;isect++){
     AliTOFRawSector* currentSector = (AliTOFRawSector*)rawsectors->UncheckedAt(isect);
     TClonesArray* rocData = (TClonesArray*)currentSector->GetRocData();

     for (Int_t digit=0; digit<ndigits; digit++){
        AliTOFdigit* currentDigit = (AliTOFdigit*)tofdigits->UncheckedAt(digit);
        Int_t sector = currentDigit->GetSector();
        if (sector==isect){
	    Int_t   pad    = currentDigit -> GetTotPad();
	    Int_t   roc    = (Int_t)(pad/fNPadXRoc)-1;
	    if (roc>=fNRoc) printf("Wrong n. of ROC ! Roc = %i",roc);
            Int_t   padRoc = (Int_t) pad%fNPadXRoc;
	    Int_t   fec    = (Int_t)(padRoc/fNFec)-1;
            Int_t   tdc    = (Int_t)(padRoc%fNFec)-1;
            Float_t time   = currentDigit->GetTdc();
            Float_t charge = currentDigit->GetAdc();
	    AliTOFRoc* currentROC = (AliTOFRoc*)rocData->UncheckedAt(roc);
	    Int_t error    = 0;
            currentROC->AddItem(fec, tdc, error, charge, time);
	}
     }
     
     UInt_t totSize=16,rocSize=0;
     UInt_t rocHead[14],rocChek[14];
     UInt_t globalCheckSum=0;

     for (UInt_t iRoc = 1; iRoc<(UInt_t)fNRoc; iRoc++){
        AliTOFRoc* currentRoc = (AliTOFRoc*)rocData->UncheckedAt(iRoc); 
	rocSize  = currentRoc->GetItems()*2+1;
	totSize += rocSize*4;
	if (rocSize>=TMath::Power(2,16)) rocSize=0;
	rocHead[iRoc]   = iRoc<<28;
	rocHead[iRoc]  += rocSize;
	rocChek[iRoc]   = currentRoc->GetCheckSum();
        Int_t headCheck = currentRoc->BitCount(rocHead[iRoc]);
	globalCheckSum += headCheck;
	globalCheckSum += rocChek[iRoc];
     }
     
     AliTOFRoc* dummyRoc = new AliTOFRoc();
     totSize *= 4;
     if (totSize>=TMath::Power(2,24)) totSize=0;
     UInt_t header = totSize;
     UInt_t sectId = ((UInt_t)isect)<<24;
     header += sectId;
     globalCheckSum += dummyRoc->BitCount(header);
     currentSector->SetGlobalCS(globalCheckSum);
     currentSector->SetHeader(header);
  }  
}
 
//____________________________________________________________________________
void AliTOF::Raw2Digits(Int_t evNumber)
{
//
//  Converts Raw Data objects into digits objects.
//  We schematize the raw data with a 
//  TClonesArray of 18 AliTOFRawSector objects
//

  TTree    *tD;
  Int_t    vol[5];
  Int_t    tracks[3];
  Float_t  digit[2];
 
  tracks[0]=0;
  tracks[1]=0;
  tracks[2]=0;
 
  Int_t nparticles = gAlice->GetEvent(evNumber); 
  if (nparticles <= 0) return;

  tD = gAlice->TreeD();
  
  TClonesArray* rawsectors = new TClonesArray("AliTOFRawSector",fNTof+2);
  
  for(Int_t nSec=1; nSec<=fNTof; nSec++){
     AliTOFRawSector* currentSector = (AliTOFRawSector*)rawsectors->UncheckedAt(nSec);
     TClonesArray* rocData = (TClonesArray*)currentSector->GetRocData();
     for(Int_t nRoc=1; nRoc<=14; nRoc++){
        AliTOFRoc* currentRoc = (AliTOFRoc*)rocData->UncheckedAt(nRoc);
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
           Int_t pad = nPad%fPadXStr;
	   if (pad==0) pad=fPadXStr;
	   Int_t nPadX=0, nPadZ=0;
	   (pad>fNpadX)? nPadX -= fNpadX : nPadX = pad ;
	   vol[3] = nPadX;
	   (pad>fNpadX)? nPadZ = 2 : nPadZ = 1 ;
	   vol[4] = nPadZ;
	   UInt_t error=0;
	   Float_t tdc = currentRoc->GetTime(item,error);
	   if (!error) digit[0]=tdc;
	   digit[1] = currentRoc->GetCharge(item);
	   AddDigit(tracks,vol,digit);
        }
     }
  }
  tD->Fill();
  tD->Write();
} 


/******************************************************************************/

ClassImp(AliTOFhit)

//____________________________________________________________________________
AliTOFhit::AliTOFhit(const AliTOFhit & hit)
{
   //
   // copy ctor for AliTOFhit object
   //
  fTrack  = hit.fTrack;  
  fX      = hit.fX;
  fY      = hit.fY;
  fZ      = hit.fZ;
  fSector = hit.fSector;
  fPlate  = hit.fPlate;
  fStrip  = hit.fStrip;
  fPadx   = hit.fPadx;
  fPadz   = hit.fPadz;
  fPx     = hit.fPx;
  fPy     = hit.fPy;
  fPz     = hit.fPz;
  fPmom   = hit.fPmom;
  fTof    = hit.fTof;
  fDx     = hit.fDx;
  fDy     = hit.fDy;
  fDz     = hit.fDz;
  fIncA   = hit.fIncA;
  fEdep   = hit.fEdep;

}
 
//______________________________________________________________________________
AliTOFhit::AliTOFhit(Int_t shunt, Int_t track, Int_t *vol,
                     Float_t *hits)
:AliHit(shunt, track)
{
//
// Constructor of hit object
//
  //
  // Hit Volume
  // 
  fSector= vol[0];
  fPlate = vol[1];
  fStrip = vol[2];
  fPadx = vol[3];
  fPadz = vol[4];
  //
  //Position of the hit
  fX = hits[0];
  fY = hits[1];
  fZ = hits[2];
  //
  // Momentum components of the particle in the ALICE frame when hit is produced
  fPx  = hits[3];
  fPy  = hits[4];
  fPz  = hits[5];
  fPmom= hits[6];
  //
  // Time Of Flight for the particle that produces hit
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
{
//
// Constructor of digit object
//

  fSector = vol[0];
  fPlate  = vol[1];
  fStrip  = vol[2];
  fPadx  = vol[3];
  fPadz  = vol[4];
  fTdc    = digit[0];
  fAdc    = digit[1];
}

//____________________________________________________________________________
AliTOFdigit::AliTOFdigit(const AliTOFdigit & digit)
{
  // 
  // copy ctor for AliTOFdigit object
  //

  Int_t i ;
  for ( i = 0; i < 3 ; i++)
    fTracks[i]  = digit.fTracks[i] ;
  fSector = digit.fSector;
  fPlate  = digit.fPlate;
  fStrip  = digit.fStrip;
  fPadx   = digit.fPadx;
  fPadz   = digit.fPadz;
  fTdc    = digit.fTdc;
  fAdc    = digit.fAdc;

}
   
//______________________________________________________________________________
void AliTOFdigit::GetLocation(Int_t *Loc)
{
//
// Get the cohordinates of the digit
// in terms of Sector - Plate - Strip - Pad
//

   Loc[0]=fSector;
   Loc[1]=fPlate;
   Loc[2]=fStrip;
   Loc[3]=fPadx;
   Loc[4]=fPadz;
}

//______________________________________________________________________________
Int_t AliTOFdigit::GetTotPad()
{
//
// Get the "total" index of the pad inside a Sector
// starting from the digits data.
//

  AliTOF* tof;
  
  if(gAlice){
     tof =(AliTOF*) gAlice->GetDetector("TOF");
  }else{
     printf("AliTOFdigit::GetTotPad - No AliRun object present, exiting");
     return 0;
  }
  
  Int_t pad = fPadx+tof->GetNpadX()*(fPadz-1);
  Int_t before=0;

  switch(fPlate){ 
  case 1: before = 0;
          break;
  case 2: before = tof->GetNStripC();
          break;
  case 3: before = tof->GetNStripB() + tof->GetNStripC();
          break;
  case 4: before = tof->GetNStripA() + tof->GetNStripB() + tof->GetNStripC();
          break;
  case 5: before = tof->GetNStripA() + 2*tof->GetNStripB() + tof->GetNStripC();
          break;
  }
  
  Int_t strip = fStrip+before;
  Int_t padTot = tof->GetPadXStr()*(strip-1)+pad;
  return padTot;
}

//______________________________________________________________________________
void AliTOFdigit::AddTrack(Int_t track)
{
//
// Add a track to the digit 
//

  if (fTracks[1]==0){
     fTracks[1] = track;
  }else if (fTracks[2]==0){
     fTracks[2] = track;
  }else{
  //   printf("AliTOFdigit::AddTrack ERROR: Too many Tracks (>3) \n");
  }
}

