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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD MCM (Multi Chip Module) simulator                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/* $Id$ */

/*

  New release on 2007/08/17

AliTRDmcmSim is now stably working and zero suppression function seems ok.
From now, the default version of raw data is set to 3 in AliTRDfeeParam.

The following internal parameters were abolished because it is useless and
made trouble:

   fColOfADCbeg
   fColOfADCend

GetCol member was modified accordingly. 

New member function DumpData was prepared for diagnostics.

ZSMapping member function was debugged. It was causing crash due to
wrong indexing in 1 dimensional numbering. Also code was shaped up better.

*/

/*Semi-final version of TRD raw data simulation code with zero suppression (ZS)
similar to TRD FEE. ZS is realized by the class group:

  AliTRDfeeParam
  AliTRDmcmSim
  AliTRDrawData

AliTRDfeeParam has been modified to have more parameters like raw data
production version and so on. AliTRDmcmSim is new class and this is the core
of MCM (PASA+TRAP) simulator. It has still very simple function and it will be
another project to improve this to make it closer to the reall FEE.
AliTRDrawData has been modified to use new class AliTRDmcmSim.

These modifications were tested on Aug. 02 HEAD version that code itself
compiles. I'm sure there must be still bugs and we need testing by as many as
possible persons now. Especially it seems HLT part is impacted by problems
because some parameters were moved from AliTRDrawData to AliTRDfeeParam (like
fRawVersion disappeared from AliTRDrawData).

In TRD definition, we have now 4 raw data versions.

  0 very old offline version (by Bogdan)
  1 test version (no zero suppression)
  2 final version (no zero suppression)
  3 test version (with zero suppression)

The default is still set to 2 in AliTRDfeeParam::fgkRAWversion and it uses
previously existing codes. If you set this to 3, AliTRDrawData changes behavior
to use AliTRDmcmSim with ZS.

Plan is after we make sure it works stably, we delete AliTRDmcm which is obsolete.
However it still take time because tracklet part is not yet touched.
The default raw version is 2.

                                                                 Ken Oyama
*/

// if no histo is drawn, these are obsolete
#include <TH1.h>
#include <TCanvas.h>

// only needed if I/O of tracklets is activated
#include <TObject.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

#include <fstream>

#include <TMath.h>

#include "AliLog.h"

#include "AliTRDmcmSim.h"
#include "AliTRDfeeParam.h"
#include "AliTRDSimParam.h"
#include "AliTRDCommonParam.h"
#include "AliTRDgeometry.h"
#include "AliTRDcalibDB.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDarrayADC.h"
// additional for new tail filter and/or tracklet
#include "AliTRDtrapAlu.h"
#include "AliTRDpadPlane.h"
#include "AliTRDtrackletMCM.h"

#include "AliRun.h"
#include "AliLoader.h"

ClassImp(AliTRDmcmSim)

//_____________________________________________________________________________
AliTRDmcmSim::AliTRDmcmSim() :TObject()
  ,fInitialized(kFALSE)
  ,fNextEvent(-1)    
  ,fMaxTracklets(-1) 
  ,fChaId(-1)
  ,fSector(-1)
  ,fStack(-1)
  ,fLayer(-1)
  ,fRobPos(-1)
  ,fMcmPos(-1)
  ,fNADC(-1)
  ,fNTimeBin(-1)
  ,fRow (-1)
  ,fADCR(NULL)
  ,fADCF(NULL)
  ,fADCT(NULL)     
  ,fPosLUT(NULL)    
  ,fMCMT(NULL)      
  ,fZSM(NULL)
  ,fZSM1Dim(NULL)
  ,fFeeParam(NULL)
  ,fSimParam(NULL)
  ,fCal(NULL)
  ,fGeo(NULL)
{
  //
  // AliTRDmcmSim default constructor
  //

  // By default, nothing is initialized.
  // It is necessary to issue Init before use.
}

//_____________________________________________________________________________
AliTRDmcmSim::AliTRDmcmSim(const AliTRDmcmSim &m) 
  :TObject(m)
  ,fInitialized(kFALSE) 
  ,fNextEvent(-1)    
  ,fMaxTracklets(-1) 
  ,fChaId(-1)
  ,fSector(-1)
  ,fStack(-1)
  ,fLayer(-1)
  ,fRobPos(-1)
  ,fMcmPos(-1)
  ,fNADC(-1)
  ,fNTimeBin(-1)
  ,fRow(-1)
  ,fADCR(NULL)
  ,fADCF(NULL)
  ,fADCT(NULL)      
  ,fPosLUT(NULL)    
  ,fMCMT(NULL)      
  ,fZSM(NULL)
  ,fZSM1Dim(NULL)
  ,fFeeParam(NULL)
  ,fSimParam(NULL)
  ,fCal(NULL)
  ,fGeo(NULL)
  
{
  //
  // AliTRDmcmSim copy constructor
  //

  // By default, nothing is initialized.
  // It is necessary to issue Init before use.
}

//_____________________________________________________________________________
AliTRDmcmSim::~AliTRDmcmSim() 
{
  //
  // AliTRDmcmSim destructor
  //

  if( fADCR != NULL ) {
    for( Int_t iadc = 0 ; iadc < fNADC; iadc++ ) {
      delete [] fADCR[iadc];
      delete [] fADCF[iadc];
      delete [] fADCT[iadc];
      delete [] fZSM [iadc];
    }
    delete [] fADCR;
    delete [] fADCF;
    delete [] fADCT;
    delete [] fZSM;
    delete [] fZSM1Dim;
  }
 
  if(fInitialized){
    delete [] fPosLUT;
    delete [] fMCMT;
  }

  delete fGeo;

}

//_____________________________________________________________________________
AliTRDmcmSim &AliTRDmcmSim::operator=(const AliTRDmcmSim &m)
{
  //
  // Assignment operator
  //

  if (this != &m) {
    ((AliTRDmcmSim &) m).Copy(*this);
  }
  return *this;

}

//_____________________________________________________________________________
void AliTRDmcmSim::Copy(TObject &m) const
{
  //
  // Copy function
  //
  ((AliTRDmcmSim &) m).fNextEvent     = 0; //new
  ((AliTRDmcmSim &) m).fMaxTracklets  = 0; //new
  ((AliTRDmcmSim &) m).fInitialized   = 0;
  ((AliTRDmcmSim &) m).fChaId         = 0;
  ((AliTRDmcmSim &) m).fSector        = 0;
  ((AliTRDmcmSim &) m).fStack         = 0;
  ((AliTRDmcmSim &) m).fLayer         = 0;
  ((AliTRDmcmSim &) m).fRobPos        = 0;
  ((AliTRDmcmSim &) m).fMcmPos        = 0;
  ((AliTRDmcmSim &) m).fNADC          = 0;
  ((AliTRDmcmSim &) m).fNTimeBin      = 0;
  ((AliTRDmcmSim &) m).fRow           = 0;
  ((AliTRDmcmSim &) m).fADCR          = 0;
  ((AliTRDmcmSim &) m).fADCF          = 0;
  ((AliTRDmcmSim &) m).fADCT          = 0; //new
  ((AliTRDmcmSim &) m).fPosLUT        = 0; //new
  ((AliTRDmcmSim &) m).fMCMT          = 0; //new
  ((AliTRDmcmSim &) m).fZSM           = 0;
  ((AliTRDmcmSim &) m).fZSM1Dim       = 0;
  ((AliTRDmcmSim &) m).fFeeParam      = 0;
  ((AliTRDmcmSim &) m).fSimParam      = 0;
  ((AliTRDmcmSim &) m).fCal           = 0;
  ((AliTRDmcmSim &) m).fGeo           = 0;

}

//_____________________________________________________________________________

//void AliTRDmcmSim::Init( Int_t chaId, Int_t robPos, Int_t mcmPos ) 
void AliTRDmcmSim::Init( Int_t chaId, Int_t robPos, Int_t mcmPos, Bool_t newEvent = kFALSE ) // only for readout tree (new event)
{
  //
  // Initialize the class with new geometry information
  // fADC array will be reused with filled by zero
  //
   
  fNextEvent     = 0; 
  fFeeParam      = AliTRDfeeParam::Instance();
  fSimParam      = AliTRDSimParam::Instance();
  fCal           = AliTRDcalibDB::Instance();
  fGeo           = new AliTRDgeometry();
  fChaId         = chaId;
  fSector        = fGeo->GetSector( fChaId );
  fStack         = fGeo->GetStack( fChaId );
  fLayer         = fGeo->GetLayer( fChaId );
  fRobPos        = robPos;
  fMcmPos        = mcmPos;
  fNADC          = fFeeParam->GetNadcMcm();
  fNTimeBin      = fCal->GetNumberOfTimeBins();
  fRow           = fFeeParam->GetPadRowFromMCM( fRobPos, fMcmPos );

  fMaxTracklets  = fFeeParam->GetMaxNrOfTracklets();

 

  
  if (newEvent == kTRUE) {
      fNextEvent = 1;
  }



  // Allocate ADC data memory if not yet done
  if( fADCR == NULL ) {
    fADCR    = new Int_t *[fNADC];
    fADCF    = new Int_t *[fNADC];
    fADCT    = new Int_t *[fNADC]; //new
    fZSM     = new Int_t *[fNADC];
    fZSM1Dim = new Int_t  [fNADC];
    for( Int_t iadc = 0 ; iadc < fNADC; iadc++ ) {
      fADCR[iadc] = new Int_t[fNTimeBin];
      fADCF[iadc] = new Int_t[fNTimeBin];
      fADCT[iadc] = new Int_t[fNTimeBin]; //new
      fZSM [iadc] = new Int_t[fNTimeBin];
    }
  }

  // Initialize ADC data
  for( Int_t iadc = 0 ; iadc < fNADC; iadc++ ) {
    for( Int_t it = 0 ; it < fNTimeBin ; it++ ) {
      fADCR[iadc][it] = 0;
      fADCF[iadc][it] = 0;
      fADCT[iadc][it] = -1;  //new
      fZSM [iadc][it] = 1;   // Default unread = 1
    }
    fZSM1Dim[iadc] = 1;      // Default unread = 1
  }
  
  //new:
  fPosLUT = new Int_t[128];
  for(Int_t i = 0; i<128; i++){
    fPosLUT[i] = 0;
  }
  
  fMCMT = new UInt_t[fMaxTracklets];
  for(Int_t i = 0; i < fMaxTracklets; i++) {
    fMCMT[i] = 0;
  }


  fInitialized = kTRUE;
}

//_____________________________________________________________________________
Bool_t AliTRDmcmSim::CheckInitialized()
{
  //
  // Check whether object is initialized
  //

  if( ! fInitialized ) {
    AliDebug(2, Form ("AliTRDmcmSim is not initialized but function other than Init() is called."));
  }
  return fInitialized;
}

//_____________________________________________________________________________


void AliTRDmcmSim::SetPosLUT() {
  Double_t iHi  = (Double_t)fCal->GetPRFhi();
  Double_t iLo  = (Double_t)fCal->GetPRFlo();
  Int_t   nBin  = fCal->GetPRFbin();
  Int_t   iOff  = fLayer * nBin;
  Int_t kNlayer = fGeo->Nlayer();

  Float_t  *sPRFsmp   = new Float_t[nBin*kNlayer];
  Double_t *sPRFlayer = new Double_t[nBin];
  
  
  for(Int_t i = 0; i<nBin*kNlayer; i++){
    
    //printf("%f\n",fCal->GetSampledPRF()[i]);
    sPRFsmp[i] = fCal->GetSampledPRF()[i]; 
  
  }

  Double_t sWidth = (iHi-iLo)/((Double_t) nBin);
  Int_t   sPad    = (Int_t) (1.0/sWidth);
  
  // get the PRF for actual layer (interpolated to ibin data-points; 61 measured)
  for(Int_t iBin = 0; iBin < nBin; iBin++){
    sPRFlayer[iBin] = (Double_t)sPRFsmp[iOff+iBin];
  }

  Int_t bin0 = (Int_t)(-iLo / sWidth - 0.5);                           // bin-nr. for pad-position 0
  
  Int_t bin1 = (Int_t)((Double_t)(0.5 - iLo) / sWidth - 0.5);          // bin-nr. for pad-position 0.5
  bin1 = bin1 + 1;
  bin0 = bin0 + 1;  //avoid negative values in aYest (start right of symmetry center)
  while (bin0-sPad<0) {
    bin0 = bin0 + 1;
  }
  while (bin1+sPad>=nBin) {
    bin1 = bin1 - 1;
  }
  
  Double_t* aYest = new Double_t[bin1-bin0+1];

  /*TH1F* hist1 = new TH1F("h1","yest(y)",128,0,0.5);
  TH1F* hist2 = new TH1F("h2","y(yest)",128,0,0.5);
  TH1F* hist3 = new TH1F("h3","y(yest)-yest",128,0,0.5);
  TH1F* hist4 = new TH1F("h4","y(yest)-yest,discrete",128,0,0.5);
 
  TCanvas *c1 = new TCanvas("c1","c1",800,1000);
  hist1->Draw();
  TCanvas *c2 = new TCanvas("c2","c2",800,1000);
  hist2->Draw();
  TCanvas *c3 = new TCanvas("c3","c3",800,1000);
  hist3->Draw();
  TCanvas *c4 = new TCanvas("c4","c4",800,1000);
  hist4->Draw();*/
  
  for(Int_t iBin = bin0; iBin <= bin1; iBin++){
    aYest[iBin-bin0] = 0.5*(sPRFlayer[iBin-sPad] - sPRFlayer[iBin+sPad])/(sPRFlayer[iBin]); // estimated position from PRF; between 0 and 1
    //Double_t position = ((Double_t)(iBin)+0.5)*sWidth+iLo;
    //  hist1->Fill(position,aYest[iBin-bin0]);
  }
  


  Double_t aY[128]; // reversed function

  AliTRDtrapAlu a;
  a.Init(1,8,0,31);
  
  for(Int_t j = 0; j<128; j++) { // loop over all Yest; LUT has 128 entries; 
    Double_t yest = ((Double_t)j)/256; 
    
    Int_t iBin = 0;
    while (yest>aYest[iBin] && iBin<(bin1-bin0)) {
      iBin = iBin+1;
    }
    if((iBin == bin1 - bin0)&&(yest>aYest[iBin])) {
      aY[j] = 0.5;                      // yest too big
      //hist2->Fill(yest,aY[j]);
      
    }
    else {
      Int_t bin_d = iBin + bin0 - 1;
      Int_t bin_u = iBin + bin0;
      Double_t y_d = ((Double_t)bin_d + 0.5)*sWidth + iLo; // lower y
      Double_t y_u = ((Double_t)bin_u + 0.5)*sWidth + iLo; // upper y
      Double_t yest_d = aYest[iBin-1];                     // lower estimated y
      Double_t yest_u = aYest[iBin];                       // upper estimated y
      
      aY[j] = ((yest-yest_d)/(yest_u-yest_d))*(y_u-y_d) + y_d;
      //hist2->Fill(yest,aY[j]);
     
    }
    aY[j] = aY[j] - yest;
    //hist3->Fill(yest,aY[j]);
    // formatting
    a.AssignDouble(aY[j]);
    //a.WriteWord();
    fPosLUT[j] = a.GetValue(); // 1+8Bit value;128 entries;LUT is steered by abs(Q(i+1)-Q(i-1))/Q(i)=COG and gives the correction to COG/2
    //hist4->Fill(yest,fPosLUT[j]);
    
  }
 
   
  
  delete [] sPRFsmp;
  delete [] sPRFlayer;
  delete [] aYest;
  
}


//_____________________________________________________________________________
Int_t* AliTRDmcmSim::GetPosLUT(){
  return fPosLUT;
}



void AliTRDmcmSim::SetData( Int_t iadc, Int_t *adc )
{
  //
  // Store ADC data into array of raw data
  //

  if( !CheckInitialized() ) return;

  if( iadc < 0 || iadc >= fNADC ) {
    //Log (Form ("Error: iadc is out of range (should be 0 to %d).", fNADC-1));
    return;
  }

  for( int it = 0 ;  it < fNTimeBin ; it++ ) {
    fADCR[iadc][it] = (Int_t)(adc[it]);
  }
}

//_____________________________________________________________________________
void AliTRDmcmSim::SetData( Int_t iadc, Int_t it, Int_t adc )
{
  //
  // Store ADC data into array of raw data
  //

  if( !CheckInitialized() ) return;

  if( iadc < 0 || iadc >= fNADC ) {
    //Log (Form ("Error: iadc is out of range (should be 0 to %d).", fNADC-1));
    return;
  }

  fADCR[iadc][it] = adc;
}

//_____________________________________________________________________________
void AliTRDmcmSim::SetDataPedestal( Int_t iadc )
{
  //
  // Store ADC data into array of raw data
  //

  if( !CheckInitialized() ) return;

  if( iadc < 0 || iadc >= fNADC ) {
    //Log (Form ("Error: iadc is out of range (should be 0 to %d).", fNADC-1));
    return;
  }

  for( Int_t it = 0 ; it < fNTimeBin ; it++ ) {
    fADCR[iadc][it] = fSimParam->GetADCbaseline();
  }
}

//_____________________________________________________________________________
Int_t AliTRDmcmSim::GetCol( Int_t iadc )
{
  //
  // Return column id of the pad for the given ADC channel
  //

  if( !CheckInitialized() ) return -1;

  return fFeeParam->GetPadColFromADC(fRobPos, fMcmPos, iadc);
}

//_____________________________________________________________________________
Int_t AliTRDmcmSim::ProduceRawStream( UInt_t *buf, Int_t maxSize )
{
  //
  // Produce raw data stream from this MCM and put in buf
  // Returns number of words filled, or negative value 
  // with -1 * number of overflowed words
  //

  UInt_t  x;
  UInt_t  iEv = 0;
  Int_t   nw  = 0;  // Number of written words
  Int_t   of  = 0;  // Number of overflowed words
  Int_t   rawVer   = fFeeParam->GetRAWversion();
  Int_t **adc;

  if( !CheckInitialized() ) return 0;

  if( fFeeParam->GetRAWstoreRaw() ) {
    adc = fADCR;
  } else {
    adc = fADCF;
  }

  // Produce MCM header
  x = ((fRobPos * fFeeParam->GetNmcmRob() + fMcmPos) << 24) | ((iEv % 0x100000) << 4) | 0xC;
  if (nw < maxSize) {
    buf[nw++] = x;
  }
  else {
    of++;
  }

  // Produce ADC mask
  if( rawVer >= 3 ) {
    x = 0;
    for( Int_t iAdc = 0 ; iAdc < fNADC ; iAdc++ ) {
      if( fZSM1Dim[iAdc] == 0 ) { //  0 means not suppressed
	x = x | (1 << iAdc);
      }
    }
    if (nw < maxSize) {
      buf[nw++] = x;
    }
    else {
      of++;
    }
  }

  // Produce ADC data. 3 timebins are packed into one 32 bits word
  // In this version, different ADC channel will NOT share the same word

  UInt_t aa=0, a1=0, a2=0, a3=0;

  for (Int_t iAdc = 0; iAdc < 21; iAdc++ ) {
    if( rawVer>= 3 && fZSM1Dim[iAdc] != 0 ) continue; // suppressed
    aa = !(iAdc & 1) + 2;
    for (Int_t iT = 0; iT < fNTimeBin; iT+=3 ) {
      a1 = ((iT    ) < fNTimeBin ) ? adc[iAdc][iT  ] : 0;
      a2 = ((iT + 1) < fNTimeBin ) ? adc[iAdc][iT+1] : 0;
      a3 = ((iT + 2) < fNTimeBin ) ? adc[iAdc][iT+2] : 0;
      x = (a3 << 22) | (a2 << 12) | (a1 << 2) | aa;
      if (nw < maxSize) {
	buf[nw++] = x;
      }
      else {
	of++;
      }
    }
  }

  if( of != 0 ) return -of; else return nw;
}

//_____________________________________________________________________________
Int_t AliTRDmcmSim::ProduceRawStreamV2( UInt_t *buf, Int_t maxSize, UInt_t iEv )
{
  //
  // Produce raw data stream from this MCM and put in buf
  // Returns number of words filled, or negative value 
  // with -1 * number of overflowed words
  //

  UInt_t  x;
  //UInt_t  iEv = 0;
  Int_t   nw  = 0;  // Number of written words
  Int_t   of  = 0;  // Number of overflowed words
  Int_t   rawVer   = fFeeParam->GetRAWversion();
  Int_t **adc;
  Int_t   nActiveADC = 0;	// number of activated ADC bits in a word

  if( !CheckInitialized() ) return 0;

  if( fFeeParam->GetRAWstoreRaw() ) {
    adc = fADCR;
  } else {
    adc = fADCF;
  }

  // Produce MCM header : xrrr mmmm eeee eeee eeee eeee eeee 1100
  //                      x : 0 before , 1 since 10.2007
  //                      r : Readout board position (Alice numbering)
  //                      m : MCM posi
  //                      e : Event counter from 1
  //x = (1<<31) | ((fRobPos * fFeeParam->GetNmcmRob() + fMcmPos) << 24) | ((iEv % 0x100000) << 4) | 0xC;
  x = (1<<31) | (fRobPos << 28) | (fMcmPos << 24) | ((iEv % 0x100000) << 4) | 0xC;
  if (nw < maxSize) {
    buf[nw++] = x;
	//printf("\nMCM header: %X ",x);
  }
  else {
    of++;
  }

  // Produce ADC mask : nncc cccm mmmm mmmm mmmm mmmm mmmm 1100
  // 				n : unused , c : ADC count, m : selected ADCs
  if( rawVer >= 3 ) {
    x = 0;
    for( Int_t iAdc = 0 ; iAdc < fNADC ; iAdc++ ) {
      if( fZSM1Dim[iAdc] == 0 ) { //  0 means not suppressed
		x = x | (1 << (iAdc+4) );	// last 4 digit reserved for 1100=0xc
		nActiveADC++;		// number of 1 in mmm....m
      }
    }
	x = x | (1 << 30) | ( ( 0x3FFFFFFC ) & (~(nActiveADC) << 25) ) | 0xC;	// nn = 01, ccccc are inverted, 0xc=1100
	//printf("nActiveADC=%d=%08X, inverted=%X ",nActiveADC,nActiveADC,x );

    if (nw < maxSize) {
      buf[nw++] = x;
	  //printf("ADC mask: %X nMask=%d ADC data: ",x,nActiveADC);
    }
    else {
      of++;
    }
  }

  // Produce ADC data. 3 timebins are packed into one 32 bits word
  // In this version, different ADC channel will NOT share the same word

  UInt_t aa=0, a1=0, a2=0, a3=0;

  for (Int_t iAdc = 0; iAdc < 21; iAdc++ ) {
    if( rawVer>= 3 && fZSM1Dim[iAdc] != 0 ) continue; // Zero Suppression, 0 means not suppressed
    aa = !(iAdc & 1) + 2;	// 3 for the even ADC channel , 2 for the odd ADC channel
    for (Int_t iT = 0; iT < fNTimeBin; iT+=3 ) {
      a1 = ((iT    ) < fNTimeBin ) ? adc[iAdc][iT  ] : 0;
      a2 = ((iT + 1) < fNTimeBin ) ? adc[iAdc][iT+1] : 0;
      a3 = ((iT + 2) < fNTimeBin ) ? adc[iAdc][iT+2] : 0;
      x = (a3 << 22) | (a2 << 12) | (a1 << 2) | aa;
      if (nw < maxSize) {
	buf[nw++] = x;
	//printf("%08X ",x);
      }
      else {
	of++;
      }
    }
  }

  if( of != 0 ) return -of; else return nw;
}

//_____________________________________________________________________________
Int_t AliTRDmcmSim::ProduceTrackletStream( UInt_t *buf, Int_t maxSize )
{
  //
  // Produce tracklet data stream from this MCM and put in buf
  // Returns number of words filled, or negative value 
  // with -1 * number of overflowed words
  //

  UInt_t  x;
  Int_t   nw  = 0;  // Number of written words
  Int_t   of  = 0;  // Number of overflowed words
    
  if( !CheckInitialized() ) return 0;

  // Produce tracklet data. A maximum of four 32 Bit words will be written per MCM 
  // fMCMT is filled continuously until no more tracklet words available

  Int_t wd = 0;
  while ( (wd < fMaxTracklets) && (fMCMT[wd] > 0) ){
      x = fMCMT[wd];
      if (nw < maxSize) {
	buf[nw++] = x;
      }
      else {
	of++;
      }
      wd++;
  }
  
  if( of != 0 ) return -of; else return nw;
}


//_____________________________________________________________________________
void AliTRDmcmSim::Filter()
{
  //
  // Apply digital filter
  //

  if( !CheckInitialized() ) return;

  // Initialize filtered data array with raw data
  for( Int_t iadc = 0 ; iadc < fNADC; iadc++ ) {
    for( Int_t it = 0 ; it < fNTimeBin ; it++ ) {
      fADCF[iadc][it] = fADCR[iadc][it]; 
    }
  }

  // Then apply fileters one by one to filtered data array
  if( fFeeParam->IsPFon() ) FilterPedestal();
  if( fFeeParam->IsGFon() ) FilterGain();
  if( fFeeParam->IsTFon() ) FilterTail();
}

//_____________________________________________________________________________
void AliTRDmcmSim::FilterPedestal()
{

     
  //
  // Apply pedestal filter
  //

  Int_t ap = fSimParam->GetADCbaseline();      // ADC instrinsic pedestal
  Int_t ep = fFeeParam->GetPFeffectPedestal(); // effective pedestal
  //Int_t tc = fFeeParam->GetPFtimeConstant();   // this makes no sense yet

  for( Int_t iadc = 0 ; iadc < fNADC; iadc++ ) {
    for( Int_t it = 0 ; it < fNTimeBin ; it++ ) {
      fADCF[iadc][it] = fADCF[iadc][it] - ap + ep;
    }
  }
}

//_____________________________________________________________________________
void AliTRDmcmSim::FilterGain()
{
  //
  // Apply gain filter (not implemented)
  // Later it will be implemented because gain digital filiter will
  // increase noise level.
  //

}

//_____________________________________________________________________________
void AliTRDmcmSim::FilterTail()
{
  //
  // Apply exponential tail filter (Bogdan's version)
  //

  Double_t *dtarg  = new Double_t[fNTimeBin];
  Int_t    *itarg  = new Int_t[fNTimeBin];
  Int_t     nexp   = fFeeParam->GetTFnExp();
  Int_t     tftype = fFeeParam->GetTFtype();

  switch( tftype ) {
    
  case 0: // Exponential Filter Analog Bogdan
    for (Int_t iCol = 0; iCol < fNADC; iCol++) {
      FilterSimDeConvExpA( fADCF[iCol], dtarg, fNTimeBin, nexp);
      for (Int_t iTime = 0; iTime < fNTimeBin; iTime++) {
	fADCF[iCol][iTime] = (Int_t) TMath::Max(0.0,dtarg[iTime]);
      }
    }
    break;

  case 1: // Exponential filter digital Bogdan
    for (Int_t iCol = 0; iCol < fNADC; iCol++) {
      FilterSimDeConvExpD( fADCF[iCol], itarg, fNTimeBin, nexp);
      for (Int_t iTime = 0; iTime < fNTimeBin; iTime++) {
	fADCF[iCol][iTime] = itarg[iTime];
      }
    }
    break;
    
  case 2: // Exponential filter Marian special
    for (Int_t iCol = 0; iCol < fNADC; iCol++) {
      FilterSimDeConvExpMI( fADCF[iCol], dtarg, fNTimeBin);
      for (Int_t iTime = 0; iTime < fNTimeBin; iTime++) {
	fADCF[iCol][iTime] = (Int_t) TMath::Max(0.0,dtarg[iTime]);
      }
    }
    break;

    //new
  case 3: // Exponential filter using AliTRDtrapAlu class
    for (Int_t iCol = 0; iCol < fNADC; iCol++) {
      FilterSimDeConvExpEl( fADCF[iCol], itarg, fNTimeBin, nexp);
      for (Int_t iTime = 0; iTime < fNTimeBin; iTime++) {
	fADCF[iCol][iTime] = itarg[iTime]>>2; // to be used for raw-data
	fADCT[iCol][iTime] = itarg[iTime];    // 12bits; to be used for tracklet; tracklet will have own container; 
      }
    }
    break;

    
  default:
    AliError(Form("Invalid filter type %d ! \n", tftype ));
    break;
  }

  delete [] dtarg;
  delete [] itarg;

}

//_____________________________________________________________________________
void AliTRDmcmSim::ZSMapping()
{
  //
  // Zero Suppression Mapping implemented in TRAP chip
  //
  // See detail TRAP manual "Data Indication" section:
  // http://www.kip.uni-heidelberg.de/ti/TRD/doc/trap/TRAP-UserManual.pdf
  //

  Int_t eBIS = fFeeParam->GetEBsglIndThr();       // TRAP default = 0x4  (Tis=4)
  Int_t eBIT = fFeeParam->GetEBsumIndThr();       // TRAP default = 0x28 (Tit=40)
  Int_t eBIL = fFeeParam->GetEBindLUT();          // TRAP default = 0xf0
                                                  // (lookup table accept (I2,I1,I0)=(111)
                                                  // or (110) or (101) or (100))
  Int_t eBIN = fFeeParam->GetEBignoreNeighbour(); // TRAP default = 1 (no neighbor sensitivity)
  Int_t ep   = AliTRDfeeParam::GetPFeffectPedestal();

  if( !CheckInitialized() ) return;

  for( Int_t iadc = 1 ; iadc < fNADC-1; iadc++ ) {
    for( Int_t it = 0 ; it < fNTimeBin ; it++ ) {

      // Get ADC data currently in filter buffer
      Int_t ap = fADCF[iadc-1][it] - ep; // previous
      Int_t ac = fADCF[iadc  ][it] - ep; // current
      Int_t an = fADCF[iadc+1][it] - ep; // next

      // evaluate three conditions
      Int_t i0 = ( ac >=  ap && ac >=  an ) ? 0 : 1; // peak center detection
      Int_t i1 = ( ap + ac + an > eBIT )    ? 0 : 1; // cluster
      Int_t i2 = ( ac > eBIS )              ? 0 : 1; // absolute large peak

      Int_t i = i2 * 4 + i1 * 2 + i0;    // Bit position in lookup table
      Int_t d = (eBIL >> i) & 1;         // Looking up  (here d=0 means true
                                         // and d=1 means false according to TRAP manual)

      fZSM[iadc][it] &= d;
      if( eBIN == 0 ) {  // turn on neighboring ADCs
	fZSM[iadc-1][it] &= d;
	fZSM[iadc+1][it] &= d;
      }

    }
  }

  // do 1 dim projection
  for( Int_t iadc = 0 ; iadc < fNADC; iadc++ ) {
    for( Int_t it = 0 ; it < fNTimeBin ; it++ ) {
      fZSM1Dim[iadc] &= fZSM[iadc][it];
    }
  }

}

//_____________________________________________________________________________
void AliTRDmcmSim::DumpData( char *f, char *target )
{
  //
  // Dump data stored (for debugging).
  // target should contain one or multiple of the following characters
  //   R   for raw data
  //   F   for filtered data
  //   Z   for zero suppression map
  //   S   Raw dat astream
  // other characters are simply ignored
  //

  UInt_t tempbuf[1024];

  if( !CheckInitialized() ) return;

  std::ofstream of( f, std::ios::out | std::ios::app );
  of << Form("AliTRDmcmSim::DumpData det=%03d sm=%02d stack=%d layer=%d rob=%d mcm=%02d\n",
	     fChaId, fSector, fStack, fLayer, fRobPos, fMcmPos );

  for( int t=0 ; target[t] != 0 ; t++ ) {
    switch( target[t] ) {
    case 'R' :
    case 'r' :
      of << Form("fADCR (raw ADC data)\n");
      for( Int_t iadc = 0 ; iadc < fNADC; iadc++ ) {
	of << Form("  ADC %02d: ", iadc);
	for( Int_t it = 0 ; it < fNTimeBin ; it++ ) {
	  of << Form("% 4d",  fADCR[iadc][it]);
	}
	of << Form("\n");
      }
      break;
    case 'F' :
    case 'f' :
      of << Form("fADCF (filtered ADC data)\n");
      for( Int_t iadc = 0 ; iadc < fNADC; iadc++ ) {
	of << Form("  ADC %02d: ", iadc);
	for( Int_t it = 0 ; it < fNTimeBin ; it++ ) {
	  of << Form("% 4d",  fADCF[iadc][it]);
	}
	of << Form("\n");
      }
      break;
    case 'Z' :
    case 'z' :
      of << Form("fZSM and fZSM1Dim (Zero Suppression Map)\n");
      for( Int_t iadc = 0 ; iadc < fNADC; iadc++ ) {
	of << Form("  ADC %02d: ", iadc);
	if( fZSM1Dim[iadc] == 0 ) { of << " R   " ; } else { of << " .   "; } // R:read .:suppressed
	for( Int_t it = 0 ; it < fNTimeBin ; it++ ) {
	  if( fZSM[iadc][it] == 0 ) { of << " R"; } else { of << " ."; } // R:read .:suppressed
	}
	of << Form("\n");
      }
      break;
    case 'S' :
    case 's' :
      Int_t s = ProduceRawStream( tempbuf, 1024 ); 
      of << Form("Stream for Raw Simulation size=%d rawver=%d\n", s, fFeeParam->GetRAWversion());
      of << Form("  address  data\n");
      for( int i = 0 ; i < s ; i++ ) {
	of << Form("  %04x     %08x\n", i, tempbuf[i]);
      }
    }
  }
}

//_____________________________________________________________________________
void AliTRDmcmSim::FilterSimDeConvExpA(Int_t *source, Double_t *target
                                     , Int_t n, Int_t nexp) 
{
  //
  // Exponential filter "analog"
  // source will not be changed
  //

  Int_t    i = 0;
  Int_t    k = 0;
  Double_t reminder[2];
  Double_t correction;
  Double_t result;
  Double_t rates[2];
  Double_t coefficients[2];

  // Initialize (coefficient = alpha, rates = lambda)
  // FilterOpt.C (aliroot@pel:/homel/aliroot/root/work/beamt/CERN02)

  Double_t r1 = (Double_t)fFeeParam->GetTFr1();
  Double_t r2 = (Double_t)fFeeParam->GetTFr2();
  Double_t c1 = (Double_t)fFeeParam->GetTFc1();
  Double_t c2 = (Double_t)fFeeParam->GetTFc2();
  
  coefficients[0] = c1;
  coefficients[1] = c2;

  Double_t dt = 0.1;
  rates[0] = TMath::Exp(-dt/(r1));
  rates[1] = TMath::Exp(-dt/(r2));

  // Attention: computation order is important
  correction = 0.0;
  for (k = 0; k < nexp; k++) {
    reminder[k] = 0.0;
  }
    
  for (i = 0; i < n; i++) {

    result    = ((Double_t)source[i] - correction);    // no rescaling
    target[i] = result;
    
    for (k = 0; k < nexp; k++) {
      reminder[k] = rates[k] * (reminder[k] + coefficients[k] * result);
    }
      
    correction = 0.0;
    for (k = 0; k < nexp; k++) {
      correction += reminder[k];
    }
  }
}

//_____________________________________________________________________________
void AliTRDmcmSim::FilterSimDeConvExpD(Int_t *source, Int_t *target, Int_t n
                                     , Int_t nexp) 
{
  //
  // Exponential filter "digital"
  // source will not be changed
  //

  Int_t i        = 0;
  Int_t fAlphaL  = 0;
  Int_t fAlphaS  = 0;
  Int_t fTailPed = 0;
  Int_t iAlphaL  = 0;
  Int_t iAlphaS  = 0;

  // FilterOpt.C (aliroot@pel:/homel/aliroot/root/work/beamt/CERN02)
  // initialize (coefficient = alpha, rates = lambda)

  Double_t dt = 0.1;
  Double_t r1 = (Double_t)fFeeParam->GetTFr1();
  Double_t r2 = (Double_t)fFeeParam->GetTFr2();
  Double_t c1 = (Double_t)fFeeParam->GetTFc1();
  Double_t c2 = (Double_t)fFeeParam->GetTFc2();

  Int_t fLambdaL = (Int_t)((TMath::Exp(-dt/r1) - 0.75) * 2048.0);
  Int_t fLambdaS = (Int_t)((TMath::Exp(-dt/r2) - 0.25) * 2048.0);
  Int_t iLambdaL = fLambdaL & 0x01FF; iLambdaL |= 0x0600; //  9 bit paramter + fixed bits
  Int_t iLambdaS = fLambdaS & 0x01FF; iLambdaS |= 0x0200; //  9 bit paramter + fixed bits

  if (nexp == 1) {
    fAlphaL = (Int_t) (c1 * 2048.0);
    iAlphaL = fAlphaL & 0x03FF;				// 10 bit paramter
  }
  if (nexp == 2) {
    fAlphaL = (Int_t) (c1 * 2048.0);
    fAlphaS = (Int_t) ((c2 - 0.5) * 2048.0);
    iAlphaL = fAlphaL & 0x03FF;				// 10 bit paramter
    iAlphaS = fAlphaS & 0x03FF; iAlphaS |= 0x0400;	        // 10 bit paramter + fixed bits
  }
  
  Double_t iAl = iAlphaL  / 2048.0;	       // alpha  L: correspondence to floating point numbers
  Double_t iAs = iAlphaS  / 2048.0;	       // alpha  S: correspondence to floating point numbers
  Double_t iLl = iLambdaL / 2048.0;	       // lambda L: correspondence to floating point numbers
  Double_t iLs = iLambdaS / 2048.0;	       // lambda S: correspondence to floating point numbers

  Int_t h1;
  Int_t h2;
  Int_t rem1;
  Int_t rem2;
  Int_t correction;
  Int_t result;
  Int_t iFactor = ((Int_t) fFeeParam->GetPFeffectPedestal() ) << 2;

  Double_t xi = 1 - (iLl*iAs + iLs*iAl);	     // Calculation of equilibrium values of the
  rem1 = (Int_t) ((iFactor/xi) * ((1-iLs)*iLl*iAl)); // Internal registers to prevent switch on effects.
  rem2 = (Int_t) ((iFactor/xi) * ((1-iLl)*iLs*iAs));
  
  // further initialization
  if ((rem1 + rem2) > 0x0FFF) {
    correction = 0x0FFF;
  } 
  else {
    correction = (rem1 + rem2) & 0x0FFF;
  }

  fTailPed = iFactor - correction;

  for (i = 0; i < n; i++) {

    result = (source[i]  - correction);
    if (result < 0) { // Too much undershoot
      result = 0;
    }

    target[i] = result;
                                                        
    h1 = (rem1 + ((iAlphaL * result) >> 11));
    if (h1 > 0x0FFF) {
      h1 = 0x0FFF;
    } 
    else {
      h1 &= 0x0FFF;
    }

    h2 = (rem2 + ((iAlphaS * result) >> 11));
    if (h2 > 0x0FFF) {
      h2 = 0x0FFF;
    } 
    else {
      h2 &= 0x0FFF;
    }
  
    rem1 = (iLambdaL * h1 ) >> 11;
    rem2 = (iLambdaS * h2 ) >> 11;
    
    if ((rem1 + rem2) > 0x0FFF) {
      correction = 0x0FFF;
    } 
    else {
      correction = (rem1 + rem2) & 0x0FFF;
    }

  }

}

//_____________________________________________________________________________
void AliTRDmcmSim::FilterSimDeConvExpMI(Int_t *source, Double_t *target
                                      , Int_t n) 
{
  //
  // Exponential filter (M. Ivanov)
  // source will not be changed
  //

  Int_t i = 0;
  Double_t sig1[100];
  Double_t sig2[100];
  Double_t sig3[100];

  for (i = 0; i < n; i++) {
    sig1[i] = (Double_t)source[i];
  }

  Float_t dt      = 0.1;
  Float_t lambda0 = (1.0 / fFeeParam->GetTFr2()) * dt;
  Float_t lambda1 = (1.0 / fFeeParam->GetTFr1()) * dt;

  FilterSimTailMakerSpline( sig1, sig2, lambda0, n);
  FilterSimTailCancelationMI( sig2, sig3, 0.7, lambda1, n);

  for (i = 0; i < n; i++) {
    target[i] = sig3[i];
  }

}

//______________________________________________________________________________
void AliTRDmcmSim::FilterSimTailMakerSpline(Double_t *ampin, Double_t *ampout
                                          , Double_t lambda, Int_t n) 
{
  //
  // Special filter (M. Ivanov)
  //

  Int_t    i = 0;
  Double_t l = TMath::Exp(-lambda*0.5);
  Double_t in[1000];
  Double_t out[1000];

  // Initialize in[] and out[] goes 0 ... 2*n+19
  for (i = 0; i < n*2+20; i++) {
    in[i] = out[i] = 0;
  }

  // in[] goes 0, 1
  in[0] = ampin[0];
  in[1] = (ampin[0] + ampin[1]) * 0.5;
   
  // Add charge to the end
  for (i = 0; i < 22; i++) {
    in[2*(n-1)+i] = ampin[n-1]; // in[] goes 2*n-2, 2*n-1, ... , 2*n+19 
  }

  // Use arithmetic mean
  for (i = 1; i < n-1; i++) {
    in[2*i]   = ampin[i];    // in[] goes 2, 3, ... , 2*n-4, 2*n-3
    in[2*i+1] = ((ampin[i]+ampin[i+1]))/2.;
  }

  Double_t temp;
  out[2*n]    = in[2*n];
  temp        = 0;
  for (i = 2*n; i >= 0; i--) {
    out[i]    = in[i] + temp;
    temp      = l*(temp+in[i]);
  }

  for (i = 0; i < n; i++){
    //ampout[i] = out[2*i+1];  // org
    ampout[i] = out[2*i];
  }

}

//______________________________________________________________________________
void AliTRDmcmSim::FilterSimTailCancelationMI(Double_t *ampin, Double_t *ampout
                                            , Double_t norm, Double_t lambda
                                            , Int_t n) 
{
  //
  // Special filter (M. Ivanov)
  //

  Int_t    i = 0;

  Double_t l = TMath::Exp(-lambda*0.5);
  Double_t k = l*(1.0 - norm*lambda*0.5);
  Double_t in[1000];
  Double_t out[1000];

  // Initialize in[] and out[] goes 0 ... 2*n+19
  for (i = 0; i < n*2+20; i++) {
    in[i] = out[i] = 0;
  }

  // in[] goes 0, 1
  in[0] = ampin[0];
  in[1] = (ampin[0]+ampin[1])*0.5;

  // Add charge to the end
  for (i =-2; i < 22; i++) {
    // in[] goes 2*n-4, 2*n-3, ... , 2*n+19 
    in[2*(n-1)+i] = ampin[n-1];
  }

  for (i = 1; i < n-2; i++) {
    // in[] goes 2, 3, ... , 2*n-6, 2*n-5
    in[2*i]    = ampin[i];
    in[2*i+1]  = (9.0 * (ampin[i]+ampin[i+1]) - (ampin[i-1]+ampin[i+2])) / 16.0;
    //in[2*i+1]  = ((ampin[i]+ampin[i+1]))/2.0;
  }

  Double_t temp;
  out[0] = in[0];
  temp   = in[0];
  for (i = 1; i <= 2*n; i++) {
    out[i] = in[i] + (k-l)*temp;
    temp   = in[i] +  k   *temp;
  }

  for (i = 0; i < n; i++) {
    //ampout[i] = out[2*i+1];  // org
    //ampout[i] = TMath::Max(out[2*i+1],0.0);  // org
    ampout[i] = TMath::Max(out[2*i],0.0);
  }
}


//_____________________________________________________________________________________
//the following filter uses AliTRDtrapAlu-class

void AliTRDmcmSim::FilterSimDeConvExpEl(Int_t *source, Int_t *target, Int_t n, Int_t nexp) {
  //static Int_t count = 0;
 
  Double_t dt = 0.1;
  Double_t r1 = (Double_t)fFeeParam->GetTFr1();
  Double_t r2 = (Double_t)fFeeParam->GetTFr2();
  Double_t c1 = (Double_t)fFeeParam->GetTFc1();
  Double_t c2 = (Double_t)fFeeParam->GetTFc2();
  
  nexp = 1;

  //it is assumed that r1,r2,c1,c2 are given such, that the configuration values are in the ranges according to TRAP-manual
  //parameters need to be adjusted
  AliTRDtrapAlu lambdaL;
  AliTRDtrapAlu lambdaS;
  AliTRDtrapAlu alphaL;
  AliTRDtrapAlu alphaS;
  
  AliTRDtrapAlu correction;
  AliTRDtrapAlu result;
  AliTRDtrapAlu bufL;
  AliTRDtrapAlu bufS;
 
  AliTRDtrapAlu bSource;
  
  lambdaL.Init(1,11);
  lambdaS.Init(1,11);
  alphaL.Init(1,11);
  alphaS.Init(1,11);
  
  //count=count+1;

  lambdaL.AssignDouble(TMath::Exp(-dt/r1));
  lambdaS.AssignDouble(TMath::Exp(-dt/r2));
  alphaL.AssignDouble(c1); // in AliTRDfeeParam the number of exponentials is set and also the according time constants
  alphaS.AssignDouble(c2); // later it should be: alphaS=1-alphaL
  
  //data is enlarged to 12 bits, including 2 bits after the comma; class AliTRDtrapAlu is used to handle arithmetics correctly
  correction.Init(10,2);
  result.Init(10,2);
  bufL.Init(10,2);
  bufS.Init(10,2);
  bSource.Init(10,2);
  
  for(Int_t i = 0; i < n; i++) {
    bSource.AssignInt(source[i]);
    result = bSource - correction; // subtraction can produce an underflow
    if(result.GetSign() == kTRUE) {
      result.AssignInt(0);
    }
    
    //target[i] = result.GetValuePre();  // later, target and source should become AliTRDtrapAlu,too in order to simulate the 10+2Bits through the filter properly
    
    target[i] = result.GetValue(); // 12 bit-value; to get the corresponding integer value, target must be shifted: target>>2 

    //printf("target-Wert zur Zeit %d : %d",i,target[i]);
    //printf("\n");
    
    bufL  =  bufL + (result * alphaL);
    bufL  =  bufL * lambdaL; 
    
    bufS  =  bufS + (result * alphaS);
    bufS  =  bufS * lambdaS;  // eventually this should look like:
    // bufS = (bufS + (result - result * alphaL)) * lambdaS // alphaS=1-alphaL; then alphaS-variable is not needed any more

    correction = bufL + bufS; //check for overflow intrinsic; if overflowed, correction is set to 0x03FF
  }
 
  
}







//__________________________________________________________________________________


// in order to use the Tracklets, please first 
// -- set AliTRDfeeParam::fgkTracklet to kTRUE, in order to switch on Tracklet-calculation
// -- set AliTRDfeeParam::fgkTFtype   to 3, in order to use the new tail cancellation filter
//    currently tracklets from filtered digits are only given when setting fgkTFtype (AliTRDfeeParam) to 3
// -- set AliTRDfeeParam::fgkMCTrackletOutput to kTRUE, if you want to use the Tracklet output container with 		information about the Tracklet position (MCM, channel number)

// The code is designed such that the less possible calculations with AliTRDtrapAlu class-objects are performed; whenever possible calculations are done with doubles or integers and the results are transformed into the right format

void AliTRDmcmSim::Tracklet(){
    // tracklet calculation
    // if you use this code after a simulation, please make sure the same filter-settings as in the simulation are set in AliTRDfeeParam

  if(!CheckInitialized()){ return; }
  
  Bool_t filtered = kTRUE;
  
  
  
  AliTRDtrapAlu data;
  data.Init(10,2);
  if(fADCT[0][0]==-1){                      // check if filter was applied
    filtered = kFALSE;
    for( Int_t iadc = 0 ; iadc < fNADC; iadc++ ) {
      for( Int_t iT = 0 ; iT < fNTimeBin ; iT++ ) {
	data.AssignInt(fADCR[iadc][iT]);
	fADCT[iadc][iT] = data.GetValue(); // all incoming values are positive 10+2 bit values; if el.filter was called, this is done correctly
      }
    }
   
  }
  
  // the online ordering of mcm's is reverse to the TRAP-manual-ordering! reverse fADCT (to be consistent to TRAP), then do all calculations
  // reverse fADCT:
  Int_t** rev0 = new Int_t *[fNADC];
  Int_t** rev1 = new Int_t *[fNADC];
  
  for( Int_t iadc = 0 ; iadc < fNADC; iadc++ ) {
    rev0[iadc] = new Int_t[fNTimeBin];
    rev1[iadc] = new Int_t[fNTimeBin];
    for( Int_t iT = 0; iT < fNTimeBin; iT++) {
      if( iadc <= fNADC-iadc-1 ) {
	rev0[iadc][iT]  = fADCT[fNADC-iadc-1][iT];
	rev1[iadc][iT]  = fADCT[iadc][iT];
	fADCT[iadc][iT] = rev0[iadc][iT];
      }
      else {
	rev0[iadc][iT]  = rev1[fNADC-iadc-1][iT];
	fADCT[iadc][iT] = rev0[iadc][iT];
      }
    }
  }
  for( Int_t iadc = 0 ; iadc < fNADC; iadc++ ) {
    delete[] rev0[iadc];
    delete[] rev1[iadc];
  }
  
  delete[] rev0;
  delete[] rev1;
  
  rev0 = NULL;
  rev1 = NULL;
    
  // get the filtered pedestal; supports only electronic tail-cancellation filter
  AliTRDtrapAlu filPed;
  Int_t ep = 0;
  Int_t *ieffped = new Int_t[fNTimeBin];
  for(Int_t iT = 0; iT < fNTimeBin; iT++){
    ieffped[iT] = ep; 
  }
 
  if( filtered == kTRUE ) {
    if( fFeeParam->IsPFon() ){
      ep = fFeeParam->GetPFeffectPedestal();
    }
    Int_t      nexp  = fFeeParam->GetTFnExp();
    Int_t  *isource  = new Int_t[fNTimeBin];
    filPed.Init(10,2);
    filPed.AssignInt(ep);           
    Int_t epf = filPed.GetValue();  
    for(Int_t iT = 0; iT < fNTimeBin; iT++){
      isource[iT] = ep;                  
      ieffped[iT] = epf;
    }
 
    if( fFeeParam->IsTFon() ) {
      FilterSimDeConvExpEl( isource, ieffped, fNTimeBin, nexp);
    }
  
    delete[] isource;
  }
  
  //the following values should go to AliTRDfeeParam once they are defined; then they have to be read in properly
  //naming follows conventions in TRAP-manual
  
  
  Bool_t bVBY = kTRUE;                         // cluster-verification bypass

  Double_t cQTParam = 0;                      // cluster quality threshold; granularity 2^-10; range: 0<=cQT/2^-10<=2^-4 - 2^-10
  AliTRDtrapAlu cQTAlu; 
  cQTAlu.Init(1,10,0,63);
  cQTAlu.AssignDouble(cQTParam);
  Int_t cQT = cQTAlu.GetValue();

  // linear fit 
  Int_t tFS = fFeeParam->GetLinearFitStart();  // linear fit start
  Int_t tFE = fFeeParam->GetLinearFitEnd();    // linear fit stop
   
  // charge accumulators
  Int_t tQS0 = fFeeParam->GetQacc0Start();     // start-time for charge-accumulator 0
  Int_t tQE0 = fFeeParam->GetQacc0End();       // stop-time for charge-accumulator 0
  Int_t tQS1 = fFeeParam->GetQacc1Start();     // start-time for charge-accumulator 1 
  Int_t tQE1 = fFeeParam->GetQacc1End();       // stop-time for charge-accumulator 1
  // values set such that tQS0=tFS; tQE0=tQS1-1; tFE=tQE1; want to do (QS0+QS1)/N
 
  Double_t cTHParam = (Double_t)fFeeParam->GetMinClusterCharge(); // cluster charge threshold
  AliTRDtrapAlu cTHAlu;  
  cTHAlu.Init(12,2);
  cTHAlu.AssignDouble(cTHParam);
  Int_t cTH = cTHAlu.GetValue();                                 // cTH used for comparison

  struct List_t {
    List_t *next;
    Int_t iadc;
    Int_t value;
  };
  
  List_t selection[7];            // list with 7 elements
  List_t *list = NULL;
  List_t *listLeft = NULL;
    
  Int_t* qsum = new Int_t[fNADC];
   
  // fit sums
  AliTRDtrapAlu qsumAlu;
  qsumAlu.Init(12,2);           // charge sum will be 12+2 bits
  AliTRDtrapAlu dCOGAlu; 
  dCOGAlu.Init(1,7,0,127);      // COG will be 1+7 Bits; maximum 1 - 2^-7 for LUT
  AliTRDtrapAlu yrawAlu;
  yrawAlu.Init(1,8,-1,255);
  AliTRDtrapAlu yAlu;
  yAlu.Init(1,16,-1,0xFF00);    // only first 8 past-comma bits filled;additional 8 bits for accuracy;maximum 1 - 2^-8; sign is given by + or -
  AliTRDtrapAlu xAlu;
  xAlu.Init(5,8);               // 8 past-comma bits because value will be added/multiplied to another value with this accuracy
  AliTRDtrapAlu xxAlu;
  xxAlu.Init(10,0);            
  AliTRDtrapAlu yyAlu;
  yyAlu.Init(1,16,0,0xFFFF);    // maximum is 2^16-1; 16Bit for past-commas
  AliTRDtrapAlu xyAlu;
  xyAlu.Init(6,8);
  AliTRDtrapAlu XAlu;
  XAlu.Init(9,0);
  AliTRDtrapAlu XXAlu;
  XXAlu.Init(14,0);
  AliTRDtrapAlu YAlu;
  YAlu.Init(5,8);               // 14 bit, 1 is sign-bit; therefore only 13 bit 
  AliTRDtrapAlu YYAlu;
  YYAlu.Init(5,16);
  AliTRDtrapAlu XYAlu;
  XYAlu.Init(8,8);              // 17 bit, 1 is sign-bit; therefore only 16 bit        
  AliTRDtrapAlu qtruncAlu;
  qtruncAlu.Init(12,0);
  AliTRDtrapAlu QT0Alu;
  QT0Alu.Init(15,0);
  AliTRDtrapAlu QT1Alu;
  QT1Alu.Init(16,0);

  AliTRDtrapAlu oneAlu;
  oneAlu.Init(1,8);
 
  
  AliTRDtrapAlu inverseNAlu;
  inverseNAlu.Init(1,8);        // simulates the LUT for 1/N
  AliTRDtrapAlu MeanChargeAlu;  // mean charge in ADC counts
  MeanChargeAlu.Init(8,0);
  AliTRDtrapAlu TotalChargeAlu;
  TotalChargeAlu.Init(17,8);
  //nr of post comma bits should be the same for inverseN and TotalCharge
  
  
  SetPosLUT();                    // initialize the position correction LUT for this MCM;


  // fit-sums; remapping!; 0,1,2->0; 1,2,3->1; ... 18,19,20->18
  Int_t *X   = new Int_t[fNADC-2];
  Int_t *XX  = new Int_t[fNADC-2];
  Int_t *Y   = new Int_t[fNADC-2];
  Int_t *YY  = new Int_t[fNADC-2];
  Int_t *XY  = new Int_t[fNADC-2];
  Int_t *N   = new Int_t[fNADC-2];
  Int_t *QT0 = new Int_t[fNADC-2]; // accumulated charge
  Int_t *QT1 = new Int_t[fNADC-2]; // accumulated charge
  
  for (Int_t iCol = 0; iCol < fNADC-2; iCol++) { 
      
      // initialize fit-sums 
      X[iCol]   = 0;
      XX[iCol]  = 0;
      Y[iCol]   = 0;
      YY[iCol]  = 0;
      XY[iCol]  = 0;
      N[iCol]   = 0;
      QT0[iCol] = 0;
      QT1[iCol] = 0;
  }
  

  filPed.Init(7,2);                         // convert filtered pedestal into 7+2Bits
  
  for(Int_t iT = 0; iT < fNTimeBin; iT++){
    
    if(iT<tFS || iT>=tFE) continue;         // linear fit yes/no? 

    // reset
    Int_t portChannel[4]   = {-1,-1,-1,-1};   
    Int_t clusterCharge[4] = {0,0,0,0};
    Int_t leftCharge[4]    = {0,0,0,0};
    Int_t centerCharge[4]  = {0,0,0,0}; 
    Int_t rightCharge[4]   = {0,0,0,0};
    
    Int_t mark = 0;
    
    filPed.AssignFormatted(ieffped[iT]);   // no size-checking when using AssignFormatted; ieffped>=0
    filPed = filPed;                       // this checks the size
    
    ieffped[iT] = filPed.GetValue();
        
    for(Int_t i = 0; i<7; i++){
      selection[i].next       = NULL;
      selection[i].iadc       =   -1;     // value of -1: invalid adc
      selection[i].value      =    0;
   
    }
    // selection[0] is starting list-element; just for pointing

    // loop over inner adc's 
    for (Int_t iCol = 1; iCol < fNADC-1; iCol++) { 
      
      Int_t left   = fADCT[iCol-1][iT]; 
      Int_t center = fADCT[iCol][iT];
      Int_t right  = fADCT[iCol+1][iT];  

      Int_t sum = left + center + right;            // cluster charge sum
      qsumAlu.AssignFormatted(sum);    
      qsumAlu = qsumAlu;                        // size-checking; redundant
 
      qsum[iCol] = qsumAlu.GetValue(); 
      
      //hit detection and masking
      if(center>=left){
	if(center>right){
	  if(qsum[iCol]>=(cTH + 3*ieffped[iT])){    // effective pedestal of all three channels must be added to cTH(+20); this is not parallel to TRAP manual; maybe cTH has to be adjusted in fFeeParam; therefore channels are not yet reduced by their pedestal
	    mark |= 1;                              // marker
	  }
	}
      }
      mark = mark<<1;                
    }
    mark = mark>>1;

       
    // get selection of 6 adc's and sort,starting with greatest values

    //read three from right side and sort (primitive sorting algorithm)
    Int_t i = 0; // adc number
    Int_t j = 1; // selection number
    while(i<fNADC-2 && j<=3){
      i = i + 1;
      if( ((mark>>(i-1)) & 1) == 1) {
	selection[j].iadc  = fNADC-1-i;
	selection[j].value = qsum[fNADC-1-i]>>6;   // for hit-selection only the first 8 out of the 14 Bits are used for comparison
	
	// insert into sorted list
	listLeft = &selection[0];
	list = listLeft->next;
	
	if(list!=NULL) {
	  while((list->next != NULL) && (selection[j].value <= list->value)){
	    listLeft = list;
	    list = list->next;
	  }
	  
	  if(selection[j].value<=list->value){
	    selection[j].next = list->next;
	    list->next = &selection[j];
	  }
	  else {
	    listLeft->next = &selection[j];
	    selection[j].next = list;
	  }
	}
	else{
	  listLeft->next = &selection[j];
	  selection[j].next = list;
	}
	
	j = j + 1;
      }
    }


    // read three from left side
    Int_t k = fNADC-2;
    while(k>i && j<=6) {
      if( ((mark>>(k-1)) & 1) == 1) {
	selection[j].iadc  = fNADC-1-k;
	selection[j].value = qsum[fNADC-1-k]>>6;
	
	listLeft = &selection[0];
	list = listLeft->next;
	
	if(list!=NULL){
	  while((list->next != NULL) && (selection[j].value <= list->value)){
	    listLeft = list;
	    list = list->next;
	  }
	
	  if(selection[j].value<=list->value){
	    selection[j].next = list->next;
	    list->next = &selection[j];
	  }
	  else {
	    listLeft->next = &selection[j];
	    selection[j].next = list;
	  }
	}
	else{
	  listLeft->next = &selection[j];
	  selection[j].next = list;
	}

	j = j + 1;
      }
      k = k - 1;
    }

    // get the four with greatest charge-sum
    list = &selection[0];
    for(i = 0; i<4; i++){
      if(list->next == NULL) continue;
      list = list->next;
      if(list->iadc == -1) continue;
      Int_t adc = list->iadc;                              // channel number with selected hit
      
      // the following arrays contain the four chosen channels in 1 time-bin
      portChannel[i]   = adc; 
      clusterCharge[i] = qsum[adc];
      leftCharge[i]    = fADCT[adc-1][iT] - ieffped[iT]; // reduce by filtered pedestal (pedestal is part of the signal)
      centerCharge[i]  = fADCT[adc][iT] - ieffped[iT];           
      rightCharge[i]   = fADCT[adc+1][iT] - ieffped[iT];         
    }

    // arithmetic unit
    
    // cluster verification
    if(!bVBY){
      for(i = 0; i<4; i++){
	Int_t lr = leftCharge[i]*rightCharge[i]*1024;
	Int_t cc = centerCharge[i]*centerCharge[i]*cQT;
	if (lr>=cc){
	  portChannel[i]   = -1;                                 // set to invalid address 
	  clusterCharge[i] = 0;
	}
      }
    }

    // fit-sums of valid channels
    // local hit position
    for(i = 0; i<4; i++){
      if (centerCharge[i] ==  0) {
	portChannel[i] = -1; 
      }// prevent division by 0
      
      if (portChannel[i]  == -1) continue;
      
      Double_t dCOG = (Double_t)(rightCharge[i]-leftCharge[i])/centerCharge[i];
       
      Bool_t sign = (dCOG>=0.0) ? kFALSE : kTRUE;
      dCOG = (sign == kFALSE) ? dCOG : -dCOG;     // AssignDouble doesn't allow for signed doubles
      dCOGAlu.AssignDouble(dCOG);
      Int_t iLUTpos = dCOGAlu.GetValue();       // steers position in LUT
            
      dCOG = dCOG/2;
      yrawAlu.AssignDouble(dCOG);
      Int_t iCOG = yrawAlu.GetValue();
      Int_t y = iCOG + fPosLUT[iLUTpos % 128];    // local position in pad-units
      yrawAlu.AssignFormatted(y);               // 0<y<1           
      yAlu  = yrawAlu;                        // convert to 16 past-comma bits
      
      if(sign == kTRUE) yAlu.SetSign(-1);       // buffer width of 9 bits; sign on real (not estimated) position
      xAlu.AssignInt(iT);                       // buffer width of 5 bits 
      

      xxAlu = xAlu * xAlu;                  // buffer width of 10 bits -> fulfilled by x*x       
      
      yyAlu = yAlu * yAlu;                  // buffer width of 16 bits
   
      xyAlu = xAlu * yAlu;                  // buffer width of 14 bits
                  
      Int_t adc = portChannel[i]-1;              // remapping! port-channel contains channel-nr. of inner adc's (1..19; mapped to 0..18)

      // calculate fit-sums recursively
      // interpretation of their bit-length is given as comment
      
      // be aware that the accuracy of the result of a calculation is always determined by the accuracy of the less accurate value

      XAlu.AssignFormatted(X[adc]);
      XAlu = XAlu + xAlu;                   // buffer width of 9 bits 
      X[adc] = XAlu.GetValue();
             
      XXAlu.AssignFormatted(XX[adc]);
      XXAlu = XXAlu + xxAlu;                // buffer width of 14 bits    
      XX[adc] = XXAlu.GetValue();

      if (Y[adc] < 0) {
	YAlu.AssignFormatted(-Y[adc]);          // make sure that only positive values are assigned; sign-setting must be done by hand
	YAlu.SetSign(-1);
      }
      else {
	YAlu.AssignFormatted(Y[adc]);
	YAlu.SetSign(1);
      }
	
      YAlu = YAlu + yAlu;                   // buffer width of 14 bits (8 past-comma);     
      Y[adc] = YAlu.GetSignedValue();
            
      YYAlu.AssignFormatted(YY[adc]);
      YYAlu = YYAlu + yyAlu;                // buffer width of 21 bits (16 past-comma) 
      YY[adc] = YYAlu.GetValue();
           
      if (XY[adc] < 0) {
	XYAlu.AssignFormatted(-XY[adc]);
	XYAlu.SetSign(-1);
      }
      else {
	XYAlu.AssignFormatted(XY[adc]);
	XYAlu.SetSign(1);
      }

      XYAlu = XYAlu + xyAlu;                // buffer allows 17 bits (8 past-comma) 
      XY[adc] = XYAlu.GetSignedValue();
            
      N[adc]  = N[adc] + 1;
   

      // accumulated charge
      qsumAlu.AssignFormatted(qsum[adc+1]); // qsum was not remapped!
      qtruncAlu = qsumAlu;

      if(iT>=tQS0 && iT<=tQE0){
	QT0Alu.AssignFormatted(QT0[adc]);
	QT0Alu = QT0Alu + qtruncAlu;
	QT0[adc] = QT0Alu.GetValue();
	//interpretation of QT0 as 12bit-value (all pre-comma); is this as it should be done?; buffer allows 15 Bit
      }
      
      if(iT>=tQS1 && iT<=tQE1){
	QT1Alu.AssignFormatted(QT1[adc]);
	QT1Alu = QT1Alu + qtruncAlu;
	QT1[adc] = QT1Alu.GetValue();
	//interpretation of QT1 as 12bit-value; buffer allows 16 Bit
      }
    }// i
      
    // remapping is done!!
     
  }//iT
 
  
    
  // tracklet-assembly
  
  // put into AliTRDfeeParam and take care that values are in proper range
  const Int_t cTCL = 1;      // left adc: number of hits; 8<=TCL<=31 (?? 1<=cTCL<+8 ??) 
  const Int_t cTCT = 8;      // joint number of hits;     8<=TCT<=31; note that according to TRAP manual this number cannot be lower than 8; however it should be adjustable to the number of hits in the fit time range (40%)
  
  Int_t mPair   = 0;         // marker for possible tracklet pairs
  Int_t* hitSum = new Int_t[fNADC-3];
  // hitSum[0] means: hit sum of remapped channels 0 and 1; hitSum[17]: 17 and 18; 
  
  // check for all possible tracklet-pairs of adjacent channels (two are merged); mark the left channel of the chosen pairs
  for (Int_t iCol = 0; iCol < fNADC-3; iCol++) {
    hitSum[iCol] = N[iCol] + N[iCol+1];
    if ((N[iCol]>=cTCL) && (hitSum[iCol]>=cTCT)) {
	mPair |= 1;         // mark as possible channel-pair
     
    }
    mPair = mPair<<1;
  }
  mPair = mPair>>1;
  
  List_t* selectPair = new List_t[fNADC-2];      // list with 18 elements (0..18) containing the left channel-nr and hit sums
                                                 // selectPair[18] is starting list-element just for pointing
  for(Int_t k = 0; k<fNADC-2; k++){
      selectPair[k].next       = NULL;
      selectPair[k].iadc       =   -1;           // invalid adc
      selectPair[k].value      =    0;
   
    }

 list = NULL;
 listLeft = NULL;
  
  // read marker and sort according to hit-sum
  
  Int_t adcL  = 0;            // left adc-channel-number (remapped)
  Int_t selNr = 0;            // current number in list
  
  // insert marked channels into list and sort according to hit-sum
  while(adcL < fNADC-3 && selNr < fNADC-3){
     
    if( ((mPair>>((fNADC-4)-(adcL))) & 1) == 1) {
      selectPair[selNr].iadc  = adcL;
      selectPair[selNr].value = hitSum[adcL];   
      
      listLeft = &selectPair[fNADC-3];
      list = listLeft->next;
	
      if(list!=NULL) {
	while((list->next != NULL) && (selectPair[selNr].value <= list->value)){
	  listLeft = list;
	  list = list->next;
	}
	
	if(selectPair[selNr].value <= list->value){
	  selectPair[selNr].next = list->next;
	  list->next = &selectPair[selNr];
	}
	else {
	  listLeft->next = &selectPair[selNr];
	  selectPair[selNr].next = list;
	}
	
      }
      else{
	listLeft->next = &selectPair[selNr];
	selectPair[selNr].next = list;
      }
      
      selNr = selNr + 1;
    }
    adcL = adcL + 1;
  }
  
  //select up to 4 channels with maximum number of hits
  Int_t lpairChannel[4] = {-1,-1,-1,-1}; // save the left channel-numbers of pairs with most hit-sum
  Int_t rpairChannel[4] = {-1,-1,-1,-1}; // save the right channel, too; needed for detecting double tracklets
  list = &selectPair[fNADC-3];
  
  for (Int_t i = 0; i<4; i++) {
    if(list->next == NULL) continue;
    list = list->next;
    if(list->iadc == -1) continue;
    lpairChannel[i] = list->iadc;        // channel number with selected hit
    rpairChannel[i] = lpairChannel[i]+1;
  }
  
  // avoid submission of double tracklets  
  for (Int_t i = 3; i>0; i--) {
    for (Int_t j = i-1; j>-1; j--) {
      if(lpairChannel[i] == rpairChannel[j]) {
	lpairChannel[i] = -1;
	rpairChannel[i] = -1;
	break;
      }
      /* if(rpairChannel[i] == lpairChannel[j]) {
	lpairChannel[i] = -1;
	rpairChannel[i] = -1;
	break;
	}*/
    }
  }
  
  // merging of the fit-sums of the remainig channels
  // assume same data-word-width as for fit-sums for 1 channel
  // relative scales!
  Int_t mADC[4];                      
  Int_t mN[4];
  Int_t mQT0[4];
  Int_t mQT1[4];
  Int_t mX[4];
  Int_t mXX[4];
  Int_t mY[4];
  Int_t mYY[4];
  Int_t mXY[4];
  Int_t mOffset[4];
  Int_t mSlope[4];
  Int_t mMeanCharge[4]; 
  Int_t inverseN = 0;
  Double_t invN = 0;
  Int_t one = 0;

  for (Int_t i = 0; i<4; i++){
    mADC[i] = -1;                        // set to invalid number
    mN[i]   =  0;
    mQT0[i] =  0;
    mQT1[i] =  0;
    mX[i]   =  0;
    mXX[i]  =  0;
    mY[i]   =  0;
    mYY[i]  =  0;
    mXY[i]  =  0;
    mOffset[i] = 0;
    mSlope[i]  = 0;
    mMeanCharge[i] = 0;
  }
  
  oneAlu.AssignInt(1);
  one = oneAlu.GetValue();              // one with 8 past comma bits
 
  for (Int_t i = 0; i<4; i++){
          

    mADC[i] = lpairChannel[i];          // mapping of merged sums to left channel nr. (0,1->0; 1,2->1; ... 17,18->17)
                                        // the adc and pad-mapping should now be one to one: adc i is linked to pad i; TRAP-numbering
    Int_t madc = mADC[i];
    if (madc == -1) continue;
    
    YAlu.AssignInt(N[rpairChannel[i]]);
    Int_t wpad  = YAlu.GetValue();       // enlarge hit counter of right channel by 8 past-comma bits; YAlu can have 5 pre-comma bits (values up to 63); hit counter<=nr of time bins (24)

    mN[i]    = hitSum[madc];
  
    // don't merge fit sums in case of a stand-alone tracklet (consisting of only 1 channel); in that case only left channel makes up the fit sums
    if (N[madc+1] == 0) {
	mQT0[i] = QT0[madc];
	mQT1[i] = QT1[madc];
	
    }
    else {

	// is it ok to do the size-checking for the merged fit-sums with the same format as for single-channel fit-sums?
	
	mQT0[i]   = QT0[madc] + QT0[madc+1];
	QT0Alu.AssignFormatted(mQT0[i]);   
	QT0Alu  = QT0Alu;                // size-check
	mQT0[i]   = QT0Alu.GetValue();     // write back
	
	mQT1[i]   = QT1[madc] + QT1[madc+1];
	QT1Alu.AssignFormatted(mQT1[i]);
	QT1Alu  = QT1Alu;
	mQT1[i]   = QT1Alu.GetValue();
    }
    
    // calculate the mean charge in adc values; later to be replaced by electron likelihood
    mMeanCharge[i] = mQT0[i] + mQT1[i]; // total charge
    mMeanCharge[i] = mMeanCharge[i]>>2; // losing of accuracy; accounts for high mean charge
    // simulate LUT for 1/N; LUT is fed with the double-accurate pre-calculated value of 1/N; accuracy of entries has to be adjusted to real TRAP
    invN = 1.0/(mN[i]);
    inverseNAlu.AssignDouble(invN);
    inverseN = inverseNAlu.GetValue();
    mMeanCharge[i] = mMeanCharge[i] * inverseN;  // now to be interpreted with 8 past-comma bits
    TotalChargeAlu.AssignFormatted(mMeanCharge[i]);
    TotalChargeAlu = TotalChargeAlu;
    MeanChargeAlu = TotalChargeAlu;
    mMeanCharge[i] = MeanChargeAlu.GetValue();
    
    // this check is not necessary; it is just for efficiency reasons
    if (N[madc+1] == 0) {
	mX[i]     =   X[madc];
	mXX[i]    =  XX[madc];
	mY[i]     =   Y[madc];
	mXY[i]    =  XY[madc];
	mYY[i]    =  YY[madc];
    }
    else {
	
	mX[i]     =   X[madc] +  X[madc+1];
	XAlu.AssignFormatted(mX[i]);
	XAlu      = XAlu;
	mX[i]     = XAlu.GetValue();
	
	mXX[i]    =  XX[madc] + XX[madc+1];
	XXAlu.AssignFormatted(mXX[i]);
	XXAlu     = XXAlu;
	mXX[i]    = XXAlu.GetValue();
 
    
	mY[i]     =   Y[madc] + Y[madc+1] + wpad;
	if (mY[i] < 0) {
	    YAlu.AssignFormatted(-mY[i]);
	    YAlu.SetSign(-1);
	}
	else {
	    YAlu.AssignFormatted(mY[i]);
	    YAlu.SetSign(1);
	}
	YAlu    = YAlu;
	mY[i]     = YAlu.GetSignedValue();
	
	mXY[i]    = XY[madc] + XY[madc+1] + X[madc+1]*one;    // multiplication by one to maintain the data format
	
	if (mXY[i] < 0) {
	    XYAlu.AssignFormatted(-mXY[i]);
	    XYAlu.SetSign(-1);
	}
	else {
	    XYAlu.AssignFormatted(mXY[i]);
	    XYAlu.SetSign(1);
	}
	XYAlu   = XYAlu;
	mXY[i]    = XYAlu.GetSignedValue();
    
	mYY[i]    = YY[madc] + YY[madc+1] + 2*Y[madc+1]*one+ wpad*one;
	if (mYY[i] < 0) {
	    YYAlu.AssignFormatted(-mYY[i]);
	    YYAlu.SetSign(-1);
	}
	else {
	    YYAlu.AssignFormatted(mYY[i]);
	    YYAlu.SetSign(1);
	}
	
	YYAlu   = YYAlu;
	mYY[i]    = YYAlu.GetSignedValue();
    }
  
  }
    
  // calculation of offset and slope from the merged fit-sums; 
  // YY is needed for some error measure only; still to be done
  // be aware that all values are relative values (scale: timebin-width; pad-width) and are integer values on special scale
  
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!important note: the offset is calculated from hits in the time bin range between tFS and tFE; it corresponds to the value at the height of the time bin tFS which does NOT need to correspond to the upper side of the drift   !!
  // !!volume (cathode wire plane). The offset cannot be rescaled as long as it is unknown which is the first time bin that contains hits from the drift region and thus to which distance from the cathode plane tFS corresponds.    !!
  // !!This has to be taken into account by the GTU. Furthermore a Lorentz correction might have to be applied to the offset (see below).                                                                                             !!
  // !!In this implementation it is assumed that no miscalibration containing changing drift velocities in the amplification region is used.                                                                                          !!
  // !!The corrections to the offset (e.g. no ExB correction applied as offset is supposed to be on top of drift region; however not at anode wire, so some inclination of drifting clusters due to Lorentz angle exists) are only    !!
  // !!valid (in approximation) if tFS is close to the beginning of the drift region.                                                                                                                                                 !!
  // !!The slope however can be converted to a deflection length between electrode and cathode wire plane as it is clear that the drift region is sampled 20 times                                                                    !!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  // which formats should be chosen?
  AliTRDtrapAlu denomAlu;
  denomAlu.Init(20,8);       
  AliTRDtrapAlu numAlu;
  numAlu.Init(20,8);     
  // is this enough pre-comma place? covers the range of the 13 bit-word of the transmitted offset
  // offset measured in coord. of left channel must be between -0.5 and 1.5; 14 pre-comma bits because numerator can be big

  for (Int_t i = 0; i<4; i++) {
    if (mADC[i] == -1) continue;
      
    Int_t num0  = (mN[i]*mXX[i]-mX[i]*mX[i]);
    if (num0 < 0) {
      denomAlu.AssignInt(-num0);    // num0 does not have to be interpreted as having past-comma bits -> AssignInt
      denomAlu.SetSign(-1);
    }
    else {
      denomAlu.AssignInt(num0);
      denomAlu.SetSign(1);
    }
    
    Int_t num1  = mN[i]*mXY[i] - mX[i]*mY[i];
    if (num1 < 0) {
      numAlu.AssignFormatted(-num1); // value of num1 is already formatted to have 8 past-comma bits
      numAlu.SetSign(-1);
    }
    else {
      numAlu.AssignFormatted(num1);
      numAlu.SetSign(1);
    }
    numAlu    = numAlu/denomAlu;
    mSlope[i]   = numAlu.GetSignedValue();
   
    Int_t num2  = mXX[i]*mY[i] - mX[i]*mXY[i];
   
    if (num2 < 0) {
      numAlu.AssignFormatted(-num2);
      numAlu.SetSign(-1);
    }
    else {
      numAlu.AssignFormatted(num2);
      numAlu.SetSign(1);
    }
   
    numAlu    = numAlu/denomAlu;
   
    
    mOffset[i]  = numAlu.GetSignedValue();
    numAlu.SetSign(1);
    denomAlu.SetSign(1);
       
                                 
    //numAlu.AssignInt(mADC[i]+1);   // according to TRAP-manual but trafo not to middle of chamber (0.5 channels away)             
    numAlu.AssignDouble((Double_t)mADC[i] + 1.5);      // numAlu has enough pre-comma place for that; correct trafo, best values
    mOffset[i]  = mOffset[i] + numAlu.GetValue();      // transform offset to a coord.system relative to chip; +1 to avoid neg. values 
    
    // up to here: adc-mapping according to TRAP manual and in line with pad-col mapping
    // reverse adc-counting to be again in line with the online mapping
    mADC[i]     = fNADC - 4 - mADC[i];                 // fNADC-4-mADC[i]: 0..17; remapping necessary;
    mADC[i]     = mADC[i] + 2; 
    // +2: mapping onto original ADC-online-counting: inner adc's corresponding to a chip's pasa: number 2..19
  }

  // adc-counting is corresponding to online mapping; use AliTRDfeeParam::GetPadColFromADC to get the pad to which adc is connected; 
  // pad-column mapping is reverse to adc-online mapping; TRAP adc-mapping is in line with pad-mapping (increase in same direction);
  
  // transform parameters to the local coordinate-system of a stack (used by GTU)
  AliTRDpadPlane* padPlane = fGeo->CreatePadPlane(fLayer,fStack);
  
  Double_t padWidthI = padPlane->GetWidthIPad()*10.0; // get values in cm; want them in mm
  //Double_t padWidthO = padPlane->GetWidthOPad()*10; // difference between outer pad-widths not included; in real TRAP??
  
  // difference between width of inner and outer pads of a row is not accounted for;
  
  Double_t magField = 0.4;                           // z-component of magnetic field in Tesla; adjust to current simulation!!; magnetic field can hardly be evaluated for the position of each mcm 
  Double_t eCharge  = 0.3;                           // unit charge in (GeV/c)/m*T
  Double_t ptMin   = 2.3;                            // minimum transverse momentum (GeV/c); to be adjusted(?)
  
  Double_t granularityOffset = 0.160;                // granularity for offset in mm
  Double_t granularitySlope  = 0.140;                // granularity for slope  in mm     
    
  // get the coordinates in SM-system; parameters: 
  
  Double_t zPos       =  (padPlane->GetRowPos(fRow))*10.0;  // z-position of the MCM; fRow is counted on a chamber; SM consists of 5 
  // zPos is position of pad-borders;
  Double_t zOffset = 0.0;
  if ( fRow == 0 || fRow == 15 ) {
      zOffset = padPlane->GetLengthOPad();
  }
  else {
      zOffset = padPlane->GetLengthIPad();
  }
  zOffset = (-1.0) * zOffset/2.0;
  // turn zPos to be z-coordinate at middle of pad-row
  zPos = zPos + zOffset;

      
  Double_t xPos       =  0.0;                               // x-position of the upper border of the drift-chamber of actual layer
  Int_t    icol       =  0;                                 // column-number of adc-channel
  Double_t yPos[4];                                         // y-position of the pad to which ADC is connected
  Double_t dx         = 30.0;                               // height of drift-chamber in mm; maybe retrieve from AliTRDGeometry
  Double_t freqSample = fFeeParam->GetSamplingFrequency();  // retrieve the sampling frequency (10.019750 MHz)
  Double_t vdrift     = fCal->GetVdriftAverage(fChaId);     // averaged drift velocity for this detector (1.500000 cm/us)
  Int_t    nrOfDriftTimeBins = Int_t(dx/10.0*freqSample/vdrift); // the number of time bins in the drift region (20)
  Int_t    nrOfAmplTimeBins  = 2;                           // the number of time bins between anode wire and cathode wires in ampl.region (3.5mm)(guess)(suppose v_drift+3.5cm/us there=>all clusters arrive at anode wire within one time bin (100ns))
  Int_t    nrOfOffsetCorrTimeBins = tFS - nrOfAmplTimeBins - 1; // -1 is  to be conservative; offset correction will not remove the shift but is supposed to improve it; if tFS = 5, 2 drift time bins before tFS are assumed
  if(nrOfOffsetCorrTimeBins < 0) nrOfOffsetCorrTimeBins = 0;// don't apply offset correction if no drift time bins before tFS can be assumed 
  Double_t lorTan     = AliTRDCommonParam::Instance()->GetOmegaTau(vdrift); // tan of the Lorentz-angle for this detector; could be evaluated and set as a parameter for each mcm
  //Double_t lorAngle   =  7.0;                             // Lorentz-angle in degrees
  Double_t tiltAngle  = padPlane->GetTiltingAngle();        // sign-respecting tilting angle of pads in actual layer
  Double_t tiltTan    = TMath::Tan(TMath::Pi()/180.0 * tiltAngle);
  //Double_t lorTan     = TMath::Tan(TMath::Pi()/180.0 * lorAngle);

  Double_t alphaMax[4];                            // maximum deflection from the direction to the primary vertex; granularity of hit pads
  Double_t slopeMin[4];                            // local limits for the deflection
  Double_t slopeMax[4];
  Int_t   mslopeMin[4];                            // in granularity units; to be compared to mSlope[i]
  Int_t   mslopeMax[4];


  // x coord. of upper side of drift chambers in local SM-system (in mm)
  // obtained by evaluating the x-range of the hits; should be crosschecked; only drift, not amplification region taken into account (30mm);
  // the y-deflection is given as difference of y between lower and upper side of drift-chamber, not pad-plane;
  switch(fLayer) 
    {
    case 0: 
      xPos = 3003.0;
      break;
    case 1:
      xPos = 3129.0;
      break;
    case 2:
      xPos = 3255.0;
      break;
    case 3:
      xPos = 3381.0;
      break;
    case 4:
      xPos = 3507.0;
      break;
    case 5:
      xPos = 3633.0;
      break;
    }
 
  // calculation of offset-correction n: 

  Int_t nCorrectOffset = (fRobPos % 2 == 0) ? ((fMcmPos % 4)) : ( 4 + (fMcmPos % 4));  
 
  nCorrectOffset = (nCorrectOffset - 4)*18 - 1;
  if (nCorrectOffset < 0) {
    numAlu.AssignInt(-nCorrectOffset);
    numAlu.SetSign(-1);
  }
  else {
    numAlu.AssignInt(nCorrectOffset);
    numAlu.SetSign(1);
  }
  nCorrectOffset = numAlu.GetSignedValue();   

  // the Lorentz correction to the offset
  Double_t lorCorrectOffset = lorTan *(Double_t)nrOfOffsetCorrTimeBins*vdrift*10.0/freqSample; // Lorentz offset correction in mm
  

  lorCorrectOffset = lorCorrectOffset/padWidthI; // Lorentz correction in pad width units
  
  if(lorCorrectOffset < 0) {
      numAlu.AssignDouble(-lorCorrectOffset);
      numAlu.SetSign(-1);
  }
  else{
      numAlu.AssignDouble(lorCorrectOffset);
      numAlu.SetSign(1);
  }
  
  Int_t mlorCorrectOffset = numAlu.GetSignedValue();
  
  
  Double_t mCorrectOffset = padWidthI/granularityOffset; // >= 0.0
 
  // calculation of slope-correction

  // this is only true for tracks coming (approx.) from primary vertex
  // everything is evaluated for a tracklet covering the whole drift chamber
  Double_t cCorrectSlope = (-lorTan*dx + zPos/xPos*dx*tiltTan)/granularitySlope;
  // Double_t cCorrectSlope =  zPos/xPos*dx*tiltTan/granularitySlope;
  // zPos can be negative! for track from primary vertex: zOut-zIn > 0 <=> zPos > 0
  
  if (cCorrectSlope < 0) {
      numAlu.AssignDouble(-cCorrectSlope);
      numAlu.SetSign(-1);
  }
  else {
      numAlu.AssignDouble(cCorrectSlope);
      numAlu.SetSign(1);
  }
  cCorrectSlope = numAlu.GetSignedValue();
 
  // convert slope to deflection between upper and lower drift-chamber position (slope is given in pad-unit/time-bins)
  // different pad-width of outer pads of a pad-plane not taken into account
  // note that the fit was only done in the range tFS to tFE, however this range does not need to cover the whole drift region (neither start nor end of it)
  // however the tracklets are supposed to be a fit in the drift region thus the linear function is stretched to fit the drift region of 30 mm
  
  
  Double_t mCorrectSlope = (Double_t)(nrOfDriftTimeBins)*padWidthI/granularitySlope;  // >= 0.0

  AliTRDtrapAlu correctAlu;
  correctAlu.Init(20,8);
  
  AliTRDtrapAlu offsetAlu;
  offsetAlu.Init(13,0,-0x1000,0x0FFF);          // 13 bit-word; 2-complement (1 sign-bit); asymmetric range
  
  AliTRDtrapAlu slopeAlu;
  slopeAlu.Init(7,0,-0x40,0x3F);                // 7 bit-word;  2-complement (1 sign-bit);

  for (Int_t i = 0; i<4; i++) {
    
    if (mADC[i] == -1) continue;
    
    icol = fFeeParam->GetPadColFromADC(fRobPos,fMcmPos,mADC[i]); // be aware that mADC[i] contains the ADC-number according to online-mapping
    yPos[i]   = (padPlane->GetColPos(icol))*10.0;
    
    
    // offset:
    
    correctAlu.AssignDouble(mCorrectOffset);     // done because max. accuracy is 8 bit
    mCorrectOffset = correctAlu.GetValueWhole(); // cut offset correction to 8 past-comma bit accuracy
    mOffset[i]  = (Int_t)((mCorrectOffset)*(Double_t)(mOffset[i] + nCorrectOffset - mlorCorrectOffset)); 
    //mOffset[i]  = mOffset[i]*(-1);                   // adjust to direction of y-axes in online simulation
    
    if (mOffset[i] < 0) {
      numAlu.AssignFormatted(-mOffset[i]);
      numAlu.SetSign(-1);
    }
    else {
      numAlu.AssignFormatted(mOffset[i]);
      numAlu.SetSign(1);
    }

    offsetAlu = numAlu; 
    mOffset[i]  = offsetAlu.GetSignedValue();  

    
    // slope:
    
    correctAlu.AssignDouble(mCorrectSlope);
    mCorrectSlope = correctAlu.GetValueWhole();
    
    mSlope[i]   = (Int_t)((mCorrectSlope*(Double_t)mSlope[i]) + cCorrectSlope);

    if (mSlope[i] < 0) {
      numAlu.AssignFormatted(-mSlope[i]);
      numAlu.SetSign(-1);
    }
    else {
      numAlu.AssignFormatted(mSlope[i]);
      numAlu.SetSign(1);
    }

    slopeAlu  = numAlu;     // here all past-comma values are cut, not rounded; alternatively add +0.5 before cutting (means rounding)
    mSlope[i]   = slopeAlu.GetSignedValue(); 
       
    // local (LTU) limits for the deflection 
    // ATan returns angles in radian
    alphaMax[i]  = TMath::ASin(eCharge*magField/(2.0*ptMin)*(TMath::Sqrt(xPos*xPos + yPos[i]*yPos[i]))/1000.0); // /1000: mm->m
    slopeMin[i]  = dx*(TMath::Tan(TMath::ATan(yPos[i]/xPos) - alphaMax[i]))/granularitySlope;
    slopeMax[i]  = dx*(TMath::Tan(TMath::ATan(yPos[i]/xPos) + alphaMax[i]))/granularitySlope;
    
    if (slopeMin[i] < 0) {
      slopeAlu.AssignDouble(-slopeMin[i]);
      slopeAlu.SetSign(-1);
    }
    else { 
      slopeAlu.AssignDouble(slopeMin[i]);
      slopeAlu.SetSign(1);
    }
    mslopeMin[i] = slopeAlu.GetSignedValue();  // the borders should lie inside the range of mSlope -> usage of slopeAlu again
   
    if (slopeMax[i] < 0) {
      slopeAlu.AssignDouble(-slopeMax[i]);
      slopeAlu.SetSign(-1);
    }
    else {
      slopeAlu.AssignDouble(slopeMax[i]);
      slopeAlu.SetSign(1);
    }
    mslopeMax[i] = slopeAlu.GetSignedValue();
  }

  // suppress submission of tracks with low stiffness
  // put parameters in 32bit-word and submit (write to file as root-file; sort after SM, stack, layer, chamber) 

  // sort tracklet-words in ascending y-order according to the offset (according to mADC would also be possible)
  // up to now they are sorted according to maximum hit sum
  // is the sorting really done in the TRAP-chip?
  
  Int_t order[4] = {-1,-1,-1,-1};
  Int_t wordnr = 0;   // number of tracklet-words
  
  for(Int_t j = 0; j < fMaxTracklets; j++) {
      //if( mADC[j] == -1) continue; 
      if( (mADC[j] == -1) || (mSlope[j] < mslopeMin[j]) || (mSlope[j] > mslopeMax[j])) continue; // this applies a pt-cut
      wordnr++;
      if( wordnr-1 == 0) {
	  order[0] = j;
	  continue;
      }
      // wordnr-1>0, wordnr-1<4
      order[wordnr-1] = j;
      for( Int_t k = 0; k < wordnr-1; k++) {
	  if( mOffset[j] < mOffset[order[k]] ) {
	      for( Int_t l = wordnr-1; l > k; l-- ) {
		  order[l] = order[l-1];
	      }
	      order[k] = j;
	      break;
	  }
	  
      }
  }
        
  // fill the bit-words in ascending order and without gaps
  UInt_t bitWord[4] = {0,0,0,0};                 // attention: unsigned int to have real 32 bits (no 2-complement)
  for(Int_t j = 0; j < wordnr; j++) { // only "wordnr" tracklet-words
      //Bool_t rem1 = kTRUE;
    
    Int_t i = order[j];
    //bit-word is 2-complement and therefore without sign
    bitWord[j] =   1; // this is the starting 1 of the bit-word (at 33rd position); the 1 must be ignored
    //printf("\n");
    UInt_t shift  = 0;
    UInt_t shift2 = 0;
	
	


    /*printf("mean charge: %d\n",mMeanCharge[i]);
    printf("row: %d\n",fRow);
    printf("slope: %d\n",mSlope[i]);
    printf("pad position: %d\n",mOffset[i]);
    printf("channel: %d\n",mADC[i]);*/

    // electron probability (currently not implemented; the mean charge is just scaled)
    shift = (UInt_t)mMeanCharge[i];
    for(Int_t iBit = 0; iBit < 8; iBit++) {
      bitWord[j]  = bitWord[j]<<1;
      bitWord[j] |= (shift>>(7-iBit))&1;               
      //printf("0");
    }

    // pad row
    shift = (UInt_t)fRow;
    for(Int_t iBit = 0; iBit < 4; iBit++) {
      bitWord[j]  = bitWord[j]<<1;
      bitWord[j] |= (shift>>(3-iBit))&1;
      //printf("%d", (fRow>>(3-iBit))&1);
    }
    
    // deflection length
    if(mSlope[i] < 0) {
	shift = (UInt_t)(-mSlope[i]);
	// shift2 is 2-complement of shift
	shift2 = 1;
	for(Int_t iBit = 1; iBit < 7; iBit++) {
	    shift2  = shift2<<1;
	    shift2 |= (1- (((shift)>>(6-iBit))&1) );
	    //printf("%d",(1-((-mSlope[i])>>(6-iBit))&1));
	}
	shift2 = shift2 + 1;
	//printf("1");
	for(Int_t iBit = 0; iBit < 7; iBit++) {
	    bitWord[j]  = bitWord[j]<<1;
	    bitWord[j] |= (shift2>>(6-iBit))&1;
	    //printf("%d",(1-((-mSlope[i])>>(6-iBit))&1));
	}
    }
    else {
	shift = (UInt_t)(mSlope[i]);
	bitWord[j]  = bitWord[j]<<1;
	bitWord[j]   |= 0;
	//printf("0");
	for(Int_t iBit = 1; iBit < 7; iBit++) {
	    bitWord[j]  = bitWord[j]<<1;
	    bitWord[j] |= (shift>>(6-iBit))&1;
	    //printf("%d",(mSlope[i]>>(6-iBit))&1);
	}
    }

    // pad position
    if(mOffset[i] < 0) {
	shift = (UInt_t)(-mOffset[i]);
	shift2 = 1;
	for(Int_t iBit = 1; iBit < 13; iBit++) {
	    shift2  = shift2<<1;
	    shift2 |= (1-(((shift)>>(12-iBit))&1));
	    //printf("%d",(1-((-mOffset[i])>>(12-iBit))&1));
	}
	shift2 = shift2 + 1;
	//printf("1");
	for(Int_t iBit = 0; iBit < 13; iBit++) {
	    bitWord[j]  = bitWord[j]<<1;
	    bitWord[j] |= (shift2>>(12-iBit))&1;
	    //printf("%d",(1-((-mSlope[i])>>(6-iBit))&1));
	}
    }
    else {
	shift = (UInt_t)mOffset[i];
	bitWord[j] = bitWord[j]<<1;
	bitWord[j]   |= 0; 
	//printf("0");
	for(Int_t iBit = 1; iBit < 13; iBit++) {
	    bitWord[j]  = bitWord[j]<<1;
	    bitWord[j] |= (shift>>(12-iBit))&1;
	    //printf("%d",(mOffset[i]>>(12-iBit))&1);
	}
    }


        
    //printf("bitWord: %u\n",bitWord[j]);
    //printf("adc: %d\n",mADC[i]);
    fMCMT[j] = bitWord[j];
  }
    
  //printf("\n");

  
  delete [] qsum;
  delete [] ieffped;

  delete [] X;
  delete [] XX;
  delete [] Y;
  delete [] YY;
  delete [] XY;
  delete [] N;
  delete [] QT0;
  delete [] QT1;

  delete [] hitSum;
  delete [] selectPair;

  delete padPlane;

//if you want to activate the MC tracklet output, set fgkMCTrackletOutput=kTRUE in AliTRDfeeParam
	
  if (!fFeeParam->GetMCTrackletOutput()) 
      return;
 
  AliLog::SetClassDebugLevel("AliTRDmcmSim", 10);
  AliLog::SetFileOutput("../log/tracklet.log");
  
  // testing for wordnr in order to speed up the simulation
  if (wordnr == 0) 
    return;
   
  UInt_t 	*trackletWord = new UInt_t[fMaxTracklets];
  Int_t 	*adcChannel   = new Int_t[fMaxTracklets];
  Int_t 	*trackRef     = new Int_t[fMaxTracklets];

  Int_t u = 0;

  AliTRDdigitsManager *digman = new AliTRDdigitsManager();
  digman->ReadDigits(AliRunLoader::Instance()->GetLoader("TRDLoader")->TreeD());
  digman->SetUseDictionaries(kTRUE);
  AliTRDfeeParam *feeParam = AliTRDfeeParam::Instance();

  for (Int_t j = 0; j < fMaxTracklets; j++) {
      Int_t i = order[j];
      trackletWord[j] = 0;
      adcChannel[j] = -1;
      if (bitWord[j]!=0) {
	  trackletWord[u] = bitWord[j];
	  adcChannel[u]   = mADC[i];   // mapping onto the original adc-array to be in line with the digits-adc-ordering (21 channels in total on 1 mcm, 18 belonging to pads); mADC[i] should be >-1 in case bitWord[i]>0

// Finding label of MC track
	  TH1F *hTrkRef = new TH1F("trackref", "trackref", 100000, 0, 100000);
	  Int_t track[3];
	  Int_t padcol = feeParam->GetPadColFromADC(fRobPos, fMcmPos, adcChannel[u]);
	  Int_t padcol_ngb = feeParam->GetPadColFromADC(fRobPos, fMcmPos, adcChannel[u] - 1);
	  Int_t padrow = 4 * (fRobPos / 2) + fMcmPos / 4;
	  Int_t det = 30 * fSector + 6 * fStack + fLayer;
	  for(Int_t iTimebin = feeParam->GetLinearFitStart(); iTimebin < feeParam->GetLinearFitEnd(); iTimebin++) {
	      track[0] = digman->GetTrack(0, padrow, padcol, iTimebin, det);
	      track[1] = digman->GetTrack(1, padrow, padcol, iTimebin, det);
	      track[2] = digman->GetTrack(2, padrow, padcol, iTimebin, det);
	      hTrkRef->Fill(track[0]);
	      if (track[1] != track[0] && track[1] != -1)
		  hTrkRef->Fill(track[1]);
	      if (track[2] != track[0] && track[2] != track[1] && track[2] != -1)
		  hTrkRef->Fill(track[2]);
	      if (padcol_ngb >= 0) {
		  track[0] = digman->GetTrack(0, padrow, padcol, iTimebin, det);
		  track[1] = digman->GetTrack(1, padrow, padcol, iTimebin, det);
		  track[2] = digman->GetTrack(2, padrow, padcol, iTimebin, det);
		  hTrkRef->Fill(track[0]);
		  if (track[1] != track[0] && track[1] != -1)
		      hTrkRef->Fill(track[1]);
		  if (track[2] != track[0] && track[2] != track[1] && track[2] != -1)
		      hTrkRef->Fill(track[2]);
	      }
	  }
	  trackRef[u] = hTrkRef->GetMaximumBin() - 1;
	  delete hTrkRef;
	  u = u + 1;
      }
  }

  AliDataLoader *dl = AliRunLoader::Instance()->GetLoader("TRDLoader")->GetDataLoader("tracklets");
  if (!dl) {
    AliError("Could not get the tracklets data loader!");
  }
  else {
    TTree *trackletTree = dl->Tree();
    if (!trackletTree)
      dl->MakeTree();
    trackletTree = dl->Tree();

   AliTRDtrackletMCM *trkl = new AliTRDtrackletMCM(); 
   TBranch *trkbranch = trackletTree->GetBranch("mcmtrklbranch");
   if (!trkbranch)
       trkbranch = trackletTree->Branch("mcmtrklbranch", "AliTRDtrackletMCM", &trkl, 32000);
    trkbranch->SetAddress(&trkl);

    for (Int_t iTracklet = 0; iTracklet < fMaxTracklets; iTracklet++) {
	if (trackletWord[iTracklet] == 0)
	    continue;
	trkl->SetTrackletWord(trackletWord[iTracklet]);
	trkl->SetDetector(30*fSector + 6*fStack + fLayer);
	trkl->SetROB(fRobPos);
	trkl->SetMCM(fMcmPos);
	trkl->SetLabel(trackRef[iTracklet]);
	trackletTree->Fill();
    }
    delete trkl;
    dl->WriteData("OVERWRITE");
  }

  delete [] trackletWord;
  delete [] adcChannel; 
  delete [] trackRef;
  delete digman;

  // to be done:
  // error measure for quality of fit (not necessarily needed for the trigger)
  // cluster quality threshold (not yet set)
  // electron probability
}
//_____________________________________________________________________________________
void AliTRDmcmSim::GeneratefZSM1Dim()
{
  //
  // Generate the array fZSM1Dim necessary
  // for the method ProduceRawStream
  //

  // Fill the mapping
  // Supressed zeros indicated by -1 in digits array
  for( Int_t iadc = 1 ; iadc < fNADC-1; iadc++ ) 
    {
      for( Int_t it = 0 ; it < fNTimeBin ; it++ ) 
	{
	  
	  if(fADCF[iadc][it]==-1)  // If is a supressed value
	    {
	      fZSM[iadc][it]=1;
	    }
	  else                    // Not suppressed
	    {
	      fZSM[iadc][it]=0;
	    }
	}
    }

  // Make the 1 dim projection
  for( Int_t iadc = 0 ; iadc < fNADC; iadc++ ) 
    {
      for( Int_t it = 0 ; it < fNTimeBin ; it++ ) 
	{
	  fZSM1Dim[iadc] &= fZSM[iadc][it];
	}
    }
}
//_______________________________________________________________________________________
void AliTRDmcmSim::CopyArrays()
{
  //
  // Initialize filtered data array with raw data
  // Method added for internal consistency
  //

  for( Int_t iadc = 0 ; iadc < fNADC; iadc++ ) 
    {
      for( Int_t it = 0 ; it < fNTimeBin ; it++ ) 
	{
	  fADCF[iadc][it] = fADCR[iadc][it]; 
	}
    }
}
//_______________________________________________________________________________________
void AliTRDmcmSim::StartfastZS(Int_t pads, Int_t timebins)
{
  //
  // Initialize just the necessary elements to perform
  // the zero suppression in the digitizer
  //
   
  fFeeParam  = AliTRDfeeParam::Instance();
  fSimParam  = AliTRDSimParam::Instance();
  fNADC      = pads;      
  fNTimeBin  = timebins; 

  if( fADCR == NULL ) 
    {
      fADCR    = new Int_t *[fNADC];
      fADCF    = new Int_t *[fNADC];
      fADCT    = new Int_t *[fNADC]; 
      fZSM     = new Int_t *[fNADC];
      fZSM1Dim = new Int_t  [fNADC];
    for( Int_t iadc = 0 ; iadc < fNADC; iadc++ )
      {
	fADCR[iadc] = new Int_t[fNTimeBin];
	fADCF[iadc] = new Int_t[fNTimeBin];
	fADCT[iadc] = new Int_t[fNTimeBin]; 
	fZSM [iadc] = new Int_t[fNTimeBin];
      }
    }

  for( Int_t iadc = 0 ; iadc < fNADC; iadc++ ) 
    {
      for( Int_t it = 0 ; it < fNTimeBin ; it++ ) 
	{
	  fADCR[iadc][it] =  0;
	  fADCF[iadc][it] =  0;
	  fADCT[iadc][it] = -1;  
	  fZSM [iadc][it] =  1;   
	}
      fZSM1Dim[iadc] = 1;      
    }
  
  fInitialized = kTRUE;
}
//_______________________________________________________________________________________
void AliTRDmcmSim::FlagDigitsArray(AliTRDarrayADC *tempdigs, Int_t valrow)
{
  //
  // Modify the digits array to flag suppressed values
  //

  for( Int_t iadc = 1 ; iadc < fNADC-1; iadc++ ) 
    {
      for( Int_t it = 0 ; it < fNTimeBin ; it++ ) 
	{
	  if(fZSM[iadc][it]==1)
	    {
	      tempdigs->SetData(valrow,iadc,it,-1);
	    }
	}
    }
}
//_______________________________________________________________________________________
void AliTRDmcmSim::RestoreZeros()
{
  //
  // Restore the zero-suppressed values (set as -1) to the value 0
  //

  for( Int_t iadc = 1 ; iadc < fNADC-1; iadc++ ) 
    {
      for( Int_t it = 0 ; it < fNTimeBin ; it++ ) 
	{
	  
	  if(fADCF[iadc][it]==-1)  //if is a supressed zero, reset to zero
	    {
	      fADCF[iadc][it]=0;
	      fADCR[iadc][it]=0;
	    }	  
	}
    }

}

