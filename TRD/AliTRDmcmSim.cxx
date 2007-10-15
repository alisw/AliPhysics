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

#include <fstream>

#include <TMath.h>

#include "AliLog.h"

#include "AliTRDmcmSim.h"
#include "AliTRDfeeParam.h"
#include "AliTRDSimParam.h"
#include "AliTRDgeometry.h"
#include "AliTRDcalibDB.h"

ClassImp(AliTRDmcmSim)

//_____________________________________________________________________________
AliTRDmcmSim::AliTRDmcmSim() :TObject()
  ,fInitialized(kFALSE)
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
      delete [] fZSM [iadc];
    }
    delete [] fADCR;
    delete [] fADCF;
    delete [] fZSM;
    delete [] fZSM1Dim;
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
  ((AliTRDmcmSim &) m).fZSM           = 0;
  ((AliTRDmcmSim &) m).fZSM1Dim       = 0;
  ((AliTRDmcmSim &) m).fFeeParam      = 0;
  ((AliTRDmcmSim &) m).fSimParam      = 0;
  ((AliTRDmcmSim &) m).fCal           = 0;
  ((AliTRDmcmSim &) m).fGeo           = 0;

}

//_____________________________________________________________________________
void AliTRDmcmSim::Init( Int_t chaId, Int_t robPos, Int_t mcmPos )
{
  //
  // Initialize the class with new geometry information
  // fADC array will be reused with filled by zero
  //

  fFeeParam      = AliTRDfeeParam::Instance();
  fSimParam      = AliTRDSimParam::Instance();
  fCal           = AliTRDcalibDB::Instance();
  fGeo           = new AliTRDgeometry();
  fChaId         = chaId;
  fSector        = fGeo->GetSector( fChaId );
  fStack         = fGeo->GetChamber( fChaId );
  fLayer         = fGeo->GetPlane( fChaId );
  fRobPos        = robPos;
  fMcmPos        = mcmPos;
  fNADC          = fFeeParam->GetNadcMcm();
  fNTimeBin      = fCal->GetNumberOfTimeBins();
  fRow           = fFeeParam->GetPadRowFromMCM( fRobPos, fMcmPos );

  // Allocate ADC data memory if not yet done
  if( fADCR == NULL ) {
    fADCR    = new Int_t *[fNADC];
    fADCF    = new Int_t *[fNADC];
    fZSM     = new Int_t *[fNADC];
    fZSM1Dim = new Int_t  [fNADC];
    for( Int_t iadc = 0 ; iadc < fNADC; iadc++ ) {
      fADCR[iadc] = new Int_t[fNTimeBin];
      fADCF[iadc] = new Int_t[fNTimeBin];
      fZSM [iadc] = new Int_t[fNTimeBin];
    }
  }

  // Initialize ADC data
  for( Int_t iadc = 0 ; iadc < fNADC; iadc++ ) {
    for( Int_t it = 0 ; it < fNTimeBin ; it++ ) {
      fADCR[iadc][it] = 0;
      fADCF[iadc][it] = 0;
      fZSM [iadc][it] = 1;   // Default unread = 1
    }
    fZSM1Dim[iadc] = 1;      // Default unread = 1
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
  if( fFeeParam->isPFon() ) FilterPedestal();
  if( fFeeParam->isGFon() ) FilterGain();
  if( fFeeParam->isTFon() ) FilterTail();
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
    
  default:
    AliError(Form("Invalid filter type %d ! \n", tftype ));
    break;
  }

  delete dtarg;
  delete itarg;
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

