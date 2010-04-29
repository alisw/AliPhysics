
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

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  TRD front end electronics parameters class                            //
//  Contains all FEE (MCM, TRAP, PASA) related                            //
//  parameters, constants, and mapping.                                   //
//                                                                        //
//  New release on 2007/08/17:                                            //
//   The default raw data version (now fRAWversion ) is set to 3          //
//   in the constructor because version 3 raw data read and write         //
//   are fully debugged.                                                  //
//                                                                        //
//  Author:                                                               //
//    Ken Oyama (oyama@physi.uni-heidelberg.de)                           //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

//#include <TMath.h>

#include "AliLog.h"

#include "AliTRDfeeParam.h"
//#include "AliTRDgeometry.h"
#include "AliTRDCommonParam.h"

ClassImp(AliTRDfeeParam)

AliTRDfeeParam *AliTRDfeeParam::fgInstance   = 0;
Bool_t          AliTRDfeeParam::fgTerminated = kFALSE;
Bool_t          AliTRDfeeParam::fgTracklet = kFALSE;

//_____________________________________________________________________________
AliTRDfeeParam* AliTRDfeeParam::Instance()
{
  //
  // Instance constructor
  //

  if (fgTerminated != kFALSE) {
    return 0;
  }

  if (fgInstance == 0) {
    fgInstance = new AliTRDfeeParam();
  }  

  return fgInstance;

}

//_____________________________________________________________________________
void AliTRDfeeParam::Terminate()
{
  //
  // Terminate the class and release memory
  //
  
  fgTerminated = kTRUE;

  if (fgInstance != 0) {
    delete fgInstance;
    fgInstance = 0;
  }

}

//_____________________________________________________________________________
AliTRDfeeParam::AliTRDfeeParam()
  :TObject()
  ,fCP(0)
  ,fTFnExp(1)
  ,fTFr1(0)
  ,fTFr2(0)
  ,fTFc1(0)
  ,fTFc2(0)
  ,fEBsglIndThr(5)
  ,fEBsumIndThr(5)
  ,fEBindLUT(0xF0)
  ,fEBignoreNeighbour(0)
  ,fRAWversion(3)
  ,fRAWstoreRaw(kTRUE)
{
  //
  // Default constructor
  //
  
  // PASA V.4
  if      (fTFnExp == 1) {
    fTFr1 = 1.1563;
    fTFr2 = 0.1299;
    fTFc1 = 0.0657;
    fTFc2 = 0.0000;
  }
  else if (fTFnExp == 2) {
    fTFr1 = 1.1563;
    fTFr2 = 0.1299;
    fTFc1 = 0.1141;
    fTFc2 = 0.6241;
  }

  fCP  = AliTRDCommonParam::Instance();

}

//_____________________________________________________________________________
AliTRDfeeParam::AliTRDfeeParam(TRootIoCtor *)
  :TObject()
  ,fCP(0)
  ,fTFnExp(1)
  ,fTFr1(0)
  ,fTFr2(0)
  ,fTFc1(0)
  ,fTFc2(0)
  ,fEBsglIndThr(0)
  ,fEBsumIndThr(0)
  ,fEBindLUT(0)
  ,fEBignoreNeighbour(0)
  ,fRAWversion(0)
  ,fRAWstoreRaw(0)
{
  //
  // IO constructor
  //

}

//_____________________________________________________________________________
AliTRDfeeParam::AliTRDfeeParam(const AliTRDfeeParam &p)
  :TObject(p)
  ,fCP(p.fCP)
  ,fTFnExp(p.fTFnExp)
  ,fTFr1(p.fTFr1)
  ,fTFr2(p.fTFr2)
  ,fTFc1(p.fTFc1)
  ,fTFc2(p.fTFc2)
  ,fEBsglIndThr(p.fEBsglIndThr)
  ,fEBsumIndThr(p.fEBsumIndThr)
  ,fEBindLUT(p.fEBindLUT)
  ,fEBignoreNeighbour (p.fEBignoreNeighbour)
  ,fRAWversion(p.fRAWversion)
  ,fRAWstoreRaw(p.fRAWstoreRaw)
{
  //
  // AliTRDfeeParam copy constructor
  //

}

//_____________________________________________________________________________
AliTRDfeeParam::~AliTRDfeeParam()
{
  //
  // AliTRDfeeParam destructor
  //

}

//_____________________________________________________________________________
AliTRDfeeParam &AliTRDfeeParam::operator=(const AliTRDfeeParam &p)
{
  //
  // Assignment operator
  //

  if (this != &p) {
    ((AliTRDfeeParam &) p).Copy(*this);
  }

  return *this;

}

//_____________________________________________________________________________
void AliTRDfeeParam::Copy(TObject &p) const
{
  //
  // Copy function
  //

  ((AliTRDfeeParam &) p).fCP          = fCP;
  ((AliTRDfeeParam &) p).fTFnExp      = fTFnExp;
  ((AliTRDfeeParam &) p).fTFr1        = fTFr1;
  ((AliTRDfeeParam &) p).fTFr2        = fTFr2;
  ((AliTRDfeeParam &) p).fTFc1        = fTFc1;
  ((AliTRDfeeParam &) p).fTFc2        = fTFc2;
  ((AliTRDfeeParam &) p).fEBsglIndThr = fEBsglIndThr;
  ((AliTRDfeeParam &) p).fEBsumIndThr = fEBsumIndThr;
  ((AliTRDfeeParam &) p).fEBindLUT    = fEBindLUT;
  ((AliTRDfeeParam &) p).fEBignoreNeighbour = fEBignoreNeighbour;
  ((AliTRDfeeParam &) p).fRAWversion  = fRAWversion;
  ((AliTRDfeeParam &) p).fRAWstoreRaw = fRAWstoreRaw;
  
  TObject::Copy(p);

}

//_____________________________________________________________________________
Int_t AliTRDfeeParam::GetPadRowFromMCM(Int_t irob, Int_t imcm) const
{
  //
  // Return on which pad row this mcm sits
  //
  
  return fgkNmcmRobInRow*(irob/2) + imcm/fgkNmcmRobInCol;

}

//_____________________________________________________________________________
Int_t AliTRDfeeParam::GetPadColFromADC(Int_t irob, Int_t imcm, Int_t iadc) const
{
  //
  // Return which pad is connected to this adc channel.
  //
  // Return virtual pad number even if ADC is outside chamber
  // to keep compatibility of data processing at the edge MCM.
  // User has to check that this is in the chamber if it is essential.
  // Return -100 if iadc is invalid.
  //
  // Caution: ADC ordering in the online data is opposite to the pad column ordering.
  // And it is not one-by-one correspondence. Precise drawing can be found in:
  // http://wiki.kip.uni-heidelberg.de/ti/TRD/index.php/Image:ROB_MCM_numbering.pdf
  //

  if (iadc < 0 || iadc > fgkNadcMcm ) return -100;
  Int_t mcmcol = imcm%fgkNmcmRobInCol + GetRobSide(irob)*fgkNmcmRobInCol;  // MCM column number on ROC [0..7]
  Int_t padcol = mcmcol*fgkNcolMcm + fgkNcolMcm + 1 - iadc;
  if( padcol < 0 || padcol >= fgkNcol ) return -1;   // this is commented because of reason above OK

  return padcol;

}

//_____________________________________________________________________________
Int_t AliTRDfeeParam::GetExtendedPadColFromADC(Int_t irob, Int_t imcm, Int_t iadc) const
{     
  //
  // Return which pad coresponds to the extended digit container pad numbering
  // Extended digit container is designed to store all pad data including shared pad, 
  // so we have to introduce new virtual pad numbering scheme for this purpose. 
  //
    
  if (iadc < 0 || iadc > fgkNadcMcm ) return -100;
  Int_t mcmcol = imcm%fgkNmcmRobInCol + GetRobSide(irob)*fgkNmcmRobInCol;  // MCM column number on ROC [0..7]
  Int_t padcol = mcmcol*fgkNadcMcm + fgkNcolMcm + 2 - iadc;

  return padcol;

}

//_____________________________________________________________________________
Int_t AliTRDfeeParam::GetMCMfromPad(Int_t irow, Int_t icol) const
{
  //
  // Return on which MCM this pad is directry connected.
  // Return -1 for error.
  //

  if ( irow < 0 || icol < 0 || irow > fgkNrowC1 || icol > fgkNcol ) return -1;

  return (icol%(fgkNcol/2))/fgkNcolMcm + fgkNmcmRobInCol*(irow%fgkNmcmRobInRow);

}

//_____________________________________________________________________________
Int_t AliTRDfeeParam::GetMCMfromSharedPad(Int_t irow, Int_t icol) const
{
  //
  // Return on which MCM this pad is directry connected.
  // Return -1 for error.
  //
  
  if ( irow < 0 || icol < 0 || irow > fgkNrowC1 || icol > fgkNcol+8*3 ) return -1;

  Int_t adc = 20 - (icol%18) -1;
  switch(adc) {
    case 2:  icol += 5; break;
    case 18: icol -= 5; break;
    case 19: icol -= 5; break;
    default: icol += 0; break;
  }

  return (icol%(fgkNcol/2))/fgkNcolMcm + fgkNmcmRobInCol*(irow%fgkNmcmRobInRow);

}

//_____________________________________________________________________________
Int_t AliTRDfeeParam::GetROBfromPad(Int_t irow, Int_t icol) const
{
  //
  // Return on which rob this pad is
  //

  return (irow/fgkNmcmRobInRow)*2 + GetColSide(icol);

}

//_____________________________________________________________________________
Int_t AliTRDfeeParam::GetROBfromSharedPad(Int_t irow, Int_t icol) const
{
  //
  // Return on which rob this pad is for shared pads
  //

  if(icol<72) return (irow/fgkNmcmRobInRow)*2 + GetColSide(icol+5);
  else return (irow/fgkNmcmRobInRow)*2 + GetColSide(icol-5);

}

//_____________________________________________________________________________
Int_t AliTRDfeeParam::GetRobSide(Int_t irob) const
{
  //
  // Return on which side this rob sits (A side = 0, B side = 1)
  //

  if ( irob < 0 || irob >= fgkNrobC1 ) return -1;

  return irob%2;

}

//_____________________________________________________________________________
Int_t AliTRDfeeParam::GetColSide(Int_t icol) const
{
  //
  // Return on which side this column sits (A side = 0, B side = 1)
  //

  if ( icol < 0 || icol >= fgkNcol ) return -1;

  return icol/(fgkNcol/2);

}



UInt_t AliTRDfeeParam::AliToExtAli(Int_t rob, Int_t aliid)
{
   if(aliid!= 127)
      return ( (1 << 10) | (rob << 7) | aliid);

   return 127;
}


Int_t AliTRDfeeParam::ExtAliToAli(UInt_t dest, UShort_t linkpair, UShort_t rocType, Int_t *mcmList, Int_t listSize)
{
   // Converts an extended ALICE ID which identifies a single MCM or a group of MCMs to
   // the corresponding list of MCMs. Only broadcasts (127) are encoded as 127 
   // The return value is the number of MCMs in the list

  mcmList[0]=-1;

  Short_t nmcm = 0;
  UInt_t mcm, rob, robAB;
  UInt_t cmA = 0, cmB = 0;  // Chipmask for each A and B side
  
  // Default chipmask for 4 linkpairs (each bit correponds each alice-mcm)
  static const UInt_t gkChipmaskDefLp[4] = { 0x1FFFF, 0x1FFFF, 0x3FFFF, 0x1FFFF };
  
  rob = dest >> 7;                              // Extract ROB pattern from dest.
  mcm = dest & 0x07F;                           // Extract MCM pattern from dest.
  robAB = GetRobAB( rob, linkpair ); // Get which ROB sides are selected.
  
  // Abort if no ROB is selected
  if( robAB == 0 ) {
    return 0;
  }
  
  // Special case
  if( mcm == 127 ) {
    if( robAB == 3 ) {      // This is very special 127 can stay only if two ROBs are selected
      mcmList[0]=127;      // broadcase to ALL
      mcmList[1]=-1;
      return 1;
    }
    cmA = cmB = 0x3FFFF;
  } else if( (mcm & 0x40) != 0 ) { // If top bit is 1 but not 127, this is chip group.
    if( (mcm & 0x01) != 0 )                  { cmA |= 0x04444; cmB |= 0x04444; } // chip_cmrg
    if( (mcm & 0x02) != 0 )                  { cmA |= 0x10000; cmB |= 0x10000; } // chip_bmrg
    if( (mcm & 0x04) != 0 && rocType == 0 ) { cmA |= 0x20000; cmB |= 0x20000; } // chip_hm3
    if( (mcm & 0x08) != 0 && rocType == 1 ) { cmA |= 0x20000; cmB |= 0x20000; } // chip_hm4
    if( (mcm & 0x10) != 0 )                  { cmA |= 0x01111; cmB |= 0x08888; } // chip_edge
    if( (mcm & 0x20) != 0 )                  { cmA |= 0x0aaaa; cmB |= 0x03333; } // chip_norm
  } else { // Otherwise, this is normal chip ID, turn on only one chip.
    cmA = 1 << mcm;
    cmB = 1 << mcm;
  }
  
  // Mask non-existing MCMs
  cmA &= gkChipmaskDefLp[linkpair];
  cmB &= gkChipmaskDefLp[linkpair];
  // Remove if only one side is selected
  if( robAB == 1 ) 
    cmB = 0;
  if( robAB == 2 ) 
    cmA = 0;
  if( robAB == 4 && linkpair != 2 ) 
    cmA = cmB = 0; // Restrict to only T3A and T3B
  
  // Finally convert chipmask to list of slaves
  nmcm = ChipmaskToMCMlist( cmA, cmB, linkpair, mcmList, listSize);
  
  return nmcm;
}


Short_t AliTRDfeeParam::GetRobAB( UShort_t robsel, UShort_t linkpair )
{
  // Converts the ROB part of the extended ALICE ID to robs

  if( (robsel & 0x8) != 0 ) { // 1000 .. direct ROB selection. Only one of the 8 ROBs are used.
    robsel = robsel & 7;
    if( (robsel % 2) == 0 && (robsel / 2) == linkpair ) 
      return 1;  // Even means A side (position 0,2,4,6)
    if( (robsel % 2) == 1 && (robsel / 2) == linkpair ) 
      return 2;  // Odd  means B side (position 1,3,5,7)
    return 0;
  }
  
  // ROB group
  if( robsel == 0 ) { return 3; } // Both   ROB
  if( robsel == 1 ) { return 1; } // A-side ROB
  if( robsel == 2 ) { return 2; } // B-side ROB
  if( robsel == 3 ) { return 3; } // Both   ROB
  if( robsel == 4 ) { return 4; } // Only T3A and T3B
  // Other number 5 to 7 are ignored (not defined) 
  
  return 0;
}


Short_t AliTRDfeeParam::ChipmaskToMCMlist( UInt_t cmA, UInt_t cmB, UShort_t linkpair, Int_t *mcmList, Int_t listSize )
{
  // Converts the chipmask to a list of MCMs 
  
  Short_t nmcm = 0;
  Short_t i;
  for( i = 0 ; i < listSize ; i++ ) {
     if( (cmA & (1 << i)) != 0 ) {
        mcmList[nmcm] = ((linkpair*2) << 7) | i;
	++nmcm;
    }
    if( (cmB & (1 << i)) != 0 ) {
       mcmList[nmcm] = ((linkpair*2+1) << 7) | i;
       ++nmcm;
    }
  }

  mcmList[nmcm]=-1;
  return nmcm;
}




//_____________________________________________________________________________
//void AliTRDfeeParam::GetFilterParam( Float_t &r1, Float_t &r2, Float_t &c1
//                                   , Float_t &c2, Float_t &ped ) const
//{
  //
  // Return current filter parameter
  //

  //  r1            = fR1;
  //r2            = fR2;
  //c1            = fC1;
  //c2            = fC2;
  //ped           = fPedestal;
//};

//_____________________________________________________________________________
void AliTRDfeeParam::SetEBsglIndThr(Int_t val)
{
  //
  // Set Event Buffer Sngle Indicator Threshold (EBIS in TRAP conf).
  // Timebin is indicated if ADC value >= val.
  //

  if( val >= 0 && val <= 1023 ) { 
    fEBsglIndThr = val;
  } else {
    AliError(Form("EBsglIndThr value %d is out of range, keep previously set value (%d).",
		  val, fEBsglIndThr));
  }

}

//_____________________________________________________________________________
void AliTRDfeeParam::SetEBsumIndThr(Int_t val)
{
  //
  // Set Event Buffer Sum Indicator Threshold (EBIT in TRAP conf).
  // Timebin is indicated if ADC sum value >= val.
  //

  if( val >= 0 && val <= 4095 ) { 
    fEBsumIndThr = val;
  } 
  else {
    AliError(Form("EBsumIndThr value %d is out of range, keep previously set value (%d).",
		  val, fEBsumIndThr));
  }

}

//_____________________________________________________________________________
void AliTRDfeeParam::SetEBindLUT(Int_t val)
{
  //
  // Set Event Buffer Indicator Look-Up Table (EBIL in TRAP conf).
  // 8 bits value forms lookup table for combination of three criterions.
  //

  if( val >= 0 && val <= 255 ) {
    fEBindLUT = val;
  } 
  else {
    AliError(Form("EBindLUT value %d is out of range, keep previously set value (%d).",
		  val, fEBindLUT));
  }

}

//_____________________________________________________________________________
void AliTRDfeeParam::SetEBignoreNeighbour(Int_t val)
{
  //
  // Set Event Buffer Indicator Neighbor Sensitivity. (EBIN in TRAP conf).
  // If 0, take account of neigbor's values.
  //

  if( val >= 0 && val <= 1 ) {
    fEBignoreNeighbour = val;
  } 
  else {
    AliError(Form("EBignoreNeighbour value %d is out of range, keep previously set value (%d).",
		  val, fEBignoreNeighbour));
  }
}

//_____________________________________________________________________________
void AliTRDfeeParam::SetRAWversion( Int_t rawver )
{
  //
  // Set raw data version (major number only)
  // Maximum available number is preset in fgkMaxRAWversion
  //

  if( rawver >= 0 && rawver <= fgkMaxRAWversion ) {
    fRAWversion = rawver;
  } 
  else {
    AliError(Form("Raw version is out of range: %d",rawver));
  }

}

//_____________________________________________________________________________
void AliTRDfeeParam::SetXenon()
{
  //
  // Sets the filter parameters for the Xenon gas mixture
  //

  fTFnExp = 1;
  fTFr1   = 1.1563;
  fTFr2   = 0.1299;
  fTFc1   = 0.0657;
  fTFc2   = 0.0000;

}

//_____________________________________________________________________________
void AliTRDfeeParam::SetArgon()
{
  //
  // Sets the filter parameters for the Argon gas mixture
  //

  fTFnExp = 2;
  fTFr1   = 6.0;
  fTFr2   = 0.62;
  fTFc1   = 0.0087;
  fTFc2   = 0.07;

}


