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
#include <Riostream.h>
#include <TMath.h>
#include <TRandom.h>

#include "AliITSresponseSPDdubna.h"
 
//___________________________________________
ClassImp(AliITSresponseSPDdubna)	

AliITSresponseSPDdubna::AliITSresponseSPDdubna() : AliITSresponse(){
   // Default constructor
   // Inputs:
   //    none.
   // Outputs:
   //    none.
   // Return:
   //    A default constructed AliITSresponseSPD class

   SetNoiseParam();  // fNoise, fBaseline
   SetThreshold();   // fThreshold
   SetCouplings();   // fCouplCol, fCouplRow
   SetFractionDeadPixels(); // fDeadPixels
   SetDataType();    // fDataType
}
//_________________________________________________________________________
Bool_t AliITSresponseSPDdubna::IsPixelDead(Int_t mod,Int_t ix,Int_t iz) const {
  // Returns kTRUE if pixel is dead
  // Inputs:
  //    Int_t mod      module number
  //    Int_t ix       x pixel number
  //    Int_t iz       z pixel number
  // Outputs:
  //    none.
  // Return:
  //    kFALSE if pixel is alive, or kTRUE if pixel is dead.
  Bool_t  dead = kFALSE;
  Int_t   seed;
  static TRandom ran; // don't use gRandom. This must not be a true randome
  // sequence. These sequence must be random one and then fully repetable.

  seed = mod*256*256+iz*256+ix;
  ran.SetSeed(seed);
  if(ran.Rndm(0)<fDeadPixels) dead = kTRUE;
  return dead;
}
//----------------------------------------------------------------------
void AliITSresponseSPDdubna::Print(ostream *os){
    // Standard output format for this class.
    // Inputs:
    //    ostream *os  Pointer to the output stream
    // Outputs:
    //    none:
    // Return:
    //    none.
#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif

    AliITSresponse::Print(os);
    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os << fNoise << " " << fBaseline << " " << fCouplCol << " ";
    *os << fCouplRow << " "<< fThreshold << " " << fDeadPixels << " ";
    *os << fDataType;
//    *os << " " << endl;
    os->flags(fmt); // reset back to old formating.
    return;
}
//----------------------------------------------------------------------
void AliITSresponseSPDdubna::Read(istream *is){
    // Standard input format for this class.
    // Inputs:
    //    ostream *os  Pointer to the output stream
    // Outputs:
    //    none:
    // Return:
    //    none.

    AliITSresponseSPDdubna::Read(is);
    *is >> fNoise >> fBaseline >> fCouplCol >> fCouplRow;
    *is >> fThreshold >> fDeadPixels >> fDataType;
    return;
}
//----------------------------------------------------------------------
ostream &operator<<(ostream &os,AliITSresponseSPDdubna &p){
    // Standard output streaming function.
    // Inputs:
    //    ostream *os  Pointer to the output stream
    // Outputs:
    //    none:
    // Return:
    //    none.

    p.Print(&os);
    return os;
}
//----------------------------------------------------------------------
istream &operator>>(istream &is,AliITSresponseSPDdubna &r){
    // Standard input streaming function.
    // Inputs:
    //    ostream *os  Pointer to the output stream
    // Outputs:
    //    none:
    // Return:
    //    none.

    r.Read(&is);
    return is;
}
//----------------------------------------------------------------------

