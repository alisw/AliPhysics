/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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
///////////////////////////////////////////////////////////////////////////
//  Base Plane Efficiency class for ITS            
//  Specific subdetector implementation is done in  
//  AliITSPlaneEffSPD                               
//  AliITSPlaneEffSDD                               
//  AliITSPlaneEffSSD                               
//
//  Author: G.E. Bruno 
//          giuseppe.bruno@ba.infn.it
//
///////////////////////////////////////////////////////////////////////////

/* $Id$ */

#include <TMath.h>
#include "AliITSPlaneEff.h"
#include "AliLog.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"

ClassImp(AliITSPlaneEff)
//______________________________________________________________________
AliITSPlaneEff::AliITSPlaneEff(): AliPlaneEff(),
fRunNumber(0), 
fCDBUri(""),
fInitCDBCalled(kFALSE),
fHis(kFALSE)
{
    // Default constructor
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    a default constructed AliITSPlaneEff class
 InitCDB();
}
//______________________________________________________________________
AliITSPlaneEff::AliITSPlaneEff(const AliITSPlaneEff &s) : AliPlaneEff(s),
fRunNumber(s.fRunNumber),
fCDBUri(s.fCDBUri),
fInitCDBCalled(s.fInitCDBCalled),
fHis(s.fHis)
{
    //     Copy Constructor
    // Inputs:
    //    const AliITSPlaneEff &s  simulation class to copy from
    // Outputs:
    //    none.
    // Return:
    //    a standard constructed AliITSPlaneEff class with values the same
    //    as that of s.

}
//_________________________________________________________________________
//void AliITSPlaneEff::operator+=(const AliITSPlaneEff &add){
    //    Add to me operator
    // Inputs:
    //    const AliITSPlaneEff &add  simulation class to be added 
    // Outputs:
    //    none.
    // Return:
    //    none

//    return;
//}
//_________________________________________________________________________
AliITSPlaneEff&  AliITSPlaneEff::operator=(const AliITSPlaneEff &source){
    //    Assignment operator
    // Inputs:
    //    const AliITSPlaneEff &source  simulation class to copy from
    // Outputs:
    //    none.
    // Return:
    //    a standard constructed AliITSPlaneEff class with values the same
    //    as that of s.
    if(this != &source){
       source.Copy(*this);
    }
    return *this;
}
//_________________________________________________________________________
/*
AliPlaneEff&  AliITSPlaneEff::operator=(const
                                           AliPlaneEff &s){
    //    Assignment operator
    // Inputs:
    //    AliPlaneEff &s The original class for which
    //                          this class is a copy of
    // Outputs:
    //    none.
    // Return: 

    if(&s == this) return *this;
    AliWarning("AliITSPlaneEff Not allowed to make a = Using default creator instead");
    return *this;
}
*/
//_________________________________________________________________________
void AliITSPlaneEff::Copy(TObject &obj) const {
  // copy this to obj
  ((AliITSPlaneEff& ) obj).fRunNumber		= fRunNumber;
  ((AliITSPlaneEff& ) obj).fCDBUri		= fCDBUri;
  ((AliITSPlaneEff& ) obj).fInitCDBCalled	= fInitCDBCalled;
  ((AliITSPlaneEff& ) obj).fHis			= fHis;
}
//_________________________________________________________________________
Double_t AliITSPlaneEff::PlaneEff(Int_t nf,Int_t nt) const {
   // Compute the efficiency for a basic block, 
    // Inputs:
    //        number of associated cluslters (nf) 
    //        number of used tracks (nt)
    // Outputs:
    //    none.
    // Return:
    //        the efficiency 
if(nf<0 || nt<=0 || nt<nf) {
   AliInfo(Form("AliITSPlaneEff::PlaneEff: nfound= %i, ntried= %i",nf,nt)); 
   return -1.;}
 Double_t eff=nf;
 return eff/=nt;
}
//_________________________________________________________________________
Double_t AliITSPlaneEff::ErrPlaneEff(Int_t nf,Int_t nt) const{
    // Compute the statistical error on efficiency for a basic block,
    // using binomial statistics 
    // Inputs:
    //        number of associated cluslters (nf)
    //        number of used tracks (nt)
    // Outputs:
    //    none.
    // Return:
    //        the error on the efficiency 
if(nf<0 || nt<=0 || nt<nf) {
   AliInfo(Form("AliITSPlaneEff::ErrPlaneEff: nfound= %i, ntried= %i",nf,nt)); 
   return -1.;}
 Double_t err=TMath::Sqrt((Double_t)nf*(1.-(Double_t)nf/(Double_t)nt));
 return err/=(Double_t)nt;
}
//______________________________________________________________________
Int_t AliITSPlaneEff::GetNTracksForGivenEff(Double_t eff, Double_t RelErr) const {
    // Estimate of the number of tracks needed for measuring efficiency 
    // with a given relative error using binomial statistics
    // Inputs:
    //        exspected efficiency eff (e.g. from previous measurements)
    //        wished relative error RelErr
    // Outputs: 
    //    none.
    // Return: number of tracks given as the nearest integer 
if(eff<=0 || eff>1 || RelErr<=0 ) {
   AliInfo(Form("AliITSPlaneEff::GetNTracksForGivenEff: eff= %f, RelErr= %f",eff,RelErr));
   return -1;}
return TMath::Nint((1-eff)/(eff*RelErr*RelErr));
}
//________________________________________________________________________
void AliITSPlaneEff::InitCDB() 
{
// activate a default CDB storage
// First check if we have any CDB storage set, because it is used
// to retrieve the calibration and alignment constants

  if (fInitCDBCalled) return;
  fInitCDBCalled = kTRUE;

  AliCDBManager* man = AliCDBManager::Instance();
  if (man->IsDefaultStorageSet())
  {
    AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    AliWarning("Default CDB storage has been already set !");
    AliWarning(Form("Ignoring the default storage declared in AliITSPlaneEff: %s",fCDBUri.Data()));
    AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    fCDBUri = man->GetDefaultStorage()->GetURI();
  }
  else {
    if (fCDBUri.Length() > 0)
    {
        AliDebug(2,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        AliDebug(2, Form("Default CDB storage is set to: %s", fCDBUri.Data()));
        AliDebug(2, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    } else {
        fCDBUri="local://$ALICE_ROOT/OCDB";
        AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        AliWarning("Default CDB storage not yet set !!!!");
        AliWarning(Form("Setting it now to: %s", fCDBUri.Data()));
        AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");

    }
    man->SetDefaultStorage(fCDBUri);
  }
}
//_____________________________________________________________________________
void AliITSPlaneEff::SetDefaultStorage(const char* uri) {
// Store the desired default CDB storage location
// Activate it later within the Run() method

  fCDBUri = uri;

}

