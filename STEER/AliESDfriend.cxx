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

//-------------------------------------------------------------------------
//               Implementation of the AliESDfriend class
//  This class contains some additional to the ESD information like
//  the clusters associated to tracks.
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-------------------------------------------------------------------------

#include "AliESDfriend.h"
#include "AliESDVZEROfriend.h"
#include "AliESDTZEROfriend.h"

ClassImp(AliESDfriend)

AliESDfriend::AliESDfriend(): TObject(), fTracks("AliESDfriendTrack",15000),
  fESDVZEROfriend(NULL),
  fESDTZEROfriend(NULL)

{
 //
 // Default constructor
 //
}

AliESDfriend::AliESDfriend(const AliESDfriend &f) :
  TObject(f),
  fTracks(f.fTracks),
  fESDVZEROfriend(f.fESDVZEROfriend ? new AliESDVZEROfriend(*f.fESDVZEROfriend) : NULL),
  fESDTZEROfriend(f.fESDTZEROfriend ? new AliESDTZEROfriend(*f.fESDTZEROfriend) : NULL)
{
 //
 // Copy constructor
 //
}

AliESDfriend& AliESDfriend::operator=(const AliESDfriend& esd)
{
    
    // Assignment operator
    if(&esd == this) return *this;
    TObject::operator=(esd);
    fTracks = esd.fTracks;

    delete fESDVZEROfriend; fESDVZEROfriend = NULL;
    if (!esd.fESDVZEROfriend) fESDVZEROfriend = new AliESDVZEROfriend(*esd.fESDVZEROfriend);

    delete fESDTZEROfriend; fESDTZEROfriend = NULL;
    if (!esd.fESDTZEROfriend) fESDTZEROfriend = new AliESDTZEROfriend(*esd.fESDTZEROfriend);
 
 
 
    return *this;
}



AliESDfriend::~AliESDfriend() {
  //
  // Destructor
  //
  fTracks.Delete();
  delete fESDVZEROfriend;
  delete fESDTZEROfriend;
}



void AliESDfriend::SetVZEROfriend(AliESDVZEROfriend * obj)
{
  //
  // Set the VZERO friend data object
  // (complete raw data)
  if (!fESDVZEROfriend) fESDVZEROfriend = new AliESDVZEROfriend();
  if (obj) *fESDVZEROfriend = *obj;
}
void AliESDfriend::SetTZEROfriend(AliESDTZEROfriend * obj)
{
  //
  // Set the TZERO friend data object
  // (complete raw data)
  if (!fESDTZEROfriend) fESDTZEROfriend = new AliESDTZEROfriend();
  if (obj) *fESDTZEROfriend = *obj;
}
