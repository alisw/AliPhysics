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

// $Id$
// $MpId: AliMpMotifType.cxx,v 1.10 2006/05/24 13:58:41 ivana Exp $
// Category: motif

//-----------------------------------------------------------------------------
// Class AliMpMotifType
// --------------------
// Class that defines the motif properties.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay
//-----------------------------------------------------------------------------

#include "AliMpMotifType.h"
#include "AliMpExMapIterator.h"
#include "AliMpMotifTypePadIterator.h"
#include "AliMpConnection.h"
#include "AliMpConstants.h"
#include "AliMpFiles.h"
#include "AliMpEncodePair.h"

#include "AliLog.h"

#include <TSystem.h>
#include <Riostream.h>

#include <cstdlib>

/// \cond CLASSIMP
ClassImp(AliMpMotifType)
/// \endcond

const Int_t AliMpMotifType::fgkPadNumForA = 65;

//______________________________________________________________________________
AliMpMotifType::AliMpMotifType(const TString &id) 
: TObject(),
fID(id),
fNofPadsX(0),   
fNofPadsY(0),
fNofPads(0),
fMaxNofPads(AliMpConstants::ManuNofChannels()),
fConnectionsByLocalIndices(fMaxNofPads*fMaxNofPads),
fConnectionsByManuChannel(fMaxNofPads)
{
  /// Standard constructor                                                   \n
  /// Please note that id should be of the form %s for station 1,2,
  //  %s-%e-%e for station345 and %sx%e for stationTrigger

  fConnectionsByLocalIndices.SetOwner(kTRUE);
  fConnectionsByManuChannel.SetOwner(kFALSE);
  AliDebug(1,Form("this=%p id=%s",this,id.Data()));
}

//______________________________________________________________________________
AliMpMotifType::AliMpMotifType(TRootIOCtor*) 
: TObject(),
fID(""),
fNofPadsX(0),   
fNofPadsY(0),
fNofPads(0),
fMaxNofPads(0),
fConnectionsByLocalIndices(),
fConnectionsByManuChannel()
{
  /// Default constructor
  AliDebug(1,Form("this=%p",this));
}

//______________________________________________________________________________
AliMpMotifType::AliMpMotifType(const AliMpMotifType& rhs)
: TObject(),
fID(""),
fNofPadsX(0),   
fNofPadsY(0),
fNofPads(0),
fMaxNofPads(0),
fConnectionsByLocalIndices(),
fConnectionsByManuChannel()
{
  /// Copy constructor
  AliDebug(1,Form("this=%p (copy ctor)",this));
  rhs.Copy(*this);
}

//______________________________________________________________________________
AliMpMotifType&
AliMpMotifType::operator=(const AliMpMotifType& rhs)
{
  /// Assignment operator
  
  TObject::operator=(rhs);
  rhs.Copy(*this);
  return *this;  
}

//______________________________________________________________________________
TObject*
AliMpMotifType::Clone(const char* /*newname*/) const 
{
  /// Returns a full copy of this object
  return new AliMpMotifType(*this);
}

//______________________________________________________________________________
void
AliMpMotifType::Copy(TObject& object) const
{
  /// Copy object

  TObject::Copy(object);
  AliMpMotifType& mt = static_cast<AliMpMotifType&>(object);
  mt.fID = fID;
  mt.fNofPadsX = fNofPadsX;
  mt.fNofPadsY = fNofPadsY;
  mt.fNofPads = fNofPads;
  mt.fMaxNofPads = fMaxNofPads;
  mt.fConnectionsByLocalIndices = fConnectionsByLocalIndices;
  mt.fConnectionsByManuChannel = fConnectionsByManuChannel;  
}

//______________________________________________________________________________
AliMpMotifType::~AliMpMotifType() 
{
  /// Destructor

  AliDebug(1,Form("this=%p",this));
}

//______________________________________________________________________________
AliMpVPadIterator* AliMpMotifType::CreateIterator() const
{
/// Create new motif type iterator

  return new AliMpMotifTypePadIterator(this);
}

//______________________________________________________________________________
void AliMpMotifType::SetNofPads(Int_t nofPadsX, Int_t nofPadsY)
{
  /// Change the number of pads in this motif

  fNofPadsX = nofPadsX;
  fNofPadsY = nofPadsY;
}


//______________________________________________________________________________
Int_t AliMpMotifType::PadNum(const TString &padName) const
{
  /// Transform a pad name into the equivalent pad number

  if ( (padName[0]>='A') && (padName[0]<='Z') )
    return fgkPadNumForA+padName[0]-'A';
  else
    return atoi(padName.Data());
}

//______________________________________________________________________________
TString AliMpMotifType::PadName(Int_t padNum) const
{
  /// Transform a pad number into its equivalent pad name

  if (padNum<fgkPadNumForA)
    return Form("%d",padNum);
  else
    return char('A'+padNum-fgkPadNumForA);
}

//______________________________________________________________________________
Bool_t 
AliMpMotifType::AddConnection(AliMpConnection* connection)
{
  /// Add the connection to the map
  
  if (!connection) return kFALSE;
  
  Int_t ix = connection->GetLocalIx();
  Int_t iy = connection->GetLocalIy();
  
  Int_t manuChannel = connection->GetManuChannel();
  
  if ( ix >=0 && ix < fMaxNofPads &&
      iy >=0 && iy < fMaxNofPads && 
      manuChannel >= 0 && manuChannel < AliMpConstants::ManuNofChannels())
  {
  
    Int_t index = ix + iy*AliMpConstants::ManuNofChannels();
    
    AliMpConnection* c = FindConnectionByLocalIndices(
                             connection->GetLocalIndices());
    
    if (c)
    {
      AliError(Form("Connection already exists for ix=%d iy=%d",ix,iy));
      return kFALSE;
    }
    
    ++fNofPads;

    fConnectionsByLocalIndices[index] = connection;
    fConnectionsByManuChannel[manuChannel] = connection;
    
    connection->SetOwner(this);
    
    return kTRUE;
  
  }
  return kFALSE;
}  

//______________________________________________________________________________
AliMpConnection*
AliMpMotifType::FindConnectionByPadNum(Int_t padNum) const
{
  /// Retrieve the AliMpConnection pointer from its pad num
  /// This method is quite inefficient as we're looping over all connections
  
  TIter next(&fConnectionsByManuChannel);
  AliMpConnection* connection;
  
  while ( ( connection = static_cast<AliMpConnection*>(next()) ) )
  {
    if (connection->GetPadNum()==padNum) return connection;
  }    
  return 0x0;
}

//______________________________________________________________________________
AliMpConnection*
AliMpMotifType::FindConnectionByLocalIndices(MpPair_t localIndices) const
{
  /// Retrieve the AliMpConnection pointer from its position (in pad unit)

  return FindConnectionByLocalIndices(AliMp::PairFirst(localIndices),
                                      AliMp::PairSecond(localIndices));
}

//______________________________________________________________________________
AliMpConnection*
AliMpMotifType::FindConnectionByLocalIndices(Int_t ix, Int_t iy) const
{
  /// Retrieve the AliMpConnection pointer from its position (in pad unit)

  if ( ix < fNofPadsX && iy < fNofPadsY && ix >= 0 && iy >= 0 )
  {  
    Int_t index = ix + iy*fMaxNofPads;

    return static_cast<AliMpConnection*>(fConnectionsByLocalIndices.UncheckedAt(index));
  }
  else
  {
    return 0x0;
  }
}

//______________________________________________________________________________
AliMpConnection*
AliMpMotifType::FindConnectionByGassiNum(Int_t gassiNum) const
{
  /// Return the connection for the given gassiplex number
  
  if ( gassiNum >=0 && gassiNum < fMaxNofPads ) 
  {
    return static_cast<AliMpConnection*>(fConnectionsByManuChannel.UncheckedAt(gassiNum));
  }
  
  return 0x0;
}

//______________________________________________________________________________
AliMpConnection*
AliMpMotifType::FindConnectionByKaptonNum(Int_t kaptonNum) const
{
  /// Give the connection related to the given kapton number
  /// Inefficient method as we loop over connections to find the right one
  
  TIter next(&fConnectionsByManuChannel);
  AliMpConnection* connection;
  
  while ( ( connection = static_cast<AliMpConnection*>(next()) ) )
  {
    if ( connection && connection->GetKaptonNum()==kaptonNum) return connection;
  }
  return 0x0;
}

//______________________________________________________________________________
AliMpConnection*
AliMpMotifType::FindConnectionByBergNum(Int_t bergNum) const
{
  /// Retrieve the connection from a Berg connector number
  /// Inefficient method as we loop over connections to find the right one
  
  TIter next(&fConnectionsByManuChannel);
  AliMpConnection* connection;
  
  while ( ( connection = static_cast<AliMpConnection*>(next()) ) )
  {
    if ( connection && connection->GetBergNum()==bergNum) return connection;
  }
  return 0x0;
}


//______________________________________________________________________________
MpPair_t AliMpMotifType::FindLocalIndicesByConnection(const AliMpConnection* connection) const
{
  /// Reurn the pad position from the connection pointer.

  return connection->GetLocalIndices();
}

//______________________________________________________________________________
MpPair_t AliMpMotifType::FindLocalIndicesByPadNum(Int_t padNum) const
{
  /// Retrieve the AliMpConnection pointer from its pad num
  
  AliMpConnection* connection = FindConnectionByPadNum(padNum);
  
  if ( ! connection) return -1;
  
  return connection->GetLocalIndices();
}

//______________________________________________________________________________
MpPair_t AliMpMotifType::FindLocalIndicesByGassiNum(Int_t gassiNum) const
{
  /// Return the connection for the given gassiplex number
  
  AliMpConnection* connection = FindConnectionByGassiNum(gassiNum);
  
  if ( ! connection) return -1;

  return connection->GetLocalIndices();
}

//______________________________________________________________________________
MpPair_t AliMpMotifType::FindLocalIndicesByKaptonNum(Int_t kaptonNum) const
{
  /// Give the connection related to the given kapton number

  AliMpConnection* connection = FindConnectionByKaptonNum(kaptonNum);
  
  if ( ! connection) return -1;

  return connection->GetLocalIndices();
}

//______________________________________________________________________________
MpPair_t AliMpMotifType::FindLocalIndicesByBergNum(Int_t bergNum) const
{
  /// Retrieve the connection from a Berg connector number
  
  AliMpConnection* connection = FindConnectionByBergNum(bergNum);
  
  if ( ! connection) return -1;

  return connection->GetLocalIndices();
}

//______________________________________________________________________________
Bool_t 
AliMpMotifType::HasPadByLocalIndices(MpPair_t localIndices) const
{
  /// Return true if the pad indexed by \a localIndices has a connection
    
  return ( FindConnectionByLocalIndices(localIndices) != 0x0 );
}

//______________________________________________________________________________
Bool_t 
AliMpMotifType::HasPadByLocalIndices(Int_t localIx, Int_t localIy) const
{
  /// Return true if the pad indexed by \a localIndices has a connection
    
  return ( FindConnectionByLocalIndices(localIx, localIy) != 0x0 );
}

//______________________________________________________________________________
Bool_t 
AliMpMotifType::HasPadByManuChannel(Int_t manuChannel) const
{
  /// Return true if the pad indexed by \a localIndices has a connection
  
//  if ( manuChannel >= fNofPads ) return kFALSE;
  
  return ( FindConnectionByGassiNum(manuChannel) != 0x0 );
}

//______________________________________________________________________________
void AliMpMotifType::Print(Option_t *option) const
{
  /// Print the map of the motif. In each cell, the value
  /// printed depends of option, as the following:
  /// - option="N" the "name" of the pad is written
  /// - option="K" the Kapton connect. number attached to the pad is written
  /// - option="B" the Berg connect. number attached to the pad is written
  /// - option="G" the Gassiplex channel number attached to the pad is written
  /// otherwise the number of the pad is written
  ///
  /// NOTE : this method is really not optimized, in case 'N' or '',
  /// but the Print() this should not be very important in a Print() method

  switch (option[0]){
  case 'N':cout<<"Name mapping";
    break;
  case 'K':cout<<"Kapton mapping";
    break;
  case 'B':cout<<"Berg mapping";
    break;
  case 'G':cout<<"Gassiplex number mapping";
    break;
  default:cout<<"Pad mapping";
  }
  cout<<" in the motif "<<fID<<endl;
  cout<<"-----------------------------------"<<endl;

  for (Int_t j=fNofPadsY-1;j>=0;j--){
    for (Int_t i=0;i<fNofPadsX;i++){
      AliMpConnection *connexion = FindConnectionByLocalIndices(i,j);
      TString str;
      if (connexion){
        AliDebug(1,Form("i,j=%2d,%2d connexion=%p",i,j,connexion));
        
        switch (option[0]){
          case 'N':str=PadName(connexion->GetPadNum());
            break;
          case 'K':str=Form("%d",connexion->GetKaptonNum());
            break;
          case 'B':str=Form("%d",connexion->GetBergNum());
            break;
          case 'G':str=Form("%d",connexion->GetManuChannel());
            break;
          default:str= Form("%d",connexion->GetPadNum());
        }
        cout<<setw(2)<<str;
      } else cout<<setw(2)<<"--";
      cout<<" ";
    }
    cout<<endl;
  }
}

//_____________________________________________________________________________
Bool_t
AliMpMotifType::Save() const
{
/// Save this motif type

  return Save(fID.Data());
}

//_____________________________________________________________________________
Bool_t
AliMpMotifType::Save(const char* motifName) const
{
  /// Generate the 2 files needed to describe the motif
  
  TString padPosFileName(AliMpFiles::PadPosFileName(motifName));
  
  TString motifTypeFileName(AliMpFiles::MotifFileName(motifName));

  // first a protection : do not allow overwriting existing files...
  Bool_t test = gSystem->AccessPathName(padPosFileName.Data());
  if (test==kFALSE) // AccessPathName has a strange return value convention...
  {
    AliError("Cannot overwrite existing padPos file");
    return kFALSE;
  }
  test = gSystem->AccessPathName(motifTypeFileName.Data());
  if (test==kFALSE)
  {
    AliError("Cannot overwrite existing motifType file");
    return kFALSE;    
  }
  
  ofstream padPosFile(padPosFileName.Data());
  ofstream motifFile(motifTypeFileName.Data());
  
  motifFile <<  "# Motif " << motifName << endl
    << "#" << endl
    << "#connecteur_berg kapton padname not_used" << endl
    << "#for slats there's no kapton connector, so it's always 1" 
    << " (zero make the reader" << endl
    << "#exit, so it's not a valid value here)." << endl
    << "#" << endl;
  
  for ( Int_t ix = 0; ix < GetNofPadsX(); ++ix ) 
  {
    for ( Int_t iy = 0; iy < GetNofPadsY(); ++iy ) 
    {
      AliMpConnection* con = FindConnectionByLocalIndices(ix,iy);
      if (con)
      {
        motifFile << con->GetBergNum() << "\t1\t" << con->GetPadNum() << "\t-" << endl;
        padPosFile << con->GetPadNum() << "\t" << ix << "\t" << iy << endl;
      }
    }
  }
  
  padPosFile.close();
  motifFile.close();
  
  return kTRUE;
}



