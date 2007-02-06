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
//
// Class AliMpMotifType
// --------------------
// Class that defines the motif properties.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpMotifType.h"
#include "AliMpMotifTypePadIterator.h"
#include "AliMpConnection.h"

#include "AliLog.h"
#include "AliMpFiles.h"
#include "TSystem.h"

#include <Riostream.h>

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
    fVerboseLevel(0),
#ifdef WITH_STL
    fConnections()
#endif
#ifdef WITH_ROOT
    fConnections(true)
#endif
{
  /// Standard constructor
      AliDebug(1,Form("this=%p id=%s",this,id.Data()));
}

//______________________________________________________________________________
AliMpMotifType::AliMpMotifType() 
  : TObject(),
    fID(""),
    fNofPadsX(0),   
    fNofPadsY(0),
    fVerboseLevel(0),
    fConnections()
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
  fVerboseLevel(0),
  fConnections()
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
  mt.fVerboseLevel = fVerboseLevel;
  mt.fConnections = fConnections;
}

//______________________________________________________________________________
AliMpMotifType::~AliMpMotifType() 
{
/// Destructor

#ifdef WITH_STL
 for(ConnectionMapCIterator i = fConnections.begin();
  i!=fConnections.end();++i)
   delete i->second;

  fConnections.erase(fConnections.begin(),fConnections.end());
#endif  
  
  AliDebug(1,Form("this=%p",this));
//  StdoutToAliDebug(1,this->Print(););
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
void AliMpMotifType::AddConnection(const AliMpIntPair &localIndices, 
                               AliMpConnection* connection)
{
  /// Add the connection to the map
  
#ifdef WITH_STL
  fConnections[localIndices]=connection;
#endif

#ifdef WITH_ROOT
  fConnections.Add(localIndices, connection);
#endif   

  connection->SetOwner(this);
}  

//______________________________________________________________________________
AliMpConnection *AliMpMotifType::FindConnectionByPadNum(Int_t padNum) const
{
  /// Retrieve the AliMpConnection pointer from its pad num
  
#ifdef WITH_STL
 for(ConnectionMapCIterator i = fConnections.begin();
  i!=fConnections.end();++i)
   if (i->second->GetPadNum()==padNum) return i->second;
 return 0;
#endif

#ifdef WITH_ROOT
  TExMapIter i = fConnections.GetIterator();
  Long_t key, value;
  while ( i.Next(key, value) ) {
    AliMpConnection* connection = (AliMpConnection*)value;
    if (connection->GetPadNum()==padNum) return connection;
  }  
 return 0;
#endif
}

//______________________________________________________________________________
AliMpConnection *AliMpMotifType::FindConnectionByLocalIndices(
                                       const AliMpIntPair& localIndices) const
{
  /// Retrieve the AliMpConnection pointer from its position (in pad unit)
  
  if (!localIndices.IsValid()) return 0;

#ifdef WITH_STL
  ConnectionMapCIterator i = fConnections.find(localIndices);
 if (i != fConnections.end())
   return i->second;
 else return 0;
#endif

#ifdef WITH_ROOT
  return (AliMpConnection*)fConnections.GetValue(localIndices);
#endif
}

//______________________________________________________________________________
AliMpConnection *AliMpMotifType::FindConnectionByGassiNum(Int_t gassiNum) const
{
  /// Return the connection for the given gassiplex number
  
#ifdef WITH_STL
 for(ConnectionMapCIterator i = fConnections.begin();
  i!=fConnections.end();++i)
   if (i->second->GetGassiNum()==gassiNum) return i->second;
 return 0;
#endif

#ifdef WITH_ROOT
  TExMapIter i = fConnections.GetIterator();
  Long_t key, value;
  while ( i.Next(key, value) ) {
    AliMpConnection* connection = (AliMpConnection*)value;
    if (connection->GetGassiNum()==gassiNum) return connection;
  }  
 return 0;
#endif
}

//______________________________________________________________________________
AliMpConnection *AliMpMotifType::FindConnectionByKaptonNum(Int_t kaptonNum) const
{
  /// Give the connection related to the given kapton number
  
#ifdef WITH_STL
 for(ConnectionMapCIterator i = fConnections.begin();
  i!=fConnections.end();++i)
   if (i->second->GetKaptonNum()==kaptonNum) return i->second;
 return 0;
#endif

#ifdef WITH_ROOT
  TExMapIter i = fConnections.GetIterator();
  Long_t key, value;
  while ( i.Next(key, value) ) {
    AliMpConnection* connection = (AliMpConnection*)value;
    if (connection->GetKaptonNum()==kaptonNum) return connection;
  }  
 return 0;
#endif
}
//______________________________________________________________________________
AliMpConnection *AliMpMotifType::FindConnectionByBergNum(Int_t bergNum) const
{
  /// Retrieve the connection from a Berg connector number
  
#ifdef WITH_STL
 for(ConnectionMapCIterator i = fConnections.begin();
  i!=fConnections.end();++i)
   if (i->second->GetBergNum()==bergNum) return i->second;
 return 0;
#endif

#ifdef WITH_ROOT
  TExMapIter i = fConnections.GetIterator();
  Long_t key, value;
  while ( i.Next(key, value) ) {
    AliMpConnection* connection = (AliMpConnection*)value;
    if (connection->GetBergNum()==bergNum) return connection;
  }  
  return 0;
#endif
}


//______________________________________________________________________________
AliMpIntPair AliMpMotifType::FindLocalIndicesByConnection(
                                 const AliMpConnection* connection) const
{
  /// Retrieve the pad position from the connection pointer.
  /// Not to be used widely, since it use a search in the
  /// connection list...

#ifdef WITH_STL
 for(ConnectionMapCIterator i = fConnections.begin();
  i!=fConnections.end();++i)
   if (i->second==connection) return i->first;
#endif

#ifdef WITH_ROOT
  TExMapIter i = fConnections.GetIterator();
  Long_t key, value;
  while ( i.Next(key, value) ) {
    AliMpConnection* aConnection = (AliMpConnection*)value;
    if (aConnection == connection) return AliMpExMap::GetPair(key);
  }  
#endif

  return AliMpIntPair::Invalid();
}

//______________________________________________________________________________
AliMpIntPair AliMpMotifType::FindLocalIndicesByPadNum(Int_t padNum) const
{
  /// Retrieve the AliMpConnection pointer from its pad num
  
#ifdef WITH_STL
 for(ConnectionMapCIterator i = fConnections.begin();
  i!=fConnections.end();++i)
   if (i->second->GetPadNum()==padNum) return i->first;
#endif
   
#ifdef WITH_ROOT
  TExMapIter i = fConnections.GetIterator();
  Long_t key, value;
  while ( i.Next(key, value) ) {
    AliMpConnection* connection = (AliMpConnection*)value;
    if (connection->GetPadNum() == padNum) return AliMpExMap::GetPair(key);
  }  
#endif
 return AliMpIntPair::Invalid();
}

//______________________________________________________________________________
AliMpIntPair AliMpMotifType::FindLocalIndicesByGassiNum(Int_t gassiNum) const
{
  /// Return the connection for the given gassiplex number
  
#ifdef WITH_STL
 for(ConnectionMapCIterator i = fConnections.begin();
  i!=fConnections.end();++i)
   if (i->second->GetGassiNum()==gassiNum) return i->first;
#endif
   
#ifdef WITH_ROOT
  TExMapIter i = fConnections.GetIterator();
  Long_t key, value;
  while ( i.Next(key, value) ) {
    AliMpConnection* connection = (AliMpConnection*)value;
    if (connection->GetGassiNum()==gassiNum) return AliMpExMap::GetPair(key);
  }  
#endif
   
 return AliMpIntPair::Invalid();
}

//______________________________________________________________________________
AliMpIntPair AliMpMotifType::FindLocalIndicesByKaptonNum(Int_t kaptonNum) const
{
  /// Give the connection related to the given kapton number
  
#ifdef WITH_STL
 for(ConnectionMapCIterator i = fConnections.begin();
  i!=fConnections.end();++i)
   if (i->second->GetKaptonNum()==kaptonNum) return i->first;
#endif
   
#ifdef WITH_ROOT
  TExMapIter i = fConnections.GetIterator();
  Long_t key, value;
  while ( i.Next(key, value) ) {
    AliMpConnection* connection = (AliMpConnection*)value;
    if (connection->GetKaptonNum()==kaptonNum) return AliMpExMap::GetPair(key);
  }  
#endif
   
 return AliMpIntPair::Invalid();
}

//______________________________________________________________________________
AliMpIntPair AliMpMotifType::FindLocalIndicesByBergNum(Int_t bergNum) const
{
  /// Retrieve the connection from a Berg connector number
  
#ifdef WITH_STL
 for(ConnectionMapCIterator i = fConnections.begin();
  i!=fConnections.end();++i)
   if (i->second->GetBergNum()==bergNum) return i->first;
#endif
   
#ifdef WITH_ROOT
  TExMapIter i = fConnections.GetIterator();
  Long_t key, value;
  while ( i.Next(key, value) ) {
    AliMpConnection* connection = (AliMpConnection*)value;
    if (connection->GetBergNum()==bergNum) return AliMpExMap::GetPair(key);
  }  
#endif
   
 return AliMpIntPair::Invalid();
}

//______________________________________________________________________________
Int_t  AliMpMotifType::GetNofPads() const   
{
/// Return the number of pads

#ifdef WITH_STL
  return fConnections.size();
#endif
   
#ifdef WITH_ROOT
  return fConnections.GetSize();
#endif
}

//______________________________________________________________________________
Bool_t AliMpMotifType::HasPad(const AliMpIntPair& localIndices) const
{
  /// Return true if the pad indexed by \a localIndices has a connection
  
  if (!localIndices.IsValid()) return false;

#ifdef WITH_STL
  return fConnections.find(localIndices)!=fConnections.end();
#endif

#ifdef WITH_ROOT
  TObject* value = fConnections.GetValue(localIndices);
  return value!=0;
#endif
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
      AliMpConnection *connexion = FindConnectionByLocalIndices(AliMpIntPair(i,j));
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
          case 'G':str=Form("%d",connexion->GetGassiNum());
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
    << "#abort, so it's not a valid value here)." << endl
    << "#" << endl;
  
  for ( Int_t ix = 0; ix < GetNofPadsX(); ++ix ) 
  {
    for ( Int_t iy = 0; iy < GetNofPadsY(); ++iy ) 
    {
      AliMpConnection* con = FindConnectionByLocalIndices(AliMpIntPair(ix,iy));
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



