// $Id$
// Category: motif
//
// Class AliMpMotifType
// --------------------
// Class that defines the motif properties.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <Riostream.h>

#include "AliMpMotifType.h"
#include "AliMpMotifTypePadIterator.h"
#include "AliMpConnection.h"

ClassImp(AliMpMotifType)

const Int_t AliMpMotifType::fgkPadNumForA = 65;



//______________________________________________________________________________
AliMpMotifType::AliMpMotifType(const TString &id) 
  : TObject(),
    fID(id),
    fNofPadsX(0),   
    fNofPadsY(0),
    fVerboseLevel(0),
    fConnections()
{
  // Constructor
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
  // Default constructor (dummy)
}

//______________________________________________________________________________
AliMpMotifType::~AliMpMotifType() {
// Destructor

 for(ConnectionMapCIterator i = fConnections.begin();
  i!=fConnections.end();++i)
   delete i->second;

  fConnections.erase(fConnections.begin(),fConnections.end());
}

//______________________________________________________________________________
AliMpVPadIterator* AliMpMotifType::CreateIterator() const
{
  return new AliMpMotifTypePadIterator(this);
}

//______________________________________________________________________________
void AliMpMotifType::SetNofPads(Int_t nofPadsX, Int_t nofPadsY)
{
  // Change the number of pads in this motif

  fNofPadsX = nofPadsX;
  fNofPadsY = nofPadsY;
}


//______________________________________________________________________________
Int_t AliMpMotifType::PadNum(const TString &padName) const
{
  // Transform a pad name into the equivalent pad number
  if ( (padName[0]>='A') && (padName[0]<='Z') )
    return fgkPadNumForA+padName[0]-'A';
  else
    return atoi(padName.Data());
}

//______________________________________________________________________________
TString AliMpMotifType::PadName(Int_t padNum) const
{
  // Transform a pad number into its equivalent pad name
  if (padNum<fgkPadNumForA)
    return Form("%d",padNum);
  else
    return char('A'+padNum-fgkPadNumForA);
}

//______________________________________________________________________________
void AliMpMotifType::AddConnection(const AliMpIntPair &localIndices, 
                               AliMpConnection* connection)
{
  // Add the connection to the map
  
  fConnections[localIndices]=connection;
  connection->SetOwner(this);
}  
//______________________________________________________________________________
AliMpConnection *AliMpMotifType::FindConnectionByPadNum(Int_t padNum) const
{
  // Retrive the AliMpConnection pointer from its pad num
 for(ConnectionMapCIterator i = fConnections.begin();
  i!=fConnections.end();++i)
   if (i->second->GetPadNum()==padNum) return i->second;
 return 0;
}

//______________________________________________________________________________
AliMpConnection *AliMpMotifType::FindConnectionByLocalIndices(
                                       const AliMpIntPair& localIndices) const
{
  if (!localIndices.IsValid()) return 0;

  // Retrive the AliMpConnection pointer from its position (in pad unit)
  ConnectionMapCIterator i = fConnections.find(localIndices);
 if (i != fConnections.end())
   return i->second;
 else return 0;
}

//______________________________________________________________________________
AliMpConnection *AliMpMotifType::FindConnectionByGassiNum(Int_t gassiNum) const
{
  // return the connection for the given gassiplex number
 for(ConnectionMapCIterator i = fConnections.begin();
  i!=fConnections.end();++i)
   if (i->second->GetGassiNum()==gassiNum) return i->second;
 return 0;
}
//______________________________________________________________________________
AliMpConnection *AliMpMotifType::FindConnectionByKaptonNum(Int_t kaptonNum) const
{
  // Gives the connection related to the given kapton number
 for(ConnectionMapCIterator i = fConnections.begin();
  i!=fConnections.end();++i)
   if (i->second->GetKaptonNum()==kaptonNum) return i->second;
 return 0;
}
//______________________________________________________________________________
AliMpConnection *AliMpMotifType::FindConnectionByBergNum(Int_t bergNum) const
{
  // Retrieve the connection from a Berg connector number
 for(ConnectionMapCIterator i = fConnections.begin();
  i!=fConnections.end();++i)
   if (i->second->GetBergNum()==bergNum) return i->second;
 return 0;
}


//______________________________________________________________________________
AliMpIntPair AliMpMotifType::FindLocalIndicesByConnection(
                                         const AliMpConnection* connection)
{
  // Retrieve the pad position from the connection pointer.
  // Not to be used widely, since it use a search in the
  // connection list...

 for(ConnectionMapCIterator i = fConnections.begin();
  i!=fConnections.end();++i)
   if (i->second==connection) return i->first;

 return AliMpIntPair::Invalid();
}

//______________________________________________________________________________
AliMpIntPair AliMpMotifType::FindLocalIndicesByPadNum(Int_t padNum) const
{
  // Retrive the AliMpConnection pointer from its pad num
 for(ConnectionMapCIterator i = fConnections.begin();
  i!=fConnections.end();++i)
   if (i->second->GetPadNum()==padNum) return i->first;
   
 return AliMpIntPair::Invalid();
}

//______________________________________________________________________________
AliMpIntPair AliMpMotifType::FindLocalIndicesByGassiNum(Int_t gassiNum) const
{
  // return the connection for the given gassiplex number
 for(ConnectionMapCIterator i = fConnections.begin();
  i!=fConnections.end();++i)
   if (i->second->GetGassiNum()==gassiNum) return i->first;
   
 return AliMpIntPair::Invalid();
}

//______________________________________________________________________________
AliMpIntPair AliMpMotifType::FindLocalIndicesByKaptonNum(Int_t kaptonNum) const
{
  // Gives the connection related to the given kapton number
 for(ConnectionMapCIterator i = fConnections.begin();
  i!=fConnections.end();++i)
   if (i->second->GetKaptonNum()==kaptonNum) return i->first;
   
 return AliMpIntPair::Invalid();
}

//______________________________________________________________________________
AliMpIntPair AliMpMotifType::FindLocalIndicesByBergNum(Int_t bergNum) const
{
  // Retrieve the connection from a Berg connector number
 for(ConnectionMapCIterator i = fConnections.begin();
  i!=fConnections.end();++i)
   if (i->second->GetBergNum()==bergNum) return i->first;
   
 return AliMpIntPair::Invalid();
}

//______________________________________________________________________________
Bool_t AliMpMotifType::HasPad(const AliMpIntPair& localIndices) const
{
  if (!localIndices.IsValid()) return false;

  // return true if the pad indexed by <localIndices> has a connection
  return fConnections.find(localIndices)!=fConnections.end();

}

//______________________________________________________________________________
void AliMpMotifType::Print(Option_t *option) const
{
  // Print the map of the motif. In each cel, the value
  // printed depends of option, as the following:
  // option="N" the "name" of the pad is written
  // option="K" the Kapton connect. number attached to the pad is written
  // option="B" the Berg connect. number attached to the pad is written
  // option="G" the Gassiplex channel number attached to the pad is written
  // otherwise the number of the pad is written

  // NOTE : this method is really not optimized, in case 'N' or '',
  // but the Print() this should not be very important in a Print() method

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
