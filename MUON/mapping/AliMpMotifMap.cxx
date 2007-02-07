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
// $MpId: AliMpMotifMap.cxx,v 1.16 2006/05/24 13:58:41 ivana Exp $
// Category: motif
// -------------------
// Class AliMpMotifMap
// -------------------
// Class describing the motif map container, where motifs are
// mapped to their string IDs.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpMotifMap.h"
#include "AliMpVMotif.h"
#include "AliMpMotif.h"
#include "AliMpMotifSpecial.h"
#include "AliMpMotifType.h"
#include "AliMpMotifPosition.h"

#include "AliLog.h"

#include <Riostream.h>
#include <TVector2.h>
#include <TArrayI.h>

/// \cond CLASSIMP
ClassImp(AliMpMotifMap)
/// \endcond

//_____________________________________________________________________________
AliMpMotifMap::AliMpMotifMap(Bool_t /*standardConstructor*/) 
  : TObject()
#ifdef WITH_ROOT
    ,fMotifs(true),
     fMotifTypes(true),
     fMotifPositions(true),
     fMotifPositions2(true)
#endif 
{
/// Standard constructor
  
  //fMotifPositions2.SetOwner(false);
}

//_____________________________________________________________________________
AliMpMotifMap::AliMpMotifMap() 
  : TObject(),
    fMotifs(),
    fMotifTypes(),
    fMotifPositions(),
    fMotifPositions2()
{
/// Default constructor
}

//_____________________________________________________________________________
AliMpMotifMap::~AliMpMotifMap() 
{
/// Destructor  

  // Delete all registered motifs, motif types, motif positions
  
#ifdef WITH_STL
  for (MotifMapIterator im=fMotifs.begin(); im != fMotifs.end(); im++) {
    delete im->second;
  }  

  for (MotifTypeMapIterator it=fMotifTypes.begin(); 
       it != fMotifTypes.end(); it++) {
    delete it->second;
  }

  for (MotifPositionMapIterator ip=fMotifPositions.begin(); 
       ip != fMotifPositions.end(); ip++) {
    delete ip->second;
  }  
#endif  
}

// 
// private methods
//

//_____________________________________________________________________________
void  AliMpMotifMap::PrintMotif(const AliMpVMotif* motif) const
{
/// Print the motif.

  cout << motif->GetID().Data() << "  "
       << motif->GetMotifType()->GetID() << "    "
       << motif->Dimensions().X() << " "
       << motif->Dimensions().Y();
}

//_____________________________________________________________________________
void  AliMpMotifMap::PrintMotifType(const AliMpMotifType* motifType) const
{
/// Print the motif type.

  cout << motifType->GetID().Data() << "  "
       << motifType->GetNofPadsX() << "  " 
       << motifType->GetNofPadsY() << "  ";
}

//_____________________________________________________________________________
void  AliMpMotifMap::PrintMotifPosition(
                          const AliMpMotifPosition* motifPosition) const
{
/// Print the motif position.

  cout << " ID " << motifPosition->GetID() << "  "
       << " Motif ID " << motifPosition->GetMotif()->GetID() << "  " 
       << " Pos (X,Y) = (" << motifPosition->Position().X() << ","
       << motifPosition->Position().Y() << ")";
}

//_____________________________________________________________________________
void  AliMpMotifMap::PrintMotifPosition2(
                          const AliMpMotifPosition* motifPosition) const
{
/// Print the motif position.

  cout << setw(3) << motifPosition->GetLowIndicesLimit().GetFirst() << "  "
       << setw(3) << motifPosition->GetLowIndicesLimit().GetSecond() << "  "
       << setw(3) << motifPosition->GetHighIndicesLimit().GetFirst()  << " " 
       << setw(3) << motifPosition->GetHighIndicesLimit().GetSecond()  << " "
       << motifPosition->GetID() << "  ";
}

//_____________________________________________________________________________
void  AliMpMotifMap::PrintMotifs() const
{
/// Print all the motifs and their motif types 
/// for all motifs in the motifs map.

#ifdef WITH_STL
  if (fMotifs.size()) {
    cout << "Dump of Motif Map - " << fMotifs.size() << " entries:" << endl;
    Int_t counter = 0;        
    for (MotifMapIterator i=fMotifs.begin(); i != fMotifs.end(); i++) {
      const TString& id  = (*i).first;
      cout << "Map element " 
           << setw(3) << counter++ << "   " 
           << id.Data() << "   " ;
      PrintMotif((*i).second);	   
      cout << endl;
    }
    cout << endl;
  }
#endif

#ifdef WITH_ROOT
  if (fMotifs.GetSize()) {
    cout << "Dump of Motif Map - " << fMotifs.GetSize() << " entries:" << endl;
    Int_t counter = 0;        
    TExMapIter i = fMotifs.GetIterator();
    Long_t key, value;
    while ( i.Next(key, value) ) {
      TString id  = fMotifs.AliMpExMap::GetString(key);
      AliMpVMotif* motif = (AliMpVMotif*)value;
      cout << "Map element " 
           << setw(3) << counter++ << "   " 
           << id.Data() << "   " ;
      PrintMotif(motif);	   
      cout << endl;
    }
    cout << endl;
  }
#endif  
}

//_____________________________________________________________________________
void  AliMpMotifMap::PrintMotifTypes() const
{
/// Print all the the motifs types and their motif dimensions
/// for all motif types in the motif types map.

#ifdef WITH_STL
  if (fMotifTypes.size()) {
    cout << "Dump of Motif Type Map - " << fMotifTypes.size() << " entries:" << endl;
    Int_t counter = 0;        
    for (MotifTypeMapIterator i=fMotifTypes.begin(); i != fMotifTypes.end(); i++) {
      const TString& id  = (*i).first;
      cout << "Map element " 
           << setw(3) << counter++ << "   " 
           << id.Data() << "   ";
      PrintMotifType((*i).second);	   
      cout << endl;
    }
    cout << endl;
  }
#endif  

#ifdef WITH_ROOT
  if (fMotifTypes.GetSize()) {
    cout << "Dump of Motif Type Map - " << fMotifTypes.GetSize() << " entries:" << endl;
    Int_t counter = 0;        
    TExMapIter i = fMotifTypes.GetIterator();
    Long_t key, value;
    while ( i.Next(key, value) ) {
      TString id  = AliMpExMap::GetString(key);
      AliMpMotifType* motifType = (AliMpMotifType*)value;
      cout << "Map element " 
           << setw(3) << counter++ << "   " 
           << id.Data() << "   " ;
      PrintMotifType(motifType);	   
      cout << endl;
    }
    cout << endl;
  }
#endif  
}

//_____________________________________________________________________________
void 
AliMpMotifMap::GetAllMotifPositionsIDs(TArrayI& ecn) const
{
/// Fill the given array with all motif positions IDs (electronic card numbers)
/// defined in the map

#ifdef WITH_STL
  ecn.Set(fMotifPositions.size());  
  Int_t i(0);
  MotifPositionMapIterator it;
  for (it=fMotifPositions.begin(); it != fMotifPositions.end(); it++) {
    AliMpMotifPosition* motifPosition = (*it).second;
    ecn[i++] = motifPosition->GetID();
  }
#endif
  
#ifdef WITH_ROOT  
  ecn.Set(fMotifPositions.GetSize());
  TExMapIter it = fMotifPositions.GetIterator();
  Long_t key, value;
  Int_t i(0);
  
  while ( it.Next(key, value) ) 
  {
    AliMpMotifPosition* motifPosition = reinterpret_cast<AliMpMotifPosition*>(value);
    ecn[i] = motifPosition->GetID();
    ++i;
  }
  
#endif  
}

//_____________________________________________________________________________
UInt_t  AliMpMotifMap::GetNofMotifPositions() const
{
/// Return the number of all motif positions IDs (electronic card numbers)

#ifdef WITH_STL
  return fMotifPositions.size();  
#endif
  
#ifdef WITH_ROOT  
  return fMotifPositions.GetSize();
#endif 
} 

//_____________________________________________________________________________
AliMpMotifPosition* AliMpMotifMap::GetMotifPosition(UInt_t index) const
{
/// Return the motif position which is in the map on the index-th position

  if ( index >= GetNofMotifPositions() ) {
    AliErrorStream() << "Index " << index << " outside limits." << endl;
    return 0;
  }   

#ifdef WITH_STL
  MotifPositionMapIterator it = fMotifPositions.begin();
  std::advance(it, index);
  return it->second;
#endif
  
#ifdef WITH_ROOT  
  return (AliMpMotifPosition*)fMotifPositions.GetObject(index);
#endif 
}

//_____________________________________________________________________________
Int_t AliMpMotifMap::CalculateNofPads() const 
{
/// Calculate total number of pads in the map

  Int_t nofPads = 0;

#ifdef WITH_STL
  MotifPositionMapIterator it;
  for (it=fMotifPositions.begin(); it != fMotifPositions.end(); it++) {
    AliMpMotifPosition* motifPosition = (*it).second;
    nofPads += motifPosition->GetMotif()->GetMotifType()->GetNofPads();
  }
#endif
  
#ifdef WITH_ROOT  
  TExMapIter it = fMotifPositions.GetIterator();
  Long_t key, value;
  
  while ( it.Next(key, value) ) {
    AliMpMotifPosition* motifPosition = reinterpret_cast<AliMpMotifPosition*>(value);
    nofPads += motifPosition->GetMotif()->GetMotifType()->GetNofPads();
  }
#endif  

  return nofPads;
}

//_____________________________________________________________________________
void  AliMpMotifMap::PrintMotifPositions() const
{
/// Print all motif positions.

#ifdef WITH_STL
  if (fMotifPositions.size()) {
    cout << "Dump of Motif Position Map - " << fMotifPositions.size() << " entries:" << endl;
    Int_t counter = 0;        
    for (MotifPositionMapIterator i=fMotifPositions.begin(); 
                                  i != fMotifPositions.end(); i++) {

      cout << "Map element " 
           << setw(3) << counter++ << "   "; 
      PrintMotifPosition((*i).second);	   
      cout << endl;
    }
    cout << endl;
  }
#endif  

#ifdef WITH_ROOT
  if (fMotifPositions.GetSize()) {
    cout << "Dump of Motif Position Map - " << fMotifPositions.GetSize() << " entries:" << endl;
    Int_t counter = 0;        
    TExMapIter i = fMotifPositions.GetIterator();
    Long_t key, value;
    while ( i.Next(key, value) ) {
      AliMpMotifPosition* motifPosition = (AliMpMotifPosition*)value;
      cout << "Map element " 
           << setw(3) << counter++ << "   "; 
      PrintMotifPosition(motifPosition);	   
      cout << endl;
    }
    cout << endl;
  }
#endif  
}

//_____________________________________________________________________________
void  AliMpMotifMap::PrintMotifPositions2() const
{
/// Print all motif positions from the second map
/// (by global indices)

#ifdef WITH_STL
  if (fMotifPositions2.size()) {
    cout << "Dump of Motif Position Map 2 - " << fMotifPositions2.size() << " entries:" << endl;
    Int_t counter = 0;        
    for (MotifPositionMap2Iterator i=fMotifPositions2.begin(); 
                                   i != fMotifPositions2.end(); i++) {

      cout << "Map element " 
           << setw(3) << counter++ << "   "; 
      PrintMotifPosition2((*i).second);  
      cout << endl;
    }
    cout << endl;
  }
#endif  

#ifdef WITH_ROOT
  if (fMotifPositions2.GetSize()) {
    cout << "Dump of Motif Position Map 2 - " << fMotifPositions2.GetSize() << " entries:" << endl;
    Int_t counter = 0;        
    TExMapIter i = fMotifPositions2.GetIterator();
    Long_t key, value;
    while ( i.Next(key, value) ) {
      AliMpMotifPosition* motifPosition = (AliMpMotifPosition*)value;
      cout << "Map element " 
           << setw(3) << counter++ << "   "; 
      PrintMotifPosition2(motifPosition);	   
      cout << endl;
    }
    cout << endl;
  }
#endif  
}

//
// public methods
//

//_____________________________________________________________________________
Bool_t AliMpMotifMap::AddMotif(AliMpVMotif* motif, Bool_t warn)
{
/// Add the specified motif 
/// if the motif with this ID is not yet present.

  AliMpVMotif* found = FindMotif(motif->GetID());
  if (found) {    
    if (warn && found == motif) 
      AliWarningStream() << "The motif is already in map." << endl;

    if (warn && found != motif) {
      AliWarningStream() 
        << "Another motif with the same ID is already in map." << endl; 
    }	     
    return false;
  }  

#ifdef WITH_STL
  fMotifs[motif->GetID()] = motif;
#endif

#ifdef WITH_ROOT
  fMotifs.Add(motif->GetID(), motif);
#endif

  return true;
}

//_____________________________________________________________________________
Bool_t AliMpMotifMap::AddMotifType(AliMpMotifType* motifType, Bool_t warn)
{
/// Add the specified motif type
/// if the motif with this ID is not yet present.

  AliMpMotifType* found = FindMotifType(motifType->GetID());
  if (found) {    
    if (warn && found == motifType) 
      AliWarningStream() << "The motif type is already in map." << endl;
      
    if (warn && found != motifType) { 
      AliWarningStream() 
        << "Another motif type with the same ID is already in map." << endl;
    }	     
    return false;
  }  

#ifdef WITH_STL
  fMotifTypes[motifType->GetID()] = motifType;
#endif

#ifdef WITH_ROOT
  fMotifTypes.Add(motifType->GetID(), motifType);
#endif

  return true;
}

//_____________________________________________________________________________
Bool_t AliMpMotifMap::AddMotifPosition(AliMpMotifPosition* motifPosition, Bool_t warn)
{
/// Add the specified motif position
/// if this position is not yet present.

  AliMpMotifPosition* found = FindMotifPosition(motifPosition->GetID());
  if (found) { 
    if (warn && found == motifPosition) {
      AliWarningStream()
           << "ID: " << motifPosition->GetID() 
           << "  found: " << found 
	   << "  new:   " << motifPosition << endl
	   << "This motif position is already in map." << endl;
    }  
    
    if (warn && found != motifPosition) { 
      AliWarningStream()
           << "ID: " << motifPosition->GetID() 
           << "  found: " << found 
	   << "  new:   " << motifPosition << endl
	   << "Another motif position with the same ID is already in map."
	   << endl;
    }
    	            
    return false;
  }  

#ifdef WITH_STL
  fMotifPositions[motifPosition->GetID()] = motifPosition;
#endif

#ifdef WITH_ROOT
  fMotifPositions.Add(motifPosition->GetID(), motifPosition);
#endif

  return true;
}

//_____________________________________________________________________________
void AliMpMotifMap::FillMotifPositionMap2()
{
/// Fill the second map (by global indices) of motif positions.

#ifdef WITH_STL
  if (fMotifPositions2.size() > 0 ) {
    AliWarningStream() << "Map has been already filled." << endl;
    return;
  }  

  for (MotifPositionMapIterator ip=fMotifPositions.begin(); 
       ip != fMotifPositions.end(); ip++) {

    fMotifPositions2[(*ip).second->GetLowIndicesLimit()] = (*ip).second;
  }  
#endif

#ifdef WITH_ROOT
  if (fMotifPositions2.GetSize() > 0 ) {
    AliWarningStream() <<"Map has been already filled." << endl;
    return;
  }  

  TExMapIter i = fMotifPositions.GetIterator();
  Long_t key, value;
  while ( i.Next(key, value) ) {
    AliMpMotifPosition* motifPosition = (AliMpMotifPosition*)value;
    fMotifPositions2.Add(motifPosition->GetLowIndicesLimit(), motifPosition);
  }
#endif

}

//_____________________________________________________________________________
void  AliMpMotifMap::Print(const char* opt) const
{
/// Print the motifs and motif types maps.

  TString sopt(opt);
  
  sopt.ToUpper();
  
  if ( sopt.Contains("MOTIFS") || sopt == "ALL" ) PrintMotifs();
  if ( sopt.Contains("MOTIFTYPES") || sopt == "ALL" ) PrintMotifTypes();
  if ( sopt.Contains("MOTIFPOSITIONS") || sopt == "ALL" ) PrintMotifPositions();
  if ( sopt.Contains("MOTIFPOSITIONS2") || sopt == "ALL" ) PrintMotifPositions2();
}

//_____________________________________________________________________________
void  AliMpMotifMap::PrintGlobalIndices(const char* fileName) const
{
/// Print all motif positions and their global indices.

  ofstream out(fileName, ios::out);

#ifdef WITH_STL
  if (fMotifPositions.size()) {
    for (MotifPositionMapIterator i=fMotifPositions.begin(); 
                                   i != fMotifPositions.end(); i++) {

      AliMpMotifPosition* motifPosition = (*i).second;
      out << setw(5) << motifPosition->GetID() << "     "
	  << setw(3) << motifPosition->GetLowIndicesLimit().GetFirst()  << " " 
	  << setw(3) << motifPosition->GetLowIndicesLimit().GetSecond() 
         << endl;
    }
    out << endl;
  }
#endif

#ifdef WITH_ROOT
  if (fMotifPositions.GetSize()) {
    TExMapIter i = fMotifPositions.GetIterator();
    Long_t key, value;
    while ( i.Next(key, value) ) {
      AliMpMotifPosition* motifPosition = (AliMpMotifPosition*)value;
      out << setw(5) << motifPosition->GetID() << "     "
	  << setw(3) << motifPosition->GetLowIndicesLimit().GetFirst()  << " " 
	  << setw(3) << motifPosition->GetLowIndicesLimit().GetSecond() 
         << endl;
    }
    out << endl;
  }
#endif
}

//_____________________________________________________________________________
void  AliMpMotifMap::UpdateGlobalIndices(const char* fileName)
{
/// Update the motif positions global indices from the file.

  ifstream in(fileName, ios::in);

  Int_t motifPositionId, offx, offy;
    
  do {
    in >> motifPositionId >> offx >> offy;
    
    if (in.eof()) {
      FillMotifPositionMap2();
      return;
    }  
    
    AliMpMotifPosition* motifPosition = FindMotifPosition(motifPositionId);
	  
    if (motifPosition) {
       AliDebugStream(1) 
            << "Processing " 
            << motifPosition->GetID() << " " << offx << " " << offy << endl; 

       motifPosition->SetLowIndicesLimit(AliMpIntPair(offx, offy));
       
       Int_t offx2 
         = offx + motifPosition->GetMotif()->GetMotifType()->GetNofPadsX() - 1;
	 
       Int_t offy2 
         = offy + motifPosition->GetMotif()->GetMotifType()->GetNofPadsY() - 1;
       
       motifPosition->SetHighIndicesLimit(AliMpIntPair(offx2, offy2));
    }
    else {   
       AliWarningStream()
         << "Motif position " << motifPositionId << " not found" << endl;
    }
  }    
  while (!in.eof());
}


//_____________________________________________________________________________
AliMpVMotif* AliMpMotifMap::FindMotif(const TString& motifID) const
{
/// Find the motif with the specified ID.
  
#ifdef WITH_STL
  MotifMapIterator i = fMotifs.find(motifID);
  if (i != fMotifs.end()) 
    return (*i).second;
  else                 
    return 0;
#endif

#ifdef WITH_ROOT
  return (AliMpVMotif*)fMotifs.GetValue(motifID);
#endif
}

//_____________________________________________________________________________
AliMpVMotif* AliMpMotifMap::FindMotif(const TString& motifID, 
                                      const TString& motifTypeID,
			              const TVector2& padDimensions ) const
{
/// Find the motif with the specified ID and returns it
/// only if its motif type and motif dimensions agree
/// with the given motifTypeID and motifDimensions.
/// Disagreement causes fatal error.
 
  AliMpVMotif* motif = FindMotif(motifID);

  if (motif && motif->GetMotifType()->GetID() != motifTypeID) {
      Fatal("FindMotif", 
            "Motif has been already defined with a different type.");
      return 0;	    
  }

  // check pad dimension in case of a normal motif
  if (motif && 
      dynamic_cast<AliMpMotif*>(motif) && 
      ( motif->GetPadDimensions(0).X() != padDimensions.X() ||
        motif->GetPadDimensions(0).Y() != padDimensions.Y())) { 
      
      Fatal("FindMotifType", 
            "Motif type has been already defined with different dimensions.");
      return 0;

  } 

  // check case of a special motif
  if (motif && 
      (padDimensions.X() == 0. && padDimensions.Y() == 0.) &&
      !dynamic_cast<AliMpMotifSpecial*>(motif)) {

      Fatal("FindMotifType", 
            "Motif type has been already defined with different dimensions.");
      return 0;

  } 
  
  return motif;
}

//_____________________________________________________________________________
AliMpMotifType* AliMpMotifMap::FindMotifType(const TString& motifTypeID) const
{
/// Find the motif type with the specified motif type ID.
  
#ifdef WITH_STL
  MotifTypeMapIterator i = fMotifTypes.find(motifTypeID);
  if (i != fMotifTypes.end()) 
    return (*i).second;
  else                 
    return 0;
#endif

#ifdef WITH_ROOT
  return (AliMpMotifType*)fMotifTypes.GetValue(motifTypeID);
#endif
}

//_____________________________________________________________________________
AliMpMotifPosition* 
AliMpMotifMap::FindMotifPosition(Int_t motifPositionID) const
{
/// Find the motif position with the specified motif position ID.
  
#ifdef WITH_STL
  MotifPositionMapIterator i = fMotifPositions.find(motifPositionID);
  if (i != fMotifPositions.end()) 
    return (*i).second;
  else                 
    return 0;
#endif

#ifdef WITH_ROOT
  return (AliMpMotifPosition*)fMotifPositions.GetValue(motifPositionID);
#endif
}

/*
//_____________________________________________________________________________
AliMpMotifPosition* 
AliMpMotifMap::FindMotifPosition(const AliMpIntPair& indices) const
{
/// Find the last motif position which has the global indices (low limit)
/// less then the indices specified.

#ifdef WITH_STL
  MotifPositionMap2Iterator found 
    = fMotifPositions2.lower_bound(indices);
  
  if (found == fMotifPositions2.end()) found--; 

  MotifPositionMap2Iterator i=found;
  do {
    AliMpIntPair low = (*i).second->GetLowIndicesLimit();
    AliMpIntPair up = (*i).second->GetHighIndicesLimit();
    
    if ( indices.GetFirst()  >= low.GetFirst() &&
         indices.GetSecond() >= low.GetSecond() &&
	 indices.GetFirst()  <= up.GetFirst() &&
         indices.GetSecond() <= up.GetSecond())
	 
	 return (*i).second;    	 	 
  }
  while ( i-- != fMotifPositions2.begin());
  
  return 0;
#endif

#ifdef WITH_ROOT
  // HOW TO DO THIS WITH ROOT ????
  // Fortunately it seems not to be used anywhere
  Fatal("FindMotifPosition", "Difficult in Root to do this.");
  return 0;
#endif
}
*/
