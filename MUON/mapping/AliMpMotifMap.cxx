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
// $MpId: AliMpMotifMap.cxx,v 1.9 2005/09/26 16:11:20 ivana Exp $
// Category: motif
//
// Class AliMpMotifMap
// -------------------
// Class describing the motif map container, where motifs are
// mapped to their string IDs.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <Riostream.h>
#include <TVector2.h>

#include "AliMpMotifMap.h"
#include "AliMpVMotif.h"
#include "AliMpMotif.h"
#include "AliMpMotifSpecial.h"
#include "AliMpMotifType.h"
#include "AliMpMotifPosition.h"

ClassImp(AliMpMotifMap)

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
// ---

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

  cout << motifPosition->GetID() << "  "
       << motifPosition->GetMotif()->GetID() << "  " 
       << motifPosition->Position().X() << "  "
       << motifPosition->Position().Y() << "    ";
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
void  AliMpMotifMap::PrintMotifPositions() const
{
/// Print all the the motifs positions.

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
/// Print all the the motifs positions from the second map
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
      Warning("AddMotif", "The motif is already in map.");
    if (warn && found != motif) 
      Warning("AddMotif", "Another motif with the same ID is already in map.");      
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
      Warning("AddMotifType", "The motif type is already in map.");
    if (warn && found != motifType) 
      Warning("AddMotifType", 
              "Another motif type with the same ID is already in map.");      
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
      cerr << "ID: " << motifPosition->GetID() 
           << "  found: " << found 
	   << "  new:   " << motifPosition << endl;   
      Warning("AddMotifPosition", "This motif position is already in map.");
    }  
    if (warn && found != motifPosition) { 
      cerr << "ID: " << motifPosition->GetID() 
           << "  found: " << found 
	   << "  new:   " << motifPosition << endl;   
      Warning("AddMotifposition", 
              "Another motif position with the same ID is already in map.");
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
    Warning("FillMotifPositionMap2", "Map has been already filled.");
    return;
  }  

  for (MotifPositionMapIterator ip=fMotifPositions.begin(); 
       ip != fMotifPositions.end(); ip++) {

    fMotifPositions2[(*ip).second->GetLowIndicesLimit()] = (*ip).second;
  }  
#endif

#ifdef WITH_ROOT
  if (fMotifPositions2.GetSize() > 0 ) {
    Warning("FillMotifPositionMap2", "Map has been already filled.");
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
void  AliMpMotifMap::Print(const char* /*option*/) const
{
/// Print the motifs and motif types maps.

  PrintMotifs();
  PrintMotifTypes();
  PrintMotifPositions();
  PrintMotifPositions2();
}

//_____________________________________________________________________________
void  AliMpMotifMap::PrintGlobalIndices(const char* fileName) const
{
/// Print all the motifs positions and their global indices.

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
/// Updates the motifs positions global indices
/// from the file.

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
       cout << "Processing " 
            << motifPosition->GetID() << " " << offx << " " << offy << endl; 

       motifPosition->SetLowIndicesLimit(AliMpIntPair(offx, offy));
       
       Int_t offx2 
         = offx + motifPosition->GetMotif()->GetMotifType()->GetNofPadsX() - 1;
	 
       Int_t offy2 
         = offy + motifPosition->GetMotif()->GetMotifType()->GetNofPadsY() - 1;
       
       motifPosition->SetHighIndicesLimit(AliMpIntPair(offx2, offy2));
    }
    else {   
       cerr <<"Motif position " << motifPositionId << endl;
       Warning("UpdateGlobalIndices", "Motif position not found !!!");
    }
  }    
  while (!in.eof());
}


//_____________________________________________________________________________
AliMpVMotif* AliMpMotifMap::FindMotif(const TString& motifID) const
{
/// Finds the motif with the specified ID.
  
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
/// Finds the motif with the specified ID and returns it
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
