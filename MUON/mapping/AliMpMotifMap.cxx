// $Id$
// Category: motif
//
// Class AliMpMotifMap
// -------------------
// Class describing the motif map container, where motifs are
// mapped to their string IDs.
//
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
AliMpMotifMap::AliMpMotifMap() 
  : TObject()
{
//
}

//_____________________________________________________________________________
AliMpMotifMap::~AliMpMotifMap() {
//  

  // Delete all registered motifs, motif types, motif positions
  
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
}

// 
// private methods
//

//_____________________________________________________________________________
void  AliMpMotifMap::PrintMotifs() const
{
// Prints all the motifs and their motif types 
// for all motifs in the motifs map.
// ---

  if (fMotifs.size()) {
    cout << "Dump of Motif Map - " << fMotifs.size() << " entries:" << endl;
    Int_t counter = 0;        
    for (MotifMapIterator i=fMotifs.begin(); i != fMotifs.end(); i++) {
      const TString& id  = (*i).first;
      AliMpVMotif* motif = (*i).second;
      cout << "Map element " 
           << setw(3) << counter++ << "   " 
           << id.Data() << "   " 
	   << motif->GetID().Data() << "  "
	   << motif->GetMotifType()->GetID() << "    "
	   << motif->Dimensions().X() << " "
	   << motif->Dimensions().Y()
	   << endl;
    }
    cout << endl;
  }
}

//_____________________________________________________________________________
void  AliMpMotifMap::PrintMotifTypes() const
{
// Prints all the the motifs types and their motif dimensions
// for all motif types in the motif types map.
// ---

  if (fMotifTypes.size()) {
    cout << "Dump of Motif Type Map - " << fMotifTypes.size() << " entries:" << endl;
    Int_t counter = 0;        
    for (MotifTypeMapIterator i=fMotifTypes.begin(); i != fMotifTypes.end(); i++) {
      const TString& id  = (*i).first;
      AliMpMotifType* motifType = (*i).second;
      cout << "Map element " 
           << setw(3) << counter++ << "   " 
           << id.Data() << "   " 
	   << motifType->GetID().Data() << "  "
	   << motifType->GetNofPadsX() << "  " 
	   << motifType->GetNofPadsY() << "  "
	   << endl;
    }
    cout << endl;
  }
}

//_____________________________________________________________________________
void  AliMpMotifMap::PrintMotifPositions() const
{
// Prints all the the motifs positions.
// ---

  if (fMotifPositions.size()) {
    cout << "Dump of Motif Position Map - " << fMotifPositions.size() << " entries:" << endl;
    Int_t counter = 0;        
    for (MotifPositionMapIterator i=fMotifPositions.begin(); 
                                  i != fMotifPositions.end(); i++) {

      AliMpMotifPosition* motifPosition = (*i).second;
      cout << "Map element " 
           << setw(3) << counter++ << "   " 
	   << motifPosition->GetID() << "  "
	   << motifPosition->GetMotif()->GetID() << "  " 
	   << motifPosition->Position().X() << "  "
	   << motifPosition->Position().Y() << "    "
	   << endl;
    }
    cout << endl;
  }
}

//_____________________________________________________________________________
void  AliMpMotifMap::PrintMotifPositions2() const
{
// Prints all the the motifs positions from the second map
// (by global indices)
// ---

  if (fMotifPositions2.size()) {
    cout << "Dump of Motif Position Map 2 - " << fMotifPositions2.size() << " entries:" << endl;
    Int_t counter = 0;        
    for (MotifPositionMap2Iterator i=fMotifPositions2.begin(); 
                                   i != fMotifPositions2.end(); i++) {

      AliMpMotifPosition* motifPosition = (*i).second;
      cout << "Map element " 
           << setw(3) << counter++ << "   " 
	   << setw(3) << motifPosition->GetLowIndicesLimit().GetFirst() << "  "
	   << setw(3) << motifPosition->GetLowIndicesLimit().GetSecond() << "  "
 	   << setw(3) << motifPosition->GetHighIndicesLimit().GetFirst()  << " " 
	   << setw(3) << motifPosition->GetHighIndicesLimit().GetSecond()  << " "
	   << motifPosition->GetID() << "  "
	   << endl;
    }
    cout << endl;
  }
}

//
// public methods
//

//_____________________________________________________________________________
Bool_t AliMpMotifMap::AddMotif(AliMpVMotif* motif, Bool_t warn)
{
// Adds the specified motif 
// if the motif with this ID is not yet present.
// ---

  AliMpVMotif* found = FindMotif(motif->GetID());
  if (found) {    
    if (warn && found == motif) 
      Warning("AddMotif", "The motif is already in map.");
    if (warn && found != motif) 
      Warning("AddMotif", "Another motif with the same ID is already in map.");      
    return false;
  }  

  fMotifs[motif->GetID()] = motif;
  return true;
}

//_____________________________________________________________________________
Bool_t AliMpMotifMap::AddMotifType(AliMpMotifType* motifType, Bool_t warn)
{
// Adds the specified motif type
// if the motif with this ID is not yet present.
// ---

  AliMpMotifType* found = FindMotifType(motifType->GetID());
  if (found) {    
    if (warn && found == motifType) 
      Warning("AddMotifType", "The motif type is already in map.");
    if (warn && found != motifType) 
      Warning("AddMotifType", 
              "Another motif type with the same ID is already in map.");      
    return false;
  }  

  fMotifTypes[motifType->GetID()] = motifType;
  return true;
}

//_____________________________________________________________________________
Bool_t AliMpMotifMap::AddMotifPosition(AliMpMotifPosition* motifPosition, Bool_t warn)
{
// Adds the specified motif position
// if this position is not yet present.
// ---

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

  fMotifPositions[motifPosition->GetID()] = motifPosition;
  return true;
}

//_____________________________________________________________________________
void AliMpMotifMap::FillMotifPositionMap2()
{
// Fills the second map (by global indices) of motif positions.
// ---

  if (fMotifPositions2.size() > 0 ) {
    Warning("FillMotifPositionMap2", "Map has been already filled.");
    return;
  }  

  for (MotifPositionMapIterator ip=fMotifPositions.begin(); 
       ip != fMotifPositions.end(); ip++) {

    fMotifPositions2[(*ip).second->GetLowIndicesLimit()] = (*ip).second;
  }  

}

//_____________________________________________________________________________
void  AliMpMotifMap::Print(const char* /*option*/) const
{
// Prints the motifs and motif types maps.
// ---

  PrintMotifs();
  PrintMotifTypes();
  PrintMotifPositions();
  PrintMotifPositions2();
}

//_____________________________________________________________________________
void  AliMpMotifMap::PrintGlobalIndices(const char* fileName) const
{
// Prints all the motifs positions and their global indices.
// ---

  ofstream out(fileName, ios::out);

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
}

//_____________________________________________________________________________
void  AliMpMotifMap::UpdateGlobalIndices(const char* fileName)
{
// Updates the motifs positions global indices
// from the file.
// ---

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
// Finds the motif with the specified ID.
// ---
  
  MotifMapIterator i = fMotifs.find(motifID);
  
  if (i != fMotifs.end()) 
    return (*i).second;
  else                 
    return 0;
}

//_____________________________________________________________________________
AliMpVMotif* AliMpMotifMap::FindMotif(const TString& motifID, 
                                      const TString& motifTypeID,
			              const TVector2& padDimensions ) const
{
// Finds the motif with the specified ID and returns it
// only if its motif type and motif dimensions agree
// with the given motifTypeID and motifDimensions.
// Disagreement causes fatal error.
// ---
  
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
// Finds the motif type with the specified motif type ID.
// ---
  
  MotifTypeMapIterator i = fMotifTypes.find(motifTypeID);
  
  if (i != fMotifTypes.end()) 
    return (*i).second;
  else                 
    return 0;
}

//_____________________________________________________________________________
AliMpMotifPosition* 
AliMpMotifMap::FindMotifPosition(Int_t motifPositionID) const
{
// Finds the motif position with the specified motif position ID.
// ---
  
  MotifPositionMapIterator i = fMotifPositions.find(motifPositionID);
  
  if (i != fMotifPositions.end()) 
    return (*i).second;
  else                 
    return 0;
}

//_____________________________________________________________________________
AliMpMotifPosition* 
AliMpMotifMap::FindMotifPosition(const AliMpIntPair& indices) const
{
// Finds the last motif position which has the global indices (low limit)
// less then the indices specified.
// ---

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
}
