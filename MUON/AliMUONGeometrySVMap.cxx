// $Id$
//
// Class AliMUONGeometrySVMap
// ------------------------------------ 
// As the detection element frame is different from the
// frame of the sensitive volume(s) defined in Geant,
// the sensitive volumes have to be mapped to the detection 
// elements. In the map, fSVMap, the sensitive voolumes are specified
// by the full path in the volume hierarchy, defined as:
//  /volname.copyNo/volName.copyNo1/...
//
// The array of global positions of sensitive volumes fSVPositions
// is included to make easier the verification of the assignements 
// in the fSVMap.
//
// Author: Ivana Hrivnacova, IPN Orsay

#include <Riostream.h>
#include <TGeoMatrix.h>
#include <TObjString.h>

#include "AliMUONGeometrySVMap.h"
#include "AliLog.h"

ClassImp(AliMUONGeometrySVMap)

//
// Class AliMUONStringIntMap
//

//______________________________________________________________________________
AliMUONStringIntMap::AliMUONStringIntMap()
 : TObject(),
   fNofItems(0),
   fFirstArray(100),
   fSecondArray(100)
{
// Standard constructor

  fFirstArray.SetOwner(true);
}

//______________________________________________________________________________
AliMUONStringIntMap::AliMUONStringIntMap(const AliMUONStringIntMap& rhs)
  : TObject(rhs)
{
  AliFatal("Copy constructor is not implemented.");
}

//______________________________________________________________________________
AliMUONStringIntMap& 
AliMUONStringIntMap::operator = (const AliMUONStringIntMap& rhs) 
{
  // check assignement to self
  if (this == &rhs) return *this;

  AliFatal("Assignment operator is not implemented.");
    
  return *this;  
}

//______________________________________________________________________________
AliMUONStringIntMap::~AliMUONStringIntMap()
{
// Destructor

  fFirstArray.Delete();
}


//______________________________________________________________________________
Bool_t  AliMUONStringIntMap::Add(const TString& first, Int_t second)
{
// Add map element if first not yet present
// ---

  Int_t second2 = Get(first);
  if ( second2 > 0 ) {
    AliError(Form("%s is already present in the map", first.Data()));
    return false;
  }
  
  // Resize TArrayI if needed
  if (fSecondArray.GetSize() == fNofItems) fSecondArray.Set(2*fNofItems);
  
  fFirstArray.Add(new TObjString(first)); 
  fSecondArray.AddAt(second, fNofItems);
  fNofItems++;
   
  return true;
}  

//______________________________________________________________________________
Int_t  AliMUONStringIntMap::Get(const TString& first) const
{
// Find the element with specified key (first)
// ---
  
  for (Int_t i=0; i<fNofItems; i++) {
    if ( ((TObjString*)fFirstArray.At(i))->GetString() == first )
      return fSecondArray.At(i);
  }
  
  return 0;
}      

//______________________________________________________________________________
Int_t  AliMUONStringIntMap::GetNofItems() const
{
// Returns the number of elements
// ---

  return fNofItems;
}  

//______________________________________________________________________________
void  AliMUONStringIntMap::Clear()
{
// Deletes the elements
// ---

  cout << "######### clearing map " << endl;

  fNofItems = 0;
  fFirstArray.Delete();
  fSecondArray.Reset();

  cout << "######### clearing map done " << endl;
}  
    
//______________________________________________________________________________
void AliMUONStringIntMap::Print(const char* /*option*/) const
{
// Prints the map elements

  for (Int_t i=0; i<fNofItems; i++) {
    cout << setw(4)
         << i << "  "
         << ((TObjString*)fFirstArray.At(i))->GetString()
	 << "  "
	 << setw(5)
	 << fSecondArray.At(i)
	 << endl;
  }
}  	 

//______________________________________________________________________________
void AliMUONStringIntMap::Print(const TString& key, ofstream& out) const
{
// Prints the map elements

  for (Int_t i=0; i<fNofItems; i++) {
    out  << key << "  "
         << ((TObjString*)fFirstArray.At(i))->GetString()
	 << "  "
	 << setw(5)
	 << fSecondArray.At(i)
	 << endl;
  }
}  	 


//
// Class AliMUONGeometrySVMap
//

//______________________________________________________________________________
AliMUONGeometrySVMap::AliMUONGeometrySVMap(Int_t initSize)
 : TObject(),
   fSVMap(),
   fSVPositions(initSize)
{ 
// Standard constructor
  
  fSVPositions.SetOwner(true);
}

//______________________________________________________________________________
AliMUONGeometrySVMap::AliMUONGeometrySVMap()
 : TObject(),
   fSVMap(),
   fSVPositions()
{
// Default constructor
}

//______________________________________________________________________________
AliMUONGeometrySVMap::AliMUONGeometrySVMap(const AliMUONGeometrySVMap& rhs)
  : TObject(rhs)
{
  AliFatal("Copy constructor is not implemented.");
}

//______________________________________________________________________________
AliMUONGeometrySVMap::~AliMUONGeometrySVMap() {
//
  fSVPositions.Delete();
}

//______________________________________________________________________________
AliMUONGeometrySVMap& 
AliMUONGeometrySVMap::operator = (const AliMUONGeometrySVMap& rhs) 
{
  // check assignement to self
  if (this == &rhs) return *this;

  AliFatal("Assignment operator is not implemented.");
    
  return *this;  
}

//
// private methods
//

//______________________________________________________________________________
const TGeoCombiTrans* 
AliMUONGeometrySVMap::FindByName(const TString& name) const
{
// Finds TGeoCombiTrans in the array of positions by name 
// ---

  for (Int_t i=0; i<fSVPositions.GetEntriesFast(); i++) { 
     TGeoCombiTrans* transform = (TGeoCombiTrans*) fSVPositions.At(i);
     if ( transform && TString(transform->GetTitle()) == name )
       return transform;
  }     
       
  return 0;
}  


//
// public methods
//

//______________________________________________________________________________
void AliMUONGeometrySVMap::Add(const TString& volumePath, 
                               Int_t detElemId)
{
// Add the specified sensitive volume path and the detElemId 
// to the map
// ---
 
  fSVMap.Add(volumePath, detElemId);
}		          
    
//______________________________________________________________________________
void AliMUONGeometrySVMap::AddPosition(const TString& volumePath, 
                              const TGeoTranslation& globalPosition)
{
// Add global position for the sensitive volume specified by volumePath  
// in the array of transformations if this volumePath is not yet present. 
// ---
 
  TGeoTranslation* newTransform = new TGeoTranslation(globalPosition);
  Int_t detElemId = fSVMap.Get(volumePath);

  TString detElemIdString("");
  detElemIdString += detElemId;

  newTransform->SetName(detElemIdString);
  newTransform->SetTitle(volumePath);
  
  // cout << ".. adding " << volumePath << "  " << detElemId << endl;

  // Add to the map  
  if ( !FindByName(volumePath )) {
  
    newTransform->SetUniqueID(detElemId);
      // Set detector element id as unique id
 
    fSVPositions.Add(newTransform);
  } 
}		      
    
//______________________________________________________________________________
void AliMUONGeometrySVMap::Clear()
{
// Clears the sensitive volumes map

  fSVMap.Clear();
}  

//______________________________________________________________________________
void AliMUONGeometrySVMap::ClearPositions()
{
// Clears the array of transformations

  fSVPositions.Delete();
}  

//______________________________________________________________________________
void AliMUONGeometrySVMap::SortPositions()
{
// Sort the array of positions by names.
// ---

  fSVPositions.Sort(fSVPositions.GetEntriesFast());
}
  
//______________________________________________________________________________
void  AliMUONGeometrySVMap::Print(const char* option) const
{    
// Prints the map of sensitive volumes and detector elements 
// ---

  fSVMap.Print(option);
}  

//______________________________________________________________________________
void  AliMUONGeometrySVMap::PrintPositions() const
{
// Prints the sensitive volumes global positions
// ---

  for (Int_t i=0; i<fSVPositions.GetEntriesFast(); i++) {
    
    TGeoTranslation* matrix = (TGeoTranslation*)fSVPositions.At(i);

    cout << "DetElemId: " << matrix->GetUniqueID();
    cout << "  name: " << matrix->GetTitle() << endl;

    const double* translation = matrix->GetTranslation();
    cout << "   translation: "
         << std::fixed
         << std::setw(7) << std::setprecision(4) << translation[0] << ", " 
         << std::setw(7) << std::setprecision(4) << translation[1] << ", "
         << std::setw(7) << std::setprecision(4) << translation[2] << endl;
  }
}     

//______________________________________________________________________________
void  AliMUONGeometrySVMap::WriteMap(ofstream& out) const
{    
// Prints the map of sensitive volumes and detector elements 
// into specified stream
// ---

  fSVMap.Print("SV", out);
}  

//______________________________________________________________________________
Int_t  AliMUONGeometrySVMap::GetDetElemId(const TString& volumePath) const
{
// Returns detection element Id for the sensitive volume specified by path
// ---

  return fSVMap.Get(volumePath);
}  
