// $Id$
// Category: plane
//
// Class AliMpPlaneSegmentation
// ----------------------------
// Class describing the segmentation of the plane.
//
// Transformation of pad characteristics according to sectors:
//   I.  ( posId,  Guassi ), ( i, j), ( x, y)         II. |  I.
//  II.  ( posId', Guassi'), (-i, j), (-x, y)       _____ | ____
// III.  (-posId,  Guassi),  (-i,-j), (-x,-y)             |
//  IV.  (-posId', Guassi'), ( i,-j), ( x,-y)        III. |  IV.
//   
// Where (posId', Guassi') is the location of the pad
// in the clipped sector.
//
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <Riostream.h>
#include <TMath.h>
#include <TError.h>

#include "AliMpPlaneSegmentation.h"
#include "AliMpPlaneAreaPadIterator.h"
#include "AliMpPlane.h"
#include "AliMpSectorPosition.h"
#include "AliMpSectorSegmentation.h"

ClassImp(AliMpPlaneSegmentation)

//_____________________________________________________________________________
AliMpPlaneSegmentation::AliMpPlaneSegmentation(const AliMpPlane* plane) 
  : AliMpVSegmentation(),
    fkPlane(plane),
    fFrontSectorSegmentation(0),
    fBackSectorSegmentation(0)
{
//
  fFrontSectorSegmentation = new AliMpSectorSegmentation(plane->GetFrontSector());
  fBackSectorSegmentation = new AliMpSectorSegmentation(plane->GetBackSector());
  
  for (Int_t i=0; i<fkPlane->GetNofSectorPositions(); i++) {

#ifdef WITH_STL
    fTransformers.push_back(
      new AliMpTransformer(fkPlane->GetSectorPosition(i)->GetOffset(),
                           fkPlane->GetSectorPosition(i)->GetScale()));
#endif

#ifdef WITH_ROOT
    fTransformers.Add(
      new AliMpTransformer(fkPlane->GetSectorPosition(i)->GetOffset(),
                           fkPlane->GetSectorPosition(i)->GetScale()));
#endif
 }		       
}

///_____________________________________________________________________________
AliMpPlaneSegmentation::AliMpPlaneSegmentation() 
  : AliMpVSegmentation(),
    fkPlane(0),
    fFrontSectorSegmentation(0),
    fBackSectorSegmentation(0)
{
//
}

//_____________________________________________________________________________
AliMpPlaneSegmentation::AliMpPlaneSegmentation(
                                  const AliMpPlaneSegmentation& right) 
  : AliMpVSegmentation(right) {
// 
  Fatal("AliMpPlaneSegmentation", "Copy constructor not provided.");
}

//_____________________________________________________________________________
AliMpPlaneSegmentation::~AliMpPlaneSegmentation() {
// 
  delete fFrontSectorSegmentation;
  delete fBackSectorSegmentation;

  for (Int_t i=0; i<GetNofTransformers(); i++) 
    delete GetTransformer(i);
}

//
// operators
//

//_____________________________________________________________________________
AliMpPlaneSegmentation& 
AliMpPlaneSegmentation::operator=(const AliMpPlaneSegmentation& right)
{
  // check assignement to self
  if (this == &right) return *this;

  Fatal("operator =", "Assignement operator not provided.");
    
  return *this;  
}    

//
// private methods
//

//_____________________________________________________________________________
const AliMpTransformer* 
AliMpPlaneSegmentation::GetTransformer(const AliMpIntPair& scale) const
{
// Returns the sector position specified by scale.
// ---

  for (Int_t i=0; i<GetNofTransformers(); i++) 
    if (GetTransformer(i)->GetScale() == scale) return GetTransformer(i);

  Fatal("GetTransformer", "Wrong scale");
  return 0; 
}

//_____________________________________________________________________________
AliMpIntPair AliMpPlaneSegmentation::GetScale(const AliMpIntPair& pair) const
{
// Returns pair of the signs of the values of the given pair.
// ---

  AliMpIntPair scale(1, 1);
  
  if (pair.GetFirst() < 0)  scale.SetFirst(-1);
  if (pair.GetSecond() < 0) scale.SetSecond(-1);
  
  return scale;
}

//_____________________________________________________________________________
AliMpIntPair AliMpPlaneSegmentation::GetScale(const TVector2& vector) const
{
// Returns pair of the signs of the values of the given vector.
// ---

  AliMpIntPair scale(1, 1);
  
  if (vector.X() < 0) scale.SetFirst(-1);
  if (vector.Y() < 0) scale.SetSecond(-1);
  
  return scale;
}

//_____________________________________________________________________________
AliMpIntPair 
AliMpPlaneSegmentation::GetLocationScale(const AliMpIntPair& location) const
{
// Returns the scale transformation of the specified location. 
// ---

  // Find the sector
  Bool_t inFront;
  if (fFrontSectorSegmentation
        ->HasMotifPosition(TMath::Abs(location.GetFirst())))
    inFront = true;
  else if (fBackSectorSegmentation
             ->HasMotifPosition(TMath::Abs(location.GetFirst())))  
    inFront = false;
  else {
    Fatal("GetLocationScale", "Motif position not found.");
    return AliMpIntPair();
  }  
    
  if      (inFront  && location.GetFirst() > 0) return  AliMpIntPair(1, 1); 
  else if (inFront  && location.GetFirst() < 0) return  AliMpIntPair(-1, -1);
  else if (!inFront && location.GetFirst() > 0) return  AliMpIntPair(-1, 1);
  else if (!inFront && location.GetFirst() < 0) return  AliMpIntPair( 1,-1);  

  // cannot get there
  Fatal("GetLocationScale", "Condition failed.");
  return AliMpIntPair();
}


//_____________________________________________________________________________
AliMpSectorSegmentation* 
AliMpPlaneSegmentation::GetSectorSegmentation(const AliMpIntPair& scale) const
{    
// Returns front sector or back sector segmentation
// according to quadrant specified by scale.
// ---

  if (scale.GetFirst()*scale.GetSecond() > 0) {
    // quadrant I or III
    return fFrontSectorSegmentation;
  }  
  else  {
    // quadrant II or IV
    return fBackSectorSegmentation;
  }  
}

//_____________________________________________________________________________
AliMpSectorSegmentation* 
AliMpPlaneSegmentation::GetSectorSegmentation(Int_t motifPositionId) const
{    
// Returns front sector or back sector segmentation
// according to specified motifPositionId
// ---

  if (fFrontSectorSegmentation->HasMotifPosition(motifPositionId))
    return fFrontSectorSegmentation;
  else if (fBackSectorSegmentation->HasMotifPosition(motifPositionId)) 
    return fBackSectorSegmentation;
  else {
    Fatal("GetSectorSegmentation", "Motif position not found.");
    return 0;
  }
}

//
// public methods
//

//_____________________________________________________________________________
AliMpVPadIterator* 
AliMpPlaneSegmentation::CreateIterator(const AliMpArea& area) const
{
// Creates the are iterator. 
// (The inherited method cannot be used)
// ---

  return new AliMpPlaneAreaPadIterator(this, area);  
}   
  
//______________________________________________________________________________
AliMpPad AliMpPlaneSegmentation::PadByLocation(const AliMpIntPair& location, 
                                       Bool_t warning) const
{
// Find the pad which corresponds to the given location
// ---

  // Get segmentation
  AliMpSectorSegmentation* segmentation 
    = GetSectorSegmentation(TMath::Abs(location.GetFirst()));

  // Get pad in the segmentation
  AliMpPad pad
     = segmentation
       ->PadByLocation(
           AliMpIntPair(TMath::Abs(location.GetFirst()),location.GetSecond()), 
	                warning);

  // Get transformation
  AliMpIntPair scale  = GetLocationScale(location);
  const AliMpTransformer* kTransformer = GetTransformer(scale);
  
  // Transform pad characteristics
  return kTransformer->Transform(pad);	      
}

//______________________________________________________________________________
AliMpPad AliMpPlaneSegmentation::PadByIndices (const AliMpIntPair& indices,
                                               Bool_t warning ) const
{
// Find the pad which corresponds to the given indices  
//

  AliMpIntPair scale = GetScale(indices);
  const AliMpTransformer* kTransformer = GetTransformer(scale);

  AliMpIntPair scaledIndices = kTransformer->Scale(indices);
  AliMpPad pad 
    = GetSectorSegmentation(scale)->PadByIndices(scaledIndices, warning);
    
  return kTransformer->Transform(pad);	      
}

//_____________________________________________________________________________
AliMpPad AliMpPlaneSegmentation::PadByPosition(const TVector2& position,
                                               Bool_t warning) const
{
// Find the pad which corresponds to the given position
// ---

  AliMpIntPair scale = GetScale(position);
  const AliMpTransformer* kTransformer = GetTransformer(scale);

  TVector2 scaledPosition = kTransformer->ITransform(position);  
  AliMpPad pad 
    = GetSectorSegmentation(scale)->PadByPosition(scaledPosition, warning);
  
  return kTransformer->Transform(pad);	      
}

//_____________________________________________________________________________
Bool_t AliMpPlaneSegmentation::HasPad(const AliMpIntPair& indices) const
{
// Does the pad located by <indices> exists ?
// ---

  AliMpIntPair scale = GetScale(indices);
  const AliMpTransformer* kTransformer = GetTransformer(scale);

  AliMpIntPair scaledIndices = kTransformer->Scale(indices);

  return GetSectorSegmentation(scale)->HasPad(scaledIndices);
}

//_____________________________________________________________________________
Int_t AliMpPlaneSegmentation::Zone(const AliMpPad& pad, Bool_t warning) const
{
// Returns the zone index of the zone containing the specified pad.
// This zone index is different from the zone ID,
// as it is unique for each pad dimensions.
// It is composed in this way:
//   sectorID*100 + zoneID*10 + specific index 
// Where sectorID = 0,1 for front/back sector.
// Specific index is present only for zones containing special motifs.
// ---

  if (!pad.IsValid()) {
    if (warning) Warning("Zone(AliMpPad)", "Invalid pad");
    return 0;
  }  

  AliMpIntPair scale = GetScale(pad.GetIndices());  
  const AliMpTransformer* kTransformer = GetTransformer(scale);

  AliMpPad scaledPad = kTransformer->ITransform(pad);
  
  AliMpSectorSegmentation* segmentation = GetSectorSegmentation(scale);
  Int_t zoneID = segmentation->Zone(scaledPad, warning);
  
  // Distinguish zones from front/back sector
  // For back sector - add 10
  if (segmentation == fBackSectorSegmentation) zoneID += 100;

  return zoneID;
}  

//_____________________________________________________________________________
TVector2 
AliMpPlaneSegmentation::PadDimensions(Int_t zone, Bool_t warning) const
{
// Returns the pad dimensions for the zone with the specified zone index.
// ---

  if (zone < 100)
    return fFrontSectorSegmentation->PadDimensions(zone, warning);
  else  
    return fBackSectorSegmentation->PadDimensions(zone - 100, warning);
}  

//_____________________________________________________________________________
Bool_t AliMpPlaneSegmentation::CircleTest(const AliMpIntPair& indices) const
{
// Verifies that all methods for retrieving pads are consistents between them.
// Returns true if the pad with specified indices was found and verified,
// false otherwise.
// ---

  if (!HasPad(indices)) return false;

  // Verify the indice->location->position->indice way
  AliMpIntPair location = PadByIndices(indices).GetLocation();
  TVector2 position = PadByLocation(location).Position();
  AliMpIntPair retIndices = PadByPosition(position).GetIndices();
    
  if (retIndices != indices) {
    cout << "Pad " << indices << " lead to inconsistency" << endl;
    cout << "in indice->location->position->indice way..." << endl;
    cout << "starting from " << indices << "-->" << location << "-->" 
         << '(' << position.X() << ',' << position.Y() << ')'
         << " and retIndices: " << retIndices << endl;
  }
    
    
  // Verify the indice->position->location->indice way    
  position = PadByIndices(indices).Position();
  location = PadByPosition(position).GetLocation();
  retIndices = PadByLocation(location).GetIndices();

  if (retIndices != indices) {
    cout << "Pad " << indices << " lead to inconsistency" << endl;
    cout << "in indice->position->location->indice way..." <<endl;
    cout << "starting from " << indices 
         << " and retIndices: " << retIndices << endl;
  }
  
  return true;

}

//_____________________________________________________________________________
Int_t AliMpPlaneSegmentation::GetNofTransformers() const
{
// Returns number of transformers.
// ---

#ifdef WITH_STL
  return fTransformers.size();
#endif

#ifdef WITH_ROOT
  return fTransformers.GetEntriesFast();
#endif
}  


//_____________________________________________________________________________
AliMpTransformer* AliMpPlaneSegmentation::GetTransformer(Int_t i) const
{
// Returns i-th transformer.
// ---
 
#ifdef WITH_STL
  return  fTransformers[i];
#endif

#ifdef WITH_ROOT
  return  (AliMpTransformer*)fTransformers[i];
#endif
}     


