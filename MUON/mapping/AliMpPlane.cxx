// $Id$
// Category: plane
//
// Class AliMpPlane
// ----------------
// Class represents the plane composed of 4 sector positions:
//   I.  FS                             II. |  I.
//  II.  BS inverted in x             _____ | ____
// III.  FS inverted in x, y                |
//  IV.  BS inverted in y              III. |  IV.
//   
// FS - front sector
// BS - back sector    
//
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <Riostream.h>
#include <TError.h>

#include "AliMpPlane.h"
#include "AliMpReader.h"
#include "AliMpSector.h"
#include "AliMpSectorPosition.h"
#include "AliMpIntPair.h"

ClassImp(AliMpPlane)

//
// static methods
//

//______________________________________________________________________________
AliMpPlane* 
AliMpPlane::Create(AliMpStationType station, AliMpPlaneType type, 
                   const TVector2& q1Position, const TVector2& q2Position,
	           const TVector2& q3Position, const TVector2& q4Position)
{
// Factory method for creating planes.
// ---

  // Build sectors
  AliMpReader bReader(station, kBendingPlane);
  // bReader.SetVerboseLevel(1);
  AliMpSector* bSector = bReader.BuildSector();
  cout << "bending sector is built" << endl;

  AliMpReader nbReader(station, kNonBendingPlane);
  // nbReader.SetVerboseLevel(1);
  AliMpSector* nbSector = nbReader.BuildSector();
  cout << "non-bending sector is built" << endl;

  if (type == kBendingPlane)
    return new AliMpPlane(bSector, nbSector,
                      q1Position, q2Position, q3Position, q4Position );
  else  
    return new AliMpPlane(nbSector, bSector,
                      q1Position, q2Position, q3Position, q4Position );
}

//______________________________________________________________________________
AliMpPlane* AliMpPlane::Create(AliMpStationType station, AliMpPlaneType type) 
{
// Factory method for creating planes with 
// not shifted qudrants.
// ---

  return Create(station, type, TVector2(), TVector2(), TVector2(), TVector2());
}

//
// constructors, destructors
//

//______________________________________________________________________________
AliMpPlane::AliMpPlane(AliMpSector* frontSector, AliMpSector* backSector,
                       const TVector2& q1Position, const TVector2& q2Position,
	               const TVector2& q3Position, const TVector2& q4Position)
  : TObject(),
    fkFrontSector(frontSector),
    fkBackSector(backSector),
    fSectorPositions()
{
//
  
#ifdef WITH_STL
  // Create sector positions
  fSectorPositions.push_back(
    new AliMpSectorPosition(fkFrontSector, q1Position, AliMpIntPair( 1, 1)));
  fSectorPositions.push_back(
    new AliMpSectorPosition(fkBackSector,  q2Position, AliMpIntPair(-1, 1)));
  fSectorPositions.push_back(
    new AliMpSectorPosition(fkFrontSector, q3Position, AliMpIntPair(-1,-1)));
  fSectorPositions.push_back(
    new AliMpSectorPosition(fkBackSector,  q4Position, AliMpIntPair( 1,-1)));
#endif

#ifdef WITH_ROOT
  // Create sector positions
  fSectorPositions.Add(
    new AliMpSectorPosition(fkFrontSector, q1Position, AliMpIntPair( 1, 1)));
  fSectorPositions.Add(
    new AliMpSectorPosition(fkBackSector,  q2Position, AliMpIntPair(-1, 1)));
  fSectorPositions.Add(
    new AliMpSectorPosition(fkFrontSector, q3Position, AliMpIntPair(-1,-1)));
  fSectorPositions.Add(
    new AliMpSectorPosition(fkBackSector,  q4Position, AliMpIntPair( 1,-1)));
#endif
}

//_____________________________________________________________________________
AliMpPlane::AliMpPlane(const AliMpPlane& right) 
  : TObject(right) {
// 
  Fatal("AliMpPlane", "Copy constructor not provided.");
}

//______________________________________________________________________________
AliMpPlane::AliMpPlane() 
  : TObject(),
    fkFrontSector(0),
    fkBackSector(0),
    fSectorPositions()
{
//
}

//______________________________________________________________________________
AliMpPlane::~AliMpPlane() {
// 

  delete fkFrontSector; 
  delete fkBackSector; 

  for (Int_t i=0; i<GetNofSectorPositions(); i++) 
    delete GetSectorPosition(i);    
}

//
// operators
//

//_____________________________________________________________________________
AliMpPlane& AliMpPlane::operator=(const AliMpPlane& right)
{
  // check assignement to self
  if (this == &right) return *this;

  Fatal("operator =", "Assignement operator not provided.");
    
  return *this;  
}    

//
// public methods
//

//______________________________________________________________________________
const AliMpSectorPosition* 
AliMpPlane::SectorPosition(const AliMpIntPair& scale) const
{
// Returns the sector position specified by scale.
// ---

  for (Int_t i=0; i<GetNofSectorPositions(); i++) 
    if (GetSectorPosition(i)->GetScale() == scale) return GetSectorPosition(i);

  Fatal("SectorPosition", "Wrong scale");
  return 0; 
}

//______________________________________________________________________________
Int_t AliMpPlane::GetNofSectorPositions() const
{
// Returns number of sector positions.
// ---

#ifdef WITH_STL
  return fSectorPositions.size();
#endif

#ifdef WITH_ROOT
  return fSectorPositions.GetEntriesFast();
#endif
}  


//______________________________________________________________________________
AliMpSectorPosition* AliMpPlane::GetSectorPosition(Int_t i) const
{
// Returns i-th sector position.
// ---
 
#ifdef WITH_STL
  return  fSectorPositions[i];
#endif

#ifdef WITH_ROOT
  return  (AliMpSectorPosition*)fSectorPositions[i];
#endif
}     


