// $Id$
// Category: plane
//
// Class AliMpPlane
// ----------------
// Class represents the plane composed of 4 sector positions:
// 
//   I.  FS                             II. |  I.
//  II.  BS inverted in x             _____ | ____
// III.  FS inverted in x, y                |
//  IV.  BS inverted in y              III. |  IV.
//   
// FS - front sector
// BS - back sector    
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <Riostream.h>

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
  
  // Create sector positions
  fSectorPositions.push_back(
    new AliMpSectorPosition(fkFrontSector, q1Position, AliMpIntPair( 1, 1)));
  fSectorPositions.push_back(
    new AliMpSectorPosition(fkBackSector,  q2Position, AliMpIntPair(-1, 1)));
  fSectorPositions.push_back(
    new AliMpSectorPosition(fkFrontSector, q3Position, AliMpIntPair(-1,-1)));
  fSectorPositions.push_back(
    new AliMpSectorPosition(fkBackSector,  q4Position, AliMpIntPair( 1,-1)));
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
// public methods
//

//______________________________________________________________________________
const AliMpSectorPosition* 
AliMpPlane::SectorPosition(const AliMpIntPair& scale) const
{
// Returns the sector position specified by scale.
// ---

  for (UInt_t i=0; i<fSectorPositions.size(); i++) 
    if (fSectorPositions[i]->GetScale() == scale) return GetSectorPosition(i);

  Fatal("SectorPosition", "Wrong scale");
  return 0; 
}

//______________________________________________________________________________
Int_t AliMpPlane::GetNofSectorPositions() const
{
// Returns number of sector positions.
// ---

  return fSectorPositions.size();
}  


//______________________________________________________________________________
AliMpSectorPosition* AliMpPlane::GetSectorPosition(Int_t i) const
{
// Returns i-th sector position.
// ---
 
  return  fSectorPositions[i];
}     


