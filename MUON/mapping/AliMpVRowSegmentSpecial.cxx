// $Id$
// Category: sector
//
// Class AliMpVRowSegmentSpecial
// ----------------------------
// Class describing a special row segment composed of the 
// pad rows.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <TError.h>
#include <Riostream.h>

#include "AliMpVRowSegmentSpecial.h"
#include "AliMpRow.h"
#include "AliMpPadRow.h"
#include "AliMpVPadRowSegment.h"
#include "AliMpMotif.h"
#include "AliMpMotifType.h"
#include "AliMpMotifMap.h"
#include "AliMpMotifPosition.h"
#include "AliMpConstants.h"

ClassImp(AliMpVRowSegmentSpecial)

#ifdef WITH_ROOT
const Int_t  AliMpVRowSegmentSpecial::fgkMaxNofMotifPositionIds = 20;
#endif    

//______________________________________________________________________________
AliMpVRowSegmentSpecial::AliMpVRowSegmentSpecial(AliMpRow* row, Double_t offsetX)
  : AliMpVRowSegment(),
    fRow(row),
    fOffsetX(offsetX),
    fPadRows(),
    fMotifs(),
    fMotifPositionIds()
#ifdef WITH_ROOT
    ,fNofMotifPositionIds(0)
#endif    
{
// 
}

//______________________________________________________________________________
AliMpVRowSegmentSpecial::AliMpVRowSegmentSpecial() 
  : AliMpVRowSegment(),
    fRow(0),
    fOffsetX(0.),
    fPadRows(),
    fMotifs(),
    fMotifPositionIds()
#ifdef WITH_ROOT
    ,fNofMotifPositionIds(0)
#endif    
{
//
#ifdef WITH_ROOT
   fMotifPositionIds.Set(fgkMaxNofMotifPositionIds);
#endif    
}

//_____________________________________________________________________________
AliMpVRowSegmentSpecial::AliMpVRowSegmentSpecial(
                                  const AliMpVRowSegmentSpecial& right) 
  : AliMpVRowSegment(right) {
// 
  Fatal("AliMpVRowSegmentSpecial", "Copy constructor not provided.");
}

//______________________________________________________________________________
AliMpVRowSegmentSpecial::~AliMpVRowSegmentSpecial() 
{
//  
  for (Int_t i=0; i<GetNofPadRows(); i++)
    delete fPadRows[i];
}

//
// operators
//

//_____________________________________________________________________________
AliMpVRowSegmentSpecial& 
AliMpVRowSegmentSpecial::operator=(const AliMpVRowSegmentSpecial& right)
{
  // check assignement to self
  if (this == &right) return *this;

  Fatal("operator =", "Assignement operator not provided.");
    
  return *this;  
}    

//
// protected methods  
//

//______________________________________________________________________________
AliMpPadRow*  AliMpVRowSegmentSpecial::FindPadRow(Double_t y) const
{
// Finds the pad row in the given y coordinate.
// ---

  Double_t lowBorder =  fRow->LowBorderY();
  Double_t highBorder = fRow->LowBorderY();
  
  for (Int_t i=0; i<GetNofPadRows(); i++) {    

    AliMpPadRow* padRow = GetPadRow(i);
    highBorder += 2.*padRow->HalfSizeY();

    if ( y >= lowBorder &&  y <= highBorder)
      return padRow;

    lowBorder = highBorder;
  }
  
  return 0;   	
}

//______________________________________________________________________________
AliMpVPadRowSegment*  
AliMpVRowSegmentSpecial::FindPadRowSegment(Int_t motifPositionId) const
{
// Find the most down pad row segment with this motifPositionId.
// ---

  for (Int_t i=0; i<GetNofPadRows(); i++) {
    AliMpPadRow* padRow = GetPadRow(i);    

    for (Int_t j=0; j<padRow->GetNofPadRowSegments(); j++) { 
      AliMpVPadRowSegment* padRowSegment = padRow->GetPadRowSegment(j);

      if (padRowSegment->GetMotifPositionId() == motifPositionId) 
        return padRowSegment;
    }
  }
  return 0;   	
}

//______________________________________________________________________________
AliMpIntPair 
AliMpVRowSegmentSpecial::FindRelativeLowIndicesOf(Int_t motifPositionId) const 
{ 
// Returns the lowest pad indices where the motif of the given position ID
// exists in this segment.
// ---

  AliMpIntPair ans(0,1000);
  AliMpIntPair ans0 = ans;
  Int_t maxNofPadsX=0;
  
  for (Int_t i=0; i<GetNofPadRows(); i++) {
    AliMpPadRow* padRow = GetPadRow(i);

    Int_t nofPadsX=0;
    for (Int_t j=0; j<padRow->GetNofPadRowSegments(); j++) {
      AliMpVPadRowSegment* padRowSegment = padRow->GetPadRowSegment(j);
      nofPadsX += padRowSegment->GetNofPads();
      if (padRowSegment->GetMotifPositionId() == motifPositionId) {
         if (ans.GetFirst() < nofPadsX) ans.SetFirst(nofPadsX);
         if (ans.GetSecond()>i) ans.SetSecond(i);
                  // ans.First = max (nof pads of this pos ID)
                  // ans.Second = min of pad row number
      }
    }  
    if (nofPadsX > maxNofPadsX) maxNofPadsX = nofPadsX;
  }    
  if (ans == ans0) return AliMpIntPair::Invalid();
  
  return AliMpIntPair(maxNofPadsX-ans.GetFirst(), ans.GetSecond());
}
 
//______________________________________________________________________________
Int_t  AliMpVRowSegmentSpecial::MaxNofPadsInRow() const 
{ 
// Returns the maximum number of pads in this row segment along the X direction
// ---

  Int_t maxNofPads = 0;    

  for (Int_t i=0; i<GetNofPadRows(); i++){
    Int_t nofPads = GetPadRow(i)->GetNofPads(); 

    // Find maximum
    if (nofPads > maxNofPads) maxNofPads = nofPads;
  }
    
  return maxNofPads;
}

//______________________________________________________________________________
Bool_t AliMpVRowSegmentSpecial::HasMotif(const AliMpVMotif* motif) const
{
// Returns true if the specified motif is already in fMotifs vector,
// returns false otherwise.
// ---

#ifdef WITH_STL
  for (UInt_t i=0; i<fMotifs.size(); i++)
    if (fMotifs[i] == motif) return true;
#endif

#ifdef WITH_ROOT
  for (Int_t i=0; i<fMotifs.GetEntriesFast(); i++)
    if (fMotifs[i] == (const TObject*)motif) return true;
#endif

  return false;	 
}

//______________________________________________________________________________
Int_t AliMpVRowSegmentSpecial::GetNofPadRows() const
{
// Returns number of pad rows.
// ---

#ifdef WITH_STL
  return fPadRows.size();
#endif

#ifdef WITH_ROOT
  return fPadRows.GetEntriesFast();
#endif
}  

//______________________________________________________________________________
AliMpPadRow* AliMpVRowSegmentSpecial::GetPadRow(Int_t i) const
{
// Returns number of pad rows.
// ---

#ifdef WITH_STL
  return fPadRows[i];
#endif

#ifdef WITH_ROOT
  return (AliMpPadRow*)fPadRows[i];
#endif
}  

//
// public methods  
//

//______________________________________________________________________________
void  AliMpVRowSegmentSpecial::AddPadRow(AliMpPadRow* padRow)
{
// Adds a pad row.
// ---

  padRow->SetOffsetX(fOffsetX);
  padRow->SetID(GetNofPadRows());

#ifdef WITH_STL
  fPadRows.push_back(padRow);
#endif

#ifdef WITH_ROOT
  fPadRows.Add(padRow);
#endif
}  

//______________________________________________________________________________
void AliMpVRowSegmentSpecial::UpdateMotifVector()
{
// Add motifs associated with the pad row segments in the specified
// pad row in the fMotifs vector.
// ---

  for (Int_t i=0; i<GetNofPadRows(); i++) {
    AliMpPadRow* padRow = GetPadRow(i);
 
    for (Int_t j=0; j<padRow->GetNofPadRowSegments(); j++) {
      AliMpVMotif* motif = padRow->GetPadRowSegment(j)->GetMotif();            

      if (!HasMotif(motif)) {
#ifdef WITH_STL
        fMotifs.push_back(motif);	 
        fMotifPositionIds.push_back(
          padRow->GetPadRowSegment(j)->GetMotifPositionId());
#endif
#ifdef WITH_ROOT
        fMotifs.Add(motif);
	
	// resize array if needed
	if (fNofMotifPositionIds<fgkMaxNofMotifPositionIds)
	  fMotifPositionIds.Set(fMotifPositionIds.GetSize()+
	                        fgkMaxNofMotifPositionIds);	 
        fMotifPositionIds.AddAt(
          padRow->GetPadRowSegment(j)->GetMotifPositionId(),
	  fNofMotifPositionIds++);
#endif
      }
    }  
  }
}

//______________________________________________________________________________
Double_t  AliMpVRowSegmentSpecial::HalfSizeY() const
{
// Returns the size in y of this row segment.
// ---

  Double_t halfSizeY = 0.;
  for (Int_t i=0; i<GetNofPadRows(); i++) {
    halfSizeY += GetPadRow(i)->HalfSizeY();
  }  
  
  return halfSizeY;
}

//______________________________________________________________________________
AliMpVMotif*  AliMpVRowSegmentSpecial::FindMotif(const TVector2& position) const
{
// Returns the motif of this row; 
// ---

  AliMpPadRow* padRow 
    = FindPadRow(position.Y());
  
  if (!padRow) return 0;

  AliMpVPadRowSegment* padRowSegment 
    = padRow->FindPadRowSegment(position.X());
    
  if (!padRowSegment) return 0;

  return padRowSegment->GetMotif();
}  

//______________________________________________________________________________
Int_t AliMpVRowSegmentSpecial::FindMotifPositionId(const TVector2& position) const
{
// Returns the motif position identified for the given
// geometric position.
// ---

  AliMpPadRow* padRow 
    = FindPadRow(position.Y());
  
  if (!padRow) return 0;

  AliMpVPadRowSegment* padRowSegment 
    = padRow->FindPadRowSegment(position.X());
    
  if (!padRowSegment) return 0;

  return padRowSegment->GetMotifPositionId();
}

//______________________________________________________________________________
Bool_t AliMpVRowSegmentSpecial::HasMotifPosition(Int_t motifPositionId) const
{
// Returns true if the motif specified with the given position identifier
// is in this segment.
// ---

  if (FindPadRowSegment(motifPositionId))
    return true;
  else  
    return false;   	
}

//______________________________________________________________________________
TVector2 AliMpVRowSegmentSpecial::MotifCenter(Int_t motifPositionId) const
{
// Returns the coordinates of the motif specified with
// the given position identifier.
// ---

  // Try to get the motif position from the motif map first
  AliMpMotifPosition* motifPosition
    = GetRow()->GetMotifMap()->FindMotifPosition(motifPositionId);
  if (motifPosition) return motifPosition->Position();

  // Use slow method otherwise
  return MotifCenterSlow(motifPositionId);
}

//______________________________________________________________________________
TVector2 AliMpVRowSegmentSpecial::Dimensions() const
{
// Returns the halflengths in x, y of the row segment rectangular envelope.
// ---

  Double_t x = 0.;		    
  Double_t y = 0.;  
  for (Int_t i=0; i<GetNofPadRows(); i++) {    
    AliMpPadRow* padRow = GetPadRow(i); 
    
    // Add all pad rows y halfsizes   
    y += padRow->HalfSizeY();

    // Find the biggest pad rows x halfsize
    Double_t xx 
      = (padRow->GetPadRowSegment(0)->RightBorderX() -
         padRow->GetPadRowSegment(padRow->GetNofPadRowSegments()-1)->LeftBorderX())/2.;
    if (xx > x) x = xx;		   
  }                  
    
  return TVector2(x, y);   
}

//______________________________________________________________________________
void AliMpVRowSegmentSpecial::SetGlobalIndices()
{
// Sets indices limits.
// ---

  AliMpMotifPosition* firstPos = 0;
  AliMpMotifPosition* lastPos = 0;
	
  for (Int_t i=0;i<GetNofMotifs();i++) {
    AliMpMotifPosition* mPos
      = GetRow()->GetMotifMap()
        ->FindMotifPosition(GetMotifPositionId(i));
	
    if (!firstPos || 
        mPos->GetLowIndicesLimit().GetFirst()
	< firstPos->GetLowIndicesLimit().GetFirst())
       firstPos = mPos;
       
    if (!lastPos || 
        mPos->GetHighIndicesLimit().GetFirst()
	>lastPos->GetHighIndicesLimit().GetFirst())
       lastPos = mPos;
  }

  // Check if the motif positions has the limits set
  if ( !firstPos->HasValidIndices() || !lastPos->HasValidIndices())
    Fatal("SetGlobalIndices", "Indices of motif positions have to be set first.");

  SetLowIndicesLimit(firstPos->GetLowIndicesLimit());
  SetHighIndicesLimit(lastPos->GetHighIndicesLimit());
}  

//______________________________________________________________________________
AliMpRow*  AliMpVRowSegmentSpecial::GetRow() const
{
// Returns the row.which this row segment belongs to.
// ---

  return fRow;
}  

//______________________________________________________________________________
Int_t  AliMpVRowSegmentSpecial::GetNofMotifs() const 
{ 
// Returns the number of different motifs present in this row segment.
// ---

#ifdef WITH_STL
  return fMotifs.size();
#endif
#ifdef WITH_ROOT
  return fMotifs.GetEntriesFast();
#endif
}  

//______________________________________________________________________________
AliMpVMotif* AliMpVRowSegmentSpecial::GetMotif(Int_t i) const  
{
// Returns the i-th motif present in this row segment.
// ---

#ifdef WITH_STL
   return fMotifs[i]; 
#endif
#ifdef WITH_ROOT
   return (AliMpVMotif*)fMotifs[i]; 
#endif
}

//______________________________________________________________________________
Int_t  AliMpVRowSegmentSpecial::GetMotifPositionId(Int_t i) const 
{ 
// Returns the i-th motif position Id present in this row segment.
// ---

   return fMotifPositionIds[i]; 
} 

