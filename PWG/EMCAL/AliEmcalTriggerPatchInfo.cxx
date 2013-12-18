// $Id$
//
// Emcal particle trigger class, which can contain either
//
// Author: J.Kral

#include "AliEmcalTriggerPatchInfo.h"
#include "AliLog.h"
#include "AliEMCALGeometry.h"
#include "TArrayI.h"

//_________________________________________________________________________________________________
AliEmcalTriggerPatchInfo::AliEmcalTriggerPatchInfo() :
  TObject(),
  fCenterGeo(),
  fCenterMass(),
  fEdge1(),
  fEdge2(),
  fADCAmp(0),
  fTriggerBits(0),
  fOffSet(kTriggerTypeEnd)
{
  // Default constructor.
  fEdgeCell[0] = -1;
  fEdgeCell[1] = -1;
}

  
//_________________________________________________________________________________________________
AliEmcalTriggerPatchInfo::AliEmcalTriggerPatchInfo(const AliEmcalTriggerPatchInfo &p) :
  TObject(p),
  fCenterGeo(p.fCenterGeo),
  fCenterMass(p.fCenterMass),
  fEdge1(p.fEdge1),
  fEdge2(p.fEdge2),
  fADCAmp(p.fADCAmp),
  fTriggerBits(p.fTriggerBits),
  fOffSet(p.fOffSet)
{
  // Copy constructor.
  fEdgeCell[0] = p.fEdgeCell[0];
  fEdgeCell[1] = p.fEdgeCell[1];
}

//_________________________________________________________________________________________________
AliEmcalTriggerPatchInfo::~AliEmcalTriggerPatchInfo()
{
  // Destructor.
}

//_________________________________________________________________________________________________
AliEmcalTriggerPatchInfo &AliEmcalTriggerPatchInfo::operator=(const AliEmcalTriggerPatchInfo &p)
{
  // Assignment operator.

  if (this != &p) {
    fCenterGeo = p.fCenterGeo;
    fCenterMass = p.fCenterMass;
    fEdge1 = p.fEdge1;
    fEdge2 = p.fEdge2;
    fADCAmp = p.fADCAmp;
    fTriggerBits = p.fTriggerBits;
    fEdgeCell[0] = p.fEdgeCell[0];
    fEdgeCell[1] = p.fEdgeCell[1];
  }

  return *this;
}

//_________________________________________________________________________________________________
void AliEmcalTriggerPatchInfo::GetCellIndices( AliEMCALGeometry *geom, TArrayI *cells ){

	// return cell indices of the given patch in hte cell array
	Int_t globCol, globRow, i, j, k, absId, cellAbsId[4];;

	cells->Set( 1024 );
	
	// get corner, convert from cells to trigger channels
	globCol = GetEdgeCellX() / 2;
	globRow = GetEdgeCellY() / 2;
	
	// get the absolute trigger ID
	geom->GetAbsFastORIndexFromPositionInEMCAL( globCol, globRow, absId );
	// convert to the 4 absId of the cells composing the trigger channel
	geom->GetCellIndexFromFastORIndex( absId, cellAbsId );
	
	// sum the available energy in the 32/32 window of cells
	// step over trigger channels and get all the corresponding cells
	for( i = 0; i < 16; i++ ){
		for( j = 0; j < 16; j++ ){
			// get the 4 cells composing the trigger channel
			geom->GetAbsFastORIndexFromPositionInEMCAL( globCol+i, globRow+j, absId );
			geom->GetCellIndexFromFastORIndex( absId, cellAbsId );
			// add amplitudes and find patch edges
			for( k = 0; k < 4; k++ ){
				cells->SetAt( cellAbsId[k], i*16*4+j*4+k );
			}
		}
	} // 32x32 cell window

	
}


//_________________________________________________________________________________________________
void AliEmcalTriggerPatchInfo::SetLorentzVector( TLorentzVector &lv, TVector3 &v, Double_t e ){
  // sets the vector
  Double_t r = TMath::Sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2] ) ; 
  
  lv.SetPxPyPzE( e*v[0]/r,  e*v[1]/r,  e*v[2]/r,  e) ;   
}

