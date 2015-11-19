/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved. *
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
#include "AliEmcalTriggerPatchInfoAPV1.h"
#include "AliLog.h"
#include "AliEMCALGeometry.h"
#include "TArrayI.h"

/// \cond CLASSIMP
ClassImp(AliEmcalTriggerPatchInfoAPV1)
/// \endcond

/**
 * Default constructor
 */
AliEmcalTriggerPatchInfoAPV1::AliEmcalTriggerPatchInfoAPV1() :
  TObject(),
  fCenterGeo(),
  fCenterMass(),
  fEdge1(),
  fEdge2(),
  fADCAmp(0),
  fADCOfflineAmp(0),
  fTriggerBits(0),
  fOffSet(0),            // To be set explictly by the trigger maker in order to avoid hard coding
  fRow0(-1),
  fCol0(-1),
  fPatchSize(0),
  fDetectorType(kEMCALdet),
  fTriggerBitConfig()
{
  fEdgeCell[0] = -1;
  fEdgeCell[1] = -1;
}

/**
 * Copy constructor
 *
 * @param p Reference for the copy
 */
AliEmcalTriggerPatchInfoAPV1::AliEmcalTriggerPatchInfoAPV1(const AliEmcalTriggerPatchInfoAPV1 &p) :
  TObject(p),
  fCenterGeo(p.fCenterGeo),
  fCenterMass(p.fCenterMass),
  fEdge1(p.fEdge1),
  fEdge2(p.fEdge2),
  fADCAmp(p.fADCAmp),
  fADCOfflineAmp(p.fADCOfflineAmp),
  fTriggerBits(p.fTriggerBits),
  fOffSet(p.fOffSet),
  fRow0(p.fRow0),
  fCol0(p.fCol0),
  fPatchSize(p.fPatchSize),
  fDetectorType(p.fDetectorType),
  fTriggerBitConfig(p.fTriggerBitConfig)
{
  // .
  fEdgeCell[0] = p.fEdgeCell[0];
  fEdgeCell[1] = p.fEdgeCell[1];
}

/**
 * Destructor
 */
AliEmcalTriggerPatchInfoAPV1::~AliEmcalTriggerPatchInfoAPV1()
{
}

/**
 * Assignment operator
 *
 * @param p Reference for assignment
 * @return This object after assignment
 */
AliEmcalTriggerPatchInfoAPV1 &AliEmcalTriggerPatchInfoAPV1::operator=(const AliEmcalTriggerPatchInfoAPV1 &p)
{
  if (this != &p) {
    fCenterGeo = p.fCenterGeo;
    fCenterMass = p.fCenterMass;
    fEdge1 = p.fEdge1;
    fEdge2 = p.fEdge2;
    fADCAmp = p.fADCAmp;
    fADCOfflineAmp = p.fADCOfflineAmp;
    fTriggerBits = p.fTriggerBits;
    fEdgeCell[0] = p.fEdgeCell[0];
    fEdgeCell[1] = p.fEdgeCell[1];
  }

  return *this;
}
void AliEmcalTriggerPatchInfoAPV1::Initialize(UChar_t col0, UChar_t row0, UChar_t size, UInt_t adc, UInt_t offlineAdc, Double_t patchE, UInt_t bitmask, const TVector3& vertex, const AliEMCALGeometry* geom)
{
  fCol0 = col0;
  fRow0 = row0;
  fPatchSize = size;
  fADCAmp = adc;
  fADCOfflineAmp = offlineAdc;
  fTriggerBits = bitmask;
  SetEdgeCell(col0*2, row0*2);
  RecalculateKinematics(patchE, vertex, geom);
}

AliEmcalTriggerPatchInfoAPV1* AliEmcalTriggerPatchInfoAPV1::CreateAndInitialize(UChar_t col0, UChar_t row0, UChar_t size, UInt_t adc, UInt_t offlineAdc, Double_t patchE, UInt_t bitmask, const TVector3& vertex, const AliEMCALGeometry* geom)
{
  AliEmcalTriggerPatchInfoAPV1* patch = new AliEmcalTriggerPatchInfoAPV1;
  patch->Initialize(col0, row0, size, adc, offlineAdc, patchE, bitmask, vertex, geom);
  return patch;
}

void AliEmcalTriggerPatchInfoAPV1::RecalculateKinematics(Double_t patchE, const TVector3& vertex, const AliEMCALGeometry* geom)
{
  if (!geom) {
    AliError("EMCal geometry pointer not set! Unable to recalculate the trigger patch kinematics!");
    return;
  }

  // get the absolute trigger ID
  Int_t absId=-1;
  geom->GetAbsFastORIndexFromPositionInEMCAL(fCol0, fRow0, absId);
  // convert to the 4 absId of the cells composing the trigger channel
  Int_t cellAbsId[4]={-1,-1,-1,-1};
  geom->GetCellIndexFromFastORIndex(absId, cellAbsId);

  // get low left edge (eta max, phi min)
  TVector3 edge1;
  geom->GetGlobal(cellAbsId[0], edge1);
  Int_t colEdge1 = fCol0, rowEdge1 = fRow0, absIdEdge1 = absId, cellIdEdge1 = cellAbsId[0]; // Used in warning for invalid patch position

  // get up right edge (eta min, phi max)
  // get the absolute trigger ID
  Int_t posOffset = fPatchSize - 1;

  geom->GetAbsFastORIndexFromPositionInEMCAL(fCol0+posOffset, fRow0+posOffset, absId);
  geom->GetCellIndexFromFastORIndex(absId, cellAbsId);
  TVector3 edge2;
  geom->GetGlobal(cellAbsId[3], edge2);
  Int_t colEdge2 = fCol0+posOffset, rowEdge2 = fRow0+posOffset, absIdEdge2 = absId, cellIdEdge2 = cellAbsId[3]; // Used in warning for invalid patch position

  // get the geometrical center as an average of two diagonally
  // adjacent patches in the center
  // picking two diagonally closest cells from the patches
  posOffset = fPatchSize / 2 - 1;

  geom->GetAbsFastORIndexFromPositionInEMCAL(fCol0+posOffset, fRow0+posOffset, absId);
  geom->GetCellIndexFromFastORIndex(absId, cellAbsId);
  TVector3 center1;
  geom->GetGlobal(cellAbsId[3], center1);

  posOffset = fPatchSize / 2;

  geom->GetAbsFastORIndexFromPositionInEMCAL(fCol0+posOffset, fRow0+posOffset, absId);
  geom->GetCellIndexFromFastORIndex(absId, cellAbsId);
  TVector3 center2;
  geom->GetGlobal(cellAbsId[0], center2);

  TVector3 centerGeo(center1);
  centerGeo += center2;
  centerGeo *= 0.5;

  // relate all to primary vertex
  TVector3 edge1tmp = edge1, edge2tmp = edge2; // Used in warning for invalid patch position
  centerGeo -= vertex;
  edge1 -= vertex;
  edge2 -= vertex;
  // Check for invalid patch positions
  if(!(edge1[0] || edge1[1] || edge1[2])){
    AliWarning(Form("Inconsistency in patch position for edge1: [%f|%f|%f]", edge1[0], edge1[1], edge1[2]));
    AliWarning("Original vectors:");
    AliWarning(Form("edge1: [%f|%f|%f]", edge1tmp[0], edge1tmp[1], edge1tmp[2]));
    AliWarning(Form("vertex: [%f|%f|%f]", vertex[0], vertex[1], vertex[2]));
    AliWarning(Form("Col: %d, Row: %d, FABSID: %d, Cell: %d", colEdge1, rowEdge1, absIdEdge1, cellIdEdge1));
  }
  if(!(edge2[0] || edge2[1] || edge2[2])){
    AliWarning(Form("Inconsistency in patch position for edge2: [%f|%f|%f]", edge2[0], edge2[1], edge2[2]));
    AliWarning("Original vectors:");
    AliWarning(Form("edge2: [%f|%f|%f]", edge2tmp[0], edge2tmp[1], edge2tmp[2]));
    AliWarning(Form("vertex: [%f|%f|%f]", vertex[0], vertex[1], vertex[2]));
    AliWarning(Form("Col: %d, Row: %d, FABSID: %d, Cell: %d", colEdge2, rowEdge2, absIdEdge2, cellIdEdge2));
  }

  fCenterMass.SetPxPyPzE(0,0,0,0);

  SetCenterGeo(centerGeo, patchE);
  SetEdge1(edge1, patchE);
  SetEdge2(edge2, patchE);
}


/**
 * Return cell indices of the given patch in the cell array
 * @param geom EMCAL Geometry used in the run where the trigger patch was created from
 * @param cells Output array of cell indices corresponding to the given trigger patch
 */
void AliEmcalTriggerPatchInfoAPV1::GetCellIndices( AliEMCALGeometry *geom, TArrayI *cells ){

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


/**
 * Define Lorentz vector of the given trigger patch
 * @param lv Lorentz vector to be defined
 * @param v Patch vector position
 * @param e Patch energy
 */
void AliEmcalTriggerPatchInfoAPV1::SetLorentzVector( TLorentzVector &lv, const TVector3 &v, Double_t e ){
  // sets the vector
  Double_t r = TMath::Sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2] ) ; 
  
  lv.SetPxPyPzE( e*v[0]/r,  e*v[1]/r,  e*v[2]/r,  e) ;   
}

