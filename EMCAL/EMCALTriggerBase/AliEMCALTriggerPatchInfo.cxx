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
#include "AliEMCALTriggerPatchInfo.h"
#include "AliLog.h"
#include "AliEMCALGeometry.h"
#include "TArrayI.h"
#include "TMath.h"

/// \cond CLASSIMP
ClassImp(AliEMCALTriggerPatchInfo)
/// \endcond

AliEMCALTriggerPatchInfo::AliEMCALTriggerPatchInfo() :
  TObject(),
  fCenterGeo(),
  fCenterMass(),
  fEdge1(),
  fEdge2(),
  fADCAmp(0),
  fADCOfflineAmp(0),
  fTriggerBits(0),
  fOffSet(0),            // To be set explicitly by the trigger maker in order to avoid hard coding
  fRow0(-1),
  fCol0(-1),
  fPatchSize(0),
  fDetectorType(kEMCALdet),
  fTriggerBitConfig()
{
  fEdgeCell[0] = -1;
  fEdgeCell[1] = -1;
}

AliEMCALTriggerPatchInfo::AliEMCALTriggerPatchInfo(const AliEMCALTriggerPatchInfo &p) :
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

AliEMCALTriggerPatchInfo::~AliEMCALTriggerPatchInfo()
{
}

AliEMCALTriggerPatchInfo &AliEMCALTriggerPatchInfo::operator=(const AliEMCALTriggerPatchInfo &p)
{
  if (this != &p) {
    new (this) AliEMCALTriggerPatchInfo(p);
  }

  return *this;
}

void AliEMCALTriggerPatchInfo::Reset()
{
  new (this) AliEMCALTriggerPatchInfo();
}

void AliEMCALTriggerPatchInfo::Initialize(UChar_t col0, UChar_t row0, UChar_t size, UInt_t adc, UInt_t offlineAdc, Double_t patchE, UInt_t bitmask, const TVector3& vertex, const AliEMCALGeometry* geom)
{
  fCol0 = col0;
  fRow0 = row0;
  fPatchSize = size;
  fADCAmp = adc;
  fADCOfflineAmp = offlineAdc;
  fTriggerBits = bitmask;
  SetEdgeCell(col0*2, row0*2);

  if (geom) {
    Int_t absId=-1;
    geom->GetAbsFastORIndexFromPositionInEMCAL(fCol0, fRow0, absId);
    Int_t iSM = -1, iEta = -1, iPhi = -1;
    geom->GetPositionInSMFromAbsFastORIndex(absId, iSM, iEta, iPhi);
    if (geom->IsDCALSM(iSM)) {
      SetDetectorType(kDCALPHOSdet);
    }
    else {
      SetDetectorType(kEMCALdet);
    }
  }
  RecalculateKinematics(patchE, vertex, geom);
}

void AliEMCALTriggerPatchInfo::Initialize(UChar_t col0, UChar_t row0, UChar_t size, UInt_t adc, UInt_t offlineAdc, UInt_t bitmask, const AliEMCALGeometry* geom)
{
  fCol0 = col0;
  fRow0 = row0;
  fPatchSize = size;
  fADCAmp = adc;
  fADCOfflineAmp = offlineAdc;
  fTriggerBits = bitmask;
  SetEdgeCell(col0*2, row0*2);

  if (geom) {
    Int_t absId=-1;
    geom->GetAbsFastORIndexFromPositionInEMCAL(fCol0, fRow0, absId);
    Int_t iSM = -1, iEta = -1, iPhi = -1;
    geom->GetPositionInSMFromAbsFastORIndex(absId, iSM, iEta, iPhi);
    if (geom->IsDCALSM(iSM)) {
      SetDetectorType(kDCALPHOSdet);
    }
    else {
      SetDetectorType(kEMCALdet);
    }
  }
}

AliEMCALTriggerPatchInfo* AliEMCALTriggerPatchInfo::CreateAndInitialize(UChar_t col0, UChar_t row0, UChar_t size, UInt_t adc, UInt_t offlineAdc, Double_t patchE, UInt_t bitmask, const TVector3& vertex, const AliEMCALGeometry* geom)
{
  AliEMCALTriggerPatchInfo* patch = new AliEMCALTriggerPatchInfo;
  patch->Initialize(col0, row0, size, adc, offlineAdc, patchE, bitmask, vertex, geom);
  return patch;
}

AliEMCALTriggerPatchInfo* AliEMCALTriggerPatchInfo::CreateAndInitialize(UChar_t col0, UChar_t row0, UChar_t size, UInt_t adc, UInt_t offlineAdc, UInt_t bitmask, const AliEMCALGeometry* geom)
{
  AliEMCALTriggerPatchInfo* patch = new AliEMCALTriggerPatchInfo;
  patch->Initialize(col0, row0, size, adc, offlineAdc, bitmask, geom);
  return patch;
}

void AliEMCALTriggerPatchInfo::RecalculateKinematics(Double_t patchE, const TVector3& vertex, const AliEMCALGeometry* geom)
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
  if(cellAbsId[0] < 0 || cellAbsId[1] < 0 || cellAbsId[2] < 0 || cellAbsId[3] < 0){
    AliWarning(Form("Invalid cell ID [%d|%d|%d|%d]\n", cellAbsId[0], cellAbsId[1], cellAbsId[2], cellAbsId[3]));
  }

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

void AliEMCALTriggerPatchInfo::GetCellIndices( AliEMCALGeometry *geom, TArrayI *cells ){

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

void AliEMCALTriggerPatchInfo::SetLorentzVector( TLorentzVector &lv, const TVector3 &v, Double_t e ){
  // sets the vector
  Double_t r = TMath::Sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2] ) ; 
  
  lv.SetPxPyPzE( e*v[0]/r,  e*v[1]/r,  e*v[2]/r,  e) ;   
}

Double_t AliEMCALTriggerPatchInfo::GetPhiTransform(Double_t phiin) const{
  if(phiin < 0) return phiin + TMath::TwoPi();
  return phiin;
}

Double_t AliEMCALTriggerPatchInfo::GetET(Double_t energy) const {
  TLorentzVector en(fCenterGeo.Px(), fCenterGeo.Py(), fCenterGeo.Pz(), energy);
  return en.Et();
}
