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
#include <AliEMCALTriggerPatchInfo.h>
#include "AliLog.h"
#include "AliEMCALGeometry.h"
#include "TArrayI.h"

const Double_t AliEMCALTriggerPatchInfo::fgkEMCL1ADCtoGeV = 0.07874;

ClassImp(AliEMCALTriggerPatchInfo)

AliEMCALTriggerPatchInfo::AliEMCALTriggerPatchInfo() :
  TObject(),
  fCenterGeo(),
  fCenterMass(),
  fEdge1(),
  fEdge2(),
  fADCAmp(0),
  fADCOfflineAmp(0),
  fTriggerBits(0),
  fOffSet(0),            // To be set explictly by the trigger maker in order to avoid hard coding
  fRow0(0),
  fCol0(0),
  fPatchSize(0),
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

