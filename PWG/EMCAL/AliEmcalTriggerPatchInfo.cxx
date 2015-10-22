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
#include "AliEmcalTriggerPatchInfo.h"
#include "AliLog.h"
#include "AliEMCALGeometry.h"
#include "TArrayI.h"

/**
 * Default constructor
 */
AliEmcalTriggerPatchInfo::AliEmcalTriggerPatchInfo() :
  TObject(),
  fCenterGeo(),
  fCenterMass(),
  fEdge1(),
  fEdge2(),
  fADCAmp(0),
  fADCOfflineAmp(0),
  fTriggerBits(0),
  fOffSet(0),            // To be set explictly by the trigger maker in order to avoid hard coding
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
AliEmcalTriggerPatchInfo::AliEmcalTriggerPatchInfo(const AliEmcalTriggerPatchInfo &p) :
  TObject(p),
  fCenterGeo(p.fCenterGeo),
  fCenterMass(p.fCenterMass),
  fEdge1(p.fEdge1),
  fEdge2(p.fEdge2),
  fADCAmp(p.fADCAmp),
  fADCOfflineAmp(p.fADCOfflineAmp),
  fTriggerBits(p.fTriggerBits),
  fOffSet(p.fOffSet),
  fTriggerBitConfig(p.fTriggerBitConfig)
{
  // .
  fEdgeCell[0] = p.fEdgeCell[0];
  fEdgeCell[1] = p.fEdgeCell[1];
}

/**
 * Destructor
 */
AliEmcalTriggerPatchInfo::~AliEmcalTriggerPatchInfo()
{
}

/**
 * Assignment operator
 *
 * @param p Reference for assignment
 * @return This object after assignment
 */
AliEmcalTriggerPatchInfo &AliEmcalTriggerPatchInfo::operator=(const AliEmcalTriggerPatchInfo &p)
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

/**
 * Return cell indices of the given patch in the cell array
 * @param geom EMCAL Geometry used in the run where the trigger patch was created from
 * @param cells Output array of cell indices corresponding to the given trigger patch
 */
void AliEmcalTriggerPatchInfo::GetCellIndices( AliEMCALGeometry *geom, TArrayI *cells ){

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
void AliEmcalTriggerPatchInfo::SetLorentzVector( TLorentzVector &lv, const TVector3 &v, Double_t e ){
  // sets the vector
  Double_t r = TMath::Sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2] ) ; 
  
  lv.SetPxPyPzE( e*v[0]/r,  e*v[1]/r,  e*v[2]/r,  e) ;   
}

