/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id: AliMUONAlignment.cxx 51000 2011-08-08 17:58:17Z ivana $ */

//-----------------------------------------------------------------------------
/// \class AliMUONAlignment
/// Alignment class for the ALICE DiMuon spectrometer
///
/// MUON specific alignment class which interface to AliMillepede.
/// For each track ProcessTrack calculates the local and global derivatives
/// at each cluster and fill the corresponding local equations. Provide methods
/// for fixing or constraining detection elements for best results.
///
/// \author Bruce Becker, Javier Castillo
//-----------------------------------------------------------------------------

#include "AliMUONAlignment.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONVCluster.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryModuleTransformer.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMUONGeometryBuilder.h"
#include "AliMillePede2.h"

#include "AliMpExMap.h"
#include "AliMpExMapIterator.h"

#include "AliAlignObjMatrix.h"
#include "AliLog.h"

#include <TMath.h>
#include <TMatrixDSym.h>
#include <TClonesArray.h>
#include <TGraphErrors.h>

/// \cond CLASSIMP
ClassImp(AliMUONAlignment)
/// \endcond

//_____________________________________________________________________
// static variables
const Int_t AliMUONAlignment::fgNDetElemCh[AliMUONAlignment::fgNCh] = { 4, 4, 4, 4, 18, 18, 26, 26, 26, 26 };
const Int_t AliMUONAlignment::fgSNDetElemCh[AliMUONAlignment::fgNCh+1] = { 0, 4, 8, 12, 16, 34, 52, 78, 104, 130, 156 };

//_____________________________________________________________________
/// self initialized array, used for adding constraints
class Array
{

  public:

  /// contructor
  Array( void )
  {
    for( Int_t i=0; i < AliMUONAlignment::fNGlobal; ++i )
    { values[i] = 0; }
  }

  /// array
  Double_t values[AliMUONAlignment::fNGlobal];

  private:

  /// Not implemented
  Array(const Array& );

  /// Not implemented
  Array&  operator = (const Array& );

};

//_____________________________________________________________________
AliMUONAlignment::AliMUONAlignment()
  : TObject(),
    fInitialized( kFALSE ),
    fRunNumber( 0 ),
    fBFieldOn( kTRUE ),
    fStartFac( 256 ),
    fResCutInitial( 100 ),
    fResCut( 100 ),
    fMillepede( 0L ),
    fCluster( 0L ),
    fNStdDev( 3 ),
    fDetElemNumber( 0 ),
    fTrackRecord(),
    fTransform( 0 ),
    fGeoCombiTransInverse()
{
  /// constructor
  fSigma[0] = 1.5e-1;
  fSigma[1] = 1.0e-2;

  // default allowed variations
  fAllowVar[0] = 0.5;  // x
  fAllowVar[1] = 0.5;  // y
  fAllowVar[2] = 0.01; // phi_z
  fAllowVar[3] = 5;    // z

  // initialize millepede
  fMillepede = new AliMillePede2();

  // initialize degrees of freedom
  // by default all parameters are free
  for( Int_t iPar = 0; iPar < fNGlobal; ++iPar )
  { fGlobalParameterStatus[iPar] = kFreeParId; }

  // initialize local equations
  for(int i=0; i<fNLocal; ++i )
  { fLocalDerivatives[i] = 0.0; }

  for(int i=0; i<fNGlobal; ++i )
  { fGlobalDerivatives[i] = 0.0; }

}

//_____________________________________________________________________
AliMUONAlignment::~AliMUONAlignment()
{
  /// destructor
}

//_____________________________________________________________________
void AliMUONAlignment::Init( void )
{

  /// initialize
  /**
  initialize millipede
  must be called after necessary detectors have been fixed,
  but before constrains are added and before global parameters initial value are set
  */
  if( fInitialized )
  { AliFatal( "Millepede already initialized" ); }

  // assign proper groupID to free parameters
  Int_t nGlobal = 0;
  for( Int_t iPar = 0; iPar < fNGlobal; ++iPar )
  {

    if( fGlobalParameterStatus[iPar] == kFixedParId )
    {
      // fixed parameters are left unchanged
      continue;

    } else if( fGlobalParameterStatus[iPar] == kFreeParId || fGlobalParameterStatus[iPar] == kGroupBaseId ) {

      // free parameters or first element of group are assigned a new group id
      fGlobalParameterStatus[iPar] = nGlobal++;
      continue;

    } else if( fGlobalParameterStatus[iPar] < kGroupBaseId ) {

      // get detector element id from status, get chamber parameter id
      const Int_t iDeBase( kGroupBaseId - 1 - fGlobalParameterStatus[iPar] );
      const Int_t iParBase = iPar%fgNParCh;

      // check
      if( iDeBase < 0 || iDeBase >= iPar/fgNParCh )
      { AliFatal( Form( "Group for parameter index %i has wrong base detector element: %i", iPar, iDeBase ) ); }

      // assign identical group id to current
      fGlobalParameterStatus[iPar] = fGlobalParameterStatus[iDeBase*fgNParCh + iParBase];
      AliInfo( Form( "Parameter %i grouped to detector %i (%s)", iPar, iDeBase, GetParameterMaskString( 1<<iParBase ).Data() ) );

    } else AliFatal( Form( "Unrecognized parameter status for index %i: %i", iPar, fGlobalParameterStatus[iPar] ) );

  }

  AliInfo( Form( "Free Parameters: %i out of %i", nGlobal, fNGlobal ) );

  // initialize millepede
  fMillepede->InitMille( fNGlobal, fNLocal, fNStdDev, fResCut, fResCutInitial, fGlobalParameterStatus );
  fInitialized = kTRUE;

  // some debug output
  for( Int_t iPar = 0; iPar < fgNParCh; ++iPar )
  {  AliInfo( Form( "fAllowVar[%i]= %f", iPar, fAllowVar[iPar] ) ); }

  // set allowed variations for all parameters
  for( Int_t iDet = 0; iDet < fgNDetElem; ++iDet )
  {
    for( Int_t iPar = 0; iPar < fgNParCh; ++iPar )
    { fMillepede->SetParSigma( iDet*fgNParCh + iPar, fAllowVar[iPar] ); }
  }

  // Set iterations
  if (fStartFac>1) fMillepede->SetIterations(fStartFac);

}

//_____________________________________________________
AliMillePedeRecord* AliMUONAlignment::ProcessTrack( AliMUONTrack* track, Bool_t doAlignment, Double_t weight )
{
  /// process track for alignment minimization
  /**
  returns the alignment records for this track.
  They can be stored in some output for later reprocessing.
  */

  // reset track records
  fTrackRecord.Reset();
  if( fMillepede->GetRecord() ) fMillepede->GetRecord()->Reset();

  // get number of track parameters
  Int_t nTrackParam = track->GetTrackParamAtCluster()->GetEntries();

  Bool_t first( kTRUE );
  for( Int_t iTrackParam = 0; iTrackParam < nTrackParam; ++iTrackParam )
  {

    // get new pointers
    AliMUONTrackParam* trackParam( (AliMUONTrackParam *) track->GetTrackParamAtCluster()->At(iTrackParam) );
    if( !trackParam ) continue;

    AliMUONVCluster* cluster = trackParam->GetClusterPtr();
    if( !cluster ) continue;

    // fill local variables for this position --> one measurement
    FillDetElemData( cluster );
    FillRecPointData( cluster );
    FillTrackParamData( trackParam );

    if( first )
    {

      // for first valid cluster, save track position as "starting" values
      first = kFALSE;

      fTrackPos0[0] = fTrackPos[0];
      fTrackPos0[1] = fTrackPos[1];
      fTrackPos0[2] = fTrackPos[2];
      fTrackSlope0[0] = fTrackSlope[0];
      fTrackSlope0[1] = fTrackSlope[1];

    }

    // 'inverse' (GlobalToLocal) rotation matrix
    const Double_t* r( fGeoCombiTransInverse.GetRotationMatrix() );

    // calculate measurements
    if( fBFieldOn )
    {

      // use residuals (cluster - track) for measurement
      fMeas[0] = r[0]*(fClustPos[0] - fTrackPos[0]) + r[1]*(fClustPos[1] - fTrackPos[1]);
      fMeas[1] = r[3]*(fClustPos[0] - fTrackPos[0]) + r[4]*(fClustPos[1] - fTrackPos[1]);

    } else {

      // use cluster position for measurement
      fMeas[0] = ( r[0]*fClustPos[0] + r[1]*fClustPos[1] );
      fMeas[1] = ( r[3]*fClustPos[0] + r[4]*fClustPos[1] );

    }

    // Set local equations
    LocalEquationX();
    LocalEquationY();

  }

  // copy track record
  fMillepede->SetRecordRun(fRunNumber);
  fMillepede->SetRecordWeight(weight);
  fTrackRecord = *fMillepede->GetRecord();

  // save record data
  if( doAlignment )
  {
    fMillepede->SaveRecordData();
  }

  // return record
  return &fTrackRecord;

}

//______________________________________________________________________________
void AliMUONAlignment::ProcessTrack( AliMillePedeRecord* trackRecord )
{
  /// process track record
  if( !trackRecord ) return;

  // make sure record storage is initialized
  if( !fMillepede->GetRecord() ) fMillepede->InitDataRecStorage();

  // copy content
  *fMillepede->GetRecord() = *trackRecord;

  // save record
  fMillepede->SaveRecordData();

  return;

}

//_____________________________________________________________________
void AliMUONAlignment::FixAll( UInt_t mask )
{
  /// fix parameters matching mask, for all chambers
  AliInfo( Form( "Fixing %s for all detector elements", GetParameterMaskString( mask ).Data() ) );

  // fix all stations
  for( Int_t i = 0; i < fgNDetElem; ++i )
  {
    if( mask & ParX )  FixParameter(i, 0);
    if( mask & ParY )  FixParameter(i, 1);
    if( mask & ParTZ ) FixParameter(i, 2);
    if( mask & ParZ )  FixParameter(i, 3);
  }

}

//_____________________________________________________________________
void AliMUONAlignment::FixChamber( Int_t iCh, UInt_t mask )
{
  /// fix parameters matching mask, for all detector elements in a given chamber, counting from 1

  // check boundaries
  if( iCh < 1 || iCh > 10 )
  { AliFatal( Form( "Invalid chamber index %i", iCh ) ); }

  // get first and last element
  const Int_t iDetElemFirst = fgSNDetElemCh[iCh-1];
  const Int_t iDetElemLast = fgSNDetElemCh[iCh];
  for( Int_t i = iDetElemFirst; i < iDetElemLast; ++i )
  {

    AliInfo( Form( "Fixing %s for detector element %i", GetParameterMaskString(mask).Data(), i ) );

    if( mask & ParX )  FixParameter(i, 0);
    if( mask & ParY )  FixParameter(i, 1);
    if( mask & ParTZ ) FixParameter(i, 2);
    if( mask & ParZ )  FixParameter(i, 3);

  }
}

//_____________________________________________________________________
void AliMUONAlignment::FixDetElem( Int_t iDetElemId, UInt_t mask )
{
  /// fix parameters matching mask, for a given detector element, counting from 0
  const Int_t iDet( GetDetElemNumber( iDetElemId ) );
  if ( mask & ParX )  FixParameter(iDet, 0);
  if ( mask & ParY )  FixParameter(iDet, 1);
  if ( mask & ParTZ ) FixParameter(iDet, 2);
  if ( mask & ParZ )  FixParameter(iDet, 3);

}

//_____________________________________________________________________
void AliMUONAlignment::FixHalfSpectrometer( const Bool_t *lChOnOff, UInt_t sidesMask, UInt_t mask )
{

  /// Fix parameters matching mask for all detectors in selected chambers and selected sides of the spectrometer
  for( Int_t i = 0; i < fgNDetElem; ++i )
  {

    // get chamber matching detector
    const Int_t iCh( GetChamberId(i) );
    if( !lChOnOff[iCh-1] ) continue;

    // get detector element in chamber
    Int_t lDetElemNumber = i-fgSNDetElemCh[iCh-1];

    // skip detector if its side is off
    // stations 1 and 2
    if( iCh>=1 && iCh<=4 )
    {
      if( lDetElemNumber == 0 && !( sidesMask & SideTopRight ) ) continue;
      if( lDetElemNumber == 1 && !( sidesMask & SideTopLeft ) ) continue;
      if( lDetElemNumber == 2 && !( sidesMask & SideBottomLeft ) ) continue;
      if( lDetElemNumber == 3 && !( sidesMask & SideBottomRight ) ) continue;
    }

    // station 3
    if (iCh>=5 && iCh<=6)
    {
      if( lDetElemNumber >= 0 && lDetElemNumber <= 4 && !( sidesMask & SideTopRight ) ) continue;
      if( lDetElemNumber >= 5 && lDetElemNumber <= 10 && !( sidesMask & SideTopLeft ) ) continue;
      if( lDetElemNumber >= 11 && lDetElemNumber <= 13 && !( sidesMask & SideBottomLeft ) ) continue;
      if( lDetElemNumber >= 14 && lDetElemNumber <= 17 && !( sidesMask & SideBottomRight ) ) continue;
    }

    // stations 4 and 5
    if (iCh>=7 && iCh<=10)
    {
      if( lDetElemNumber >= 0 && lDetElemNumber <= 6 && !( sidesMask & SideTopRight ) ) continue;
      if( lDetElemNumber >= 7 && lDetElemNumber <= 13 && !( sidesMask & SideTopLeft ) ) continue;
      if( lDetElemNumber >= 14 && lDetElemNumber <= 19 && !( sidesMask & SideBottomLeft ) ) continue;
      if( lDetElemNumber >= 20 && lDetElemNumber <= 25 && !( sidesMask & SideBottomRight ) ) continue;
    }

    // detector is accepted, fix it
    FixDetElem( i, mask );

  }

}

//______________________________________________________________________
void AliMUONAlignment::FixParameter( Int_t iPar )
{

  /// fix a given parameter, counting from 0
  if( fInitialized )
  { AliFatal( "Millepede already initialized" ); }

  fGlobalParameterStatus[iPar] = kFixedParId;

}


//_____________________________________________________________________
void AliMUONAlignment::ReleaseChamber( Int_t iCh, UInt_t mask )
{
  /// release parameters matching mask, for all detector elements in a given chamber, counting from 1

  // check boundaries
  if( iCh < 1 || iCh > 10 )
  { AliFatal( Form( "Invalid chamber index %i", iCh ) ); }

  // get first and last element
  const Int_t iDetElemFirst = fgSNDetElemCh[iCh-1];
  const Int_t iDetElemLast = fgSNDetElemCh[iCh];
  for( Int_t i = iDetElemFirst; i < iDetElemLast; ++i )
  {

    AliInfo( Form( "Releasing %s for detector element %i", GetParameterMaskString(mask).Data(), i ) );

    if( mask & ParX )  ReleaseParameter(i, 0);
    if( mask & ParY )  ReleaseParameter(i, 1);
    if( mask & ParTZ ) ReleaseParameter(i, 2);
    if( mask & ParZ )  ReleaseParameter(i, 3);

  }
}

//_____________________________________________________________________
void AliMUONAlignment::ReleaseDetElem( Int_t iDetElemId, UInt_t mask )
{
  /// release parameters matching mask, for a given detector element, counting from 0
  const Int_t iDet( GetDetElemNumber( iDetElemId ) );
  if ( mask & ParX )  ReleaseParameter(iDet, 0);
  if ( mask & ParY )  ReleaseParameter(iDet, 1);
  if ( mask & ParTZ ) ReleaseParameter(iDet, 2);
  if ( mask & ParZ )  ReleaseParameter(iDet, 3);

}

//______________________________________________________________________
void AliMUONAlignment::ReleaseParameter( Int_t iPar )
{

  /// release a given parameter, counting from 0
  if( fInitialized )
  { AliFatal( "Millepede already initialized" ); }

  fGlobalParameterStatus[iPar] = kFreeParId;

}

//_____________________________________________________________________
void AliMUONAlignment::GroupChamber( Int_t iCh, UInt_t mask )
{
  /// group parameters matching mask for all detector elements in a given chamber, counting from 1
  if( iCh < 1 || iCh > 10 )
  { AliFatal( Form( "Invalid chamber index %i", iCh ) ); }

  const Int_t detElemMin = 100*iCh;
  const Int_t detElemMax = 100*iCh + fgNDetElemCh[iCh]-1;
  GroupDetElems( detElemMin, detElemMax, mask );

}

//_____________________________________________________________________
void AliMUONAlignment::GroupDetElems( Int_t detElemMin, Int_t detElemMax, UInt_t mask )
{
  /// group parameters matching mask for all detector elements between min and max
  // check number of detector elements
  const Int_t nDetElem = detElemMax - detElemMin + 1;
  if( nDetElem<2 )
  { AliFatal( Form( "Requested group of DEs %d-%d contains less than 2 DE's", detElemMin, detElemMax ) ); }

  // create list
  Int_t* detElemList = new int[nDetElem];
  for( Int_t i = 0; i < nDetElem; ++i )
  { detElemList[i] = detElemMin+i; }

  // group
  GroupDetElems( detElemList, nDetElem, mask );
  delete[] detElemList;

}

//_____________________________________________________________________
void AliMUONAlignment::GroupDetElems( Int_t* detElemList, Int_t nDetElem, UInt_t mask )
{
  /// group parameters matching mask for all detector elements in list
  if( fInitialized )
  { AliFatal( "Millepede already initialized" ); }

  const Int_t iDeBase( GetDetElemNumber( detElemList[0] ) );
  for( Int_t i = 0; i < nDetElem; ++i )
  {
    const Int_t iDeCurrent( GetDetElemNumber( detElemList[i] ) );
    if( mask & ParX ) fGlobalParameterStatus[iDeCurrent*fgNParCh + 0] = (i==0) ?  kGroupBaseId : (kGroupBaseId-iDeBase-1);
    if( mask & ParY ) fGlobalParameterStatus[iDeCurrent*fgNParCh + 1] = (i==0) ?  kGroupBaseId : (kGroupBaseId-iDeBase-1);
    if( mask & ParTZ ) fGlobalParameterStatus[iDeCurrent*fgNParCh + 2] = (i==0) ?  kGroupBaseId : (kGroupBaseId-iDeBase-1);
    if( mask & ParZ ) fGlobalParameterStatus[iDeCurrent*fgNParCh + 3] = (i==0) ?  kGroupBaseId : (kGroupBaseId-iDeBase-1);

    if( i== 0 ) AliInfo( Form( "Creating new group for detector %i and variable %s", detElemList[i], GetParameterMaskString( mask ).Data() ) );
    else AliInfo( Form( "Adding detector element %i to current group", detElemList[i] ) );
  }

}

//______________________________________________________________________
void AliMUONAlignment::SetChamberNonLinear( Int_t iCh, UInt_t mask )
{
  /// Set parameters matching mask as non linear, for all detector elements in a given chamber, counting from 1
  const Int_t iDetElemFirst = fgSNDetElemCh[iCh-1];
  const Int_t iDetElemLast = fgSNDetElemCh[iCh];
  for( Int_t i = iDetElemFirst; i < iDetElemLast; ++i )
  {

      if( mask & ParX ) SetParameterNonLinear(i, 0);
      if( mask & ParY ) SetParameterNonLinear(i, 1);
      if( mask & ParTZ ) SetParameterNonLinear(i, 2);
      if( mask & ParZ ) SetParameterNonLinear(i, 3);

  }

}

//_____________________________________________________________________
void AliMUONAlignment::SetDetElemNonLinear( Int_t iDetElemId, UInt_t mask )
{
  /// Set parameters matching mask as non linear, for a given detector element, counting from 0
  const Int_t iDet( GetDetElemNumber( iDetElemId ) );
  if ( mask & ParX )  SetParameterNonLinear(iDet, 0);
  if ( mask & ParY )  SetParameterNonLinear(iDet, 1);
  if ( mask & ParTZ ) SetParameterNonLinear(iDet, 2);
  if ( mask & ParZ )  SetParameterNonLinear(iDet, 3);

}

//______________________________________________________________________
void AliMUONAlignment::SetParameterNonLinear( Int_t iPar )
{
  /// Set nonlinear flag for parameter iPar
  if( !fInitialized )
  { AliFatal( "Millepede not initialized" ); }

  fMillepede->SetNonLinear( iPar );
  AliInfo( Form( "Parameter %i set to non linear", iPar ) );
}

//______________________________________________________________________
void AliMUONAlignment::AddConstraints( const Bool_t *lChOnOff, UInt_t mask )
{
  /// Add constraint equations for selected chambers and degrees of freedom

  Array fConstraintX;
  Array fConstraintY;
  Array fConstraintTZ;
  Array fConstraintZ;

  for( Int_t i = 0; i < fgNDetElem; ++i )
  {

    // get chamber matching detector
    const Int_t iCh( GetChamberId(i) );
    if (lChOnOff[iCh-1])
    {

      if( mask & ParX ) fConstraintX.values[i*fgNParCh+0]=1.0;
      if( mask & ParY ) fConstraintY.values[i*fgNParCh+1]=1.0;
      if( mask & ParTZ ) fConstraintTZ.values[i*fgNParCh+2]=1.0;
      if( mask & ParZ ) fConstraintTZ.values[i*fgNParCh+3]=1.0;

    }
  }

  if( mask & ParX ) AddConstraint(fConstraintX.values,0.0);
  if( mask & ParY ) AddConstraint(fConstraintY.values,0.0);
  if( mask & ParTZ ) AddConstraint(fConstraintTZ.values,0.0);
  if( mask & ParZ ) AddConstraint(fConstraintZ.values,0.0);

}

//______________________________________________________________________
void AliMUONAlignment::AddConstraints(const Bool_t *lChOnOff, const Bool_t *lVarXYT, UInt_t sidesMask )
{
  /*
  questions:
  - is there not redundancy/inconsistency between lDetTLBR and lSpecLROnOff ? shouldn't we use only lDetTLBR ?
  - why is weight ignored for ConstrainT and ConstrainB
  - why is there no constrain on z
  */

  /// Add constraint equations for selected chambers, degrees of freedom and detector half
  Double_t lMeanY = 0.;
  Double_t lSigmaY = 0.;
  Double_t lMeanZ = 0.;
  Double_t lSigmaZ = 0.;
  Int_t lNDetElem = 0;

  for( Int_t i = 0; i < fgNDetElem; ++i )
  {

    // get chamber matching detector
    const Int_t iCh( GetChamberId(i) );

    // skip detector if chamber is off
    if( lChOnOff[iCh-1] ) continue;

    // get detector element id from detector element number
    const Int_t lDetElemNumber = i-fgSNDetElemCh[iCh-1];
    const Int_t lDetElemId = iCh*100+lDetElemNumber;

    // skip detector if its side is off
    // stations 1 and 2
    if( iCh>=1 && iCh<=4 )
    {
      if( lDetElemNumber == 0 && !( sidesMask & SideTopRight ) ) continue;
      if( lDetElemNumber == 1 && !( sidesMask & SideTopLeft ) ) continue;
      if( lDetElemNumber == 2 && !( sidesMask & SideBottomLeft ) ) continue;
      if( lDetElemNumber == 3 && !( sidesMask & SideBottomRight ) ) continue;
    }

    // station 3
    if (iCh>=5 && iCh<=6)
    {
      if( lDetElemNumber >= 0 && lDetElemNumber <= 4 && !( sidesMask & SideTopRight ) ) continue;
      if( lDetElemNumber >= 5 && lDetElemNumber <= 10 && !( sidesMask & SideTopLeft ) ) continue;
      if( lDetElemNumber >= 11 && lDetElemNumber <= 13 && !( sidesMask & SideBottomLeft ) ) continue;
      if( lDetElemNumber >= 14 && lDetElemNumber <= 17 && !( sidesMask & SideBottomRight ) ) continue;
    }

    // stations 4 and 5
    if (iCh>=7 && iCh<=10)
    {
      if( lDetElemNumber >= 0 && lDetElemNumber <= 6 && !( sidesMask & SideTopRight ) ) continue;
      if( lDetElemNumber >= 7 && lDetElemNumber <= 13 && !( sidesMask & SideTopLeft ) ) continue;
      if( lDetElemNumber >= 14 && lDetElemNumber <= 19 && !( sidesMask & SideBottomLeft ) ) continue;
      if( lDetElemNumber >= 20 && lDetElemNumber <= 25 && !( sidesMask & SideBottomRight ) ) continue;
    }

    // get global x, y and z position
    Double_t lDetElemGloX = 0.;
    Double_t lDetElemGloY = 0.;
    Double_t lDetElemGloZ = 0.;
    fTransform->Local2Global( lDetElemId, 0, 0, 0, lDetElemGloX, lDetElemGloY, lDetElemGloZ );

    // increment mean Y, mean Z, sigmas and number of accepted detectors
    lMeanY += lDetElemGloY;
    lSigmaY += lDetElemGloY*lDetElemGloY;
    lMeanZ += lDetElemGloZ;
    lSigmaZ += lDetElemGloZ*lDetElemGloZ;
    lNDetElem++;

  }

  // calculate mean values
  lMeanY /= lNDetElem;
  lSigmaY /= lNDetElem;
  lSigmaY = TMath::Sqrt(lSigmaY-lMeanY*lMeanY);
  lMeanZ /= lNDetElem;
  lSigmaZ /= lNDetElem;
  lSigmaZ = TMath::Sqrt(lSigmaZ-lMeanZ*lMeanZ);
  AliInfo( Form( "Used %i DetElem, MeanZ= %f , SigmaZ= %f", lNDetElem,lMeanZ,lSigmaZ ) );

  // create all possible arrays
  Array fConstraintX[4];  //Array for constraint equation X
  Array fConstraintY[4];  //Array for constraint equation Y
  Array fConstraintP[4];  //Array for constraint equation P
  Array fConstraintXZ[4];  //Array for constraint equation X vs Z
  Array fConstraintYZ[4];  //Array for constraint equation Y vs Z
  Array fConstraintPZ[4];  //Array for constraint equation P vs Z

  // do we really need these ?
  Array fConstraintXY[4];  //Array for constraint equation X vs Y
  Array fConstraintYY[4];  //Array for constraint equation Y vs Y
  Array fConstraintPY[4];  //Array for constraint equation P vs Y

  // fill Bool_t sides array based on masks, for convenience
  Bool_t lDetTLBR[4];
  lDetTLBR[0] = sidesMask & SideTop;
  lDetTLBR[1] = sidesMask & SideLeft;
  lDetTLBR[2] = sidesMask & SideBottom;
  lDetTLBR[3] = sidesMask & SideRight;

  for( Int_t i = 0; i < fgNDetElem; ++i )
  {

    // get chamber matching detector
    const Int_t iCh( GetChamberId(i) );

    // skip detector if chamber is off
    if( !lChOnOff[iCh-1] ) continue;

    // get detector element id from detector element number
    const Int_t lDetElemNumber = i-fgSNDetElemCh[iCh-1];
    const Int_t lDetElemId = iCh*100+lDetElemNumber;

    // get global x, y and z position
    Double_t lDetElemGloX = 0.;
    Double_t lDetElemGloY = 0.;
    Double_t lDetElemGloZ = 0.;
    fTransform->Local2Global( lDetElemId, 0, 0, 0, lDetElemGloX, lDetElemGloY, lDetElemGloZ );

    // loop over sides
    for( Int_t iSide = 0; iSide < 4; iSide++ )
    {

      // skip if side is not selected
      if( !lDetTLBR[iSide] ) continue;

      // skip detector if it is not in the selected side
      // stations 1 and 2
      if( iCh>=1 && iCh<=4 )
      {
        if( lDetElemNumber == 0 && !(iSide == 0 || iSide == 3) ) continue; // top-right
        if( lDetElemNumber == 1 && !(iSide == 0 || iSide == 1) ) continue; // top-left
        if( lDetElemNumber == 2 && !(iSide == 2 || iSide == 1) ) continue; // bottom-left
        if( lDetElemNumber == 3 && !(iSide == 2 || iSide == 3) ) continue; // bottom-right
      }

      // station 3
      if (iCh>=5 && iCh<=6)
      {
        if( lDetElemNumber >= 0 && lDetElemNumber <= 4 && !(iSide == 0 || iSide == 3) ) continue; // top-right
        if( lDetElemNumber >= 5 && lDetElemNumber <= 9 && !(iSide == 0 || iSide == 1) ) continue; // top-left
        if( lDetElemNumber >= 10 && lDetElemNumber <= 13 && !(iSide == 2 || iSide == 1) ) continue; // bottom-left
        if( lDetElemNumber >= 14 && lDetElemNumber <= 17 && !(iSide == 2 || iSide == 3) ) continue; // bottom-right
      }

      // stations 4 and 5
      if (iCh>=7 && iCh<=10)
      {
        if( lDetElemNumber >= 0 && lDetElemNumber <= 6 && !(iSide == 0 || iSide == 3) ) continue; // top-right
        if( lDetElemNumber >= 7 && lDetElemNumber <= 13 && !(iSide == 0 || iSide == 1) ) continue; // top-left
        if( lDetElemNumber >= 14 && lDetElemNumber <= 19 && !(iSide == 2 || iSide == 1) ) continue; // bottom-left
        if( lDetElemNumber >= 20 && lDetElemNumber <= 25 && !(iSide == 2 || iSide == 3) ) continue; // bottom-right
      }

      // constrain x
      if( lVarXYT[0] ) fConstraintX[iSide].values[i*fgNParCh+0] = 1;

      // constrain y
      if( lVarXYT[1] ) fConstraintY[iSide].values[i*fgNParCh+1] = 1;

      // constrain phi (rotation around z)
      if( lVarXYT[2] ) fConstraintP[iSide].values[i*fgNParCh+2] = 1;

      // x-z shearing
      if( lVarXYT[3] ) fConstraintXZ[iSide].values[i*fgNParCh+0] = (lDetElemGloZ-lMeanZ)/lSigmaZ;

      // y-z shearing
      if( lVarXYT[4] ) fConstraintYZ[iSide].values[i*fgNParCh+1] = (lDetElemGloZ-lMeanZ)/lSigmaZ;

      // phi-z shearing
      if( lVarXYT[5] ) fConstraintPZ[iSide].values[i*fgNParCh+2] = (lDetElemGloZ-lMeanZ)/lSigmaZ;

      // x-y shearing
      if( lVarXYT[6] ) fConstraintXY[iSide].values[i*fgNParCh+0] = (lDetElemGloY-lMeanY)/lSigmaY;

      // y-y shearing
      if( lVarXYT[7] ) fConstraintYY[iSide].values[i*fgNParCh+1] = (lDetElemGloY-lMeanY)/lSigmaY;

      // phi-y shearing
      if( lVarXYT[8] ) fConstraintPY[iSide].values[i*fgNParCh+2] = (lDetElemGloY-lMeanY)/lSigmaY;

    }

  }

  // pass constraints to millepede
  for( Int_t iSide = 0; iSide < 4; iSide++ )
  {
    // skip if side is not selected
    if( !lDetTLBR[iSide] ) continue;

    if( lVarXYT[0] ) AddConstraint(fConstraintX[iSide].values,0.0);
    if( lVarXYT[1] ) AddConstraint(fConstraintY[iSide].values,0.0);
    if( lVarXYT[2] ) AddConstraint(fConstraintP[iSide].values,0.0);
    if( lVarXYT[3] ) AddConstraint(fConstraintXZ[iSide].values,0.0);
    if( lVarXYT[4] ) AddConstraint(fConstraintYZ[iSide].values,0.0);
    if( lVarXYT[5] ) AddConstraint(fConstraintPZ[iSide].values,0.0);
    if( lVarXYT[6] ) AddConstraint(fConstraintXY[iSide].values,0.0);
    if( lVarXYT[7] ) AddConstraint(fConstraintYY[iSide].values,0.0);
    if( lVarXYT[8] ) AddConstraint(fConstraintPY[iSide].values,0.0);
  }

}

//______________________________________________________________________
void AliMUONAlignment::InitGlobalParameters(Double_t *par)
{
  /// Initialize global parameters with par array
  if( !fInitialized )
  { AliFatal( "Millepede is not initialized" ); }

  fMillepede->SetGlobalParameters(par);
}

//______________________________________________________________________
void AliMUONAlignment::SetAllowedVariation( Int_t iPar, Double_t value )
{
  /// "Encouraged" variation for degrees of freedom
  // check initialization
  if( fInitialized )
  { AliFatal( "Millepede already initialized" ); }

  // check initialization
  if( !(iPar >= 0 && iPar < fgNParCh ) )
  { AliFatal( Form( "Invalid index: %i", iPar ) ); }

  fAllowVar[iPar] = value;
}

//______________________________________________________________________
void AliMUONAlignment::SetSigmaXY(Double_t sigmaX, Double_t sigmaY)
{

  /// Set expected measurement resolution
  fSigma[0] = sigmaX;
  fSigma[1] = sigmaY;

  // print
  for( Int_t i=0; i<2; ++i )
  { AliInfo( Form( "fSigma[%i]=%f", i, fSigma[i] ) ); }

}

//_____________________________________________________
void AliMUONAlignment::GlobalFit( Double_t *parameters, Double_t *errors, Double_t *pulls )
{

  /// Call global fit; Global parameters are stored in parameters
  fMillepede->GlobalFit( parameters, errors, pulls );

  AliInfo( "Done fitting global parameters" );
  for( int iDet=0; iDet<fgNDetElem; ++iDet )
  {
    AliInfo( Form( "%d\t %f\t %f\t %f\t %f",
      iDet,
      parameters[iDet*fgNParCh+0],parameters[iDet*fgNParCh+1],
      parameters[iDet*fgNParCh+3],parameters[iDet*fgNParCh+2]
      ) );
  }

}

//_____________________________________________________
void AliMUONAlignment::PrintGlobalParameters() const
{ fMillepede->PrintGlobalParameters(); }

//_____________________________________________________
Double_t AliMUONAlignment::GetParError(Int_t iPar) const
{ return fMillepede->GetParError(iPar); }

//______________________________________________________________________
AliMUONGeometryTransformer* AliMUONAlignment::ReAlign(
  const AliMUONGeometryTransformer * transformer,
  const double *misAlignments, Bool_t )
{

  /// Returns a new AliMUONGeometryTransformer with the found misalignments
  /// applied.

  // Takes the internal geometry module transformers, copies them
  // and gets the Detection Elements from them.
  // Takes misalignment parameters and applies these
  // to the local transform of the Detection Element
  // Obtains the global transform by multiplying the module transformer
  // transformation with the local transformation
  // Applies the global transform to a new detection element
  // Adds the new detection element to a new module transformer
  // Adds the new module transformer to a new geometry transformer
  // Returns the new geometry transformer

  Double_t lModuleMisAlignment[fgNParCh] = {0};
  Double_t lDetElemMisAlignment[fgNParCh] = {0};
  const TClonesArray* oldMisAlignArray( transformer->GetMisAlignmentData() );

  AliMUONGeometryTransformer *newGeometryTransformer = new AliMUONGeometryTransformer();
  for( Int_t iMt = 0; iMt < transformer->GetNofModuleTransformers(); ++iMt )
  {

    // module transformers
    const AliMUONGeometryModuleTransformer *kModuleTransformer = transformer->GetModuleTransformer(iMt, kTRUE);

    AliMUONGeometryModuleTransformer *newModuleTransformer = new AliMUONGeometryModuleTransformer(iMt);
    newGeometryTransformer->AddModuleTransformer(newModuleTransformer);

    // get transformation
    TGeoHMatrix deltaModuleTransform( DeltaTransform( lModuleMisAlignment ) );

    // update module
    TGeoHMatrix moduleTransform( *kModuleTransformer->GetTransformation() );
    TGeoHMatrix newModuleTransform( AliMUONGeometryBuilder::Multiply( deltaModuleTransform, moduleTransform ) );
    newModuleTransformer->SetTransformation(newModuleTransform);

    // Get matching old alignment and update current matrix accordingly
    if( oldMisAlignArray )
    {

      const AliAlignObjMatrix* oldAlignObj(0);
      const Int_t moduleId( kModuleTransformer->GetModuleId() );
      const Int_t volId = AliGeomManager::LayerToVolUID(AliGeomManager::kMUON, moduleId );
      for( Int_t pos = 0; pos < oldMisAlignArray->GetEntriesFast(); ++pos )
      {

        const AliAlignObjMatrix* localAlignObj( dynamic_cast<const AliAlignObjMatrix*>(oldMisAlignArray->At( pos ) ) );
        if( localAlignObj && localAlignObj->GetVolUID() == volId )
        {
          oldAlignObj = localAlignObj;
          break;
        }

      }

      // multiply
      if( oldAlignObj )
      {

        TGeoHMatrix oldMatrix;
        oldAlignObj->GetMatrix( oldMatrix );
        deltaModuleTransform.Multiply( &oldMatrix );

      }

    }

    // Create module mis alignment matrix
    newGeometryTransformer ->AddMisAlignModule(kModuleTransformer->GetModuleId(), deltaModuleTransform);

    AliMpExMap *detElements = kModuleTransformer->GetDetElementStore();

    TIter next(detElements->CreateIterator());
    AliMUONGeometryDetElement* detElement;
    Int_t iDe(-1);
    while ( ( detElement = static_cast<AliMUONGeometryDetElement*>(next()) ) )
    {
      ++iDe;
      // make a new detection element
      AliMUONGeometryDetElement *newDetElement = new AliMUONGeometryDetElement(detElement->GetId(), detElement->GetVolumePath());
      TString lDetElemName(detElement->GetDEName());
      lDetElemName.ReplaceAll("DE","");

      // store detector element id and number
      const Int_t iDetElemId = lDetElemName.Atoi();
      if( !DetElemIsValid( iDetElemId ) )
      {
        AliInfo( Form( "Skipping invalid detector element %i", iDetElemId ) );
        continue;
      }

      const Int_t iDetElemNumber( GetDetElemNumber( iDetElemId ) );

      for( int i=0; i<fgNParCh; ++i )
      {
        lDetElemMisAlignment[i] = 0.0;
        if( iMt<fgNTrkMod ) { lDetElemMisAlignment[i] =  misAlignments[iDetElemNumber*fgNParCh+i]; }
      }

      // get transformation
      TGeoHMatrix deltaGlobalTransform( DeltaTransform( lDetElemMisAlignment ) );

      // update module
      TGeoHMatrix globalTransform( *detElement->GetGlobalTransformation() );
      TGeoHMatrix newGlobalTransform( AliMUONGeometryBuilder::Multiply( deltaGlobalTransform, globalTransform ) );
      newDetElement->SetGlobalTransformation( newGlobalTransform );
      newModuleTransformer->GetDetElementStore()->Add(newDetElement->GetId(), newDetElement);

      // Get matching old alignment and update current matrix accordingly
      if( oldMisAlignArray )
      {

        const AliAlignObjMatrix* oldAlignObj(0);
        const int detElemId( detElement->GetId() );
        const Int_t volId = AliGeomManager::LayerToVolUID(AliGeomManager::kMUON, detElemId );
        for( Int_t pos = 0; pos < oldMisAlignArray->GetEntriesFast(); ++pos )
        {

          const AliAlignObjMatrix* localAlignObj( dynamic_cast<const AliAlignObjMatrix*>(oldMisAlignArray->At( pos ) ) );
          if( localAlignObj && localAlignObj->GetVolUID() == volId )
          {
            oldAlignObj = localAlignObj;
            break;
          }

        }

        // multiply
        if( oldAlignObj )
        {

          TGeoHMatrix oldMatrix;
          oldAlignObj->GetMatrix( oldMatrix );
          deltaGlobalTransform.Multiply( &oldMatrix );

        }

      }

      // Create misalignment matrix
      newGeometryTransformer->AddMisAlignDetElement(detElement->GetId(), deltaGlobalTransform);

    }

    newGeometryTransformer->AddModuleTransformer(newModuleTransformer);
  }

  return newGeometryTransformer;

}

//______________________________________________________________________
void AliMUONAlignment::SetAlignmentResolution( const TClonesArray* misAlignArray, Int_t rChId, Double_t chResX, Double_t chResY, Double_t deResX, Double_t deResY )
{

  /// Set alignment resolution to misalign objects to be stored in CDB
  /// if rChId is > 0 set parameters for this chamber only, counting from 1
  TMatrixDSym mChCorrMatrix(6);
  mChCorrMatrix[0][0]=chResX*chResX;
  mChCorrMatrix[1][1]=chResY*chResY;

  TMatrixDSym mDECorrMatrix(6);
  mDECorrMatrix[0][0]=deResX*deResX;
  mDECorrMatrix[1][1]=deResY*deResY;

  AliAlignObjMatrix *alignMat = 0x0;

  for( Int_t chId = 0; chId <= 9; ++chId )
  {

    // skip chamber if selection is valid, and does not match
    if( rChId > 0 && chId+1 != rChId ) continue;

    TString chName1;
    TString chName2;
    if (chId<4)
    {

      chName1 = Form("GM%d",chId);
      chName2 = Form("GM%d",chId);

    } else {

      chName1 = Form("GM%d",4+(chId-4)*2);
      chName2 = Form("GM%d",4+(chId-4)*2+1);

    }

    for( int i=0; i<misAlignArray->GetEntries(); ++i )
    {

      alignMat = (AliAlignObjMatrix*)misAlignArray->At(i);
      TString volName(alignMat->GetSymName());
      if((volName.Contains(chName1)&&
        ((volName.Last('/')==volName.Index(chName1)+chName1.Length())||
        (volName.Length()==volName.Index(chName1)+chName1.Length())))||
        (volName.Contains(chName2)&&
        ((volName.Last('/')==volName.Index(chName2)+chName2.Length())||
        (volName.Length()==volName.Index(chName2)+chName2.Length()))))
      {

        volName.Remove(0,volName.Last('/')+1);
        if (volName.Contains("GM")) alignMat->SetCorrMatrix(mChCorrMatrix);
        else if (volName.Contains("DE")) alignMat->SetCorrMatrix(mDECorrMatrix);

      }

    }

  }

}


//_____________________________________________________
void AliMUONAlignment::FillDetElemData( AliMUONVCluster* cluster )
{

  /// Get information of current detection element
  // get detector element number from Alice ID
  const Int_t detElemId = cluster->GetDetElemId();
  fDetElemNumber = GetDetElemNumber( detElemId );

  // get detector element
  const AliMUONGeometryDetElement* detElement = fTransform->GetDetElement( detElemId );

  /*
  get the global transformation matrix and store its inverse, in order to manually perform
  the global to Local transformations needed to calculate the derivatives
  */
  fGeoCombiTransInverse = detElement->GetGlobalTransformation()->Inverse();

}

//______________________________________________________________________
void AliMUONAlignment::FillRecPointData( AliMUONVCluster* cluster )
{

  /// Get information of current cluster
  fClustPos[0] = cluster->GetX();
  fClustPos[1] = cluster->GetY();
  fClustPos[2] = cluster->GetZ();

}

//______________________________________________________________________
void AliMUONAlignment::FillTrackParamData( AliMUONTrackParam* trackParam )
{

  /// Get information of current track at current cluster
  fTrackPos[0] = trackParam->GetNonBendingCoor();
  fTrackPos[1] = trackParam->GetBendingCoor();
  fTrackPos[2] = trackParam->GetZ();
  fTrackSlope[0] = trackParam->GetNonBendingSlope();
  fTrackSlope[1] = trackParam->GetBendingSlope();

}

//______________________________________________________________________
void AliMUONAlignment::LocalEquationX( void )
{
  /// local equation along X

  // 'inverse' (GlobalToLocal) rotation matrix
  const Double_t* r( fGeoCombiTransInverse.GetRotationMatrix() );

  // local derivatives
  SetLocalDerivative( 0, r[0] );
  SetLocalDerivative( 1, r[0]*(fTrackPos[2] - fTrackPos0[2]) );
  SetLocalDerivative( 2, r[1] );
  SetLocalDerivative( 3, r[1]*(fTrackPos[2] - fTrackPos0[2]) );

  // global derivatives
  /*
  alignment parameters are
  0: delta_x
  1: delta_y
  2: delta_phiz
  3: delta_z
  */

  SetGlobalDerivative( fDetElemNumber*fgNParCh + 0, -r[0] );
  SetGlobalDerivative( fDetElemNumber*fgNParCh + 1, -r[1] );

  if( fBFieldOn )
  {

    // use local position for derivatives vs 'delta_phi_z'
    SetGlobalDerivative( fDetElemNumber*fgNParCh + 2, -r[1]*fTrackPos[0] + r[0]*fTrackPos[1] );

    // use local slopes for derivatives vs 'delta_z'
    SetGlobalDerivative( fDetElemNumber*fgNParCh + 3, r[0]*fTrackSlope[0] + r[1]*fTrackSlope[1] );

  } else {

    // local copy of extrapolated track positions
    const Double_t trackPosX = fTrackPos0[0]+fTrackSlope0[0]*( fTrackPos[2]-fTrackPos0[2] );
    const Double_t trackPosY = fTrackPos0[1]+fTrackSlope0[1]*( fTrackPos[2]-fTrackPos0[2] );

    // use properly extrapolated position for derivatives vs 'delta_phi_z'
    SetGlobalDerivative( fDetElemNumber*fgNParCh + 2, -r[1]*trackPosX + r[0]*trackPosY );

    // use slopes at origin for derivatives vs 'delta_z'
    SetGlobalDerivative( fDetElemNumber*fgNParCh + 3, r[0]*fTrackSlope0[0] + r[1]*fTrackSlope0[1] );

  }

  // store local equation
  fMillepede->SetLocalEquation( fGlobalDerivatives, fLocalDerivatives, fMeas[0], fSigma[0] );

}

//______________________________________________________________________
void AliMUONAlignment::LocalEquationY( void )
{
  /// local equation along Y

  // 'inverse' (GlobalToLocal) rotation matrix
  const Double_t* r( fGeoCombiTransInverse.GetRotationMatrix() );

  // store local derivatives
  SetLocalDerivative( 0, r[3] );
  SetLocalDerivative( 1, r[3]*(fTrackPos[2] - fTrackPos0[2] ) );
  SetLocalDerivative( 2, r[4] );
  SetLocalDerivative( 3, r[4]*(fTrackPos[2] - fTrackPos0[2] ) );

  // set global derivatives
  SetGlobalDerivative( fDetElemNumber*fgNParCh + 0, -r[3]);
  SetGlobalDerivative( fDetElemNumber*fgNParCh + 1, -r[4]);

  if( fBFieldOn )
  {

    // use local position for derivatives vs 'delta_phi'
    SetGlobalDerivative( fDetElemNumber*fgNParCh + 2, -r[4]*fTrackPos[0] + r[3]*fTrackPos[1]);

    // use local slopes for derivatives vs 'delta_z'
    SetGlobalDerivative( fDetElemNumber*fgNParCh + 3, r[3]*fTrackSlope[0]+r[4]*fTrackSlope[1] );

  } else {

    // local copy of extrapolated track positions
    const Double_t trackPosX = fTrackPos0[0]+fTrackSlope0[0]*( fTrackPos[2]-fTrackPos0[2] );
    const Double_t trackPosY = fTrackPos0[1]+fTrackSlope0[1]*( fTrackPos[2]-fTrackPos0[2] );

    // use properly extrapolated position for derivatives vs 'delta_phi'
    SetGlobalDerivative( fDetElemNumber*fgNParCh + 2, -r[4]*trackPosX + r[3]*trackPosY );

    // use slopes at origin for derivatives vs 'delta_z'
    SetGlobalDerivative( fDetElemNumber*fgNParCh + 3, r[3]*fTrackSlope0[0]+r[4]*fTrackSlope0[1] );

  }

  // store local equation
  fMillepede->SetLocalEquation( fGlobalDerivatives, fLocalDerivatives, fMeas[1], fSigma[1] );

}

//_________________________________________________________________________
TGeoCombiTrans AliMUONAlignment::DeltaTransform( const double *lMisAlignment) const
{
  /// Get Delta Transformation, based on alignment parameters

  // translation
  const TGeoTranslation deltaTrans( lMisAlignment[0], lMisAlignment[1], lMisAlignment[3]);

  // rotation
  TGeoRotation deltaRot;
  deltaRot.RotateZ(lMisAlignment[2]*180./TMath::Pi());

  // combined rotation and translation.
  return TGeoCombiTrans(deltaTrans,deltaRot);

}

//______________________________________________________________________
void AliMUONAlignment::AddConstraint(Double_t *par, Double_t value)
{
  /// Constrain equation defined by par to value
  if( !fInitialized )
  { AliFatal( "Millepede is not initialized" ); }

  fMillepede->SetGlobalConstraint(par, value);
}

//______________________________________________________________________
Bool_t AliMUONAlignment::DetElemIsValid( Int_t iDetElemId ) const
{
  /// return true if given detector element is valid (and belongs to muon tracker)
  const Int_t iCh = iDetElemId/100;
  const Int_t iDet = iDetElemId%100;
  return ( iCh > 0 && iCh <= fgNCh && iDet < fgNDetElemCh[iCh-1] );
}

//______________________________________________________________________
Int_t AliMUONAlignment::GetDetElemNumber( Int_t iDetElemId ) const
{
  /// get det element number from ID
  // get chamber and element number in chamber
  const Int_t iCh = iDetElemId/100;
  const Int_t iDet = iDetElemId%100;

  // make sure detector index is valid
  if( !( iCh > 0 && iCh <= fgNCh && iDet < fgNDetElemCh[iCh-1] ) )
  { AliFatal( Form( "Invalid detector element id: %i", iDetElemId ) ); }

  // add number of detectors up to this chamber
  return iDet + fgSNDetElemCh[iCh-1];

}

//______________________________________________________________________
Int_t AliMUONAlignment::GetChamberId( Int_t iDetElemNumber ) const
{
  /// get chamber (counting from 1) matching a given detector element id
  Int_t iCh( 0 );
  for( iCh=0; iCh<fgNCh; iCh++ )
  { if( iDetElemNumber < fgSNDetElemCh[iCh] ) break; }

  return iCh;
}

//______________________________________________________________________
TString AliMUONAlignment::GetParameterMaskString( UInt_t mask ) const
{
  TString out;
  if( mask & ParX ) out += "X";
  if( mask & ParY ) out += "Y";
  if( mask & ParZ ) out += "Z";
  if( mask & ParTZ ) out += "T";
  return out;
}

//______________________________________________________________________
TString AliMUONAlignment::GetSidesMaskString( UInt_t mask ) const
{
  TString out;
  if( mask & SideTop ) out += "T";
  if( mask & SideLeft ) out += "L";
  if( mask & SideBottom ) out += "B";
  if( mask & SideRight ) out += "R";
  return out;
}
