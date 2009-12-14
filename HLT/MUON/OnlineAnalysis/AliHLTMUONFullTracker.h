#ifndef ALIHLTMUONFULLTRACKER_H
#define ALIHLTMUONFULLTRACKER_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */
/**********************************************************************
 Created on : 08/12/2009
 Purpose    : First version implementation of the Full tracker for dHLT.
 Author     : Indranil Das, HEP Division, SINP
 Email      : indra.das@saha.ac.in | indra.ehep@gmail.com
**********************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <iostream>

#include "TMap.h"
#include "TString.h"
#include "TVector3.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TGeoGlobalMagField.h"

#include "AliMagF.h"
#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "AliGRPObject.h"

#include "AliMUONTrackExtrap.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONConstants.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONTrackParam.h"

#include "AliMpDDLStore.h"
#include "AliMpSegmentation.h"
#include "AliMpCDB.h"

#include "AliCDBEntry.h"
#include "AliLog.h"

#include "AliRunInfo.h"

#include "AliHLTLogging.h"

#include "AliHLTMUONConstants.h"
#include "AliHLTMUONDataTypes.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONDataBlockReader.h"
#include "AliHLTMUONDataBlockWriter.h"
#endif

using namespace std;

class AliHLTMUONConstants;

class AliHLTMUONFullTracker : public AliHLTLogging
{

 public :
  ///Constructor
  AliHLTMUONFullTracker() ;
  ///Destructor
  ~AliHLTMUONFullTracker();

  ///Print message
  void Print();
  ///Set the input of trigrec blocks
  Bool_t SetInput(AliHLTInt32_t ddl, const AliHLTMUONTriggerRecordStruct  *data, AliHLTInt32_t size);
  ///Set the input of rechit blocks
  Bool_t SetInput(AliHLTInt32_t ddl, const AliHLTMUONRecHitStruct  *data, AliHLTInt32_t size);
  ///Main method to run and compute the tracks
  Bool_t Run(int iEvent,AliHLTMUONMansoTrackStruct *data, AliHLTUInt32_t& size);
  ///To be called once from DoInit method of component
  Bool_t Init();
  
 protected:

  /// copy constructor
  AliHLTMUONFullTracker(const AliHLTMUONFullTracker& rhs); 
  /// assignment operator
  AliHLTMUONFullTracker& operator=(const AliHLTMUONFullTracker& rhs); 

 private :

  /// intger pair needed for QuadTrackSeg method
  struct IntPair{
    Int_t fFirst,fSecond;
  };
  ///Structure for internal track segments
  struct TrackSeg{
    Int_t fIndex[4];
    AliHLTInt32_t fTrigRec;
  };

  ///Sructure for clusters
  struct Cluster{
    Float_t fX,fY,fZ;
    Float_t fErrX2,fErrY2;
  };

  
  static const Float_t TrackDetCoordinate[3]; /// set the constant value for third station position and size
  
  static const Double_t fgkAbsoedge[4] ;     /// edge of the absorber
  static const Double_t fgkRadLen[3] ;       /// radiation length of the main three matirials of the front absorber
  static const Double_t fgkRho[3] ;          /// density of the main three matirials of the front absorber
  static const Double_t fgkAtomicZ[3] ;      /// atomic number the main three matirials of the front absorber
  static const Double_t fgkAtomicA[3] ;      /// atomic mass of the main three matirials of the front absorber

  static const Int_t fgkMaxNofCellsPerCh ;      /// maximum number of cell are allowed to create
  static const Int_t fgkMaxNofPointsPerCh ;     /// maximim number of points per chamber
  static const Int_t fgkMaxNofCh ;              /// maximum number of chambrs
  static const Int_t fgkMaxNofTracks;           /// maximum number of allowed tracks
  static const Int_t fgkMaxNofConnectedTracks;  /// maximum number of back to front connected tracks
  
  /// Pointer to AliMUONGeometryTransformer
  AliMUONGeometryTransformer *fChamberGeometryTransformer;
  
  /// array of pointer to rechit data
  AliHLTMUONRecHitStruct ***fChPoint;
  ///array of pointer to trigrec data
  AliHLTMUONTriggerRecordStruct **fChPoint11;
  /// track segments at the rear part of the spectrometer
  TrackSeg *fBackTrackSeg;
  /// track segments close the part of interaction point  of ALICE
  TrackSeg *fFrontTrackSeg;
  
  /// Extrapolated x position in third station
  Float_t *fExtrapSt3X ;
  /// Extrapolated y position in third station
  Float_t *fExtrapSt3Y ;
  /// values of inclination angle of back track segments
  Float_t *fInclinationBack;

  /// nof connected tracks in front direction for each back track segments
  Int_t *fNofConnectedfrontTrackSeg ;
  /// Pointer to back to front segment mapping
  Int_t **fBackToFront;
  /// Charge of the tracks
  Float_t *fCharge;
  /// Number of points for each stations
  Int_t *fNofPoints ;
  /// track parameters;
  AliMUONTrackParam *fTrackParam ;

  /// Total number of points received from all rechit source
  Int_t fTotNofPoints;
  /// Total number of track segments
  Int_t fTotTrackSeg;
  /// Number of cells in QuadSeg
  Int_t fNofCells[2]; // nof cell count per station
  /// Check if overflowed
  Bool_t fOverflowed;
  /// number of back track segments
  Int_t fNofbackTrackSeg;
  /// number of front track segments
  Int_t fNoffrontTrackSeg;
  /// number of connected track segments
  Int_t fNofConnected ;

  /// Slat Track segments 
  Bool_t SlatTrackSeg();
  /// Quad Track segments 
  Bool_t QuadTrackSeg();
  /// Kalman Chi2 test
  Bool_t KalmanChi2Test();
  /// track extrapolation through  dipole magnet to connect front and back track seg
  Bool_t SelectFront();
  /// Propagate tracks
  void PropagateTracks(Double_t charge, Float_t& px, Float_t& py, Float_t& pz, 
		       Float_t& xr, Float_t& yr, Float_t& zr, Float_t zprop);
  /// extrapolate to origin
  Bool_t ExtrapolateToOrigin(Bool_t);
  /// Clean after each run
  Bool_t Clear();


  /// Angle calculate
  inline Double_t Angle(AliHLTMUONRecHitStruct *v1, AliHLTMUONRecHitStruct *v2);
  /// Subtracktion of two point
  inline void Sub(AliHLTMUONRecHitStruct *v1, AliHLTMUONRecHitStruct *v2, AliHLTMUONRecHitStruct *v3);
  /// Kalman Filter
  inline Double_t KalmanFilter(AliMUONTrackParam &trackParamAtCluster, Cluster *cluster);
  /// Try onecluster
  inline Double_t TryOneCluster(const AliMUONTrackParam &trackParam, Cluster* cluster,
				AliMUONTrackParam &trackParamAtCluster, Bool_t updatePropagator);
  inline Bool_t TryOneClusterFast(const AliMUONTrackParam &trackParam, Cluster* cluster);

  /// MCS effect correction
  inline void CorrectMCSEffectInAbsorber(AliMUONTrackParam* param,
					 Double_t xVtx, Double_t yVtx, Double_t zVtx,
					 Double_t absZBeg, 
					 Double_t f1, Double_t f2);
  /// Covariant handling function
  inline void Cov2CovP(const TMatrixD &param, TMatrixD &cov);
  /// Covariant handling function
  inline void CovP2Cov(const TMatrixD &param, TMatrixD &covP);
  /// Energy loss coreection in front absorber
  inline void CorrectELossEffectInAbsorber(AliMUONTrackParam* param, Double_t eLoss);
  /// Linear Extrapolation to Z position
  inline void LinearExtrapToZ(AliMUONTrackParam* trackParam, Double_t zEnd);
  /// Energy loss
  inline Double_t EnergyLossFluctuation2(Double_t pTotal, Double_t pathLength, Double_t rho, Double_t atomicA, Double_t atomicZ);
  /// Bethe Bloch formula of enrgy loss
  inline Double_t BetheBloch(Double_t pTotal, Double_t pathLength, Double_t rho, Double_t atomicA, Double_t atomicZ);
  
  /// Runge Kutta method of track extrapolation through mag field
  inline void OneStepRungekutta(Double_t charge, Double_t step,Double_t* vect, Double_t* vout);
  /// Helix3 method of track extrapolation through mag field
  inline void OneStepHelix3(Double_t field, Double_t step, Double_t *vect, Double_t *vout);				  
  /// Initialise GRP when running without reconstruction chain
  Bool_t InitGRP();
  /// Fill the tracks to output pointer
  Bool_t FillOutData(AliHLTMUONMansoTrackStruct *data, AliHLTUInt32_t& size);
  
};
#endif // ALIHLTMUONMANSOTRACKERFSM_H
