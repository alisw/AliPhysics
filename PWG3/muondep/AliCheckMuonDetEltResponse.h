#ifndef ALICHECKMUONDETELTRESPONSE_H 
#define ALICHECKMUONDETELTRESPONSE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup base
/// \class AliCheckMuonDetEltResponse
/// \brief tracking chamber efficiency from data
//Author: Nicolas LE BRIS - SUBATECH Nantes

#include <TH2F.h>
#include <TObject.h>
#include <TClonesArray.h>

class AliMUONTrackParam;
class AliMUONTrack;
class AliMUONVCluster;
class AliMUONGeometryTransformer;
class AliESDEvent;

class AliCheckMuonDetEltResponse : public TObject
{
public:
  AliCheckMuonDetEltResponse();
  AliCheckMuonDetEltResponse(const AliCheckMuonDetEltResponse& rhs);
  AliCheckMuonDetEltResponse& operator=(const AliCheckMuonDetEltResponse& rhs);
//Constructor:
  AliCheckMuonDetEltResponse(const AliMUONGeometryTransformer* transformer,
			     AliESDEvent* esd,
			     TClonesArray* detEltTDHistList,
			     TClonesArray* detEltTTHistList,
			     TClonesArray* chamberTDHistList,
			     TClonesArray* chamberTTHistList,
			     Bool_t isCosmic = kFALSE);

//Destructor:
  virtual ~AliCheckMuonDetEltResponse();


  void CheckDetEltResponse ();
  void TrackLoop ();
  void TrackParamLoop ();

  void SetCosmic(Bool_t isCosmic) {fIsCosmicData = isCosmic;};
  Bool_t IsCosmic() {return fIsCosmicData;};
 
private:
  
  void FillTDHistos (Int_t chamber, Int_t detElt,
		     Double_t posXL, Double_t posYL);

  void FillTTHistos (Int_t chamber, Int_t detElt,
		     Double_t posXL, Double_t posYL);

  void FindAndFillMissedDetElt (AliMUONTrackParam* extrapTrackParam,
				Int_t firstMissCh, Int_t lastChamber);

  void CoordinatesOfMissingCluster(Double_t x1, Double_t y1, Double_t z1,
				   Double_t x2, Double_t y2, Double_t z2,
				   Double_t& x, Double_t& y);

  Bool_t CoordinatesInDetEltSt12(Int_t DeId, Double_t x, Double_t y);
  Bool_t CoordinatesInDetEltSt345(Int_t DeId, Double_t x, Double_t y);


  Int_t fNCh; //!<Number of tracking chamber.
  Int_t fNSt; //!<Number of tracking station.
  Int_t fNDE; //!<Number of detection element in the tracking system.

  Int_t FromDetElt2iDet (Int_t chamber, Int_t detElt);
  Int_t FromDetElt2LocalId (Int_t chamber, Int_t detElt);

  const AliMUONGeometryTransformer* fTransformer; //!<Geometry transformer

  AliESDEvent* fESD;

  Int_t fNbrClustersCh[10];      //!<Number of clusters in the chamber [fChamberNbr].
  Int_t fTracksTotalNbr;         //!<Total number of tracks in the event.
  Int_t fTrackFilter[10];        //!<To select track for the efficiency calculation.
  Bool_t fIsCosmicData;          ///< Check if the data are cosmic rays (used only to cut cosmic shower at hte trigger level if true)

  TClonesArray     * fTrackParams;
  AliMUONTrackParam* fTrackParam;
  AliMUONVCluster  * fCluster;

  TClonesArray* fDetEltTDHistList; //!<List of histograms of the tracks detected in the detection elements 
  TClonesArray* fDetEltTTHistList; //!<List of histograms of the tracks which have passed through the detection elements
  TClonesArray* fChamberTDHistList; //!<List of histograms of the tracks detected in the chambers 
  TClonesArray* fChamberTTHistList; //!<List of histograms of the tracks which have passed through the chambers

  static const Int_t fNbrOfChamber;            ///< The total number of chamber in the tracking system.
  static const Int_t fNbrOfStation;            ///< The total number of station in the tracking system.
  static const Int_t fNbrOfDetectionElt[10];   ///< The total number of detection element in each chamber.
  static const Int_t fFirstDetectionElt[10];   ///< The Id of the first detection element of each chamber.
  static const Int_t fOffset;                  ///< fFirstDetectionElt[iChamber] = fOffset * (iChamber+1).
  static const Int_t fOverlapSize;             ///< Average size (in cm) of the overlap area between two detection eltement.
  static const Int_t fYSlatSize;               ///< Average size (in cm) of the overlap area between two detection eltement.

  ClassDef(AliCheckMuonDetEltResponse, 0)
};
#endif
