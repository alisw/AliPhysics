#ifndef ALICHECKMUONDETELTRESPONSE_H 
#define ALICHECKMUONDETELTRESPONSE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup base
/// \class AliCheckMuonDetEltResponse
/// \brief tracking chamber efficiency from data
//Author: Nicolas LE BRIS - SUBATECH Nantes

#include <TObject.h>

class AliESDEvent;
class AliMUONTrackParam;
class AliMUONTrack;
class AliMUONVCluster;
class AliMUONGeometryTransformer;
class TClonesArray;

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

  Int_t GetNbrUsableTracks() const {return fNbrUsableTracks;};
  void SetNbrUsableTracks(Int_t nUsable){fNbrUsableTracks = nUsable;};
  Bool_t IsCosmic() const {return fIsCosmicData;};
  void SetCosmic(Bool_t isCosmic) {fIsCosmicData = isCosmic;};
 
private:
  
  void FillTDHistos (Int_t chamber, Int_t detElt,
		     Double_t posXL, Double_t posYL);

  void FillTTHistos (Int_t chamber, Int_t detElt,
		     Double_t posXL, Double_t posYL);

  void FindAndFillMissedDetElt (AliMUONTrackParam* extrapTrackParam,
				Int_t firstMissCh, Int_t lastChamber);

  void CoordinatesOfMissingCluster(Double_t x1, Double_t y1, Double_t z1,
				   Double_t x2, Double_t y2, Double_t z2,
				   Double_t& x, Double_t& y) const;

  Bool_t CoordinatesInDetEltSt12(Int_t DeId, Double_t x, Double_t y);
  Bool_t CoordinatesInDetEltSt345(Int_t DeId, Double_t x, Double_t y);

  Int_t FromDetElt2iDet (Int_t chamber, Int_t detElt) const;
  Int_t FromDetElt2LocalId (Int_t chamber, Int_t detElt) const;


  const AliMUONGeometryTransformer* fkTransformer; //!<Geometry transformer

  AliESDEvent* fESD;             //<!Current event

  Int_t fNbrClustersCh[10];      //!<Number of clusters in the chamber [fChamberNbr].
  Int_t fTracksTotalNbr;         //!<Total number of tracks in the event.
  Int_t fTrackFilter[10];        //!<To select track for the efficiency calculation.
  Bool_t fIsCosmicData;          //!<Check if the data are cosmic rays (used only to cut cosmic shower at the trigger level if true)
  Int_t fNbrUsableTracks;        //!<Number of usable tracks (matches trigger and contains traker data, plus a trigger condition for cosmic)

  TClonesArray     * fTrackParams;  //!<Array of track param
  AliMUONTrackParam* fTrackParam;   //!<Current track param
  AliMUONVCluster  * fCluster;      //!<Current cluster

  TClonesArray* fDetEltTDHistList; //!<List of histograms of the tracks detected in the detection elements 
  TClonesArray* fDetEltTTHistList; //!<List of histograms of the tracks which have passed through the detection elements
  TClonesArray* fChamberTDHistList; //!<List of histograms of the tracks detected in the chambers 
  TClonesArray* fChamberTTHistList; //!<List of histograms of the tracks which have passed through the chambers

  static const Int_t fgkNCh;                     ///< The total number of chamber in the tracking system.
  static const Int_t fgkNSt;                     ///< The total number of station in the tracking system.
  static const Int_t fgkNDE;                     ///< Number of detection element in the tracking system.
  static const Int_t fgkNbrOfDetectionElt[10];   ///< The total number of detection element in each chamber.
  static const Int_t fgkFirstDetectionElt[10];   ///< The Id of the first detection element of each chamber.
  static const Int_t fgkOffset;                  ///< fFirstDetectionElt[iChamber] = fOffset * (iChamber+1).

  ClassDef(AliCheckMuonDetEltResponse, 0)
};
#endif
