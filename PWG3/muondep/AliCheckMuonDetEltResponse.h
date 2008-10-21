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
			     TClonesArray* detEltTTHistList);

//Destructor:
  virtual ~AliCheckMuonDetEltResponse();


  void CheckDetEltResponse ();
  void TrackLoop ();
  void TrackParamLoop ();
 
private:
  
  void FillDetEltTDHisto (Int_t chamber, Int_t detElt,
			  Double_t posXL, Double_t posYL);

  void FillDetEltTTHisto (Int_t chamber, Int_t detElt,
			  Double_t posXG, Double_t posYG, Double_t posZG,
			  Double_t posXL, Double_t posYL, Double_t posZL);

  void CalculMissClusterParam (AliMUONTrackParam* extrapTrackParam,
			       Int_t firstMissCh, Int_t nbrOfMissCh);

  void GetDetEltFromPosition (Int_t chamber,
			      Double_t posX, Double_t posY, Double_t posZ);


  Int_t fNCh; //!<Number of tracking chamber.
  Int_t fNSt; //!<Number of tracking station.
  Int_t fNDE; //!<Number of detection element in the tracking system.

  Int_t FromDetElt2iDet (Int_t chamber, Int_t detElt);

  const AliMUONGeometryTransformer* fTransformer; //!<Geometry transformer

  AliESDEvent* fESD;

  Int_t fNbrClustersCh[10];      //!<Number of clusters in the chamber [fChamberNbr].
  Int_t fTracksTotalNbr;         //!<Total number of tracks in the event.
  Int_t fGetDetElt[2];           //!<Detection elemtents obtained from position(X,Y,Z). fGetDetElt[1] is for the case where there is an overlap.
  Int_t fTrackFilter[10];        //!<To select track for the efficiency calculation.

  TClonesArray     * fTrackParams;
  AliMUONTrackParam* fTrackParam;
  AliMUONVCluster  * fCluster;

  TClonesArray* fDetEltTDHistList; //!<List of histograms of the tracks detected in the detection elements 
  TClonesArray* fDetEltTTHistList; //!<List of histograms of the tracks which have passed through the detection elements

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
