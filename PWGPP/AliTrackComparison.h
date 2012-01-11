#ifndef ALITRACKCOMPARISON_H
#define ALITRACKCOMPARISON_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//-------------------------------------------------------------------
// Base class to test the extrapolation performance from TPC to outer 
// detectors. Several member functions AddTracks() with different
// arguments can be called to execute extrapolation and several 
// THnSparse are filled with residuls and basic track information
//
// Anthor: R.Ma, M.Ivanov 02/04/2011
//--------------------------------------------------------------------

#include "TNamed.h"
#include "TMatrixDfwd.h"
class THnSparse;
class TCollection;
class AliExternalTrackParam;
class AliTrackPointArray;
class AliTrackPoint;
class TParticle;
class AliTrackReference;
class TObjArray;
class TTreeSRedirector;
 
class AliTrackComparison:public TNamed {
public:
  AliTrackComparison();
  AliTrackComparison(const Text_t *name, const Text_t *title);
  AliTrackComparison(const AliTrackComparison& comp);
  AliTrackComparison& operator=(const AliTrackComparison& comp);
  virtual ~AliTrackComparison();
  //Make histograms
  void Init();

  //Main process functions  
  Int_t AddTracks(AliExternalTrackParam *param0,  AliExternalTrackParam *param1, Double_t mass);
  Int_t AddTracks(AliTrackReference *ref0,  AliTrackReference *ref1, Double_t mass, Int_t charge);
  Int_t AddTracks(AliExternalTrackParam *param0,  AliTrackPoint *point1, Double_t mass, Double_t energy, Double_t *vxyz);
  AliExternalTrackParam *MakeTrack(const AliTrackReference *ref, Int_t charge);
  Bool_t PropagateToPoint(AliExternalTrackParam *param0,AliExternalTrackParam *param1, Double_t mass);
  
  //Set the step used in extrapolatoin
  void SetStep(Double_t step){fStep=step;}

  //Set the range of the first bins of the THnSparse
  void SetRangeDY(Double_t lowBinDY, Double_t upBinDY) {fLowBinDY=lowBinDY,fUpBinDY=upBinDY;}
  void SetRangeDZ(Double_t lowBinDZ, Double_t upBinDZ) {fLowBinDZ=lowBinDZ,fUpBinDZ=upBinDZ;}
  void SetRangeDSnp(Double_t lowBinDSnp, Double_t upBinDSnp){fLowBinDSnp=lowBinDSnp;fUpBinDSnp=upBinDSnp;}
  void SetRangeDTheta(Double_t lowBinDTheta, Double_t upBinDTheta) {fLowBinDTheta=lowBinDTheta;fUpBinDTheta=upBinDTheta;}
  void SetRange1Pt(Double_t lowBin1Pt, Double_t upBin1Pt) {fLowBin1Pt=lowBin1Pt;fUpBin1Pt=upBin1Pt;}
  void SetRange1PtLoss(Double_t lowBin1PtLoss, Double_t upBin1PtLoss) {fLowBin1PtLoss=lowBin1PtLoss;fUpBin1PtLoss=upBin1PtLoss;}


  //Set number of bins 
  void SetNBins(Int_t nBinsDY, Int_t nBinsDZ, Int_t nBinsDSnp, Int_t nBinsDTheta, Int_t nBins1Pt, Int_t nBins1PtLoss)
  {
    fNBinsDY=nBinsDY;
    fNBinsDZ=nBinsDZ;
    fNBinsDSnp=nBinsDSnp;
    fNBinsDTheta=nBinsDTheta;
    fNBins1Pt=nBins1Pt;
    fNBins1PtLoss=nBins1PtLoss;
  }


  //Set the layer ID of the detectors. Refer to AliGeomManager.h
  void SetLayerID(Double_t layerID) {fLayerID=layerID;}

  //Set the flag to fill all the THnSparse's. By default it is kTURE
  void SetFillAll(Bool_t fillAll) {fFillAll=fillAll;}
  
  //Fill the THnSparse
  void FillHistos(AliExternalTrackParam *param0, AliExternalTrackParam *param1, Double_t tr1Pt);

  //Getter
  THnSparse *GetHnSparse(Int_t index) {return fResolHisto[index];}

  //Used to redirect the output THnSparse into TTreeSRedirector
  void MakeDistortionMap(Double_t refX, Int_t type);

  //Set number of bins used to combine statistics
  void SetNCombineBin(Int_t nCombineBin) {fNCombineBin=nCombineBin;}

  virtual Long64_t       Merge(TCollection *const li);
  virtual void           Add(AliTrackComparison *const comp);
  void Analyze(); //Not implemented

protected:
  void    MakeHistos();     //Initialize all the THnSparse
  Double_t fStep;           //Step used in extrapolation
  Double_t fLowBinDY;       //First bin of residual in Y
  Double_t fUpBinDY;        //Last bin of residual in Y
  Double_t fLowBinDZ;       //First bin of residual in Z
  Double_t fUpBinDZ;        //Last bin of residual in Z
  Double_t fLowBinDSnp;     //First bin of residual in sin(phi)
  Double_t fUpBinDSnp;      //Last bin of residual in sin(phi) 
  Double_t fLowBinDTheta;   //First bin of residual in Theta
  Double_t fUpBinDTheta;    //Last bin of residual in Theta
  Double_t fLowBin1Pt;      //First bin of residual in 1/pT
  Double_t fUpBin1Pt;       //Last bin of residual in 1/pT
  Double_t fLowBin1PtLoss;  //First bin of corrected energy loss
  Double_t fUpBin1PtLoss;   //Last bin of corrected energy loss
  Int_t fNBinsDY;           //NBins of residual in Y
  Int_t fNBinsDZ;           //NBins of residual in Z
  Int_t fNBinsDSnp;         //NBins of residual in sin(phi)
  Int_t fNBinsDTheta;       //NBins of residual in Theta
  Int_t fNBins1Pt;          //NBins of residual in 1/pT
  Int_t fNBins1PtLoss;      //NBins of corrected energy loss
  Double_t fLayerID;        //Detector layer ID
  Bool_t fFillAll;          //Flag to fill all the THnSparse
  Int_t fNCombineBin;       //Number of bins used to combine statistics in MakeDistortionMap()

  THnSparse *fResolHisto[6];//Output THnSparse

  ClassDef(AliTrackComparison, 4); 
};

#endif


