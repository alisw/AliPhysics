#ifndef ALICFEVENTRECCUTS_H
#define ALICFEVENTRECCUTS_H
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
// Cut on the Event at reconstructed level: for the moment 
// just the requirements on the number of charged tracks and on 
// the vertex position and resolution are implemented
// The argument of IsSelected member function (passed object) is cast into 
// an AliESDEvent. In the future may be modified to use AliVEvent interface
// and include more cut variables.
// The class derives from AliCFCutBase
// Author:S.Arcelli Silvia.Arcelli@cern.ch


#include "AliCFCutBase.h"
class TH1F;
class TBits;
//_____________________________________________________________________________
class AliCFEventRecCuts: public AliCFCutBase 
{
 public :
  AliCFEventRecCuts() ;
  AliCFEventRecCuts(Char_t* name, Char_t* title) ;
  AliCFEventRecCuts(const AliCFEventRecCuts& c) ;
  AliCFEventRecCuts& operator=(const AliCFEventRecCuts& c) ;
  ~AliCFEventRecCuts();
  void GetBitMap(TObject *obj, TBits*bitmap);
  Bool_t IsSelected(TObject* obj);
  void Init();
  void AddQAHistograms(TList *list) const;

  void SetNTracksCut(Int_t xMin=-1, Int_t xMax=1000000) {fNTracksMin=xMin; fNTracksMax=xMax;} // cut values setter

  void SetRequireVtxCuts(Bool_t vtx=kFALSE) {fRequireVtxCuts=vtx;} // cut values setter
  void SetVertexXCut(Double_t xMin=-1.e99, Double_t xMax=1.e99) { fVtxXMin=xMin; fVtxXMax=xMax;} // cut values setter
  void SetVertexYCut(Double_t yMin=-1.e99, Double_t yMax=1.e99) { fVtxYMin=yMin; fVtxYMax=yMax;} // cut values setter
  void SetVertexZCut(Double_t zMin=-1.e99, Double_t zMax=1.e99) { fVtxZMin=zMin; fVtxZMax=zMax;} // cut values setter

  void SetVertexXResCut(Double_t xMax=1.e99) {fVtxXResMax=xMax;} // cut values setter
  void SetVertexYResCut(Double_t yMax=1.e99){fVtxYResMax=yMax;} // cut values setter
  void SetVertexZResCut(Double_t zMax=1.e99){fVtxZResMax=zMax;} // cut values setter

  Int_t    GetNTracksMin() const {return fNTracksMin;} // cut values getter
  Int_t    GetNTracksMax() const {return fNTracksMax;} // cut values getter
  Bool_t   GetRequireVtxCuts() const {return fRequireVtxCuts;} // cut value getter
  Double_t GetVertexXMax() const {return fVtxXMax;} // cut values getter
  Double_t GetVertexYMax() const {return fVtxYMax;} // cut values getter
  Double_t GetVertexZMax() const {return fVtxZMax;} // cut values getter
  Double_t GetVertexXMin() const {return fVtxXMin;} // cut values getter
  Double_t GetVertexYMin() const {return fVtxYMin;} // cut values getter
  Double_t GetVertexZMin() const {return fVtxZMin;} // cut values getter
  Double_t GetVertexXResMax() const {return fVtxXResMax;} // cut values getter
  Double_t GetVertexYResMax() const {return fVtxYResMax;} // cut values getter
  Double_t GetVertexZResMax() const {return fVtxZResMax;} // cut values getter

  // QA histograms
  void FillHistogramsBeforeCuts(TObject* obj) {return FillHistograms(obj,kFALSE);}
  void FillHistogramsAfterCuts(TObject* obj)  {return FillHistograms(obj,kTRUE);}
  // QA histogram setter
  // please use indices from the enumeration below
  void SetHistogramBins(Int_t index, Int_t nbins, Double_t *bins);
  void SetHistogramBins(Int_t index, Int_t nbins, Double_t xmin, Double_t xmax);
  enum{kNTracks=0,
	 kVtxPosX,
	 kVtxPosY,
	 kVtxPosZ,
	 kVtxResX,
	 kVtxResY,
	 kVtxResZ,
	 kNCuts=7,
	 kNStepQA=2
	 };
  
 private:
  TBits *SelectionBitMap(TObject* obj);
  void DefineHistograms(); 		// books histograms 
  void Initialise();			// sets everything to 0
  void FillHistograms(TObject* obj, Bool_t b);
  Int_t fNTracksMin; //minimum number of esd tracks
  Int_t fNTracksMax; //maximum number of esd tracks
  Bool_t fRequireVtxCuts ; //The type of trigger to be checked
  Double_t fVtxXMax ; //X vertex position, maximum value
  Double_t fVtxYMax ; //Y vertex position, maximum value 
  Double_t fVtxZMax ; //Z vertex position, maximum value
  Double_t fVtxXMin ; //X vertex position, minimum value
  Double_t fVtxYMin ; //Y vertex position, minimum value
  Double_t fVtxZMin ; //Z vertex position, minimum value
  Double_t fVtxXResMax ;//Maximum value of sigma_vtx in X
  Double_t fVtxYResMax ;//Maximum value of sigma_vtx in X
  Double_t fVtxZResMax ;//Maximum value of sigma_vtx in X


  TBits *fBitMap ; //cut mask


  TH1F* fhQA[kNCuts][kNStepQA];		// QA Histograms
  //QA Histogram parameters
  Int_t fhNBinsNTracks;//size of array of bin limits, N Tracks ESD
  Double_t *fhBinLimNTracks;//[fhNBinsNTracks] bin limits, N Tracks ESD
  Int_t fhNBinsVtxPosX;//size of array of bin limits, Vtx Pos X
  Double_t *fhBinLimVtxPosX;//[fhNBinsVtxPosX] bin limits, Vtx Pos X
  Int_t fhNBinsVtxPosY;//size of array of bin limits, Vtx Pos Y
  Double_t *fhBinLimVtxPosY;//[fhNBinsVtxPosY] bin limits, Vtx Pos Y
  Int_t fhNBinsVtxPosZ;//size of array of bin limits, Vtx Pos Z
  Double_t *fhBinLimVtxPosZ;//[fhNBinsVtxPosZ] bin limits, Vtx Pos Z
  Int_t fhNBinsVtxResX;//size of array of bin limits, Vtx Res X
  Double_t *fhBinLimVtxResX;//[fhNBinsVtxResX] bin limits, Vtx Res X
  Int_t fhNBinsVtxResY;//size of array of bin limits, Vtx Res Y
  Double_t *fhBinLimVtxResY;//[fhNBinsVtxResY] bin limits, Vtx Res Y
  Int_t fhNBinsVtxResZ;//size of array of bin limits, Vtx Res Z
  Double_t *fhBinLimVtxResZ;//[fhNBinsVtxResZ] bin limits, Vtx Res Z

  ClassDef(AliCFEventRecCuts,2);
};

#endif
