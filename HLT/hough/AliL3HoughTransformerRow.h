// @(#) $Id$

#ifndef ALIL3HOUGHTRANSFORMERROW_H
#define ALIL3HOUGHTRANSFORMERROW_H

#include "AliL3RootTypes.h"
#include "AliL3HoughBaseTransformer.h"

#define MAX_N_GAPS 5
#define MIN_TRACK_LENGTH 70

struct AliL3EtaRow {
  UChar_t fStartPad; //First pad in the cluster
  UChar_t fEndPad; //Last pad in the cluster
  Bool_t fIsFound; //Is the cluster already found
  Float_t fStartY; //Y position of the first pad in the cluster
#ifdef do_mc
  Int_t fMcLabels[MaxTrack]; //Array to store mc labels inside cluster
#endif
};

struct AliL3PadHoughParams {
  // Parameters which represent given pad in the hough space
  // Used in order to avoid as much as possible floating
  // point operations during the hough transform
  Float_t fAlpha; // Starting value for the hough parameter alpha1
  Float_t fDeltaAlpha; // Slope of alpha1
  Int_t fFirstBin; // First alpha2 bin to be filled 
  Int_t fLastBin; // Last alpha2 bin to be filled
};

class AliL3DigitData;
class AliL3Histogram;

class AliL3HoughTransformerRow : public AliL3HoughBaseTransformer {

 public:
  AliL3HoughTransformerRow(); 
  AliL3HoughTransformerRow(Int_t slice,Int_t patch,Int_t netasegments,Bool_t DoMC=kFALSE,Float_t zvertex=0.0);
  virtual ~AliL3HoughTransformerRow();

  //void CreateHistograms(Float_t ptmin,Float_t ptmax,Float_t ptres,Int_t nybin,Float_t psi);
  void CreateHistograms(Int_t /*nxbin*/,Float_t /*ptmin*/,Int_t /*nybin*/,Float_t /*phimin*/,Float_t /*phimax*/)
  {STDCERR<<"This method for creation of parameter space histograms is not supported for this Transformer!"<<STDENDL;}
  void CreateHistograms(Int_t nxbin,Float_t xmin,Float_t xmax,
			Int_t nybin,Float_t ymin,Float_t ymax);
  void Reset();
  void TransformCircle();

  Int_t GetEtaIndex(Double_t eta) const;
  AliL3Histogram *GetHistogram(Int_t etaindex);
  Double_t GetEta(Int_t etaindex,Int_t slice) const;
  Int_t GetTrackID(Int_t etaindex,Double_t alpha1,Double_t alpha2) const;
  UChar_t *GetGapCount(Int_t etaindex) const { return fGapCount[etaindex]; }
  UChar_t *GetCurrentRowCount(Int_t etaindex) const { return fCurrentRowCount[etaindex]; }
  UChar_t *GetPrevBin(Int_t etaindex) const { return fPrevBin[etaindex]; }
  UChar_t *GetNextBin(Int_t etaindex) const { return fNextBin[etaindex]; }
  UChar_t *GetNextRow(Int_t etaindex) const { return fNextRow[etaindex]; }
  UChar_t *GetTrackNRows() const { return fTrackNRows; }
  UChar_t *GetTrackFirstRow() const { return fTrackFirstRow; }
  UChar_t *GetTrackLastRow() const { return fTrackLastRow; }
  static Float_t GetBeta1() {return fgBeta1;}
  static Float_t GetBeta2() {return fgBeta2;}

  void SetTPCRawStream(AliTPCRawStream *rawstream) {fTPCRawStream=rawstream;}

 private:

  UChar_t **fGapCount; //!
  UChar_t **fCurrentRowCount; //!
#ifdef do_mc
  AliL3TrackIndex **fTrackID; //!
#endif

  UChar_t *fTrackNRows; //!
  UChar_t *fTrackFirstRow; //!
  UChar_t *fTrackLastRow; //!
  UChar_t *fInitialGapCount; //!

  UChar_t **fPrevBin; //!
  UChar_t **fNextBin; //!
  UChar_t **fNextRow; //!

  AliL3PadHoughParams **fStartPadParams; //!
  AliL3PadHoughParams **fEndPadParams; //!
  Float_t **fLUTr2; //!

  Float_t *fLUTforwardZ; //!
  Float_t *fLUTforwardZ2; //!
  Float_t *fLUTbackwardZ; //!
  Float_t *fLUTbackwardZ2; //!

  AliL3Histogram **fParamSpace; //!

  void TransformCircleFromDigitArray();
  void TransformCircleFromRawStream();

  void DeleteHistograms(); //Method to clean up the histograms containing Hough space

  inline void FillClusterRow(UChar_t i,Int_t binx1,Int_t binx2,UChar_t *ngaps2,UChar_t *currentrow2,UChar_t *lastrow2
#ifdef do_mc
			     ,AliL3EtaRow etaclust,AliL3TrackIndex *trackid
#endif
			     );
  inline void FillCluster(UChar_t i,Int_t etaindex,AliL3EtaRow *etaclust,Int_t ilastpatch,Int_t firstbinx,Int_t lastbinx,Int_t nbinx,Int_t firstbiny);
#ifdef do_mc
  inline void FillClusterMCLabels(AliL3DigitData digpt,AliL3EtaRow *etaclust);
#endif

  void SetTransformerArrays(AliL3HoughTransformerRow *tr);

  static Float_t fgBeta1,fgBeta2; // Two curves which define the Hough space

  AliTPCRawStream *fTPCRawStream; // Pointer to the raw stream in case of fast reading of the raw data (fast_raw flag)

  ClassDef(AliL3HoughTransformerRow,1) //TPC Rows Hough transformation class

};

#endif




