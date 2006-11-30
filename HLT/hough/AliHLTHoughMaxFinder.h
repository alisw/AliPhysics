// @(#) $Id$

#ifndef ALIL3HOUGHMAXFINDER_H
#define ALIL3HOUGHMAXFINDER_H

#include "AliHLTRootTypes.h"
#include "AliHLTStandardIncludes.h"

class AliHLTHistogram;
class AliHLTTrackArray;
class AliHLTHoughTrack;
class TNtuple;

struct AliHLTAxisWindow
{
  Int_t fYmin; // min Y
  Int_t fYmax; // max Y
  Int_t fXbin; // X bin
  Int_t fWeight; // weight
};

struct AliHLTPre2DPeak
{
  Float_t fX; // X coordinate of the preak
  Float_t fY; // Y coordinate of the preak
  Float_t fSizeX; // Size of the peak
  Float_t fSizeY; // Size of the peak
  Int_t fStartX; // Start position of the peak
  Int_t fStartY; // Start position of the peak
  Int_t fEndX; // End position of the peak
  Int_t fEndY; // End position of the peak
  Float_t fWeight; // Weight assigned to the peak
};

class AliHLTHoughMaxFinder {

 public:
  AliHLTHoughMaxFinder(); 
  AliHLTHoughMaxFinder(Char_t *histotype,Int_t nmax,AliHLTHistogram *hist=0);
  virtual ~AliHLTHoughMaxFinder();
  void Reset();

  void CreateNtuppel();
  void WriteNtuppel(Char_t *filename);

  //Simple maxima finders:
  void FindAbsMaxima();
  void FindBigMaxima();
  void FindMaxima(Int_t threshold=0);
  void FindAdaptedPeaks(Int_t nkappawindow,Float_t cutratio);
  //Peak finder for HoughTransformerRow
  void FindAdaptedRowPeaks(Int_t kappawindow,Int_t xsize,Int_t ysize);
  //More sophisticated peak finders:
  void FindPeak(Int_t t1,Double_t t2,Int_t t3);
  void FindPeak1(Int_t ywindow=2,Int_t xbinsides=1);
  void SortPeaks(struct AliHLTAxisWindow **a,Int_t first,Int_t last);
  Int_t PeakCompare(struct AliHLTAxisWindow *a,struct AliHLTAxisWindow *b) const;
  
  //Setters:
  void SetGradient(Float_t x,Float_t y) {fGradX=x; fGradY=y;}
  void SetThreshold(Int_t f) {fThreshold = f;}
  void SetHistogram(AliHLTHistogram *hist) {fCurrentHisto = hist;}
  void SetTrackLUTs(UChar_t *tracknrows, UChar_t *trackfirstrow, UChar_t *tracklastrow, UChar_t *nextrow) {fTrackNRows = tracknrows; fTrackFirstRow = trackfirstrow; fTrackLastRow = tracklastrow; fNextRow = nextrow;}
  void SetEtaSlice(Int_t etaslice) {fCurrentEtaSlice = etaslice;}
  
  //Getters:
  Float_t GetXPeak(Int_t i) const;
  Float_t GetYPeak(Int_t i) const;
  Float_t GetXPeakSize(Int_t i) const;
  Float_t GetYPeakSize(Int_t i) const;
  Int_t GetWeight(Int_t i) const;
  Int_t GetStartEta(Int_t i) const;
  Int_t GetEndEta(Int_t i) const;
  Int_t GetEntries() const {return fNPeaks;}

  //Method for merging of peaks produced by AliHLTHoughTransfromerRow
  Bool_t MergeRowPeaks(AliHLTPre2DPeak *maxima1, AliHLTPre2DPeak *maxima2,Float_t distance);
  
 private:

  Int_t fThreshold; // Threshold for Peak Finder
  Int_t fCurrentEtaSlice; // Current eta slice being processed
  AliHLTHistogram *fCurrentHisto;  //!

  UChar_t *fTrackNRows; //!
  UChar_t *fTrackFirstRow; //!
  UChar_t *fTrackLastRow; //!
  UChar_t *fNextRow; //!
  
  Float_t fGradX; // Gradient threshold inside Peak Finder 
  Float_t fGradY; // Gradient threshold inside Peak Finder 
  Float_t *fXPeaks; //!
  Float_t *fYPeaks; //!
  Int_t *fSTARTXPeaks; //!
  Int_t *fSTARTYPeaks; //!
  Int_t *fENDXPeaks; //!
  Int_t *fENDYPeaks; //!
  Int_t *fSTARTETAPeaks; //!
  Int_t *fENDETAPeaks; //!
  Int_t *fWeight;   //!
  Int_t fN1PeaksPrevEtaSlice; // Index of the first peak in the previous eta slice
  Int_t fN2PeaksPrevEtaSlice; // Index of the  last peak in the previous eta slice
  Int_t fNPeaks; // Index of the last accumulated peak
  Int_t fNMax; // Maximum allowed number of peaks
  
  Char_t fHistoType; // Histogram type

#ifndef no_root
  TNtuple *fNtuppel; //!
#endif

  ClassDef(AliHLTHoughMaxFinder,1) //Maximum finder class

};

typedef AliHLTHoughMaxFinder AliL3HoughMaxFinder; // for backward compatibility

inline Float_t AliHLTHoughMaxFinder::GetXPeak(Int_t i) const
{
  if(i<0 || i>fNMax)
    {
      STDCERR<<"AliHLTHoughMaxFinder::GetXPeak : Invalid index "<<i<<STDENDL;
      return 0;
    }
  return fXPeaks[i];
}

inline Float_t AliHLTHoughMaxFinder::GetYPeak(Int_t i) const
{
  if(i<0 || i>fNMax)
    {
      STDCERR<<"AliHLTHoughMaxFinder::GetYPeak : Invalid index "<<i<<STDENDL;
      return 0;
    }
  return fYPeaks[i];

}

inline Int_t AliHLTHoughMaxFinder::GetWeight(Int_t i) const
{
  if(i<0 || i>fNMax)
    {
      STDCERR<<"AliHLTHoughMaxFinder::GetWeight : Invalid index "<<i<<STDENDL;
      return 0;
    }
  return fWeight[i];
}

inline Int_t AliHLTHoughMaxFinder::GetStartEta(Int_t i) const
{
  if(i<0 || i>fNMax)
    {
      STDCERR<<"AliHLTHoughMaxFinder::GetStartEta : Invalid index "<<i<<STDENDL;
      return 0;
    }
  return fSTARTETAPeaks[i];
}

inline Int_t AliHLTHoughMaxFinder::GetEndEta(Int_t i) const
{
  if(i<0 || i>fNMax)
    {
      STDCERR<<"AliHLTHoughMaxFinder::GetStartEta : Invalid index "<<i<<STDENDL;
      return 0;
    }
  return fENDETAPeaks[i];
}

#endif

