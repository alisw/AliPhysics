// @(#) $Id$

#ifndef ALIL3HOUGHMAXFINDER_H
#define ALIL3HOUGHMAXFINDER_H

#include "AliL3RootTypes.h"
#include "AliL3StandardIncludes.h"

class AliL3Histogram;
class AliL3TrackArray;
class AliL3HoughTrack;
class TNtuple;

struct AliL3AxisWindow
{
  Int_t fYmin; // min Y
  Int_t fYmax; // max Y
  Int_t fXbin; // X bin
  Int_t fWeight; // weight
};

class AliL3HoughMaxFinder {

 public:
  AliL3HoughMaxFinder(); 
  AliL3HoughMaxFinder(Char_t *histotype,Int_t nmax,AliL3Histogram *hist=0);
  virtual ~AliL3HoughMaxFinder();
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
  void SortPeaks(struct AliL3AxisWindow **a,Int_t first,Int_t last);
  Int_t PeakCompare(struct AliL3AxisWindow *a,struct AliL3AxisWindow *b) const;
  
  //Setters:
  void SetGradient(Float_t x,Float_t y) {fGradX=x; fGradY=y;}
  void SetThreshold(Int_t f) {fThreshold = f;}
  void SetHistogram(AliL3Histogram *hist) {fCurrentHisto = hist;}
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
  
 private:

  Int_t fThreshold; // Threshold for Peak Finder
  Int_t fCurrentEtaSlice; // Current eta slice being processed
  AliL3Histogram *fCurrentHisto;  //!
  
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

  ClassDef(AliL3HoughMaxFinder,1) //Maximum finder class

};

inline Float_t AliL3HoughMaxFinder::GetXPeak(Int_t i) const
{
  if(i<0 || i>fNMax)
    {
      STDCERR<<"AliL3HoughMaxFinder::GetXPeak : Invalid index "<<i<<STDENDL;
      return 0;
    }
  return fXPeaks[i];
}

inline Float_t AliL3HoughMaxFinder::GetYPeak(Int_t i) const
{
  if(i<0 || i>fNMax)
    {
      STDCERR<<"AliL3HoughMaxFinder::GetYPeak : Invalid index "<<i<<STDENDL;
      return 0;
    }
  return fYPeaks[i];

}

inline Int_t AliL3HoughMaxFinder::GetWeight(Int_t i) const
{
  if(i<0 || i>fNMax)
    {
      STDCERR<<"AliL3HoughMaxFinder::GetWeight : Invalid index "<<i<<STDENDL;
      return 0;
    }
  return fWeight[i];
}

inline Int_t AliL3HoughMaxFinder::GetStartEta(Int_t i) const
{
  if(i<0 || i>fNMax)
    {
      STDCERR<<"AliL3HoughMaxFinder::GetStartEta : Invalid index "<<i<<STDENDL;
      return 0;
    }
  return fSTARTETAPeaks[i];
}

inline Int_t AliL3HoughMaxFinder::GetEndEta(Int_t i) const
{
  if(i<0 || i>fNMax)
    {
      STDCERR<<"AliL3HoughMaxFinder::GetStartEta : Invalid index "<<i<<STDENDL;
      return 0;
    }
  return fENDETAPeaks[i];
}

#endif

