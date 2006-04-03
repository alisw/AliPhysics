#ifndef ALITRDMCMTRACKLET_H
#define ALITRDMCMTRACKLET_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////
//  Tracklet object (MCM/TRIGGER)                    //
///////////////////////////////////////////////////////

#include <TObject.h>

const Int_t kNclsPads  =  3;
const Int_t kNtimeBins = 30;
const Int_t kNdict     =  3;

class TGraph;

class AliTRDgeometry;

class AliTRDmcmTracklet : public TObject {

 public:

  AliTRDmcmTracklet();
  AliTRDmcmTracklet(Int_t det, Int_t row, Int_t n);
  virtual ~AliTRDmcmTracklet();

  void AddCluster(Int_t icol, Int_t itb, Float_t *adc, Int_t *track) { 
    if (fNclusters >= kNtimeBins) return;
    for (Int_t icl = 0; icl < kNclsPads; icl++) {
      //fADC[fNclusters][icl] = (Int_t)adc[icl]; 
      fADC[fNclusters][icl] = adc[icl]; 
    }
    fTrack[fNclusters][0] = track[0];
    fTrack[fNclusters][1] = track[1];
    fTrack[fNclusters][2] = track[2];
    fTime[fNclusters] = itb;
    fCol[fNclusters] = icol;
    fNclusters++;
  };

  Float_t  *GetClusterADC(Int_t icl)      { return fADC[icl];  };
  Int_t     GetClusterTime(Int_t icl)     { return fTime[icl]; };
  Int_t     GetClusterCol(Int_t icl)      { return fCol[icl];  };

  TGraph   *GetTrackletGraph() { return fGPos; };
  TGraph   *GetClusAmpGraph()  { return fGAmp; };

  virtual void MakeTrackletGraph(AliTRDgeometry *geo = 0, Float_t field = 0);
  virtual void MakeClusAmpGraph();

  Float_t   GetClusY(Float_t *adc, Int_t pla);

  Int_t     GetNumber() { return fN; };

  void      CookLabel(Float_t frac);
  Int_t     GetLabel() { return fTrackLabel; };

  Int_t     GetNclusters()    { return fNclusters; };
  Int_t     GetDetector()     { return fDetector; };
  Int_t     GetRow()          { return fRow; };
  Float_t   GetOffset() { return fOffset; };
  Float_t   GetSlope()  { return fSlope; };
  Float_t   GetTime0()  { return fTime0; };
  Float_t   GetRowz()   { return fRowz; };
  Float_t   GetPt()     { return fPt; };
  Float_t   GetdQdl()   { return fdQdl; };

  Float_t   GetOmegaTau(Float_t vdrift, Float_t field);

 protected:

  Int_t   fDetector;                     //  TRD detector number (0 ... 539)
  Int_t   fRow;                          //  Row number in the detector
  Float_t fADC[kNtimeBins][kNclsPads];   //  Array of ADC values in a pad group
  Int_t   fTrack[kNtimeBins][kNdict];    //! Array of track dictionary values
  Int_t   fTrackLabel;                   //  Cooked track label
  Int_t   fTime[kNtimeBins];             //  Array of time bin values
  Int_t   fCol[kNtimeBins];              //  Array of pad column values
  Int_t   fNclusters;                    //  Number of clusters in the tracklet

  Int_t   fN;                            // Tracklet number

  TGraph *fGPos;                         //! Positions
  TGraph *fGAmp;                         //! Amplitudes

  Float_t fTime0;                        // X position at the entrance window
  Float_t fRowz;                         // Z position of the row center
  Float_t fSlope;                        // Slope [deg]
  Float_t fOffset;                       // Offset
  Float_t fPt;                           // Transverse momentum
  Float_t fdQdl;                         // Charge per unit length

  ClassDef(AliTRDmcmTracklet,1)          // Track segment for the TRD (Tracklet)

};

#endif
