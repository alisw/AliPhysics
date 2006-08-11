#ifndef ALITRDMCMTRACKLET_H
#define ALITRDMCMTRACKLET_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////
//                                                   //
//  Tracklet object (MCM/LTU)                        //
//                                                   //
///////////////////////////////////////////////////////

#include <TObject.h>

class TGraph;

class AliTRDgeometry;

class AliTRDmcmTracklet : public TObject {

 public:

  enum { kNclsPads = 3, kNtimeBins = 30, kNdict = 3 };

  AliTRDmcmTracklet();
  AliTRDmcmTracklet(Int_t det, Int_t row, Int_t n);
  AliTRDmcmTracklet(const AliTRDmcmTracklet &t);
  virtual           ~AliTRDmcmTracklet();
  AliTRDmcmTracklet &operator=(const AliTRDmcmTracklet &t);
  virtual void       Copy(TObject &t) const;

  void      Reset();

  void      AddCluster(Int_t icol, Int_t itb, Float_t *adc, Int_t *track);
  void      MakeTrackletGraph(AliTRDgeometry *geo = 0, Float_t field = 0);
  void      MakeClusAmpGraph();
  void      CookLabel(Float_t frac);

  void      SetRow(Int_t row)               { fRow      = row;    };
  void      SetDetector(Int_t det)          { fDetector = det;    };
  void      SetN(Int_t n)                   { fN        = n;      };

  Int_t     GetNclusters() const            { return fNclusters;  };
  Int_t     GetDetector()  const            { return fDetector;   };
  Int_t     GetRow()       const            { return fRow;        };
  Float_t   GetOffset()    const            { return fOffset;     };
  Float_t   GetSlope()     const            { return fSlope;      };
  Float_t   GetTime0()     const            { return fTime0;      };
  Float_t   GetRowz()      const            { return fRowz;       };
  Float_t   GetPt()        const            { return fPt;         };
  Float_t   GetdQdl()      const            { return fdQdl;       };
  Int_t     GetLabel()     const            { return fTrackLabel; };
  Int_t     GetNumber()    const            { return fN;          };
  Float_t   GetOmegaTau(Float_t vdrift, Float_t field) const;
  Float_t   GetClusY(Float_t *adc, Int_t pla) const;
  TGraph   *GetTrackletGraph() const        { return fGPos;       };
  TGraph   *GetClusAmpGraph()  const        { return fGAmp;       };
  Float_t  *GetClusterADC(Int_t icl)        { return fADC[icl];   };
  Int_t     GetClusterTime(Int_t icl) const { return fTime[icl];  };
  Int_t     GetClusterCol(Int_t icl)  const { return fCol[icl];   };

 protected:

  Int_t   fDetector;                     //  TRD detector number (0 ... 539)
  Int_t   fRow;                          //  Row number in the detector
  Float_t fADC[kNtimeBins][kNclsPads];   //  Array of ADC values in a pad group
  Int_t   fTrack[kNtimeBins][kNdict];    //! Array of track dictionary values
  Int_t   fTrackLabel;                   //  Cooked track label
  Int_t   fTime[kNtimeBins];             //  Array of time bin values
  Int_t   fCol[kNtimeBins];              //  Array of pad column values
  Int_t   fNclusters;                    //  Number of clusters in the tracklet

  Int_t   fN;                            //  Tracklet number

  TGraph *fGPos;                         //! Positions
  TGraph *fGAmp;                         //! Amplitudes

  Float_t fTime0;                        //  X position at the entrance window
  Float_t fRowz;                         //  Z position of the row center
  Float_t fSlope;                        //  Slope [deg]
  Float_t fOffset;                       //  Offset
  Float_t fPt;                           //  Transverse momentum
  Float_t fdQdl;                         //  Charge per unit length

  ClassDef(AliTRDmcmTracklet,2)          //  Track segment for the TRD (Tracklet)

};

#endif
