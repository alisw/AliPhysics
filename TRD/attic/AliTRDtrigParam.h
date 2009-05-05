#ifndef ALITRDTRIGPARAM_H
#define ALITRDTRIGPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD trigger parameters class                                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"

class AliTRDtrigParam : public TObject {

 public:

  AliTRDtrigParam(const AliTRDtrigParam &p);   
  virtual         ~AliTRDtrigParam();
  AliTRDtrigParam &operator=(const AliTRDtrigParam &p); 

  static AliTRDtrigParam *Instance();
  static  void     Terminate();

  virtual void     Copy(TObject &p) const;

          void     Init();

          Int_t    GetTime1() const                        { return fTime1;         };
          Int_t    GetTime2() const                        { return fTime2;         };
          Float_t  GetClusThr() const                      { return fClusThr;       };
          Float_t  GetPadThr() const                       { return fPadThr;        };
          Int_t    GetSum10() const                        { return fSum10;         };
          Int_t    GetSum12() const                        { return fSum12;         };
          Int_t    GetTailCancelation() const              { return fTCOn;          };
          Int_t    GetNexponential() const                 { return fTCnexp;        };
          Int_t    GetFilterType() const                   { return fFilterType;    };
          void     GetFilterParam(Float_t &r1, Float_t &r2, Float_t &c1, Float_t &c2, Float_t &ped) const 
                                                           { r1            = fR1; 
                                                             r2            = fR2; 
                                                             c1            = fC1; 
                                                             c2            = fC2; 
                                                             ped           = fPedestal; };
          Float_t  GetADCnoise() const                     { return fADCnoise;      };
          Float_t  GetDeltaY() const                       { return fDeltaY;        };
          Float_t  GetDeltaS() const                       { return fDeltaS;        }; 
          Float_t  GetXprojPlane() const                   { return fXprojPlane;    };
          Float_t  GetLtuPtCut() const                     { return fLtuPtCut;      };
          Float_t  GetGtuPtCut() const                     { return fGtuPtCut;      };
          Float_t  GetHighPt() const                       { return fHighPt;        };
          Int_t    GetNPartJetLow() const                  { return fNPartJetLow;   };
          Int_t    GetNPartJetHigh() const                 { return fNPartJetHigh;  };
          Float_t  GetJetLowPt() const                     { return fJetLowPt;      };
          Float_t  GetJetHighPt() const                    { return fJetHighPt;     };

          void     SetTimeRange(Int_t time1, Int_t time2)  { fTime1        = time1; 
                                                             fTime2        = time2; };
          void     SetClusThr(Float_t clth)                { fClusThr      = clth;  };
          void     SetPadThr(Float_t path)                 { fPadThr       = path;  };
          void     SetSum10(Int_t sum)                     { fSum10        = sum;   };
          void     SetSum12(Int_t sum)                     { fSum12        = sum;   }; 
          void     SetTailCancelation(Int_t tcOn = 0)      { fTCOn         = tcOn;  };
          void     SetNexponential(Int_t nexp = 1)         { fTCnexp       = nexp;  };
          void     SetFilterType(Int_t ftype = 0)          { fFilterType   = ftype; };
          void     SetFilterParam(Float_t r1, Float_t r2, Float_t c1, Float_t c2, Float_t ped) 
                                                           { fR1           = r1; 
                                                             fR2           = r2; 
                                                             fC1           = c1; 
                                                             fC2           = c2; 
                                                             fPedestal     = ped;   };
          void     SetADCnoise(Float_t adcn)               { fADCnoise     = adcn;  };
          void     SetDeltaY(Float_t dy)                   { fDeltaY       = dy;    };
          void     SetDeltaS(Float_t ds)                   { fDeltaS       = ds;    };
          void     SetLtuPtCut(Float_t ptcut)              { fLtuPtCut     = ptcut; };
          void     SetGtuPtCut(Float_t ptcut)              { fGtuPtCut     = ptcut; };
          void     SetHighPt(Float_t hpt)                  { fHighPt       = hpt;   };
          void     SetNPartJetLow(Int_t npj)               { fNPartJetLow  = npj;   };
          void     SetNPartJetHigh(Int_t npj)              { fNPartJetHigh = npj;   };
          void     SetJetLowPt(Float_t thr)                { fJetLowPt     = thr;   };
          void     SetJetHighPt(Float_t thr)               { fJetHighPt    = thr;   };

 protected:

  static  AliTRDtrigParam *fgInstance;             // Instance of this class (singleton implementation)
  static  Bool_t           fgTerminated;           // Defines if this class has already been terminated

          Int_t    fTime1;                         // First time bin for tracking (incl.)
          Int_t    fTime2;                         // Last  time bin for tracking (incl.)
          Float_t  fClusThr;                       // Cluster threshold
          Float_t  fPadThr;                        // Pad threshold
          Int_t    fSum10;                         // MCM CreateSeeds: Min_Thr_Left_Neighbour
          Int_t    fSum12;                         // MCM CreateSeeds: Min_Sum_From_Two_Neighbours
          Int_t    fTCOn;                          // Tail cancelation flag
          Int_t    fTCnexp;                        // Number of exp in filter
          Int_t    fFilterType;                    // Filter type (0=A - analog, 1=D - digital)

          // Filter parameters (1 = long, 2 = short component)
          Float_t  fR1;                            // Time constant [microseconds]
          Float_t  fR2;                            // Time constant [microseconds]
          Float_t  fC1;                            // Weight
          Float_t  fC2;                            // Weight
          Float_t  fPedestal;                      // ADC baseline
          Float_t  fADCnoise;                      // ADC noise (not contained in the digitizer)

          Float_t  fDeltaY;                        // Y (offset) matching window in the GTU
          Float_t  fDeltaS;                        // Slope matching window in the GTU

          Float_t  fXprojPlane;                    // Projection plane (X) for GTU matching

          Float_t  fLtuPtCut;                      // Local pt cut
          Float_t  fGtuPtCut;                      // Global pt cut

          Float_t  fHighPt;                        // High pt selection

          Int_t    fNPartJetLow;                   // Number of tracks for jet (low)
          Int_t    fNPartJetHigh;                  // Number of tracks for jet (high)
          Float_t  fJetLowPt;                      // Low pt threshold for jet particles
          Float_t  fJetHighPt;                     // High pt threshold for jet particles

 private:

  // This is a singleton, constructor is private!
  AliTRDtrigParam();

  ClassDef(AliTRDtrigParam,3)                      // TRD trigger parameter class

};

#endif
