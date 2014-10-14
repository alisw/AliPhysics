#ifndef ALITRDCALDCSv2_H
#define ALITRDCALDCSv2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalDCSv2.h 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class v2 for TRD DCS parameters                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"
#include "TObjArray.h"

class TString;

class AliTRDCalDCSFEEv2;
class AliTRDCalDCSPTR;
class AliTRDCalDCSGTU;

class AliTRDCalDCSv2 : public TNamed {

 public:

  AliTRDCalDCSv2();
  AliTRDCalDCSv2(const Text_t *name, const Text_t *title);
  AliTRDCalDCSv2(const AliTRDCalDCSv2 &cd);
  AliTRDCalDCSv2 &operator=(const AliTRDCalDCSv2 &cd);
  virtual ~AliTRDCalDCSv2() { };

  void    EvaluateGlobalParameters();
  void    SetFEEArr(TObjArray * const fa)      { fFEEArr              = fa;    }
  void    SetPTRArr(TObjArray * const pa)      { fPTRArr              = pa;    }
  void    SetGTUObj(AliTRDCalDCSGTU *go)       { fGTUObj              = go;    }
  void    SetRunType(TString rt)               { fRunType = rt;                }
  void    SetStartTime(UInt_t st)              { fStartTime = st;              }
  void    SetEndTime(UInt_t et)                { fEndTime = et;                }
  
  Int_t   GetGlobalNumberOfTimeBins() const    { return fGNumberOfTimeBins;    }
  Int_t   GetGlobalConfigTag() const           { return fGConfigTag;           }
  Int_t   GetGlobalSingleHitThres() const      { return fGSingleHitThres;      }
  Int_t   GetGlobalThreePadClustThres() const  { return fGThreePadClustThres;  }
  Int_t   GetGlobalSelectiveNoZS() const       { return fGSelNoZS;             }
  Int_t   GetGlobalTCFilterWeight() const      { return fGTCFilterWeight;      }
  Int_t   GetGlobalTCFilterShortDecPar() const { return fGTCFilterShortDecPar; }
  Int_t   GetGlobalTCFilterLongDecPar() const  { return fGTCFilterLongDecPar;  }
  Int_t   GetGlobalModeFastStatNoise() const   { return fGFastStatNoise;       }
  TString GetGlobalConfigVersion() const       { return fGConfigVersion;       }
  TString GetGlobalConfigName() const          { return fGConfigName;          }
  TString GetGlobalFilterType() const          { return fGFilterType;          }
  TString GetGlobalReadoutParam() const        { return fGReadoutParam;        }
  TString GetGlobalTestPattern() const         { return fGTestPattern;         }
  TString GetGlobalTrackletMode() const        { return fGTrackletMode;        }
  TString GetGlobalTrackletDef() const         { return fGTrackletDef;         }
  TString GetGlobalTriggerSetup() const        { return fGTriggerSetup;        }
  TString GetGlobalAddOptions() const          { return fGAddOptions;          }
  TString GetRunType() const                   { return fRunType;              }
  UInt_t  GetStartTime() const                 { return fStartTime;            }
  UInt_t  GetEndTime() const                   { return fEndTime;              }
  TObjArray*       GetFEEArr() const           { return fFEEArr;               }
  TObjArray*       GetPTRArr() const           { return fPTRArr;               }
  AliTRDCalDCSFEEv2* GetCalDCSFEEObj(Int_t det) const
  		  	          { return (AliTRDCalDCSFEEv2*)fFEEArr->At(det); }
  AliTRDCalDCSPTR* GetCalDCSPTRObj(Int_t det) const
  			          { return (AliTRDCalDCSPTR*)fPTRArr->At(det); }
  AliTRDCalDCSGTU* GetGTUObj() const
           		          { return (AliTRDCalDCSGTU*)fGTUObj;          }

 protected:

  // global configuration parameters
  Int_t   fGNumberOfTimeBins;    // Number of timebins (-1 if diverse)
  Int_t   fGConfigTag;           // Configuration Tag (-1 if diverse)
  Int_t   fGSingleHitThres;      // thres. of single hits (arg of readout param) (-1 if diverse)
  Int_t   fGThreePadClustThres;  // thres. of 3-pad clusters (arg of readout param) (-1 if diverse)
  Int_t   fGSelNoZS;             // write every fGSelNoZS'th event without ZS (-1 if diverse)
  Int_t   fGTCFilterWeight;      // tail cancellation filter weight (-1 if diverse)
  Int_t   fGTCFilterShortDecPar; // tail cancellation filter short decay parameter (-1 if diverse)
  Int_t   fGTCFilterLongDecPar;  // tail cancellation filter long decay parameter (-1 if diverse)
  Int_t   fGFastStatNoise;       // collect stat. f. fast noise mode (0: no, 1: yes, -1: diverse)
  TString fGConfigVersion;       // Configuration version (empty if diverse)
  TString fGConfigName;          // Configuration name (empty if diverse)
  TString fGFilterType;          // filter type (p, pgt, nf) (empty if diverse)
  TString fGReadoutParam;        // readout parameter (zs, nozs, testpattern) (empty if diverse)
  TString fGTestPattern;         // value of testpattern (for readout param) (empty if diverse)
  TString fGTrackletMode;        // tracklet mode (trk, csmtrk, notrk) (empty if diverse)
  TString fGTrackletDef;         // definition for tracklet mode trk (empty if diverse)
  TString fGTriggerSetup;        // trigger setup (ptrg, autotrg, autol0) (empty if diverse)
  TString fGAddOptions;          // additional options (nopm, nion) (empty if diverse)
  TString fRunType;              // the type of run (physics, pedestal, ...)
  UInt_t  fStartTime;            // value from GetStartTimeDCSQuery
  UInt_t  fEndTime;              // value from GetiEndTimeDCSQuery
  
  // individual configuration parameters
  TObjArray *fFEEArr;            // config param of the individual chambers
  TObjArray *fPTRArr;            // config param of the pretrigger

  AliTRDCalDCSGTU *fGTUObj;      // GTU object

  ClassDef(AliTRDCalDCSv2,1)     // TRD calibration class v2 for TRD DCS parameters

};
#endif

