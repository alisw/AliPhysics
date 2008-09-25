#ifndef AliTRDCALDCS_H
#define AliTRDCALDCS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalDCS.h 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD DCS parameters                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"
#include "TObjArray.h"

class TString;

class AliTRDCalDCSFEE;
class AliTRDCalDCSPTR;
class AliTRDCalDCSGTU;

class AliTRDCalDCS : public TNamed {

 public:

  AliTRDCalDCS();
  AliTRDCalDCS(const Text_t *name, const Text_t *title);
  AliTRDCalDCS(const AliTRDCalDCS &cd);
  AliTRDCalDCS &operator=(const AliTRDCalDCS &cd);
  virtual ~AliTRDCalDCS() { };

  void    EvaluateGlobalParameters();
  void    SetFEEArr(TObjArray *fa)             { fFEEArr              = fa;    }
  void    SetPTRArr(TObjArray *pa)             { fPTRArr              = pa;    }
  void    SetGTUArr(TObjArray *ga)             { fGTUArr              = ga;    }
  
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
  TString GetGlobalAddOptions() const          { return fGAddOptions;           } 
  TObjArray*       GetFEEArr() const           { return fFEEArr;               }
  TObjArray*       GetPTRArr() const           { return fPTRArr;               }
  TObjArray*       GetGTUArr() const           { return fGTUArr;               }
  AliTRDCalDCSFEE* GetCalDCSFEEObj(Int_t det) 
  		  	          { return (AliTRDCalDCSFEE*)fFEEArr->At(det); }
  AliTRDCalDCSPTR* GetCalDCSPTRObj(Int_t det) 
  			          { return (AliTRDCalDCSPTR*)fPTRArr->At(det); }
  AliTRDCalDCSGTU* GetCalDCSGTUObj(Int_t det) 
           		          { return (AliTRDCalDCSGTU*)fGTUArr->At(det); }

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
  
  //individual configuration parameters
  TObjArray *fFEEArr; // config param of the individual chambers
  TObjArray *fPTRArr; // config param of the pretrigger
  TObjArray *fGTUArr; // config param of the GTU

  ClassDef(AliTRDCalDCS,2)         //  TRD calibration class for TRD DCS parameters

};
#endif
