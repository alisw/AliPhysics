#ifndef AliTRDCALDCSFEE_H
#define AliTRDCALDCSFEE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalDCSFEE.h 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD FEE configuration parameters               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class TString;

class AliTRDCalDCSFEE : public TNamed {

 public:

  AliTRDCalDCSFEE();
  AliTRDCalDCSFEE(const char *name, const char *title);
  virtual ~AliTRDCalDCSFEE() { };

  void    SetStatusBit(Int_t stbit)          { fStatusBit           = stbit; }  
  void    SetDCSid(Int_t dcsid)              { fDCSID               = dcsid; }  
  void    SetSM(Int_t smid)                  { fSM                  = smid;  }
  void    SetStack(Int_t stid)               { fStack               = stid;  }
  void    SetLayer(Int_t lyid)               { fLayer               = lyid;  }
  void    SetNumberOfTimeBins(Int_t value)   { fNumberOfTimeBins    = value; }
  void    SetPedestal(Int_t ped)             { fPedestal            = ped;   }
  void    SetConfigTag(Int_t cfgt)           { fConfigTag           = cfgt;  }
  void    SetSingleHitThres(Int_t sht)       { fSingleHitThres      = sht;   }
  void    SetThreePadClustThres(Int_t tpct)  { fThrPdClsThres       = tpct;  }
  void    SetSelectiveNoZS(Int_t snzs)       { fSelNoZS             = snzs;  }
  void    SetFastStatNoise(Int_t fstn)       { fFastStatNoise       = fstn;  }
  void    SetTCFilterWeight(Int_t tcfw)      { fTCFilterWeight      = tcfw;  }
  void    SetTCFilterShortDecPar(Int_t sdp)  { fTCFilterShortDecPar = sdp;   }
  void    SetTCFilterLongDecPar(Int_t ldp)   { fTCFilterLongDecPar  = ldp;   }
  void    SetFilterType(TString fity)        { fFilterType          = fity;  }
  void    SetReadoutParam(TString rpar)      { fReadoutParam        = rpar;  }
  void    SetTestPattern(TString tpat)       { fTestPattern         = tpat;  }
  void    SetTrackletMode(TString tmde)      { fTrackletMode        = tmde;  }
  void    SetTrackletDef(TString tdef)       { fTrackletDef         = tdef;  }
  void    SetTriggerSetup(TString trse)      { fTriggerSetup        = trse;  }
  void    SetAddOptions(TString adop)        { fAddOptions          = adop;  }
  void    SetConfigName(TString cfgn)        { fConfigName          = cfgn;  }
  void    SetConfigVersion(TString cfgv)     { fConfigVersion       = cfgv;  }
  void    SetGainTableID(TString id)         { fGainTableID         = id;    }

  Int_t   GetStatusBit() const               { return fStatusBit;            }
  Int_t   GetDCSid() const                   { return fDCSID;                }
  Int_t   GetSM() const                      { return fSM;                   }
  Int_t   GetStack() const                   { return fStack;                }
  Int_t   GetLayer() const                   { return fLayer;                }
  Int_t   GetNumberOfTimeBins() const        { return fNumberOfTimeBins;     }
  Int_t   GetPedestal() const                { return fPedestal;             }
  Int_t   GetConfigTag() const               { return fConfigTag;            }
  Int_t   GetSingleHitThres() const          { return fSingleHitThres;       }
  Int_t   GetThreePadClustThres() const      { return fThrPdClsThres;        }
  Int_t   GetSelectiveNoZS() const           { return fSelNoZS;              }
  Int_t   GetFastStatNoise() const           { return fFastStatNoise;        }
  Int_t   GetTCFilterWeight() const          { return fTCFilterWeight;       }
  Int_t   GetTCFilterShortDecPar() const     { return fTCFilterShortDecPar;  }
  Int_t   GetTCFilterLongDecPar() const      { return fTCFilterLongDecPar;   }
  TString GetFilterType() const              { return fFilterType;           }
  TString GetReadoutParam() const            { return fReadoutParam;         }
  TString GetTestPattern() const             { return fTestPattern;          }
  TString GetTrackletMode() const            { return fTrackletMode;         }
  TString GetTrackletDef() const             { return fTrackletDef;          }
  TString GetTriggerSetup() const            { return fTriggerSetup;         }
  TString GetAddOptions() const              { return fAddOptions;           }
  TString GetConfigName() const              { return fConfigName;           }
  TString GetConfigVersion() const           { return fConfigVersion;        }
  TString GetGainTableID() const             { return fGainTableID;          }

 protected:
  
  Int_t   fStatusBit;              // 0 if everything is OK, otherwise !=0 (see impl. file)
  Int_t   fDCSID;                  // ID of the DCS-Board
  Int_t   fSM;                     // the number of the supermode 0..17
  Int_t   fStack;                  // the number of the stack 0..4
  Int_t   fLayer;                  // the number of the layer 0..5
  Int_t   fNumberOfTimeBins;       // Number of timebins  
  Int_t   fPedestal;               // Pedestal
  Int_t   fConfigTag;              // Configuration tag
  Int_t   fSingleHitThres;         // threshold of single hits (arg of readout param)
  Int_t   fThrPdClsThres;          // threshold of 3-pad clusters (arg of readout param)
  Int_t   fSelNoZS;                // write every NNNth event without ZS
  Int_t   fFastStatNoise;          // collect statistics for fast noise mode
  Int_t   fTCFilterWeight;         // tail cancellation filter weight
  Int_t   fTCFilterShortDecPar;    // tail cancellation filter short decay parameter
  Int_t   fTCFilterLongDecPar;     // tail cancellation filter short decay parameter

  TString fFilterType;             // filter type (p, pgt, nf)
  TString fReadoutParam;           // zs, nozs, testpattern
  TString fTestPattern;            // value of testpattern (for readout param)
  TString fTrackletMode;           // trk, csmtrk, notrk
  TString fTrackletDef;            // definition for tracklet mode trk
  TString fTriggerSetup;           // ptrg, autotrg, autol0
  TString fAddOptions;             // additional options, like nopm, nion
  TString fConfigName;             // Configuration name
  TString fConfigVersion;          // Configuration version
  TString fGainTableID;            // Gain table ID

  ClassDef(AliTRDCalDCSFEE,2)      // TRD calibration class for TRD FEE parameters
};
#endif
