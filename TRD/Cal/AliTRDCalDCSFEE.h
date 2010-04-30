#ifndef ALITRDCALDCSFEE_H
#define ALITRDCALDCSFEE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalDCSFEE.h 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for FEE configuration parameters                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TString.h"

//class TString;

class AliTRDCalDCSFEE : public TObject {

 public:

  AliTRDCalDCSFEE();
  AliTRDCalDCSFEE(const AliTRDCalDCSFEE &c);
  virtual ~AliTRDCalDCSFEE() { };
  AliTRDCalDCSFEE &operator=(const AliTRDCalDCSFEE &c);

  void    SetStatusBit(Int_t stbit)                  { fStatusBit           = (Char_t)stbit;        }
  void    SetDCSid(Int_t dcsid)                      { fDCSID               = (Char_t)dcsid;        }
  void    SetSM(Int_t smid)                          { fSM                  = (Char_t)smid;         }
  void    SetStack(Int_t stid)                       { fStack               = (Char_t)stid;         }
  void    SetLayer(Int_t lyid)                       { fLayer               = (Char_t)lyid;         }
  void    SetNumberOfTimeBins(Int_t value)           { fNumberOfTimeBins    = (UShort_t)value;      }
  void    SetConfigTag(Int_t cfgt)                   { fConfigTag           = (UShort_t)cfgt;       }
  void    SetSingleHitThres(Int_t sht)               { fSingleHitThres      = (Short_t)sht;         }
  void    SetThreePadClustThres(Int_t tpct)          { fThrPdClsThres       = (Short_t)tpct;        }
  void    SetSelectiveNoZS(Int_t snzs)               { fSelNoZS             = (Short_t)snzs;        }
  void    SetFastStatNoise(Int_t fstn)               { fFastStatNoise       = (Short_t)fstn;        }
  void    SetTCFilterWeight(Int_t tcfw)              { fTCFilterWeight      = (Short_t)tcfw;        }
  void    SetTCFilterShortDecPar(Int_t sdp)          { fTCFilterShortDecPar = (Short_t)sdp;         }
  void    SetTCFilterLongDecPar(Int_t ldp)           { fTCFilterLongDecPar  = (Short_t)ldp;         }
  void    SetGainTableRocSerial(Int_t gts)           { fGainTableRocSerial  = (UChar_t)gts;         }
  void    SetFilterType(TString fity)                { fFilterType          = fity;                 }
  void    SetReadoutParam(TString rpar)              { fReadoutParam        = rpar;                 }
  void    SetTestPattern(TString tpat)               { fTestPattern         = tpat;                 }
  void    SetTrackletMode(TString tmde)              { fTrackletMode        = tmde;                 }
  void    SetTrackletDef(TString tdef)               { fTrackletDef         = tdef;                 }
  void    SetTriggerSetup(TString trse)              { fTriggerSetup        = trse;                 }
  void    SetAddOptions(TString adop)                { fAddOptions          = adop;                 }
  void    SetConfigName(TString cfgn)                { fConfigName          = cfgn;                 }
  void    SetConfigVersion(TString cfgv)             { fConfigVersion       = cfgv;                 }
  void    SetGainTableName(TString gt)               { fGainTableName       = gt;                   }
  void    SetGainTableDesc(TString gd)               { fGainTableDesc       = gd;                   }
  void    SetGainTableRocType(TString gr)            { fGainTableRocType    = gr;                   }
  void    SetMCMGlobalState(Int_t r,Int_t m,Int_t g) { fRStateGSM[r][m]     = g;                    }
  void    SetMCMStateNI(Int_t r,Int_t m,Int_t v)     { fRStateNI[r][m]      = v;                    }
  void    SetMCMEventCnt(Int_t r,Int_t m,Int_t v)    { fRStateEV[r][m]      = v;                    }
  void    SetMCMPtCnt(Int_t r,Int_t m,Int_t v)       { fRStatePTRG[r][m]    = v;                    }
  void    SetGainTableAdcdac(Int_t r,Int_t m,Int_t v){ fGainTableAdcdac[r][m]         = (Char_t)v;  }
  void    SetGainTableFgfn(Int_t r,Int_t m,Int_t a,Int_t v) { fGainTableFgfn[r][m][a] = (Short_t)v; }
  void    SetGainTableFgan(Int_t r,Int_t m,Int_t a,Int_t v) { fGainTableFgan[r][m][a] = (Char_t)v;  }

  Int_t   GetStatusBit() const                       { return (Int_t)fStatusBit;                    }
  Int_t   GetDCSid() const                           { return (Int_t)fDCSID;                        }
  Int_t   GetSM() const                              { return (Int_t)fSM;                           }
  Int_t   GetStack() const                           { return (Int_t)fStack;                        }
  Int_t   GetLayer() const                           { return (Int_t)fLayer;                        }
  Int_t   GetNumberOfTimeBins() const                { return (Int_t)fNumberOfTimeBins;             }
  Int_t   GetConfigTag() const                       { return (Int_t)fConfigTag;                    }
  Int_t   GetSingleHitThres() const                  { return (Int_t)fSingleHitThres;               }
  Int_t   GetThreePadClustThres() const              { return (Int_t)fThrPdClsThres;                }
  Int_t   GetSelectiveNoZS() const                   { return (Int_t)fSelNoZS;                      }
  Int_t   GetTCFilterWeight() const                  { return (Int_t)fTCFilterWeight;               }
  Int_t   GetTCFilterShortDecPar() const             { return (Int_t)fTCFilterShortDecPar;          }
  Int_t   GetTCFilterLongDecPar() const              { return (Int_t)fTCFilterLongDecPar;           }
  Int_t   GetFastStatNoise() const                   { return (Int_t)fFastStatNoise;                }
  Int_t   GetGainTableRocSerial() const              { return (Int_t)fGainTableRocSerial;           }
  TString GetFilterType() const                      { return fFilterType;                          }
  TString GetReadoutParam() const                    { return fReadoutParam;                        }
  TString GetTestPattern() const                     { return fTestPattern;                         }
  TString GetTrackletMode() const                    { return fTrackletMode;                        }
  TString GetTrackletDef() const                     { return fTrackletDef;                         }
  TString GetTriggerSetup() const                    { return fTriggerSetup;                        }
  TString GetAddOptions() const                      { return fAddOptions;                          }
  TString GetConfigName() const                      { return fConfigName;                          }
  TString GetConfigVersion() const                   { return fConfigVersion;                       }
  TString GetGainTableName() const                   { return fGainTableName;                       }
  TString GetGainTableDesc() const                   { return fGainTableDesc;                       }
  TString GetGainTableRocType() const                { return fGainTableRocType;                    }
  Int_t   GetMCMGlobalState(Int_t r,Int_t m) const   { return (UChar_t)fRStateGSM[r][m];            }
  Int_t   GetMCMStateNI(Int_t r,Int_t m) const       { return (UChar_t)fRStateNI[r][m];             }
  Int_t   GetMCMEventCnt(Int_t r,Int_t m) const      { return fRStateEV[r][m];                      }
  Int_t   GetMCMPtCnt(Int_t r,Int_t m) const         { return fRStatePTRG[r][m];                    }
  Int_t   GetGainTableAdcdac(Int_t r,Int_t m) const  { return (Int_t)fGainTableAdcdac[r][m];        }
  Int_t   GetGainTableFgfn(Int_t r,Int_t m,Int_t a) const  { return (Int_t)fGainTableFgfn[r][m][a]; }
  Int_t   GetGainTableFgan(Int_t r,Int_t m,Int_t a) const  { return (Int_t)fGainTableFgan[r][m][a]; }

 protected:

  static const Char_t fgkROB = 8;       // Number of readout boards
  static const Char_t fgkMCM = 18;      // Number of MCMs
  static const Char_t fgkADC = 21;      // Number of ADC channels
  
  Char_t   fStatusBit;                  // 0 if everything is OK, otherwise !=0 (see impl. file)
  Char_t   fSM;                         // the number of the supermode 0..17
  Char_t   fStack;                      // the number of the stack 0..4
  Char_t   fLayer;                      // the number of the layer 0..5
  Char_t   fGainTableFgan[(Int_t)fgkROB][(Int_t)fgkMCM][(Int_t)fgkADC]; // array of gain table fgan values
  Char_t   fGainTableAdcdac[(Int_t)fgkROB][(Int_t)fgkMCM]; // array of gain table adcdac values
  UChar_t  fRStateGSM[(Int_t)fgkROB][(Int_t)fgkMCM];  // array of the global states of the MCMs
  UChar_t  fRStateNI[(Int_t)fgkROB][(Int_t)fgkMCM];   // array of the network interface states of the MCMs
  UChar_t  fGainTableRocSerial;         // the roc serial of the chamber from the gain table
  UShort_t fDCSID;                      // ID of the DCS-Board
  UShort_t fNumberOfTimeBins;           // Number of timebins  
  UShort_t fConfigTag;                  // Configuration tag
  Short_t  fSingleHitThres;             // threshold of single hits (arg of readout param)
  Short_t  fThrPdClsThres;              // threshold of 3-pad clusters (arg of readout param)
  Short_t  fSelNoZS;                    // write every fSelNoZS'th event without ZS
  Short_t  fTCFilterWeight;             // tail cancellation filter weight
  Short_t  fTCFilterShortDecPar;        // tail cancellation filter short decay parameter
  Short_t  fTCFilterLongDecPar;         // tail cancellation filter long decay parameter
  Short_t  fFastStatNoise;              // collect statistics for fast noise mode
  Short_t  fGainTableFgfn[(Int_t)fgkROB][(Int_t)fgkMCM][(Int_t)fgkADC]; // array of gain table fgfn values
  Int_t    fRStateEV[(Int_t)fgkROB][(Int_t)fgkMCM];   // array of the event counters of the MCMs
  Int_t    fRStatePTRG[(Int_t)fgkROB][(Int_t)fgkMCM]; // array of the pretrigger counters of the MCMs
  TString  fGainTableRocType;           // the roc type from the gain table
  TString  fFilterType;                 // filter type (p, pgt, nf)
  TString  fReadoutParam;               // readout parameter (zs, nozs, testpattern)
  TString  fTestPattern;                // value of testpattern (for readout param)
  TString  fTrackletMode;               // tracklet mode (trk, csmtrk, notrk)
  TString  fTrackletDef;                // definition for tracklet mode trk
  TString  fTriggerSetup;               // trigger setup (ptrg, autotrg, autol0)
  TString  fAddOptions;                 // additional options (nopm, nion)
  TString  fConfigName;                 // Configuration name
  TString  fConfigVersion;              // Configuration version
  TString  fGainTableName;              // the name of the gain table
  TString  fGainTableDesc;              // the description of the gain table

  ClassDef(AliTRDCalDCSFEE,5)          // TRD calibration class for TRD FEE parameters
};
#endif
