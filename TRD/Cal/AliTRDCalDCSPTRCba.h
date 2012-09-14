#ifndef ALITRDCALDCSPTRCBA_H
#define ALITRDCALDCSPTRCBA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalDCSPTRCba.h 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD GTU configuration parameters               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class TString;

class AliTRDCalDCSPTRCba : public TNamed {

 public:

  AliTRDCalDCSPTRCba();
  AliTRDCalDCSPTRCba(const char *name, const char *title);
  AliTRDCalDCSPTRCba(const AliTRDCalDCSPTRCba &);
  virtual ~AliTRDCalDCSPTRCba() { };

  
  UInt_t  GetChDelayT0() const                        { return fChDelayT0;                    }
  UInt_t  GetChDelayV0() const                        { return fChDelayV0;                    }
  UInt_t  GetChDelayV1() const                        { return fChDelayV1;                    }
  UInt_t  GetChDelayV2() const                        { return fChDelayV2;                    }
  UInt_t  GetChDelayV3() const                        { return fChDelayV3;                    }
  UInt_t  GetChDisableT0() const                      { return fChDisableT0;                  }
  UInt_t  GetChDisableV0() const                      { return fChDisableV0;                  }
  UInt_t  GetChDisableV1() const                      { return fChDisableV1;                  }
  UInt_t  GetChDisableV2() const                      { return fChDisableV2;                  }
  UInt_t  GetChDisableV3() const                      { return fChDisableV3;                  }
  
  UInt_t  GetTo27ParralelLb() const                   { return fTo27ParralelLb;               }
  UInt_t  GetTo27ParralelHb() const                   { return fTo27ParralelHb;               }
  UInt_t  GetTo28ParralelLb() const                   { return fTo28ParralelLb;               }
  UInt_t  GetTo28ParralelHb() const                   { return fTo28ParralelHb;               }
  UInt_t  GetTo29ParralelLb() const                   { return fTo29ParralelLb;               }
  UInt_t  GetTo29ParralelHb() const                   { return fTo29ParralelHb;               }
  UInt_t  GetTo30ParralelLb() const                   { return fTo30ParralelLb;               }
  UInt_t  GetTo30ParralelHb() const                   { return fTo30ParralelHb;               }
  UInt_t  GetTo31ParralelLb() const                   { return fTo31ParralelLb;               }
  UInt_t  GetTo31ParralelHb() const                   { return fTo31ParralelHb;               }
  UInt_t  GetTo32ParralelLb() const                   { return fTo32ParralelLb;               }
  UInt_t  GetTo32ParralelHb() const                   { return fTo32ParralelHb;               }
  UInt_t  GetTo33ParralelLb() const                   { return fTo33ParralelLb;               }
  UInt_t  GetTo33ParralelHb() const                   { return fTo33ParralelHb;               }
  UInt_t  GetTo34ParralelLb() const                   { return fTo34ParralelLb;               }
  UInt_t  GetTo34ParralelHb() const                   { return fTo34ParralelHb;               }
  UInt_t  GetTo35ParralelLb() const                   { return fTo35ParralelLb;               }
  UInt_t  GetTo35ParralelHb() const                   { return fTo35ParralelHb;               }
  UInt_t  GetTo36ParralelLb() const                   { return fTo36ParralelLb;               }
  UInt_t  GetTo36ParralelHb() const                   { return fTo36ParralelHb;               }
  
  UInt_t  GetClkLb() const                            { return fClkLb;                        }
  UInt_t  GetClkHb() const                            { return fClkHb;                        }
  
  UInt_t  GetBitsToCbB42Lb() const                    { return fBitsToCbB42Lb;                }
  UInt_t  GetBitsToCbB42Hb() const                    { return fBitsToCbB42Hb;                }
  UInt_t  GetBitsToCbB43Lb() const                    { return fBitsToCbB43Lb;                }
  UInt_t  GetBitsToCbB43Hb() const                    { return fBitsToCbB43Hb;                }
  UInt_t  GetBitsToCbB44Lb() const                    { return fBitsToCbB44Lb;                }
  UInt_t  GetBitsToCbB44Hb() const                    { return fBitsToCbB44Hb;                }
  UInt_t  GetBitsToCbB45Lb() const                    { return fBitsToCbB45Lb;                }
  UInt_t  GetBitsToCbB45Hb() const                    { return fBitsToCbB45Hb;                }
  
  TString GetControlBoxSide() const                   { return fSide;                         }
  Int_t   GetControlBoxPrimary() const                { return fPrimary;                      }

  void    SetControlBoxSide(TString bs)               { fSide = bs;                           }
  void    SetControlBoxPrimary(Int_t bp)              { fPrimary = bp;                        }
  
  void    SetClkLb(UInt_t cl)                         { fClkLb = cl;                          }
  void    SetClkHb(UInt_t ch)                         { fClkHb = ch;                          }
 
  void    SetTo27ParralelHb(UInt_t tp)                { fTo27ParralelHb = tp;                 }                                            
  void    SetTo27ParralelLb(UInt_t tp)                { fTo27ParralelLb = tp;                 }                                            
  void    SetTo28ParralelHb(UInt_t tp)                { fTo28ParralelHb = tp;                 }                                            
  void    SetTo28ParralelLb(UInt_t tp)                { fTo28ParralelLb = tp;                 }                                            
  void    SetTo29ParralelHb(UInt_t tp)                { fTo29ParralelHb = tp;                 }                                            
  void    SetTo29ParralelLb(UInt_t tp)                { fTo29ParralelLb = tp;                 }                                            
  void    SetTo30ParralelHb(UInt_t tp)                { fTo30ParralelHb = tp;                 }                                            
  void    SetTo30ParralelLb(UInt_t tp)                { fTo30ParralelLb = tp;                 }                                            
  void    SetTo31ParralelHb(UInt_t tp)                { fTo31ParralelHb = tp;                 }                                            
  void    SetTo31ParralelLb(UInt_t tp)                { fTo31ParralelLb = tp;                 }
  void    SetTo32ParralelHb(UInt_t tp)                { fTo32ParralelHb = tp;                 }
  void    SetTo32ParralelLb(UInt_t tp)                { fTo32ParralelLb = tp;                 }
  void    SetTo33ParralelHb(UInt_t tp)                { fTo33ParralelHb = tp;                 }
  void    SetTo33ParralelLb(UInt_t tp)                { fTo33ParralelLb = tp;                 }
  void    SetTo34ParralelHb(UInt_t tp)                { fTo34ParralelHb = tp;                 }
  void    SetTo34ParralelLb(UInt_t tp)                { fTo34ParralelLb = tp;                 }
  void    SetTo35ParralelHb(UInt_t tp)                { fTo35ParralelHb = tp;                 }
  void    SetTo35ParralelLb(UInt_t tp)                { fTo35ParralelLb = tp;                 }
  void    SetTo36ParralelHb(UInt_t tp)                { fTo36ParralelHb = tp;                 }
  void    SetTo36ParralelLb(UInt_t tp)                { fTo36ParralelLb = tp;                 }

  void    SetBitsToCbB42Hb(UInt_t bc)                 { fBitsToCbB42Hb = bc;                  }
  void    SetBitsToCbB42Lb(UInt_t bc)                 { fBitsToCbB42Lb = bc;                  }
  void    SetBitsToCbB43Hb(UInt_t bc)                 { fBitsToCbB43Hb = bc;                  }
  void    SetBitsToCbB43Lb(UInt_t bc)                 { fBitsToCbB43Lb = bc;                  }
  void    SetBitsToCbB44Hb(UInt_t bc)                 { fBitsToCbB44Hb = bc;                  }
  void    SetBitsToCbB44Lb(UInt_t bc)                 { fBitsToCbB44Lb = bc;                  }
  void    SetBitsToCbB45Hb(UInt_t bc)                 { fBitsToCbB45Hb = bc;                  }
  void    SetBitsToCbB45Lb(UInt_t bc)                 { fBitsToCbB45Lb = bc;                  }

  void    SetChDelayT0(UInt_t cd)                     { fChDelayT0 = cd;                      }
  void    SetChDelayV0(UInt_t cd)                     { fChDelayV0 = cd;                      }
  void    SetChDelayV1(UInt_t cd)                     { fChDelayV1 = cd;                      }
  void    SetChDelayV2(UInt_t cd)                     { fChDelayV2 = cd;                      }
  void    SetChDelayV3(UInt_t cd)                     { fChDelayV3 = cd;                      }
  void    SetChDisableT0(UInt_t cd)                   { fChDisableT0 = cd;                    }
  void    SetChDisableV0(UInt_t cd)                   { fChDisableV0 = cd;                    }
  void    SetChDisableV1(UInt_t cd)                   { fChDisableV1 = cd;                    }
  void    SetChDisableV2(UInt_t cd)                   { fChDisableV2 = cd;                    }
  void    SetChDisableV3(UInt_t cd)                   { fChDisableV3 = cd;                    }

 protected:
  TString fSide; // side of the control box, either A, B or C 
  Int_t   fPrimary; // 1 if its the primary control box, 2 for backup

  UInt_t  fChDelayT0; // value from the ChDelayT0 tag in the pt box type CbA
  UInt_t  fChDelayV0; // value from the ChDelayV0 tag in the pt box type CbA
  UInt_t  fChDelayV1; // value from the ChDelayV1 tag in the pt box type CbA
  UInt_t  fChDelayV2; // value from the ChDelayV2 tag in the pt box type CbA
  UInt_t  fChDelayV3; // value from the ChDelayV3 tag in the pt box type CbA
  UInt_t  fChDisableT0; // value from the ChDisableT0 tag in the pt box type CbA
  UInt_t  fChDisableV0; // value from the ChDisableV0 tag in the pt box type CbA
  UInt_t  fChDisableV1; // value from the ChDisableV1 tag in the pt box type CbA
  UInt_t  fChDisableV2; // value from the ChDisableV2 tag in the pt box type CbA
  UInt_t  fChDisableV3; // value from the ChDisableV3 tag in the pt box type CbA
  
  UInt_t  fTo27ParralelLb; // value from the To27Parralel low bit tag in the pt box type CbA
  UInt_t  fTo27ParralelHb; // value from the To27Parralel high bit tag in the pt box type CbA
  UInt_t  fTo28ParralelLb; // value from the To28Parralel low bit tag in the pt box type CbA
  UInt_t  fTo28ParralelHb; // value from the To28Parralel high bit tag in the pt box type CbA
  UInt_t  fTo29ParralelLb; // value from the To29Parralel low bit tag in the pt box type CbA
  UInt_t  fTo29ParralelHb; // value from the To29Parralel high bit tag in the pt box type CbA
  UInt_t  fTo30ParralelLb; // value from the To30Parralel low bit tag in the pt box type CbA
  UInt_t  fTo30ParralelHb; // value from the To30Parralel high bit tag in the pt box type CbA
  UInt_t  fTo31ParralelLb; // value from the To31Parralel low bit tag in the pt box type CbA
  UInt_t  fTo31ParralelHb; // value from the To31Parralel high bit tag in the pt box type CbA
  UInt_t  fTo32ParralelLb; // value from the To32Parralel low bit tag in the pt box type CbA
  UInt_t  fTo32ParralelHb; // value from the To32Parralel high bit tag in the pt box type CbA
  UInt_t  fTo33ParralelLb; // value from the To33Parralel low bit tag in the pt box type CbA
  UInt_t  fTo33ParralelHb; // value from the To33Parralel high bit tag in the pt box type CbA
  UInt_t  fTo34ParralelLb; // value from the To34Parralel low bit tag in the pt box type CbA
  UInt_t  fTo34ParralelHb; // value from the To34Parralel high bit tag in the pt box type CbA
  UInt_t  fTo35ParralelLb; // value from the To35Parralel low bit tag in the pt box type CbA
  UInt_t  fTo35ParralelHb; // value from the To35Parralel high bit tag in the pt box type CbA
  UInt_t  fTo36ParralelLb; // value from the To36Parralel low bit tag in the pt box type CbA
  UInt_t  fTo36ParralelHb; // value from the To36Parralel high bit tag in the pt box type CbA
  
  UInt_t  fClkLb; // value from the Clk low bit tag in the pt box type CbA
  UInt_t  fClkHb; // value from the Clk high bit tag in the pt box type CbA
  
  UInt_t  fBitsToCbB42Lb; // value from the BitsToCbB42 low bit tag in the pt box type CbA
  UInt_t  fBitsToCbB42Hb; // value from the BitsToCbB42 high bit tag in the pt box type CbA
  UInt_t  fBitsToCbB43Lb; // value from the BitsToCbB43 low bit tag in the pt box type CbA
  UInt_t  fBitsToCbB43Hb; // value from the BitsToCbB43 high bit tag in the pt box type CbA
  UInt_t  fBitsToCbB44Lb; // value from the BitsToCbB44 low bit tag in the pt box type CbA
  UInt_t  fBitsToCbB44Hb; // value from the BitsToCbB44 high bit tag in the pt box type CbA
  UInt_t  fBitsToCbB45Lb; // value from the BitsToCbB45 low bit tag in the pt box type CbA
  UInt_t  fBitsToCbB45Hb; // value from the BitsToCbB45 high bit tag in the pt box type CbA


  ClassDef(AliTRDCalDCSPTRCba,1)      //  TRD calibration class for TRD GTU parameters

};
#endif
