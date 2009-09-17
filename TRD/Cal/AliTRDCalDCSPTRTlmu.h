#ifndef AliTRDCALDCSPTRTlmu_H
#define AliTRDCALDCSPTRTlmu_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalDCSPTRTlmu.h 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD GTU configuration parameters               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class TString;

class AliTRDCalDCSPTRTlmu : public TNamed {

 public:

  AliTRDCalDCSPTRTlmu();
  AliTRDCalDCSPTRTlmu(const char *name, const char *title);
  AliTRDCalDCSPTRTlmu(const AliTRDCalDCSPTRTlmu &);
  virtual ~AliTRDCalDCSPTRTlmu() { };

  Int_t   GetRunNumber()                              { return fRunNumber;                    }
  Int_t   GetSORFlag()                                { return fSORFlag;                      }
  Int_t   GetSerial()                                 { return fSerial;                       }
  Int_t   GetDNR()                                    { return fDNR;                          }

  void    SetRunNumber(Int_t rn)                      { fRunNumber = rn;                      }
  void    SetSORFlag(Int_t fg)                        { fSORFlag = fg;                        }
  void    SetSerial(Int_t se)                         { fSerial = se;                         }
  void    SetDNR(Int_t dn)                            { fDNR = dn;                            }

  TObjArray* GetSegmentArray() const                  { return fSegmentsArr;                  }
  void SetSegmentArray(TObjArray *sa)                 { fSegmentsArr = sa;                    }

  AliTRDCalDCSGTUTgu* GetTgu() const                  { return fTgu;                          }
  void SetTgu(AliTRDCalDCSGTUTgu* tg)                 { fTgu = tg;                            }

  UInt_t  GetTo27ParralelLb()                         { return fTo27ParralelLb;               }
  UInt_t  GetTo27ParralelHb()                         { return fTo27ParralelHb;               }            
  UInt_t  GetTo28ParralelLb()                         { return fTo28ParralelLb;               }
  UInt_t  GetTo28ParralelHb()                         { return fTo28ParralelHb;               }            
  UInt_t  GetTo29ParralelLb()                         { return fTo29ParralelLb;               }
  UInt_t  GetTo29ParralelHb()                         { return fTo29ParralelHb;               }            
  UInt_t  GetTo30ParralelLb()                         { return fTo30ParralelLb;               }
  UInt_t  GetTo30ParralelHb()                         { return fTo30ParralelHb;               }            
  UInt_t  GetTo31ParralelLb()                         { return fTo31ParralelLb;               }
  UInt_t  GetTo31ParralelHb()                         { return fTo31ParralelHb;               }
  UInt_t  GetTo32ParralelLb()                         { return fTo32ParralelLb;               }
  UInt_t  GetTo32ParralelHb()                         { return fTo32ParralelHb;               }
  UInt_t  GetTo33ParralelLb()                         { return fTo33ParralelLb;               }
  UInt_t  GetTo33ParralelHb()                         { return fTo33ParralelHb;               }
  UInt_t  GetTo34ParralelLb()                         { return fTo34ParralelLb;               }
  UInt_t  GetTo34ParralelHb()                         { return fTo34ParralelHb;               }
  UInt_t  GetTo35ParralelLb()                         { return fTo35ParralelLb;               }
  UInt_t  GetTo35ParralelHb()                         { return fTo35ParralelHb;               }
  UInt_t  GetTo36ParralelLb()                         { return fTo36ParralelLb;               }
  UInt_t  GetTo36ParralelHb()                         { return fTo36ParralelHb;               }
 
  void    SetTo27ParralelHb(Int_t tp)                 { fTo27ParralelHb = tp;                 }                                            
  void    SetTo27ParralelLb(Int_t tp)                 { fTo27ParralelLb = tp;                 }                                            
  void    SetTo28ParralelHb(Int_t tp)                 { fTo28ParralelHb = tp;                 }                                            
  void    SetTo28ParralelLb(Int_t tp)                 { fTo28ParralelLb = tp;                 }                                            
  void    SetTo29ParralelHb(Int_t tp)                 { fTo29ParralelHb = tp;                 }                                            
  void    SetTo29ParralelLb(Int_t tp)                 { fTo29ParralelLb = tp;                 }                                            
  void    SetTo30ParralelHb(Int_t tp)                 { fTo30ParralelHb = tp;                 }                                            
  void    SetTo30ParralelLb(Int_t tp)                 { fTo30ParralelLb = tp;                 }                                            
  void    SetTo31ParralelHb(Int_t tp)                 { fTo31ParralelHb = tp;                 }                                            
  void    SetTo31ParralelLb(Int_t tp)                 { fTo31ParralelLb = tp;                 }
  void    SetTo32ParralelHb(Int_t tp)                 { fTo32ParralelHb = tp;                 }
  void    SetTo32ParralelLb(Int_t tp)                 { fTo32ParralelLb = tp;                 }
  void    SetTo33ParralelHb(Int_t tp)                 { fTo33ParralelHb = tp;                 }
  void    SetTo33ParralelLb(Int_t tp)                 { fTo33ParralelLb = tp;                 }
  void    SetTo34ParralelHb(Int_t tp)                 { fTo34ParralelHb = tp;                 }
  void    SetTo34ParralelLb(Int_t tp)                 { fTo34ParralelLb = tp;                 }
  void    SetTo35ParralelHb(Int_t tp)                 { fTo35ParralelHb = tp;                 }
  void    SetTo35ParralelLb(Int_t tp)                 { fTo35ParralelLb = tp;                 }
  void    SetTo36ParralelHb(Int_t tp)                 { fTo36ParralelHb = tp;                 }
  void    SetTo36ParralelLb(Int_t tp)                 { fTo36ParralelLb = tp;                 }

  UInt_t  GetBitsToCbB42Lb()                          { return fBitsToCbB42Lb;                }                                        
  UInt_t  GetBitsToCbB42Hb()                          { return fBitsToCbB42Hb;                }                                                   
  UInt_t  GetBitsToCbB43Lb()                          { return fBitsToCbB43Lb;                }                                        
  UInt_t  GetBitsToCbB43Hb()                          { return fBitsToCbB43Hb;                }                                                   
  UInt_t  GetBitsToCbB44Lb()                          { return fBitsToCbB44Lb;                }                                        
  UInt_t  GetBitsToCbB44Hb()                          { return fBitsToCbB44Hb;                }
  UInt_t  GetBitsToCbB45Lb()                          { return fBitsToCbB45Lb;                }
  UInt_t  GetBitsToCbB45Hb()                          { return fBitsToCbB45Hb;                }

  void    SetBitsToCbB42Hb(Int_t bc)                  { fBitsToCbB42Hb = bc;                  }
  void    SetBitsToCbB42Lb(Int_t bc)                  { fBitsToCbB42Lb = bc;                  }
  void    SetBitsToCbB43Hb(Int_t bc)                  { fBitsToCbB43Hb = bc;                  }
  void    SetBitsToCbB43Lb(Int_t bc)                  { fBitsToCbB43Lb = bc;                  }
  void    SetBitsToCbB44Hb(Int_t bc)                  { fBitsToCbB44Hb = bc;                  }
  void    SetBitsToCbB44Lb(Int_t bc)                  { fBitsToCbB44Lb = bc;                  }
  void    SetBitsToCbB45Hb(Int_t bc)                  { fBitsToCbB45Hb = bc;                  }
  void    SetBitsToCbB45Lb(Int_t bc)                  { fBitsToCbB45Lb = bc;                  }

  UInt_t  GetChDelayT0()                              { return fChDelayT0;                    }  
  UInt_t  GetChDelayV0()                              { return fChDelayV0;                    }
  UInt_t  GetChDelayV1()                              { return fChDelayV1;                    }
  UInt_t  GetChDelayV2()                              { return fChDelayV2;                    }
  UInt_t  GetChDelayV3()                              { return fChDelayV3;                    }
  UInt_t  GetChDisableT0()                            { return fChDisableT0;                  }
  UInt_t  GetChDisableV0()                            { return fChDisableV0;                  }
  UInt_t  GetChDisableV1()                            { return fChDisableV1;                  }
  UInt_t  GetChDisableV2()                            { return fChDisableV2;                  }
  UInt_t  GetChDisableV3()                            { return fChDisableV3;                  }

  void    SetChDelayT0(Int_t cd)                      { fChDelayT0 = cd;                      }
  void    SetChDelayV0(Int_t cd)                      { fChDelayV0 = cd;                      }
  void    SetChDelayV1(Int_t cd)                      { fChDelayV1 = cd;                      }
  void    SetChDelayV2(Int_t cd)                      { fChDelayV2 = cd;                      }
  void    SetChDelayV3(Int_t cd)                      { fChDelayV3 = cd;                      }
  void    SetChDisableT0(Int_t cd)                    { fChDisableT0 = cd;                    }
  void    SetChDisableV0(Int_t cd)                    { fChDisableV0 = cd;                    }
  void    SetChDisableV1(Int_t cd)                    { fChDisableV1 = cd;                    }
  void    SetChDisableV2(Int_t cd)                    { fChDisableV2 = cd;                    }
  void    SetChDisableV3(Int_t cd)                    { fChDisableV3 = cd;                    }

 protected:
  TString fSide; // side of the control box, either A, B or C 
  TInt_t  fPrimary; // 1 if its the primary control box, 2 for backup

  UInt_t  fClkLb;
  UInt_t  fClkHb;

  UInt_t  fTo27ParralelLb;
  UInt_t  fTo27ParralelLb;
  UInt_t  fTo28ParralelLb;                                                                                                                
  UInt_t  fTo28ParralelHb;                                                                                                                
  UInt_t  fTo29ParralelLb;                                                                                                                
  UInt_t  fTo29ParralelHb;                                                                                                                
  UInt_t  fTo30ParralelLb;                                                                                                                
  UInt_t  fTo30ParralelHb;                                                                                                                
  UInt_t  fTo31ParralelLb;                                                                                                                
  UInt_t  fTo31ParralelHb;
  UInt_t  fTo32ParralelLb;
  UInt_t  fTo32ParralelHb;
  UInt_t  fTo33ParralelLb;
  UInt_t  fTo33ParralelHb;
  UInt_t  fTo34ParralelLb;
  UInt_t  fTo34ParralelHb;
  UInt_t  fTo35ParralelLb;
  UInt_t  fTo35ParralelHb;
  UInt_t  fTo36ParralelLb;
  UInt_t  fTo36ParralelHb;

  UInt_t  fBitsToCbB42Lb;
  UInt_t  fBitsToCbB42Hb;
  UInt_t  fBitsToCbB43Lb;
  UInt_t  fBitsToCbB43Hb;
  UInt_t  fBitsToCbB44Lb;
  UInt_t  fBitsToCbB44Hb;
  UInt_t  fBitsToCbB45Lb;
  UInt_t  fBitsToCbB45Hb;

  UInt_t  fChDelayT0;
  UInt_t  fChDelayV0;
  UInt_t  fChDelayV1;
  UInt_t  fChDelayV2;
  UInt_t  fChDelayV3;
  UInt_t  fChDisableT0;
  UInt_t  fChDisableV0;
  UInt_t  fChDisableV1;
  UInt_t  fChDisableV2;
  UInt_t  fChDisableV3;

  ClassDef(AliTRDCalDCSPTRTlmu,1)      //  TRD calibration class for TRD GTU parameters

};
#endif
