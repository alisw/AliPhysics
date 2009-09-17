#ifndef AliTRDCALDCSPTRFeb_H
#define AliTRDCALDCSPTRFeb_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalDCSPTRFeb.h 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD GTU configuration parameters               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class TString;

class AliTRDCalDCSPTRFeb : public TNamed {

 public:

  AliTRDCalDCSPTRFeb();
  AliTRDCalDCSPTRFeb(const char *name, const char *title);
  AliTRDCalDCSPTRFeb(const AliTRDCalDCSPTRFeb &);
  virtual ~AliTRDCalDCSPTRFeb() { };

  TString GetControlBoxSide()                         { return fSide;                         }
  TString GetDetectorName()                           { return fName;                         }
  Int_t   GetControlBoxPrimary()                      { return fPrimary;                      }
  UInt_t  GetControlBoxSide()                         { return fSide;                         }
  UInt_t  GetControlBoxPrimary()                      { return fPrimary;                      }

  void    SetControlBoxSide(TString bs)               { fSide = bs;                           }
  void    SetDetectorName(TString bs)                 { fName = bs;                           }
  void    SetControlBoxPrimary(Int_t bp)              { fPrimary = bp;                        }
  void    SetClkLb(UInt_t cl)                         { fClkLb = cl;                          }
  void    SetClkHb(UInt_t ch)                         { fClkHb = ch;                          }
  
  UInt_t  GetCh0CountLb()                             { return fCh0CountLb;                   }
  UInt_t  GetCh0CountHb()                             { return fCh0CountHb;                   }
  UInt_t  GetCh1CountLb()                             { return fCh1CountLb;                   }
  UInt_t  GetCh1CountHb()                             { return fCh1CountHb;                   }
  UInt_t  GetCh2CountLb()                             { return fCh2CountLb;                   }
  UInt_t  GetCh2CountHb()                             { return fCh2CountHb;                   }
  UInt_t  GetCh3CountLb()                             { return fCh3CountLb;                   }
  UInt_t  GetCh3CountHb()                             { return fCh3CountHb;                   }
  UInt_t  GetCh4CountLb()                             { return fCh4CountLb;                   }
  UInt_t  GetCh4CountHb()                             { return fCh4CountHb;                   }
  UInt_t  GetCh5CountLb()                             { return fCh5CountLb;                   }
  UInt_t  GetCh5CountHb()                             { return fCh5CountHb;                   }
  UInt_t  GetCh6CountLb()                             { return fCh6CountLb;                   }
  UInt_t  GetCh6CountHb()                             { return fCh6CountHb;                   }
  UInt_t  GetCh7CountLb()                             { return fCh7CountLb;                   }
  UInt_t  GetCh7CountHb()                             { return fCh7CountHb;                   }
  UInt_t  GetCh8CountLb()                             { return fCh8CountLb;                   }
  UInt_t  GetCh8CountHb()                             { return fCh8CountHb;                   }
  UInt_t  GetCh9CountLb()                             { return fCh9CountLb;                   }
  UInt_t  GetCh9CountHb()                             { return fCh9CountHb;                   }
  UInt_t  GetCh10CountLb()                            { return fCh10CountLb;                  }
  UInt_t  GetCh10CountHb()                            { return fCh10CountHb;                  }
  UInt_t  GetCh11CountLb()                            { return fCh11CountLb;                  }
  UInt_t  GetCh11CountHb()                            { return fCh11CountHb;                  }
  UInt_t  GetTrigParallel0Lb()                        { return fTrigParallel0Lb;              }
  UInt_t  GetTrigParallel0Hb()                        { return fTrigParallel0Hb;              }
  UInt_t  GetTrigParallel1Lb()                        { return fTrigParallel1Lb;              }
  UInt_t  GetTrigParallel1Hb()                        { return fTrigParallel1Hb;              }
  UInt_t  GetPtChDelay0()                             { return fPtChDelay0;                   }
  UInt_t  GetPtChDelay1()                             { return fPtChDelay1;                   }
  UInt_t  GetPtChDelay2()                             { return fPtChDelay2;                   }
  UInt_t  GetPtChDelay3()                             { return fPtChDelay3;                   }
  UInt_t  GetPtChDelay4()                             { return fPtChDelay4;                   }
  UInt_t  GetPtChDelay5()                             { return fPtChDelay5;                   }
  UInt_t  GetPtChDelay6()                             { return fPtChDelay6;                   }
  UInt_t  GetPtChDelay7()                             { return fPtChDelay7;                   }
  UInt_t  GetPtChDelay8()                             { return fPtChDelay8;                   }
  UInt_t  GetPtChDelay9()                             { return fPtChDelay9;                   }
  UInt_t  GetPtChDelay10()                            { return fPtChDelay10;                  }
  UInt_t  GetPtChDelay11()                            { return fPtChDelay11;                  }
  UInt_t  GetPtChThr0()                               { return fPtChThr0;                     }
  UInt_t  GetPtChThr1()                               { return fPtChThr1;                     }
  UInt_t  GetPtChThr2()                               { return fPtChThr2;                     }
  UInt_t  GetPtChThr3()                               { return fPtChThr3;                     }
  UInt_t  GetPtChThr4()                               { return fPtChThr4;                     }
  UInt_t  GetPtChThr5()                               { return fPtChThr5;                     }
  UInt_t  GetPtChThr6()                               { return fPtChThr6;                     }
  UInt_t  GetPtChThr7()                               { return fPtChThr7;                     }
  UInt_t  GetPtChThr8()                               { return fPtChThr8;                     }
  UInt_t  GetPtChThr9()                               { return fPtChThr9;                     }
  UInt_t  GetPtChThr10()                              { return fPtChThr10;                    }
  UInt_t  GetPtChThr11()                              { return fPtChThr11;                    }


  void    SetCh0CountLb(UInt_ ar)                     { fCh0CountLb = ar;                     }
  void    SetCh0CountHb(UInt_ ar)                     { fCh0CountHb = ar;                     }
  void    SetCh1CountLb(UInt_ ar)                     { fCh1CountLb = ar;                     }
  void    SetCh1CountHb(UInt_ ar)                     { fCh1CountHb = ar;                     }
  void    SetCh2CountLb(UInt_ ar)                     { fCh2CountLb = ar;                     }
  void    SetCh2CountHb(UInt_ ar)                     { fCh2CountHb = ar;                     }
  void    SetCh3CountLb(UInt_ ar)                     { fCh3CountLb = ar;                     }
  void    SetCh3CountHb(UInt_ ar)                     { fCh3CountHb = ar;                     }
  void    SetCh4CountLb(UInt_ ar)                     { fCh4CountLb = ar;                     }
  void    SetCh4CountHb(UInt_ ar)                     { fCh4CountHb = ar;                     }
  void    SetCh5CountLb(UInt_ ar)                     { fCh5CountLb = ar;                     }
  void    SetCh5CountHb(UInt_ ar)                     { fCh5CountHb = ar;                     }
  void    SetCh6CountLb(UInt_ ar)                     { fCh6CountLb = ar;                     }
  void    SetCh6CountHb(UInt_ ar)                     { fCh6CountHb = ar;                     }
  void    SetCh7CountLb(UInt_ ar)                     { fCh7CountLb = ar;                     }
  void    SetCh7CountHb(UInt_ ar)                     { fCh7CountHb = ar;                     }
  void    SetCh8CountLb(UInt_ ar)                     { fCh8CountLb = ar;                     }
  void    SetCh8CountHb(UInt_ ar)                     { fCh8CountHb = ar;                     }
  void    SetCh9CountLb(UInt_ ar)                     { fCh9CountLb = ar;                     }
  void    SetCh9CountHb(UInt_ ar)                     { fCh9CountHb = ar;                     }
  void    SetCh10CountLb(UInt_ ar)                    { fCh10CountLb = ar;                    }
  void    SetCh10CountHb(UInt_ ar)                    { fCh10CountHb = ar;                    }
  void    SetCh11CountLb(UInt_ ar)                    { fCh11CountLb = ar;                    }
  void    SetCh11CountHb(UInt_ ar)                    { fCh11CountHb = ar;                    }
  void    SetTrigParallel0Lb(UInt_ ar)                { fTrigParallel0Lb = ar;                }
  void    SetTrigParallel0Hb(UInt_ ar)                { fTrigParallel0Hb = ar;                }
  void    SetTrigParallel1Lb(UInt_ ar)                { fTrigParallel1Lb = ar;                }
  void    SetTrigParallel1Hb(UInt_ ar)                { fTrigParallel1Hb = ar;                }
  void    SetPtChDelay0(UInt_ ar)                     { fPtChDelay0 = ar;                     }
  void    SetPtChDelay1(UInt_ ar)                     { fPtChDelay1 = ar;                     }
  void    SetPtChDelay2(UInt_ ar)                     { fPtChDelay2 = ar;                     }
  void    SetPtChDelay3(UInt_ ar)                     { fPtChDelay3 = ar;                     }
  void    SetPtChDelay4(UInt_ ar)                     { fPtChDelay4 = ar;                     }
  void    SetPtChDelay5(UInt_ ar)                     { fPtChDelay5 = ar;                     }
  void    SetPtChDelay6(UInt_ ar)                     { fPtChDelay6 = ar;                     }
  void    SetPtChDelay7(UInt_ ar)                     { fPtChDelay7 = ar;                     }
  void    SetPtChDelay8(UInt_ ar)                     { fPtChDelay8 = ar;                     }
  void    SetPtChDelay9(UInt_ ar)                     { fPtChDelay9 = ar;                     }
  void    SetPtChDelay10(UInt_ ar)                    { fPtChDelay10 = ar;                    }
  void    SetPtChDelay11(UInt_ ar)                    { fPtChDelay11 = ar;                    }
  void    SetPtChThr0(UInt_ ar)                       { fPtChThr0 = ar;                       }
  void    SetPtChThr1(UInt_ ar)                       { fPtChThr1 = ar;                       }
  void    SetPtChThr2(UInt_ ar)                       { fPtChThr2 = ar;                       }
  void    SetPtChThr3(UInt_ ar)                       { fPtChThr3 = ar;                       }
  void    SetPtChThr4(UInt_ ar)                       { fPtChThr4 = ar;                       }
  void    SetPtChThr5(UInt_ ar)                       { fPtChThr5 = ar;                       }
  void    SetPtChThr6(UInt_ ar)                       { fPtChThr6 = ar;                       }
  void    SetPtChThr7(UInt_ ar)                       { fPtChThr7 = ar;                       }
  void    SetPtChThr8(UInt_ ar)                       { fPtChThr8 = ar;                       }
  void    SetPtChThr9(UInt_ ar)                       { fPtChThr9 = ar;                       }
  void    SetPtChThr10(UInt_ ar)                      { fPtChThr10 = ar;                      }
  void    SetPtChThr11(UInt_ ar)                      { fPtChThr11 = ar;                      }


 protected:
  TString fSide; // side of the control box, either A, B or C 
  TString fName; // 1 if its the primary control box, 2 for backup
  Int_t   fPrimary; // 1 if its the primary control box, 2 for backup

  UInt_t  fClkLb;
  UInt_t  fClkHb;
  
  UInt_t  fCh0CountLb;
  UInt_t  fCh0CountHb;
  UInt_t  fCh1CountLb;
  UInt_t  fCh1CountHb;
  UInt_t  fCh2CountLb;
  UInt_t  fCh2CountHb;
  UInt_t  fCh3CountLb;
  UInt_t  fCh3CountHb;
  UInt_t  fCh4CountLb;
  UInt_t  fCh4CountHb;
  UInt_t  fCh5CountLb;
  UInt_t  fCh5CountHb;
  UInt_t  fCh6CountLb;
  UInt_t  fCh6CountHb;
  UInt_t  fCh7CountLb;
  UInt_t  fCh7CountHb;
  UInt_t  fCh8CountLb;
  UInt_t  fCh8CountHb;
  UInt_t  fCh9CountLb;
  UInt_t  fCh9CountHb;
  UInt_t  fCh10CountLb;
  UInt_t  fCh10CountHb;
  UInt_t  fCh11CountLb;
  UInt_t  fCh11CountHb;
  UInt_t  fTrigParallel0Lb;
  UInt_t  fTrigParallel0Hb;
  UInt_t  fTrigParallel1Lb;
  UInt_t  fTrigParallel1Hb;
  UInt_t  fPtChDelay0;
  UInt_t  fPtChDelay1;
  UInt_t  fPtChDelay2;
  UInt_t  fPtChDelay3;
  UInt_t  fPtChDelay4;
  UInt_t  fPtChDelay5;
  UInt_t  fPtChDelay6;
  UInt_t  fPtChDelay7;
  UInt_t  fPtChDelay8;
  UInt_t  fPtChDelay9;
  UInt_t  fPtChDelay10;
  UInt_t  fPtChDelay11;
  UInt_t  fPtChThr0;
  UInt_t  fPtChThr1;
  UInt_t  fPtChThr2;
  UInt_t  fPtChThr3;
  UInt_t  fPtChThr4;
  UInt_t  fPtChThr5;
  UInt_t  fPtChThr6;
  UInt_t  fPtChThr7;
  UInt_t  fPtChThr8;
  UInt_t  fPtChThr9;
  UInt_t  fPtChThr10;
  UInt_t  fPtChThr11;

  ClassDef(AliTRDCalDCSPTRFeb,1)      //  TRD calibration class for TRD GTU parameters

};
#endif
