#ifndef AliTRDCALDCSPTRCbb_H
#define AliTRDCALDCSPTRCbb_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalDCSPTRCbb.h 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD GTU configuration parameters               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class TString;

class AliTRDCalDCSPTRCbb : public TNamed {

 public:

  AliTRDCalDCSPTRCbb();
  AliTRDCalDCSPTRCbb(const char *name, const char *title);
  AliTRDCalDCSPTRCbb(const AliTRDCalDCSPTRCbb &);
  virtual ~AliTRDCalDCSPTRCbb() { };

  TString GetControlBoxSide()                         { return fSide;                         }
  Int_t   GetControlBoxPrimary()                      { return fPrimary;                      }
  UInt_t  GetControlBoxSide()                         { return fSide;                         }
  UInt_t  GetControlBoxPrimary()                      { return fPrimary;                      }

  void    SetControlBoxSide(TString bs)               { fSide = bs;                           }
  void    SetControlBoxPrimary(Int_t bp)              { fPrimary = bp;                        }
  void    SetClkLb(UInt_t cl)                         { fClkLb = cl;                          }
  void    SetClkHb(UInt_t ch)                         { fClkHb = ch;                          }

  UInt_t  GetPulseToSmLb()                            { return fPulseToSmLb;                  }
  UInt_t  GetPulseToSmHb()                            { return fPulseToSmHb;                  }
  UInt_t  GetDeadTimeLb()                             { return fDeadTimeLb;                   }
  UInt_t  GetDeadTimeHb()                             { return fDeadTimeHb;                   }
  UInt_t  GetPrePulsesLb()                            { return fPrePulsesLb;                  }
  UInt_t  GetPrePulsesHb()                            { return fPrePulsesHb;                  }
  UInt_t  GetL0CtpLb()                                { return fL0CtpLb;                      }
  UInt_t  GetL0CtpHb()                                { return fL0CtpHb;                      }
  UInt_t  GetL1CtpLb()                                { return fL1CtpLb;                      }
  UInt_t  GetL1CtpHb()                                { return fL1CtpHb;                      }
  UInt_t  GetL0LutLb()                                { return fL0LutLb;                      }
  UInt_t  GetL0LutHb()                                { return fL0LutHb;                      }
  UInt_t  GetNotTinToCtpPreLb()                       { return fNotTinToCtpPreLb;             }
  UInt_t  GetNotTinToCtpPreHb()                       { return fNotTinToCtpPreHb;             }
  UInt_t  GetPatternMatchCbALb()                      { return fPatternMatchCbALb;            }
  UInt_t  GetPatternMatchCbAHb()                      { return fPatternMatchCbAHb;            }
  UInt_t  GetPatternMatchCbCLb()                      { return fPatternMatchCbCLb;            }
  UInt_t  GetPatternMatchCbCHb()                      { return fPatternMatchCbCHb;            }
  UInt_t  GetPatternMatchCbTlmuLb()                   { return fPatternMatchCbTlmuLb;         }
  UInt_t  GetPatternMatchCbTlmuHb()                   { return fPatternMatchCbTlmuHb;         }
  UInt_t  GetTriggerSynchrLb()                        { return fTriggerSynchrLb;              }
  UInt_t  GetTriggerSynchrHb()                        { return fTriggerSynchrHb;              }
  UInt_t  GetTriggerToSmLb()                          { return fTriggerToSmLb;                }
  UInt_t  GetTriggerToSmHb()                          { return fTriggerToSmHb;                }
  UInt_t  GetPreLutLb()                               { return fPreLutLb;                     }
  UInt_t  GetPreLutHb()                               { return fPreLutHb;                     }
  UInt_t  GetMissingPreLb()                           { return fMissingPreLb;                 }
  UInt_t  GetMissingPreHb()                           { return fMissingPreHb;                 }
  UInt_t  GetUnnecessaryPreLb()                       { return fUnnecessaryPreLb;             }
  UInt_t  GetUnnecessaryPreHb()                       { return fUnnecessaryPreHb;             }
  UInt_t  GetPreToNextCycles()                        { return fPreToNextCycles;              }
  UInt_t  GetL0ToNextCycles()                         { return fL0ToNextCycles;               }
  UInt_t  GetL1ToNextCycles()                         { return fL1ToNextCycles;               }
  UInt_t  GetClkTtcexShifts()                         { return fClkTtcexShifts;               }
  UInt_t  GetDisableGtuBusy()                         { return fDisableGtuBusy;               }
  UInt_t  GetDisableSorEorBusy()                      { return fDisableSorEorBusy;            }
  UInt_t  GetScintEn()                                { return fScintEn;                      }
  UInt_t  GetLutAsPre()                               { return fLutAsPre;                     }
  UInt_t  GetCtpAsPre()                               { return fCtpAsPre;                     }
  UInt_t  GetL0En()                                   { return fL0En;                         }
  UInt_t  GetL1En()                                   { return fL1En;                         }
  UInt_t  GetLutL0ToCtpEn()                           { return fLutL0ToCtpEn;                 }
  UInt_t  GetLutPreToCtpEn()                          { return fLutPreToCtpEn;                }
  UInt_t  GetSoftwTrigToCtpEn()                       { return fSoftwTrigToCtpEn;             }
  UInt_t  GetL0Autogenerate()                         { return fL0Autogenerate;               }
  UInt_t  GetL1Autogenerate()                         { return fL1Autogenerate;               }
  UInt_t  GetChannelBDisable()                        { return fChannelBDisable;              }
  UInt_t  GetTtcexClkDisable()                        { return fTtcexClkDisable;              }
  UInt_t  GetPreAsL0En()                              { return fPreAsL0En;                    }
  UInt_t  GetNoTriggerToSm()                          { return fNoTriggerToSm;                }
  UInt_t  GetSoftwTrigEn()                            { return fSoftwTrigEn;                  }
  UInt_t  GetSoftwTrigPattern()                       { return fSoftwTrigPattern;             }
  UInt_t  GetSychrInput()                             { return fSychrInput;                   }
  UInt_t  GetRandomTriggerThr()                       { return fRandomTriggerThr;             }
  UInt_t  GetPtChDelayCbA()                           { return fPtChDelayCbA;                 }
  UInt_t  GetPtChDelayCbC()                           { return fPtChDelayCbC;                 }
  UInt_t  GetPtChDelayGtu()                           { return fPtChDelayGtu;                 }
  UInt_t  GetPtChDelayT0()                            { return fPtChDelayT0;                  }
  UInt_t  GetPtChDelayT0()                            { return fPtChDelayT0;                  }
  UInt_t  GetPtChDelayT1()                            { return fPtChDelayT1;                  }
  UInt_t  GetPtChDelayT2()                            { return fPtChDelayT2;                  }
  UInt_t  GetPtChDelayT3()                            { return fPtChDelayT3;                  }
  UInt_t  GetPtChDelayT4()                            { return fPtChDelayT4;                  }
  UInt_t  GetPtChDelayT5()                            { return fPtChDelayT5;                  }
  UInt_t  GetPtChDelayT6()                            { return fPtChDelayT6;                  }
  UInt_t  GetPtChDelayT7()                            { return fPtChDelayT7;                  }


  void    SetPulseToSmLb(UInt_ ar)                    { fPulseToSmLb = ar;                    }
  void    SetPulseToSmHb(UInt_ ar)                    { fPulseToSmHb = ar;                    }
  void    SetDeadTimeLb(UInt_ ar)                     { fDeadTimeLb = ar;                     }
  void    SetDeadTimeHb(UInt_ ar)                     { fDeadTimeHb = ar;                     }
  void    SetPrePulsesLb(UInt_ ar)                    { fPrePulsesLb = ar;                    }
  void    SetPrePulsesHb(UInt_ ar)                    { fPrePulsesHb = ar;                    }
  void    SetL0CtpLb(UInt_ ar)                        { fL0CtpLb = ar;                        }
  void    SetL0CtpHb(UInt_ ar)                        { fL0CtpHb = ar;                        }
  void    SetL1CtpLb(UInt_ ar)                        { fL1CtpLb = ar;                        }
  void    SetL1CtpHb(UInt_ ar)                        { fL1CtpHb = ar;                        }
  void    SetL0LutLb(UInt_ ar)                        { fL0LutLb = ar;                        }
  void    SetL0LutHb(UInt_ ar)                        { fL0LutHb = ar;                        }
  void    SetNotTinToCtpPreLb(UInt_ ar)               { fNotTinToCtpPreLb = ar;               }
  void    SetNotTinToCtpPreHb(UInt_ ar)               { fNotTinToCtpPreHb = ar;               }
  void    SetPatternMatchCbALb(UInt_ ar)              { fPatternMatchCbALb = ar;              }
  void    SetPatternMatchCbAHb(UInt_ ar)              { fPatternMatchCbAHb = ar;              }
  void    SetPatternMatchCbCLb(UInt_ ar)              { fPatternMatchCbCLb = ar;              }
  void    SetPatternMatchCbCHb(UInt_ ar)              { fPatternMatchCbCHb = ar;              }
  void    SetPatternMatchCbTlmuLb(UInt_ ar)           { fPatternMatchCbTlmuLb = ar;           }
  void    SetPatternMatchCbTlmuHb(UInt_ ar)           { fPatternMatchCbTlmuHb = ar;           }
  void    SetTriggerSynchrLb(UInt_ ar)                { fTriggerSynchrLb = ar;                }
  void    SetTriggerSynchrHb(UInt_ ar)                { fTriggerSynchrHb = ar;                }
  void    SetTriggerToSmLb(UInt_ ar)                  { fTriggerToSmLb = ar;                  }
  void    SetTriggerToSmHb(UInt_ ar)                  { fTriggerToSmHb = ar;                  }
  void    SetPreLutLb(UInt_ ar)                       { fPreLutLb = ar;                       }
  void    SetPreLutHb(UInt_ ar)                       { fPreLutHb = ar;                       }
  void    SetMissingPreLb(UInt_ ar)                   { fMissingPreLb = ar;                   }
  void    SetMissingPreHb(UInt_ ar)                   { fMissingPreHb = ar;                   }
  void    SetUnnecessaryPreLb(UInt_ ar)               { fUnnecessaryPreLb = ar;               }
  void    SetUnnecessaryPreHb(UInt_ ar)               { fUnnecessaryPreHb = ar;               }
  void    SetPreToNextCycles(UInt_ ar)                { fPreToNextCycles = ar;                }
  void    SetL0ToNextCycles(UInt_ ar)                 { fL0ToNextCycles = ar;                 }
  void    SetL1ToNextCycles(UInt_ ar)                 { fL1ToNextCycles = ar;                 }
  void    SetClkTtcexShifts(UInt_ ar)                 { fClkTtcexShifts = ar;                 }
  void    SetDisableGtuBusy(UInt_ ar)                 { fDisableGtuBusy = ar;                 }
  void    SetDisableSorEorBusy(UInt_ ar)              { fDisableSorEorBusy = ar;              }
  void    SetScintEn(UInt_ ar)                        { fScintEn = ar;                        }
  void    SetLutAsPre(UInt_ ar)                       { fLutAsPre = ar;                       }
  void    SetCtpAsPre(UInt_ ar)                       { fCtpAsPre = ar;                       }
  void    SetL0En(UInt_ ar)                           { fL0En = ar;                           }
  void    SetL1En(UInt_ ar)                           { fL1En = ar;                           }
  void    SetLutL0ToCtpEn(UInt_ ar)                   { fLutL0ToCtpEn = ar;                   }
  void    SetLutPreToCtpEn(UInt_ ar)                  { fLutPreToCtpEn = ar;                  }
  void    SetSoftwTrigToCtpEn(UInt_ ar)               { fSoftwTrigToCtpEn = ar;               }
  void    SetL0Autogenerate(UInt_ ar)                 { fL0Autogenerate = ar;                 }
  void    SetL1Autogenerate(UInt_ ar)                 { fL1Autogenerate = ar;                 }
  void    SetChannelBDisable(UInt_ ar)                { fChannelBDisable = ar;                }
  void    SetTtcexClkDisable(UInt_ ar)                { fTtcexClkDisable = ar;                }
  void    SetPreAsL0En(UInt_ ar)                      { fPreAsL0En = ar;                      }
  void    SetNoTriggerToSm(UInt_ ar)                  { fNoTriggerToSm = ar;                  }
  void    SetSoftwTrigEn(UInt_ ar)                    { fSoftwTrigEn = ar;                    }
  void    SetSoftwTrigPattern(UInt_ ar)               { fSoftwTrigPattern = ar;               }
  void    SetSychrInput(UInt_ ar)                     { fSychrInput = ar;                     }
  void    SetRandomTriggerThr(UInt_ ar)               { fRandomTriggerThr = ar;               }
  void    SetPtChDelayCbA(UInt_ ar)                   { fPtChDelayCbA = ar;                   }
  void    SetPtChDelayCbC(UInt_ ar)                   { fPtChDelayCbC = ar;                   }
  void    SetPtChDelayGtu(UInt_ ar)                   { fPtChDelayGtu = ar;                   }
  void    SetPtChDelayT0(UInt_ ar)                    { fPtChDelayT0 = ar;                    }
  void    SetPtChDelayT0(UInt_ ar)                    { fPtChDelayT0 = ar;                    }
  void    SetPtChDelayT1(UInt_ ar)                    { fPtChDelayT1 = ar;                    }
  void    SetPtChDelayT2(UInt_ ar)                    { fPtChDelayT2 = ar;                    }
  void    SetPtChDelayT3(UInt_ ar)                    { fPtChDelayT3 = ar;                    }
  void    SetPtChDelayT4(UInt_ ar)                    { fPtChDelayT4 = ar;                    }
  void    SetPtChDelayT5(UInt_ ar)                    { fPtChDelayT5 = ar;                    }
  void    SetPtChDelayT6(UInt_ ar)                    { fPtChDelayT6 = ar;                    }
  void    SetPtChDelayT7(UInt_ ar)                    { fPtChDelayT7 = ar;                    }


 protected:
  TString fSide; // side of the control box, either A, B or C 
  TInt_t  fPrimary; // 1 if its the primary control box, 2 for backup

  UInt_t  fClkLb;
  UInt_t  fClkHb;

  UInt_t  fPulseToSmLb;
  UInt_t  fPulseToSmHb;
  UInt_t  fDeadTimeLb;
  UInt_t  fDeadTimeHb;
  UInt_t  fPrePulsesLb;
  UInt_t  fPrePulsesHb;
  UInt_t  fL0CtpLb;
  UInt_t  fL0CtpHb;
  UInt_t  fL1CtpLb;
  UInt_t  fL1CtpHb;
  UInt_t  fL0LutLb;
  UInt_t  fL0LutHb;
  UInt_t  fNotTinToCtpPreLb;
  UInt_t  fNotTinToCtpPreHb;
  UInt_t  fPatternMatchCbALb;
  UInt_t  fPatternMatchCbAHb;
  UInt_t  fPatternMatchCbCLb;
  UInt_t  fPatternMatchCbCHb;
  UInt_t  fPatternMatchCbTlmuLb;
  UInt_t  fPatternMatchCbTlmuHb;
  UInt_t  fTriggerSynchrLb;
  UInt_t  fTriggerSynchrHb;
  UInt_t  fTriggerToSmLb;
  UInt_t  fTriggerToSmHb;
  UInt_t  fPreLutLb;
  UInt_t  fPreLutHb;
  UInt_t  fMissingPreLb;
  UInt_t  fMissingPreHb;
  UInt_t  fUnnecessaryPreLb;
  UInt_t  fUnnecessaryPreHb;
  UInt_t  fPreToNextCycles;
  UInt_t  fL0ToNextCycles;
  UInt_t  fL1ToNextCycles;
  UInt_t  fClkTtcexShifts;
  UInt_t  fDisableGtuBusy;
  UInt_t  fDisableSorEorBusy;
  UInt_t  fScintEn;
  UInt_t  fLutAsPre;
  UInt_t  fCtpAsPre;
  UInt_t  fL0En;
  UInt_t  fL1En;
  UInt_t  fLutL0ToCtpEn;
  UInt_t  fLutPreToCtpEn;
  UInt_t  fSoftwTrigToCtpEn;
  UInt_t  fL0Autogenerate;
  UInt_t  fL1Autogenerate;
  UInt_t  fChannelBDisable;
  UInt_t  fTtcexClkDisable;
  UInt_t  fPreAsL0En;
  UInt_t  fNoTriggerToSm;
  UInt_t  fSoftwTrigEn;
  UInt_t  fSoftwTrigPattern;
  UInt_t  fSychrInput;
  UInt_t  fRandomTriggerThr;
  UInt_t  fPtChDelayCbA;
  UInt_t  fPtChDelayCbC;
  UInt_t  fPtChDelayGtu;
  UInt_t  fPtChDelayT0;
  UInt_t  fPtChDelayT0;
  UInt_t  fPtChDelayT1;
  UInt_t  fPtChDelayT2;
  UInt_t  fPtChDelayT3;
  UInt_t  fPtChDelayT4;
  UInt_t  fPtChDelayT5;
  UInt_t  fPtChDelayT6;
  UInt_t  fPtChDelayT7;

  ClassDef(AliTRDCalDCSPTRCbb,1)      //  TRD calibration class for TRD GTU parameters

};
#endif
