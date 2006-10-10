#ifndef ALITRDSEED_H
#define ALITRDSEED_H   

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */ 

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  The TRD track seed                                                    //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TObject.h" 

class AliTRDcluster;

class AliTRDseed : public TObject {

 public:

  AliTRDseed(); 
  AliTRDseed(const AliTRDseed &s);
  ~AliTRDseed() {};                 

  AliTRDseed      &operator=(const AliTRDseed &/*s*/)       { return *this;   } 

  static  void     EvaluateUni(Int_t nvectors, Double_t *data, Double_t &mean, Double_t &sigma, Int_t hh);
  static  Float_t  FitRiemanTilt(AliTRDseed *seed, Bool_t error);
          void     UseClusters();
          void     Update();
          void     CookLabels();
          void     UpdateUsed();
          void     Reset();

          Bool_t   IsOK() const                             { return fN2 > 8;        }
          Bool_t   IsUsable(Int_t i) const                  { return fUsable[i];     }
 
          Float_t  GetTilt() const                          { return fTilt;          }
          Float_t  GetPadLength() const                     { return fPadLength;     }
          Float_t  GetX0() const                            { return fX0;            }
          Float_t  GetX(Int_t i) const                      { return fX[i];          }
          Float_t  GetY(Int_t i) const                      { return fY[i];          }
          Float_t  GetZ(Int_t i) const                      { return fZ[i];          }
          Int_t    GetIndexes(Int_t i) const                { return fIndexes[i];    }
  AliTRDcluster   *GetClusters(Int_t i) const               { return fClusters[i];   }
          Float_t  GetYref(Int_t i) const                   { return fYref[i];       }
          Float_t  GetZref(Int_t i) const                   { return fZref[i];       }
          Float_t  GetYfit(Int_t i) const                   { return fYfit[i];       }
          Float_t  GetYfitR(Int_t i) const                  { return fYfitR[i];      }
          Float_t  GetZfit(Int_t i) const                   { return fZfit[i];       }
          Float_t  GetZfitR(Int_t i) const                  { return fZfitR[i];      }
          Float_t  GetSigmaY() const                        { return fSigmaY;        }
          Float_t  GetSigmaY2() const                       { return fSigmaY2;       }
          Float_t  GetMeanz() const                         { return fMeanz;         }
          Float_t  GetZProb() const                         { return fZProb;         }
          Int_t    GetLabels(Int_t i) const                 { return fLabels[i];     }
          Int_t    GetN2() const                            { return fN2;            }
          Int_t    GetNUsed() const                         { return fNUsed;         }
          Int_t    GetFreq() const                          { return fFreq;          }
          Float_t  GetC() const                             { return fC;             }
          Float_t  GetCC() const                            { return fCC;            }
          Float_t  GetChi2() const                          { return fChi2;          }
          Float_t  GetChi2Z() const                         { return fChi2Z;         }

          void     SetTilt(Float_t tilt)                    { fTilt        = tilt;   }
          void     SetPadLength(Float_t len)                { fPadLength   = len;    }
          void     SetX0(Float_t x0)                        { fX0          = x0;     }
          void     SetX(Int_t i, Float_t x)                 { fX[i]        = x;      } 
          void     SetY(Int_t i, Float_t y)                 { fY[i]        = y;      }
          void     SetZ(Int_t i, Float_t z)                 { fZ[i]        = z;      }
          void     SetIndexes(Int_t i, Int_t idx)           { fIndexes[i]  = idx;    }
          void     SetClusters(Int_t i, AliTRDcluster *c)   { fClusters[i] = c;      }
          void     SetUsable(Int_t i, Bool_t usable)        { fUsable[i]   = usable; }
          void     SetYref(Int_t i, Float_t yref)           { fYref[i]     = yref;   }
          void     SetZref(Int_t i, Float_t zref)           { fZref[i]     = zref;   }
          void     SetYfit(Int_t i, Float_t yfit)           { fYfit[i]     = yfit;   }
          void     SetYfitR(Int_t i, Float_t yfitr)         { fYfitR[i]    = yfitr;  }
          void     SetZfit(Int_t i, Float_t zfit)           { fZfit[i]     = zfit;   }
          void     SetZfitR(Int_t i, Float_t zfitr)         { fZfitR[i]    = zfitr;  }
          void     SetSigmaY(Float_t sigmay)                { fSigmaY      = sigmay; }
          void     SetSigmaY2(Float_t sigmay)               { fSigmaY2     = sigmay; }
          void     SetMeanz(Float_t meanz)                  { fMeanz       = meanz;  }
          void     SetZProb(Float_t zprob)                  { fZProb       = zprob;  }
          void     SetLabels(Int_t i, Int_t label)          { fLabels[i]   = label;  }
          void     SetN2(Int_t n2)                          { fN2          = n2;     }
          void     SetNUsed(Int_t nused)                    { fNUsed       = nused;  }
          void     SetFreq(Int_t freq)                      { fFreq        = freq;   }
          void     SetC(Float_t c)                          { fC           = c;      }
          void     SetCC(Float_t cc)                        { fCC          = cc;     }
          void     SetChi2(Float_t chi2)                    { fChi2        = chi2;   }
          void     SetChi2Z(Float_t chi2z)                  { fChi2Z       = chi2z;  }

 private:

          Float_t  fTilt;               //  Tilting angle
          Float_t  fPadLength;          //  Pad length
          Float_t  fX0;                 //  X0 position
          Float_t  fX[25];              //! X position
          Float_t  fY[25];              //! Y position
          Float_t  fZ[25];              //! Z position
          Int_t    fIndexes[25];        //! Indexes
          AliTRDcluster *fClusters[25]; //! Clusters
          Bool_t   fUsable[25];         //! Indication  - usable cluster
          Float_t  fYref[2];            //  Reference y
          Float_t  fZref[2];            //  Reference z
          Float_t  fYfit[2];            //  Y fit position +derivation
          Float_t  fYfitR[2];           //  Y fit position +derivation
          Float_t  fZfit[2];            //  Z fit position
          Float_t  fZfitR[2];           //  Z fit position
          Float_t  fSigmaY;             //  "Robust" sigma in Y - constant fit
          Float_t  fSigmaY2;            //  "Robust" sigma in Y - line fit
          Float_t  fMeanz;              //  Mean vaue of z
          Float_t  fZProb;              //  Max probbable z
          Int_t    fLabels[2];          //  Labels
          Int_t    fN;                  //  Number of associated clusters
          Int_t    fN2;                 //  Number of not crossed
          Int_t    fNUsed;              //  Number of used clusters
          Int_t    fFreq;               //  Frequency
          Int_t    fNChange;            //  Change z counter
          Float_t  fMPads;              //  Mean number of pads per cluster

          Float_t  fC;                  //  Curvature
          Float_t  fCC;                 //  Curvature with constrain
          Float_t  fChi2;               //  Global chi2
          Float_t  fChi2Z;              //  Global chi2

  ClassDef(AliTRDseed,1)                //  Seed for a local TRD track

};

#endif 
