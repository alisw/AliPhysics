#ifndef ALITRDPIDLQ_H
#define ALITRDPIDLQ_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//   The TRD particle identification class                                   //
//                                                                           //
//   Its main purposes are:                                                  //
//      - Creation and bookkeeping of the propability distributions          //
//      - Assignment of a e/pi propability to a given track based on         //
//        the LQ method                                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDpid.h"

class AliTRDpidLQ : public AliTRDpid {

 public:

  AliTRDpidLQ();
  AliTRDpidLQ(const char* name, const char* title);
  AliTRDpidLQ(const AliTRDpidLQ &p);
  virtual ~AliTRDpidLQ();
  AliTRDpidLQ &operator=(const AliTRDpidLQ &p);

  virtual void          Copy(TObject &p) const;
  virtual Bool_t        Init();
  virtual Bool_t        AssignLikelihood()            { return 0; };
  virtual Bool_t        AssignLikelihood(TObjArray *) { return 0; };
  virtual Bool_t        AssignLikelihood(AliTRDtrack *t);
  virtual Bool_t        CreateHistograms(Int_t nmom, Float_t minmom, Float_t maxmom);
  virtual Bool_t        FillSpectra()                 { return 0; };
  virtual Bool_t        FillSpectra(TObjArray*)       { return 0; };
  virtual Bool_t        FillSpectra(const AliTRDtrack *t);

          Int_t         GetIndex(const AliTRDtrack *t);
          Int_t         GetIndex(Int_t imom, Int_t ipid);
          Int_t         GetIndex(Float_t mom, Int_t ipid);

          TObjArray*    GetHist() const                     { return fHist;        };

          Float_t       GetChargeMin() const                { return fChargeMin;   };
          Int_t         GetNClusterMin() const              { return fNClusterMin; };

          Int_t         GetNLh() const                      { return fNLh;         };
          Float_t       GetMinLh() const                    { return fMinLh;       };
          Float_t       GetMaxLh() const                    { return fMaxLh;       };

          void          SetChargeMin(Float_t min)     { fChargeMin   = min;  };
          void          SetNClusterMin(Int_t min)     { fNClusterMin = min;  };

          void          SetNLh(Int_t n)               { fNLh         = n;    };
          void          SetMinLh(Float_t min)         { fMinLh       = min;  };
          void          SetMaxLh(Float_t max)         { fMaxLh       = max;  };

 protected:

  Int_t           fNMom;             //  Number of momentum bins
  Float_t         fMinMom;           //  Lower momentum
  Float_t         fMaxMom;           //  Upper momentum
  Float_t         fWidMom;           //  Width of the momentum bins
  TObjArray      *fHist;             //  Array of histograms

  Int_t           fNLh;              //  Number of bins of the likelihood spectra
  Float_t         fMinLh;            //  Lower range of the likelihood spectra
  Float_t         fMaxLh;            //  Upper range of the likelihood spectra

  Float_t         fChargeMin;        //  Minimum charge required in one plane
  Int_t           fNClusterMin;      //  Minimum number of clusters required in one plane

  ClassDef(AliTRDpidLQ,1)            //  Assigns e/pi propability to the tracks based on LQ method 

};
#endif
