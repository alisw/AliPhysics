#ifndef ALITRDPIDLQ_H
#define ALITRDPIDLQ_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliTRDpid.h"

class AliTRDpidLQ : public AliTRDpid {

 public:

  AliTRDpidLQ();
  AliTRDpidLQ(const char* name, const char* title);
  AliTRDpidLQ(const AliTRDpidLQ &p);
  virtual ~AliTRDpidLQ();
  AliTRDpidLQ &operator=(const AliTRDpidLQ &p);

  virtual void          Copy(TObject &p);
  virtual Bool_t        Init();
  virtual Bool_t        AssignLikelihood(AliTRDtrack *t);
  virtual Bool_t        CreateHistograms(const Int_t   nmom
                                       , const Float_t minmom
                                       , const Float_t maxmom);
  virtual Bool_t        FillSpectra(const AliTRDtrack *t);

  Int_t         GetIndex(const AliTRDtrack *t);
  Int_t         GetIndex(const Int_t imom, const Int_t ipid);
  Int_t         GetIndex(const Float_t mom, const Int_t ipid);

          TObjArray*    GetHist() const                     { return fHist;        };

          Float_t       GetChargeMin() const                { return fChargeMin;   };
          Int_t         GetNClusterMin() const              { return fNClusterMin; };

          Int_t         GetNLh() const                      { return fNLh;         };
          Float_t       GetMinLh() const                    { return fMinLh;       };
          Float_t       GetMaxLh() const                    { return fMaxLh;       };

          void          SetChargeMin(const Float_t min)     { fChargeMin   = min;  };
          void          SetNClusterMin(const Int_t min)     { fNClusterMin = min;  };

          void          SetNLh(const Int_t n)               { fNLh         = n;    };
          void          SetMinLh(const Float_t min)         { fMinLh       = min;  };
          void          SetMaxLh(const Float_t max)         { fMaxLh       = max;  };

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
