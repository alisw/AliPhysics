#ifndef ALITRDMCM_H
#define ALITRDMCM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDmcm.h 22653 2007-11-29 21:18:30Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Multi Chip Module object                                              //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliTRDmcm : public TObject {

 public:

  enum { kMaxTrackletsPerMCM = 4
       , kMcmCol             = 21
       , kMcmTBmax           = 60
       , kSelClus            = 6
       , kMaxClus            = 4 };

  AliTRDmcm();
  AliTRDmcm(const AliTRDmcm &m);
  AliTRDmcm(Int_t id);
  virtual         ~AliTRDmcm();
  AliTRDmcm       &operator=(const AliTRDmcm &m);

  virtual void     Copy(TObject &m) const;

          void     SetRobId(Int_t id)                          { fRobId = id; };
          void     SetChaId(Int_t id)                          { fChaId = id; };
          void     SetRow(Int_t row)                           { fRow = row;  };
          void     SetColRange(Int_t colf, Int_t coll)         { fColFirst = colf; fColLast = coll; };
          void     SetADC(Int_t icol, Int_t itb, Float_t adc)  { fADC[icol][itb] = adc;       };
          void     SetCluster(Int_t icol, Int_t itb)           { fIsClus[icol][itb] = kTRUE;  };
          void     UnSetCluster(Int_t icol, Int_t itb)         { fIsClus[icol][itb] = kFALSE; };

          Int_t   *GetTrkIndex()                               { return &fTrkIndex[0];        };
          Int_t    GetRobId() const { return fRobId; };
          Int_t    GetChaId() const { return fChaId; };
          Int_t    GetRow() const { return fRow; };
          void     GetColRange(Int_t &colf, Int_t &coll) const { colf = fColFirst; coll = fColLast; };
          Float_t  GetADC(Int_t icol, Int_t itb) const         { return fADC[icol][itb];     };
          Int_t    GetNtrkSeeds() const                        { return fNtrkSeeds;          };
          Int_t   *GetSeedCol()                                { return &fSeedCol[0];        };
          Int_t   *GetPadHits()                                { return &fPadHits[0];        };
          Bool_t   IsCluster(Int_t icol, Int_t itim) const     { return fIsClus[icol][itim]; };
          Int_t    Ntrk() const { return fNtrk; };
          Int_t    Id() const { return fId; };

          void     AddTrk(Int_t id);
          void     Reset();
          Bool_t   Run();
          Bool_t   IsCluster(Float_t amp[3]) const;
          void     AddTimeBin(Int_t itime);
          Int_t    CreateSeeds();
          void     Sort(Int_t nel, Int_t *x1, Int_t *x2, Int_t dir) const;

          void     Filter(Int_t nexp, Int_t ftype = 0);
          void     DeConvExpA(Double_t *source, Double_t *target, Int_t n, Int_t nexp);
          void     DeConvExpD(Double_t *source, Int_t    *target, Int_t n, Int_t nexp);
          void     DeConvExpMI(Double_t *source, Double_t *target, Int_t n);
          void     TailMakerSpline(Double_t *ampin, Double_t *ampout, Double_t lambda, Int_t n);
          void     TailCancelationMI(Double_t *ampin, Double_t *ampout, Double_t norm, Double_t lambda, Int_t n);

 protected:

          Int_t    fNtrk;                              //  Number of found tracklets
          Int_t    fTrkIndex[kMaxTrackletsPerMCM];     //  Index of found tracklets
          Int_t    fRobId;                             //  ROB id
          Int_t    fChaId;                             //  Chamber id
          Int_t    fRow;                               //  Pad row number (0-11 or 0-15)
          Int_t    fColFirst;                          //  First pad column
          Int_t    fColLast;                           //  Last pad column (<)
          Float_t  fADC[kMcmCol][kMcmTBmax];           //! Array with MCM ADC values
          Bool_t   fIsClus[kMcmCol][kMcmTBmax];        //! Flag of a cluster maximum
          Int_t    fTime1;                             //  First time bin for tracking (incl.)
          Int_t    fTime2;                             //  Last  time bin for tracking (incl.)
          Float_t  fClusThr;                           //  Cluster threshold
          Float_t  fPadThr;                            //  Pad threshold
          Int_t    fPadHits[kMcmCol];                  //  Hit counter in pads
          Int_t    fNtrkSeeds;                         //  Number of found seeds
          Int_t    fSeedCol[kMaxTrackletsPerMCM];      //  Column number of found tracklet seeds
                                                       //  Filter parameters (1 = long, 2 = short component)
          Float_t  fR1;                                //  Time constant [microseconds] 
          Float_t  fR2;                                //  Time constant [microseconds]
          Float_t  fC1;                                //  Weight
          Float_t  fC2;                                //  Weight
          Float_t  fPedestal;                          //  ADC baseline (pedestal)

          Int_t    fId;                                //  Dummy id
 
  ClassDef(AliTRDmcm,3)                                //  TRD MCM class

};

#endif
