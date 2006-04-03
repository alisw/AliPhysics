#ifndef ALITRDMCM_H
#define ALITRDMCM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////
//  Multi Chip Module object                         //
///////////////////////////////////////////////////////

#include <TObject.h>

const Int_t kMaxTrackletsPerMCM = 4;
const Int_t kMcmCol = 21;
const Int_t kMcmTBmax = 60;
const Int_t kSelClus = 6;
const Int_t kMaxClus = 4;

class AliTRDtrigParam;

class AliTRDmcm : public TObject {

 public:

  AliTRDmcm();
  AliTRDmcm(AliTRDtrigParam *trigp, const Int_t id);
  virtual ~AliTRDmcm();

  Int_t   Ntrk()                                          { return fNtrk; };
  Int_t  *GetTrkIndex()                                   { return &fTrkIndex[0]; };
  Int_t   GetRobId()                                      { return fRobId; };
  Int_t   GetChaId()                                      { return fChaId; };
  void    SetRobId(const Int_t id)                        { fRobId = id; };
  void    SetChaId(const Int_t id)                        { fChaId = id; };
  void    SetRow(const Int_t row)                         { fRow = row; };
  Int_t   GetRow()                                        { return fRow; };
  void    SetColRange(const Int_t colf, const Int_t coll) { fColFirst = colf; fColLast = coll; };
  void    GetColRange(Int_t &colf, Int_t &coll)           { colf = fColFirst; coll = fColLast; };
  void    AddTrk(const Int_t id);
  void    SetADC(const Int_t icol, const Int_t itb, const Float_t adc) 
    { fADC[icol][itb] = adc; };
  Float_t GetADC(const Int_t icol, const Int_t itb)       { return fADC[icol][itb]; };
  void    Reset();
  Bool_t  Run();
  void    SetCluster(const Int_t icol, const Int_t itb)      { fIsClus[icol][itb] = kTRUE; };
  void    UnSetCluster(const Int_t icol, const Int_t itb)    { fIsClus[icol][itb] = kFALSE; };
  Bool_t  IsCluster(Float_t amp[3]);
  void    AddTimeBin(const Int_t itime);
  Int_t   CreateSeeds();
  void    Sort(const Int_t nel, Int_t *x1, Int_t *x2, Int_t dir);
  Int_t   GetNtrkSeeds()                           { return fNtrkSeeds; };
  Int_t  *GetSeedCol()                             { return &fSeedCol[0]; };
  Int_t  *GetPadHits()                             { return &fPadHits[0]; };
  Bool_t  IsCluster(const Int_t icol, const Int_t itim) { return fIsClus[icol][itim]; };

  void Filter(Int_t nexp, Int_t ftype = 0);
  void DeConvExpA(Double_t *source, Double_t *target, Int_t n, Int_t nexp);
  void DeConvExpD(Double_t *source, Int_t    *target, Int_t n, Int_t nexp);
  void DeConvExpMI(Double_t *source, Double_t *target, Int_t n);
  void TailMakerSpline(Double_t *ampin, Double_t *ampout, Double_t lambda, Int_t n);
  void TailCancelationMI(Double_t *ampin, Double_t *ampout, Double_t norm, Double_t lambda, Int_t n);

  Int_t  Id()                          { return fId; };

 protected:

  AliTRDtrigParam *fTrigParam;                //! pointer to the trigger parameters class

  Int_t   fNtrk;                              //  number of found tracklets
  Int_t   fTrkIndex[kMaxTrackletsPerMCM];     //  index of found tracklets
  Int_t   fRobId;                             //  ROB id
  Int_t   fChaId;                             //  Chamber id
  Int_t   fRow;                               //  pad row number (0-11 or 0-15)
  Int_t   fColFirst;                          //  first pad column
  Int_t   fColLast;                           //  last pad column (<)
  Float_t fADC[kMcmCol][kMcmTBmax];           //! array with MCM ADC values
  Bool_t  fIsClus[kMcmCol][kMcmTBmax];        //! flag of a cluster maximum
  Int_t   fTime1;                             //  first time bin for tracking (incl.)
  Int_t   fTime2;                             //  last  time bin for tracking (incl.)
  Float_t fClusThr;                           //  cluster threshold
  Float_t fPadThr;                            //  pad threshold
  Int_t   fPadHits[kMcmCol];                  //  hit counter in pads
  Int_t   fNtrkSeeds;                         //  number of found seeds
  Int_t   fSeedCol[kMaxTrackletsPerMCM];      //  column number of found tracklet seeds
  Float_t fR1;                                // filter parameters (1 = long, 2 = short component)
  Float_t fR2;                                //
  Float_t fC1;                                //
  Float_t fC2;                                //
  Float_t fPedestal;                          //

  Int_t fId;                                  // dummy id

  ClassDef(AliTRDmcm,1)                       // TRD mcm

};

#endif
