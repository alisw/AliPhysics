#ifndef ALITRDRECONSTRUCTOR_H
#define ALITRDRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class for TRD reconstruction                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliReconstructor.h"
#include "AliDetectorRecoParam.h"
#include "AliTRDrecoParam.h"

class TClonesArray;
class AliRawReader;
class AliTRDReconstructor: public AliReconstructor 
{
public:
  enum AliTRDsteerParam {
    kWriteClusters   = BIT(0)
    ,kSeeding        = BIT(1)
    ,kSteerPID       = BIT(2)
    ,kWriteTracklets = BIT(3)
    ,kDriftGas       = BIT(4)
    ,kHLT            = BIT(5)
  };
  enum AliTRDReconstructorTask {
    kClusterizer = 0
    ,kTracker    = 1
    ,kPID        = 2
  };
  enum AliTRDpidMethod {
    kLQPID = 0,
    kNNPID = 1
  };
  enum AliTRDdriftGas {
    kXe = 0,
    kAr = 1
  };

  AliTRDReconstructor();
  AliTRDReconstructor(const AliTRDReconstructor &r);
  virtual ~AliTRDReconstructor();
  AliTRDReconstructor& operator = (const AliTRDReconstructor&)          { return *this;}
	
	virtual void        Init();

  virtual void        ConvertDigits(AliRawReader *rawReader, TTree *digitsTree) const;
  virtual AliTracker* CreateTracker() const;

  virtual void        FillESD(AliRawReader *, TTree *clusterTree, AliESDEvent *esd) const { FillESD((TTree * )NULL, clusterTree, esd);                    }
  virtual void        FillESD(TTree *digitsTree, TTree *clusterTree, AliESDEvent *esd) const;
  static TClonesArray* GetClusters() {return fgClusters;}
  Int_t               GetNdEdxSlices() const     { return GetPIDMethod() == kNNPID ? kNNslices : kLQslices;}
  AliTRDdriftGas      GetDriftGas() const        { return fSteerParam&kDriftGas ? kAr : kXe;}
  AliTRDpidMethod     GetPIDMethod() const       { return fSteerParam&kSteerPID ? kNNPID : kLQPID;}
  static const AliTRDrecoParam* GetRecoParam() { return dynamic_cast<const AliTRDrecoParam*>(AliReconstructor::GetRecoParam(2)); }
  Int_t               GetStreamLevel(AliTRDReconstructorTask task) const    { return fStreamLevel[task];} 
  inline void         GetTCParams(Double_t *par) const;
  virtual Bool_t      HasDigitConversion() const                   { return kFALSE;           };
  Bool_t              IsWritingClusters() const  { return fSteerParam&kWriteClusters;}
  Bool_t              IsWritingTracklets() const { return fSteerParam&kWriteTracklets;}
  Bool_t              IsHLT() const              { return fSteerParam&kHLT;}
  Bool_t              IsSeeding() const          { return fSteerParam&kSeeding;}

  virtual void        Reconstruct(AliRawReader *rawReader, TTree *clusterTree) const;
  virtual void        Reconstruct(TTree *digitsTree, TTree *clusterTree) const;

  static void         SetClusters(TClonesArray *clusters) {fgClusters = clusters;}
  void	              SetOption(Option_t *opt);
  inline void         SetTCParams(Double_t *par);
  void                SetStreamLevel(Int_t level, AliTRDReconstructorTask task= kTracker);

private:
  enum{
    kNNslices = 8
   ,kLQslices = 3
  };
  UChar_t       fStreamLevel[5];      // stream level for each reconstruction task         
  UInt_t        fSteerParam;          // steering flags
  Double_t      fTCParams[8];         // Tail Cancellation parameters for drift gases 
 
  static TClonesArray *fgClusters;    // list of clusters for local reconstructor

  ClassDef(AliTRDReconstructor, 1)    //  Class for the TRD reconstruction

};

//___________________________________________________
inline void AliTRDReconstructor::GetTCParams(Double_t *par) const
{
  if(!par) return;
  if(GetDriftGas()==kAr) memcpy(par, &fTCParams[4], 4*sizeof(Double_t));
  else memcpy(par, &fTCParams[0], 4*sizeof(Double_t));
}

//___________________________________________________
inline void AliTRDReconstructor::SetTCParams(Double_t *par)
{
  if(!par) return;
  memcpy(fTCParams, par, 8*sizeof(Double_t));
}

#endif
