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
#include "AliTRDpidUtil.h"

class TClonesArray;
class TTreeSRedirector;
class AliRawReader;
class AliTRDReconstructor: public AliReconstructor 
{
public:
  enum ETRDReconstructorSteer {
    kDigitsConversion= BIT(0)
    ,kTC             = BIT(1) // tail cancelation
    ,kLUT            = BIT(2) // look up table for cluster position determination 
    ,kGAUS           = BIT(3) // look up table for cluster position determination 
    ,kClusterSharing = BIT(4) // Toggle cluster sharing
    ,kSteerPID       = BIT(5)
    ,kEightSlices    = BIT(6)
    ,kWriteClusters  = BIT(7)
    ,kWriteTracklets = BIT(8)
    ,kDriftGas       = BIT(9)
    ,kSeeding        = BIT(10)
    ,kVertexConstrained = BIT(11) // Perform vertex constrained fit
    ,kImproveTracklet   = BIT(12) // Improve tracklet in the SA TRD track finder 
    ,kHLT            = BIT(13)
    ,kCosmic         = BIT(14)
    ,kOwner          = BIT(15)
    ,kNsteer         = 15       // number of tasks
  };
  enum ETRDReconstructorTask {
    kRawReader    = 0
    ,kClusterizer = 1
    ,kTracker     = 2
    ,kPID         = 3
    ,kNtasks      = 4  // number of reconsruction tasks
  };
  enum ETRDReconstructorGas {
    kXe = 0
    ,kAr = 1
  };

  AliTRDReconstructor();
  AliTRDReconstructor(const AliTRDReconstructor &r);
  virtual ~AliTRDReconstructor();
  AliTRDReconstructor& operator = (const AliTRDReconstructor&)          { return *this;}
	
	virtual void        Init();

  virtual void        ConvertDigits(AliRawReader *rawReader, TTree *digitsTree) const;
  virtual AliTracker* CreateTracker() const;
  TTreeSRedirector*   GetDebugStream(ETRDReconstructorTask task) const { return task < kNtasks ? fDebugStream[task] : 0x0; }

  virtual void        FillESD(AliRawReader *, TTree *clusterTree, AliESDEvent *esd) const { FillESD((TTree * )NULL, clusterTree, esd);                    }
  virtual void        FillESD(TTree *digitsTree, TTree *clusterTree, AliESDEvent *esd) const;
  static TClonesArray* GetClusters() {return fgClusters;}
  Int_t               GetNdEdxSlices() const     { return (Int_t)AliTRDpidUtil::GetNdEdxSlices(GetPIDMethod());}
  ETRDReconstructorGas GetDriftGas() const        { return fSteerParam&kDriftGas ? kAr : kXe;}
  AliTRDpidUtil::ETRDPIDMethod       GetPIDMethod() const       { return fSteerParam&kSteerPID ? AliTRDpidUtil::kNN : AliTRDpidUtil::kLQ;}
  static const AliTRDrecoParam* GetRecoParam() { return dynamic_cast<const AliTRDrecoParam*>(AliReconstructor::GetRecoParam(2)); }
  Int_t               GetStreamLevel(ETRDReconstructorTask task) const    { return fStreamLevel[task];} 
  inline void         GetTCParams(Double_t *par) const;
  virtual Bool_t      HasDigitConversion() const { return fSteerParam&kDigitsConversion;  };
  Bool_t              HasVertexConstrained() const { return fSteerParam&kVertexConstrained; }
  Bool_t              HasImproveTracklets() const  { return fSteerParam&kImproveTracklet; }
  Bool_t              IsWritingClusters() const  { return fSteerParam&kWriteClusters;}
  Bool_t              IsWritingTracklets() const { return fSteerParam&kWriteTracklets;}
  Bool_t              IsHLT() const              { return fSteerParam&kHLT;}
  Bool_t              IsSeeding() const          { return fSteerParam&kSeeding;}
  Bool_t              IsCosmic() const           { return fSteerParam&kCosmic;}
  Bool_t              IsEightSlices() const      { return fSteerParam&kEightSlices;}
  Bool_t              UseClusterSharing() const  { return fSteerParam&kClusterSharing;}
  Bool_t              UseLUT() const             { return fSteerParam&kLUT;}
  Bool_t              UseGAUS() const             { return fSteerParam&kGAUS;}
  Bool_t              UseTailCancelation() const { return fSteerParam&kTC;}
  
  static void         Options(UInt_t steer=0, UChar_t *stream=0x0);
  virtual void        Reconstruct(AliRawReader *rawReader, TTree *clusterTree) const;
  virtual void        Reconstruct(TTree *digitsTree, TTree *clusterTree) const;

  static void         SetClusters(TClonesArray *clusters) {fgClusters = clusters;}
  void	              SetOption(Option_t *opt);
  inline void         SetTCParams(Double_t *par);
  void                SetStreamLevel(Int_t level, ETRDReconstructorTask task= kTracker);

private:
  static Char_t    *fgSteerNames[kNsteer];//! steering names
  static Char_t    *fgSteerFlags[kNsteer];//! steering flags
  static Char_t    *fgTaskNames[kNtasks]; //! tasks names
  static Char_t    *fgTaskFlags[kNtasks]; //! tasks flags
  UChar_t           fStreamLevel[kNtasks];// stream level for each reconstruction task
  UInt_t            fSteerParam;          // steering bits
  Double_t          fTCParams[8];         // Tail Cancellation parameters for drift gases 
  TTreeSRedirector *fDebugStream[kNtasks];// Debug Streamer container;
 
  static TClonesArray *fgClusters;    // list of clusters for local reconstructor

  ClassDef(AliTRDReconstructor, 2)    //  Class for the TRD reconstruction

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
