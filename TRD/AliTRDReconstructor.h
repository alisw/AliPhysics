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
#include "AliRecoParam.h"
#include "AliDetectorRecoParam.h"
#include "AliTRDpidUtil.h"
#include "AliTRDrecoParam.h"

class TClonesArray;
class TTreeSRedirector;
class AliRawReader;
class AliTRDReconstructor: public AliReconstructor 
{
public:
  enum ETRDReconstructorSteer {
    kDigitsConversion= BIT(0)
    ,kWriteClusters  = BIT(1)
    ,kWriteTracklets = BIT(2)
    ,kSeeding        = BIT(3)
    ,kHLT            = BIT(4)
    ,kProcTracklets  = BIT(5) // process online tracklets
    ,kDebug          = BIT(6)
    ,kOwner          = BIT(7)
    ,kNsteer         = 7       // number of tasks
  };
  AliTRDReconstructor();
  virtual ~AliTRDReconstructor();
	
	virtual void        Init();

  virtual void        ConvertDigits(AliRawReader *rawReader, TTree *digitsTree) const;
  virtual AliTracker* CreateTracker() const;
  TTreeSRedirector*   GetDebugStream(AliTRDrecoParam::ETRDReconstructionTask task) const { return task < AliTRDrecoParam::kTRDreconstructionTasks ? fDebugStream[task] : 0x0; }

  virtual void        FillESD(AliRawReader *, TTree *clusterTree, AliESDEvent *esd) const { FillESD((TTree * )NULL, clusterTree, esd);                    }
  virtual void        FillESD(TTree *digitsTree, TTree *clusterTree, AliESDEvent *esd) const;
  static TClonesArray* GetClusters() {return fgClusters;}
  static TClonesArray* GetTracklets() {return fgTracklets;}
  Int_t               GetNdEdxSlices() const     { return (Int_t)AliTRDpidUtil::GetNdEdxSlices(GetPIDMethod());}
  AliTRDpidUtil::ETRDPIDMethod       GetPIDMethod() const       { return GetRecoParam()->IsPIDNeuralNetwork() ? AliTRDpidUtil::kNN : AliTRDpidUtil::kLQ;}
  static const AliTRDrecoParam* GetRecoParam() { return dynamic_cast<const AliTRDrecoParam*>(AliReconstructor::GetRecoParam(2)); }
  virtual Bool_t      HasDigitConversion() const { return fSteerParam&kDigitsConversion;  };
  Bool_t              IsCosmic() const { return GetRecoParam()->GetEventSpecie() & AliRecoParam::kCosmic;}
  Bool_t              IsWritingClusters() const  { return fSteerParam&kWriteClusters;}
  Bool_t              IsWritingTracklets() const { return fSteerParam&kWriteTracklets;}
  Bool_t              IsHLT() const              { return fSteerParam&kHLT;}
  Bool_t              IsSeeding() const          { return fSteerParam&kSeeding;}
  Bool_t              IsProcessingTracklets() const { return fSteerParam&kProcTracklets;}
  Bool_t              IsDebugStreaming() const { return fSteerParam&kDebug;}

  static void         Options(UInt_t steer=0);
  virtual void        Reconstruct(AliRawReader *rawReader, TTree *clusterTree) const;
  virtual void        Reconstruct(TTree *digitsTree, TTree *clusterTree) const;

  static void         SetClusters(TClonesArray *clusters) {fgClusters = clusters;} 
  static void         SetTracklets(TClonesArray *tracklets) {fgTracklets = tracklets;}
  void	              SetOption(Option_t *opt);

private:
  AliTRDReconstructor(const AliTRDReconstructor &r); //Not implemented
  AliTRDReconstructor& operator = (const AliTRDReconstructor&); //Not implemented

  static Char_t    *fgSteerNames[kNsteer];//! steering names
  static Char_t    *fgSteerFlags[kNsteer];//! steering flags
  static Char_t    *fgTaskNames[AliTRDrecoParam::kTRDreconstructionTasks]; //! tasks names
  static Char_t    *fgTaskFlags[AliTRDrecoParam::kTRDreconstructionTasks]; //! tasks flags
  UInt_t            fSteerParam;          // steering bits
  TTreeSRedirector *fDebugStream[AliTRDrecoParam::kTRDreconstructionTasks];// Debug Streamer container;
 
  static TClonesArray *fgClusters;    // list of clusters for local reconstructor
  static TClonesArray *fgTracklets;   // list of online tracklets for local reconstructor

  ClassDef(AliTRDReconstructor, 2)    //  Class for the TRD reconstruction

};



#endif
