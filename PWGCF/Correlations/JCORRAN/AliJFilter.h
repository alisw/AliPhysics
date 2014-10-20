// $Id: AliJFilter.h,v 1.5 2012/04/19 15:19:52 jkral Exp $

//______________________________________________________________________________
// Analysis task for high pt particle correlations 
// author: R.Diaz, J. Rak,  D.J. Kim
// ALICE Group University of Jyvaskyla 
// Finland 
//
// Fill the analysis containers for ESD or AOD
// Note: Adapted for AliAnalysisTaskSE
//////////////////////////////////////////////////////////////////////////////

#ifndef ALIJFILTER_H
#define ALIJFILTER_H

#include "TNamed.h"
#include "AliJRunHeader.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAnalysisTaskSE.h"
#include <iostream>

#include <AliJConst.h>
#include <TVectorT.h>

//==============================================================

#ifndef AliJMaxDimBuffer
#define AliJMaxDimBuffer
const int kMaxDimBuffer = 300;//max length of a line read to a buffe
#endif

class AliJEventHeader;
class AliJRunHeader;
class AliJTrack;
class AliAnalysisTaskSE;

class TH1D;
class TH2D;
class TNtuple;
class TList;
class TTree;
class TFormula;
class TRefArray;
class TArrayI;

class AliMCEvent; 
class AliESDEvent; 
class AliESDtrack;
class AliAODEvent; 
class AliAODTrack; 
class AliESDtrackCuts;
class AliESDVZERO;
class AliESDCentrality;
class AliVCluster;
class AliVCaloCells;
class AliVEvent;

class AliEMCALGeometry;
class AliEMCALGeoUtils;
class AliEMCALRecoUtils;
class AliPHOSGeoUtils;

class AliMCEvent;
class AliAnalysisFilter;

class AliESDTZERO;
class AliESDZDC;
class AliPIDResponse;
class AliPIDResponse;
class AliPIDCombined;
class AliESDTZERO;
class AliAnalysisUtils;

using namespace std;

class AliJFilter : public TNamed  {

 public:
  AliJFilter();
  AliJFilter(const char *name,  AliAnalysisTaskSE *task);
  AliJFilter(const AliJFilter& ap);   
  AliJFilter& operator = (const AliJFilter& ap);
  virtual ~AliJFilter();

  // methods to fill from AliAnalysisTaskSE
  virtual void UserCreateOutputObjects(); 
  virtual void Init();   
  virtual void LocalInit() { Init(); }
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t * opt = "");


  void SetTrackFilter( AliAnalysisFilter * filter ){ fESDFilter = filter; }
  bool GetEventSuccess() const { return fEventSuccess; }
  //are ESDs from real  exp or MonteCarlo 
  //  AliEMCALGeoUtils* GetEMCALGeoUtils (bool doDelete=kFALSE);
  void SetClusterThreshold(double t) { fClusterThreshold = t; }
  void SetTrackThreshold(double t) { fTrackThreshold = t; }
  void SetMyTask( AliAnalysisTaskSE *t ) { fMyTask = t; }
  TClonesArray *GetTrackList() const { return fTrackList; }
  TClonesArray *GetPhotonList() const { return fPhotonList; }
  TClonesArray *GetCaloCellList() const { return fCaloCellList; }
  TClonesArray *GetMCTrackList() const { return fMCTrackList; }
  TClonesArray *GetHeaderList() const { return fHeaderList; }
  TList        *GetRunInfoList() const { return fRunInfoList; }

  AliESDVZERO  *GetESDVZERO() const { return fVZEROData; }
  AliESDTZERO  *GetESDTZERO() const { return fTZEROData; }
  //  AliESDFMD*          fFMDData;
  AliESDZDC* GetESDZDC() const { return fZDCData; }
 
  TClonesArray **GetTrackListP()  { return &fTrackList; }
  TClonesArray **GetPhotonListP()  { return &fPhotonList; }
  TClonesArray **GetCaloCellListP()  { return &fCaloCellList; }
  TClonesArray **GetMCTrackListP()  { return &fMCTrackList; }
  TClonesArray **GetHeaderListP() { return &fHeaderList; }
  TList **GetRunInfoListP()  { return &fRunInfoList; }

  AliESDVZERO** GetESDVZEROP() { return &fVZEROData; }
  AliESDTZERO** GetESDTZEROP() { return &fTZEROData; }
  //  AliESDFMD*          fFMDData;
  AliESDZDC** GetESDZDCP() { return &fZDCData; }

  Bool_t      GetStoreEventPlaneSource(){ return fAliJRunHeader->GetStoreEventPlaneSource(); }
  Bool_t      GetStoreEMCalInfo(){ return fAliJRunHeader->GetStoreEMCalInfo(); }
  AliESDEvent * ESDEvent(){ return FromESD()? dynamic_cast<AliESDEvent*>(Event()):NULL;}
  AliAODEvent * AODEvent(){ return FromAOD()? dynamic_cast<AliAODEvent*>(Event()):NULL;}
  AliVEvent   * Event(){ return fMyTask->InputEvent(); }
  AliMCEvent  * MCEvent(){ return IsMC()?fMyTask->MCEvent():NULL; }

  Bool_t       IsMC(){ return fAliJRunHeader->IsMC(); }
  Bool_t       FromESD(){ return fAliJRunHeader->FromESD(); }
  Bool_t       FromAOD(){ return fAliJRunHeader->FromAOD(); }
 
  AliJRunHeader* GetAliJRunHeader() const { return fAliJRunHeader; }
  void    SetAliJRunHeader( AliJRunHeader* header ){ fAliJRunHeader=header; } 
 private:

  Int_t        DebugLevel(){ return fMyTask->DebugLevel(); }
  inline void   DEBUG(int level, int type, TString msg1, TString msg2=""){
    if(DebugLevel()>level) std::cout<<type<<"\t"<<msg1<<" : "<<msg2<<std::endl;
  }

  AliJEventHeader* ReadCommonHeader(AliVEvent *event);
  // methods to read data from ESD
  void ReadESDTracks(AliESDEvent* esd);
  void ConvertESDTPCOnlyTracks(AliESDEvent* esd, int iTrack, AliJTrack * ctrack, double ptMin, double ptMax);
  void ConvertESDGCGTracks(AliESDEvent* esd, int iTrack, AliJTrack * ctrack, double ptMin, double ptMax);
  void ReadESDCaloClusters(const AliESDEvent* esd);
  void ReadESDCaloCells(const AliESDEvent* esd);
  void ReadESDHeader(AliESDEvent* esd);
  void ReadESDPID(AliESDtrack* track, AliJTrack* ctrack);
  // methods to read data from AOD
  Bool_t ReadAODTracks(const AliAODEvent* aod);
  Bool_t ReadAODCaloClusters(const AliAODEvent* aod);
  void ReadAODCaloCells(const AliAODEvent* aod);
  void ReadAODHeader(AliAODEvent* aod);
  void ReadFilter();
  void ReadMCTracksFromESD();
  void ReadMCTracksFromAOD();
	void RemapMCLabels();

  Int_t GetSuperModuleNumber(bool isemcal, AliVCluster *cluster, AliVCaloCells *cells, Int_t absId);
  Double_t* GetCellsAmplitude( bool isemcal, AliVCluster *cluster, AliVCaloCells *emCells, AliVCaloCells *phoCells );

  UInt_t ConvertTriggerMask();//Converts alice trigger mask to JCorran trigger mask
  //functions used for event selction:
  bool AcceptAODTrack(AliAODTrack* aodTrack);
  void SetOADBPath(const char* path) {fOADBPath=path;}
  const char* GetOADBPath() const { return fOADBPath.Data(); }

  // method to fill jcorran
  void PrintOut() const;
  
  // UTILS
  void AddList(const char* aname, const char* cname, TClonesArray **obj, int nlist);

  // d a t a     m e m b e r s
  AliESDtrackCuts* fEsdTrackCuts; // track selection cuts
  AliAnalysisFilter * fESDFilter; // filter set of track selection BS
  TVectorT<double>  fIsRealOrMC; // flags if the input are real (0) ESDs or MonteCarlo ESDs (1)
  TString fActiveTriggers[kRangeTriggerTableAlice]; // alice table mapping trigger bit to trigger name
  TString fTriggerTableJCorran[kRangeTriggerTableJCorran]; // JCorran trigger table TBit 0 =MinBias
  Bool_t fStoreEventPlaneSource; // store event plane
  TString fOADBPath; // oadb path
  TRefArray *fCaloClustersArr; //! calo cluster array
  Double_t fClusterThreshold; // for event tropping
  Double_t fTrackThreshold; // for event tropping
  Bool_t fEventSuccess; //! if filter was successful with current event

  TArrayI *fMcMap; //! mc index map

  // jcorran output objects
  TClonesArray *    fTrackList;   //! list of charged track objects
  TClonesArray *    fMCTrackList; //! list of charged track objects
  TClonesArray *    fPhotonList;  //! list of photons objects
  TClonesArray *    fCaloCellList;  //! list of calo cells
  TClonesArray *    fHeaderList;  //! event details
  TList *        fRunInfoList; //! run details
  AliPIDResponse  *fPIDResponse; //! PID response object
  AliPIDCombined  *fPIDCombined; //! PID response object

  AliESDVZERO*        fVZEROData;  //!
  AliESDTZERO*        fTZEROData;  //!
  //  AliESDFMD*          fFMDData;
  AliESDZDC*          fZDCData;  //!

  vector<Int_t>       fEMCLabels; //! EMCal hit labels
  vector<Int_t>       fEMCTreeLabels; //! cascades for EMCal hits

  AliJRunHeader*      fAliJRunHeader; //!  run details (mg field, trigger mask,etc...)
  AliEMCALGeometry * fEMCALGeometry;  //! emcal geometry
  AliEMCALRecoUtils * fEMCALRecoUtils;  //! reco utils
  AliPHOSGeoUtils  * fPHOSGeom; //! phos geometry matrix 
  AliAnalysisUtils  *fAnaUtils; //! analysis utils ALICE
  AliAnalysisTaskSE *fMyTask; //! task pointer


  ClassDef(AliJFilter, 1); 
};
#endif // AliJFilter_H
