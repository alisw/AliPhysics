#ifndef ALITRDCALIBRATION_H
#define ALITRDCALIBRATION_H

// macro for extremely simple analysis

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDrecoTask.h"

class TList;
class TObject;
class TH1F;
class TProfile2D;
class TGraph;
class TH2I;
class TH2D;
class TTree;
class TObjArray;
class AliTRDtrackV1;
class AliTRDCalibraFillHisto;
class AliTRDCalibraVector;
class AliTRDCalibraMode;
class AliTRDcluster;
class AliTRDtrackInfo;

class AliTRDcalibration : public AliTRDrecoTask 
{
 
public:

  // Plots registered for this task
  enum{
    kNbTrack        =  0     // Nb tracks per event
      ,kNbTracklet     =  1     // Nb of traklets per detector
      ,kNbTimeBin      =  2     // Nb of clusters per timebin
      ,kNbClusters     =  3     // Nb of clusters per tracklet
      ,kPHSum          =  4     // All integrated PH
      ,kCHSum          =  5     // All integrated CH
      ,kPH2D           =  6     // PH2D
      ,kCH2D           =  7     // CH2D
      ,kPRF2D          =  8     // PRF2D 
      ,kPH2DVector     =  9     // PH2D
      ,kCH2DVector     =  10     // CH2D
      ,kPRF2DVector    =  11     // PRF2D 
      ,kLinearFitter   =  12    // For the one with the most stats 
      ,kGainFactor     =  13    // Gain factor
      ,kVdriftT0Factor   = 14   // VdriftT0 average pulse height
      ,kVdriftLorentzAngleFactor = 15  // VdriftLorentzAngle
      ,kPRFFactor = 16                 //PRF Factor
      };
  

  AliTRDcalibration();
  AliTRDcalibration(char* name);
  virtual ~AliTRDcalibration();
  
  virtual void    UserCreateOutputObjects();
  virtual void    UserExec(Option_t *option);
  virtual void    Terminate(Option_t *);
  virtual Bool_t  GetRefFigure(Int_t ifig);
  virtual Bool_t  PostProcess();

  Bool_t FillGraphIndex(const TObjArray *vectora, TGraph *graph) const;
  Bool_t AddStatsPerDetector(const TH2I *ch);
  Bool_t AddStatsPerDetector(const TProfile2D *ph, Int_t i);

  Bool_t SetNrphiFromTObject(const char *name, Int_t i, AliTRDCalibraMode *calibMode) const;
  Bool_t SetNzFromTObject(const char *name, Int_t i, AliTRDCalibraMode *calibMode) const;

  Int_t  GetNumberOfGroupsPRF(const char* nametitle) const;
  TH2I  *GetSumCH() const                                           { return fCHSum; };
  TH2D  *GetSumDet() const                                          { return fDetSum;};
  AliTRDCalibraVector  *GetSumDetVector() const                     { return fDetSumVector;};
  TObjArray *GetArrayCalib() const                                  { return fArrayCalib;  };


  void SetHisto2d(Bool_t histo2d)                                   {fHisto2d=histo2d;};
  void SetVector2d(Bool_t vector2d)                                 {fVector2d=vector2d;};
  void SetVdriftLinear(Bool_t vdriftLinear)                         {fVdriftLinear = vdriftLinear;};
  void SetNbTimeBins(Int_t nbTimeBins)                              {fNbTimeBins=nbTimeBins;};
  void SetLow(Int_t low)                                            {flow=low;};
  void SetHigh(Int_t high)                                          {fhigh=high;};
  void SetNz(Short_t nz, Int_t i)                                   {fNz[i]=nz;};
  void SetNrphi(Short_t nrphi, Int_t i)                             {fNrphi[i]=nrphi;};
  void SetFillZero(Bool_t fillZero)                                 {ffillZero =  fillZero;};
  void SetNormalizeNbOfCluster(Bool_t normalizeNbOfCluster)         {fnormalizeNbOfCluster = normalizeNbOfCluster;};
  void SetMaxCluster(Float_t maxcluster)                            {fmaxCluster =  maxcluster; }; 
  void SetOfflineTracks()                                           {fOfflineTracks=kTRUE; fStandaloneTracks=kFALSE; };
  void SetStandaloneTracks()                                        {fStandaloneTracks=kTRUE; fOfflineTracks=kFALSE; };
  void SetCompressPerDetector(Bool_t compressPerDetector=kTRUE)     {fCompressPerDetector=compressPerDetector; };
   
private:
  AliTRDtrackInfo  *fTrackInfo;                  //track info

  AliTRDtrackV1 *ftrdTrack;                      //trdtrack
  AliTRDcluster *fcl;                            //cluster
  
  AliTRDCalibraFillHisto *fTRDCalibraFillHisto;  //! calibration analyse object
  TH1F        *fNbTRDTrack;                      //! nb ESD tracks used
  TH1F        *fNbTRDTrackOffline;               //! nb ESD tracks offline used
  TH1F        *fNbTRDTrackStandalone;            //! nb ESD tracks standalone used
  TH1F        *fNbTRDTracklet;                   //! nb tracklets used
  TH1F        *fNbTRDTrackletOffline;            //! nb tracklets offline used
  TH1F        *fNbTRDTrackletStandalone;         //! nb tracklets standalone used
  TH1F        *fNbTimeBin;                       //! nb Time Bin
  TH1F        *fNbTimeBinOffline;                //! nb Time Bin offline
  TH1F        *fNbTimeBinStandalone;             //! nb Time Bin standalone
  TH1F        *fNbClusters;                      //! nb Clusters
  TH1F        *fNbClustersOffline;               //! nb Clusters offline
  TH1F        *fNbClustersStandalone;            //! nb Clusters standalone
  TProfile2D  *fPHSM;                            //! Per SM
  TH2I        *fCHSM;                            //! Per SM

  TProfile2D  *fPHSum;                           //! All together PH
  TH2I        *fCHSum;                           //! All together CH during the task, used also if want to use AddStatsPerDetector
  TH2D        *fDetSum;                          //! AddStatsPerDetector
  AliTRDCalibraVector *fDetSumVector;            //! AddStatsPerDetector
  
  Bool_t      fHisto2d;                          // histo
  Bool_t      fVector2d;                         // vector
  Bool_t      fVdriftLinear;                     // vdrift Linear
  
  Int_t       flow;                              // lower limit nb of clusters
  Int_t       fhigh;                             // higher limit nb of clusters
  Int_t       fNbTimeBins;                       // number of timebins 
  Bool_t      ffillZero;                         // fill zero
  Short_t     fNz[3];                            // Nz mode 
  Short_t     fNrphi[3];                         // Nrphi mode
  Bool_t      fnormalizeNbOfCluster;             // normalize with number of clusters
  Float_t     fmaxCluster;                       // maxcluster (noise at the end)
  Bool_t      fOfflineTracks;                    // Offline refited tracks
  Bool_t      fStandaloneTracks;                 // Take only standalone tracks

  Bool_t      fCompressPerDetector;              //! Compress per detector 

  TObjArray  *fGraph;                            //! array of graphs filled in PostProcess
  TObjArray  *fArrayCalib;                       //! array of calibration objects filled in PostProcess
  Bool_t      fPostProcess;                      //Post process 

  AliTRDcalibration(const AliTRDcalibration&); 
  AliTRDcalibration& operator=(const AliTRDcalibration&); 

  ClassDef(AliTRDcalibration, 1) // calibration task
};
#endif
