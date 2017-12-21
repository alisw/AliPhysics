#ifndef AliPHOSEmbeddingRun2_h
#define AliPHOSEmbeddingRun2_h

// Class to perform embedding on the AOD level
// Author: D.Peressounko

class TChain ;
class TClonesArray ;
class TH2F ;

class AliPHOSClusterizerv1 ;
class AliPHOSReconstructor ;
class AliAODEvent ;
class AliESDEvent ;
class AliESDtrack ;
class AliESDCaloCells ;

#include "AliAnalysisTaskESDfilter.h"

class AliPHOSEmbeddingRun2 : public AliAnalysisTaskESDfilter {
public:
  AliPHOSEmbeddingRun2(const char *name = "AliPHOSEmbeddingRun2");
  virtual ~AliPHOSEmbeddingRun2() {}
  
  //Standard methods
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *){}
  
  //Chain with signal AOD for embedding
  void SetSignalChains(TChain * signal1,TChain * signal2,TChain * signal3) {fAODChain1 =signal1; fAODChain2 =signal2; fAODChain3 =signal3;}
  void SetEmbeddedFlag(Bool_t embpi0, Bool_t embeta,Bool_t embgamma) {fIsPi0Embedded = embpi0; fIsEtaEmbedded =embeta; fIsGammaEmbedded =embgamma;}
  void SetPrivateOADBPath(TString path) {fPathPrivateOADB=path;}//simply copy from AliPhysics which you want to use for analysis is enough.
  void SetSignalCalibration(Double_t corr) {fSignalECorrection=corr;}

private:
  AliPHOSEmbeddingRun2(const AliPHOSEmbeddingRun2&); // not implemented
  AliPHOSEmbeddingRun2& operator=(const AliPHOSEmbeddingRun2&); // not implemented

  void Init() ;
  void InitMF() ; //Mag.Field initialization for track matching
  void InitGeometry() ;
  
  void GetNextSignalEvent(void) ;
  void GetNextSignalEventPi0(void) ;
  void GetNextSignalEventEta(void) ;
  void GetNextSignalEventGamma(void) ;

  void CopyRecalibrateDigits(void) ;
  void MakeEmbedding(AliESDEvent * event, AliAODEvent * signal) ;
  void MakeDigits(AliAODEvent* signal) ;  
  
  Double_t DecalibrateSignal(Double_t cellAmplitude,Int_t cellNumber) ;
 
  //Add new branch
  void ConvertEmbeddedClusters(const AliESDEvent *esd,Int_t what) ;
  void ConvertEmbeddedCells(const AliESDEvent *esd,Int_t what) ;
  void ConvertMCParticles(const AliAODEvent *aod,Int_t what) ;

  Float_t TestCPV(Double_t dx, Double_t dz, Double_t pt, Int_t charge) ;
  Float_t TestCPVRun2(Double_t dx, Double_t dz, Double_t pt, Int_t charge) ;
  

  TChain * fAODChain1 ; //Signal1 (pi0)
  TChain * fAODChain2 ; //Signal (eta)
  TChain * fAODChain3 ; //Signal (gamma)
  
  AliAODEvent * fSignal1 ; //! pi0 signal event  
  AliAODEvent * fSignal2 ; //! eta signal event  
  AliAODEvent * fSignal3 ; //! gamma signal event  

  TTree * fDigitsTree ;  //! Digits
  TTree * fClustersTree; //! Clusters
  TTree * fTreeOut;      //Output AOD
  TClonesArray * fDigitsArr ; //!

  TClonesArray * fEmbeddedClusters1 ; //!
  TClonesArray * fEmbeddedClusters2 ; //!
  TClonesArray * fEmbeddedClusters3 ; //!
  AliAODCaloCells * fEmbeddedCells1 ; //!
  AliAODCaloCells * fEmbeddedCells2 ; //!
  AliAODCaloCells * fEmbeddedCells3 ; //!
  AliESDCaloCells * fCellsPHOS ; //! Old PHOS cells

  AliPHOSClusterizerv1 * fClusterizer ; //!
  AliPHOSReconstructor * fPHOSReconstructor ; //!
  
  TH2F * fOldPHOSCalibration[5] ; //! Calibration coeff. used in ESD production
  AliPHOSCalibData * fSignalCalibData ; //! Decalibration of signal, inverse to OADB. //new memeber variable for an additional calibration.
  TString fPathPrivateOADB; //path to private OADB.
  Double_t           fSignalECorrection;  //! Correction for the Signal clibration
  Int_t fNSignal ; // Number of signal evetns processed  
  Int_t fNSignalPi0 ; // Number of signal evetns processed  
  Int_t fNSignalEta ; // Number of signal evetns processed  
  Int_t fNSignalGamma ; // Number of signal evetns processed  
  Int_t fNCaloClustersOld ; //Number of CaloClusters already in ESD
  Bool_t fInitialized ; //!

	Bool_t fIsPi0Embedded;
	Bool_t fIsEtaEmbedded;
	Bool_t fIsGammaEmbedded;

  ClassDef(AliPHOSEmbeddingRun2, 4); // PHOS analysis task
};

#endif
