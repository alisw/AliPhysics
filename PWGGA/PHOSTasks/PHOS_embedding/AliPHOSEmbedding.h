#ifndef AliPHOSEmbedding_h
#define AliPHOSEmbedding_h

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

#include "AliAnalysisTaskSE.h"

class AliPHOSEmbedding : public AliAnalysisTaskSE {
public:
  AliPHOSEmbedding(const char *name = "AliPHOSEmbedding");
  virtual ~AliPHOSEmbedding() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *){}
  

  void SetSignalChain(TChain * signal){fAODChain =signal;}

  void SetOldCalibration(TH2F **calib) ; 
  //Calibration used in reconstruction of real data (ESDs)
  //If not set, assume same calibration as set by default

private:
  AliPHOSEmbedding(const AliPHOSEmbedding&); // not implemented
  AliPHOSEmbedding& operator=(const AliPHOSEmbedding&); // not implemented

  void Init() ;
  void InitMF() ; //Mag.Field initialization for track matching
  void InitGeometry() ;
  
  AliAODEvent * GetNextSignalEvent(void) ;

  void MakeEmbedding(AliESDEvent * data, AliAODEvent * signal) ;
  void MakeDigits(AliAODEvent* signal) ;  

  void ConvertESDtoAOD(AliESDEvent *esd) ;
  void ConvertHeader(AliESDEvent &esd) ;
  void ConvertPrimaryVertices(const AliESDEvent &esd) ;
  void ConvertCaloClusters(const AliESDEvent &esd) ;
  void ConvertEMCALCells(const AliESDEvent &esd) ;
  void ConvertPHOSCells(const AliESDEvent &esd) ;
  
  void ConvertEmbeddedClusters(const AliESDEvent *esd) ;
  void ConvertEmbeddedCells(const AliESDEvent *esd) ;
  void ConvertMCParticles(const AliAODEvent *aod) ;

  Double_t TPCrp(const AliESDEvent * event) ;
  Bool_t SelectTrack(AliESDtrack * t) ;  
  Float_t TestCPV(Double_t dx, Double_t dz, Double_t pt, Int_t charge) ;
  

  TChain * fAODChain ; //Signal

  TTree * fDigitsTree ;  //! Digits
  TTree * fClustersTree; //! Clusters
  TTree * fTreeOut; //Output AOD
  TClonesArray * fDigitsArr ; //!

  TClonesArray * fEmbeddedClusters ; //!
  AliAODCaloCells * fEmbeddedCells ; //!
  AliESDCaloCells * fCellsPHOS ; //! Old PHOS cells

  AliPHOSClusterizerv1 * fClusterizer ; //!
  AliPHOSReconstructor * fPHOSReconstructor ; //!
  
  TH2F * fOldPHOSCalibration[5] ; //! Calibration coeff. used in ESD production

  Int_t fNSignal ; // Number of signal evetns processed  
  Int_t fNCaloClustersOld ; //Number of CaloClusters already in ESD
  Bool_t fInitialized ; //!
  ClassDef(AliPHOSEmbedding, 1); // PHOS analysis task
};

#endif
