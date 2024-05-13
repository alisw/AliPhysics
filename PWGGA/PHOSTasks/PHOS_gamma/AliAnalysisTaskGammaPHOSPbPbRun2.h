#ifndef AliAnalysisTaskGammaPHOSPbPbRun2_h
#define AliAnalysisTaskGammaPHOSPbPbRun2_h

class TObjArray;
class THashList ;
class TH1F;
class TH1D ;
class TH2I;
class TH2F;
class TH3F;
class TF1 ;
class AliStack ;
class AliESDtrackCuts;
class AliPHOSGeometry;
class AliAODEvent ;
class AliPHOSCalibData;
class AliAODtrack ;
class AliAODCaloCluster ;
class AliEPFlattener ;
class AliFlowTrackCuts ;
class AliFlowEvent ;
class AliFlowVector ;
class AliAODMCParticle ;
class AliCaloPhoton;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskGammaPHOSPbPbRun2 : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskGammaPHOSPbPbRun2(const char *name = "AliAnalysisTaskGammaPHOSPbPbRun2");
  virtual ~AliAnalysisTaskGammaPHOSPbPbRun2() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  void SetHarmonics(Int_t n=2){fHarmonics=n ;}
  void SetDistCut(Bool_t on=kTRUE){fDistCut=on ;}
  void SetNCenBins(Int_t nbins) {fNCenBins = nbins;}
  void SetCentralityIntervals(TString mode);
  
private:
  AliAnalysisTaskGammaPHOSPbPbRun2(const AliAnalysisTaskGammaPHOSPbPbRun2&); // not implemented
  AliAnalysisTaskGammaPHOSPbPbRun2& operator=(const AliAnalysisTaskGammaPHOSPbPbRun2&); // not implemented

  void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z, Double_t w) const ; //Fill 3D histogram witn name key

  Int_t ConvertRunNumber(Int_t run) ; 

  void EvalV0ReactionPlane(AliAODEvent * event) ;
  void ApplyFinalFlattening() ;//apply final fine flattening
  void ApplyFinalQFlattening() ;//apply final fine flattening
  Bool_t GetTPCEventPlane(Double_t &epAngle, Double_t &qsubRes) ;
  TObjArray* GetEventPlaneTracks(Int_t &maxID) ;
  Double_t GetWeight(TObject* track1) ;
  Double_t GetPhiWeight(TObject* track1) ;

  void EvalResolution() ;
  void EvalQResolution() ;
  Double_t PHOSMultiplicity();
  Double_t CentralityWeight(Double_t c) ;

  void AddQAHistograms();
  void AddSinglePhotonHistograms();
  void AddEventPlaneHistograms();
  void AddMCHistograms();
  void AddDistBadHistograms();
  void AddPhiTitleHistograms();

  void TestMatchingTrackPID(AliCaloPhoton *ph, Bool_t mix);
  
  static  Bool_t PythiaInfoFromFile(TString currFile, Float_t & xsec, Float_t & trials) ;


private:
  THashList * fOutputContainer;    //final histogram container
  THashList * fOutputContainer2;    //final histogram container
  THashList * fOutputContainer3;    //final histogram container
  THashList * fOutputContainer4;    //final histogram container
  THashList * fOutputContainer5;    //final histogram container
  THashList * fOutputContainer6;    //final histogram container

  AliAODEvent      *fEvent; //! input event 
				    //
  TList * fPHOSEvents[1][10][11] ; // Containers for events with PHOS photons
  TClonesArray * fPHOSEvent ;      // PHOS photons in current event
  AliPIDResponse *     fPIDResponse;    // Pid response
 
  //Vector value
  Double_t fQV0A ;   //Lengths of flow vectors
  Double_t fQV0C ;
  Double_t fQTPC ;
 
  //Reaction plain for v2
  Float_t fRP ;        //! Reaction plane calculated with full TPC 
  Float_t fRPV0A ;     //! Reaction plain calculated with A-side TPC: eta>0.15 
  Float_t fRPV0C ;     //! Reaction plain calculated with C-side TPC: eta<-0.15
  Float_t fRPQ ;       //! Reaction plane calculated with full TPC 
  Float_t fRPQV0A ;    //! Reaction plain calculated with A-side TPC: eta>0.15 
  Float_t fRPQV0C ;    //! Reaction plain calculated with C-side TPC: eta<-0.15
  Double_t fV0Ares ;   //! resolution of EP
  Double_t fV0Cres ;   //! resolution of EP
  Double_t fTPCres ;   //! resolution of EP
  Double_t fV0AQres ;  //! resolution of EP
  Double_t fV0CQres ;  //! resolution of EP
  Double_t fTPCQres ;  //! resolution of EP
  Bool_t fHaveTPCRP ;  //! Is TPC RP defined?
  TH1F * fPhiDist ;    //!
  AliEPFlattener * fV0AFlat ; //!
  AliEPFlattener * fV0CFlat ; //!
  AliEPFlattener * fTPCFlat ; //!
  AliEPFlattener * fV0AQFlat ; //!
  AliEPFlattener * fV0CQFlat ; //!
  AliEPFlattener * fTPCQFlat ; //!
  
  Int_t  fHarmonics ; //Harminic in use
  Bool_t fDistCut ;

  Int_t   fRunNumber ;    //Current run number
  Float_t fCentrality ; //!Centrality of the currecnt event

  Int_t fCenBin ;       //! Current centrality bin

  AliPHOSGeometry  *fPHOSGeo;  //! PHOS geometry
  Int_t fEventCounter;         // number of analyzed events

  Int_t fInPHOS;

  TH1D*  fTPCfinalC2 ;   //! HIstos with flattening parameters
  TH1D*  fTPCfinalS2 ;   //!
  TH1D*  fTPCfinalC4 ;   //!
  TH1D*  fTPCfinalS4 ;   //!

  TH1D*  fV0AfinalC2 ;   //!
  TH1D*  fV0AfinalS2 ;   //!
  TH1D*  fV0AfinalC4 ;   //!
  TH1D*  fV0AfinalS4 ;   //!
    
  TH1D*  fV0CfinalC2 ;   //!
  TH1D*  fV0CfinalS2 ;   //!
  TH1D*  fV0CfinalC4 ;   //!
  TH1D*  fV0CfinalS4 ;   //!

  TH1D*  fTPCfinalQC2 ;   //! HIstos with flattening parameters
  TH1D*  fTPCfinalQS2 ;   //!
  TH1D*  fTPCfinalQC4 ;   //!
  TH1D*  fTPCfinalQS4 ;   //!

  TH1D*  fV0AfinalQC2 ;   //!
  TH1D*  fV0AfinalQS2 ;   //!
  TH1D*  fV0AfinalQC4 ;   //!
  TH1D*  fV0AfinalQS4 ;   //!
    
  TH1D*  fV0CfinalQC2 ;   //!
  TH1D*  fV0CfinalQS2 ;   //!
  TH1D*  fV0CfinalQC4 ;   //!
  TH1D*  fV0CfinalQS4 ;    //!
			   
  AliFlowTrackCuts * fCutsV0 ; //! 
  AliFlowTrackCuts * fCutsTPC ; //!
  AliFlowEvent     * fFlowEvent ; //!

  TClonesArray *       fMCArray;  // MC array
  Bool_t fIsMC; //MC flag

  std::vector<std::pair<TString, TString>> fPidCuts; // names and titles for pid cuts
  std::vector<Double_t> fCenBinEdges;

  Double_t fTOF; //time of flight
  Int_t fNCenBins;

  TString              fCurrFileName;      // current file path name
  Bool_t               fCheckMCCrossSection; // retrieve from the pyxsec.root file only if requested

  Bool_t Notify();

  TH1F *               fh1Xsec ;          //! Xsec pythia
  TH1F *               fh1Trials ;        //! trials pythia
  Float_t              fAvgTrials;         // avg trials
					
  Int_t GetPrimaryLabel(AliVCluster *clu);
  Int_t GetPrimaryLabelAtVertex(AliVCluster *clu);
  Int_t FindTrackMatching(Int_t mod,TVector3 *locpos, Double_t &dx, Double_t &dz, Double_t &pt,Int_t &charge);
  Double_t TestCPV(Double_t dx, Double_t dz, Double_t pt, Int_t charge);

  Double_t fVtx0[3] = {0., 0., 0}, 
	   fVtx5[3] = {0., 0., 0.};

  void ProcessMC();

  Double_t getR(AliAODMCParticle *particle);
  Int_t Pdg2Index(Int_t pdg);


  
  ClassDef(AliAnalysisTaskGammaPHOSPbPbRun2, 3); // PHOS analysis task
};

#endif
