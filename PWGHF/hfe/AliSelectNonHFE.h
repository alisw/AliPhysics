#ifndef ALISELECTNONHFE_H
#define ALISELECTNONHFE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  	Class for the Selection of Non-Heavy-Flavour-Electrons 	      //
//                                                                    //
//  		Author: Elienos Pereira de Oliveira Filho 				  //
//					(epereira@cern.ch)								  //
//					University of São Paulo							  //	
//                                                                    //
////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif

class TH1F;
class TH2F;
class AliVEvent;
class AliVParticle;
class AliESDtrackCuts;
class AliPIDResponse;

class AliSelectNonHFE : public TNamed {
 public:
  AliSelectNonHFE();
  AliSelectNonHFE(const char *name, const Char_t *title);
  virtual ~AliSelectNonHFE();
  Int_t GetNLS() const {return fNLS;};
  Int_t GetNULS() const {return fNULS;};
  Int_t* GetPartnersLS() const {return fLSPartner;};
  Int_t* GetPartnersULS() const {return fULSPartner;};
  Bool_t IsLS() const {return fIsLS;};
  Bool_t IsULS() const {return fIsULS;};
  void FindNonHFE(Int_t iTrack1, AliVParticle *Vtrack1, AliVEvent *fVevent);
  void SetAlgorithm(TString Algorithm) {fAlgorithm = Algorithm;};
	
  void SetAdditionalCuts(Double_t PtMin, Int_t TpcNcls) {fPtMin = PtMin; fTpcNcls = TpcNcls; fHasPtCut=kTRUE; };
  void SetAODanalysis(Bool_t IsAOD) {fIsAOD = IsAOD;};
	
  void SetChi2OverNDFCut(Double_t Chi2OverNDFCut) {fChi2OverNDFCut = Chi2OverNDFCut;};
  void SetDCACut(Double_t dcaCut) {fdcaCut = dcaCut;};
  void SetHistAngleBack(TH1F *HistAngleBack) {fHistAngleBack = HistAngleBack;};
  void SetHistAngle(TH1F *HistAngle) {fHistAngle = HistAngle;};
  void SetHistDCABack(TH1F *HistDCABack) {fHistDCABack = HistDCABack;};
  void SetHistDCA(TH1F *HistDCA) {fHistDCA = HistDCA;};
	
  void SetHistMassBack(TH1F *HistMassBack, Double_t weight=1) {fHistMassBack = HistMassBack; fweight2=weight;};
  void SetHistMass(TH1F *HistMass, Double_t weight=1) {fHistMass = HistMass; fweight1=weight;};
	
  void SetInvariantMassCut(Double_t MassCut) {fMassCut = MassCut;};
  void SetOpeningAngleCut(Double_t AngleCut) {fAngleCut = AngleCut;};
  void SetPIDresponse(AliPIDResponse* PIDResponse) {fPIDResponse = PIDResponse;};
    
  //New cuts to match GSI: Eta, ITS NClust, TPC NCluster for PID,  Global Tracks, DCA
  void SetEtaCutForPart(Double_t EtaMin, Double_t EtaMax) {fUseEtaCutForPart = kTRUE; fEtaCutMin = EtaMin; fEtaCutMax = EtaMax; } ;
    void SetTPCNclsForPID(Int_t NClustPID) {fRequireTPCNclusForPID = kTRUE; fTpcNclsPID = NClustPID; };
  void SetNClustITS(Int_t nclus) {fNClusITS = nclus; fRequirePointOnITS = kTRUE; };
  void SetUseGlobalTracks() {fUseGlobalTracks = kTRUE;};
  void SetDCAPartnerCuts(Double_t xy, Double_t z) {fUseDCAPartnerCut = kTRUE; fDCAcutxyPartner = xy; fDCAcutzPartner = z;};
  void SetUseITSTPCRefit(Bool_t Use) {fRequireITSAndTPCRefit = Use;};
	
	void SetTrackCuts(Double_t TPCnSigmaMin, Double_t TPCnSigmaMax) 
	{
		fTPCnSigmaMin = TPCnSigmaMin; 
		fTPCnSigmaMax = TPCnSigmaMax; 
	};
	
	void SetTrackCuts(Double_t TPCnSigmaMin, Double_t TPCnSigmaMax, AliESDtrackCuts* TrackCuts) 
	{
		fTPCnSigmaMin = TPCnSigmaMin; 
		fTPCnSigmaMax = TPCnSigmaMax; 
		fTrackCuts = TrackCuts;
	};
  
 private:
  AliESDtrackCuts	*fTrackCuts;		//! Track quality
  TString		fAlgorithm;       		//Algorithm choice: "MA" (Manual Algorithm) or "KF" (Kalman Filter Algorithm)
  Double_t 		fAngleCut; 				//Maximum opening angle between the tracks
  Double_t 		fdcaCut;   				//Maximum dca between the tracks
  Double_t 		fTPCnSigmaMin; 			//Minimum partner TPCnSigma value
  Double_t 		fTPCnSigmaMax;  		//Maximum partner TPCnSigma value
  Double_t 		fMassCut;				//Maximum Invariant Mass Value for Non-HF-Electrons
  Double_t		fChi2OverNDFCut;        //Maximum value of Chi2 over NDF in the case of KF Algorithm
  Double_t		fPtMin;					//Minimum pT value for the associated electrons
  Bool_t		fIsLS;					//Unlike signal pairs Flag
  Bool_t		fIsULS;					//like signal pairs Flag
  Bool_t		fIsAOD;					//Flag for AOD analysis
  Bool_t		fHasPtCut;				//Flag for additional pT cut
  Int_t			fNLS;					//Number of Unlike signal pairs
  Int_t			fNULS;					//Number of like signal pairs
  Int_t			fTpcNcls;				//Minimum Number of clusters in the TPC
  Double_t      fweight1;				//weight to fill the invariant mass plot
  Double_t      fweight2;				//weight to fill the invariant mass plot
    
  Bool_t fUseGlobalTracks; //Use Global tracks for partner
    
  Bool_t fRequirePointOnITS; // Require a specific number of clusters for partners on ITS
  Int_t	        fNClusITS; //Number of Clusters required
 
  Bool_t fUseDCAPartnerCut; //flag for cut
  Double_t fDCAcutxyPartner; //DCA cut partner xy
  Double_t fDCAcutzPartner; // DCA cut partner z
  
  Bool_t fRequireTPCNclusForPID; //require a TPCNClus in PID
  Int_t fTpcNclsPID; // Number of TPC NClus
    
  Bool_t fUseEtaCutForPart; // Use eta cut in partners
  Double_t fEtaCutMin; // min eta cut
  Double_t fEtaCutMax; // max eta cut
   
  Bool_t fRequireITSAndTPCRefit;

	
  Int_t			*fLSPartner;	        //! Pointer for the LS partners index
  Int_t			*fULSPartner;	        //! Pointer for the ULS partners index
  TH1F			*fHistMass;				//! Invariant mass histogram for Unlike sign pairs
  TH1F			*fHistMassBack;         //! Invariant mass histogram for like sign pairs
  TH1F			*fHistDCA;				//! DCA histogram for Unlike sign pairs
  TH1F			*fHistDCABack;	        //! DCA histogram for like sign pairs
  TH1F			*fHistAngle;	        //! Opening Angle histogram for Unlike sign pairs
  TH1F			*fHistAngleBack;        //! Opening Angle histogram for like sign pairs
  AliPIDResponse *fPIDResponse;     	//! PID response object
  
  AliSelectNonHFE(const AliSelectNonHFE&); // not implemented
  AliSelectNonHFE& operator=(const AliSelectNonHFE&); // not implemented
  
  ClassDef(AliSelectNonHFE, 2); //!example of analysis
};

#endif
