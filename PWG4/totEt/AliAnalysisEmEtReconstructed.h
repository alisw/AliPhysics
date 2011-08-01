#ifndef ALIANALYSISEMETRECONSTRUCTED_H
#define ALIANALYSISEMETRECONSTRUCTED_H
//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for MC analysis
//  - MC output
//
//*-- Author: Marcelo G. Munhoz (USP)
//_________________________________________________________________________

#include "AliAnalysisEtReconstructed.h"
class TParticle;
class TParticlePDG;
class AliESDEvent;
class AliESDtrack;
class AliEMCALTrack;
class TVector3;
class AliEMCALGeometry;
class AliExternalTrackParam;
class AliStack;

class AliAnalysisEmEtReconstructed : public AliAnalysisEtReconstructed
{

public:
   
  AliAnalysisEmEtReconstructed();
  virtual ~AliAnalysisEmEtReconstructed();

    virtual Int_t AnalyseEvent(AliVEvent* event);

    virtual void Init();
    virtual void ResetEventValues();
    virtual void CreateHistograms();
    virtual void FillOutputList(TList* list);

protected:
	
	AliESDtrack* FindMatch(const AliESDCaloCluster *caloCluster, Double_t& Res);
	Double_t GetTrackPID(const AliESDtrack *track) const;
	
	virtual Bool_t GetTrackProjection(AliExternalTrackParam *trackParam, TVector3 &trackPos); // project to a radius
	virtual Bool_t GetTrackProjection(AliEMCALTrack* emcTrack, TVector3 &trackPos, TVector3 clusPos); // project to a point

protected:

	Double_t fResCut;//Marcelo please add comment
	
	Double_t fAllRectotETDep;//Marcelo please add comment
	Double_t fElectronMatchtotETDep;//Marcelo please add comment
	Double_t fNeutralRectotET;//Marcelo please add comment
	Double_t fTotEMRectotET;//Marcelo please add comment
	Double_t fMuonMatchtotETDep, fPionMatchtotETDep, fKaonMatchtotETDep, fProtonMatchtotETDep;//Marcelo please add comment
	Double_t fTotChargedMatchtotETDep;//Marcelo please add comment
	Double_t fTotalRectotETDep;//Marcelo please add comment
	
	AliESDEvent *fESD;//Marcelo please add comment
	AliEMCALGeometry *fGeoUt;//Marcelo please add comment

	// *******************
	// all ET
	// *******************
	//TH2F *fHistAllRecEtaEDepETDep;//Marcelo please add comment 
	//TH2F *fHistAllRecEtaETDep;//Marcelo please add comment 
	
	THnSparseD* fHistAllRecETDep;//Marcelo please add comment
	THnSparseD* fHistAllRec;//Marcelo please add comment
	TH1F *fHistAllRectotETDep;//Marcelo please add comment
	
	// *******************
	// electron ET reconstructed in EMCal
	// *******************
	/*
	TH2F *fHistElectronMatchEtaEDepETDep;//Marcelo please add comment 
	TH2F *fHistElectronMatchEtaPtETDep;//Marcelo please add comment 
	TH2F *fHistElectronMatchEtaETDep;//Marcelo please add comment 
	TH2F *fHistElectronMatchEtaPt;//Marcelo please add comment 
	
	TH2F *fHistElectronRec_ResEDep_ETDep;//Marcelo please add comment 
	TH2F *fHistElectronRec_ResPt_ETDep;//Marcelo please add comment 
	TH2F *fHistElectronRec_ResEDep;//Marcelo please add comment 
	TH2F *fHistElectronRec_ResPt;//Marcelo please add comment 
	*/
	
	THnSparseD* fHistElectronRecETDep;//Marcelo please add comment
	THnSparseD* fHistElectronRec;//Marcelo please add comment
	TH1F *fHistElectronMatchtotETDep;//Marcelo please add comment 
	
	TH2F *fHistElectronRecdEdxP;//Marcelo please add comment

	// *******************
	// Neutral ET reconstructed in EMCal
	// *******************
	/*
	TH2F *fHistNeutralRec_EtaE_ET;//Marcelo please add comment  
	TH2F *fHistNeutralRec_EtaPt_ET;//Marcelo please add comment  
	TH2F *fHistNeutralRec_EtaET;//Marcelo please add comment  
	TH2F *fHistNeutralRec_EtaE;//Marcelo please add comment  
	TH2F *fHistNeutralRec_EtaPt;//Marcelo please add comment  
	*/
	
	TH1F *fHistNeutralRectotET;//Marcelo please add comment  

	// *******************
	// total EM ET reconstructed in EMCal
	// *******************
	TH1F *fHistTotEMRectotET;//Marcelo please add comment

	// *******************
	// muon ET (+ and -)
	// *******************
	/*
	TH2F *fHistMuonMatchEtaEDepETDep;//Marcelo please add comment 
	TH2F *fHistMuonMatchEtaPtETDep;//Marcelo please add comment 
	TH2F *fHistMuonMatchEtaETDep;//Marcelo please add comment 
	TH2F *fHistMuonMatchEtaPt;//Marcelo please add comment 
	
	TH2F *fHistMuonRecResEDepETDep;//Marcelo please add comment 
	TH2F *fHistMuonRecResPtETDep;//Marcelo please add comment 
	TH2F *fHistMuonRecResEDep;//Marcelo please add comment 
	TH2F *fHistMuonRecResPt;//Marcelo please add comment 
	*/
	
	THnSparseD* fHistMuonRecETDep;//Marcelo please add comment
	THnSparseD* fHistMuonRec;//Marcelo please add comment
	TH1F *fHistMuonMatchtotETDep;//Marcelo please add comment 

	TH2F *fHistMuonRecdEdxP;//Marcelo please add comment
	
	// *******************
	// pion ET (+ and -)
	// *******************
	/*
	TH2F *fHistPionMatchEtaEDepETDep;//Marcelo please add comment 
	TH2F *fHistPionMatchEtaPtETDep;//Marcelo please add comment 
	TH2F *fHistPionMatchEtaETDep;//Marcelo please add comment 
	TH2F *fHistPionMatchEtaPt;//Marcelo please add comment 
	
	TH2F *fHistPionRecResEDepETDep;//Marcelo please add comment 
	TH2F *fHistPionRecResPtETDep;//Marcelo please add comment 
	TH2F *fHistPionRecResEDep;//Marcelo please add comment 
	TH2F *fHistPionRecResPt;//Marcelo please add comment 
	*/
	
	THnSparseD* fHistPionRecETDep;//Marcelo please add comment
	THnSparseD* fHistPionRec;//Marcelo please add comment
	TH1F *fHistPionMatchtotETDep;//Marcelo please add comment 

	TH2F *fHistPionRecdEdxP;//Marcelo please add comment

	// *******************
	// charged kaon (+ and -) ET
	// *******************
	/*
	TH2F *fHistKaonMatchEtaEDepETDep;//Marcelo please add comment 
	TH2F *fHistKaonMatchEtaPtETDep;//Marcelo please add comment 
	TH2F *fHistKaonMatchEtaETDep;//Marcelo please add comment 
	TH2F *fHistKaonMatchEtaPt;//Marcelo please add comment 
	
	TH2F *fHistKaonRecResEDepETDep;//Marcelo please add comment 
	TH2F *fHistKaonRecResPtETDep;//Marcelo please add comment 	
	TH2F *fHistKaonRecResEDep;//Marcelo please add comment 
	TH2F *fHistKaonRecResPt;//Marcelo please add comment 
	*/
	
	THnSparseD* fHistKaonRecETDep;//Marcelo please add comment
	THnSparseD* fHistKaonRec;//Marcelo please add comment
	TH1F *fHistKaonMatchtotETDep;//Marcelo please add comment 

	TH2F *fHistKaonRecdEdxP;//Marcelo please add comment
	
	// *******************
	// proton (anti) ET
	// *******************
	/*
	TH2F *fHistProtonMatchEtaEDepETDep;//Marcelo please add comment 
	TH2F *fHistProtonMatchEtaPtETDep;//Marcelo please add comment 
	TH2F *fHistProtonMatchEtaETDep;//Marcelo please add comment 
	TH2F *fHistProtonMatchEtaPt;//Marcelo please add comment 
	
	TH2F *fHistProtonRecResEDepETDep;//Marcelo please add comment 
	TH2F *fHistProtonRecResPtETDep;//Marcelo please add comment 
	TH2F *fHistProtonRecResEDep;//Marcelo please add comment 
	TH2F *fHistProtonRecResPt;//Marcelo please add comment 
	*/
	
	THnSparseD* fHistProtonRecETDep;//Marcelo please add comment
	THnSparseD* fHistProtonRec;//Marcelo please add comment
	TH1F *fHistProtonMatchtotETDep;//Marcelo please add comment 

	TH2F *fHistProtonRecdEdxP;//Marcelo please add comment
	
	// *******************
	// total charged ET
	// *******************
	TH1F *fHistTotChargedMatchtotETDep;//Marcelo please add comment
	
	// *******************
	// total ET
	// *******************
	TH1F *fHistTotalRectotETDep;//Marcelo please add comment
	
	//few checks
	TH2F *fHistDeltaRZ;//Marcelo please add comment
	
 private:
  //Declare it private to avoid compilation warning
    AliAnalysisEmEtReconstructed & operator = (const AliAnalysisEmEtReconstructed & g) ;//cpy assignment
    AliAnalysisEmEtReconstructed(const AliAnalysisEmEtReconstructed & g) ; // cpy ctor
    ClassDef(AliAnalysisEmEtReconstructed, 1);
};

#endif // ALIANALYSISEMETRECONSTRUCTED
