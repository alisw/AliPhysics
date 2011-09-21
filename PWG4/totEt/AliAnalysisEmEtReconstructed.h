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

	Double_t fResCut;//track matching residual cut
	
	Double_t fAllRectotETDep;//total ET deposited - all particles
	Double_t fElectronMatchtotETDep;// total Et deposited - matched electrons
	Double_t fNeutralRectotET;// total ET - neutral particles
	Double_t fTotEMRectotET;// total electromagnetic ET
	Double_t fMuonMatchtotETDep, fPionMatchtotETDep, fKaonMatchtotETDep, fProtonMatchtotETDep;// total ET deposited - hadrons
	Double_t fTotChargedMatchtotETDep;// total Et deposited - matched chagred particles
	Double_t fTotalRectotETDep;// total ET deposited
	
	AliESDEvent *fESD;// ESD object
	AliEMCALGeometry *fGeoUt;// EMCal geometry object

	// *******************
	// all ET
	// *******************
	THnSparseD* fHistAllRecETDep;// ET deposited - all particles
	THnSparseD* fHistAllRec;// mutliplicity - all particles
	TH1F *fHistAllRectotETDep;// total ET deposited - all particles
	
	// *******************
	// electron ET reconstructed in EMCal
	// *******************
	THnSparseD* fHistElectronRecETDep;// Et deposited - matched electrons
	THnSparseD* fHistElectronRec;// multiplicity - matched electrons
	TH1F *fHistElectronMatchtotETDep;// total Et deposited - matched electrons
	
	TH2F *fHistElectronRecdEdxP;// electron dEdx vs p

	// *******************
	// Neutral ET reconstructed in EMCal
	// *******************
	TH1F *fHistNeutralRectotET;// total ET - neutral particles 

	// *******************
	// total EM ET reconstructed in EMCal
	// *******************
	TH1F *fHistTotEMRectotET;// total electromagnetic ET

	// *******************
	// muon ET (+ and -)
	// *******************
	THnSparseD* fHistMuonRecETDep;// Et deposited
	THnSparseD* fHistMuonRec;// multiplicity
	TH1F *fHistMuonMatchtotETDep;// total Et deposited

	TH2F *fHistMuonRecdEdxP;// dEdx vs p
	
	// *******************
	// pion ET (+ and -)
	// *******************
	THnSparseD* fHistPionRecETDep;// Et deposited 
	THnSparseD* fHistPionRec;// multiplicity
	TH1F *fHistPionMatchtotETDep;// total Et deposited

	TH2F *fHistPionRecdEdxP;// dEdx vs p

	// *******************
	// charged kaon (+ and -) ET
	// *******************
	THnSparseD* fHistKaonRecETDep;// Et deposited
	THnSparseD* fHistKaonRec;// multiplicity
	TH1F *fHistKaonMatchtotETDep;// total Et deposited

	TH2F *fHistKaonRecdEdxP;// dEdx vs p
	
	// *******************
	// proton (anti) ET
	// *******************
	THnSparseD* fHistProtonRecETDep;// Et deposited
	THnSparseD* fHistProtonRec;// multiplicity
	TH1F *fHistProtonMatchtotETDep;// total Et deposited

	TH2F *fHistProtonRecdEdxP;// dEdx vs p
	
	// *******************
	// total charged ET
	// *******************
	TH1F *fHistTotChargedMatchtotETDep;// total Et deposited - all charged particles 
	
	// *******************
	// total ET
	// *******************
	TH1F *fHistTotalRectotETDep;// total Et deposited - all particles
	
	//few checks
	TH2F *fHistDeltaRZ;// track-cluster matching residual
	
 private:
  //Declare it private to avoid compilation warning
    AliAnalysisEmEtReconstructed & operator = (const AliAnalysisEmEtReconstructed & g) ;//cpy assignment
    AliAnalysisEmEtReconstructed(const AliAnalysisEmEtReconstructed & g) ; // cpy ctor
    ClassDef(AliAnalysisEmEtReconstructed, 1);
};

#endif // ALIANALYSISEMETRECONSTRUCTED
