#ifndef ALIANALYSISTASKJETIPQA_H
#define ALIANALYSISTASKJETIPQA_H
#include "AliHFJetsTagging.h"
#include "AliAnalysisTaskEmcalJet.h"
class AliEmcalJet;
class AliRDHFJetsCuts;
class AliAODVertex;
class AliAODTrack;
class TList;
class TH1D;
class TH2D;
class AliHFJetsTagging;
class TParticle;
class TClonesArray;
class AliAODMCParticle;
class AliMCEvent;
class AliESDEvent;
class AliESDtrack;
class AliAnalysisUtils;
#include "AliESDtrackCuts.h"



class AliAnalysisTaskHFJetIPQA: public AliAnalysisTaskEmcalJet
{
public:
	enum heavyType {kCharm=4, kBeauty=5, kOthers=6, kElectronPDG=11};
	enum SourceType {
		kDirectCharm=1,
		kDirectBeauty=2,
		kBeautyCharm=3,
		kGamma=4,
		kPi0=5,
		kElse=6,
		kMisID=7,
		kEta=8,
		kOmega=9,
		kPhi=10,
		kEtaPrime=11,
		kRho0=12,
		kGammaPi0=13,
		kGammaEta=14,
		kGammaOmega=15,
		kGammaPhi=16,
		kGammaEtaPrime=17,
		kGammaRho0=18,
		kJpsi=19,
		kB2Jpsi=20,
		kKe3=21,
		kGammaB2M=22,
		kGammaD2M=23,
		kGammaM2M=24,
		kB2M=25,
		kD2M=26,
		kM2M=27,
		kK0L=28,
		kK0s2P=31,
		kK0l2P=32,
		kLamda2P=33,
		kSigma2P=34,
		kGammaK0s2P=35,
		kGammaK0l2P=36,
		kGammaLamda2P=37,
		kGammaSigma2P=38,
		kK2P=39,
		kGammaK2P=40
	};


	enum particletype
	{
	bPi0=111,
	bEta=221,
	bEtaPrime=331,
	bPhi=333,
	bRho=113,
	bOmega=223,
	bSigma0=3212,
	bK0s=310,
	bLambda=3122,
	bPi=211,
	bProton=2212,
	bKaon=321,
	bOmegaBaryon=3334,
	bAntiOmegaBaryon=-3334,
	bXiBaryon=3312,
	bAntiXiBaryon=-3312,
	bD0=411,
	bDPlus=421,
	bDStarPlus=413,
	bDSPlus=431,
	bK0l=431,
	bSigmaPlus = 3222,
	bRhoPlus=213,
	bBPlus = 521,
	bB0 = 511,
	bLambdaB =5122,
	bLambdaC=4122,
	bBStarPlus=523
	};

	enum particlearraxidx
	{
	bIdxPi0=0,
	bIdxEta=1,
	bIdxEtaPrime=2,
	bIdxPhi=3,
	bIdxRho=4,
	bIdxOmega=5,
	bIdxK0s=6,
	bIdxLambda=7,
	bIdxPi=8,
	bIdxProton=9,
	bIdxKaon=10,
	bIdxD0=11,
	bIdxDPlus=12,
	bIdxDStarPlus=13,
	bIdxDSPlus=14,
	bIdxLambdaC=15,
	bIdxBPlus = 16,
	bIdxB0 = 17,
	bIdxLambdaB = 18,
	bIdxBStarPlus=19
	};


	struct myvaluetuple {
		myvaluetuple(double v1,double v2, bool b,bool c){
			first = v1;
			second = v2;
			is_electron = b;
			is_fromB = c;

		};
		double first; // to be compatible with std::pair
		double second;// to be compatible with std::pair
		bool   is_electron; // added for electron contribution check
		bool   is_fromB; // added for electron contribution check
	};



	Int_t fNparents; // number of heavy hadrons to be considered
	Int_t fParentSelect[2][7]; // heavy hadron species
	static const Int_t fgkMaxGener=10; // ancester level wanted to be checked
	static const Int_t fgkMaxIter=100; // number of iteration to find out matching particle


	enum TTypeImpPar {kXY,kXYSig,kXYZ,kXYZSig,kXYZSigmaOnly};
	AliAnalysisTaskHFJetIPQA();
	AliAnalysisTaskHFJetIPQA(const char *name);
	virtual ~AliAnalysisTaskHFJetIPQA(){;}
	virtual void  UserCreateOutputObjects();
	virtual Bool_t Run();

	void SetESDCuts (AliESDtrackCuts  *cuts =NULL){fESDTrackCut =  new AliESDtrackCuts(*cuts);};
	void SetRunESD (Bool_t val = kTRUE){fESD = val;};


	/*
	AliAnalysisTaskHFJetIPQA(const AliAnalysisTaskHFJetIPQA&);
	AliAnalysisTaskHFJetIPQA& operator=(const AliAnalysisTaskHFJetIPQA&);*/
	virtual AliRDHFJetsCuts* GetJetCutsHF(){return fJetCutsHF;};
	Bool_t IsEventSelected();
	  enum EPileup {kNoPileupSelection,kRejectPileupEvent,kRejectTracksFromPileupVertex};
	  enum ERejBits {kNotSelTrigger,kNoVertex,kTooFewVtxContrib,kZVtxOutFid,kPileupSPD,kOutsideCentrality,kVertexZContrib,kPhysicsSelection,kNoContributors,kDeltaVertexZ,kNoVertexTracks,kVertexZResolution,kMVPileup,kSPDClusterCut,kZVtxSPDOutFid,kCentralityFlattening};
	Bool_t IsSelected(AliVEvent *event, Int_t &WhyRejected,ULong_t &RejectionBits);

	void SetUseMonteCarloWeighingLinus(
			TH1F *Pi0 ,
			TH1F *Eta,
			TH1F *EtaP,
			TH1F *Rho,
			TH1F *Phi,
			TH1F *Omega,
			TH1F *K0s,
			TH1F *Lambda,
			TH1F *ChargedPi,
			TH1F *ChargedKaon,
			TH1F *Proton,
			TH1F *D0,
			TH1F *DPlus,
			TH1F *DStarPlus,
			TH1F *DSPlus,
			TH1F *LambdaC,
			TH1F *BPlus,
			TH1F *B0,
			TH1F *LambdaB,
			TH1F *BStarPlus
 //
	)
	{
		for(int i =1 ; i< Pi0->GetNbinsX()+1;++i){
			fBackgroundFactorLinus[bIdxPi0][i-1] =Pi0->GetBinContent(i);
			fBackgroundFactorLinus[bIdxEta][i-1] =Eta->GetBinContent(i);
			fBackgroundFactorLinus[bIdxEtaPrime][i-1] =EtaP->GetBinContent(i);
			fBackgroundFactorLinus[bIdxRho][i-1] =Rho->GetBinContent(i);
			fBackgroundFactorLinus[bIdxPhi][i-1] =Phi->GetBinContent(i);
			fBackgroundFactorLinus[bIdxOmega][i-1] =Omega->GetBinContent(i);
			fBackgroundFactorLinus[bIdxK0s][i-1] =K0s->GetBinContent(i);
			fBackgroundFactorLinus[bIdxLambda][i-1] =Lambda->GetBinContent(i);
			fBackgroundFactorLinus[bIdxPi][i-1] =ChargedPi->GetBinContent(i);
			fBackgroundFactorLinus[bIdxKaon][i-1] =ChargedKaon->GetBinContent(i);
			fBackgroundFactorLinus[bIdxProton][i-1] =Proton->GetBinContent(i);
			fBackgroundFactorLinus[bIdxD0][i-1] =D0->GetBinContent(i);
			fBackgroundFactorLinus[bIdxDPlus][i-1] =DPlus->GetBinContent(i);
			fBackgroundFactorLinus[bIdxDStarPlus][i-1] =DStarPlus->GetBinContent(i);
			fBackgroundFactorLinus[bIdxDSPlus][i-1] =DSPlus->GetBinContent(i);
			fBackgroundFactorLinus[bIdxLambdaC][i-1] =LambdaC->GetBinContent(i);
			fBackgroundFactorLinus[bIdxBPlus][i-1] =BPlus->GetBinContent(i);
			fBackgroundFactorLinus[bIdxB0][i-1] =B0->GetBinContent(i);
			fBackgroundFactorLinus[bIdxLambdaB][i-1] =LambdaB->GetBinContent(i);
			fBackgroundFactorLinus[bIdxBStarPlus][i-1] =BStarPlus->GetBinContent(i);
 //
		}
		return;
	};




	void SetUseMonteCarloWeighing(
			TH1D *PionHist ,
			TH1D *EtaHist,
			TH1D *OmegaHist,
			TH1D *PhiHist,
			TH1D *EtapHist,
			TH1D *RhoHist,
			TH1D *KaonHist,
			TH1D *K0sHist,
			TH1D *LambdaHist
	)
	{
		const Double_t binLimit[45] =  {
				0.1,0.112797,0.127231,0.143512,0.161877,0.182592,0.205957,0.232313,0.262041,
				0.295573,0.333397,0.37606,0.424183,0.478465,0.539692,0.608754,0.686654,0.774523,0.873636,0.985432,
				1.11153,1.25377,1.41421,1.59519,1.79932,2.02957,2.28928,2.58223,2.91267,3.2854,3.70582,4.18004,
				4.71494,5.3183,5.99886,6.76651,7.6324,8.60909,9.71076,10.9534,12.3551,13.9361,15.7195,17.731,20
		};
		for (int i=0; i<45;++i){
			fBackgroundFactorBins[i] = binLimit[i];
		}
		for (int i=0; i<44;++i){
			if (PionHist)fBackgroundFactor[0][i] = PionHist->GetBinContent(i+1);
			if (EtaHist)fBackgroundFactor[1][i] = EtaHist->GetBinContent(i+1);
			if (OmegaHist)fBackgroundFactor[2][i] = OmegaHist->GetBinContent(i+1);
			if (PhiHist)fBackgroundFactor[3][i] = PhiHist->GetBinContent(i+1);
			if (EtapHist)fBackgroundFactor[4][i] = EtapHist->GetBinContent(i+1);
			if (RhoHist)fBackgroundFactor[5][i] = RhoHist->GetBinContent(i+1);
			if (KaonHist)fBackgroundFactor[6][i] = KaonHist->GetBinContent(i+1);
			if (K0sHist)fBackgroundFactor[7][i] = K0sHist->GetBinContent(i+1);
			if (LambdaHist)fBackgroundFactor[8][i] = LambdaHist->GetBinContent(i+1);
		}

	};

	void UseCorrectedRhoPt(bool val = true){fUseCorrPt =val;};
	void UseGammaV0RejectionESD(bool val = true){fEnableV0GammaRejection =val;};
private:

	Double_t GetArmenteros(AliESDv0 * v0 , int pidneg,int pidpos ,double &alpha);
	Double_t GetPsiPair(AliESDv0 * v0);
	Bool_t CalculateTrackImpactParameter(AliAODTrack * track,double *impar, double * cov); // Removes track from Vertex calculation first
	Bool_t CalculateTrackImpactParameter(AliESDtrack * track,double *impar, double * cov); // Removes track from Vertex calculation first
	Bool_t CalculateTrackImpactParameterTruth(AliAODTrack * track,double *impar, double * cov); // calculates DCA on MC particle/event information
	Bool_t CalculateTrackImpactParameterTruth(AliESDtrack * track,double *impar, double * cov); // calculates DCA on MC particle/event information
	Bool_t CalculateJetSignedTrackImpactParameter(AliAODTrack * track,AliEmcalJet * jet ,double *impar, double * cov, double &sign, double &dcajetrack, double &lineardecaylength);
	Bool_t CalculateJetSignedTrackImpactParameter(AliESDtrack * track,AliEmcalJet * jet ,double *impar, double * cov, double &sign, double &dcajetrack, double &lineardecaylength);
	Double_t GetValImpactParameter(TTypeImpPar type,double *impar, double * cov);
	Bool_t IsV0PhotonFromBeamPipeDaughter(const AliAODTrack* track);
	Bool_t IsV0PhotonFromBeamPipeDaughter(const AliESDtrack* track);
	Bool_t IsTrackAccepted(AliAODTrack* track);
	Bool_t IsTrackAccepted(AliESDtrack* track);
	Bool_t MatchJetsGeometricDefault(); //jet matching function 1/4
	void   DoJetLoop(); //jet matching function 2/4
	void   SetMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, int matching=0);
	void   GetGeometricalMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d) const;
	Double_t GetMonteCarloCorrectionFactor(AliAODTrack* track,bool &ise,bool &fromB);
	Double_t GetMonteCarloCorrectionFactor(AliESDtrack* track,bool &ise,bool &fromB);
	Int_t GetElecSource(const AliAODMCParticle *  const mcpart,Double_t &mpt) const;
	Double_t GetWeightFactor(const AliAODMCParticle * const mcpart, const Int_t iBgLevel);
	Double_t GetWeightFactorLinus( AliAODMCParticle * mcpart,bool &isTrackFromPromptB);
	Double_t GetWeightFactorLinus( AliMCParticle * mcpart,bool &isTrackFromPromptB);
	Bool_t ParticleIsPossibleSource(int pdg);
	Bool_t IsSelectionParticle( AliAODMCParticle * mcpart ,int &pdg,double &pT,int &idx  );
	Bool_t IsSelectionParticle( AliMCParticle * mcpart ,int &pdg,double &pT,int &idx  );
	Bool_t IsSecondaryFromWeakDecay( AliAODMCParticle * particle ) ;
	Bool_t IsSecondaryFromWeakDecay( AliMCParticle * particle ) ;

	bool   IsPromptDMeson(AliAODMCParticle * part );
	bool   IsPromptDMeson(AliMCParticle * part );

	bool   IsPromptBMeson(AliAODMCParticle * part );
	bool   IsPromptBMeson(AliMCParticle * part );

	Bool_t GetBMesonWeight( AliAODMCParticle * mcpart ,int &pdg,double &pT,int &idx  );
	Bool_t GetBMesonWeight( AliMCParticle * mcpart ,int &pdg,double &pT,int &idx  );

	AliAODMCParticle* GetMCTrack( const AliAODTrack* _track);

	Float_t  GetRapidity(const TParticle *part);

	AliRDHFJetsCuts* fJetCutsHF;//
	AliHFJetsTagging* fHFJetUtils;//!

	Bool_t fUseCorrPt;
	Bool_t fESD ;
	Bool_t fEnableV0GammaRejection;

	TH1D * fh1dEventRejectionRDHFCuts;//! Store Rejection reasons and number of accepted events
	TH1D * fh1dVertexZ;//!
	TH1D * fh1dVertexZAccepted;//!
	TH2D * fh1dVertexR;//!
	TH2D * fh1dVertexRAccepted;//!


	TH1D * fh1dTracksAccepeted;//!

	TH1D * fh1dTracksImpParXY;//! R Impact Parameter
	TH1D * fh1dTracksImpParXYZ;//! R+z Impact Parameter
	TH1D * fh1dTracksImpParXYSignificance;//! R Impact Parameter Significance
	TH1D * fh1dTracksImpParXYZSignificance;//! R+z Impact Parameter Significance

	TH1D * fh1dTracksImpParXYTruth;//! R Impact Parameter truth
	TH1D * fh1dTracksImpParXYZTruth;//! R+z Impact Parameter truth
	TH1D * fh1dTracksImpParXYResidualTruth;//! Delta R Impact Parameter over uncertainty
	TH1D * fh1dTracksImpParXYZResidualTruth;//! Delta R+z Impact Parameter truth over uncertainty
	// inclusive  impact parameter distributions  for tracks with light meson (TODO D meson) correction in Monte Carlo
	TH1D * fh1dTracksImpParXY_McCorr;//!
	TH1D * fh1dTracksImpParXYZ_McCorr;//!
	TH1D * fh1dTracksImpParXYSignificance_McCorr;//!
	TH1D * fh1dTracksImpParXYZSignificance_McCorr;//!

	TH2D * fh2dVertexChi2NDFNESDTracks;//!

	TH1D * fh1dJetGenPt; //! Generator level jet momentum for unfolding
	TH1D * fh1dJetGenPtUnidentified; //!
	TH1D * fh1dJetGenPtudsg; //!
	TH1D * fh1dJetGenPtc; //!
	TH1D * fh1dJetGenPtb; //!

	TH1D * fh1dJetRecPt; //! Detector level jets
	TH1D * fh1dJetRecPtAccepted; //! Detector level jets accepted
	TH2D * fh1dJetRecEtaPhiAccepted; //! Detector level jets accepted

	TH1D * fh1dJetRecPtUnidentified; //!
	TH1D * fh1dJetRecPtudsg; //!
	TH1D * fh1dJetRecPtc; //!
	TH1D * fh1dJetRecPtb; //!
	TH1D * fh1dJetRecPtUnidentifiedAccepted; //!
	TH1D * fh1dJetRecPtudsgAccepted; //!
	TH1D * fh1dJetRecPtcAccepted; //!
	TH1D * fh1dJetRecPtbAccepted; //!

	TH2D * fh2dJetGenPtVsJetRecPt; //! raw momentum response matrix
	// inclusive signed impact parameter distributions
	TH2D * fh2dJetSignedImpParXY; //!
	TH2D * fh2dJetSignedImpParXYUnidentified; //!
	TH2D * fh2dJetSignedImpParXYudsg; //!
	TH2D * fh2dJetSignedImpParXYb; //!
	TH2D * fh2dJetSignedImpParXYbNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYc; //!

	TH2D * fh2dJetSignedImpParXYSignificance; //!
	TH2D * fh2dJetSignedImpParXYSignificanceUnidentified; //!
	TH2D * fh2dJetSignedImpParXYSignificanceudsg; //!
	TH2D * fh2dJetSignedImpParXYSignificanceb; //!
	TH2D * fh2dJetSignedImpParXYSignificancebNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYSignificancec; //!

	TH2D * fh2dJetSignedImpParXYZ; //!
	TH2D * fh2dJetSignedImpParXYZUnidentified; //!
	TH2D * fh2dJetSignedImpParXYZudsg; //!
	TH2D * fh2dJetSignedImpParXYZb; //!
	TH2D * fh2dJetSignedImpParXYZbNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYZc; //!

	TH2D * fh2dJetSignedImpParXYZSignificance; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceUnidentified; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceudsg; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceb; //!
	TH2D * fh2dJetSignedImpParXYZSignificancebNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYZSignificancec; //!
	// inclusive signed impact parameter distributions with light meson (TODO D meson) correction in Monte Carlo
	TH2D * fh2dJetSignedImpParXY_McCorr; //!
	TH2D * fh2dJetSignedImpParXYUnidentified_McCorr; //!
	TH2D * fh2dJetSignedImpParXYudsg_McCorr; //!
	TH2D * fh2dJetSignedImpParXYb_McCorr; //!
	TH2D * fh2dJetSignedImpParXYb_McCorrNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYc_McCorr; //!

	TH2D * fh2dJetSignedImpParXYSignificance_McCorr; //!
	TH2D * fh2dJetSignedImpParXYSignificanceUnidentified_McCorr; //!
	TH2D * fh2dJetSignedImpParXYSignificanceudsg_McCorr; //!
	TH2D * fh2dJetSignedImpParXYSignificanceb_McCorr; //!
	TH2D * fh2dJetSignedImpParXYSignificanceb_McCorrNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYSignificancec_McCorr; //!

	TH2D * fh2dJetSignedImpParXYZ_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZUnidentified_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZudsg_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZb_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZb_McCorr_McCorrNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYZc_McCorr; //!

	TH2D * fh2dJetSignedImpParXYZSignificance_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceUnidentified_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceudsg_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceb_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceb_McCorrNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYZSignificancec_McCorr; //!



	// inclusive signed impact parameter distributions
	//First
	TH2D * fh2dJetSignedImpParXYFirst; //!
	TH2D * fh2dJetSignedImpParXYUnidentifiedFirst; //!
	TH2D * fh2dJetSignedImpParXYudsgFirst; //!
	TH2D * fh2dJetSignedImpParXYbFirst; //!
	TH2D * fh2dJetSignedImpParXYbFirstNonBDecay; //
	TH2D * fh2dJetSignedImpParXYcFirst; //!

	TH2D * fh2dJetSignedImpParXYSignificanceFirst; //!
	TH2D * fh2dJetSignedImpParXYSignificanceUnidentifiedFirst; //!
	TH2D * fh2dJetSignedImpParXYSignificanceudsgFirst; //!
	TH2D * fh2dJetSignedImpParXYSignificancebFirst; //!
	TH2D * fh2dJetSignedImpParXYSignificancebFirstNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYSignificancecFirst; //!

	TH2D * fh2dJetSignedImpParXYZFirst; //!
	TH2D * fh2dJetSignedImpParXYZUnidentifiedFirst; //!
	TH2D * fh2dJetSignedImpParXYZudsgFirst; //!
	TH2D * fh2dJetSignedImpParXYZbFirst; //!
	TH2D * fh2dJetSignedImpParXYZbFirstNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYZcFirst; //!

	TH2D * fh2dJetSignedImpParXYZSignificanceFirst; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceudsgFirst; //!
	TH2D * fh2dJetSignedImpParXYZSignificancebFirst; //!
	TH2D * fh2dJetSignedImpParXYZSignificancebFirstNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYZSignificancecFirst; //!
	//Second
	TH2D * fh2dJetSignedImpParXYSecond; //!
	TH2D * fh2dJetSignedImpParXYUnidentifiedSecond; //!
	TH2D * fh2dJetSignedImpParXYudsgSecond; //!
	TH2D * fh2dJetSignedImpParXYbSecond; //!
	TH2D * fh2dJetSignedImpParXYbSecondNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYcSecond; //!

	TH2D * fh2dJetSignedImpParXYSignificanceSecond; //!
	TH2D * fh2dJetSignedImpParXYSignificanceUnidentifiedSecond; //!
	TH2D * fh2dJetSignedImpParXYSignificanceudsgSecond; //!
	TH2D * fh2dJetSignedImpParXYSignificancebSecond; //!
	TH2D * fh2dJetSignedImpParXYSignificancebSecondNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYSignificancecSecond; //!

	TH2D * fh2dJetSignedImpParXYZSecond; //!
	TH2D * fh2dJetSignedImpParXYZUnidentifiedSecond; //!
	TH2D * fh2dJetSignedImpParXYZudsgSecond; //!
	TH2D * fh2dJetSignedImpParXYZbSecond; //!
	TH2D * fh2dJetSignedImpParXYZbSecondNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYZcSecond; //!

	TH2D * fh2dJetSignedImpParXYZSignificanceSecond; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceudsgSecond; //!
	TH2D * fh2dJetSignedImpParXYZSignificancebSecond; //!
	TH2D * fh2dJetSignedImpParXYZSignificancebSecondNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYZSignificancecSecond; //!
	//Third
	TH2D * fh2dJetSignedImpParXYThird; //!
	TH2D * fh2dJetSignedImpParXYUnidentifiedThird; //!
	TH2D * fh2dJetSignedImpParXYudsgThird; //!
	TH2D * fh2dJetSignedImpParXYbThird; //!
	TH2D * fh2dJetSignedImpParXYbThirdNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYcThird; //!

	TH2D * fh2dJetSignedImpParXYSignificanceThird; //!
	TH2D * fh2dJetSignedImpParXYSignificanceUnidentifiedThird; //!
	TH2D * fh2dJetSignedImpParXYSignificanceudsgThird; //!
	TH2D * fh2dJetSignedImpParXYSignificancebThird; //!
	TH2D * fh2dJetSignedImpParXYSignificancebThirdNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYSignificancecThird; //!

	TH2D * fh2dJetSignedImpParXYZThird; //!
	TH2D * fh2dJetSignedImpParXYZUnidentifiedThird; //!
	TH2D * fh2dJetSignedImpParXYZudsgThird; //!
	TH2D * fh2dJetSignedImpParXYZbThird; //!
	TH2D * fh2dJetSignedImpParXYZbThirdNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYZcThird; //!

	TH2D * fh2dJetSignedImpParXYZSignificanceThird; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceUnidentifiedThird; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceudsgThird; //!
	TH2D * fh2dJetSignedImpParXYZSignificancebThird; //!
	TH2D * fh2dJetSignedImpParXYZSignificancebThirdNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYZSignificancecThird; //!

	// inclusive signed impact parameter distributions corrected for light mesons (TODO D mesons) in Monte Carlo
	//First
	TH2D * fh2dJetSignedImpParXYFirst_McCorr; //!
	TH2D * fh2dJetSignedImpParXYUnidentifiedFirst_McCorr; //!
	TH2D * fh2dJetSignedImpParXYudsgFirst_McCorr; //!
	TH2D * fh2dJetSignedImpParXYbFirst_McCorr; //!
	TH2D * fh2dJetSignedImpParXYbFirst_McCorrNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYcFirst_McCorr; //!

	TH2D * fh2dJetSignedImpParXYSignificanceFirst_McCorr; //!
	TH2D * fh2dJetSignedImpParXYSignificanceUnidentifiedFirst_McCorr; //!
	TH2D * fh2dJetSignedImpParXYSignificanceudsgFirst_McCorr; //!
	TH2D * fh2dJetSignedImpParXYSignificancebFirst_McCorr; //!
	TH2D * fh2dJetSignedImpParXYSignificancebFirst_McCorrNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYSignificancecFirst_McCorr; //!

	TH2D * fh2dJetSignedImpParXYZFirst_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZUnidentifiedFirst_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZudsgFirst_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZbFirst_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZbFirst_McCorrNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYZcFirst_McCorr; //!

	TH2D * fh2dJetSignedImpParXYZSignificanceFirst_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceudsgFirst_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZSignificancebFirst_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZSignificancebFirst_McCorrNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYZSignificancecFirst_McCorr; //!
	//Second
	TH2D * fh2dJetSignedImpParXYSecond_McCorr; //!
	TH2D * fh2dJetSignedImpParXYUnidentifiedSecond_McCorr; //!
	TH2D * fh2dJetSignedImpParXYudsgSecond_McCorr; //!
	TH2D * fh2dJetSignedImpParXYbSecond_McCorr; //!
	TH2D * fh2dJetSignedImpParXYbSecond_McCorrNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYcSecond_McCorr; //!

	TH2D * fh2dJetSignedImpParXYSignificanceSecond_McCorr; //!
	TH2D * fh2dJetSignedImpParXYSignificanceUnidentifiedSecond_McCorr; //!
	TH2D * fh2dJetSignedImpParXYSignificanceudsgSecond_McCorr; //!
	TH2D * fh2dJetSignedImpParXYSignificancebSecond_McCorr; //!
	TH2D * fh2dJetSignedImpParXYSignificancebSecond_McCorrNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYSignificancecSecond_McCorr; //!

	TH2D * fh2dJetSignedImpParXYZSecond_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZUnidentifiedSecond_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZudsgSecond_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZbSecond_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZbSecond_McCorrNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYZcSecond_McCorr; //!

	TH2D * fh2dJetSignedImpParXYZSignificanceSecond_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceudsgSecond_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZSignificancebSecond_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZSignificancebSecond_McCorrNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYZSignificancecSecond_McCorr; //!
	//Third
	TH2D * fh2dJetSignedImpParXYThird_McCorr; //!
	TH2D * fh2dJetSignedImpParXYUnidentifiedThird_McCorr; //!
	TH2D * fh2dJetSignedImpParXYudsgThird_McCorr; //!
	TH2D * fh2dJetSignedImpParXYbThird_McCorr; //!
	TH2D * fh2dJetSignedImpParXYbThird_McCorrNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYcThird_McCorr; //!

	TH2D * fh2dJetSignedImpParXYSignificanceThird_McCorr; //!
	TH2D * fh2dJetSignedImpParXYSignificanceUnidentifiedThird_McCorr; //!
	TH2D * fh2dJetSignedImpParXYSignificanceudsgThird_McCorr; //!
	TH2D * fh2dJetSignedImpParXYSignificancebThird_McCorr; //!
	TH2D * fh2dJetSignedImpParXYSignificancebThird_McCorrNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYSignificancecThird_McCorr; //!

	TH2D * fh2dJetSignedImpParXYZThird_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZUnidentifiedThird_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZudsgThird_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZbThird_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZbThird_McCorrNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYZcThird_McCorr; //!

	TH2D * fh2dJetSignedImpParXYZSignificanceThird_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceUnidentifiedThird_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceudsgThird_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZSignificancebThird_McCorr; //!
	TH2D * fh2dJetSignedImpParXYZSignificancebThird_McCorrNonBDecay; //!
	TH2D * fh2dJetSignedImpParXYZSignificancecThird_McCorr; //!


	///////////// Plots for electron contribution
	TH1D * fh1dTracksImpParXY_electron;//! R Impact Parameter
	TH1D * fh1dTracksImpParXYZ_electron;//! R+z Impact Parameter
	TH1D * fh1dTracksImpParXYSignificance_electron;//! R Impact Parameter Significance
	TH1D * fh1dTracksImpParXYZSignificance_electron;//! R+z Impact Parameter Significance

	TH2D * fh2dJetSignedImpParXY_electron; //!
	TH2D * fh2dJetSignedImpParXYUnidentified_electron; //!
	TH2D * fh2dJetSignedImpParXYudsg_electron; //!
	TH2D * fh2dJetSignedImpParXYb_electron; //!
	TH2D * fh2dJetSignedImpParXYc_electron; //!

	TH2D * fh2dJetSignedImpParXYSignificance_electron; //!
	TH2D * fh2dJetSignedImpParXYSignificanceUnidentified_electron; //!
	TH2D * fh2dJetSignedImpParXYSignificanceudsg_electron; //!
	TH2D * fh2dJetSignedImpParXYSignificanceb_electron; //!
	TH2D * fh2dJetSignedImpParXYSignificancec_electron; //!

	TH2D * fh2dJetSignedImpParXYZ_electron; //!
	TH2D * fh2dJetSignedImpParXYZUnidentified_electron; //!
	TH2D * fh2dJetSignedImpParXYZudsg_electron; //!
	TH2D * fh2dJetSignedImpParXYZb_electron; //!
	TH2D * fh2dJetSignedImpParXYZc_electron; //!

	TH2D * fh2dJetSignedImpParXYZSignificance_electron; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceUnidentified_electron; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceudsg_electron; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceb_electron; //!
	TH2D * fh2dJetSignedImpParXYZSignificancec_electron; //!


	// inclusive signed impact parameter distributions
	//First
	TH2D * fh2dJetSignedImpParXYFirst_electron; //!
	TH2D * fh2dJetSignedImpParXYUnidentifiedFirst_electron; //!
	TH2D * fh2dJetSignedImpParXYudsgFirst_electron; //!
	TH2D * fh2dJetSignedImpParXYbFirst_electron; //!
	TH2D * fh2dJetSignedImpParXYcFirst_electron; //!

	TH2D * fh2dJetSignedImpParXYSignificanceFirst_electron; //!
	TH2D * fh2dJetSignedImpParXYSignificanceUnidentifiedFirst_electron; //!
	TH2D * fh2dJetSignedImpParXYSignificanceudsgFirst_electron; //!
	TH2D * fh2dJetSignedImpParXYSignificancebFirst_electron; //!
	TH2D * fh2dJetSignedImpParXYSignificancecFirst_electron; //!

	TH2D * fh2dJetSignedImpParXYZFirst_electron; //!
	TH2D * fh2dJetSignedImpParXYZUnidentifiedFirst_electron; //!
	TH2D * fh2dJetSignedImpParXYZudsgFirst_electron; //!
	TH2D * fh2dJetSignedImpParXYZbFirst_electron; //!
	TH2D * fh2dJetSignedImpParXYZcFirst_electron; //!

	TH2D * fh2dJetSignedImpParXYZSignificanceFirst_electron; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst_electron; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceudsgFirst_electron; //!
	TH2D * fh2dJetSignedImpParXYZSignificancebFirst_electron; //!
	TH2D * fh2dJetSignedImpParXYZSignificancecFirst_electron; //!
	//Second
	TH2D * fh2dJetSignedImpParXYSecond_electron; //!
	TH2D * fh2dJetSignedImpParXYUnidentifiedSecond_electron; //!
	TH2D * fh2dJetSignedImpParXYudsgSecond_electron; //!
	TH2D * fh2dJetSignedImpParXYbSecond_electron; //!
	TH2D * fh2dJetSignedImpParXYcSecond_electron; //!

	TH2D * fh2dJetSignedImpParXYSignificanceSecond_electron; //!
	TH2D * fh2dJetSignedImpParXYSignificanceUnidentifiedSecond_electron; //!
	TH2D * fh2dJetSignedImpParXYSignificanceudsgSecond_electron; //!
	TH2D * fh2dJetSignedImpParXYSignificancebSecond_electron; //!
	TH2D * fh2dJetSignedImpParXYSignificancecSecond_electron; //!

	TH2D * fh2dJetSignedImpParXYZSecond_electron; //!
	TH2D * fh2dJetSignedImpParXYZUnidentifiedSecond_electron; //!
	TH2D * fh2dJetSignedImpParXYZudsgSecond_electron; //!
	TH2D * fh2dJetSignedImpParXYZbSecond_electron; //!
	TH2D * fh2dJetSignedImpParXYZcSecond_electron; //!

	TH2D * fh2dJetSignedImpParXYZSignificanceSecond_electron; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond_electron; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceudsgSecond_electron; //!
	TH2D * fh2dJetSignedImpParXYZSignificancebSecond_electron; //!
	TH2D * fh2dJetSignedImpParXYZSignificancecSecond_electron; //!
	//Third
	TH2D * fh2dJetSignedImpParXYThird_electron; //!
	TH2D * fh2dJetSignedImpParXYUnidentifiedThird_electron; //!
	TH2D * fh2dJetSignedImpParXYudsgThird_electron; //!
	TH2D * fh2dJetSignedImpParXYbThird_electron; //!
	TH2D * fh2dJetSignedImpParXYcThird_electron; //!

	TH2D * fh2dJetSignedImpParXYSignificanceThird_electron; //!
	TH2D * fh2dJetSignedImpParXYSignificanceUnidentifiedThird_electron; //!
	TH2D * fh2dJetSignedImpParXYSignificanceudsgThird_electron; //!
	TH2D * fh2dJetSignedImpParXYSignificancebThird_electron; //!
	TH2D * fh2dJetSignedImpParXYSignificancecThird_electron; //!

	TH2D * fh2dJetSignedImpParXYZThird_electron; //!
	TH2D * fh2dJetSignedImpParXYZUnidentifiedThird_electron; //!
	TH2D * fh2dJetSignedImpParXYZudsgThird_electron; //!
	TH2D * fh2dJetSignedImpParXYZbThird_electron; //!
	TH2D * fh2dJetSignedImpParXYZcThird_electron; //!

	TH2D * fh2dJetSignedImpParXYZSignificanceThird_electron; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceUnidentifiedThird_electron; //!
	TH2D * fh2dJetSignedImpParXYZSignificanceudsgThird_electron; //!
	TH2D * fh2dJetSignedImpParXYZSignificancebThird_electron; //!
	TH2D * fh2dJetSignedImpParXYZSignificancecThird_electron; //!
	TClonesArray * fMCArray;//!
	AliMCEvent   * fMCEvent;//!;
	AliESDtrackCuts  *fESDTrackCut;// copy to manager - root from AddTask
	AliAnalysisUtils *fUtils;//!
	//Monte Carlo correction factor containers

	Double_t fBackgroundFactor[9][44];//[9][44]
	Double_t fBackgroundFactorBins[45];//[45]
	Double_t fBackgroundFactorLinus[21][498]; //[21][498]FineBinned correction factors up 0.1-25 GeV/c first value below last above 0.05 binwidth
	static bool mysort(const myvaluetuple& i, const myvaluetuple& j);

	ClassDef(AliAnalysisTaskHFJetIPQA, 4	)
};
#endif




 //
