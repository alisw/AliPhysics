#ifndef AliAnalysisTaskPsi3HBTsystFitrange_cxx
#define AliAnalysisTaskPsi3HBTsystFitrange_cxx


class TH1F;
class TH2F;
class TH3F;
class AliAODEvent;

#include "AliAnalysisTaskSE.h"
#include <vector>

class AliAnalysisTaskPsi3HBTsystFitrange : public AliAnalysisTaskSE {
	public:
		AliAnalysisTaskPsi3HBTsystFitrange();
		AliAnalysisTaskPsi3HBTsystFitrange(const char *name);
		virtual ~AliAnalysisTaskPsi3HBTsystFitrange() {;}

		virtual void   UserCreateOutputObjects()		;
		virtual void   UserExec(Option_t *option)		;
		virtual void   Terminate(Option_t *)				;


		void SetCentralityCentral(Bool_t setcentral)		{ fCentralAna = setcentral; };
		void SetCollisionCandidates(Int_t settrigger)		{ fTrigSel = settrigger; };

		Int_t TPCeta(Float_t ieta)									;
		Int_t TPCeta2(Float_t ieta)									;
		Float_t PhiZDC(Int_t itow)									;
		Int_t	trig_dec(Bool_t ikMB, Bool_t ikSemi, Bool_t ikCent);
		Int_t	HistCent(Double_t fV0Cent)						;
		Int_t	CoulCent(Int_t ftV0Cent)							;
		Int_t	PiEPAngle(Double_t fPlane)						;

		Double_t	qinv_cms(double& px1, double& px2, double& py1, double& py2, double& pz1, double& pz2, double& pe1, double& pe2)			;
		Double_t	qout_cms(double& x1, double& x2, double& y1, double& y2)																															;
		Double_t	qside_cms(double& x1, double& x2, double& y1, double& y2)																														;
		Double_t	qlong_cms(double& pz1, double& pz2, double& pe1, double& pe2)																												;

		int MagDec(Float_t fMag_t)									;
		Float_t		dec_dist(Float_t fDist)						;

		//Qinv -> 3D Q conversion
		int conv_if(Int_t qoIn, Int_t qsIn, Int_t qlIn)											;
		int conv_if2(Int_t qoIn, Int_t qsIn, Int_t qlIn)										;
		int conv_over(double qoIn, double qsIn, double qlIn)	;
		int conv_qosl(Double_t qIn)																					;
		int conv_qosl2(Double_t qIn)																				;
		int conv_qosl3(Double_t qIn)																				;

		//Calculate two track resolution
		float Sharity(TBits& cl1,  TBits& cl2, TBits& sh1, TBits& sh2)	;
		float Qfac(TBits& cl1,  TBits& cl2, TBits& sh1, TBits& sh2)			;

		//PID Function
		bool IsPionNSigma(Double_t& mom, Float_t& nsigmaTPCPi, Float_t& nsigmaTOFPi)			;
		bool IsKaonNSigma(Double_t& mom, Float_t& nsigmaTPCPi, Float_t& nsigmaTOFPi)			;
		bool IsProtonNSigma(Double_t& mom, Float_t& nsigmaTPCPi, Float_t& nsigmaTOFPi)		;


		Float_t		EtaLimitFMDC(Int_t i)							;
		Float_t		EtaLimitFMDA(Int_t i)							;
		Float_t		EtaLimitFMDAC1(Int_t i)						;
		Float_t		EtaLimitFMDAC2(Int_t i)						;

		TH1F* ProjectPhiFMD(TH2D* fmd, Float_t eta1, Float_t eta2) ;   //don't forget delete projected hist

		void FillEPCorr(Int_t fV0flag, Int_t fTPCflag, Int_t fFMDflag);
		void PairCalc_Final()																					;

		static const int fCorr	= 2		;
		static const int fTrigH	= 3		;
		static const int fV0s		= 3		;
		static const int fTPCs	= 7		;
		static const int fFMDs	= 13	;
		static const int fHarm	= 5		;

		static const int fOrd		= 8		;
		static const int fDet		= 39	;
		static const int fVnDet	= 3		;
		static const int fVnCent= 18	;
		static const int fQvCent= 9		;

		static const int fChargeH	= 2		;
		static const int fCentH		= 4		;
		static const int fCentM		= 10	;
		static const int fZvtxM		= 8		;
		static const int fPsi3M		= 20	;
		static const int fMix			= 10	;
		static const int fktH			= 2		;
		static const int fPsi3H		= 8		;

		static const int fQbin		= 2		;


	private:

		AliAnalysisTaskPsi3HBTsystFitrange(const AliAnalysisTaskPsi3HBTsystFitrange&); // not implemented
		AliAnalysisTaskPsi3HBTsystFitrange& operator=(const AliAnalysisTaskPsi3HBTsystFitrange&); // not implemented

		AliAODEvent	*fAOD						; //! AOD object
		TList		*fOutputList				; //! Output list

		TList		*fHListEPCalib			;
		TFile		*fFileEPCalib				;

		Bool_t	fCentralAna					;
		Int_t		fTrigSel						;

		Int_t		fRunNumber							;
		Int_t		fRunNumberForPrevEvent	;

		Int_t		fNEntries, fNumberE			;

		Double_t	cent_V0M, cent_Calib				;
		Int_t			fV0flag, fTPCflag, fFMDflag	;
		Int_t			fMagsign										;

		TH1F				*fH1Event				;
		TH1F				*fH1Nch					;

		TH1F				*fH1Phi					;
		TH1F				*fH1Eta					;
		TH1F				*fH1Pt					;

		TProfile		*fPrfV0Sig			;
		TH2F				*fH2V0Sig				;

		TH2D				*fFMD						;
		TH2F				*fH2FMD2D				;

		TH1F				*fH1kt					;
		TH1F				*fH1kt_all			;

		TH2F				*fH2DCA					;
		TH2F				*fH2DCAxyPt			;
		TH2F				*fH2DCAzPt 			;
		TH2F				*fH2PhiEta			;

		TH2F				*fH2TOF					;

		AliAnaCoulWave	*fCoulomb		;

		TH2F				*fH2NsigTPC	[3]	;
		TH2F				*fH2NsigTOF	[3]	;

		TH1F				*fH1Cent		[4]	;
		TH1F				*fH1Zvtx		[4]	;
		TH1F				*fH1EtaD		[7]	;
		TH1F				*fH1kt_div	[4]	;
    
		TProfile		*fPrfV0Qvcos	[fCorr][fTrigH][fV0s]	[fHarm]	;
		TProfile		*fPrfV0Qvsin	[fCorr][fTrigH][fV0s]	[fHarm]	;
		TProfile		*fPrfTPCQvcos	[fCorr][fTrigH][fTPCs][fHarm]	;
		TProfile		*fPrfTPCQvsin	[fCorr][fTrigH][fTPCs][fHarm]	;
		TProfile		*fPrfFMDQvcos	[fCorr][fTrigH][fFMDs][fHarm]	;
		TProfile		*fPrfFMDQvsin	[fCorr][fTrigH][fFMDs][fHarm]	;
    
		//[correction][trigger][side][harmonics][centrality]
		TH1F		*fH1V0Qvcos		[fCorr][fTrigH][fV0s]	[fHarm][fQvCent]	;
		TH1F		*fH1V0Qvsin		[fCorr][fTrigH][fV0s]	[fHarm][fQvCent]	;
		TH1F		*fH1TPCQvcos	[fCorr][fTrigH][fTPCs][fHarm][fQvCent]	;
		TH1F		*fH1TPCQvsin	[fCorr][fTrigH][fTPCs][fHarm][fQvCent]	;
		TH1F		*fH1FMDQvcos	[fCorr][fTrigH][fFMDs][fHarm][fQvCent]	;
		TH1F		*fH1FMDQvsin	[fCorr][fTrigH][fFMDs][fHarm][fQvCent]	;
    
		//[correction][trigger][side][harmonics][cent]
		TH1F				*fH1EPV0		[fCorr][fTrigH][fV0s]	[fHarm][fQvCent]	;
		TH1F				*fH1EPTPC		[fCorr][fTrigH][fTPCs][fHarm][fQvCent]	;
		TH1F				*fH1EPFMD		[fCorr][fTrigH][fFMDs][fHarm][fQvCent]	;
    
		//[trigger][side][harmonics]
		TProfile		*fEPCalibV0Qvx	[fTrigH][fV0s]	[fHarm]		;
		TProfile		*fEPCalibV0Qvy	[fTrigH][fV0s]	[fHarm]		;
		TProfile		*fEPCalibTPCQvx	[fTrigH][fTPCs]	[fHarm]		;
		TProfile		*fEPCalibTPCQvy	[fTrigH][fTPCs]	[fHarm]		;
		TProfile		*fEPCalibFMDQvx	[fTrigH][fFMDs]	[fHarm]		;
		TProfile		*fEPCalibFMDQvy	[fTrigH][fFMDs]	[fHarm]		;
    
		TProfile		*fPrfCos	[fDet][fHarm];
		TProfile		*fPrfSin	[fDet][fHarm];

		TProfile		*fPrfvncos_cent	[fVnDet][fV0s][fHarm]					;
		TProfile		*fPrfvnsin_cent	[fVnDet][fV0s][fHarm]					;
		TProfile		*fPrfvncos_pt		[fVnDet][fV0s][fHarm][fVnCent];
		TProfile		*fPrfvnsin_pt		[fVnDet][fV0s][fHarm][fVnCent];

		//[charge][Centrality][Psi3][Qbin]
		TH1D				*fH1Qinv		[fChargeH][fCentH][fPsi3H][fQbin]		;
		TH1D				*fH1CQinv		[fChargeH][fCentH][fPsi3H][fQbin]		;
		TH1D				*fH1Qinv_mix[fChargeH][fCentH][fPsi3H][fQbin]		;
		TH3D				*fH3Q				[fChargeH][fCentH][fPsi3H][fQbin]		;
		TH3D				*fH3Q_mix		[fChargeH][fCentH][fPsi3H][fQbin]		;
		TProfile		*fPrfConv		[fChargeH][fCentH][fPsi3H][fQbin]		;

		//[trigger][side][harmonics][centrality]
		Float_t	MeanQxV0	[fTrigH][fV0s]	[fHarm][fVnCent], MeanQyV0	[fTrigH][fV0s]	[fHarm][fVnCent];
		Float_t	RMSQxV0		[fTrigH][fV0s]	[fHarm][fVnCent], RMSQyV0		[fTrigH][fV0s]	[fHarm][fVnCent];
		Float_t	MeanQxTPC	[fTrigH][fTPCs]	[fHarm][fVnCent], MeanQyTPC	[fTrigH][fTPCs]	[fHarm][fVnCent];
		Float_t	RMSQxTPC	[fTrigH][fTPCs]	[fHarm][fVnCent], RMSQyTPC	[fTrigH][fTPCs]	[fHarm][fVnCent];
		Float_t	MeanQxFMD	[fTrigH][fFMDs]	[fHarm][fVnCent], MeanQyFMD	[fTrigH][fFMDs]	[fHarm][fVnCent];
		Float_t	RMSQxFMD	[fTrigH][fFMDs]	[fHarm][fVnCent], RMSQyFMD	[fTrigH][fFMDs]	[fHarm][fVnCent];

		//[correction][side][harmonics]
		Double_t	QxV0	[fCorr][fV0s]	[fHarm],	QyV0	[fCorr][fV0s]	[fHarm];
		Double_t	QxTPC	[fCorr][fTPCs][fHarm],	QyTPC	[fCorr][fTPCs][fHarm];
		Double_t	QxFMD	[fCorr][fFMDs][fHarm],	QyFMD	[fCorr][fFMDs][fHarm];

		//													[correction][side][harmonics]
		Double_t	fV0psi	[fCorr][fV0s]	[fHarm];
		Double_t	fTPCpsi	[fCorr][fTPCs][fHarm];
		Double_t	fFMDpsi	[fCorr][fFMDs][fHarm];

		//Mixing Parameters
		/**********************************************/
		//Pion event mixing vector
		//[Centrality][Zvertex][Psi3 angle][Event Pool#]
		std::vector<double>	pion_ppx	[fCentM][fZvtxM][fPsi3M][fMix];
		std::vector<double>	pion_ppy	[fCentM][fZvtxM][fPsi3M][fMix];
		std::vector<double>	pion_ppz	[fCentM][fZvtxM][fPsi3M][fMix];
		std::vector<double>	pion_ppt	[fCentM][fZvtxM][fPsi3M][fMix];
		std::vector<double>	pion_pe		[fCentM][fZvtxM][fPsi3M][fMix];
		std::vector<double>	pion_pphi	[fCentM][fZvtxM][fPsi3M][fMix];
		std::vector<double>	pion_peta	[fCentM][fZvtxM][fPsi3M][fMix];
		std::vector<double>	pion_pax	[fCentM][fZvtxM][fPsi3M][fMix];
		std::vector<TBits>	pion_pclu	[fCentM][fZvtxM][fPsi3M][fMix];
		std::vector<TBits>	pion_psha	[fCentM][fZvtxM][fPsi3M][fMix];
		//=================================================
		std::vector<double>	pion_mpx	[fCentM][fZvtxM][fPsi3M][fMix];
		std::vector<double>	pion_mpy	[fCentM][fZvtxM][fPsi3M][fMix];
		std::vector<double>	pion_mpz	[fCentM][fZvtxM][fPsi3M][fMix];
		std::vector<double>	pion_mpt	[fCentM][fZvtxM][fPsi3M][fMix];
		std::vector<double>	pion_me		[fCentM][fZvtxM][fPsi3M][fMix];
		std::vector<double>	pion_mphi	[fCentM][fZvtxM][fPsi3M][fMix];
		std::vector<double>	pion_meta	[fCentM][fZvtxM][fPsi3M][fMix];
		std::vector<double>	pion_max	[fCentM][fZvtxM][fPsi3M][fMix];
		std::vector<TBits>	pion_mclu	[fCentM][fZvtxM][fPsi3M][fMix];
		std::vector<TBits>	pion_msha	[fCentM][fZvtxM][fPsi3M][fMix];
		/**********************************************/

		//Mixing Flag [charge][centrality][Zvertex][Psi3 angle]
		Int_t		mix_flag[fChargeH][fCentM][fZvtxM][fPsi3M];

		ClassDef(AliAnalysisTaskPsi3HBTsystFitrange, 1); // example of analysis
};

#endif
