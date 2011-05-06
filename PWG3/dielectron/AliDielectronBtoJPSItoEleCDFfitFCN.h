#ifndef ALIDIELECTRONBTOJPSITOELECDFFITFCN_H
#define ALIDIELECTRONBTOJPSITOELECDFFITFCN_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
//                      Class AliDielectronBtoJPSItoEleCDFfitFCN
//                    Definition of main function used in 
//                     unbinned log-likelihood fit for
//                 the channel B -> JPsi + X -> e+e- + X
//      
//                          Origin: C.Di Giglio
//     Contact: Carmelo.Digiglio@ba.infn.it , Giuseppe.Bruno@ba.infn.it
//_________________________________________________________________________

#include <TNamed.h>
#include <TDatabasePDG.h>
#include <TH1F.h>

class TRandom3;
class TF1;

enum IntegralType {kBkg, 
	kBkgNorm, 
	kSig, 
	kSigNorm};

enum PtBins       {kallpt, 
	kptbin1,kptbin2,kptbin3,
	kptbin4,kptbin5,kptbin6,
	kptbin7,kptbin8,kptbin9};
//_________________________________________________________________________________________________
class AliDielectronBtoJPSItoEleCDFfitFCN : public TNamed {
	public:
		//
		AliDielectronBtoJPSItoEleCDFfitFCN();
		AliDielectronBtoJPSItoEleCDFfitFCN(const AliDielectronBtoJPSItoEleCDFfitFCN& source); 
		AliDielectronBtoJPSItoEleCDFfitFCN& operator=(const AliDielectronBtoJPSItoEleCDFfitFCN& source);  
		virtual ~AliDielectronBtoJPSItoEleCDFfitFCN();

		Double_t EvaluateLikelihood(const Double_t* pseudoproperdecaytime,
				const Double_t* invariantmass, const Int_t ncand) const;

		Double_t GetResWeight()                    const { return fParameters[0]; }
		Double_t GetFPlus()                        const { return fParameters[1]; }
		Double_t GetFMinus()                       const { return fParameters[2]; }
		Double_t GetFSym()                         const { return fParameters[3]; } 
		Double_t GetLamPlus()                      const { return fParameters[4]; }
		Double_t GetLamMinus()                     const { return fParameters[5]; }
		Double_t GetLamSym()                       const { return fParameters[6]; }
		Double_t GetFractionJpsiFromBeauty()       const { return fParameters[7]; }
		Double_t GetFsig()                         const { return fParameters[8]; }
		Double_t GetCrystalBallMmean()             const { return fParameters[9]; }
		Double_t GetCrystalBallNexp()              const { return fParameters[10]; }
		Double_t GetCrystalBallSigma()             const { return fParameters[11]; }
		Double_t GetCrystalBallAlpha()             const { return fParameters[12]; }
		Double_t GetCrystalBallNorm()              const { return fParameters[13]; }
		Double_t GetBkgInvMassNorm()               const { return fParameters[14]; }
		Double_t GetBkgInvMassMean()               const { return fParameters[15]; }
		Double_t GetBkgInvMassSlope()              const { return fParameters[16]; }  
		Double_t GetBkgInvMassConst()              const { return fParameters[17]; } 
		Double_t GetNormGaus1ResFunc()             const { return fParameters[18]; }
		Double_t GetNormGaus2ResFunc()             const { return fParameters[19]; }
		Double_t GetIntegralMassSig()              const { return fintmMassSig; }
		Double_t GetIntegralMassBkg()              const { return fintmMassBkg; }
		const Double_t* GetResolutionConstants()   const { return fResolutionConstants; }
		Bool_t GetCrystalBallParam()               const { return fCrystalBallParam; }
		TH1F * GetCsiMcHisto()                     const { return fhCsiMC; }

		// return pointer to likelihood functions  
		TF1* GetCsiMC(Double_t xmin, Double_t xmax);
		TF1* GetResolutionFunc(Double_t xmin, Double_t xmax);
		TF1* GetEvaluateCDFDecayTimeBkgDistr(Double_t xmin, Double_t xmax);
		TF1* GetEvaluateCDFDecayTimeSigDistr(Double_t xmin, Double_t xmax);

		void SetResWeight (Double_t resWgt) {fParameters[0] = resWgt;}
		void SetFPlus(Double_t plus) {fParameters[1] = plus;}
		void SetFMinus(Double_t minus) {fParameters[2]  = minus;}
		void SetFSym(Double_t sym) {fParameters[3] = sym;}
		void SetLamPlus(Double_t lamplus) {fParameters[4] = lamplus;}
		void SetLamMinus(Double_t lamminus) {fParameters[5] = lamminus;}
		void SetLamSym(Double_t lamsym) {fParameters[6] = lamsym;}
		void SetFractionJpsiFromBeauty(Double_t B) {fParameters[7] = B;}
		void SetFsig(Double_t Fsig) {fParameters[8] = Fsig;}
		void SetCrystalBallMmean(Double_t CrystalBallMmean) {fParameters[9] = CrystalBallMmean;}
		void SetCrystalBallNexp(Double_t CrystalBallNexp) {fParameters[10] = CrystalBallNexp;}
		void SetCrystalBallSigma(Double_t CrystalBallSigma) {fParameters[11] = CrystalBallSigma;}
		void SetCrystalBallAlpha(Double_t CrystalBallAlpha) {fParameters[12] = CrystalBallAlpha;}
		void SetCrystalBallNorm(Double_t CrystalBallNorm) {fParameters[13] = CrystalBallNorm;}
		void SetBkgInvMassNorm(Double_t BkgInvMassNorm) {fParameters[14] = BkgInvMassNorm;}
		void SetBkgInvMassMean(Double_t BkgInvMassMean) {fParameters[15] = BkgInvMassMean;}
		void SetBkgInvMassSlope(Double_t BkgInvMassSlope) {fParameters[16] = BkgInvMassSlope;}
		void SetBkgInvMassConst(Double_t BkgInvMassConst) {fParameters[17] = BkgInvMassConst;}
		void SetNormGaus1ResFunc(Double_t norm1) {fParameters[18] = norm1;}
		void SetNormGaus2ResFunc(Double_t norm2) {fParameters[19] = norm2;}
		void SetAllParameters(const Double_t* parameters);
		void SetIntegralMassSig(Double_t integral) { fintmMassSig = integral; }
		void SetIntegralMassBkg(Double_t integral) { fintmMassBkg = integral; }
		void SetCsiMC(const TH1F* MCtemplate) {fhCsiMC = (TH1F*)MCtemplate->Clone("fhCsiMC");}

		void SetResolutionConstants(Double_t* resolutionConst);
		void SetMassWndHigh(Double_t limit) { fMassWndHigh = TDatabasePDG::Instance()->GetParticle(443)->Mass() + limit ;}
		void SetMassWndLow(Double_t limit) { fMassWndLow = TDatabasePDG::Instance()->GetParticle(443)->Mass() - limit ;}
		void SetCrystalBallFunction(Bool_t okCB) {fCrystalBallParam = okCB;}

		void ComputeMassIntegral(); 

		void ReadMCtemplates(Int_t BinNum);

		void PrintStatus();

	private:  
		//
		Double_t fParameters[20];        /*  par[0]  = weightRes;                
						     par[1]  = fPos;
						     par[2]  = fNeg;
						     par[3]  = fSym
						     par[4]  = fOneOvLamPlus;
						     par[5]  = fOneOvLamMinus;
						     par[6]  = fOneOvLamSym;
						     par[7]  = fFractionJpsiFromBeauty;
						     par[8]  = fFsig;
						     par[9]  = fCrystalBallMmean;
						     par[10] = fCrystalBallNexp;
						     par[11] = fCrystalBallSigma;
						     par[12] = fCrystalBallAlpha;
						     par[13] = fCrystalBallNorm;
						     par[14] = fBkgNorm;
						     par[15] = fBkgMean; 
						     par[16] = fBkgSlope;
						     par[17] = fBkgConst;
						     par[18] = norm1Gaus;
						     par[19] = norm2Gaus;*/


		Double_t fFPlus;                     // parameters of the log-likelihood function
		Double_t fFMinus;                    // Slopes of the x distributions of the background 
		Double_t fFSym;                      // functions 

		Double_t fintmMassSig;               // integral of invariant mass distribution for the signal
		Double_t fintmMassBkg;               // integral of invariant mass distribution for the bkg

		TH1F *fhCsiMC;                       // X distribution used as MC template for JPSI from B
		Double_t fResolutionConstants[4];    // constants for the parametrized resolution function R(X)
		Double_t fMassWndHigh;               // JPSI Mass window higher limit
		Double_t fMassWndLow;                // JPSI Mass window lower limit
		Bool_t fCrystalBallParam;            // Boolean to switch to Crystall Ball parameterisation

		////

		Double_t EvaluateCDFfunc(Double_t x, Double_t m) const ;
		Double_t EvaluateCDFfuncNorm(Double_t x, Double_t m) const ;

		////

		Double_t EvaluateCDFfuncSignalPart(Double_t x, Double_t m) const ;      // Signal part 
		Double_t EvaluateCDFDecayTimeSigDistr(Double_t x) const ;
		Double_t EvaluateCDFDecayTimeSigDistrFunc(const Double_t* x, const Double_t *par) const { par = 0x0;  return EvaluateCDFDecayTimeSigDistr(x[0]);}
		Double_t EvaluateCDFInvMassSigDistr(Double_t m) const ;
		Double_t EvaluateCDFfuncBkgPart(Double_t x,Double_t m) const ;          // Background part
		Double_t EvaluateCDFDecayTimeBkgDistr(Double_t x) const ;
		Double_t EvaluateCDFDecayTimeBkgDistrFunc(const Double_t* x, const Double_t *par) const { par = 0x0;  return EvaluateCDFDecayTimeBkgDistr(x[0]);}
		Double_t EvaluateCDFInvMassBkgDistr(Double_t m) const ;

		////

		Double_t FunB(Double_t x) const;
		Double_t FunP(Double_t x) const ;
		Double_t CsiMC(Double_t x) const;
		Double_t CsiMCfunc(const Double_t* x, const Double_t *par) const {  par = 0x0; return CsiMC(x[0]);}
		Double_t FunBkgPos(Double_t x) const ;
		Double_t FunBkgNeg(Double_t x) const ;
		Double_t FunBkgSym(Double_t x) const ;
		Double_t ResolutionFunc(Double_t x) const ;
		Double_t ResolutionFuncf(const Double_t* x, const Double_t *par) const { par = 0x0;  return ResolutionFunc(x[0]);}

		ClassDef (AliDielectronBtoJPSItoEleCDFfitFCN,1);         // Unbinned log-likelihood fit 

};

#endif
