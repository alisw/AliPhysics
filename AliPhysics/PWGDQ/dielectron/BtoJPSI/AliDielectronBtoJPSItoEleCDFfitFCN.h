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

//_________________________________________________________________________________________________
class AliDielectronBtoJPSItoEleCDFfitFCN : public TNamed {
	public:
		//
		AliDielectronBtoJPSItoEleCDFfitFCN();
		AliDielectronBtoJPSItoEleCDFfitFCN(const AliDielectronBtoJPSItoEleCDFfitFCN& source); 
		AliDielectronBtoJPSItoEleCDFfitFCN& operator=(const AliDielectronBtoJPSItoEleCDFfitFCN& source);  
		virtual ~AliDielectronBtoJPSItoEleCDFfitFCN();

		Double_t EvaluateLikelihood(const Double_t* pseudoproperdecaytime,
				const Double_t* invariantmass, const Double_t* pt, const Int_t* type, Int_t ncand) const;

                // getters for parameters
		Double_t GetResWeight()                    const { return fParameters[0]; }
		Double_t GetFPlus()                        const { return fParameters[1]; }
		Double_t GetFMinus()                       const { return fParameters[2]; }
		Double_t GetFSym()                         const { return fParameters[3]; } 
		Double_t GetFSym1()                        const { return fParameters[46]; }
		Double_t GetLamPlus()                      const { return fParameters[4]; }
		Double_t GetLamMinus()                     const { return fParameters[5]; }
		Double_t GetLamSym()                       const { return fParameters[6]; }
		Double_t GetLamSym1()                      const { return fParameters[45]; }
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
		Double_t GetNormGaus1ResFunc(Int_t type)   const { return fParameters[18+(2-type)*9]; }
		Double_t GetNormGaus2ResFunc(Int_t type)   const { return fParameters[19+(2-type)*9]; }
		Double_t GetIntegralMassSig()              const { return fintmMassSig; }
		Double_t GetIntegralMassBkg()              const { return fintmMassBkg; }
                Double_t GetResMean1(Int_t type)           const { return fParameters[20+(2-type)*9]; } 
                Double_t GetResSigma1(Int_t type)          const { return fParameters[21+(2-type)*9]; }
                Double_t GetResMean2(Int_t type)           const { return fParameters[22+(2-type)*9]; }
                Double_t GetResSigma2(Int_t type)          const { return fParameters[23+(2-type)*9]; }
                Double_t GetResAlfa(Int_t type)            const { return fParameters[24+(2-type)*9]; }   
                Double_t GetResLambda(Int_t type)          const { return fParameters[25+(2-type)*9]; } 
                Double_t GetResNormExp(Int_t type)         const { return fParameters[26+(2-type)*9]; }
                Double_t GetPolyn4()                       const { return fParameters[47]; }
                Double_t GetPolyn5()                       const { return fParameters[48]; }
                Bool_t GetCrystalBallParam()               const { return fCrystalBallParam; }
                Bool_t GetExponentialParam()               const { return fExponentialParam; }
		TH1F * GetCsiMcHisto()                     const { return fhCsiMC; }
                Double_t GetResWeight(Int_t iW)            const { return fWeightType[iW]; }

		// return pointer to likelihood functions  
		TF1* GetCsiMC(Double_t xmin, Double_t xmax,Double_t normalization);
		TF1* GetResolutionFunc(Double_t xmin, Double_t xmax,Double_t normalization, Double_t pt, Int_t type=2);
	        TF1* GetResolutionFuncAllTypes(Double_t xmin, Double_t xmax,Double_t normalization);
         	TF1* GetFunB(Double_t xmin, Double_t xmax, Double_t normalization, Double_t pt, Int_t type=2, Int_t npx = 5000);
                TF1* GetFunBAllTypes(Double_t xmin, Double_t xmax, Double_t normalization);
                TF1* GetEvaluateCDFDecayTimeBkgDistr(Double_t xmin, Double_t xmax, Double_t normalization, Int_t type = 2, Double_t mass = 3.09, Double_t pt = 200.,Int_t npx = 5000);
		TF1* GetEvaluateCDFDecayTimeBkgDistrAllTypes(Double_t xmin, Double_t xmax, Double_t normalization);
                TF1* GetEvaluateCDFDecayTimeSigDistr(Double_t xmin, Double_t xmax, Double_t normalization, Double_t type);
                TF1* GetEvaluateCDFInvMassBkgDistr(Double_t mMin, Double_t mMax, Double_t normalization);
                TF1* GetEvaluateCDFInvMassSigDistr(Double_t mMin, Double_t mMax, Double_t normalization);
                TF1* GetEvaluateCDFInvMassTotalDistr(Double_t mMin, Double_t mMax, Double_t normalization);
                TF1* GetEvaluateCDFDecayTimeTotalDistr(Double_t xMin, Double_t xMax, Double_t normalization,Double_t pt = 200., Int_t type=2);
                TF1 *GetEvaluateCDFDecayTimeTotalDistrAllTypes(Double_t xMin, Double_t xMax, Double_t normalization);

		void SetResWeight(Double_t resWgt) {fParameters[0] = resWgt;}
		void SetFPlus(Double_t plus) {fParameters[1] = plus;}
		void SetFMinus(Double_t minus) {fParameters[2]  = minus;}
		void SetFSym(Double_t sym) {fParameters[3] = sym;}
		void SetFSym1(Double_t sym) {fParameters[46] = sym;}
		void SetLamPlus(Double_t lamplus) {fParameters[4] = lamplus;}
		void SetLamMinus(Double_t lamminus) {fParameters[5] = lamminus;}
		void SetLamSym(Double_t lamsym) {fParameters[6] = lamsym;}
		void SetLamSym1(Double_t lamsym) {fParameters[45] = lamsym;}
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
                void SetBkgInvMassPolyn4(Double_t coeffPol4) {fParameters[47] = coeffPol4;}
                void SetBkgInvMassPolyn5(Double_t coeffPol5) {fParameters[48] = coeffPol5;}
		void SetNormGaus1ResFunc(Double_t norm1) {fParameters[18] = norm1;}
		void SetNormGaus2ResFunc(Double_t norm2) {fParameters[19] = norm2;}
		void SetAllParameters(const Double_t* parameters);
		void SetIntegralMassSig(Double_t integral) { fintmMassSig = integral; }
		void SetIntegralMassBkg(Double_t integral) { fintmMassBkg = integral; }
		void SetCsiMC(const TH1F* MCtemplate) {fhCsiMC = (TH1F*)MCtemplate->Clone("fhCsiMC");}
                void SetTemplateShift(Double_t shift = 0.){fShiftTemplate = shift;}
 
		void SetResolutionConstants(const Double_t* resolutionConst, Int_t type);
		void SetMassWndHigh(Double_t limit) { fMassWndHigh = TDatabasePDG::Instance()->GetParticle(443)->Mass() + limit ;}
		void SetMassWndLow(Double_t limit) { fMassWndLow = TDatabasePDG::Instance()->GetParticle(443)->Mass() - limit ;}
		void SetCrystalBallFunction(Bool_t okCB) {fCrystalBallParam = okCB;}
		void SetExponentialFunction(Bool_t okExp) {fExponentialParam = okExp;}
 
                void SetWeightType(Double_t wFF, Double_t wFS, Double_t wSS) {fWeightType[0]= wSS; fWeightType[1]= wFS; fWeightType[2]= wFF;}

                void SetChangeResolution(Double_t change){fChangeResolution = change;}   
                void SetChangeMass(Double_t change){fChangeMass = change;}   
		void ComputeMassIntegral(); 
		void ReadMCtemplates(Int_t BinNum);
		void PrintStatus();

		//specific methods for multi-variational fit
		void SetBkgWeights(Double_t ***bkgWgt){for(int kpt=0; kpt<(fPtWindows->GetSize()-1); kpt++){for(int k=0; k<(fMassWindows->GetSize()-2); k++) { for(int ktype=0; ktype<3; ktype++) fWeights[k][kpt][ktype]=bkgWgt[k][kpt][ktype]; }}}
 
		void SetBkgFunction(Int_t massRange,Int_t type, Int_t ptB, TF1 *histBkg){
                fFunBkgSaved[ptB][massRange][type] = histBkg;
                }
                TF1 *GetBkgFunction(Int_t massRange, Int_t ptB, Int_t type) const {return fFunBkgSaved[ptB][massRange][type];} 
                TF1 *GetFunBFunction(Int_t ptB, Int_t type) const {return  fFunBSaved[ptB][type];} 
                void SetFunBFunction(Int_t type, Int_t ptB, TF1 *histSec) { fFunBSaved[ptB][type] = histSec; }
                void SetBackgroundSpecificParameters(Int_t pt, Int_t mb, Int_t tp);
                void SetExtrapolationRegion(Int_t extrRegion){fSignalBinForExtrapolation = extrRegion;}   					   void SetLoadFunction(Bool_t loadFunc) { fLoadFunctions = loadFunc;}
                void SetMultivariateFit(Bool_t multVar) { fMultivariate = multVar;}
                Bool_t GetMultivariate() const { return fMultivariate;} 
                void SetFunctionsSaved(Int_t npxFunB=5000, Int_t npxFunBkg=5000, Double_t funBLimits = 20000., Double_t funBkgLimits = 40000., Int_t signalRegion=2);
                void SetResParams(Double_t ***pars){ fResParams = pars;}
                void SetBkgParams(Float_t ****pars){ fBkgParams = pars;}
                void SetMassWindows(TArrayD *msWnd){ fMassWindows = msWnd;}
                void SetPtWindows(TArrayD *ptWnd){ fPtWindows = ptWnd;}
                void InitializeFunctions(Int_t ptSize, Int_t massSize);
		Double_t EvaluateCDFInvMassSigDistr(Double_t m) const ;
		Double_t EvaluateCDFInvMassBkgDistr(Double_t m) const;
		Double_t ResolutionFunc(Double_t x, Double_t pt, Int_t type) const;
		Double_t EvaluateCDFDecayTimeBkgDistr(Double_t x, Int_t type, Double_t m=3.09, Double_t pt=200.) const ;
	

	private:  
		Double_t fParameters[49];        /*  par[0]  = weightRes;                
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
						     par[18] = norm1Gaus; // resolution param used for First-First
						     par[19] = norm2Gaus;
                                                     par[20] = fMean1ResFunc;
                                                     par[21] = fSigma1ResFunc;
                                                     par[22] = fMean2ResFunc;
                                                     par[23] = fSigma2ResFunc;
                                                     par[24] = fResAlfa;  
                                                     par[25] = fResLambda;
						     par[26] = fResNormExp;
                                                     par[27] = norm1Gaus;    // resolution param used for First-Second
                                                     par[28] = norm2Gaus;
                                                     par[29] = fMean1ResFunc;
                                                     par[30] = fSigma1ResFunc;
                                                     par[31] = fMean2ResFunc;
                                                     par[32] = fSigma2ResFunc;
                                                     par[33] = fResAlfa;  
                                                     par[34] = fResLambda;
                                                     par[35] = fResNormExp;
                                                     par[36] = norm1Gaus;    // resolution param used for Second-Second
                                                     par[37] = norm2Gaus;
                                                     par[38] = fMean1ResFunc;
                                                     par[39] = fSigma1ResFunc;
                                                     par[40] = fMean2ResFunc;
                                                     par[41] = fSigma2ResFunc;
                                                     par[42] = fResAlfa; 
                                                     par[43] = fResLambda;
                                                     par[44] = fResNormExp;   
						     ------------------------------> parameters to describe x and m à la CDF
						     // additional parameters below (added for PbPb analysis): if not used should be fixed
						     // by the FitHandler pointer
                                                     par[45] = fOneOvLamSym1; //  additional parameter for background
						     par[46] = fSym1;         //  additional parameter for background
                                                     par[47] = fPolyn4; //additional parameter for inv. mass background
                                                     par[48] = fPolyn5; //additional parameter for inv. mass background */
						   
 		Double_t fFPlus;                     // parameter of the log-likelihood function
		Double_t fFMinus;                    // slopes of the x distributions of the background 
		Double_t fFSym;                      // functions 

		Double_t fintmMassSig;               // integral of invariant mass distribution for the signal
		Double_t fintmMassBkg;               // integral of invariant mass distribution for the bkg

		TH1F *fhCsiMC;                       // X distribution used as MC template for JPSI from B

                Double_t fShiftTemplate;             // to shift the MC template
		Double_t fMassWndHigh;               // JPSI Mass window higher limit
		Double_t fMassWndLow;                // JPSI Mass window lower limit
		Bool_t fCrystalBallParam;            // Boolean to switch to Crystall Ball parameterisation

                Double_t fWeightType[3];             // vector with weights of candidates types (used to draw functions)            
                Double_t fChangeResolution; 	     // constant to change the RMS of the resolution function
                Double_t fChangeMass;		     // constant to change the RMS of the signal mass function
                
		// data members for multivariate fit
		Double_t ***fWeights;		     // weights to interpolate bkg shape in pt, mass, type bins under signal region
                Bool_t fLoadFunctions;		     // boolean to load functions saved 
                Bool_t fMultivariate;	             // switch-on multivariate fit
                TF1 ***fFunBSaved;		     // pointers to save functions for x of non-prompt J/psi in pt, mass, type bins
                TF1 ****fFunBkgSaved;		     // pointers to save functions for x of background in pt, mass, type bins
                Double_t ***fResParams;		     // resolution function parameters in pt and type bins
                Float_t ****fBkgParams;              // x background parameters in pt, mass, type bins
                TArrayD *fMassWindows;		     // limits for invariant mass bins
                TArrayD *fPtWindows;		     // limits for pt bins
		Bool_t fExponentialParam;            // Boolean to switch to Exponential parameterisation
                Double_t fSignalBinForExtrapolation; // inv. mass region in which extrapolate the background shape
						     // starting from background shapes in the near inv. mass regions


		Double_t EvaluateCDFfunc(Double_t x, Double_t m, Double_t pt, Int_t type) const ;
		Double_t EvaluateCDFfuncNorm(Double_t x, Double_t m, Double_t pt, Int_t type) const ;
		Double_t EvaluateCDFfuncSignalPart(Double_t x, Double_t m, Double_t pt, Int_t type) const ;      // Signal part 
		Double_t EvaluateCDFDecayTimeSigDistr(Double_t x, Double_t pt, Int_t type) const ;
		Double_t EvaluateCDFDecayTimeSigDistrFunc(const Double_t* x, const Double_t *par) const { return par[0]*EvaluateCDFDecayTimeSigDistr(x[0],par[1],(Int_t)par[2]);}
                Double_t EvaluateCDFInvMassSigDistrFunc(const Double_t* x, const Double_t *par) const {return par[0]*EvaluateCDFInvMassSigDistr(x[0])/fintmMassSig;}
		Double_t EvaluateCDFfuncBkgPart(Double_t x,Double_t m, Double_t pt, Int_t type) const ;          // Background part
		Double_t EvaluateCDFDecayTimeBkgDistrFunc(const Double_t* x, const Double_t *par) const { return EvaluateCDFDecayTimeBkgDistr(x[0],(Int_t)par[1],par[2],par[3])*par[0];}
                Double_t EvaluateCDFDecayTimeBkgDistrFuncAllTypes(const Double_t* x, const Double_t *par) const {return (fWeightType[2]*EvaluateCDFDecayTimeBkgDistr(x[0],2)+fWeightType[1]*EvaluateCDFDecayTimeBkgDistr(x[0],1)+fWeightType[0]*EvaluateCDFDecayTimeBkgDistr(x[0],0))*par[0];}
                Double_t EvaluateCDFInvMassBkgDistrFunc(const Double_t* x, const Double_t *par) const {return par[0]*EvaluateCDFInvMassBkgDistr(x[0])/fintmMassBkg;} 
                  
                Double_t EvaluateCDFInvMassTotalDistr(const Double_t* x, const Double_t *par) const;
	        Double_t EvaluateCDFDecayTimeTotalDistr(const Double_t* x, const Double_t *par) const;	
                ////
                Double_t EvaluateCDFDecayTimeTotalDistrAllTypes(const Double_t* x, const Double_t *par) const;

		Double_t FunB(Double_t x, Double_t pt, Int_t type) const;
		Double_t FunBfunc(const Double_t *x, const Double_t *par) const {return FunB(x[0],par[1],(Int_t)par[2])*par[0];}
                Double_t FunBfuncAllTypes(const Double_t *x, const Double_t *par) const {return (fWeightType[2]*FunB(x[0],200.,2)+fWeightType[1]*FunB(x[0],200.,1)+fWeightType[0]*FunB(x[0],200.,0))*par[0];}
                Double_t FunP(Double_t x, Double_t pt,Int_t type) const ;
		Double_t CsiMC(Double_t x) const;
		Double_t CsiMCfunc(const Double_t* x, const Double_t *par) const {  return CsiMC(x[0])*par[0];}
		Double_t FunBkgPos(Double_t x, Double_t pt, Int_t type) const ;
		Double_t FunBkgNeg(Double_t x, Double_t pt, Int_t type) const ;
		Double_t FunBkgSym(Double_t x, Double_t pt, Int_t type) const ;
		Double_t FunBkgSym1(Double_t x, Double_t pt, Int_t type) const ;
		Double_t ResolutionFuncf(const Double_t* x, const Double_t *par) const { return ResolutionFunc(x[0],par[1],(Int_t)par[2])*par[0];}
                Double_t ResolutionFuncAllTypes(const Double_t* x, const Double_t *par) const { return (fWeightType[2]*ResolutionFunc(x[0],200.,2)+fWeightType[1]*ResolutionFunc(x[0],200.,1)+fWeightType[0]*ResolutionFunc(x[0],200.,0))*par[0]; }                 

                // private methods for multivariate fit 
                Double_t EvaluateCDFDecayTimeBkgDistrSaved(Double_t x, Int_t type, Double_t m=3.09, Double_t pt = 200.) const ;
		Double_t EvaluateCDFDecayTimeBkgDistrDifferential(Double_t x, Int_t type, Double_t m=3.09, Double_t pt = 200.) const;
                Double_t FunBsaved(Double_t x, Double_t pt, Int_t type) const;

		ClassDef (AliDielectronBtoJPSItoEleCDFfitFCN,1);         // Unbinned log-likelihood fit 

};

#endif
