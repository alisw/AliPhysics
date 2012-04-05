#ifndef ALIDIELECTRONBTOJPSITOELECDFFITHANDLER_H
#define ALIDIELECTRONBTOJPSITOELECDFFITHANDLER_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//___________________________________________________________________________
//                  Class AliDielectronBtoJPSItoEleCDFfitHandler
//            Class to perform unbinned log-likelihood fit
//      
//                         Origin: C. Di Giglio
//     Contact: Carmelo.Digiglio@ba.infn.it; Giuseppe.Bruno@ba.infn.it
//____________________________________________________________________________

#include <TNamed.h>
#include <TBits.h>
#include <TGraph.h>
#include <TFitter.h>

class AliDielectronBtoJPSItoEleCDFfitFCN ;

class AliDielectronBtoJPSItoEleCDFfitHandler : public TNamed {
	public:
		//
		AliDielectronBtoJPSItoEleCDFfitHandler();
		AliDielectronBtoJPSItoEleCDFfitHandler& operator= (const  AliDielectronBtoJPSItoEleCDFfitHandler& c);
		AliDielectronBtoJPSItoEleCDFfitHandler(const AliDielectronBtoJPSItoEleCDFfitHandler& c);
		AliDielectronBtoJPSItoEleCDFfitHandler(Double_t* decaytime, Double_t* invariantmass, Int_t *type, Int_t ncand);
		~AliDielectronBtoJPSItoEleCDFfitHandler(); 
		Double_t Up() const { return fUp; }
		void SetErrorDef(Double_t up) {fUp = up;}
		void SetPrintStatus(Bool_t okPrintStatus) { fPrintStatus = okPrintStatus; } 
		Bool_t GetPrintStatus() const { return fPrintStatus ; }
		void SetParamStartValues (Double_t*);
		Double_t* GetStartParamValues() { return fParamStartValues; }
		TBits GetFixedParamList() const { return fIsParamFixed; }
		TFitter *GetFitter() const {return fitter;}
                Double_t GetParameter(Int_t numPar) const {return fitter->GetParameter(numPar);}
                Double_t GetParameterError(Int_t numPar) const {return fitter->GetParError(numPar);}   

                void FixParam(UInt_t param, Bool_t value) { fIsParamFixed.SetBitNumber(param,value); }
		void FixAllParam(Bool_t value) { for(UInt_t par=0;par<45;par++) fIsParamFixed.SetBitNumber(par,value); }
		Bool_t IsParamFixed(UInt_t param) { return fIsParamFixed.TestBitNumber(param); }
		void SetResolutionConstants(Double_t* resolutionConst, Int_t type);
		void SetCrystalBallFunction(Bool_t okCB);
		void SetMassWndHigh(Double_t limit);
		void SetMassWndLow(Double_t limit);

		Double_t operator()(const Double_t* par) const ;
		void CdfFCN(Int_t & /* npar */, Double_t * /* gin */, Double_t &f, Double_t *par, Int_t /* iflag */);

		Double_t* Decaytime() const         { return fX; }
		Double_t* InvariantMass() const     { return fM; }
                Int_t*    TypeCand() const          { return fType;}
 		AliDielectronBtoJPSItoEleCDFfitFCN* LikelihoodPointer() const { return fLikely; }
		Int_t DoMinimization(Int_t step = 0);
                

	private:
		//
		TBits fIsParamFixed;                               //array of bits: 0 = param free; 1 = param fixed;
		Bool_t fPrintStatus;                               //flag to enable the prit out of the algorithm at each step
		Double_t fParamStartValues[45];                    //array of parameters input value
		Double_t fUp;                                      //error definition 
		Double_t* fX; 	                     	           //pseudo-proper decay time X
		Double_t* fM;                                      //invariant mass M
                Int_t* fType;                                      //candidate type
		AliDielectronBtoJPSItoEleCDFfitFCN* fLikely;       //Log likelihood function
		Int_t fNcand;                                      //number of candidates
		TGraph* fContPlot1;                                //contour plot       
		TGraph* fContPlot2;                                //contour plot       
		TGraph* fContPlot3;                                //contour plot       
		TFitter *fitter;                                   //pointer to TFitter object	
		//
		ClassDef(AliDielectronBtoJPSItoEleCDFfitHandler,1);

}; 
#endif
