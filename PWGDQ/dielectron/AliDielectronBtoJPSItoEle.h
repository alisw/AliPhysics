#ifndef ALIDIELECTRONBTOJPSITOELE_H
#define ALIDIELECTRONBTOJPSITOELE_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                         Class AliDielectronBtoJPSItoEle
//                  Unbinned log-likelihood fit analysis class
//
//                             Origin: C.Di Giglio
//       Contact: Carmelo.Digiglio@ba.infn.it , giuseppe.bruno@ba.infn.it
//-------------------------------------------------------------------------

#include "TH1F.h"

class TNtuple ;
class AliDielectronBtoJPSItoEleCDFfitHandler ; 
class AliDielectronBtoJPSItoEleCDFfitFCN ; 

class AliDielectronBtoJPSItoEle : public TNamed {
	public:
		//
		AliDielectronBtoJPSItoEle();
		AliDielectronBtoJPSItoEle(const AliDielectronBtoJPSItoEle& source);
		AliDielectronBtoJPSItoEle& operator=(const AliDielectronBtoJPSItoEle& source);
		virtual ~AliDielectronBtoJPSItoEle();

		Int_t DoMinimization(Int_t step = 0);
		void ReadCandidates(TNtuple* nt, Double_t* &x, Double_t* &m, Double_t* &pt, Int_t * &typeCand, Int_t& n,Double_t massLow = -1., Double_t massUp = -1., Double_t ptLow = -1., Double_t ptUp = -1.); // primary JPSI + secondary JPSI + bkg couples

		void SetCsiMC();
		void SetFitHandler(Double_t* x /*pseudoproper*/, Double_t* m /*inv mass*/, Double_t *pt /*transverse momentum */, Int_t *type /*type*/, Int_t ncand /*candidates*/); 
		void CloneMCtemplate(const TH1F* MCtemplate) {fMCtemplate = (TH1F*)MCtemplate->Clone("fMCtemplate");}
		void SetResTypeAnalysis(TString resType){fResType = resType;}
                Double_t* GetResolutionConstants(Double_t* resolutionConst);
		AliDielectronBtoJPSItoEleCDFfitHandler* GetCDFFitHandler() const { return fFCNfunction ; }

	private:
		//
		AliDielectronBtoJPSItoEleCDFfitHandler* fFCNfunction; //! pointer to the interface class
		TH1F* fMCtemplate;			    //! template of the MC distribution for the x distribution of the secondary J/psi
                TString fResType;                           // string with candidate's types considered

		ClassDef(AliDielectronBtoJPSItoEle,1); // AliDielectronBtoJPSItoEle class
};

#endif
