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

		Int_t DoMinimization();
		void ReadCandidates(TNtuple* nt, Double_t* &x, Double_t* &m, Int_t& n); // primary JPSI + secondary JPSI + bkg couples

		void SetPtBin(Int_t BinNum) { fPtBin = BinNum ; }
		void SetCsiMC();
		void SetFitHandler(Double_t* x /*pseudoproper*/, Double_t* m /*inv mass*/, Int_t ncand /*candidates*/); 
		void CloneMCtemplate(const TH1F* MCtemplate) {fMCtemplate = (TH1F*)MCtemplate->Clone("fMCtemplate");}
		Double_t* GetResolutionConstants(Double_t* resolutionConst);
		AliDielectronBtoJPSItoEleCDFfitHandler* GetCDFFitHandler() const { return fFCNfunction ; }
		Int_t GetPtBin() const { return fPtBin ; }

	private:
		//
		AliDielectronBtoJPSItoEleCDFfitHandler* fFCNfunction; //! pointer to the interface class
		Int_t fPtBin;                               // number of pt bin in which the analysis is performes
		TH1F* fMCtemplate;			      //! template of the MC distribution for the x distribution of the secondary J/psi

		ClassDef(AliDielectronBtoJPSItoEle,1); // AliDielectronBtoJPSItoEle class
};

#endif
