#ifndef AliJFFlucAnalysis_cxx
#define AliJFFlucAnalysis_cxx

#include <vector>
#include <TVector.h>
#include <TRandom.h>
#include <TString.h>
#include <TPRegexp.h>
#include <TObjArray.h>
#include <TGrid.h>
#include "AliAnalysisTaskSE.h"
#include "AliGenEventHeader.h"
#include "AliJEfficiency.h"
#include "AliJHistManager.h"
#include "AliVVertex.h"

class TClonesArray;
class AliJBaseTrack;
class AliJEfficiency;

class AliJFFlucAnalysis : public AliAnalysisTaskSE {

	public:
		AliJFFlucAnalysis();
		AliJFFlucAnalysis(const char *name);
		AliJFFlucAnalysis(const AliJFFlucAnalysis& a); // not implemented
		AliJFFlucAnalysis& operator=(const AliJFFlucAnalysis& ap); // not implemented

		virtual ~AliJFFlucAnalysis(); 
		virtual void   UserCreateOutputObjects();
		virtual void   Init();
		virtual void   UserExec(Option_t *option);
		virtual void   Terminate(Option_t *);

		void SetInputList( TClonesArray *inputarray){fInputList = inputarray;};
		void SetEventCentrality( float cent ){fCent = cent;};
		void SetEventImpactParameter( float ip ){ fImpactParameter = ip; };
		void SetEventVertex( double *vtx ){ fVertex = vtx; };
		void SetInFileName( TString inName){ fInFileName = inName ; };
		void SetIsPhiModule( Bool_t isphi ){ IsPhiModule = isphi ; } ;

		void SetEtaRange( double eta_min, double eta_max){fEta_min = eta_min; fEta_max = eta_max; };
		void SetDebugLevel( int dblv ){fDebugLevel = dblv;};
		void SetEffConfig( int Mode, int FilterBit ){ fEffMode = Mode; fEffFilterBit = FilterBit; 
									cout << "fEffMode set = " << fEffMode << endl;};

		inline void DEBUG(int level, TString msg){if(level<fDebugLevel) std::cout<<level<<"\t"<<msg<<endl;};


		double Get_Qn_Real(double eta1, double eta2, int harmonics); 
		double Get_Qn_Img(double eta1, double eta2, int harmonics); 
		double Get_QC_Vn( double QnA_real, double QnA_img, double QnB_real, double QnB_img);

		void Fill_QA_plot(double eta1, double eta2 );

		double Complex_product_real( double QnA_real, double QnA_img, double QnB_real, double QnB_img);
		double Complex_product_img( double QnA_real, double QnA_img, double QnB_real, double QnB_img);
		double Complex_abs( double real, double img);
		double Complex_sqr_real( double real, double img);
		double Complex_sqr_img( double real, double img);
		double Get_ScaledMoments( int k, int harmonics);
		AliJEfficiency* GetAliJEfficiency() { return fEfficiency;}


	private:
//		TDirectory           *fOutput;     // Output
		Long64_t AnaEntry; 
		TClonesArray * fInputList;
		double * fVertex;
		Float_t		fCent;
		Float_t		fImpactParameter;
		int			fDebugLevel;
		int fNCent;
		int fCBin;
		double *fCentBin;
		int fNJacek;  
		double *fPttJacek;
		AliJEfficiency *fEfficiency;
		int fEffMode;
		int fEffFilterBit;
		TString fInFileName;
		Bool_t IsPhiModule;  

// Histograms
		enum{kH0, kH1, kH2, kH3, kH4, kH5, kH6, kH7, kNH}; //harmonics
		enum{kK0, kK1, kK2, kK3, kK4, nKL}; // order
		double fEta_min;
		double fEta_max;

		TH1D *h_phi_module[7][2]; // cent, isub 
		TFile *inclusFile; // pointer for root file  

		AliJHistManager * fHMG;
	

		AliJBin fBin_Subset;
		AliJBin fBin_h; 
		AliJBin fBin_k;
		AliJBin fBin_hh;
		AliJBin fBin_kk;
		AliJBin fHistCentBin;
		AliJBin fVertexBin; // x, y, z 
		AliJBin fCorrBin;

		AliJTH1D fh_cent; // for cent dist
		AliJTH1D fh_ImpactParameter; // for impact parameter for mc
		AliJTH1D fh_vertex;
		AliJTH1D fh_pt; // for pt dist of tracks 
		AliJTH1D fh_eta; // for eta dist of tracks
		AliJTH1D fh_phi; // for phi dist [ic][isub]
		AliJTH1D fh_Qvector; // for Q-Vector dist [ic][isub][ih]

		AliJTH1D fh_ntracks; // for number of tracks dist
		AliJTH1D fh_vn;  // single vn^k  array [ih][ik][iCent]
		AliJTH1D fh_vn_vn; // combination for <vn*vn> [ih][ik][ihh][ikk][iCent]

		AliJTH1D fh_correlator; // some more complex correlator


		ClassDef(AliJFFlucAnalysis, 1); // example of analysis
};

#endif
