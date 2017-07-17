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
#include <TComplex.h>

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

		void SetInputList( TClonesArray *inputarray){fInputList = inputarray;}
		void SetEventCentrality( float cent ){fCent = cent;}
		void SetEventImpactParameter( float ip ){ fImpactParameter = ip; }
		void SetEventVertex( double *vtx ){ fVertex = vtx; }
		void SetIsPhiModule( Bool_t isphi ){ IsPhiModule = isphi ; }
		void SetPhiModuleHistos( int cent, int sub, TH1D *hModuledPhi);

		void SetEtaRange( double eta_min, double eta_max){fEta_min = eta_min; fEta_max = eta_max; }
		void SetDebugLevel( int dblv ){ fDebugLevel = dblv; }
		void SetEffConfig( int Mode, int FilterBit ){ fEffMode = Mode; fEffFilterBit = FilterBit; cout << "fEffMode set = " << fEffMode << endl;}
		void SetIsSCptdep( Bool_t isSCptdep ){ IsSCptdep = isSCptdep; cout << "doing addtional loop to check SC pt dep = "<< IsSCptdep << endl; }
		void SetSCwithQC(Bool_t isSCwithQC){ IsSCwithQC = isSCwithQC; cout << "doing additinal loop for SC results with QC method = " << IsSCwithQC << endl; }
		void SetEbEWeight(Bool_t isEbEWeighted){ IsEbEWeighted = isEbEWeighted; cout << "use event weight = " << IsEbEWeighted << endl; }
		void SetQCEtaCut( Double_t QC_eta_cut_min, Double_t QC_eta_cut_max){
			fQC_eta_cut_min=QC_eta_cut_min; fQC_eta_cut_max=QC_eta_cut_max; cout<<"setting eta range for QC" << fQC_eta_cut_min << "~" << fQC_eta_cut_max << endl;
		}


		void SetEventTracksQA(float tpc, float glb){ fTPCtrks = tpc; fGlbtrks = glb;};

		inline void DEBUG(int level, TString msg){if(level<fDebugLevel) std::cout<<level<<"\t"<<msg<<endl;};

		TComplex CalculateQnSP( double eta1, double eta2, int harmonics);

		double Get_Qn_Real_pt(double eta1, double eta2, int harmonics, int ipt, double pt_min, double pt_max);
		double Get_Qn_Img_pt(double eta1, double eta2, int harmonics, int ipt, double pt_min, double pt_max);
		double Get_QC_Vn( double QnA_real, double QnA_img, double QnB_real, double QnB_img);
		void Fill_QA_plot(double eta1, double eta2 );

		double Get_ScaledMoments( int k, int harmonics);
		AliJEfficiency* GetAliJEfficiency() { return fEfficiency; }

		// new function for QC method //
		void CalculateQvectorsQC();
		TComplex Q(int n, int p);
		TComplex Two( int n1, int n2);
		TComplex Four( int n1, int n2, int n3, int n4);

		// Getter for single vn
		Double_t Get_vn( int ih, int imethod ){ return fSingleVn[ih][imethod]; } // method 0:SP, 1:QC(with eta gap), 2:QC(without eta gap)


	private:
		enum{kH0, kH1, kH2, kH3, kH4, kH5, kH6, kH7, kH8, kNH}; //harmonics // do we need vn up to v8? .. yes we need..
		enum{kK0, kK1, kK2, kK3, kK4, nKL}; // order // do we really need vn^8

//		TDirectory           *fOutput;     // Output
		Long64_t AnaEntry;
		TClonesArray * fInputList;
		double *fVertex;//!
		Float_t		fCent;
		Float_t		fImpactParameter;
		int			fDebugLevel;
		int fNCent;
		int fCBin;
		double *fCentBin;//!
		int fNJacek;
		double *fPttJacek;//!
		AliJEfficiency *fEfficiency;
		int fEffMode;
		int fEffFilterBit;
		float fTPCtrks;
		float fGlbtrks;
		TString fInFileName;
		Bool_t IsPhiModule;
		Bool_t IsSCwithQC; // flag to check SC with QC method
		Bool_t IsSCptdep;  // flag to check SC pt dep or not
		Bool_t IsEbEWeighted; // flag for ebe weight for QC method
		Double_t fSingleVn[kNH][3]; // 3 method

// Histograms
		double fEta_min;
		double fEta_max;
		double NSubTracks[2];

		Double_t fQC_eta_cut_min;
		Double_t fQC_eta_cut_max;

		TComplex QvectorQC[kNH][nKL];
		TComplex QvectorQCeta10[kNH][nKL][2]; // ksub

		TH1D *h_phi_module[7][2]; // cent, isub
		TFile *inclusFile; // pointer for root file

		AliJHistManager * fHMG;//!

		AliJBin fBin_Subset;//!
		AliJBin fBin_h;//!
		AliJBin fBin_k;//!
		AliJBin fBin_hh;//!
		AliJBin fBin_kk;//!
		AliJBin fHistCentBin;//!
		AliJBin fVertexBin;//! // x, y, z
		AliJBin fCorrBin;//!

		AliJTH1D fh_cent;//! // for cent dist
		AliJTH1D fh_ImpactParameter;//! // for impact parameter for mc
		AliJTH1D fh_vertex;//!
		AliJTH1D fh_pt;//! // for pt dist of tracks
		AliJTH1D fh_eta;//! // for eta dist of tracks
		AliJTH1D fh_phi;//! // for phi dist [ic][isub]
		AliJTH1D fh_Qvector;//! // for Q-Vector dist [ic][isub][ih]

		AliJTH1D fh_ntracks;//! // for number of tracks dist
		AliJTH1D fh_vn;//!  // single vn^k  array [ih][ik][iCent]
		AliJTH1D fh_vn_vn;//! // combination for <vn*vn> [ih][ik][ihh][ikk][iCent]

		AliJTH1D fh_correlator;//! // some more complex correlator
		AliJTH2D fh_TrkQA_TPCvsGlob;//! // QA histos
		AliJTH2D fh_TrkQA_TPCvsCent;//! // QA histos

		// addtinal variables for ptbins(Standard Candles only)
		enum{kPt0, kPt1, kPt2, kPt3, kPt4, kPt5, kPt6, kPt7, N_ptbins};
		double NSubTracks_pt[2][N_ptbins];
		AliJBin fBin_Nptbins;//!
		AliJTH1D fh_SC_ptdep_4corr;//! // for < vn^2 vm^2 >
		AliJTH1D fh_SC_ptdep_2corr;//!  // for < vn^2 >
		// additinal variables for SC with QC
		AliJTH1D fh_SC_with_QC_4corr;//! // for <vn^2 vm^2>
		AliJTH1D fh_SC_with_QC_2corr;//! // for <vn^2>
		AliJTH1D fh_SC_with_QC_2corr_eta10;//!
		AliJTH2D fh_QvectorQC;//! // check for Q-vec dist for [ic][ih]
		AliJTH1D fh_QvectorQCphi;//!
		AliJTH1D fh_evt_SP_QC_ratio_2p;//! // check SP QC evt by evt ratio
		AliJTH1D fh_evt_SP_QC_ratio_4p;//! // check SP QC evt by evt ratio
		ClassDef(AliJFFlucAnalysis, 1); // example of analysis
};

#endif
