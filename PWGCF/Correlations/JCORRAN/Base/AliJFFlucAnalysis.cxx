#include <TChain.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TVector.h>
#include <TVector3.h>
#include <TVectorT.h>
#include <TComplex.h>
#include "AliJBaseTrack.h"
#include "AliJFFlucAnalysis.h"
//#include "AliJCorrelations.h"
#include "AliAnalysisManager.h"
//#include "AliAODEvent.h"
//#include "AliAODTrack.h"
//#include "AliVVertex.h"
//#include "AliAODMCParticle.h"
//#include "AliAODMCHeader.h"
//#include "AliMCEventHandler.h"
//#include "AliMCEvent.h"
//#include "AliStack.h"
//#include "AliHeader.h"
//#include "AliGenEventHeader.h"
//#include "AliGenCocktailEventHeader.h"
//#include "AliGenPythiaEventHeader.h"
#include "AliInputEventHandler.h"
//#include "AliESDVertex.h"
//#include "AliVParticle.h"
//#include "AliCentrality.h"
//#include "AliEventplane.h"
//#include "AliJHistManager.h"
#include "TClonesArray.h"
#include "AliJEfficiency.h"
#pragma GCC diagnostic warning "-Wall"

ClassImp(AliJFFlucAnalysis)

		//________________________________________________________________________
AliJFFlucAnalysis::AliJFFlucAnalysis()
	: AliAnalysisTaskSE(),
	fInputList(0),
	//h_phi_module(),
	fEfficiency(0), // pointer to tracking efficiency
	fVertex(0),
	fCent(0),
	fNJacek(0),
	fEffMode(0),
	fEffFilterBit(0),
	fHMG(0),
	fBin_Subset(),
	fBin_h(),
	fBin_k(),
	fBin_hh(),
	fBin_kk(),
	fHistCentBin(),
	fh_cent(),
	fh_ImpactParameter(),
	fh_vertex(),
	fh_eta(),
	fh_phi(),
	//fh_Qvector(),
	fh_ntracks(),
	fh_vn(),
	fh_vna(),
	fh_vn_vn(),
	fh_cn_4c(),
	fh_cn_2c(),
	fh_cn_cn_2c(),
	fh_cn_2c_eta10(),
	fh_cn_cn_2c_eta10()
{
	const int NCent = 7;
	static Double_t CentBin[NCent+1] = {0, 5, 10, 20, 30, 40, 50, 60};
	fNCent = NCent;
	fCentBin = CentBin;


	// pt bins to check pt dist copied from AliJHistos
	const int nJacek = 73;
	static Double_t pttJacek[nJacek+1] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9,
		0.95,1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9,
		10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 45, 50, 60, 70, 80, 90, 100};

	fNJacek = nJacek;
	fPttJacek = pttJacek;

	// Constructor
}
//________________________________________________________________________
AliJFFlucAnalysis::AliJFFlucAnalysis(const char *name)
	: AliAnalysisTaskSE(name),
	fInputList(0),
	//h_phi_module(),
	fEfficiency(0),
	fVertex(0),
	fCent(0),
	fNJacek(0),
	fEffMode(0),
	fEffFilterBit(0),
	fHMG(0),
	fBin_Subset(),
	fBin_h(),
	fBin_k(),
	fBin_hh(),
	fBin_kk(),
	fHistCentBin(),
	fh_cent(),
	fh_ImpactParameter(),
	fh_vertex(),
	fh_eta(),
	fh_phi(),
	//fh_Qvector(),
	fh_ntracks(),
	fh_vn(),
	fh_vna(),
	fh_vn_vn(),
	fh_cn_4c(),
	fh_cn_2c(),
	fh_cn_cn_2c(),
	fh_cn_2c_eta10(),
	fh_cn_cn_2c_eta10()
{
	cout << "analysis task created " << endl;
	const int NCent = 7;
	static Double_t CentBin[NCent+1] = {0, 5, 10, 20, 30, 40, 50, 60};
	fNCent = NCent;
	fCentBin = CentBin;

	fDebugLevel = 0;
	fCent = -1;
	fCBin = -1;
	fEffMode = 0;
	fEffFilterBit =0;
	fInFileName ="";
	IsPhiModule = kFALSE;
	IsSCptdep = kFALSE;
	fEta_min = 0;
	fEta_max = 0;
	fQC_eta_cut_min = -0.8; // default setting
	fQC_eta_cut_max = 0.8; // default setting
	fQC_eta_gap_half = 0.5;
	fImpactParameter = -1;

	for(int icent=0; icent<NCent; icent++){
			for(int isub=0; isub<2; isub++){
					h_phi_module[icent][isub]=NULL;
			}
	}

	// pt bins to check pt dist copied from AliJHistos
	const int nJacek = 73;
	static Double_t pttJacek[nJacek+1] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9,
		0.95,1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9,
		10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 45, 50, 60, 70, 80, 90, 100};
	fNJacek = nJacek;
	fPttJacek = pttJacek;
}

//________________________________________________________________________
AliJFFlucAnalysis::AliJFFlucAnalysis(const AliJFFlucAnalysis& a):
	AliAnalysisTaskSE(a.GetName()),
	fInputList(a.fInputList),
	//h_phi_module(a.h_phi_module),
	fEfficiency(a.fEfficiency),
	fVertex(a.fVertex),
	fCent(a.fCent),
	fNJacek(a.fNJacek),
	fEffMode(a.fEffMode),
	fEffFilterBit(a.fEffFilterBit),
	fHMG(a.fHMG),
	fBin_Subset(a.fBin_Subset),
	fBin_h(a.fBin_h),
	fBin_k(a.fBin_k),
	fBin_hh(a.fBin_hh),
	fBin_kk(a.fBin_kk),
	fHistCentBin(a.fHistCentBin),
	fh_cent(a.fh_cent),
	fh_ImpactParameter(a.fh_ImpactParameter),
	fh_vertex(a.fh_vertex),
	fh_eta(a.fh_eta),
	fh_phi(a.fh_phi),
	//fh_Qvector(a.fh_Qvector),
	fh_ntracks(a.fh_ntracks),
	fh_vn(a.fh_vn),
	fh_vna(a.fh_vna),
	fh_vn_vn(a.fh_vn_vn),
	fh_cn_4c(a.fh_cn_4c),
	fh_cn_2c(a.fh_cn_2c),
	fh_cn_cn_2c(a.fh_cn_cn_2c),
	fh_cn_2c_eta10(a.fh_cn_2c_eta10),
	fh_cn_cn_2c_eta10(a.fh_cn_cn_2c_eta10)
{
	//copy constructor
	//	DefineOutput(1, TList::Class() );
}
//________________________________________________________________________
AliJFFlucAnalysis& AliJFFlucAnalysis::operator = (const AliJFFlucAnalysis& ap){
	// assignment operator
	this->~AliJFFlucAnalysis();
	new(this) AliJFFlucAnalysis(ap);
	return *this;
}
//________________________________________________________________________
void AliJFFlucAnalysis::Init(){
	//
}
//________________________________________________________________________
void AliJFFlucAnalysis::UserCreateOutputObjects(){
	fEfficiency = new AliJEfficiency();
	cout << "********" << endl;
	cout << fEffMode << endl ;
	cout << "********" << endl;
	fEfficiency->SetMode( fEffMode ) ; // 0:NoEff 1:Period 2:RunNum 3:Auto
	fEfficiency->SetDataPath( "alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data" );
	
	fHMG = new AliJHistManager("AliJFFlucHistManager","jfluc");
	// set AliJBin here //
	fBin_Subset .Set("Sub","Sub","Sub:%d", AliJBin::kSingle).SetBin(2);
	fBin_h .Set("NH","NH","NH:%d", AliJBin::kSingle).SetBin(kNH);
	fBin_k .Set("K","K","K:%d", AliJBin::kSingle).SetBin(nKL);

	fBin_hh .Set("NHH","NHH","NHH:%d", AliJBin::kSingle).SetBin(kNH);
	fBin_kk .Set("KK","KK","KK:%d", AliJBin::kSingle).SetBin(nKL);

	fHistCentBin .Set("CentBin","CentBin","Cent:%d",AliJBin::kSingle).SetBin(fNCent);
	fVertexBin .Set("Vtx","Vtx","Vtx:%d", AliJBin::kSingle).SetBin(3);
	fCorrBin .Set("C", "C","C:%d", AliJBin::kSingle).SetBin(28);

	fBin_Nptbins .Set("PtBin","PtBin", "Pt:%d", AliJBin::kSingle).SetBin(N_ptbins);

	// set AliJTH1D here //
	fh_cent
		<< TH1D("h_cent","h_cent", 400, 0, 100)
		<< "END" ;

	fh_ImpactParameter
		<< TH1D("h_IP", "h_IP", 400, -2, 20)
		<< "END" ;

	fh_TrkQA_TPCvsCent
		<< TH2D("h_trk_Cent_vs_TPC","h_trk_Cent_vs_TPC", 100, 0, 100, 100, 0, 3000)
		<< "END" ;

	fh_TrkQA_TPCvsGlob
		<< TH2D("h_trk_Glob_vs_TPC", "h_trk_Glob_vs_TPC", 100, 0, 2000, 100, 0, 3000)
		<< "END";

	fh_TrkQA_FB32_vs_FB32TOF
		<< TH2D("h_trk_FB32_vs_FB32TOF", "h_trk_FB32_vs_FB32TOF", 200, 0, 4000, 100, 0, 2000)
		<< "END";

	fh_vertex
		<< TH1D("h_vertex","h_vertex", 400, -20, 20)
		<< fVertexBin
		<< "END" ;
	fh_pt
		<< TH1D("hChargedPtJacek", "", fNJacek, fPttJacek)
		<< fHistCentBin
		<< "END" ;

	fh_eta
		<< TH1D("h_eta", "h_eta", 100, -15, 15 )
		<< fHistCentBin
		<< "END" ;
	fh_phi
		<< TH1D("h_phi", "h_phi", 100, -10, 10)
		<< fHistCentBin << fBin_Subset
		<< "END" ;
	/*fh_Qvector
		<< TH1D("h_QVector", "h_QVector", 100, -10, 10)
		<< fHistCentBin << fBin_Subset
		<< fBin_h
		<< "END" ;*/

	fh_ntracks
		<< TH1D("h_tracks", "h_tracks", 100, 0, 30000)
		<< fHistCentBin
		<< "END" ;

	fh_vn
		<< TH1D("hvn","hvn", 1024, -1.5, 1.5)
		<< fBin_h << fBin_k
		<< fHistCentBin
		<< "END";   // histogram of vn_h^k values for [ih][ik][iCent]
	fh_vna
		<< TH1D("hvna","hvna", 1024, -1.5, 1.5)
		<< fBin_h << fBin_k
		<< fHistCentBin
		<< "END";   // histogram of vn_h^k values for [ih][ik][iCent]
	fh_vn_vn
		<< TH1D("hvn_vn", "hvn_vn", 1024, -1.5, 1.5)
		<< fBin_h << fBin_k
		<< fBin_hh << fBin_kk
		<< fHistCentBin
		<< "END";  // histo of < vn * vn > for [ih][ik][ihh][ikk][iCent]
	fh_cn_4c
		<< TH1D("hcn_4c","hcn_4c", 1024, -1.5, 1.5)
		<< fBin_h << fBin_k
		<< fHistCentBin
		<< "END";
	fh_cn_2c
		<< TH1D("hcn_2c","hcn_2c", 1024, -1.5, 1.5)
		<< fBin_h << fBin_k
		<< fHistCentBin
		<< "END";
	fh_cn_cn_2c
		<< TH1D("hcn_cn_2c", "hcn_cn_2c", 1024, -1.5, 1.5)
		<< fBin_h << fBin_k
		<< fBin_hh << fBin_kk
		<< fHistCentBin
		<< "END";
	fh_cn_2c_eta10
		<< TH1D("hcn_2c_eta10","hcn_2c_eta10", 1024, -1.5, 1.5)
		<< fBin_h << fBin_k
		<< fHistCentBin
		<< "END";
	fh_cn_cn_2c_eta10
		<< TH1D("hcn_cn_2c_eta10", "hcn_cn_2c_eta10", 1024, -1.5, 1.5)
		<< fBin_h << fBin_k
		<< fBin_hh << fBin_kk
		<< fHistCentBin
		<< "END";
	fh_correlator
		<< TH1D("h_corr", "h_corr", 1024, -1.5, 1.5)
		<< fCorrBin
		<< fHistCentBin
		<< "END" ;

	fh_SC_ptdep_4corr
		<< TH1D("hvnvm_SC","hvnvm_SC", 1024, -1.5, 1.5)
		<< fBin_h << fBin_k
		<< fBin_hh << fBin_kk
		<< fHistCentBin << fBin_Nptbins
		<< "END" ;
	fh_SC_ptdep_2corr
		<< TH1D("hvn_SC","hvn_SC", 1024, -1.5, 1.5)
		<< fBin_h << fBin_k
		<< fHistCentBin << fBin_Nptbins
		<< "END" ;

	fh_SC_with_QC_4corr
		<< TH1D("hQC_SC4p", "hQC_SC4p", 1024, -1.5, 1.5)
		<< fBin_h << fBin_hh
		<< fHistCentBin
		<< "END" ;
	fh_SC_with_QC_2corr
		<< TH1D("hQC_SC2p", "hQC_SC2p", 1024, -1.5, 1.5)
		<< fBin_h
		<< fHistCentBin
		<< "END" ;
	fh_SC_with_QC_2corr_eta10
		<< TH1D("hQC_SC2p_eta10", "hQC_SC2p_eta10", 1024, -1.5, 1.5)
		<< fBin_h
		<< fHistCentBin
		<< "END" ;
	/*fh_QvectorQC
		<< TH2D("hQvecQC", "hQvecQC", 1024, -1.1 , 1.1, 1024, -1.1, 1.1 )
		<< fBin_h
		<< fHistCentBin
		<< "END" ;

	fh_QvectorQCphi
		<< TH1D("hQbecQCphi", "hQbecQCphi", 1024, -3.2 , 3.2 )
		<< fBin_h
		<< fHistCentBin
		<< "END" ;*/
	fh_evt_SP_QC_ratio_4p
		<< TH1D("hSPQCratio4p", "hSPQCratio4p", 1024, -100, 100)
		<< fBin_h
		<< fHistCentBin
		<< "END" ; // fBin_h > not stand for harmonics, just case(32, 42, 52, 53, 43)

	fh_evt_SP_QC_ratio_2p
		<< TH1D("hSPQCratio", "hSPQCratio", 1024, -100, 100)
		<< fBin_h
		<< fHistCentBin
		<< "END" ; // fBin_h > not stand for harmonics, only for v2, v3, v4, v5
	//AliJTH1D set done.

	fHMG->Print();
	fHMG->WriteConfig();

}

//________________________________________________________________________
AliJFFlucAnalysis::~AliJFFlucAnalysis() {
	delete fInputList;
	delete fHMG;
	delete fEfficiency;
}

//________________________________________________________________________
void AliJFFlucAnalysis::UserExec(Option_t *) {
	// find Centrality
	Double_t inputCent = fCent;
	fCBin = -1;
	for(int iCbin=0; iCbin<fNCent; iCbin++){
		if( inputCent > fCentBin[iCbin] && inputCent < fCentBin[iCbin+1] )
			fCBin = iCbin;
	}
	if(fCBin == -1){return;};
	DEBUG(3, "cent bin found" );
	int trk_number = fInputList->GetEntriesFast();
	DEBUG(3, Form("trk number is %d", trk_number) );
	fh_ntracks[fCBin]->Fill( trk_number ) ;
	DEBUG(3, "trk number filled into histo");
	fh_cent->Fill( inputCent ) ;
	DEBUG(3, "filled cent into histo" );
	fh_ImpactParameter->Fill( fImpactParameter);
	DEBUG(3, "impact parameter has been filled" );
	Fill_QA_plot( fEta_min, fEta_max );
	DEBUG(3, "QA Plot filled");

	enum{kSubA, kSubB, kNSub};
	enum{kMin, kMax};
	enum{kReal, kImg, kNPhase};
	Double_t Eta_config[kNSub][2];
	Eta_config[kSubA][kMin] = fEta_min;  // 0.4 min for SubA
	Eta_config[kSubA][kMax] = fEta_max;  // 0.8 max for SubA
	Eta_config[kSubB][kMin] = -fEta_max; // -0.8  min for SubB
	Eta_config[kSubB][kMax] = -fEta_min; // -0.4  max for SubB

	// use complex variable instead of doulbe Qn //
	TComplex QnA[kNH];
	TComplex QnB[kNH];
	TComplex QnA_star[kNH];
	TComplex QnB_star[kNH];

	//--------------- Calculate Qn--------------------
	for(int ih=0; ih<kNH; ih++){
		QnA[ih] = CalculateQnSP( Eta_config[kSubA][kMin], Eta_config[kSubA][kMax], ih);
		QnB[ih] = CalculateQnSP( Eta_config[kSubB][kMin], Eta_config[kSubB][kMax], ih);
		QnA_star[ih] = TComplex::Conjugate ( QnA[ih] ) ;
		QnB_star[ih] = TComplex::Conjugate ( QnB[ih] ) ;
	}
	NSubTracks[kSubA] = QnA[0].Re(); // this is number of tracks in Sub A
	NSubTracks[kSubB] = QnB[0].Re(); // this is number of tracks in Sub B
	
	// v2^2 :  k=1  /// remember QnQn = vn^(2k) not k
	// use k=0 for check v2, v3 only
	Double_t vn2[kNH][nKL];
	Double_t vn2_vn2[kNH][nKL][kNH][nKL];

	TComplex corr[kNH][nKL];
	TComplex ncorr[kNH][nKL];

	const TComplex *pQn[][2] = {
		{QnA,QnB_star},
		{QnB,QnA_star}
	};
	const double N[][2] = {
		{NSubTracks[0],NSubTracks[1]},
		{NSubTracks[1],NSubTracks[0]}
	};

	for(int i = 0; i < 2; ++i){
		Double_t ebe_2p_weight = 1.0;
		Double_t ebe_3p_weight = 1.0;
		Double_t ebe_4p_weightB = 1.0;
		if( IsEbEWeighted == kTRUE ){
			ebe_2p_weight = N[i][0]*N[i][1];//NSubTracks[kSubA] * NSubTracks[kSubB] ;
			ebe_3p_weight = ebe_2p_weight*(N[i][1]-1.0);// * (NSubTracks[kSubB]-1.0);
			ebe_4p_weightB = ebe_3p_weight*(N[i][1]-2.0);// * (NSubTracks[kSubB]-2.0);
		}
		Double_t ebe_2Np_weight[2*nKL] = {
			ebe_2p_weight,
		};
		if( IsEbEWeighted == kTRUE ){
			for(int ik=1; ik<2*nKL; ik++){
				double dk = (double)ik;
				ebe_2Np_weight[ik] = ebe_2Np_weight[ik-1]*max(N[i][0]-dk,1.0)*max(N[i][1]-dk,1.0);
			}
		}else for(int ik=1; ik<2*nKL; ik++)
			ebe_2Np_weight[ik] = 1.0;

		double mf = 1.0/((N[i][0]-1.0)*(N[i][1]-1.0));

		for(int ih=2; ih<kNH; ih++){
			corr[ih][1] = pQn[i][0][ih]*pQn[i][1][ih];//QnA[ih]*QnB_star[ih];
			for(int ik=2; ik<nKL; ik++)
				corr[ih][ik] = corr[ih][ik-1]*corr[ih][1];//TComplex::Power(corr[ih][1],ik);
			ncorr[ih][1] = corr[ih][1];
			ncorr[ih][2] = mf*(corr[ih][2]*N[i][0]*N[i][1]-pQn[i][1][2*ih]*pQn[i][0][ih]*pQn[i][0][ih]*N[i][0]-pQn[i][0][2*ih]*pQn[i][1][ih]*pQn[i][1][ih]*N[i][1]+pQn[i][1][2*ih]*pQn[i][0][2*ih]);
			for(int ik=3; ik<nKL; ik++)
				ncorr[ih][ik] = corr[ih][ik]; //for 6,8,...-particle correlations, ignore the autocorrelation for now
		}
		
		for(int ih=2; ih<kNH; ih++){
			for(int ik=1; ik<nKL; ik++){ // 2k(0) =1, 2k(1) =2, 2k(2)=4....
				vn2[ih][ik] = corr[ih][ik].Re();
				fh_vn[ih][ik][fCBin]->Fill( vn2[ih][ik] , ebe_2Np_weight[ik-1]);
				fh_vna[ih][ik][fCBin]->Fill(ncorr[ih][ik].Re(), ebe_2Np_weight[ik-1]);
				for( int ihh=2; ihh<kNH; ihh++){
					for(int ikk=1; ikk<nKL; ikk++){
						vn2_vn2[ih][ik][ihh][ikk] = (corr[ih][ik]*corr[ihh][ikk]).Re();
						fh_vn_vn[ih][ik][ihh][ikk][fCBin]->Fill( vn2_vn2[ih][ik][ihh][ikk], ebe_2Np_weight[ik+ikk-1]) ; // Fill hvn_vn
					}
				}
			}
			fSingleVn[ih][0] = TMath::Sqrt(vn2[ih][1]); // fill single vn with SP as method 0
		}

		//************************************************************************

		TComplex V4V2star_2 = pQn[i][0][4] * pQn[i][1][2] * pQn[i][1][2];
		TComplex V4V2starv2_2 =	V4V2star_2 * corr[2][1];//vn[2][1]
		TComplex V4V2starv2_4 = V4V2star_2 * corr[2][2];//vn2[2][2]
		TComplex V5V2starV3starv2_2 = pQn[i][0][5] * pQn[i][1][2] * pQn[i][1][3] * corr[2][1]; //vn2[2][1]
		TComplex V5V2starV3star = pQn[i][0][5] * pQn[i][1][2] * pQn[i][1][3] ;
		TComplex V5V2starV3startv3_2 = V5V2starV3star * corr[3][1]; //vn2[3][1]
		TComplex V6V2star_3 = pQn[i][0][6] * pQn[i][1][2] * pQn[i][1][2] * pQn[i][1][2];
		TComplex V6V3star_2 = pQn[i][0][6] * pQn[i][1][3] * pQn[i][1][3];
		TComplex V6V2starV4star = pQn[i][0][6] * pQn[i][1][2] * pQn[i][1][4];
		TComplex V7V2star_2V3star = pQn[i][0][7] * pQn[i][1][2] * pQn[i][1][2] * pQn[i][1][3];
		TComplex V7V2starV5star = pQn[i][0][7] * pQn[i][1][2] * pQn[i][1][5];
		TComplex V7V3starV4star = pQn[i][0][7] * pQn[i][1][3] * pQn[i][1][4];
		TComplex V8V2starV3star_2 = pQn[i][0][8] * pQn[i][1][2] * pQn[i][1][3] * pQn[i][1][3];
		TComplex V8V2star_4 = pQn[i][0][8] * TComplex::Power(pQn[i][1][2],4);

		// New correlators (Modified by You's correction term for self-correlations)
		double nf = 1.0/(N[i][1]-1.0);
		double ef = nf/(N[i][1]-2.0);
		TComplex nV4V2star_2 = nf*( V4V2star_2*N[i][1] - pQn[i][0][4]*pQn[i][1][4] );
		TComplex nV5V2starV3star = nf*( V5V2starV3star*N[i][1] - pQn[i][0][5]*pQn[i][1][5] );
		TComplex nV6V2star_3 = pQn[i][0][6]*ef*( pQn[i][1][2]*pQn[i][1][2]*pQn[i][1][2]*N[i][1]*N[i][1] - 3.0*pQn[i][1][2]*pQn[i][1][4]*N[i][1] + 2.0*pQn[i][1][6] );
		TComplex nV6V3star_2 = nf*(V6V3star_2*N[i][1] - pQn[i][0][6]*pQn[i][1][6]);
		TComplex nV6V2starV4star = nf*(V6V2starV4star*N[i][1] - pQn[i][0][6]*pQn[i][1][6]);
		TComplex nV7V2star_2V3star = pQn[i][0][7]*ef*( pQn[i][1][2]*pQn[i][1][2]*pQn[i][1][3]*N[i][1]*N[i][1] - 2.0*pQn[i][1][2]*pQn[i][1][5]*N[i][1] - pQn[i][1][3]*pQn[i][1][4]*N[i][1] + 2.0*pQn[i][1][7] );
		TComplex nV7V2starV5star = nf*(V7V2starV5star*N[i][1] - pQn[i][0][7]*pQn[i][1][7]);
		TComplex nV7V3starV4star = nf*(V7V3starV4star*N[i][1] - pQn[i][0][7]*pQn[i][1][7]);
		TComplex nV8V2starV3star_2 = pQn[i][0][8]*ef*( pQn[i][1][2]*pQn[i][1][3]*pQn[i][1][3]*N[i][1]*N[i][1] - 2.0*pQn[i][1][3]*pQn[i][1][5]*N[i][1] - pQn[i][1][2]*pQn[i][1][6]*N[i][1] + 2.0*pQn[i][1][8] );

		TComplex nV4V4V2V2 = (pQn[i][0][4]*pQn[i][1][4]*pQn[i][0][2]*pQn[i][1][2]) - ((1/(N[i][1]-1) * pQn[i][1][6] * pQn[i][0][4] *pQn[i][0][2] ))
			- ((1/(N[i][0]-1) * pQn[i][0][6]*pQn[i][1][4] * pQn[i][1][2])) + (1/((N[i][0]-1)*(N[i][1]-1))*pQn[i][0][6]*pQn[i][1][6] );
		TComplex nV3V3V2V2 = (pQn[i][0][3]*pQn[i][1][3]*pQn[i][0][2]*pQn[i][1][2]) - ((1/(N[i][1]-1) * pQn[i][1][5] * pQn[i][0][3] *pQn[i][0][2] ))
			- ((1/(N[i][0]-1) * pQn[i][0][5]*pQn[i][1][3] * pQn[i][1][2])) + (1/((N[i][0]-1)*(N[i][1]-1))*pQn[i][0][5]*pQn[i][1][5] );
		TComplex nV5V5V2V2 = (pQn[i][0][5]*pQn[i][1][5]*pQn[i][0][2]*pQn[i][1][2]) - ((1/(N[i][1]-1) * pQn[i][1][7] * pQn[i][0][5] *pQn[i][0][2] ))
			- ((1/(N[i][0]-1) * pQn[i][0][7]*pQn[i][1][5] * pQn[i][1][2])) + (1/((N[i][0]-1)*(N[i][1]-1))*pQn[i][0][7]*pQn[i][1][7] );
		TComplex nV5V5V3V3 = (pQn[i][0][5]*pQn[i][1][5]*pQn[i][0][3]*pQn[i][1][3]) - ((1/(N[i][1]-1) * pQn[i][1][8] * pQn[i][0][5] *pQn[i][0][3] ))
			- ((1/(N[i][0]-1) * pQn[i][0][8]*pQn[i][1][5] * pQn[i][1][3])) + (1/((N[i][0]-1)*(N[i][1]-1))*pQn[i][0][8]*pQn[i][1][8] );
		TComplex nV4V4V3V3 = (pQn[i][0][4]*pQn[i][1][4]*pQn[i][0][3]*pQn[i][1][3]) - ((1/(N[i][1]-1) * pQn[i][1][7] * pQn[i][0][4] *pQn[i][0][3] ))
			- ((1/(N[i][0]-1) * pQn[i][0][7]*pQn[i][1][4] * pQn[i][1][3])) + (1/((N[i][0]-1)*(N[i][1]-1))*pQn[i][0][7]*pQn[i][1][7] );

		fh_correlator[0][fCBin]->Fill( V4V2starv2_2.Re() );
		fh_correlator[1][fCBin]->Fill( V4V2starv2_4.Re() );
		fh_correlator[2][fCBin]->Fill( V4V2star_2.Re(),ebe_3p_weight ) ; // added 2015.3.18
		fh_correlator[3][fCBin]->Fill( V5V2starV3starv2_2.Re() );
		fh_correlator[4][fCBin]->Fill( V5V2starV3star.Re(),ebe_3p_weight );
		fh_correlator[5][fCBin]->Fill( V5V2starV3startv3_2.Re() );
		fh_correlator[6][fCBin]->Fill( V6V2star_3.Re(),ebe_4p_weightB );
		fh_correlator[7][fCBin]->Fill( V6V3star_2.Re(),ebe_3p_weight );
		fh_correlator[8][fCBin]->Fill( V7V2star_2V3star.Re(),ebe_4p_weightB ) ;

		fh_correlator[9][fCBin]->Fill( nV4V2star_2.Re(),ebe_3p_weight ); // added 2015.6.10
		fh_correlator[10][fCBin]->Fill( nV5V2starV3star.Re(),ebe_3p_weight );
		fh_correlator[11][fCBin]->Fill( nV6V3star_2.Re(),ebe_3p_weight ) ;

		// use this to avoid self-correlation 4p correlation (2 particles from A, 2 particles from B) -> MA(MA-1)MB(MB-1) : evt weight..
		fh_correlator[12][fCBin]->Fill( nV4V4V2V2.Re(),ebe_2Np_weight[1]);
		fh_correlator[13][fCBin]->Fill( nV3V3V2V2.Re(),ebe_2Np_weight[1]);

		fh_correlator[14][fCBin]->Fill( nV5V5V2V2.Re(),ebe_2Np_weight[1]);
		fh_correlator[15][fCBin]->Fill( nV5V5V3V3.Re(),ebe_2Np_weight[1]);
		fh_correlator[16][fCBin]->Fill( nV4V4V3V3.Re(),ebe_2Np_weight[1]);

		//higher order correlators, added 2017.8.10
		fh_correlator[17][fCBin]->Fill( V8V2starV3star_2.Re(),ebe_4p_weightB );
		fh_correlator[18][fCBin]->Fill( V8V2star_4.Re() ); //5p weight
		fh_correlator[19][fCBin]->Fill( nV6V2star_3.Re(),ebe_4p_weightB );
		fh_correlator[20][fCBin]->Fill( nV7V2star_2V3star.Re(),ebe_4p_weightB );
		fh_correlator[21][fCBin]->Fill( nV8V2starV3star_2.Re(),ebe_4p_weightB );

		fh_correlator[22][fCBin]->Fill( V6V2starV4star.Re(),ebe_3p_weight );
		fh_correlator[23][fCBin]->Fill( V7V2starV5star.Re(),ebe_3p_weight );
		fh_correlator[24][fCBin]->Fill( V7V3starV4star.Re(),ebe_3p_weight );
		fh_correlator[25][fCBin]->Fill( nV6V2starV4star.Re(),ebe_3p_weight );
		fh_correlator[26][fCBin]->Fill( nV7V2starV5star.Re(),ebe_3p_weight );
		fh_correlator[27][fCBin]->Fill( nV7V3starV4star.Re(),ebe_3p_weight );
	}

	CalculateQvectorsQC();

	//cumulants (no mixed harmonics)
	TComplex four[kNH];
	TComplex two[kNH];
	TComplex two_eta10[kNH];

	TComplex M = Q(0,1);
	Double_t qcn4 = (M*(M-TComplex(1,0))*(M-TComplex(2,0))*(M-TComplex(3,0)));
	Double_t qcn = (M*(M-TComplex(1,0))).Re();
	Double_t qcn_10 = (QvectorQCeta10[0][kSubA]*QvectorQCeta10[0][kSubB]).Re();
	Double_t qw1_4 = 1.0, qw1 = 1.0, qw1_10 = 1.0, qw2_10 = 1.0;
	if(IsEbEWeighted == kTRUE){
		qw1_4 = qcn4;
		qw1 = qcn;
		qw1_10 = qcn_10;
		qw2_10 = qcn_10*((QvectorQCeta10[0][kSubA]-TComplex(1,0))*(QvectorQCeta10[0][kSubB]-TComplex(1,0))).Re();
	}

	TComplex corr10[kNH][nKL];

	for(int ih=2; ih < kNH; ih++){
		four[ih] = ((Q(ih,1)*Q(ih,1)*Q(-ih,1)*Q(-ih,1)+Q(2*ih,1)*Q(-2*ih,1)-TComplex(2,0)*(Q(2*ih,1)*Q(-ih,1)*Q(-ih,1)).Re())
			-2.0*(2.0*(M-TComplex(2,0))*(Q(ih,1)*Q(-ih,1))-M*(M-TComplex(3,0))))/qcn4;
		two[ih] = (Q(ih,1)*Q(-ih,1)-M)/(M*(M-TComplex(1,0)));
		two_eta10[ih] = (QvectorQCeta10[ih][kSubA]*TComplex::Conjugate(QvectorQCeta10[ih][kSubB])) / qcn_10;

		corr[ih][1] = two[ih];
		corr10[ih][1] = two_eta10[ih];
		for(int ik=2; ik < nKL; ik++){
			corr[ih][ik] = TComplex::Power(two[ih],ik);
			corr10[ih][ik] = TComplex::Power(two_eta10[ih],ik);
		}
	}

	for(int ih=2; ih < kNH; ih++){
		for(int ik=1; ik<nKL; ik++){
			Double_t cn = TComplex::Power(four[ih],ik).Re();
			fh_cn_4c[ih][ik][fCBin]->Fill(cn,qw1_4);
			fh_cn_2c[ih][ik][fCBin]->Fill(corr[ih][ik].Re(),qw1);
			fh_cn_2c_eta10[ih][ik][fCBin]->Fill(corr10[ih][ik].Re(),qw1_10);

			for( int ihh=2; ihh<kNH; ihh++){
				for(int ikk=1; ikk<nKL; ikk++){
					Double_t cn_cn = (corr[ih][ik]*corr[ihh][ikk]).Re();//(TComplex::Power(two[ih],ik)*TComplex::Power(two[ihh],ikk)).Re();
					fh_cn_cn_2c[ih][ik][ihh][ikk][fCBin]->Fill(cn_cn,qw1_4);
					cn_cn = (corr10[ih][ik]*corr10[ihh][ikk]).Re();//(TComplex::Power(two[ih],ik)*TComplex::Power(two[ihh],ikk)).Re();//(TComplex::Power(two_eta10[ih],ik)*TComplex::Power(two_eta10[ihh],ikk)).Re();
					fh_cn_cn_2c_eta10[ih][ik][ihh][ikk][fCBin]->Fill(cn_cn,qw2_10);
				}
			}
		}
	}

	if(IsSCptdep == kTRUE){
		const int SCNH = 9; // 0, 1, 2(v2), 3(v3), 4(v4), 5(v5)
		Double_t ptbin_borders[N_ptbins+1] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 5.0};
		//init
		TComplex QnA_pt[SCNH][N_ptbins];
		TComplex QnB_pt[SCNH][N_ptbins];
		TComplex QnB_pt_star[SCNH][N_ptbins];
		for(int ih=2; ih<SCNH; ih++){
			for(int ipt=0; ipt<N_ptbins; ipt++){
				QnA_pt[ih][ipt] = TComplex(0,0);
				QnB_pt[ih][ipt] = TComplex(0,0);
				QnB_pt_star[ih][ipt] = TComplex(0,0);
			}
		}

		// calculate Qn for each pt bins
		for(int ih=2; ih<SCNH; ih++){
			for(int ipt=0; ipt<N_ptbins; ipt++){
				Double_t pt_bin_min = ptbin_borders[ipt];
				Double_t pt_bin_max = ptbin_borders[ipt+1];
				Double_t QAReal=Get_Qn_Real_pt( Eta_config[kSubA][0], Eta_config[kSubA][1], ih, ipt, pt_bin_min, pt_bin_max);
				Double_t QAImag=Get_Qn_Img_pt(  Eta_config[kSubA][0], Eta_config[kSubA][1], ih, ipt, pt_bin_min, pt_bin_max);

				Double_t QBReal=Get_Qn_Real_pt( Eta_config[kSubB][0], Eta_config[kSubB][1], ih, ipt, pt_bin_min, pt_bin_max);
				Double_t QBImag=Get_Qn_Img_pt(  Eta_config[kSubB][0], Eta_config[kSubB][1], ih, ipt, pt_bin_min, pt_bin_max);

				QnA_pt[ih][ipt]= TComplex(QAReal, QAImag);
				QnB_pt[ih][ipt]= TComplex(QBReal, QBImag);
				QnB_pt_star[ih][ipt] = TComplex::Conjugate( QnB_pt[ih][ipt] ) ;
			}
		}

		for(int ipt=0; ipt<N_ptbins; ipt++){
			for(int ih=2; ih<SCNH; ih++){
				int ik=1; // v2^2 only (k=1 means ^2)
				fh_SC_ptdep_2corr[ih][ik][fCBin][ipt]->Fill( ( QnA_pt[ih][ipt]*QnB_pt_star[ih][ipt]).Re()) ;
			}
		}


		for(int ipt=0; ipt<N_ptbins; ipt++){
			TComplex nV4V4V2V2_pt = (QnA_pt[4][ipt]*QnB_pt_star[4][ipt]*QnA_pt[2][ipt]*QnB_pt_star[2][ipt]) - ((1/(NSubTracks_pt[1][ipt]-1) * QnB_pt_star[6][ipt] * QnA_pt[4][ipt] *QnA_pt[2][ipt] ))
				- ((1/(NSubTracks_pt[0][ipt]-1) * QnA_pt[6][ipt]*QnB_pt_star[4][ipt] * QnB_pt_star[2][ipt])) + (1/((NSubTracks_pt[0][ipt]-1)*(NSubTracks_pt[1][ipt]-1))*QnA_pt[6][ipt]*QnB_pt_star[6][ipt] );
			TComplex nV3V3V2V2_pt = (QnA_pt[3][ipt]*QnB_pt_star[3][ipt]*QnA_pt[2][ipt]*QnB_pt_star[2][ipt]) - ((1/(NSubTracks_pt[1][ipt]-1) * QnB_pt_star[5][ipt] * QnA_pt[3][ipt] *QnA_pt[2][ipt] ))
				- ((1/(NSubTracks_pt[0][ipt]-1) * QnA_pt[5][ipt]*QnB_pt_star[3][ipt] * QnB_pt_star[2][ipt])) + (1/((NSubTracks_pt[0][ipt]-1)*(NSubTracks_pt[1][ipt]-1))*QnA_pt[5][ipt]*QnB_pt_star[5][ipt] );

			//add higher order 5225, 5335, 4334
			TComplex nV5V5V3V3_pt =  (QnA_pt[5][ipt]*QnB_pt_star[5][ipt]*QnA_pt[3][ipt]*QnB_pt_star[3][ipt]) - ((1/(NSubTracks_pt[1][ipt]-1) * QnB_pt_star[8][ipt] * QnA_pt[5][ipt] *QnA_pt[3][ipt] ))
				- ((1/(NSubTracks_pt[0][ipt]-1) * QnA_pt[8][ipt]*QnB_pt_star[5][ipt] * QnB_pt_star[3][ipt])) + (1/((NSubTracks_pt[0][ipt]-1)*(NSubTracks_pt[1][ipt]-1))*QnA_pt[8][ipt]*QnB_pt_star[8][ipt] );
			TComplex nV5V5V2V2_pt =  (QnA_pt[5][ipt]*QnB_pt_star[5][ipt]*QnA_pt[2][ipt]*QnB_pt_star[2][ipt]) - ((1/(NSubTracks_pt[1][ipt]-1) * QnB_pt_star[7][ipt] * QnA_pt[5][ipt] *QnA_pt[2][ipt] ))
				- ((1/(NSubTracks_pt[0][ipt]-1) * QnA_pt[7][ipt]*QnB_pt_star[5][ipt] * QnB_pt_star[2][ipt])) + (1/((NSubTracks_pt[0][ipt]-1)*(NSubTracks_pt[1][ipt]-1))*QnA_pt[7][ipt]*QnB_pt_star[7][ipt] );
			TComplex nV4V4V3V3_pt =  (QnA_pt[4][ipt]*QnB_pt_star[4][ipt]*QnA_pt[3][ipt]*QnB_pt_star[3][ipt]) - ((1/(NSubTracks_pt[1][ipt]-1) * QnB_pt_star[7][ipt] * QnA_pt[4][ipt] *QnA_pt[3][ipt] ))
				- ((1/(NSubTracks_pt[0][ipt]-1) * QnA_pt[7][ipt]*QnB_pt_star[4][ipt] * QnB_pt_star[3][ipt])) + (1/((NSubTracks_pt[0][ipt]-1)*(NSubTracks_pt[1][ipt]-1))*QnA_pt[7][ipt]*QnB_pt_star[7][ipt] );

			fh_SC_ptdep_4corr[2][1][4][1][fCBin][ipt]->Fill( nV4V4V2V2_pt.Re());
			fh_SC_ptdep_4corr[2][1][3][1][fCBin][ipt]->Fill( nV3V3V2V2_pt.Re());
			fh_SC_ptdep_4corr[2][1][5][1][fCBin][ipt]->Fill( nV5V5V2V2_pt.Re() );
			fh_SC_ptdep_4corr[3][1][4][1][fCBin][ipt]->Fill( nV4V4V3V3_pt.Re() );
			fh_SC_ptdep_4corr[3][1][5][1][fCBin][ipt]->Fill( nV5V5V3V3_pt.Re() ) ;
		}
	}//pt dep done

	if(IsSCwithQC==kTRUE){
		//cumulants (with mixed harmonics)
		Double_t QC_4p_value[kNH][kNH];
		Double_t QC_2p_value[kNH];

		Double_t event_weight_four = 1.0;
		Double_t event_weight_two = 1.0;
		Double_t event_weight_two_eta10 = 1.0;
		if(IsEbEWeighted == kTRUE){
			event_weight_four = Four(0,0,0,0).Re();
			event_weight_two = Two(0,0).Re();
			event_weight_two_eta10 = (QvectorQCeta10[0][kSubA]*QvectorQCeta10[0][kSubB]).Re();
		}

		for(int ih=2; ih < kNH; ih++){
			for(int ihh=2; ihh<ih; ihh++){
				TComplex scfour = Four( ih, ihh, -ih, -ihh ) / Four(0,0,0,0).Re();
				
				fh_SC_with_QC_4corr[ih][ihh][fCBin]->Fill( scfour.Re(), event_weight_four );
				QC_4p_value[ih][ihh] = scfour.Re();
			}

			// Finally we want 2p corr as like
			// 1/( M*(M-1) ) * [ QQ* - M ]
			// two(2,2) = Q2 Q2* - Q0 = Q2Q2* - M
			// two(0,0) = Q0 Q0* - Q0 = M^2 - M
			//two[ih] = Two(ih, -ih) / Two(0,0).Re();
			TComplex sctwo = Two(ih, -ih) / Two(0,0).Re();
			fh_SC_with_QC_2corr[ih][fCBin]->Fill( sctwo.Re(), event_weight_two );
			QC_2p_value[ih] = sctwo.Re();
			// fill single vn  with QC without EtaGap as method 2
			fSingleVn[ih][2] = TMath::Sqrt(sctwo.Re());
			
			TComplex sctwo10 = (QvectorQCeta10[ih][kSubA]*TComplex::Conjugate(QvectorQCeta10[ih][kSubB])) / (QvectorQCeta10[0][kSubA]*QvectorQCeta10[0][kSubB]).Re();
			fh_SC_with_QC_2corr_eta10[ih][fCBin]->Fill( sctwo10.Re(), event_weight_two_eta10 );
			// fill single vn with QC method with Eta Gap as method 1
			fSingleVn[ih][1] = TMath::Sqrt(sctwo10.Re());
		}
		
		//Check evt-by-evt SP/QC ratio. (term-by-term)
		// calculate  (vn^2 vm^2)_SP /  (vn^2 vm^2)_QC
		// 4p ( v3v3v2v2, v4v4v2v2, v5v5v2v2, v5v5v3v3, v4v4v3v3
#if 0
		Double_t SP_4p_value[5] = { nV3V3V2V2.Re(), nV4V4V2V2.Re(), nV5V5V2V2.Re(), nV5V5V3V3.Re(), nV4V4V3V3.Re()};
		Double_t evtSP_QC_ratio_2p = -99.;
		Double_t evtSP_QC_ratio_4p = -99.;
		int har1[5] = {3, 4, 5, 5, 4 }; // m of SC(m,n)
		int har2[5] = {2, 2, 3, 2, 3 }; // n of SC(m,n)
		for(int i=0; i<5; i++){ // i array index (for m, n)
			evtSP_QC_ratio_4p = SP_4p_value[i] / QC_4p_value[har1[i]][har2[i]] ;
			if( evtSP_QC_ratio_4p < -1 || evtSP_QC_ratio_4p > 5.)
				evtSP_QC_ratio_4p = -99;
			fh_evt_SP_QC_ratio_4p[i][fCBin]->Fill( evtSP_QC_ratio_4p );
			// fh_evt_SP_QC_ratio_4p[ ih ][fCBin] : ih is not harmonics in this histo. ( SC(m,n) case)
		}
		// 2p , v2, v3, v4, v5
		for(int i=0; i<4; i++){
			Double_t SP_2p_value = vn2[2+i][1];
			evtSP_QC_ratio_2p = SP_2p_value / QC_2p_value[i+2];
			if( evtSP_QC_ratio_2p < -1 || evtSP_QC_ratio_2p > 5.)
				evtSP_QC_ratio_2p = -99;
			fh_evt_SP_QC_ratio_2p[i][fCBin]->Fill(evtSP_QC_ratio_2p );
		}
#endif

	} // QC method done.

	//1 evt is done...
}

//________________________________________________________________________
void AliJFFlucAnalysis::Terminate(Option_t *)
{
	//	fList = dynamic_cast<TList*> (GetOutputData(1));
	//	if(!fList) { Printf("ERROR: fList not availabe"); return;};
	//
	inclusFile->Close();
	//
	cout<<"Sucessfully Finished"<<endl;
}
//________________________________________________________________________
//________________________________________________________________________
void AliJFFlucAnalysis::Fill_QA_plot( Double_t eta1, Double_t eta2 )
{
	Long64_t ntracks = fInputList->GetEntriesFast();
	for( Long64_t it=0; it< ntracks; it++){
		AliJBaseTrack *itrack = (AliJBaseTrack*)fInputList->At(it); // load track
		Double_t pt = itrack->Pt();
		Double_t effCorr = fEfficiency->GetCorrection( pt, fEffFilterBit, fCent);
		Double_t eta = itrack->Eta();
		int isub = -1;
		if( eta < 0 )
			isub = 0;
		if( eta > 0 )
			isub = 1;
		Double_t phi = itrack->Phi();
		Double_t phi_module_corr =1;
		if( IsPhiModule == kTRUE){
			phi_module_corr = h_phi_module[fCBin][isub]->GetBinContent( (h_phi_module[fCBin][isub]->GetXaxis()->FindBin( phi )) );
		}
		//
		if( TMath::Abs(eta) > eta1 && TMath::Abs(eta) < eta2 ){
			fh_eta[fCBin]->Fill(eta , 1./ effCorr );
			fh_pt[fCBin]->Fill(pt, 1./ effCorr );
			if( eta < 0 )
				fh_phi[fCBin][0]->Fill( phi_module_corr * phi, 1./effCorr) ;
			if( eta > 0 )
				fh_phi[fCBin][1]->Fill( phi_module_corr * phi, 1./effCorr) ;
		}
	}
	for(int iaxis=0; iaxis<3; iaxis++)
		fh_vertex[iaxis]->Fill( fVertex[iaxis] );

	fh_TrkQA_TPCvsCent->Fill(fCent,fTPCtrks);
	fh_TrkQA_TPCvsGlob->Fill(fGlbtrks,fTPCtrks);
	fh_TrkQA_FB32_vs_FB32TOF->Fill(fFB32trks,fFB32TOFtrks);
}

//________________________________________________________________________
TComplex AliJFFlucAnalysis::CalculateQnSP( Double_t eta1, Double_t eta2, int harmonics)
{
	int ih=harmonics;
	TComplex Qn = TComplex(0,0);
	Double_t Sub_Ntrk = 0; // number of Tracks * effCorr * phi modulation factor
	Long64_t ntracks = fInputList->GetEntriesFast();
	for(Long64_t it=0; it< ntracks; it++){
		AliJBaseTrack *itrack = (AliJBaseTrack*)fInputList->At(it); // load track
		Double_t pt = itrack->Pt();
		Double_t eta = itrack->Eta();
		Double_t phi = itrack->Phi();
		if( eta < eta1 || eta > eta2)
			continue; // eta cut

		Double_t phi_module_corr = 1;
		int isub = -1;
		if( eta < 0 )
			isub = 0;
		if( eta > 0 )
			isub = 1;
		if( IsPhiModule == kTRUE){
			phi_module_corr = h_phi_module[fCBin][isub]->GetBinContent( (h_phi_module[fCBin][isub]->GetXaxis()->FindBin( phi ) )  );
		}
		Double_t effCorr = fEfficiency->GetCorrection( pt, fEffFilterBit, fCent );

		Qn += TComplex( 1./effCorr * phi_module_corr * TMath::Cos(ih*phi), 1./effCorr * phi_module_corr * TMath::Sin(ih*phi) );
		Sub_Ntrk += 1./effCorr * phi_module_corr ;
	}

	if( ih !=0)
		Qn /= Sub_Ntrk; // Use Qn[0] as total number of tracks(*eff)

	return Qn;
}
///________________________________________________________________________
Double_t AliJFFlucAnalysis::Get_QC_Vn(Double_t QnA_real, Double_t QnA_img, Double_t QnB_real, Double_t QnB_img )
{

	Double_t QAB_real = QnA_real* QnB_real + QnA_img * QnB_img ;
	//Double_t QAB_img = QnA_img * QnB_real - QnA_real * QnB_img;

	//Double_t QAB_abs = TMath::Sqrt( QAB_real * QAB_real + QAB_img * QAB_img );
	//Double_t QC_Vn = TMath::Sqrt(QAB_abs);
	Double_t QC_Vn = TMath::Sqrt(QAB_real);
	return QC_Vn;
}
//________________________________________________________________________
Double_t AliJFFlucAnalysis::Get_Qn_Real_pt(Double_t eta1, Double_t eta2, int harmonics, int ptbin, Double_t pt_min, Double_t pt_max)
{
	if( eta1 > eta2) cout << "ERROR eta1 should be smaller than eta2!!!" << endl;
	int nh = harmonics;
	Double_t Qn_real = 0;
	Double_t Sub_Ntrk =0;

	Long64_t ntracks = fInputList->GetEntriesFast();
	for( Long64_t it=0; it< ntracks; it++){
		AliJBaseTrack *itrack = (AliJBaseTrack*)fInputList->At(it); // load track
		Double_t eta = itrack->Eta();
		Double_t pt = itrack->Pt();
		if( pt > pt_min && pt < pt_max){
			if( eta>eta1 &&  eta<eta2 ){
				int isub = -1;
				if( eta < 0 )
					isub = 0;
				if( eta > 0 )
					isub = 1;
				Double_t phi = itrack->Phi();
				Double_t phi_module_corr =1;
				if( IsPhiModule == kTRUE){
					phi_module_corr = h_phi_module[fCBin][isub]->GetBinContent( (h_phi_module[fCBin][isub]->GetXaxis()->FindBin( phi )) );
				}
				Double_t pt = itrack->Pt();
				Double_t effCorr = fEfficiency->GetCorrection( pt, fEffFilterBit, fCent);
				Qn_real += 1.0 / effCorr * phi_module_corr * TMath::Cos( nh * phi);
				Sub_Ntrk = Sub_Ntrk + 1.0 / effCorr * phi_module_corr;
			}
		}
	}
	Qn_real /= Sub_Ntrk;
	int iside = 0; // eta -
	if (eta1 > 0 )
		iside = 1; // eta +
	NSubTracks_pt[iside][ptbin] = Sub_Ntrk; // it will be overwrite for each harmonics but should be same.

	return Qn_real;
}
//________________________________________________________________________
Double_t AliJFFlucAnalysis::Get_Qn_Img_pt(Double_t eta1, Double_t eta2, int harmonics, int ptbin, Double_t pt_min, Double_t pt_max)
{
	int nh = harmonics;
	Double_t Qn_img = 0;
	Double_t Sub_Ntrk =0;

	Long64_t ntracks = fInputList->GetEntriesFast();
	for( Long64_t it=0; it< ntracks; it++){
		AliJBaseTrack *itrack = (AliJBaseTrack*)fInputList->At(it); // load track
		Double_t eta = itrack->Eta();
		Double_t pt = itrack->Pt();
		if(pt>pt_min && pt<pt_max){
			if( eta > eta1 && eta < eta2 ){
				int isub = -1;
				if( eta < 0 ) isub = 0;
				if( eta > 0 ) isub = 1;
				Double_t phi = itrack->Phi();
				Double_t phi_module_corr =1;
				if( IsPhiModule == kTRUE){
					phi_module_corr = h_phi_module[fCBin][isub]->GetBinContent( (h_phi_module[fCBin][isub]->GetXaxis()->FindBin( phi )) );
				}
				Double_t pt = itrack->Pt();
				Double_t effCorr = fEfficiency->GetCorrection( pt, fEffFilterBit, fCent);
				Qn_img += 1./ effCorr * phi_module_corr * TMath::Sin( nh * phi);
				Sub_Ntrk = Sub_Ntrk + 1./ effCorr * phi_module_corr;
			}
		}
	}
	Qn_img /= Sub_Ntrk;
	return Qn_img;
}
///________________________________________________________________________
/* new Function for QC method
   Please see Generic Framwork from Ante
   use Standalone method  */
//________________________________________________________________________
void AliJFFlucAnalysis::CalculateQvectorsQC(){
	// calcualte Q-vector for QC method ( no subgroup )
	//init
	for(int ih=0; ih<kNH; ih++){
		QvectorQC[ih] = TComplex(0, 0);
		for(int isub=0; isub<2; isub++){
			QvectorQCeta10[ih][isub] = TComplex(0, 0);
		}
	} // for max harmonics
	//Calculate Q-vector with particle loop
	Long64_t ntracks = fInputList->GetEntriesFast(); // all tracks from Task input
	for( Long64_t it=0; it<ntracks; it++){
		AliJBaseTrack *itrack = (AliJBaseTrack*)fInputList->At(it); // load track
		Double_t phi = itrack->Phi();
		Double_t eta = itrack->Eta();
		// track Eta cut Note! pt cuts already applied in AliJFFlucTask.cxx
		// Do we need arbitary Eta cut for QC method?
		// fixed eta ranged -0.8 < eta < 0.8 for QC
//				if( TMath::Abs(eta) > fEta_max || TMath::Abs(eta) < fEta_min ) continue; << this is SP cut
//				if( TMath::Abs(eta) > 0.8 ) continue;  //   << this is old QC cut
		// we need more configuration for to study eta dep of SC(m,n) with QC method.
		// eta cut among all tracks (this is not same with SC(m,n) SP method (SP method needs symmetric eta range)//
		if( eta < fQC_eta_cut_min || eta > fQC_eta_cut_max)
			continue;
		/////////////////////////////////////////////////

		for(int ih=0; ih<kNH; ih++){
			//for(int ik=0; ik<nKL; ik++){
			TComplex q = TComplex( TMath::Cos(ih*phi), TMath::Sin(ih*phi) );
			QvectorQC[ih] += q;
			if( TMath::Abs(eta) > fQC_eta_gap_half ){  // this is for normalized SC ( denominator needs an eta gap )
				int isub = 0;
				if( eta > 0 )
					isub = 1;
				QvectorQCeta10[ih][isub] += q;
			}
			//}
		}
	} // track loop done.
}
//________________________________________________________________________
TComplex AliJFFlucAnalysis::Q(int n, int p){
	// Retrun QvectorQC
	// Q{-n, p} = Q{n, p}*
	if(n >= 0)
		return QvectorQC[n];//[p];
	return TComplex::Conjugate( QvectorQC[-n] );//[p] );
}
//________________________________________________________________________
TComplex AliJFFlucAnalysis::Two(int n1, int n2 ){
	// two-particle correlation <exp[i(n1*phi1 + n2*phi2)]>
	//	cout << "TWO FUNCTION " << Q(n1,1) << "*" << Q(n2,1) << " - " << Q(n1+n2 , 2) << endl;
	TComplex two = Q(n1, 1) * Q(n2, 1) - Q( n1+n2, 2);
	return two;
}
//________________________________________________________________________
TComplex AliJFFlucAnalysis::Four( int n1, int n2, int n3, int n4){
	TComplex four =
		Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*Q(n1+n3,2)*Q(n4,1)
		- Q(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.*Q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)
		+ Q(n2+n3,2)*Q(n1+n4,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)+Q(n1+n3,2)*Q(n2+n4,2)
		+ 2.*Q(n3,1)*Q(n1+n2+n4,3)-Q(n1,1)*Q(n2,1)*Q(n3+n4,2)+Q(n1+n2,2)*Q(n3+n4,2)
		+ 2.*Q(n2,1)*Q(n1+n3+n4,3)+2.*Q(n1,1)*Q(n2+n3+n4,3)-6.*Q(n1+n2+n3+n4,4);
	return four;
}
//__________________________________________________________________________
void AliJFFlucAnalysis::SetPhiModuleHistos( int cent, int sub, TH1D *hModuledPhi){
	// hPhi histo setter
	h_phi_module[cent][sub] = hModuledPhi;
}
