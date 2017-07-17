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


ClassImp(AliJFFlucAnalysis)

		//________________________________________________________________________
AliJFFlucAnalysis::AliJFFlucAnalysis()
	: AliAnalysisTaskSE(),
	fInputList(0),
	//h_phi_module(),
	fVertex(0),
	fCent(0),
	fNJacek(0),
	fHMG(NULL),
	fEfficiency(0), // pointer to tracking efficiency
	fEffMode(0),
	fEffFilterBit(0),
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
	fh_vn_vn()
{
	const int NCent = 7;
	Double_t CentBin[NCent+1] = {0, 5, 10, 20, 30, 40, 50, 60};
	fNCent = NCent;
	fCentBin = new Double_t[fNCent+1];
	for(int ic=0; ic<=NCent; ic++){
		fCentBin[ic] = CentBin[ic];
	}


	// pt bins to check pt dist copied from AliJHistos
	const int nJacek = 73;
	Double_t pttJacek[nJacek+1] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
		1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9,
		10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 45, 50, 60, 70, 80, 90, 100};

	fNJacek = nJacek;
	fPttJacek = new Double_t[fNJacek+1] ;
	for(int i=0; i<= fNJacek; i++){
		fPttJacek[i] = pttJacek[i];
	}

	// Constructor
}
//________________________________________________________________________
AliJFFlucAnalysis::AliJFFlucAnalysis(const char *name)
	: AliAnalysisTaskSE(name),
	fInputList(0),
	//h_phi_module(),
	fVertex(0),
	fCent(0),
	fNJacek(0),
	fHMG(NULL),
	fEfficiency(0),
	fEffMode(0),
	fEffFilterBit(0),
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
	fh_vn_vn()
{
	cout << "analysis task created " << endl;
	const int NCent = 7;
	Double_t CentBin[NCent+1] = {0, 5, 10, 20, 30, 40, 50, 60};
	fNCent = NCent;
	fDebugLevel = 0;
	fCent = -1;
	fCBin = -1;
	fEffMode = 0;
	fEffFilterBit =0;
	fInFileName ="";
	Bool_t IsPhiModule = kFALSE;
	Bool_t IsSCptdep = kFALSE;
	fEta_min = 0;
	fEta_max = 0;
	fQC_eta_cut_min = -0.8; // default setting
	fQC_eta_cut_max = 0.8; // default setting
	fImpactParameter = -1;

	fCentBin = new Double_t[fNCent+1];
	for(int ic=0; ic<=NCent; ic++){
			fCentBin[ic] = CentBin[ic];
	}

	for(int icent=0; icent<NCent; icent++){
			for(int isub=0; isub<2; isub++){
					h_phi_module[icent][isub]=NULL;
			}
	}

	// pt bins to check pt dist copied from AliJHistos
	const int nJacek = 73;
	Double_t pttJacek[nJacek+1] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
			1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9,
			10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 45, 50, 60, 70, 80, 90, 100};

	fNJacek = nJacek;
	fPttJacek = new Double_t[fNJacek+1] ;
	for(int i=0; i<= fNJacek; i++){
			fPttJacek[i] = pttJacek[i];
	}



	// Constructor
	// Define input and output slots here
}

//________________________________________________________________________
AliJFFlucAnalysis::AliJFFlucAnalysis(const AliJFFlucAnalysis& a):
	AliAnalysisTaskSE(a.GetName()),
	fInputList(a.fInputList),
	//h_phi_module(a.h_phi_module),
	fVertex(a.fVertex),
	fCent(a.fCent),
	fHMG(a.fHMG),
	fEfficiency(a.fEfficiency),
	fEffMode(a.fEffMode),
	fEffFilterBit(a.fEffFilterBit),
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
	fh_vn_vn(a.fh_vn_vn)
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
	// Create histograms
	// Called once
	// need to fill to book a histo
	fHMG = new AliJHistManager("AliJFFlucHistManager","test");
	// set AliJBin here //
	fBin_Subset .Set("Sub","Sub","Sub:%d", AliJBin::kSingle).SetBin(2);
	fBin_h .Set("NH","NH","NH:%d", AliJBin::kSingle).SetBin(kNH);
	fBin_k .Set("K","K","K:%d", AliJBin::kSingle).SetBin(nKL);

	fBin_hh .Set("NHH","NHH","NHH:%d", AliJBin::kSingle).SetBin(kNH);
	fBin_kk .Set("KK","KK","KK:%d", AliJBin::kSingle).SetBin(nKL);

	fHistCentBin .Set("CentBin","CentBin","Cent:%d",AliJBin::kSingle).SetBin(fNCent);
	fVertexBin .Set("Vtx","Vtx","Vtx:%d", AliJBin::kSingle).SetBin(3);
	fCorrBin .Set("C", "C","C:%d", AliJBin::kSingle).SetBin(17);

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
	fh_Qvector
		<< TH1D("h_QVector", "h_QVector", 100, -10, 10)
		<< fHistCentBin << fBin_Subset
		<< fBin_h
		<< "END" ;

	fh_ntracks
		<< TH1D("h_tracks", "h_tracks", 100, 0, 30000)
		<< fHistCentBin
		<< "END" ;

	fh_vn
		<< TH1D("hvn","hvn", 1024, -1.5, 1.5)
		<< fBin_h << fBin_k
		<< fHistCentBin
		<< "END";   // histogram of vn_h^k values for [ih][ik][iCent]
	fh_vn_vn
		<< TH1D("hvn_vn", "hvn_vn", 1024, -1.5, 1.5)
		<< fBin_h << fBin_k
		<< fBin_hh << fBin_kk
		<< fHistCentBin
		<< "END";  // histo of < vn * vn > for [ih][ik][ihh][ikk][iCent]
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
	fh_QvectorQC
		<< TH2D("hQvecQC", "hQvecQC", 1024, -1.1 , 1.1, 1024, -1.1, 1.1 )
		<< fBin_h
		<< fHistCentBin
		<< "END" ;

	fh_QvectorQCphi
		<< TH1D("hQbecQCphi", "hQbecQCphi", 1024, -3.2 , 3.2 )
		<< fBin_h
		<< fHistCentBin
		<< "END" ;
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
	delete []fPttJacek;
	delete []fCentBin;
}

//________________________________________________________________________
void AliJFFlucAnalysis::UserExec(Option_t *) {
	// Main loop
	// init
	for(int ih=0; ih<kNH; ih++){
		for(int im=0; im<3; im++){ //method
			fSingleVn[ih][im] = -9999;
		}
	}
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
	Eta_config[kSubB][kMin] = -1*fEta_max; // -0.8  min for SubB
	Eta_config[kSubB][kMax] = -1*fEta_min; // -0.4  max for SubB

	// use complex variable instead of doulbe Qn //
	TComplex QnA[kNH];
	TComplex QnB[kNH];
	TComplex QnB_star[kNH];
	//---------------- Do initialize here -----------
	for(int ih=0; ih<kNH; ih++){
		QnA[ih]= TComplex(0,0);
		QnB[ih]= TComplex(0,0);
		QnB_star[ih] = TComplex(0,0);
	}
	//--------------- Calculate Qn--------------------
	for(int ih=0; ih<kNH; ih++){
		QnA[ih] = CalculateQnSP( Eta_config[kSubA][kMin], Eta_config[kSubA][kMax], ih);
		QnB[ih] = CalculateQnSP( Eta_config[kSubB][kMin], Eta_config[kSubB][kMax], ih);
		//fh_Qvector[fCBin][0][ih]->Fill( QnA[ih].Theta() );
		//fh_Qvector[fCBin][1][ih]->Fill( QnB[ih].Theta() );
		QnB_star[ih] = TComplex::Conjugate ( QnB[ih] ) ;
	}
	NSubTracks[kSubA] = QnA[0].Re(); // this is number of tracks in Sub A
	NSubTracks[kSubB] = QnB[0].Re(); // this is number of tracks in Sub B
	//-------------- Fill histos with below Values ----
	// v2^2 :  k=1  /// remember QnQn = vn^(2k) not k
	// use k=0 for check v2, v3 only
	Double_t vn2[kNH][nKL];
	Double_t vn2_vn2[kNH][nKL][kNH][nKL];
	//initiation
	for(int ih=0; ih<kNH; ih++){
		for(int ik=0; ik<nKL ; ik++){
			vn2[ih][ik] =    -999;
			for(int ihh=0; ihh<kNH; ihh++){
				for(int ikk=0; ikk<nKL; ikk++){
					vn2_vn2[ih][ik][ihh][ikk] = -999;
				}
			}
		}
	}
	// calculate vn^2k
	for(int ih=2; ih<kNH; ih++){
		for(int ik=0; ik<nKL; ik++){ // 2k(0) =1, 2k(1) =2, 2k(2)=4....
			if(ik==0){
				vn2[ih][ik] = TMath::Sqrt( ( ( QnA[ih] * QnB_star[ih] ).Re() ) );
				fSingleVn[ih][0] = vn2[ih][ik]; // fill single vn with SP as method 0
			}else{
				TComplex QnAk;
				TComplex QnBstark;
				QnAk = TComplex::Power( QnA[ih], ik);
				QnBstark = TComplex::Power(QnB_star[ih], ik);
				vn2[ih][ik] = ( QnAk * QnBstark ).Re();
			}
		}
	}
	// vn^2k calcualted for n.... k....
	// calculate hvn_vn (2 combination of vn)
	for( int ih=2; ih<kNH; ih++){
		for( int ik=1; ik<nKL; ik++){
			for( int ihh=2; ihh<kNH; ihh++){
				for(int ikk=1; ikk<nKL; ikk++){
					vn2_vn2[ih][ik][ihh][ikk] = (TComplex::Power( QnA[ih]*QnB_star[ih],ik)*TComplex::Power(QnA[ihh]*QnB_star[ihh],ikk) ).Re();
				}
			}
		}
	}
	//************************************************************************
	// doing this
	//Fill the Histos here
	Double_t ebe_2p_weight = 1;
	Double_t ebe_4p_weight = 1;
	if( IsEbEWeighted == kTRUE ){
		ebe_2p_weight = NSubTracks[kSubA] * NSubTracks[kSubB] ;
		ebe_4p_weight = NSubTracks[kSubA]* NSubTracks[kSubB] * (NSubTracks[kSubA]-1) * (NSubTracks[kSubB]-1) ;
	}

	for(int ih=2; ih< kNH; ih++){
		for(int ik=0; ik<nKL; ik++){
			if(vn2[ih][ik] != -999)
				fh_vn[ih][ik][fCBin]->Fill( vn2[ih][ik] , ebe_2p_weight ); // Fill hvn2
		}
	}

	for( int ih=2; ih<kNH; ih++){
		for( int ik=1; ik<nKL; ik++){
			for( int ihh=2; ihh<kNH; ihh++){
				for(int ikk=1; ikk<nKL; ikk++){
					if(vn2_vn2[ih][ik][ihh][ikk] != -999 )
						fh_vn_vn[ih][ik][ihh][ikk][fCBin]->Fill( vn2_vn2[ih][ik][ihh][ikk], ebe_4p_weight ) ; // Fill hvn_vn
				}
			}
		}
	}
	///	Fill more correlators in manualy
	TComplex V4V2starv2_2 =	QnA[4] *TComplex::Power( QnB_star[2] ,2) * vn2[2][1] ;
	TComplex V4V2starv2_4 = QnA[4] * TComplex::Power( QnB_star[2], 2) * vn2[2][2] ;
	TComplex V4V2star = QnA[4] * TComplex::Power( QnB_star[2], 2 );
	TComplex V5V2starV3starv2_2 = QnA[5] * QnB_star[2] * QnB_star[3] * vn2[2][1] ;
	TComplex V5V2starV3star = QnA[5] * QnB_star[2] * QnB_star[3] ;
	TComplex V5V2starV3startv3_2 = QnA[5] * QnB_star[2] * QnB_star[3] * vn2[3][1];
	TComplex V6V2star_3 = QnA[6] * TComplex::Power( QnB_star[2] , 3) ;
	TComplex V6V3star_2 = QnA[6] * TComplex::Power( QnB_star[3], 2) ;
	TComplex V7V2star_2V3star = QnA[7] * TComplex::Power( QnB_star[2] , 2) * QnB_star[3];


	// New correlattors (Modified by You's corretion term for self-correlations)
	TComplex nV4V2star = (QnA[4] * QnB_star[2] * QnB_star[2]) -( 1./(NSubTracks[1]-1) * QnA[4] * QnB_star[4] );
	TComplex nV5V2starV3star = (QnA[5] * QnB_star[2] * QnB_star[3])- (1/(NSubTracks[1]-1) * QnA[5] * QnB_star[5]);
	TComplex nV6V3star_2 = (QnA[6] * QnB_star[3] * QnB_star[3]) - (1/(NSubTracks[1]-1) * QnA[6] * QnB_star[6] );



	// New correlattors (Modifed by Ante's correction term for self-correlations for SC result)
	TComplex nV4V4V2V2 = (QnA[4]*QnB_star[4]*QnA[2]*QnB_star[2]) - ((1/(NSubTracks[1]-1) * QnB_star[6] * QnA[4] *QnA[2] ))
		- ((1/(NSubTracks[0]-1) * QnA[6]*QnB_star[4] * QnB_star[2])) + (1/((NSubTracks[0]-1)*(NSubTracks[1]-1))*QnA[6]*QnB_star[6] );
	TComplex nV3V3V2V2 = (QnA[3]*QnB_star[3]*QnA[2]*QnB_star[2]) - ((1/(NSubTracks[1]-1) * QnB_star[5] * QnA[3] *QnA[2] ))
		- ((1/(NSubTracks[0]-1) * QnA[5]*QnB_star[3] * QnB_star[2])) + (1/((NSubTracks[0]-1)*(NSubTracks[1]-1))*QnA[5]*QnB_star[5] );
	// add higher order SC results
	TComplex nV5V5V2V2 = (QnA[5]*QnB_star[5]*QnA[2]*QnB_star[2]) - ((1/(NSubTracks[1]-1) * QnB_star[7] * QnA[5] *QnA[2] ))
		- ((1/(NSubTracks[0]-1) * QnA[7]*QnB_star[5] * QnB_star[2])) + (1/((NSubTracks[0]-1)*(NSubTracks[1]-1))*QnA[7]*QnB_star[7] );
	TComplex nV5V5V3V3 = (QnA[5]*QnB_star[5]*QnA[3]*QnB_star[3]) - ((1/(NSubTracks[1]-1) * QnB_star[8] * QnA[5] *QnA[3] ))
		- ((1/(NSubTracks[0]-1) * QnA[8]*QnB_star[5] * QnB_star[3])) + (1/((NSubTracks[0]-1)*(NSubTracks[1]-1))*QnA[8]*QnB_star[8] );
	TComplex nV4V4V3V3 = (QnA[4]*QnB_star[4]*QnA[3]*QnB_star[3]) - ((1/(NSubTracks[1]-1) * QnB_star[7] * QnA[4] *QnA[3] ))
		- ((1/(NSubTracks[0]-1) * QnA[7]*QnB_star[4] * QnB_star[3])) + (1/((NSubTracks[0]-1)*(NSubTracks[1]-1))*QnA[7]*QnB_star[7] );



	fh_correlator[0][fCBin]->Fill( V4V2starv2_2.Re() );
	fh_correlator[1][fCBin]->Fill( V4V2starv2_4.Re() );
	fh_correlator[2][fCBin]->Fill( V4V2star.Re() ) ; // added 2015.3.18
	fh_correlator[3][fCBin]->Fill( V5V2starV3starv2_2.Re() );
	fh_correlator[4][fCBin]->Fill( V5V2starV3star.Re() );
	fh_correlator[5][fCBin]->Fill( V5V2starV3startv3_2.Re() );
	fh_correlator[6][fCBin]->Fill( V6V2star_3.Re() );
	fh_correlator[7][fCBin]->Fill( V6V3star_2.Re() );
	fh_correlator[8][fCBin]->Fill( V7V2star_2V3star.Re() ) ;

	fh_correlator[9][fCBin]->Fill( nV4V2star.Re() ); // added 2015. 6. 10
	fh_correlator[10][fCBin]->Fill( nV5V2starV3star.Re() );
	fh_correlator[11][fCBin]->Fill( nV6V3star_2.Re() ) ;

	// use this to avoid self-correlation 4p correlation (2 particles from A, 2 particles from B) -> MA(MA-1)MB(MB-1) : evt weight..
	fh_correlator[12][fCBin]->Fill( nV4V4V2V2.Re() , ebe_4p_weight);
	fh_correlator[13][fCBin]->Fill( nV3V3V2V2.Re() , ebe_4p_weight);

	fh_correlator[14][fCBin]->Fill( nV5V5V2V2.Re() , ebe_4p_weight);
	fh_correlator[15][fCBin]->Fill( nV5V5V3V3.Re() , ebe_4p_weight);
	fh_correlator[16][fCBin]->Fill( nV4V4V3V3.Re() , ebe_4p_weight);


	if(IsSCptdep == kTRUE){
		const int SCNH =6; // 0, 1, 2(v2), 3(v3), 4(v4), 5(v5)
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

			//add hihger order 5225, 5335, 4334
			TComplex nV5V5V3V3_pt =  (QnA_pt[5][ipt]*QnB_pt_star[5][ipt]*QnA_pt[3][ipt]*QnB_pt_star[3][ipt]) - ((1/(NSubTracks_pt[1][ipt]-1) * QnB_pt_star[8][ipt] * QnA_pt[5][ipt] *QnA_pt[3][ipt] ))
				- ((1/(NSubTracks_pt[0][ipt]-1) * QnA_pt[8][ipt]*QnB_pt_star[5][ipt] * QnB_pt_star[3][ipt])) + (1/((NSubTracks_pt[0][ipt]-1)*(NSubTracks_pt[1][ipt]-1))*QnA_pt[8][ipt]*QnB_pt_star[8][ipt] );
			TComplex nV5V5V2V2_pt =  (QnA_pt[5][ipt]*QnB_pt_star[5][ipt]*QnA_pt[2][ipt]*QnB_pt_star[2][ipt]) - ((1/(NSubTracks_pt[1][ipt]-1) * QnB_pt_star[7][ipt] * QnA_pt[5][ipt] *QnA_pt[2][ipt] ))
				- ((1/(NSubTracks_pt[0][ipt]-1) * QnA_pt[7][ipt]*QnB_pt_star[5][ipt] * QnB_pt_star[2][ipt])) + (1/((NSubTracks_pt[0][ipt]-1)*(NSubTracks_pt[1][ipt]-1))*QnA_pt[7][ipt]*QnB_pt_star[7][ipt] );
			TComplex nV4V4V3V3_pt =  (QnA_pt[4][ipt]*QnB_pt_star[4][ipt]*QnA_pt[3][ipt]*QnB_pt_star[3][ipt]) - ((1/(NSubTracks_pt[1][ipt]-1) * QnB_pt_star[7][ipt] * QnA_pt[4][ipt] *QnA_pt[3][ipt] ))
				- ((1/(NSubTracks_pt[0][ipt]-1) * QnA_pt[7][ipt]*QnB_pt_star[4][ipt] * QnB_pt_star[3][ipt])) + (1/((NSubTracks_pt[0][ipt]-1)*(NSubTracks_pt[1][ipt]-1))*QnA_pt[7][ipt]*QnB_pt_star[7][ipt] );


			int ik=1;
			fh_SC_ptdep_4corr[2][1][4][1][fCBin][ipt]->Fill( nV4V4V2V2_pt.Re());
			fh_SC_ptdep_4corr[2][1][3][1][fCBin][ipt]->Fill( nV3V3V2V2_pt.Re());
			fh_SC_ptdep_4corr[2][1][5][1][fCBin][ipt]->Fill( nV5V5V2V2_pt.Re() );
			fh_SC_ptdep_4corr[3][1][4][1][fCBin][ipt]->Fill( nV4V4V3V3_pt.Re() );
			fh_SC_ptdep_4corr[3][1][5][1][fCBin][ipt]->Fill( nV5V5V3V3_pt.Re() ) ;
		}
	}//pt dep done

	if(IsSCwithQC==kTRUE){
		// (a) calculate QC q-vecotr
		// (b) calculate 4p correaltion
		// (c) calculate 2p correaltion
		// (d) calculate 2p corrleaion with |dEta|>1.0 (for normalized SC)

		//(a)
		CalculateQvectorsQC();
		//(b)
		Double_t QC_4p_value[6][6] = {-99.};
		Double_t QC_2p_value[6] = {-99.};
		TComplex four; // declear here
		TComplex two;

		for(int ih=2; ih<=5; ih++){
			for(int ihh=2; ihh<ih; ihh++){
				Double_t event_weight = 1;
				if( IsEbEWeighted == kTRUE)
					event_weight = Four(0,0,0,0).Re();
				four = TComplex(-99, -99); // initialize
				four = Four( ih, ihh, -1*ih, -1*ihh ) / Four(0,0,0,0).Re();
				fh_SC_with_QC_4corr[ih][ihh][fCBin]->Fill( four.Re(), event_weight );
				QC_4p_value[ih][ihh] = four.Re();
			}
		};
		//(c)
		for(int ih=2; ih<=5; ih++){
			// Finally we want 2p corr as like
			// 1/( M*(M-1) ) * [ QQ* - M ]
			// two(2,2) = Q2 Q2* - Q0 = Q2Q2* - M
			// two(0,0) = Q0 Q0* - Q0 = M^2 - M
			Double_t event_weight = 1;
			if( IsEbEWeighted == kTRUE)
				event_weight = Two(0,0).Re();
			two = TComplex(-99, -99); // initialize
			two = Two(ih, -1*ih) / Two(0,0).Re();
			fh_SC_with_QC_2corr[ih][fCBin]->Fill( two.Re(), event_weight );
			QC_2p_value[ih] = two.Re();
			// fill single vn  with QC without EtaGap as mehtod 2
			fSingleVn[ih][2] = TMath::Sqrt(two.Re() );
		}
		//(d)
		for(int ih=2; ih<=5; ih++){
			Double_t event_weight = 1;
			if( IsEbEWeighted == kTRUE)event_weight = (QvectorQCeta10[0][1][kSubA]*QvectorQCeta10[0][1][kSubB] ).Re();
						 // of course Qvec[0] * Qvec[0].. -> real (all sin term =0) //
			TComplex two = (QvectorQCeta10[ih][1][kSubA]*TComplex::Conjugate(QvectorQCeta10[ih][1][kSubB])) / (QvectorQCeta10[0][1][kSubA]*QvectorQCeta10[0][1][kSubB]).Re();
							// is same as divided by event_weight.(number of comniations)
			fh_SC_with_QC_2corr_eta10[ih][fCBin]->Fill( two.Re(), event_weight );
			// fill single vn with QC method with Eta Gap as method 1
			fSingleVn[ih][1] = TMath::Sqrt(two.Re() );

		}
		//Check evt-by-evt SP/QC ratio. (term-by-term)
		// calculate  (vn^2 vm^2)_SP /  (vn^2 vm^2)_QC
		// 4p ( v3v3v2v2, v4v4v2v2, v5v5v2v2, v5v5v3v3, v4v4v3v3
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
	for(int iaxis=0; iaxis<3; iaxis++){
		fh_vertex[iaxis]->Fill(  fVertex[iaxis] );
	}
	fh_TrkQA_TPCvsCent->Fill( fCent, fTPCtrks);
	fh_TrkQA_TPCvsGlob->Fill( fGlbtrks, fTPCtrks);
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
	Double_t QAB_img = QnA_img * QnB_real - QnA_real * QnB_img;

	Double_t QAB_abs = TMath::Sqrt( QAB_real * QAB_real + QAB_img * QAB_img );
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
		for(int ik=0; ik<nKL; ik++){
			QvectorQC[ih][ik] = TComplex(0, 0);
			for(int isub=0; isub<2; isub++){
				QvectorQCeta10[ih][ik][isub] = TComplex(0, 0);
			}
		} // for max power
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
			for(int ik=0; ik<nKL; ik++){
				QvectorQC[ih][ik] += TComplex( TMath::Cos(ih*phi), TMath::Sin(ih*phi) );
				// this is not working (there are no eta gap for +0.6, +0.61 in this way..
				// fix this as like SP -> 2 sub event //
				if( TMath::Abs(eta) > 0.5 ){  // this is for Noramlized SC ( denominator need eta gap )
					int isub = 0;
					if( eta > 0 )
						isub = 1; // what about eta=0?
					QvectorQCeta10[ih][ik][isub] += TComplex( TMath::Cos(ih*phi), TMath::Sin(ih*phi) );
				}
			}
		}
	} // track loop done.


	// Q-vector[ih][ik] calculated doen. //(need ik??)

	/*
	// Save QA plot
	for(int ih=2; ih<kNH; ih++){
	fh_QvectorQC[ih][fCBin]->Fill( QvectorQC[ih][1].Re()/QvectorQC[0][1].Re() , QvectorQC[ih][1].Im()/QvectorQC[0][1].Re() ); // fill normalized Q vector
	fh_QvectorQCphi[ih][fCBin]->Fill( QvectorQC[ih][1].Theta() );
	}
	// Q-vector calculated
	 */
}
//________________________________________________________________________
TComplex AliJFFlucAnalysis::Q(int n, int p){
	// Retrun QvectorQC
	// Q{-n, p} = Q{n, p}*
	if(n >= 0)
		return QvectorQC[n][p];
	return TComplex::Conjugate( QvectorQC[-n][p] );
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
