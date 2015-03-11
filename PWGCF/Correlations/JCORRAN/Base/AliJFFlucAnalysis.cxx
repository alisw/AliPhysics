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


ClassImp(AliJFFlucAnalysis)

//________________________________________________________________________
AliJFFlucAnalysis::AliJFFlucAnalysis() 
	: AliAnalysisTaskSE(), 
	fInputList(0),
	fVertex(0),
	fCent(0),
	fNJacek(0),	
	fHMG(NULL),
	fBin_h(),
	fBin_k(),
	fBin_hh(),
	fBin_kk(),
	fHistCentBin(),
	fh_cent(),
	fh_vertex(),
	fh_eta(),
	fh_ntracks(),
	fh_vn(),
	fh_vn_vn(),
	fh_vn_test1(),
	fh_vn_test2(),
	fh_vn_vn_test1(),
	fh_vn_vn_test2()	
{
	const int NCent = 7;
	double CentBin[NCent+1] = {0, 5, 10, 20, 30, 40, 50, 60};
	fNCent = NCent;
	fCentBin = new double[fNCent+1];
	for(int ic=0; ic<=NCent; ic++){
		fCentBin[ic] = CentBin[ic];
	}	


	// pt bins to check pt dist copied from AliJHistos
	const int nJacek = 73;
	double pttJacek[nJacek+1] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
      1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9,
      10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 45, 50, 60, 70, 80, 90, 100};
	
	fNJacek = nJacek;
	fPttJacek = new double[fNJacek+1] ;
	for(int i=0; i<= fNJacek; i++){
		fPttJacek[i] = pttJacek[i];
	}




	// Constructor
}
//________________________________________________________________________
AliJFFlucAnalysis::AliJFFlucAnalysis(const char *name) 
	: AliAnalysisTaskSE(name), 
	fInputList(0),
	fVertex(0),
	fCent(0),
	fNJacek(0),
	fHMG(NULL),
	fBin_h(),
	fBin_k(),
	fBin_hh(),
	fBin_kk(),
	fHistCentBin(),
	fh_cent(),
	fh_vertex(),
	fh_eta(),
	fh_ntracks(),
	fh_vn(),
	fh_vn_vn(),
	fh_vn_test1(),
	fh_vn_test2(),
	fh_vn_vn_test1(),
	fh_vn_vn_test2()
{
 
	const int NCent = 7;
	double CentBin[NCent+1] = {0, 5, 10, 20, 30, 40, 50, 60};
	fNCent = NCent;
	fCentBin = new double[fNCent+1];
	for(int ic=0; ic<=NCent; ic++){
		fCentBin[ic] = CentBin[ic];
	}	


	// pt bins to check pt dist copied from AliJHistos
	const int nJacek = 73;
	double pttJacek[nJacek+1] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
      1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9,
      10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 45, 50, 60, 70, 80, 90, 100};
	
	fNJacek = nJacek;
	fPttJacek = new double[fNJacek+1] ;
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
	fVertex(a.fVertex),
	fCent(a.fCent),
	fHMG(a.fHMG),
	fBin_h(a.fBin_h),
	fBin_k(a.fBin_k),
	fBin_hh(a.fBin_hh),
	fBin_kk(a.fBin_kk),
	fHistCentBin(a.fHistCentBin),
	fh_cent(a.fh_cent),
	fh_vertex(a.fh_vertex),
	fh_eta(a.fh_eta),
	fh_ntracks(a.fh_ntracks),
	fh_vn(a.fh_vn),
	fh_vn_vn(a.fh_vn_vn),
	fh_vn_test1(a.fh_vn_test1),
	fh_vn_test2(a.fh_vn_test2),
	fh_vn_vn_test1(a.fh_vn_vn_test1),
	fh_vn_vn_test2(a.fh_vn_vn_test2)
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
}
//________________________________________________________________________
void AliJFFlucAnalysis::UserCreateOutputObjects(){
	// Create histograms
	// Called once
	AnaEntry = 0;
	// need to fill to book a histo

	fHMG = new AliJHistManager("AliJFFlucHistManager","test");


	// set AliJBin here // 
	fBin_h .Set("NH","NH","NH:%d", AliJBin::kSingle).SetBin(kNH);
	fBin_k .Set("K","K","K:%d", AliJBin::kSingle).SetBin(nKL);

	fBin_hh .Set("NHH","NHH","NHH:%d", AliJBin::kSingle).SetBin(kNH);
	fBin_kk .Set("KK","KK","KK:%d", AliJBin::kSingle).SetBin(nKL);

	fHistCentBin .Set("CentBin","CentBin","Cent:%d",AliJBin::kSingle).SetBin(fNCent);
	fVertexBin .Set("Vtx","Vtx","Vtx:%d", AliJBin::kSingle).SetBin(3);
	fCorrBin .Set("C", "C","C:%d", AliJBin::kSingle).SetBin(5);

	// set AliJTH1D here //
	fh_cent
		<< TH1D("h_cent","h_cent", 100, 0, 100) 
		<< "END" ;

	fh_vertex
		<< TH1D("h_vertex","h_vertex", 400, -20, 20)
		<< fVertexBin
		<< "END" ;
	fh_pt
		<< TH1D("hChargedPtJacek", "", fNJacek, fPttJacek)
		<< fHistCentBin 
		<< "END" ; 

	fh_eta
		<< TH1D("h_eta", "h_eta", 300, -15, 15 )
		<< fHistCentBin
		<< "END" ;

	fh_ntracks 
		<< TH1D("h_tracks", "h_tracks", 1000, 0, 30000)
		<< fHistCentBin
		<< "END" ;  
	
	fh_vn 
 		<< TH1D("hvn","hvn", 5096, -0.5, 0.5) 
 		<< fBin_h << fBin_k 
 		<< fHistCentBin
 		<< "END";   // histogram of vn_h^k values for [ih][ik][iCent] 
	fh_vn_vn
	  	<< TH1D("hvn_vn", "hvn_vn", 5096, -0.1, 0.1)
  		<< fBin_h << fBin_k 
 	 	<< fBin_hh << fBin_kk
 	 	<< fHistCentBin
 	 	<< "END";  // histo of < vn * vn > for [ih][ik][ihh][ikk][iCent] 
	fh_vn_test1 
 		<< TH1D("hvnTEST1","hvn", 5096, -0.5, 0.5) 
 		<< fBin_h << fBin_k 
 		<< fHistCentBin
 		<< "END";   // histogram of vn_h^k values for [ih][ik][iCent] 
	fh_vn_vn_test1
	  	<< TH1D("hvn_vnTEST1", "hvn_vn", 5096, -0.1, 0.1)
  		<< fBin_h << fBin_k 
 	 	<< fBin_hh << fBin_kk
 	 	<< fHistCentBin
 	 	<< "END";  // histo of < vn * vn > for [ih][ik][ihh][ikk][iCent] 
	fh_vn_test2
 		<< TH1D("hvnTEST2","hvn", 5096, -0.5, 0.5) 
 		<< fBin_h << fBin_k 
 		<< fHistCentBin
 		<< "END";   // histogram of vn_h^k values for [ih][ik][iCent] 
	fh_vn_vn_test2
	  	<< TH1D("hvn_vnTEST2", "hvn_vn", 5096, -0.1, 0.1)
  		<< fBin_h << fBin_k 
 	 	<< fBin_hh << fBin_kk
 	 	<< fHistCentBin
 	 	<< "END";  // histo of < vn * vn > for [ih][ik][ihh][ikk][iCent] 
	fh_correlator 
		<< TH1D("h_corr", "h_corr", 5096, -1, 1)
		<< fCorrBin
		<< fHistCentBin
		<< "END" ;


	//AliJTH1D set done.
	
	fHMG->Print();
	fHMG->WriteConfig();
	
}
//________________________________________________________________________
AliJFFlucAnalysis::~AliJFFlucAnalysis() {
	delete fInputList;
	delete fHMG;
	delete []fPttJacek;
	delete []fCentBin;
}

//________________________________________________________________________
void AliJFFlucAnalysis::UserExec(Option_t *) {
	// Main loop
	// Called for each event
	// DO ANALYSIS WORK HERE // 

	// find Centrality
	double inputCent = fCent;
	fCBin = -1;
	for(int iCbin=0; iCbin<fNCent; iCbin++){
		if( inputCent > fCentBin[iCbin] && inputCent < fCentBin[iCbin+1] ){fCBin = iCbin;};
	}
	if(fCBin == -1){return;};
	int trk_number = fInputList->GetEntriesFast();
	fh_ntracks[fCBin]->Fill( trk_number ) ;
	fh_cent->Fill( inputCent ) ; 
	Fill_QA_plot( -15, 15 );


	enum{kSubA, kSubB, kNSub};
	enum{kReal, kImg, kNPhase};
	double Eta_config[kNSub][2];
	Eta_config[kSubA][0] = fEta_min;
	Eta_config[kSubA][1] = fEta_max;
	Eta_config[kSubB][0] = -1*fEta_max;
	Eta_config[kSubB][1] = -1*fEta_min; 
	

	// use complex variable instead of doulbe Qn // 
	TComplex QnA[kNH];
	TComplex QnB[kNH];
	TComplex QnB_star[kNH];

	//---------------- Do initialize here -----------
	for(int ih=2; ih<kNH; ih++){
			QnA[ih]= TComplex(0,0);
			QnB[ih]= TComplex(0,0);
			QnB_star[ih] = TComplex(0,0);			
	}
	//--------------- Calculate Qn--------------------
	for(int ih=2; ih<kNH; ih++){
					double QAReal = Get_Qn_Real( Eta_config[kSubA][0], Eta_config[kSubA][1], ih );
					double QAImag = Get_Qn_Img(  Eta_config[kSubA][0], Eta_config[kSubA][1], ih );
					QnA[ih]= TComplex(QAReal, QAImag);
					double QBReal = Get_Qn_Real( Eta_config[kSubB][0], Eta_config[kSubB][1], ih );
					double QBImag = Get_Qn_Img(  Eta_config[kSubB][0], Eta_config[kSubB][1], ih );
					QnB[ih]= TComplex(QBReal, QBImag);
					QnB_star[ih] = TComplex::Conjugate ( QnB[ih] ) ;
	}
	//-------------- Fill histos with below Values ----
	// v2^2 :  k=1  /// remember QnQn = vn^(2k) not k
	// use k=0 for check v2, v3 only
	double vn2[kNH][nKL]; 
	double vn2_vn2[kNH][nKL][kNH][nKL];
	double vn2_test1[kNH][nKL]; 
	double vn2_test2[kNH][nKL]; 
	double vn2_vn2_test1[kNH][nKL][kNH][nKL];
	double vn2_vn2_test2[kNH][nKL][kNH][nKL];
	//initiation
	for(int ih=0; ih<kNH; ih++){
		for(int ik=0; ik<nKL ; ik++){
				vn2[ih][ik] = -999;
				vn2_test1[ih][ik] = -999;
				vn2_test2[ih][ik] = -999;
				for(int ihh=0; ihh<kNH; ihh++){
						for(int ikk=0; ikk<nKL; ikk++){
								vn2_vn2[ih][ik][ihh][ikk] = -999;
								vn2_vn2_test1[ih][ik][ihh][ikk] = -999;
								vn2_vn2_test2[ih][ik][ihh][ikk] = -999;
						}
				}
		}
	}
	// calculate vn^2k	
	for( int ih=2; ih<kNH; ih++){
			for(int ik=0; ik<nKL; ik++){ // 2k(0) =1, 2k(1) =2, 2k(2)=4....
					if(ik==0) vn2[ih][ik] = TMath::Sqrt( ( TComplex::Abs( QnA[ih] * QnB_star[ih] )) );
					if(ik!=0){
							TComplex QnAk;
							TComplex QnBstark;
							QnAk = TComplex::Power( QnA[ih], ik) ;
							QnBstark = TComplex::Power(QnB_star[ih], ik) ;
							vn2[ih][ik] = ( QnAk * QnBstark ).Re();  // for k=1 vn2 and vn2_test1 is same.
							vn2_test1[ih][ik] = (TComplex::Power( QnA[ih] * QnB_star[ih], ik)).Re();
							vn2_test2[ih][ik] = TMath::Power(  ( QnA[ih] * QnB_star[ih] ).Re() ,ik );
					}		
			}
	}
	// vn^2k calcualted for n.... k....
	// calculate hvn_vn (2 combination of vn) 
	for( int ih=2; ih<kNH; ih++){ 
			for( int ik=1; ik<nKL; ik++){
					for( int ihh=2; ihh<kNH; ihh++){ 
							for(int ikk=1; ikk<nKL; ikk++){
									vn2_vn2[ih][ik][ihh][ikk] = (  
													TComplex::Power( QnA[ih]*QnB_star[ih],ik)*TComplex::Power(QnA[ihh]*QnB_star[ihh],ikk) ).Re();
									vn2_vn2_test1[ih][ik][ihh][ikk] = 
											TMath::Power( (QnA[ih] * QnB_star[ih]).Re() , ik ) * 
											TMath::Power( (QnA[ihh] * QnB_star[ihh] ).Re(), ik ); 

									vn2_vn2_test2[ih][ik][ihh][ikk] = 
											( TComplex::Power( QnA[ih], ik ) * TComplex::Power( QnB_star[ih], ik ) *
											  TComplex::Power( QnA[ihh], ikk ) * TComplex::Power( QnB_star[ihh], ikk ) ).Re();

							}
					}
			}
	}
	//************************************************************************
	// doing this 
	//Fill the Histos here
	for(int ih=2; ih< kNH; ih++){
		if(vn2[ih][0] != -999)	fh_vn[ih][0][fCBin]->Fill( vn2[ih][0] );
		if(vn2[ih][0] != -999)  fh_vn_test1[ih][0][fCBin]->Fill( vn2[ih][0] );
		if(vn2[ih][0] != -999) fh_vn_test2[ih][0][fCBin]->Fill( vn2[ih][0] ); // fill k=0 for just v2, v3.. 

		for(int ik=1; ik<nKL; ik++){
			if(vn2[ih][ik] != -999)	fh_vn[ih][ik][fCBin]->Fill( vn2[ih][ik] ); // Fill hvn2
			if(vn2_test1[ih][ik] != -999 )fh_vn_test1[ih][ik][fCBin]->Fill( vn2_test1[ih][ik] ) ;
			if(vn2_test2[ih][ik] != -999 )fh_vn_test2[ih][ik][fCBin]->Fill( vn2_test2[ih][ik] ) ;
		}
	}

	for( int ih=2; ih<kNH; ih++){ 
		for( int ik=1; ik<nKL; ik++){
			for( int ihh=2; ihh<kNH; ihh++){ 
				for(int ikk=1; ikk<nKL; ikk++){
					if(vn2_vn2[ih][ik][ihh][ikk] != -999 ) fh_vn_vn[ih][ik][ihh][ikk][fCBin]->Fill( vn2_vn2[ih][ik][ihh][ikk] ) ; // Fill hvn_vn 
					if(vn2_vn2_test1[ih][ik][ihh][ikk] != -999)  fh_vn_vn_test1[ih][ik][ihh][ikk][fCBin]->Fill( vn2_vn2_test1[ih][ik][ihh][ikk] ) ;
					if(vn2_vn2_test2[ih][ik][ihh][ikk] != -999)  fh_vn_vn_test2[ih][ik][ihh][ikk][fCBin]->Fill( vn2_vn2_test2[ih][ik][ihh][ikk] ) ; 
				}
			}
		}
	}
	///	Fill more correlators in manualy
	TComplex V4V2starv2_2 =	QnA[4] *TComplex::Power( QnB_star[2] ,2) * vn2[2][1] ;
	TComplex V4V2starv2_4 = QnA[4] * TComplex::Power( QnB_star[2], 2) * vn2[2][2] ;
	TComplex V5V2starV3starv2_2 = QnA[5] * QnB_star[2] * QnB_star[3] * vn2[2][1] ;
	TComplex V5V2starV3star = QnA[5] * QnB_star[2] * QnB_star[3] ;
	TComplex V5V2starV3startv3_2 = QnA[5] * QnB_star[2] * QnB_star[3] * vn2[3][1];


	fh_correlator[0][fCBin]->Fill( V4V2starv2_2.Re() );
	fh_correlator[1][fCBin]->Fill( V4V2starv2_4.Re() );
	fh_correlator[2][fCBin]->Fill( V5V2starV3starv2_2.Re() );
	fh_correlator[3][fCBin]->Fill( V5V2starV3star.Re() );
	fh_correlator[4][fCBin]->Fill( V5V2starV3startv3_2.Re() );
	//
	//

	AnaEntry++;
	//higher debug level = detail information
	//---------------------------------------------------------------
	// check if the event was triggered or not and vertex   
}

//________________________________________________________________________
void AliJFFlucAnalysis::Terminate(Option_t *) 
{
//	fList = dynamic_cast<TList*> (GetOutputData(1));
//	if(!fList) { Printf("ERROR: fList not availabe"); return;};
	cout<<"Sucessfully Finished"<<endl;
}
//________________________________________________________________________
//________________________________________________________________________
double AliJFFlucAnalysis::Get_Qn_Real(double eta1, double eta2, int harmonics)
{
		int nh = harmonics;
		double Qn_real = 0;
		double Sub_Ntrk =0;

		Long64_t ntracks = fInputList->GetEntriesFast();
		for( Long64_t it=0; it< ntracks; it++){
			AliJBaseTrack *itrack = (AliJBaseTrack*)fInputList->At(it); // load track
			double eta = itrack->Eta();
			if( eta>eta1 &&  eta<eta2 ){ 
					Sub_Ntrk++;
					double phi = itrack->GetTwoPiPhi();
					Qn_real += TMath::Cos( nh * phi);
			}
		}
		Qn_real /= Sub_Ntrk;
		return Qn_real; 
}
//________________________________________________________________________
void AliJFFlucAnalysis::Fill_QA_plot( double eta1, double eta2 )
{
		Long64_t ntracks = fInputList->GetEntriesFast();
		for( Long64_t it=0; it< ntracks; it++){
			AliJBaseTrack *itrack = (AliJBaseTrack*)fInputList->At(it); // load track
			double eta = itrack->Eta();
			double pt = itrack->Pt();
			if( eta > eta1 && eta < eta2 ){ 
				fh_eta[fCBin]->Fill(eta);
				fh_pt[fCBin]->Fill(pt);	
			}
		}
	for(int iaxis=0; iaxis<3; iaxis++){
		fh_vertex[iaxis]->Fill(  fVertex[iaxis] );
	}
}
//________________________________________________________________________
double AliJFFlucAnalysis::Get_Qn_Img(double eta1, double eta2, int harmonics)
{
		int nh = harmonics;
		double Qn_img = 0;
		double Sub_Ntrk =0;

		Long64_t ntracks = fInputList->GetEntriesFast();
		for( Long64_t it=0; it< ntracks; it++){
			AliJBaseTrack *itrack = (AliJBaseTrack*)fInputList->At(it); // load track
			double eta = itrack->Eta();
			if( eta > eta1 && eta < eta2 ){ 
					Sub_Ntrk++;
					double phi = itrack->GetTwoPiPhi();
					Qn_img += TMath::Sin( nh * phi);
			}
		}
		Qn_img /= Sub_Ntrk;
		return Qn_img; 
}
///________________________________________________________________________
double AliJFFlucAnalysis::Get_QC_Vn(double QnA_real, double QnA_img, double QnB_real, double QnB_img )
{

		double QAB_real = QnA_real* QnB_real + QnA_img * QnB_img ;
		double QAB_img = QnA_img * QnB_real - QnA_real * QnB_img; 

		double QAB_abs = TMath::Sqrt( QAB_real * QAB_real + QAB_img * QAB_img );
		//double QC_Vn = TMath::Sqrt(QAB_abs);
		double QC_Vn = TMath::Sqrt(QAB_real);
		return QC_Vn; 
}
///________________________________________________________________________
double AliJFFlucAnalysis::Complex_product_real(double QnA_real, double QnA_img, double QnB_real, double QnB_img)
{
		// this values is  (Q_An Q*_Bn ) // output is complex
		return (QnA_real* QnB_real - QnA_img * QnB_img) ;

}
///___________________
double AliJFFlucAnalysis::Complex_product_img(double QnA_real, double QnA_img, double QnB_real, double QnB_img)
{
		// this values is  (Q_An Q*_Bn ) // output is complex
		return ( QnA_img * QnB_real + QnA_real * QnB_img ) ;

}
///___________________
double AliJFFlucAnalysis::Complex_abs(double real, double img)
{
		// this values is  (Q_An Q*_Bn ) // output is complex
		return ( TMath::Sqrt(  real* real + img * img ) ) ;

}
///___________________
double AliJFFlucAnalysis::Complex_sqr_real(double real, double img)
{
		// this value is  real-part of  (Qn)^2 
		return ( ( real* real - img * img ) ) ;
}
///___________________
double AliJFFlucAnalysis::Complex_sqr_img(double real, double img)
{
		// this value is  real-part of  (Qn)^2 
		return ( 2* real * img ) ;
}







