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
	fCent(0),
	fNJacek(0),	
	fHMG(NULL),
	fBin_h(),
	fBin_k(),
	fBin_hh(),
	fBin_kk(),
	fHistCentBin(),
	fh_cent(),
	fh_eta(),
	fh_ntracks(),
	fh_vn(),
	fh_vn_vn()	
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
	fCent(0),
	fNJacek(0),
	fHMG(NULL),
	fBin_h(),
	fBin_k(),
	fBin_hh(),
	fBin_kk(),
	fHistCentBin(),
	fh_cent(),
	fh_eta(),
	fh_ntracks(),
	fh_vn(),
	fh_vn_vn()
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
	fCent(a.fCent),
	fHMG(a.fHMG),
	fBin_h(a.fBin_h),
	fBin_k(a.fBin_k),
	fBin_hh(a.fBin_hh),
	fBin_kk(a.fBin_kk),
	fHistCentBin(a.fHistCentBin),
	fh_cent(a.fh_cent),
	fh_eta(a.fh_eta),
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

	// set AliJTH1D here //
	fh_cent
		<< TH1D("h_cent","h_cent", 100, 0, 100) 
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
 		<< TH1D("hvn","hvn", 5096, 0, 0.3) 
 		<< fBin_h << fBin_k 
 		<< fHistCentBin
 		<< "END";   // histogram of vn_h^k values for [ih][ik][iCent] 
	fh_vn_vn
	  	<< TH1D("hvn_vn", "hvn_vn", 5096, 0, 0.3)
  		<< fBin_h << fBin_k 
 	 	<< fBin_hh << fBin_kk
 	 	<< fHistCentBin
 	 	<< "END";  // histo of < vn * vn > for [ih][ik][ihh][ikk][iCent] 
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
		if( inputCent >= fCentBin[iCbin] && inputCent < fCentBin[iCbin+1] ){fCBin = iCbin;};
	}
	int trk_number = fInputList->GetEntriesFast();
	fh_ntracks[fCBin]->Fill( trk_number ) ;
	fh_cent->Fill( inputCent ) ; 



	enum{kSub1, kSub2, kNSub};
	enum{kReal, kImg, kNPhase};
	double Eta_config[kNSub][2];
	Eta_config[kSub1][0] = fEta_min;
	Eta_config[kSub1][1] = fEta_max;
	Eta_config[kSub2][0] = -1*fEta_max;
	Eta_config[kSub2][1] = -1*fEta_min; 
	
	double Qn[kNH][kNSub][kNPhase]; // Qn for 2nd harmonics

	//---------------- Do initialize here -----------
	for(int ih=2; ih<kNH; ih++){
			for(int isub=0; isub<kNSub; isub++){
					for(int iphase=0; iphase < kNPhase; iphase++){
							Qn[ih][isub][iphase] =0; // init
					}
			}
	}
	//--------------- Calculate Qn--------------------
	for(int ih=2; ih<kNH; ih++){
			for(int isub=0; isub<kNSub; isub++){
					Qn[ih][isub][kReal] = Get_Qn_Real( Eta_config[isub][0], Eta_config[isub][1], ih );
					Qn[ih][isub][kImg]  = Get_Qn_Img(  Eta_config[isub][0], Eta_config[isub][1], ih );
			}
	}

	//-------------- Fill histos with below Values ----
	double vn[kNH][nKL]; 
	double vn_vn[kNH][nKL][kNH][nKL];
	//initiation
	for(int ih=0; ih<kNH; ih++){
		for(int ik=0; ik<nKL ; ik++){
				vn[ih][ik] = -999;
			for(int ihh=0; ihh<kNH; ihh++){
				for(int ikk=0; ikk<nKL; ikk++){
					vn_vn[ih][ik][ihh][ikk] = -999;
				}
			}
		}
	}


	// calculate hvn	
	for( int ih=2; ih<kNH; ih++){
		for(int ik=1; ik<nKL; ik++){
			if(ik==1) vn[ih][ik] = Get_QC_Vn( Qn[ih][kSub1][kReal], Qn[ih][kSub1][kImg], Qn[ih][kSub2][kReal], Qn[ih][kSub2][kImg]);
			else{
				vn[ih][ik] = TMath::Power( vn[ih][1] , ik ); 
			}
		}		
	}
	// vn^k calcualted for n.... k....

	// calculate hvn_vn (2 combination of vn) 
	for( int ih=2; ih<kNH; ih++){ 
		for( int ik=1; ik<nKL; ik++){
			for( int ihh=2; ihh<kNH; ihh++){ 
				for(int ikk=1; ikk<nKL; ikk++){
					vn_vn[ih][ik][ihh][ikk] = vn[ih][ik] * vn[ihh][ikk]; 
				}
			}
		}
	}
	//************************************************************************
	// doing this 
	// check it when ik = 2, ikk=2  DO NOT BELIEVE LINE 290~309 (may be wrong )
	//
	for(int ih=2; ih<kNH; ih++){
				double Qn_Qn_real = Complex_product_real( Qn[ih][kSub1][kReal], Qn[ih][kSub1][kImg], Qn[ih][kSub2][kReal], -1*Qn[ih][kSub2][kImg]) ;
				double Qn_Qn_img = Complex_product_img( Qn[ih][kSub1][kReal], Qn[ih][kSub1][kImg], Qn[ih][kSub2][kReal], -1*Qn[ih][kSub2][kImg]) ;

				//vn[ih][2] = Complex_abs( Qn_Qn_real, Qn_Qn_img );
				vn[ih][2] = Qn_Qn_real;
				//vn[ih][4] = Complex_abs( Complex_sqr_real(Qn_Qn_real, Qn_Qn_img),Complex_sqr_img(Qn_Qn_real, Qn_Qn_img));
				vn[ih][4] = Complex_sqr_real(Qn_Qn_real, Qn_Qn_img) ;

	}

	for( int ih=2; ih<kNH; ih++){ 
			for( int ihh=2; ihh<kNH; ihh++){  // <v2v3> = < Q2 Q2* Q3 Q3* >
					double Qn_Qn_real = Complex_product_real( Qn[ih][kSub1][kReal], Qn[ih][kSub1][kImg], Qn[ih][kSub2][kReal], -1* Qn[ih][kSub2][kImg]) ;
					double Qn_Qn_img = Complex_product_img( Qn[ih][kSub1][kReal], Qn[ih][kSub1][kImg], Qn[ih][kSub2][kReal], -1* Qn[ih][kSub2][kImg]) ;
					double Qm_Qm_real = Complex_product_real( Qn[ihh][kSub1][kReal], Qn[ihh][kSub1][kImg], Qn[ihh][kSub2][kReal], -1* Qn[ihh][kSub2][kImg]) ;
					double Qm_Qm_img = Complex_product_img( Qn[ihh][kSub1][kReal], Qn[ihh][kSub1][kImg], Qn[ihh][kSub2][kReal], -1 * Qn[ihh][kSub2][kImg]) ;

					double QnQm_real = Complex_product_real( Qn_Qn_real, Qn_Qn_img, Qm_Qm_real, Qm_Qm_img ) ;
					double QnQm_img = Complex_product_img( Qn_Qn_real, Qn_Qn_img, Qm_Qm_real, Qm_Qm_img ) ;
					//vn_vn[ih][2][ihh][2] =  Complex_abs( QnQm_real, QnQm_img) ; 
					vn_vn[ih][2][ihh][2] =  QnQm_real ; 

				}
			}
	//************************************************************************

	//Fill the Histos here
	for( int ih=2; ih<kNH; ih++){ 
		for( int ik=1; ik<nKL; ik++){
			fh_vn[ih][ik][fCBin]->Fill( vn[ih][ik]  ); // Fill hvn
			for( int ihh=2; ihh<kNH; ihh++){ 
				for(int ikk=1; ikk<nKL; ikk++){
					 fh_vn_vn[ih][ik][ihh][ikk][fCBin]->Fill( vn_vn[ih][ik][ihh][ikk] ) ; // Fill hvn_vn 
				}
			}
		}
	}
	///	

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
double AliJFFlucAnalysis::Get_Qn_Img(double eta1, double eta2, int harmonics)
{
		int nh = harmonics;
		double Qn_img = 0;
		double Sub_Ntrk =0;

		Long64_t ntracks = fInputList->GetEntriesFast();
		for( Long64_t it=0; it< ntracks; it++){
			AliJBaseTrack *itrack = (AliJBaseTrack*)fInputList->At(it); // load track
			double eta = itrack->Eta();
			double pt = itrack->Pt();
			fh_eta[fCBin]->Fill(eta);
			fh_pt[fCBin]->Fill(pt);
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







