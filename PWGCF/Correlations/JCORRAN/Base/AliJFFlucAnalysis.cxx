#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF3.h>
#include <TMath.h>
#include <TComplex.h>
#include <TClonesArray.h>
#include "AliJBaseTrack.h"
#include "AliJFFlucAnalysis.h"
#pragma GCC diagnostic warning "-Wall"

//ClassImp(AliJFFlucAnalysis)

		//________________________________________________________________________
AliJFFlucAnalysis::AliJFFlucAnalysis() :
	//: AliAnalysisTaskSE(),
	fInputList(0),
	fVertex(0),
	fCent(0),
	fCBin(0),
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
	fh_phieta(),
	fh_phietaz(),
	//fh_Qvector(),
	fh_ntracks(),
	fh_vn(),
	fh_vna(),
	fh_vn_vn()
	/*fh_cn_4c(),
	fh_cn_2c(),
	fh_cn_cn_2c(),
	fh_cn_2c_eta10(),
	fh_cn_cn_2c_eta10()*/
{
	subeventMask = SUBEVENT_A|SUBEVENT_B;
	binning = BINNING_CENT_PbPb;
	flags = 0;
	fEta_min = 0;
	fEta_max = 0;
	fQC_eta_cut_min = -0.8; // default setting
	fQC_eta_cut_max = 0.8; // default setting
	fQC_eta_gap_half = 0.4;
	fImpactParameter = -1;
}

//________________________________________________________________________
AliJFFlucAnalysis::AliJFFlucAnalysis(const char *name) :
	//: AliAnalysisTaskSE(name),
	fInputList(0),
	fVertex(0),
	fCent(0),
	fCBin(0),
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
	fh_phieta(),
	fh_phietaz(),
	//fh_Qvector(),
	fh_ntracks(),
	fh_vn(),
	fh_vna(),
	fh_vn_vn()
	/*fh_cn_4c(),
	fh_cn_2c(),
	fh_cn_cn_2c(),
	fh_cn_2c_eta10(),
	fh_cn_cn_2c_eta10()*/
{
	cout << "analysis task created " << endl;

	subeventMask = SUBEVENT_A|SUBEVENT_B;
	binning = BINNING_CENT_PbPb;
	flags = 0;
	fEta_min = 0;
	fEta_max = 0;
	fQC_eta_cut_min = -0.8; // default setting
	fQC_eta_cut_max = 0.8; // default setting
	fQC_eta_gap_half = 0.4;
	fImpactParameter = -1;
}

//Double_t AliJFFlucAnalysis::CentBin[8] = {0, 5, 10, 20, 30, 40, 50, 60};
//Double_t AliJFFlucAnalysis::CentBin[CENTN_NAT+1] = {0, 1, 2, 5, 10, 20, 30, 40, 50, 60};
//UInt_t AliJFFlucAnalysis::NCentBin = sizeof(AliJFFlucAnalysis::CentBin)/sizeof(AliJFFlucAnalysis::CentBin[0])-1;
//Double_t AliJFFlucAnalysis::MultBin[MULTN+1] = {30.651,31.125,42.627,43.263,54.598,55.399,66.613,67.581,78.665,79.802,90.789,92.098,102.805,104.286,114.903,116.569};
//! Max bins 96
Double_t AliJFFlucAnalysis::CentBin_PbPb_default[][2] = {{0,1},{1,2},{2,5},{5,10},{10,20},{20,30},{30,40},{40,50},{50,60}};
Double_t AliJFFlucAnalysis::MultBin_PbPb_1[][2] = {{25.245,25.611},{31.555,31.982},{37.898,38.387},{44.251,44.803},{50.584,51.198},{56.942,57.62},{66.341,67.113},{79.03,79.93},{91.764,92.792},{104.431,105.587},{117.103,118.387},{136.063,137.541},{161.3,163.033},{186.707,188.698},{212.126,214.374},{281.085,284.042},{344.727,348.332},{408.197,412.447},{471.798,476.696},{535.041,540.583},{598.475,604.663},{662.018,668.854},{725.484,732.967},{788.636,796.763},{852.345,861.121},{915.632,925.054},{979.134,989.204},{1042.376,1053.091},{1105.751,1117.114},{1169.205,1181.215},{1232.423,1245.079},{1295.89,1309.194},{1359.578,1373.533},{1422.675,1437.274},{1485.99,1501.236},{1549.608,1565.503},{1612.809,1629.351},{1676.158,1693.347},{1739.647,1757.486},{1802.865,1821.351},{1866.328,1885.462},{1929.777,1949.561},{1993.119,2013.553},{2056.422,2077.504},{2119.716,2141.448},{2182.935,2205.317},{2246.261,2269.296},{2309.523,2333.209},{2373.038,2397.381},{2436.146,2461.145},{2499.398,2525.057},{2563.07,2589.39},{2626.482,2653.462},{2689.435,2717.075},{2753.073,2781.382},{2815.911,2844.889},{2879.821,2909.486},{2943.2,2973.535},{3006.237,3037.292},{3068.718,3100.47},{3132.669,3165.259}};
Double_t AliJFFlucAnalysis::MultBin_pPb_1[][2] = {{30.651,31.125},{42.627,43.263},{54.598,55.399},{66.613,67.581},{78.665,79.802},{90.789,92.098},{102.805,104.286},{114.903,116.569}};
Double_t (*AliJFFlucAnalysis::pBin[3])[2] = {&CentBin_PbPb_default[0],&MultBin_PbPb_1[0],&MultBin_pPb_1[0]};
UInt_t AliJFFlucAnalysis::NBin[3] = {
	sizeof(AliJFFlucAnalysis::CentBin_PbPb_default)/sizeof(AliJFFlucAnalysis::CentBin_PbPb_default[0]),
	sizeof(AliJFFlucAnalysis::MultBin_PbPb_1)/sizeof(AliJFFlucAnalysis::MultBin_PbPb_1[0]),
	sizeof(AliJFFlucAnalysis::MultBin_pPb_1)/sizeof(AliJFFlucAnalysis::MultBin_pPb_1[0])
};
//UInt_t AliJFFlucAnalysis::CentralityTranslationMap[CENTN_NAT] = {0,0,0,1,2,3,4,5,6};
Double_t AliJFFlucAnalysis::pttJacek[74] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 45, 50, 60, 70, 80, 90, 100};
UInt_t AliJFFlucAnalysis::NpttJacek = sizeof(AliJFFlucAnalysis::pttJacek)/sizeof(AliJFFlucAnalysis::pttJacek[0])-1;

//________________________________________________________________________
AliJFFlucAnalysis::AliJFFlucAnalysis(const AliJFFlucAnalysis& a):
	//AliAnalysisTaskSE(a.GetName()),
	fInputList(a.fInputList),
	fVertex(a.fVertex),
	fCent(a.fCent),
	fCBin(a.fCBin),
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
	fh_phieta(a.fh_phieta),
	fh_phietaz(a.fh_phietaz),
	//fh_Qvector(a.fh_Qvector),
	fh_ntracks(a.fh_ntracks),
	fh_vn(a.fh_vn),
	fh_vna(a.fh_vna),
	fh_vn_vn(a.fh_vn_vn)
	/*fh_cn_4c(a.fh_cn_4c),
	fh_cn_2c(a.fh_cn_2c),
	fh_cn_cn_2c(a.fh_cn_cn_2c),
	fh_cn_2c_eta10(a.fh_cn_2c_eta10),
	fh_cn_cn_2c_eta10(a.fh_cn_cn_2c_eta10)*/
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
	fHMG = new AliJHistManager("AliJFFlucHistManager","jfluc");
	// set AliJBin here //
	fBin_Subset .Set("Sub","Sub","Sub:%d", AliJBin::kSingle).SetBin(2);
	fBin_h .Set("NH","NH","NH:%d", AliJBin::kSingle).SetBin(kNH);
	fBin_k .Set("K","K","K:%d", AliJBin::kSingle).SetBin(nKL);

	fBin_hh .Set("NHH","NHH","NHH:%d", AliJBin::kSingle).SetBin(kcNH);
	fBin_kk .Set("KK","KK","KK:%d", AliJBin::kSingle).SetBin(nKL);

	//TODO: index with binning the array of pointers
	if(binning != BINNING_CENT_PbPb)
		fHistCentBin.Set("MultBin","MultBin","Cent:%d",AliJBin::kSingle).SetBin(NBin[binning]);
	else fHistCentBin.Set("CentBin","CentBin","Cent:%d",AliJBin::kSingle).SetBin(NBin[binning]);

	fVertexBin .Set("Vtx","Vtx","Vtx:%d", AliJBin::kSingle).SetBin(3);
	fCorrBin .Set("C", "C","C:%d", AliJBin::kSingle).SetBin(28);

	fBin_Nptbins .Set("PtBin","PtBin", "Pt:%d", AliJBin::kSingle).SetBin(N_ptbins);

	// set AliJTH1D here //
	fh_cent
		<< TH1D("h_cent","h_cent", 200, 0, 100)
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
		<< TH1D("hChargedPtJacek", "", AliJFFlucAnalysis::NpttJacek, AliJFFlucAnalysis::pttJacek)
		<< fHistCentBin
		<< "END" ;

	fh_eta
		<< TH1D("h_eta", "h_eta", 40, -2.0, 2.0 )
		<< fHistCentBin
		<< "END" ;
	fh_phi
		<< TH1D("h_phi", "h_phi", 50,-TMath::Pi(),TMath::Pi())
		<< fHistCentBin << fBin_Subset
		<< "END" ;
	fh_phieta
		<< TH2D("h_phieta","h_phieta",50,-TMath::Pi(),TMath::Pi(),40,-2.0,2.0)
		<< fHistCentBin
		<< "END";
	fh_phietaz
		<< TH3D("h_phietaz","h_phietaz",50,-TMath::Pi(),TMath::Pi(),40,-2.0,2.0,20,-10.0,10.0)
		<< fHistCentBin
		<< "END";
	/*fh_Qvector
		<< TH1D("h_QVector", "h_QVector", 100, -10, 10)
		<< fHistCentBin << fBin_Subset
		<< fBin_h
		<< "END" ;*/

	fh_psi_n
		<< TH1D("h_psi","h_psi",50,-0.5*TMath::Pi(),0.5*TMath::Pi())
		<< fBin_h
		<< fHistCentBin
		<< "END" ;
	
	fh_cos_n_phi
		<< TH1D("h_cos_n_phi","h_cos_n_phi",50,-1,1)
		<< fBin_h
		<< fHistCentBin
		<< "END" ;

	fh_sin_n_phi
		<< TH1D("h_sin_n_phi","h_sin_n_phi",50,-1,1)
		<< fBin_h
		<< fHistCentBin
		<< "END" ;

	fh_cos_n_psi_n
		<< TH1D("h_cos_n_psi","h_cos_n_psi",50,-1,1)
		<< fBin_h
		<< fHistCentBin
		<< "END" ;

	fh_sin_n_psi_n
		<< TH1D("h_sin_n_psi","h_sin_n_psi",50,-1,1)
		<< fBin_h
		<< fHistCentBin
		<< "END" ;

	fh_ntracks
		<< TH1D("h_tracks", "h_tracks", 100, 0, 5000)
		<< fHistCentBin
		<< "END" ;

	fh_vn
		<< TH1D("hvn","hvn", 1024, -1.0, 1.0)
		<< fBin_h << fBin_k
		<< fHistCentBin
		<< "END";   // histogram of vn_h^k values for [ih][ik][iCent]
	fh_vna
		<< TH1D("hvna","hvna", 1024, -1.0, 1.0)
		<< fBin_h << fBin_k
		<< fHistCentBin
		<< "END";   // histogram of vn_h^k values for [ih][ik][iCent]
	fh_vn_vn
		<< TH1D("hvn_vn", "hvn_vn", 1024, -1.0, 1.0)
		<< fBin_h << fBin_k
		<< fBin_hh << fBin_kk
		<< fHistCentBin
		<< "END";  // histo of < vn * vn > for [ih][ik][ihh][ikk][iCent]
	/*fh_cn_4c
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
		<< "END";*/
	fh_correlator
		<< TH1D("h_corr", "h_corr", 1024, -3.0, 3.0)
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
	delete fHMG;
}

#define A i
#define B (1-i)
#define C(u) TComplex::Conjugate(u)
//TODO: conjugate macro
inline TComplex TwoGap(const TComplex (*pQq)[AliJFFlucAnalysis::kNH][AliJFFlucAnalysis::nKL], uint i, uint a, uint b){
	return pQq[A][a][1]*C(pQq[B][b][1]);
}

inline TComplex ThreeGap(const TComplex (*pQq)[AliJFFlucAnalysis::kNH][AliJFFlucAnalysis::nKL], uint i, uint a, uint b, uint c){
	return pQq[A][a][1]*C(pQq[B][b][1]*pQq[B][c][1]-pQq[B][b+c][2]);
}

inline TComplex FourGap22(const TComplex (*pQq)[AliJFFlucAnalysis::kNH][AliJFFlucAnalysis::nKL], uint i, uint a, uint b, uint c, uint d){
	return pQq[A][a][1]*pQq[A][b][1]*C(pQq[B][c][1]*pQq[B][d][1])-pQq[A][a+b][2]*C(pQq[B][c][1]*pQq[B][d][1])-pQq[A][a][1]*pQq[A][b][1]*C(pQq[B][c+d][2])+pQq[A][a+b][2]*C(pQq[B][c+d][2]);
}

inline TComplex FourGap13(const TComplex (*pQq)[AliJFFlucAnalysis::kNH][AliJFFlucAnalysis::nKL], uint i, uint a, uint b, uint c, uint d){
	return pQq[A][a][1]*C(pQq[B][b][1]*pQq[B][c][1]*pQq[B][d][1]-pQq[B][b+c][2]*pQq[B][d][1]-pQq[B][b+d][2]*pQq[B][c][1]-pQq[B][c+d][2]*pQq[B][b][1]+2.0*pQq[B][b+c+d][3]);
}

inline TComplex SixGap33(const TComplex (*pQq)[AliJFFlucAnalysis::kNH][AliJFFlucAnalysis::nKL], uint i, uint n1, uint n2, uint n3, uint n4, uint n5, uint n6){
	return pQq[A][n1][1]*pQq[A][n2][1]*pQq[A][n3][1]*C(pQq[B][n4][1]*pQq[B][n5][1]*pQq[B][n6][1])-pQq[A][n1][1]*pQq[A][n2][1]*pQq[A][n3][1]*C(pQq[B][n4+n5][2]*pQq[B][n6][1])-pQq[A][n1][1]*pQq[A][n2][1]*pQq[A][n3][1]*C(pQq[B][n4+n6][2]*pQq[B][n5][1])-pQq[A][n1][1]*pQq[A][n2][1]*pQq[A][n3][1]*C(pQq[B][n5+n6][2]*pQq[B][n4][1])+2.0*pQq[A][n1][1]*pQq[A][n2][1]*pQq[A][n3][1]*C(pQq[B][n4+n5+n6][3])-pQq[A][n1+n2][2]*pQq[A][n3][1]*C(pQq[B][n4][1]*pQq[B][n5][1]*pQq[B][n6][1])+pQq[A][n1+n2][2]*pQq[A][n3][1]*C(pQq[B][n4+n5][2]*pQq[B][n6][1])+pQq[A][n1+n2][2]*pQq[A][n3][1]*C(pQq[B][n4+n6][2]*pQq[B][n5][1])+pQq[A][n1+n2][2]*pQq[A][n3][1]*C(pQq[B][n5+n6][2]*pQq[B][n4][1])-2.0*pQq[A][n1+n2][2]*pQq[A][n3][1]*C(pQq[B][n4+n5+n6][3])-pQq[A][n1+n3][2]*pQq[A][n2][1]*C(pQq[B][n4][1]*pQq[B][n5][1]*pQq[B][n6][1])+pQq[A][n1+n3][2]*pQq[A][n2][1]*C(pQq[B][n4+n5][2]*pQq[B][n6][1])+pQq[A][n1+n3][2]*pQq[A][n2][1]*C(pQq[B][n4+n6][2]*pQq[B][n5][1])+pQq[A][n1+n3][2]*pQq[A][n2][1]*C(pQq[B][n5+n6][2]*pQq[B][n4][1])-2.0*pQq[A][n1+n3][2]*pQq[A][n2][1]*C(pQq[B][n4+n5+n6][3])-pQq[A][n2+n3][2]*pQq[A][n1][1]*C(pQq[B][n4][1]*pQq[B][n5][1]*pQq[B][n6][1])+pQq[A][n2+n3][2]*pQq[A][n1][1]*C(pQq[B][n4+n5][2]*pQq[B][n6][1])+pQq[A][n2+n3][2]*pQq[A][n1][1]*C(pQq[B][n4+n6][2]*pQq[B][n5][1])+pQq[A][n2+n3][2]*pQq[A][n1][1]*C(pQq[B][n5+n6][2]*pQq[B][n4][1])-2.0*pQq[A][n2+n3][2]*pQq[A][n1][1]*C(pQq[B][n4+n5+n6][3])+2.0*pQq[A][n1+n2+n3][3]*C(pQq[B][n4][1]*pQq[B][n5][1]*pQq[B][n6][1])-2.0*pQq[A][n1+n2+n3][3]*C(pQq[B][n4+n5][2]*pQq[B][n6][1])-2.0*pQq[A][n1+n2+n3][3]*C(pQq[B][n4+n6][2]*pQq[B][n5][1])-2.0*pQq[A][n1+n2+n3][3]*C(pQq[B][n5+n6][2]*pQq[B][n4][1])+4.0*pQq[A][n1+n2+n3][3]*C(pQq[B][n4+n5+n6][3]);
}
#undef C

//________________________________________________________________________
void AliJFFlucAnalysis::UserExec(Option_t *) {
	// find Centrality
	int trk_number = fInputList->GetEntriesFast();
	fCBin = (binning != BINNING_CENT_PbPb)?
		GetBin((double)trk_number,binning):
		GetBin(fCent,binning); //--- similarly in Task

	if(fCBin == -1)
		return;
	
	fh_ntracks[fCBin]->Fill( trk_number ) ;
	fh_cent->Fill(fCent) ;
	fh_ImpactParameter->Fill( fImpactParameter);
	Fill_QA_plot( fEta_min, fEta_max );

	enum{kSubA, kSubB, kNSub};
	enum{kMin, kMax};
	Double_t Eta_config[kNSub][2];
	Eta_config[kSubA][kMin] = fEta_min;  // 0.4 min for SubA
	Eta_config[kSubA][kMax] = fEta_max;  // 0.8 max for SubA
	Eta_config[kSubB][kMin] = -fEta_max; // -0.8  min for SubB
	Eta_config[kSubB][kMax] = -fEta_min; // -0.4  max for SubB

	// use complex variable instead of doulbe Qn //
	/*TComplex QnA[kNH];
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
	NSubTracks[kSubB] = QnB[0].Re(); // this is number of tracks in Sub B*/
	
	CalculateQvectorsQC(fEta_min,fEta_max);

	for(int ih=2; ih<kNH; ih++){
		fh_cos_n_phi[ih][fCBin]->Fill(QvectorQC[ih][1].Re()/QvectorQC[0][1].Re());
		fh_sin_n_phi[ih][fCBin]->Fill(QvectorQC[ih][1].Im()/QvectorQC[0][1].Re());
		//
		//
		Double_t psi = QvectorQC[ih][1].Theta();
		fh_psi_n[ih][fCBin]->Fill(psi);
		fh_cos_n_psi_n[ih][fCBin]->Fill(TMath::Cos((Double_t)ih*psi));
		fh_sin_n_psi_n[ih][fCBin]->Fill(TMath::Sin((Double_t)ih*psi));
	}

	// v2^2 :  k=1  /// remember QnQn = vn^(2k) not k
	// use k=0 for check v2, v3 only
	Double_t vn2[kNH][nKL];
	Double_t vn2_vn2[kNH][nKL][kNH][nKL];

	TComplex corr[kNH][nKL];
	TComplex ncorr[kNH][nKL];
	TComplex ncorr2[kNH][nKL][kcNH][nKL];

	const TComplex (*pQq)[kNH][nKL] = QvectorQCeta10;

	for(int i = 0; i < 2; ++i){
		if((subeventMask & (1<<i)) == 0)
			continue;
		//Double_t ref_2p = N[i][0]*N[i][1];//TwoGap(pQq,i,0,0).Re();
		Double_t ref_2p = TwoGap(pQq,i,0,0).Re();
		Double_t ref_3p = ThreeGap(pQq,i,0,0,0).Re();
		Double_t ref_4p = FourGap22(pQq,i,0,0,0,0).Re();
		Double_t ref_4pB = FourGap13(pQq,i,0,0,0,0).Re();
		Double_t ref_6p = SixGap33(pQq,i,0,0,0,0,0,0).Re();

		Double_t ebe_2p_weight = 1.0;
		Double_t ebe_3p_weight = 1.0;
		Double_t ebe_4p_weight = 1.0;
		Double_t ebe_4p_weightB = 1.0;
		Double_t ebe_6p_weight = 1.0;
		if(flags & FLUC_EBE_WEIGHTING){
			ebe_2p_weight = ref_2p;//N[i][0]*N[i][1];//NSubTracks[kSubA] * NSubTracks[kSubB] ;
			ebe_3p_weight = ref_3p;//ebe_2p_weight*(N[i][1]-1.0);// * (NSubTracks[kSubB]-1.0);
			ebe_4p_weight = ref_4p;
			ebe_4p_weightB = ref_4pB;//ebe_3p_weight*(N[i][1]-2.0);// * (NSubTracks[kSubB]-2.0);
			ebe_6p_weight = ref_6p;
		}
		Double_t ref_2Np[2*nKL] = {
			ref_2p,
			ref_4p,
			ref_6p
		};
		Double_t ebe_2Np_weight[2*nKL] = {
			ebe_2p_weight,
			ebe_4p_weight,
			ebe_6p_weight
		};
		if(flags & FLUC_EBE_WEIGHTING){
			for(int ik=3; ik<2*nKL; ik++){
				double dk = (double)ik;
				ref_2Np[ik] = ref_2Np[ik-1]*max(pQq[A][0][1].Re()-dk,1.0)*max(pQq[B][0][1].Re()-dk,1.0);
				ebe_2Np_weight[ik] = ebe_2Np_weight[ik-1]*max(pQq[A][0][1].Re()-dk,1.0)*max(pQq[B][0][1].Re()-dk,1.0);
			}
		}else for(int ik=3; ik<2*nKL; ik++){
			double dk = (double)ik;
			ref_2Np[ik] = ref_2Np[ik-1]*max(pQq[A][0][1].Re()-dk,1.0)*max(pQq[B][0][1].Re()-dk,1.0);
			ebe_2Np_weight[ik] = 1.0;
		}

		for(int ih=2; ih<kNH; ih++){
			//corr[ih][1] = pQn[i][0][ih]*pQn[i][1][ih]*N[i][0]*N[i][1];//QnA[ih]*QnB_star[ih];
			corr[ih][1] = TwoGap(pQq,i,ih,ih);
			for(int ik=2; ik<nKL; ik++)
				corr[ih][ik] = corr[ih][ik-1]*corr[ih][1];//TComplex::Power(corr[ih][1],ik);
			ncorr[ih][1] = corr[ih][1];
			ncorr[ih][2] = FourGap22(pQq,i,ih,ih,ih,ih);//mf*(corr[ih][2]*N[i][0]*N[i][1]-pQn[i][1][2*ih]*pQn[i][0][ih]*pQn[i][0][ih]*N[i][0]-pQn[i][0][2*ih]*pQn[i][1][ih]*pQn[i][1][ih]*N[i][1]+pQn[i][1][2*ih]*pQn[i][0][2*ih]);
			ncorr[ih][3] = SixGap33(pQq,i,ih,ih,ih,ih,ih,ih);
			for(int ik=4; ik<nKL; ik++)
				ncorr[ih][ik] = corr[ih][ik]; //for 8,...-particle correlations, ignore the autocorrelation / weight dependency for now

			for(int ihh=2; ihh<kcNH; ihh++){
				ncorr2[ih][1][ihh][1] = FourGap22(pQq,i,ih,ihh,ih,ihh);
				ncorr2[ih][1][ihh][2] = SixGap33(pQq,i,ih,ihh,ihh,ih,ihh,ihh);
				ncorr2[ih][2][ihh][1] = SixGap33(pQq,i,ih,ih,ihh,ih,ih,ihh);
				for(int ik=2; ik<nKL; ik++)
					for(int ikk=2; ikk<nKL; ikk++)
						ncorr2[ih][ik][ihh][ikk] = ncorr[ih][ik]*ncorr[ihh][ikk];
			}
		}
		
		for(int ih=2; ih<kNH; ih++){
			for(int ik=1; ik<nKL; ik++){ // 2k(0) =1, 2k(1) =2, 2k(2)=4....
				vn2[ih][ik] = corr[ih][ik].Re()/ref_2Np[ik-1];
				fh_vn[ih][ik][fCBin]->Fill(vn2[ih][ik],ebe_2Np_weight[ik-1]);
				fh_vna[ih][ik][fCBin]->Fill(ncorr[ih][ik].Re()/ref_2Np[ik-1],ebe_2Np_weight[ik-1]);
				for(int ihh=2; ihh<kcNH; ihh++){
					for(int ikk=1; ikk<nKL; ikk++){
						vn2_vn2[ih][ik][ihh][ikk] = ncorr2[ih][ik][ihh][ikk]/ref_2Np[ik+ikk-1];//(ncorr[ih][ik]*ncorr[ihh][ikk]).Re()/ref_2Np[ik+ikk-1];
						fh_vn_vn[ih][ik][ihh][ikk][fCBin]->Fill(vn2_vn2[ih][ik][ihh][ikk],ebe_2Np_weight[ik+ikk-1]); // Fill hvn_vn
					}
				}
			}
			fSingleVn[ih][0] = TMath::Sqrt(vn2[ih][1]); // fill single vn with SP as method 0
		}

		//************************************************************************
		TComplex V4V2star_2 = pQq[A][4][1] * pQq[B][2][1] * pQq[B][2][1];
		TComplex V4V2starv2_2 =	V4V2star_2 * corr[2][1]/ref_2Np[0];//vn[2][1]
		TComplex V4V2starv2_4 = V4V2star_2 * corr[2][2]/ref_2Np[1];//vn2[2][2]
		TComplex V5V2starV3starv2_2 = pQq[A][5][1] * pQq[B][2][1] * pQq[B][3][1] * corr[2][1]/ref_2Np[0]; //vn2[2][1]
		TComplex V5V2starV3star = pQq[A][5][1] * pQq[B][2][1] * pQq[B][3][1];
		TComplex V5V2starV3startv3_2 = V5V2starV3star * corr[3][1]/ref_2Np[0]; //vn2[3][1]
		TComplex V6V2star_3 = pQq[A][6][1] * pQq[B][2][1] * pQq[B][2][1] * pQq[B][2][1];
		TComplex V6V3star_2 = pQq[A][6][1] * pQq[B][3][1] * pQq[B][3][1];
		TComplex V6V2starV4star = pQq[A][6][1] * pQq[B][2][1] * pQq[B][4][1];
		TComplex V7V2star_2V3star = pQq[A][7][1] * pQq[B][2][1] * pQq[B][2][1] * pQq[B][3][1];
		TComplex V7V2starV5star = pQq[A][7][1] * pQq[B][2][1] * pQq[B][5][1];
		TComplex V7V3starV4star = pQq[A][7][1] * pQq[B][3][1] * pQq[B][4][1];
		TComplex V8V2starV3star_2 = pQq[A][8][1] * pQq[B][2][1] * pQq[B][3][1] * pQq[B][3][1];
		TComplex V8V2star_4 = pQq[A][8][1] * TComplex::Power(pQq[B][2][1],4);

		// New correlators (Modified by You's correction term for self-correlations)
		//double nf = 1.0/(N[i][1]-1.0);
		///double ef = nf/(N[i][1]-2.0);
		TComplex nV4V2star_2 = ThreeGap(pQq,i,4,2,2)/ref_3p;//V4V2_star2-pQq[i][A][4][1]*pQq[i][B][4][2]; //nf*( V4V2star_2*N[i][1] - pQn[i][0][4]*pQn[i][1][4] );
		//TComplex nV4V2star_2 = nf*( pQn[i][0][4]*pQn[i][1][2]*pQn[i][1][2]*N[i][1] - pQn[i][0][4]*pQn[i][1][4] );//nf*( V4V2star_2*N[i][1] - pQn[i][0][4]*pQn[i][1][4] );
		TComplex nV5V2starV3star = ThreeGap(pQq,i,5,2,3)/ref_3p;//V5V2starV3star-pQq[i][A][5][1]*pQq[i][B][5][2]; //nf*( V5V2starV3star*N[i][1] - pQn[i][0][5]*pQn[i][1][5] );
		TComplex nV6V2star_3 = FourGap13(pQq,i,6,2,2,2)/ref_4pB;//V6V2star_3-pQq[i][A][6][2]*pQq[i][B][4][2]*pQq[i][B][2][1]- //pQn[i][0][6]*ef*( pQn[i][1][2]*pQn[i][1][2]*pQn[i][1][2]*N[i][1]*N[i][1] - 3.0*pQn[i][1][2]*pQn[i][1][4]*N[i][1] + 2.0*pQn[i][1][6] );
		TComplex nV6V3star_2 = ThreeGap(pQq,i,6,3,3)/ref_3p;//nf*(V6V3star_2*N[i][1] - pQn[i][0][6]*pQn[i][1][6]);
		TComplex nV6V2starV4star = ThreeGap(pQq,i,6,2,4)/ref_3p;//nf*(V6V2starV4star*N[i][1] - pQn[i][0][6]*pQn[i][1][6]);
		TComplex nV7V2star_2V3star = FourGap13(pQq,i,7,2,2,3)/ref_4pB;//pQn[i][0][7]*ef*( pQn[i][1][2]*pQn[i][1][2]*pQn[i][1][3]*N[i][1]*N[i][1] - 2.0*pQn[i][1][2]*pQn[i][1][5]*N[i][1] - pQn[i][1][3]*pQn[i][1][4]*N[i][1] + 2.0*pQn[i][1][7] );
		TComplex nV7V2starV5star = ThreeGap(pQq,i,7,2,5)/ref_3p;//nf*(V7V2starV5star*N[i][1] - pQn[i][0][7]*pQn[i][1][7]);
		TComplex nV7V3starV4star = ThreeGap(pQq,i,7,3,4)/ref_3p;//nf*(V7V3starV4star*N[i][1] - pQn[i][0][7]*pQn[i][1][7]);
		TComplex nV8V2starV3star_2 = FourGap13(pQq,i,8,2,3,3)/ref_4pB;//pQn[i][0][8]*ef*( pQn[i][1][2]*pQn[i][1][3]*pQn[i][1][3]*N[i][1]*N[i][1] - 2.0*pQn[i][1][3]*pQn[i][1][5]*N[i][1] - pQn[i][1][2]*pQn[i][1][6]*N[i][1] + 2.0*pQn[i][1][8] );

		TComplex nV4V4V2V2 = FourGap22(pQq,i,4,2,4,2)/ref_4p;//(pQn[i][0][4]*pQn[i][1][4]*pQn[i][0][2]*pQn[i][1][2]) - ((1/(N[i][1]-1) * pQn[i][1][6] * pQn[i][0][4] *pQn[i][0][2] ))
			//- ((1/(N[i][0]-1) * pQn[i][0][6]*pQn[i][1][4] * pQn[i][1][2])) + (1/((N[i][0]-1)*(N[i][1]-1))*pQn[i][0][6]*pQn[i][1][6] );
		TComplex nV3V3V2V2 = FourGap22(pQq,i,3,2,3,2)/ref_4p;//(pQn[i][0][3]*pQn[i][1][3]*pQn[i][0][2]*pQn[i][1][2]) - ((1/(N[i][1]-1) * pQn[i][1][5] * pQn[i][0][3] *pQn[i][0][2] ))
			//- ((1/(N[i][0]-1) * pQn[i][0][5]*pQn[i][1][3] * pQn[i][1][2])) + (1/((N[i][0]-1)*(N[i][1]-1))*pQn[i][0][5]*pQn[i][1][5] );
		TComplex nV5V5V2V2 = FourGap22(pQq,i,5,2,5,2)/ref_4p;//(pQn[i][0][5]*pQn[i][1][5]*pQn[i][0][2]*pQn[i][1][2]) - ((1/(N[i][1]-1) * pQn[i][1][7] * pQn[i][0][5] *pQn[i][0][2] ))
			//- ((1/(N[i][0]-1) * pQn[i][0][7]*pQn[i][1][5] * pQn[i][1][2])) + (1/((N[i][0]-1)*(N[i][1]-1))*pQn[i][0][7]*pQn[i][1][7] );
		TComplex nV5V5V3V3 = FourGap22(pQq,i,5,3,5,3)/ref_4p;//(pQn[i][0][5]*pQn[i][1][5]*pQn[i][0][3]*pQn[i][1][3]) - ((1/(N[i][1]-1) * pQn[i][1][8] * pQn[i][0][5] *pQn[i][0][3] ))
			//- ((1/(N[i][0]-1) * pQn[i][0][8]*pQn[i][1][5] * pQn[i][1][3])) + (1/((N[i][0]-1)*(N[i][1]-1))*pQn[i][0][8]*pQn[i][1][8] );
		TComplex nV4V4V3V3 = FourGap22(pQq,i,4,3,4,3)/ref_4p;//(pQn[i][0][4]*pQn[i][1][4]*pQn[i][0][3]*pQn[i][1][3]) - ((1/(N[i][1]-1) * pQn[i][1][7] * pQn[i][0][4] *pQn[i][0][3] ))
			//- ((1/(N[i][0]-1) * pQn[i][0][7]*pQn[i][1][4] * pQn[i][1][3])) + (1/((N[i][0]-1)*(N[i][1]-1))*pQn[i][0][7]*pQn[i][1][7] );

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

	Double_t event_weight_four = 1.0;
	Double_t event_weight_two = 1.0;
	Double_t event_weight_two_eta10 = 1.0;
	if(flags & FLUC_EBE_WEIGHTING){
		event_weight_four = Four(0,0,0,0).Re();
		event_weight_two = Two(0,0).Re();
		event_weight_two_eta10 = (QvectorQCeta10[kSubA][0][1]*QvectorQCeta10[kSubB][0][1]).Re();
	}

	for(int ih=2; ih < kNH; ih++){
		//for(int ihh=2; ihh<ih; ihh++){ //all SC
		for(int ihh=2, mm = (ih < kcNH?ih:kcNH); ihh<mm; ihh++){ //limited
			TComplex scfour = Four( ih, ihh, -ih, -ihh ) / Four(0,0,0,0).Re();
			
			fh_SC_with_QC_4corr[ih][ihh][fCBin]->Fill( scfour.Re(), event_weight_four );
			//QC_4p_value[ih][ihh] = scfour.Re();
		}

		// Finally we want 2p corr as like
		// 1/( M*(M-1) ) * [ QQ* - M ]
		// two(2,2) = Q2 Q2* - Q0 = Q2Q2* - M
		// two(0,0) = Q0 Q0* - Q0 = M^2 - M
		//two[ih] = Two(ih, -ih) / Two(0,0).Re();
		TComplex sctwo = Two(ih, -ih) / Two(0,0).Re();
		fh_SC_with_QC_2corr[ih][fCBin]->Fill( sctwo.Re(), event_weight_two );
		//QC_2p_value[ih] = sctwo.Re();
		// fill single vn  with QC without EtaGap as method 2
		fSingleVn[ih][2] = TMath::Sqrt(sctwo.Re());
		
		TComplex sctwo10 = (QvectorQCeta10[kSubA][ih][1]*TComplex::Conjugate(QvectorQCeta10[kSubB][ih][1])) / (QvectorQCeta10[kSubA][0][1]*QvectorQCeta10[kSubB][0][1]).Re();
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

	if(flags & FLUC_SCPT){
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

		// calculate Qn for each pt bin
		for(int ih=2; ih<SCNH; ih++){
			for(int ipt=0; ipt<N_ptbins; ipt++){
				Double_t pt_bin_min = ptbin_borders[ipt];
				Double_t pt_bin_max = ptbin_borders[ipt+1];
				QnA_pt[ih][ipt] = Get_Qn_pt(Eta_config[kSubA][0], Eta_config[kSubA][1], ih, ipt, pt_bin_min, pt_bin_max);
				QnB_pt[ih][ipt] = Get_Qn_pt(Eta_config[kSubB][0], Eta_config[kSubB][1], ih, ipt, pt_bin_min, pt_bin_max);

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
	}
}

//________________________________________________________________________
void AliJFFlucAnalysis::Terminate(Option_t *)
{
	//	fList = dynamic_cast<TList*> (GetOutputData(1));
	//	if(!fList) { Printf("ERROR: fList not availabe"); return;};
	//
	//
	cout<<"Sucessfully Finished"<<endl;
}
//________________________________________________________________________
//________________________________________________________________________
void AliJFFlucAnalysis::Fill_QA_plot( Double_t eta1, Double_t eta2 )
{
	Long64_t ntracks = fInputList->GetEntriesFast();
	for( Long64_t it=0; it < ntracks; it++){
		AliJBaseTrack *itrack = (AliJBaseTrack*)fInputList->At(it); // load track
		Double_t eta = itrack->Eta();
		Double_t phi = itrack->Phi();

		fh_phieta[fCBin]->Fill(phi,eta);
		fh_phietaz[fCBin]->Fill(phi,eta,fVertex[2]);

		if(TMath::Abs(eta) < eta1 || TMath::Abs(eta) > eta2)
			continue;

		Double_t phi_module_corr = itrack->GetWeight();// doing it in AliJFlucTask while filling track information.

		Double_t pt = itrack->Pt();
		Double_t effCorr = itrack->GetTrackEff();//fEfficiency->GetCorrection( pt, fEffFilterBit, fCent);
		Double_t effInv = 1.0/effCorr;
		fh_eta[fCBin]->Fill(eta,effInv);
		fh_pt[fCBin]->Fill(pt,effInv);
		fh_phi[fCBin][(int)(eta > 0.0)]->Fill( phi, effInv/phi_module_corr);
	}
	for(int iaxis=0; iaxis<3; iaxis++)
		fh_vertex[iaxis]->Fill( fVertex[iaxis] );

	fh_TrkQA_TPCvsCent->Fill(fCent,fTPCtrks);
	fh_TrkQA_TPCvsGlob->Fill(fGlbtrks,fTPCtrks);
	fh_TrkQA_FB32_vs_FB32TOF->Fill(fFB32trks,fFB32TOFtrks);
}

//________________________________________________________________________
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
TComplex AliJFFlucAnalysis::Get_Qn_pt(Double_t eta1, Double_t eta2, int harmonics, int ptbin, Double_t pt_min, Double_t pt_max)
{
	int nh = harmonics;
	TComplex Qn = TComplex(0,0);
	Double_t Sub_Ntrk =0;

	Long64_t ntracks = fInputList->GetEntriesFast();
	for( Long64_t it = 0; it < ntracks; it++){
		AliJBaseTrack *itrack = (AliJBaseTrack*)fInputList->At(it); // load track
		Double_t eta = itrack->Eta();
		if(eta < eta1 || eta > eta2)
			continue;
		Double_t pt = itrack->Pt();
		if(pt < pt_min || pt > pt_max)
			continue;
		Double_t phi = itrack->Phi();
		Double_t effCorr = itrack->GetTrackEff();//fEfficiency->GetCorrection( pt, fEffFilterBit, fCent);
		Double_t phi_module_corr = itrack->GetWeight();// doing it in AliJFlucTask while filling track information.
		Double_t tf = 1.0/(phi_module_corr*effCorr);
		Qn += TComplex(tf*TMath::Cos(nh*phi),tf*TMath::Sin(nh*phi));
		Sub_Ntrk += tf;
	}

	Qn /= Sub_Ntrk;
	int iside = (int)(eta1 > 0.0);
	NSubTracks_pt[iside][ptbin] = Sub_Ntrk;

	return Qn;
}
///________________________________________________________________________
/* new Function for QC method
   Please see Generic Framwork from Ante
   use Standalone method  */
//________________________________________________________________________
void AliJFFlucAnalysis::CalculateQvectorsQC(double etamin, double etamax){
	// calcualte Q-vector for QC method ( no subgroup )
	//init
	for(int ih=0; ih<kNH; ih++){
		for(int ik=0; ik<nKL; ++ik){
			QvectorQC[ih][ik] = TComplex(0,0);
			for(int isub=0; isub<2; isub++){
				QvectorQCeta10[isub][ih][ik] = TComplex(0,0);
			}
		}
	} // for max harmonics
	//Calculate Q-vector with particle loop
	Long64_t ntracks = fInputList->GetEntriesFast(); // all tracks from Task input
	for( Long64_t it=0; it<ntracks; it++){
		AliJBaseTrack *itrack = (AliJBaseTrack*)fInputList->At(it); // load track
		Double_t eta = itrack->Eta();
		// track Eta cut Note! pt cuts already applied in AliJFFlucTask.cxx
		// Do we need arbitary Eta cut for QC method?
		// fixed eta ranged -0.8 < eta < 0.8 for QC
//				if( TMath::Abs(eta) > fEta_max || TMath::Abs(eta) < fEta_min ) continue; << this is SP cut
//				if( TMath::Abs(eta) > 0.8 ) continue;  //   << this is old QC cut
		// we need more configuration for to study eta dep of SC(m,n) with QC method.
		// eta cut among all tracks (this is not same with SC(m,n) SP method (SP method needs symmetric eta range)//
		//if( eta < fQC_eta_cut_min || eta > fQC_eta_cut_max)
		if( eta < -etamax || eta > etamax)
			continue;
		/////////////////////////////////////////////////

		int isub = (int)(eta > 0.0);
		Double_t phi = itrack->Phi();
		Double_t pt = itrack->Pt();

		Double_t effCorr = itrack->GetTrackEff();//fEfficiency->GetCorrection( pt, fEffFilterBit, fCent);	
		Double_t phi_module_corr = itrack->GetWeight();

		for(int ih=0; ih<kNH; ih++){
			Double_t tf = 1.0;
			TComplex q[nKL];
			for(int ik=0; ik<nKL; ik++){
				q[ik] = TComplex(tf*TMath::Cos(ih*phi),tf*TMath::Sin(ih*phi));
				QvectorQC[ih][ik] += q[ik];

				//this is for normalized SC ( denominator needs an eta gap )
				if(TMath::Abs(eta) > etamin)//fQC_eta_gap_half)
					QvectorQCeta10[isub][ih][ik] += q[ik];

				tf *= 1.0/(phi_module_corr*effCorr);
			}
		}
	} // track loop done.
}
//________________________________________________________________________
TComplex AliJFFlucAnalysis::Q(int n, int p){
	// Return QvectorQC
	// Q{-n, p} = Q{n, p}*
	if(n >= 0)
		return QvectorQC[n][p];
	return TComplex::Conjugate(QvectorQC[-n][p]);
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
/*void AliJFFlucAnalysis::SetPhiModuleHistos( int cent, int sub, TH1D *hModuledPhi){
	// hPhi histo setter
	h_phi_module[cent][sub] = hModuledPhi;
}*/

/*int AliJFFlucAnalysis::GetCentralityClass(Double_t fCent){
	for(UInt_t iCbin = 0; iCbin < NCentBin; iCbin++){
		if(fCent > CentBin[iCbin] && fCent < CentBin[iCbin+1])
			return iCbin;
	}
	return -1;
}

int AliJFFlucAnalysis::GetMultiplicityBin(Double_t fMult, BINNING _binning){
	for(UInt_t iMbin = 0; iMbin < NMultBin[_binning]; iMbin++){
		//if(fMult > MultBin[iMbin] && fMult < MultBin[iMbin+1])
		if(fMult > pMultBin[_binning][iMbin][0] && fMult < pMultBin[_binning][iMbin][1])
			return iMbin;
	}
	return -1;
}*/

int AliJFFlucAnalysis::GetBin(Double_t fq, BINNING _binning){
	for(UInt_t iMbin = 0; iMbin < NBin[_binning]; iMbin++){
		if(fq > pBin[_binning][iMbin][0] && fq < pBin[_binning][iMbin][1])
			return iMbin;
	}
	return -1;
}

