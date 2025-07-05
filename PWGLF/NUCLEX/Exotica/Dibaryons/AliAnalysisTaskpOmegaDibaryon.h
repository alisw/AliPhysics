#ifndef ALIANALYSISTASKPOMEGADIBARYON_H
#define ALIANALYSISTASKPOMEGADIBARYON_H

// ROOT includes
#include <THashList.h>
#include <TH1.h>

// AliRoot includes
//Janik include 
//(https://github.com/alisw/AliPhysics/blob/master/PWGLF/NUCLEX/Hypernuclei/DoubleHypNuc/AliAnalysisTaskDoubleHypNucTree.h)
#define HomogeneousField ///homogenous field in z direction
#include <algorithm>
#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"
#include "AliInputEventHandler.h"
#include "AliAODTrack.h"
#include "AliPIDResponse.h"
#include "AliAODVertex.h"
#include "AliAODEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDpid.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliEventCuts.h"
#include "AliKalmanTrack.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliMCEvent.h"
#include "AliMCVertex.h"
#include "AliPID.h"
#include "AliStack.h"
#include "AliVertexerTracks.h"
#include "AliVVertex.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include <fstream>
#include "KFParticle.h"
#include "KFVertex.h"
#include "KFPVertex.h"
#include "KFPTrack.h"
#include <vector>
#include "TCanvas.h"
#include "TChain.h"
#include <TClonesArray.h>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TObject.h"
#include "TPDGCode.h"
#include "TLorentzVector.h"
#include "TRandom2.h"
#include "TTree.h"
#include "TVector3.h"

//my include
#include <AliAnalysisTaskSE.h>
#include <AliAODEvent.h>
#include <AliAnalysisUtils.h>
#include <AliPIDResponse.h>
#include <AliAODv0.h>
#include <AliAODcascade.h>
#include <AliEventPoolManager.h>

#include <AliMultSelection.h>
#include <AliEventplane.h>
#include <AliQnCorrectionsManager.h>
#include <AliAnalysisTaskFlowVectorCorrections.h>
#include <AliQnCorrectionsQnVector.h>

class TH1F;
class TH2F;
class AliESDEvent;
class AliESDpid;
class AliESDVertex;
class AliESDTrack;
class AliAODVertex;
class AliAODTrack;
class AliAODEvent;
class AliVEvent;
class KFParticle;
class KFVertex;
class AliKalmanTrack;
class AliITStrackerMI;

class KFParticleDibaryon : public KFParticle {
 public:
  Bool_t CheckDaughter(KFParticle daughter) {
    Float_t m[8], mV[36], D[3][3];
    Bool_t CheckMeasurement = KFParticleBase::GetMeasurement(daughter, m, mV, D);
    if(!CheckMeasurement) {
      if(fWarnings) AliWarning("Checking KF Daughter failed: GetMeasurement() returned false"); return kFALSE;
    }
    for(Int_t i=0;i<3;i++){
      for(Int_t j=0;j<3;j++) {
	if(TMath::Abs(D[i][j])>100000) {
	  if(fWarnings) AliWarning("Checking KF Daughter failed: Covariance limit exceeded"); return kFALSE;
	}
      }
    }
    return kTRUE;  
  }

  void ActivateWarnings() {fWarnings = kTRUE;}
 protected:
  Bool_t fWarnings = kFALSE;
  /*
    Float_t m[8], mV[36], D[3][3];
    if (KFParticleBase::GetMeasurement(daughter, m, mV, D)) return kTRUE;
    return kFALSE;
  */
};

class KFParticleBase2 : public KFParticle {
 public:
  float GetDStoPointBz( float B, const float xyz[3], float dsdr[6], const float* param) const
  { 

    cout<<"GetDStoPointBz[0]"<<endl;

    const float& x  = param[0];
    const float& y  = param[1];
    const float& z  = param[2];
    const float& px = param[3];
    const float& py = param[4];
    const float& pz = param[5];
    
    //cout<<"GetDStoPointBz[1]"<<endl;

    const float kCLight = 0.000299792458f;
    float bq = B*fQ*kCLight;
    float pt2 = px*px + py*py;
    float p2 = pt2 + pz*pz;     
    float dx = xyz[0] - x;
    float dy = xyz[1] - y; 
    float dz = xyz[2] - z; 
    //cout<<"GetDStoPointBz[1-1]==dx:px:dy:py = "<<dx<<":"<<px<<":"<<dy<<";"<<py<<endl;
    float a = dx*px+dy*py;
    float dS(0.f);

    //cout<<"GetDStoPointBz[2]"<<endl;

    //cout<<"a = "<<a<<endl;
    //cout<<"bq = "<<bq<<endl;
    float abq = bq*a;
    //cout<<"abq = "<<abq<<endl;
    //cout<<"GetDStoPointBz[3]"<<endl;
    
    const float LocalSmall = 1.e-8f;
    bool mask = ( fabs(bq)<LocalSmall );
    if(mask && p2>1.e-4f)
      {
	cout<<"GetDStoPointBz[3-1]"<<endl;
	dS = (a + dz*pz)/p2;  
	cout<<"GetDStoPointBz[3-2]"<<endl;
	dsdr[0] = -px/p2;
	cout<<"[1]"<<endl;
	dsdr[1] = -py/p2;
	cout<<"[2]"<<endl;
	dsdr[2] = -pz/p2;
	cout<<"[3]"<<endl;
	cout<<"dz = "<<dz;
	dsdr[3] = (dx*p2 - 2.f* px *(a + dz *pz))/(p2*p2);
	cout<<"[4]"<<endl;
	dsdr[4] = (dy*p2 - 2.f* py *(a + dz *pz))/(p2*p2);
	cout<<"[5]"<<endl;
	dsdr[5] = (dz*p2 - 2.f* pz *(a + dz *pz))/(p2*p2);
	cout<<"GetDStoPointBz[3-3]"<<endl;
      }
    if(mask)
      { 
	cout<<"GetDStoPointBz[3-3]"<<endl;
	return dS;
      }

    cout<<"GetDStoPointBz[4]"<<endl;    
    dS = atan2( abq, pt2 + bq*(dy*px -dx*py) )/bq; 
    float bs= bq*dS; 
    float s = sin(bs), c = cos(bs);
    
    cout<<"GetDStoPointBz[5]"<<endl;

    if(fabs(bq) < LocalSmall)
      bq = LocalSmall;
    float bbq = bq*(dx*py - dy*px) - pt2;   
    dsdr[0] = (px*bbq - py*abq)/(abq*abq + bbq*bbq);
    dsdr[1] = (px*abq + py*bbq)/(abq*abq + bbq*bbq);
    dsdr[2] = 0;
    dsdr[3] = -(dx*bbq + dy*abq + 2.f*px*a)/(abq*abq + bbq*bbq);
    dsdr[4] = (dx*abq - dy*bbq - 2.f*py*a)/(abq*abq + bbq*bbq);
    dsdr[5] = 0;
   
    float sz(0.f);
    float cCoeff =  (bbq*c - abq*s) - pz*pz ;
    if(fabs(cCoeff) > 1.e-8f)
      sz = (dS*pz - dz)*pz / cCoeff;
    
    cout<<"GetDStoPointBz[6]"<<endl;

    float dcdr[6] = {0.f};
    dcdr[0] = -bq*py*c - bbq*s*bq*dsdr[0] + px*bq*s - abq*c*bq*dsdr[0];
    dcdr[1] =  bq*px*c - bbq*s*bq*dsdr[1] + py*bq*s - abq*c*bq*dsdr[1];
    dcdr[3] = (-bq*dy-2*px)*c - bbq*s*bq*dsdr[3] - dx*bq*s - abq*c*bq*dsdr[3];
    dcdr[4] = ( bq*dx-2*py)*c - bbq*s*bq*dsdr[4] - dy*bq*s - abq*c*bq*dsdr[4];
    dcdr[5] = -2*pz;

    cout<<"GetDStoPointBz[7]"<<endl;

    for(int iP=0; iP<6; iP++)
      dsdr[iP] += pz*pz/cCoeff*dsdr[iP] - sz/cCoeff*dcdr[iP];
    dsdr[2] += pz/cCoeff;
    dsdr[5] += (2.f*pz*dS - dz)/cCoeff;
    
    dS += sz;
   
    bs= bq*dS;
    s = sin(bs), c = cos(bs);
    
    cout<<"GetDStoPointBz[8]"<<endl;

    float sB, cB;
    const float kOvSqr6 = 1.f/sqrt(float(6.f));
    
    if(LocalSmall < fabs(bs))
      {
	sB = s/bq;
	cB = (1.f-c)/bq;
      }
    else
      {
	sB = (1.f-bs*kOvSqr6)*(1.f+bs*kOvSqr6)*dS;
	cB = .5f*sB*bs;
      }
    
    cout<<"GetDStoPointBz[9]"<<endl;

    float p[5];
    p[0] = x + sB*px + cB*py;
    p[1] = y - cB*px + sB*py;
    p[2] = z +  dS*pz;
    p[3] =          c*px + s*py;
    p[4] =         -s*px + c*py;
  
    dx = xyz[0] - p[0];
    dy = xyz[1] - p[1];
    dz = xyz[2] - p[2];
    a = dx*p[3]+dy*p[4] + dz*pz;
    abq = bq*a;
    
    dS += atan2( abq, p2 + bq*(dy*p[3] -dx*p[4]) )/bq;
    cout<<"GetDStoPointBz[10]"<<endl;
 
    return dS;
  }

};

class AliAnalysisTaskpOmegaDibaryon : public AliAnalysisTaskSE {
 public:
  enum FlowMethod{
    kOFF = -1,
    kEP  = 0,
      kSP  = 1
  };

  enum QnDetector{
    kNone      = -1,
    kFullTPC   = 0,
    kTPCNegEta = 1,
    kTPCPosEta = 2,
    kFullV0    = 3,
    kV0A       = 4,
    kV0C       = 5
  };

  enum PdgCodeType_t {
    kPDGPion,
    kPDGProton,
    kPDGLambda,
    kPDGK0s,
    kPDGXi,
    kPDGOmega
  };

  static const Int_t fgkPdgCode[];

  AliAnalysisTaskpOmegaDibaryon();
  AliAnalysisTaskpOmegaDibaryon(const char *name);
  virtual ~AliAnalysisTaskpOmegaDibaryon();
  
  virtual void  SetTrigger(UInt_t ktriggerInt=AliVEvent::kINT7){ ftrigBit=ktriggerInt; }
  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t *option);

  //setter for event selection
  void SetCentralityEstimator(TString estimator) {fEstimator = estimator;}
  void SetCentralityMin(Float_t min) {fCentralityMin = min;}
  void SetCentralityMax(Float_t max) {fCentralityMax = max;}

 protected:
  Bool_t   EventSelection(AliAODEvent *data);
  Bool_t   ProtonSelection(AliAODTrack *trk);
  Double_t InvMassLambda(AliAODcascade *casc);
  Double_t InvMassAntiLambda(AliAODcascade *casc);
  Double_t InvMassXi(AliAODcascade *casc);
  Double_t InvMassOmega(AliAODcascade *casc);
  Double_t LambdaCosPointingAngle(AliAODcascade *casc,const Double_t *DecayVtx,const Float_t *point) const;
  Double_t DecayLengthXY(const Double_t *DecayVtx,const Float_t *point) const;
  Double_t xiDecayLengthXY(const Double_t *xiDecayVtx,const Float_t *point) const;
  Double_t InvMasslambda(AliAODv0 *v0);
  Double_t InvMassK0(AliAODv0 *v0);
  Double_t InvMassSelf(AliExternalTrackParam *track1, AliExternalTrackParam *track2, Double_t daughtermass1, Double_t daughtermass2);
  Double_t InvMassAntilambda(AliAODv0 *v0);
  Double_t CalculateInvMassAntilambda(AliAODv0 *v0,AliAODTrack *antiprotontrk,AliAODTrack *piontrk);
  Double_t CosPointingAngle(AliAODv0 *v0,const Double_t *DecayVtx,const Float_t *point) const;
  Double_t OpenAngle(Double_t px1,Double_t py1,Double_t pz1,
		     Double_t px2,Double_t py2,Double_t pz2);
  Double_t InvariantMass(Double_t px1,Double_t py1,Double_t pz1,
			 Double_t px2,Double_t py2,Double_t pz2,Double_t energysum);
  //Bool_t ExtractQnVector();
  Double_t ExtractQnVector();
  const AliQnCorrectionsQnVector *GetQnVectorFromList(const TList *qnlist, const char* subdetector, const char *expcorr, const char *altcorr);

  //Convert AOD to ESD Vertex                                                                                                                    
  AliESDVertex *AODToESDVertex(const AliAODVertex &AODVtx);
  KFVertex CreateKFVertex(const AliVVertex& vertex);
  KFParticle CreateKFParticle(AliExternalTrackParam& track, Double_t Mass, Int_t Charge); 

  Bool_t IsQualityTrack(AliAODTrack *track);
  Double_t GetTOFSignal(AliAODTrack& trackHe, Double_t tP);
  Double_t CalculateBetaTOF(AliAODTrack &track);
  Double_t CalculateMassSquareTOF(AliAODTrack &track);
  Double_t CalculateProtonSigma(Double_t pT, Double_t massSq, Double_t charge);
  Bool_t IsQualityProton(AliAODTrack *track);
  Double_t CalculateSigmadEdxITS(AliAODTrack &Track);

 private:
  AliVEvent       *fEvent;    //! ESD object
  AliESDEvent     *fESD;    //! ESD object
  AliAODEvent     *fAOD;    //! AOD object
  AliAODHeader    *fHeader; //! AOD header
  AliPIDResponse  *fPIDResponse;
  UInt_t          ftrigBit;
  TString fEstimator;//V0[M|A|C], ZN[A|C], CL[0|1]
  AliMultSelection *fMultSelection;
  Float_t fCentralityMain;
  Float_t fCentralityMin;
  Float_t fCentralityMax;

  Bool_t fIsFlowTask;
  AliQnCorrectionsManager *fFlowQnVectorMgr;
  THashList             *fOutputList; //! Output list

  //==================== variable
  Int_t nevt;
  //========== Event selection
  Bool_t isEvtMB;
  Bool_t isPileupEvt;
  Bool_t isGoodEvt;
  Float_t centralityV0M;
  Float_t centralityCL0;
  Float_t centralityCL1;
  Float_t centralityV0A;
  Float_t centralityV0C;
  Float_t centralityZNA;
  Float_t centralityZNC;
  Float_t centralityMain;
    
  //========== Event
  Int_t    nTracks;
  Int_t    nV0s;
  Int_t    nCascades;
  Float_t  vecTarget[3];
  Double_t kMagF;

  Double_t PrimVertexKF[3];
  const AliAODVertex *vertex;
  AliESDVertex  *primVertex;
  KFVertex       primKFVertex;

  //========== Track
  TLorentzVector *lorentzsum;
  TLorentzVector *particleProton;
  TLorentzVector *particlePion;
  TVector3 *h;

  int np;
  int nap; 
  int npi;
  int nppi;

  AliAODTrack *track;  
  Double_t dcaxy,dcaz,dca;
  Double_t dcaxyn,dcazn,dcan;
  Float_t  dDCA[2]; 
  Float_t  cDCA[3]; 
  Double_t v0Vtx[3];
  Double_t track_pT;

  Bool_t Daugproton;
  Bool_t Daugpion;
  Bool_t isProton;
  Bool_t isAntiProton;
  Bool_t isNPion;
  Bool_t isPPion;

  Bool_t isBachpion;
  Bool_t isBachkaon;

  Bool_t isLambdadecayProton;
  Bool_t isLambdadecayAntiProton;
  Bool_t isLambdadecayPPion;
  Bool_t isLambdadecayNPion;
  Bool_t isBachPPion;
  Bool_t isBachNPion;
  Bool_t isBachPKaon;
  Bool_t isBachNKaon;
  Bool_t pTCut_Proton;
  Bool_t pTCut_Pion;
  Bool_t pTCut_Kaon;

  //========== Array
  std::vector < int >   ProtonArray;
  std::vector < int >   AntiProtonArray;
  std::vector < int >   NPionArray;
  std::vector < int >   PPionArray;
  std::vector < int >   LambdadecayProtonArray;
  std::vector < int >   LambdadecayAntiProtonArray;
  std::vector < int >   LambdadecayPPionArray;
  std::vector < int >   LambdadecayNPionArray;
  std::vector < int >   XidecayBachPPionArray;
  std::vector < int >   XidecayBachNPionArray;
  std::vector < int >   XidecayBachPKaonArray;
  std::vector < int >   XidecayBachNKaonArray;

  //========== TrackParam
  AliExternalTrackParam *exProtonTrack;
  AliExternalTrackParam *exPionTrack;
  AliExternalTrackParam *exBachPionTrack;
  AliExternalTrackParam *exBachKaonTrack;

  //========== KF particle
  AliAODTrack *ProtonTrack;
  AliAODTrack *PionTrack;

  KFParticle KFProtonTrack;
  KFParticle KFPPionTrack;
  KFParticle KFPionTrack;
  KFParticle KFMother;
  KFParticle KFMother_K0s;

  Double_t SecVertexKF[3];
  Double_t SecVertexErrKF[3];
  AliESDVertex *KFsecVertex;
  
  Double_t fDCADaugter;
  Double_t fDcaSecDaughterProton;
  Double_t fDcaSecDaughterPion;
  Double_t impar[2];
  Double_t dd[3];
  Double_t fPA;
  Double_t DecLength;

  //========== V0
  AliAODv0 *v0;
  AliAODTrack* ptrack;
  AliAODTrack* ntrack;
  Double_t dcaxyp_v0;
  Double_t dcazp_v0;
  Double_t dcap_v0;
  Double_t dcaxyn_v0;
  Double_t dcazn_v0;
  Double_t dcan_v0;
  Float_t  dDCAp_v0[2];
  Float_t  dDCAn_v0[2];
  Float_t  cDCAp_v0[3];
  Float_t  cDCAn_v0[3];
  Double_t v0Vtx_v0[3];
  Double_t energy_v0;
  Double_t fcpa_v0;
  Double_t ftransradius_v0;
  Bool_t Daugproton_v0;
  Bool_t Daugantiproton_v0;
  Bool_t Daugpion_v0;
  Bool_t Daugnpion_v0;
  Bool_t lambda_v0;
  Bool_t antilambda_v0;

  //========== KF particle cascade 
  Double_t ndp;
  Double_t ndap;
  Double_t ndpp;
  Double_t ndnp;
  Double_t nbpp;
  Double_t nbnp;
  Double_t nbpk;
  Double_t nbnk;

  //Float_t DCAtoPrimVtxXY;
  //Float_t DCAtoPrimVtxZ;
  Float_t DCAtoPrimVtx;

  AliAODTrack *BachNPionTrack;
  KFParticle KFSubMother;
  KFParticle KFBachPionTrack;
  KFParticle KFBachKaonTrack;
  KFParticle KFMother_Omega;

  Double_t SecVertexKF_Xi[3];
  Double_t SecVertexErrKF_Xi[3];
  AliESDVertex *KFsecVertex_Xi;

  Double_t SecVertexKF_Omega[3];
  Double_t SecVertexErrKF_Omega[3];
  AliESDVertex *KFsecVertex_Omega;

  TLorentzVector *lorentzsum_Xi;
  TLorentzVector *particleBachPion;
  TVector3 *h_Xi;

  TLorentzVector *lorentzsum_Omega;
  TLorentzVector *particleBachKaon;

  Double_t dd_Xi[3];
  Double_t fPA_Xi;
  Double_t DecLength_Xi;

  Float_t DCAdauXY_Xi;
  Float_t DCAdauZ_Xi;
  Float_t DCAdau_Xi;

  Int_t ftrkID_daughter1;
  Int_t ftrkID_daughter2;
  Int_t ftrkID_daughter3;

  //========== cascade
  AliAODcascade *casc;
  AliAODTrack* pTrackXi;
  AliAODTrack* nTrackXi;
  AliAODTrack* bachTrackXi;
  double_t dcaxyp;
  double_t dcazp;
  double_t dcap;
  double_t dcaxyb;
  double_t dcazb;
  double_t dcab;
  Float_t  dDCAp[2];
  Float_t  cDCAp[3];
  Float_t  dDCAn[2];
  Float_t  cDCAn[3];
  Float_t  dDCAb[2];
  Float_t  cDCAb[3];
  double_t v0Vtx_casc[3];
  double_t xivertex[3];
  Double_t lambdacpa;
  Double_t lambdatransradius;
  Double_t xicpa;
  Double_t xitransradius;
  Double_t v0vtxx;
  Double_t v0vtxy;
  Double_t v0vtxz;
  Double_t xivtxx;
  Double_t xivtxy;
  Double_t xivtxz;

  //========= TTree_omegam
  TTree *fTree_omegam; //< 
  Long_t fevtID_omegam;
  Double_t fzvertex_omegam;
  Double_t fcentrality_omegam;

  //Daughter track                                                                                          
  Float_t fdca_daughter1_omegam;
  Float_t fdca_daughter2_omegam;
  Float_t fdca_daughter3_omegam;
  Int_t ftrkID_daughter1_omegam;
  Int_t ftrkID_daughter2_omegam;
  Int_t ftrkID_daughter3_omegam;

  //Lambda                                                                             
  Double_t fcpa_lambda_omegam;
  Double_t fpcaxy_lambda_omegam;
  Double_t fpca_lambda_omegam;
  Double_t fdeclength_lambda_omegam;
  Float_t  fchi2_lambda_omegam;
  Double_t fmass_lamdba_omegam;
  Double_t fdcatoPVxy_lambda_omegam;
  Double_t fdcatoPV_lambda_omegam;
  Double_t fpt_lambda_omegam;

  //Omega                                                                                  
  Double_t fcpa_omega_omegam;
  Double_t fpcaxy_omega_omegam;
  Double_t fpca_omega_omegam;
  Double_t fdeclength_omega_omegam;
  Float_t  fchi2_omega_omegam;
  Double_t fmass_omega_omegam;
  Double_t fpx_omega_omegam;
  Double_t fpy_omega_omegam;
  Double_t fpz_omega_omegam;
  Double_t fpt_omega_omegam;
  Double_t feta_omega_omegam;
  Double_t fphi_omega_omegam;
  Double_t fdcatoPVxy_omega_omegam;
  Double_t fdcatoPV_omega_omegam;
  Double_t fct_omega_omegam;
  Double_t fpx_omega_omegam_mc;
  Double_t fpy_omega_omegam_mc;
  Double_t fpz_omega_omegam_mc;
  Double_t fpt_omega_omegam_mc;

  //========= TTree_omegap
  TTree *fTree_omegap; //< 
  Long_t fevtID_omegap;
  Double_t fzvertex_omegap;
  Double_t fcentrality_omegap;

  //Daughter track                                                                                          
  Float_t fdca_daughter1_omegap;
  Float_t fdca_daughter2_omegap;
  Float_t fdca_daughter3_omegap;
  Int_t ftrkID_daughter1_omegap;
  Int_t ftrkID_daughter2_omegap;
  Int_t ftrkID_daughter3_omegap;

  //Lambda                                                                             
  Double_t fcpa_lambda_omegap;
  Double_t fpcaxy_lambda_omegap;
  Double_t fpca_lambda_omegap;
  Double_t fdeclength_lambda_omegap;
  Float_t  fchi2_lambda_omegap;
  Double_t fmass_lamdba_omegap;
  Double_t fdcatoPVxy_lambda_omegap;
  Double_t fdcatoPV_lambda_omegap;
  Double_t fpt_lambda_omegap;

  //Omega                                                                                  
  Double_t fcpa_omega_omegap;
  Double_t fpcaxy_omega_omegap;
  Double_t fpca_omega_omegap;
  Double_t fdeclength_omega_omegap;
  Float_t  fchi2_omega_omegap;
  Double_t fmass_omega_omegap;
  Double_t fpx_omega_omegap;
  Double_t fpy_omega_omegap;
  Double_t fpz_omega_omegap;
  Double_t fpt_omega_omegap;
  Double_t feta_omega_omegap;
  Double_t fphi_omega_omegap;
  Double_t fdcatoPVxy_omega_omegap;
  Double_t fdcatoPV_omega_omegap;
  Double_t fct_omega_omegap;
  Double_t fpx_omega_omegap_mc;
  Double_t fpy_omega_omegap_mc;
  Double_t fpz_omega_omegap_mc;
  Double_t fpt_omega_omegap_mc;

  //========= TTree sidebandm
  TTree *fTree_sidebandm; //< 
  Long_t fevtID_sidebandm;
  Double_t fzvertex_sidebandm;
  Double_t fcentrality_sidebandm;

  //Daughter track                                                                                                   
  Float_t fdca_daughter1_sidebandm;
  Float_t fdca_daughter2_sidebandm;
  Float_t fdca_daughter3_sidebandm;
  Int_t ftrkID_daughter1_sidebandm;
  Int_t ftrkID_daughter2_sidebandm;
  Int_t ftrkID_daughter3_sidebandm;

  //Lambda                                                                                                           
  Double_t fcpa_lambda_sidebandm;
  Double_t fpcaxy_lambda_sidebandm;
  Double_t fpca_lambda_sidebandm;
  Double_t fdeclength_lambda_sidebandm;
  Float_t  fchi2_lambda_sidebandm;
  Double_t fmass_lamdba_sidebandm;
  Double_t fdcatoPVxy_lambda_sidebandm;
  Double_t fdcatoPV_lambda_sidebandm;
  Double_t fpt_lambda_sidebandm;

  //Omega                                                                                                               
  Double_t fcpa_omega_sidebandm;
  Double_t fpcaxy_omega_sidebandm;
  Double_t fpca_omega_sidebandm;
  Double_t fdeclength_omega_sidebandm;
  Float_t  fchi2_omega_sidebandm;
  Double_t fmass_omega_sidebandm;
  Double_t fpx_omega_sidebandm;
  Double_t fpy_omega_sidebandm;
  Double_t fpz_omega_sidebandm;
  Double_t fpt_omega_sidebandm;
  Double_t feta_omega_sidebandm;
  Double_t fphi_omega_sidebandm;
  Double_t fdcatoPVxy_omega_sidebandm;
  Double_t fdcatoPV_omega_sidebandm;
  Double_t fct_omega_sidebandm;
  Double_t fpx_omega_sidebandm_mc;
  Double_t fpy_omega_sidebandm_mc;
  Double_t fpz_omega_sidebandm_mc;
  Double_t fpt_omega_sidebandm_mc;

  //========= TTree sidebandp
  TTree *fTree_sidebandp; //< 
  Long_t fevtID_sidebandp;
  Double_t fzvertex_sidebandp;
  Double_t fcentrality_sidebandp;

  //Daughter track                                                                                                   
  Float_t fdca_daughter1_sidebandp;
  Float_t fdca_daughter2_sidebandp;
  Float_t fdca_daughter3_sidebandp;
  Int_t ftrkID_daughter1_sidebandp;
  Int_t ftrkID_daughter2_sidebandp;
  Int_t ftrkID_daughter3_sidebandp;

  //Lambda                                                                                                           
  Double_t fcpa_lambda_sidebandp;
  Double_t fpcaxy_lambda_sidebandp;
  Double_t fpca_lambda_sidebandp;
  Double_t fdeclength_lambda_sidebandp;
  Float_t  fchi2_lambda_sidebandp;
  Double_t fmass_lamdba_sidebandp;
  Double_t fdcatoPVxy_lambda_sidebandp;
  Double_t fdcatoPV_lambda_sidebandp;
  Double_t fpt_lambda_sidebandp;

  //Omega                                                                                                               
  Double_t fcpa_omega_sidebandp;
  Double_t fpcaxy_omega_sidebandp;
  Double_t fpca_omega_sidebandp;
  Double_t fdeclength_omega_sidebandp;
  Float_t  fchi2_omega_sidebandp;
  Double_t fmass_omega_sidebandp;
  Double_t fpx_omega_sidebandp;
  Double_t fpy_omega_sidebandp;
  Double_t fpz_omega_sidebandp;
  Double_t fpt_omega_sidebandp;
  Double_t feta_omega_sidebandp;
  Double_t fphi_omega_sidebandp;
  Double_t fdcatoPVxy_omega_sidebandp;
  Double_t fdcatoPV_omega_sidebandp;
  Double_t fct_omega_sidebandp;
  Double_t fpx_omega_sidebandp_mc;
  Double_t fpy_omega_sidebandp_mc;
  Double_t fpz_omega_sidebandp_mc;
  Double_t fpt_omega_sidebandp_mc;

  //========= TTree Cascade loop                                                                          
  TTree *fTree_cascade; //<

  Double_t fpt_daughter1;
  Double_t fpt_daughter2;
  Double_t fpt_daughter3;

  Double_t fpt_daughter1_cascade;
  Double_t fpt_daughter2_cascade;
  Double_t fpt_daughter3_cascade;

  Double_t fcpa_lambda_cascade;
  Double_t fpca_lambda_cascade;
  Double_t fpt_lambda_cascade;

  Double_t fcpa_omega_cascade;
  Double_t fpca_omega_cascade;
  Double_t fmass_omega_cascade;
  Double_t fpt_omega_cascade;

  //========= TTree_proton                                                                                                            
  TTree *fTree_proton; //< 
  Long_t fevtID_proton;
  Double_t fzvertex_proton;
  Double_t fcentrality_proton;

  Int_t ftrkID_proton;
  Double_t fdcaxy_proton;
  Double_t fdcaz_proton;
  Double_t fpx_proton;
  Double_t fpy_proton;
  Double_t fpz_proton;
  Double_t fpt_proton;
  Double_t fpTPC_proton;
  Double_t fnSigmaTPC_proton;
  Double_t fnSigmaTOF_proton;
  Double_t fTrkMassTOF_proton;
  Double_t fnSigmaITS_proton;
  Double_t feta_proton;
  Double_t fphi_proton;

  //========= TTree_antiproton  
  TTree *fTree_antiproton; //< 
  Long_t fevtID_antiproton;
  Double_t fzvertex_antiproton;
  Double_t fcentrality_antiproton;

  Int_t ftrkID_antiproton;
  Double_t fdcaxy_antiproton;
  Double_t fdcaz_antiproton;
  Double_t fpx_antiproton;
  Double_t fpy_antiproton;
  Double_t fpz_antiproton;
  Double_t fpt_antiproton;
  Double_t fpTPC_antiproton;
  Double_t fnSigmaTPC_antiproton;
  Double_t fnSigmaTOF_antiproton;
  Double_t fTrkMassTOF_antiproton;
  Double_t fnSigmaITS_antiproton;
  Double_t feta_antiproton;
  Double_t fphi_antiproton;

  AliAnalysisTaskpOmegaDibaryon(const AliAnalysisTaskpOmegaDibaryon&); // not implemented
  AliAnalysisTaskpOmegaDibaryon& operator=(const AliAnalysisTaskpOmegaDibaryon&); // not implemented

  ClassDef(AliAnalysisTaskpOmegaDibaryon,2);

};

#endif
