#ifndef ALIFEMTOCORRFCTNPDTHE3_H
#define ALIFEMTOCORRFCTNPDTHE3_H

#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TString.h"

#include "AliFemtoCorrFctn.h"

#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"

/// \class AliFemtoCorrFctnpdtHe3
/// \brief KStar correlation function for p-d/t/He3 
///         Based on AliFemtoCorrFctnKStar.cxx!
/// \authors  Dong-Fang Wang, Fudan University, China, <dongfang.wang@cern.ch>
///
Double_t MassBandFunc(Double_t *x, Double_t *par);
class AliFemtoCorrFctnpdtHe3 : public AliFemtoCorrFctn {
    public:
        AliFemtoCorrFctnpdtHe3();
        AliFemtoCorrFctnpdtHe3(const char* title,
                        const int nbins,
                        const float KStarLo,
                        const float KStarHi);
        AliFemtoCorrFctnpdtHe3(const AliFemtoCorrFctnpdtHe3 & aCorrFctn);
        virtual ~AliFemtoCorrFctnpdtHe3();
        AliFemtoCorrFctnpdtHe3& operator=(const AliFemtoCorrFctnpdtHe3& aCorrFctn);
        
        virtual AliFemtoString Report();
        virtual TList* GetOutputList();
        virtual void Finish();
        void Write();

        virtual void AddRealPair(AliFemtoPair* aPair);
        virtual void AddMixedPair(AliFemtoPair* aPair);
        void SetIsHe3Pair(int isOrNot);
        AliFemtoPair * ChangeP2Mom(AliFemtoPair* aPair);
	void SetfUseVelGate(int aUseVelGate);
	void SetVelGateInit(bool aUseVelGate);
        int ReVelocityGate(AliFemtoPair* aPair);
        void SetP1AndP2Mass(float p1Mass,float p2Mass);

	//\ 2022.2.15
	void SetfHighCF(bool aHighCF);
	void SetHighCFInit(bool aHighCF);

	void SetfSideBand(bool aSideBand);
	void SetSideBandTF1Init(bool aSideBand);
	void SetSideBandHistInit(bool aSideBand);

	void SetTF1ParaUp3Sigma(float *para);
	void SetTF1ParaLow3Sigma(float *para);
	void SetTF1PareSideBandUp(float *para);
	void SetTF1PareSideBandLow(float *para);
	
	
	void FillSideBandNum(AliFemtoPair* aPair);
	void FillSideBandDum(AliFemtoPair* aPair);

	void SetfUsePt(int aUsePt);
	void SetfUseDPhiDEtaQA(int aUseDPhiDEtaQA);
	void SetDPhiDEtaQAInit(bool aDPhiDEtaQA);
        float ReAvgDphi(AliFemtoPair* aPair);
	//\ 2022.9.5
	void SetUseStavinskyMethod(int aUse);
	void SetStaSkyBkgInit(bool aInit);
	AliFemtoPair *InversePair(AliFemtoPair* aPair);
	void SetUse2DpTvsKStar(int aUse);
        void Set2DpTvsKStarInit(bool aInit,int nbinsks,float lowks,float upks,int nbinspT,float lowpT,float uppT);

	//\ 2022.9.27
	bool PairEtaPhiSelect(AliFemtoPair* aPair);
	void SetUsePairCutEtaPhi(int aUsePairCutEtaPhi);
	void SetPairCutEtaPhi(float aEtaCut,float aPhiCut);
	
	void SetPassAllPair(int aUse);
	void SetGobalVelGate(int aUse);

	void SetIsSameParticlePair(int aUse);
	float CalcMt(const AliFemtoPair* aPair);
	void SetUse3DkTvsKStarvsmT(int aUse);
	void Set3DkTvsKStarvsmTInit(bool aInit,
		int nbinskT,float lowkT,float upkT,
		int nbinsks,float lowks,float upks,
		int nbinsmT,float lowmT,float upmT);


	void Set2DKstarVsmT(int aUse);
	void Set2DkStarVsmTInit(bool aInit,
		int nbinsks,float lowks,float upks,
                int nbinsmT,float lowmT,float upmT);

void SetdBumpCheck(int aUse);
void SetdBumpCheckInit(bool aInit,
 int nbinsks,float lowks,float upks,
int nbinspT,float lowpT,float uppT,
 int nbinsMass,float lowMass,float upMass);


void SetmTLimitCheck(int aUse);

void SetmTCheckInit(bool aInit,
 int nbinsks,float lowks,float upks,
int nbinspT,float lowpT,float uppT,
int nbinsmT,float lowmT,float upmT);

        virtual AliFemtoCorrFctnpdtHe3* Clone() const  { return new AliFemtoCorrFctnpdtHe3(*this); }
    protected:
        int isHe3Pair;
        TString fTitle;
        int fNbinsKStar;
        double fKStarLow, fKStarHigh;
        float fP1Mass, fP2Mass;
        TH1D* fNumerator;          // numerator - real pairs
        TH1D* fDenominator;        // denominator - mixed pairs

        // velocity gate 11.23
	int fUseVelGate;
        TH1D *fP1EarlierP2Num;
        TH1D *fP1EarlierP2Dum;
        TH1D *fP2EarlierP1Num;
        TH1D *fP2EarlierP1Dum;
 
        TH3F *fNumHigh3F;
        TH3F *fDenHigh3F;
	// control highD hist & mass side band
	bool fHighCF;
	bool fSideBand;
	// TF1 for TOF mass cut band!
	
	TF1 *p2Up3Sigma;
	TF1 *p2Low3Sigma;

	TF1 *SideBandUp;
	TF1 *SideBandLow;
	
	// 8 area for mass side band!
	//\----------------
	// num!
	//\----------------
	// cut by 3 sigma
	TH1D *A1SideBandNum;
	TH1D *S1SideBandNum;
	TH1D *A2SideBandNum;
	// merge
	TH1D *A1A2SideBandNum;
	TH1D *SignalAndSideCFNum;
	

	//\----------------	
	//\ dum!
	//\----------------
	// cut by 3 sigma
	TH1D *A1SideBandDum;
	TH1D *S1SideBandDum;
	TH1D *A2SideBandDum;
	// merge
	TH1D *A1A2SideBandDum;
	TH1D *SignalAndSideCFDum;

	int fUsePt;
	int fUseDPhiDEtaQA;
	//TH2F *fNumDPhiDEtaQAFailCut;
	//TH2F *fDumDPhiDEtaQAFailCut;
	TH2F *fNumDPhiDEtaAvgQA;
	TH2F *fDumDPhiDEtaAvgQA;
		
	int fUseStavinskyMethod;
	TH1F *fStaSkyBkg;
	

	int fUsePairCutEtaPhi;
	float fPairCut_eta;
	float fPairCut_phi;
	
	int fPassAllPair;

	int fUseGobalVelGate;

	int fUse2DpTvsKStar;
	int IsSameParticlePair;
	TH2F *KStarVspT_P1Hist;
	TH2F *KStarVspT_P2Hist;

	int fUse3DkTvsKStarvsmT;
	TH3F *fNum3DkTvsKStarvsmT;
	TH3F *fDum3DkTvsKStarvsmT;
	
	int fUse2DkStarVsmT;
	TH2F *fNum2DkStarVsmT;
	TH2F *fDum2DkStarVsmT;

	int fUseBumpC;
	TH2F *f2DkSVspT;
	TH2F *f2DkSVsMass;
	
	int fUsemTCheck;
	TH3F *f3DmTDepkSVspT;
};

#endif


