///
/// \file AliFemtoCorrFctnpdtHe3.cxx
///

#include "AliFemtoCorrFctnpdtHe3.h"
const float ProtonMass = 0.9383;
const float DeuteronMass= 1.8756;
const float TritonMass = 2.8089;
const float He3Mass = 2.8084;
#define PIH 1.57079632679489656
#define PIT 6.28318530717958623
#define fphiL -1.4137167	//default 20bin for phi&eta
#define fphiT 4.8694686
static float TPCradii[9] = { 0.85, 1.05, 1.25, 1.45, 1.65, 1.85, 2.05, 2.25, 2.45 };
/*
AliFemtoCorrFctnpdtHe3::AliFemtoCorrFctnpdtHe3():
    AliFemtoCorrFctnpdtHe3("CorrFctnKStar", 200, 0, 1)
{
// no-op
}
*/
AliFemtoCorrFctnpdtHe3::AliFemtoCorrFctnpdtHe3(const char* title,
                                             const int nbins,
                                             const float KStarLo,
                                             const float KStarHi):
    AliFemtoCorrFctn(),
    isHe3Pair(0),
    fTitle(title),
    fNbinsKStar(nbins),
    fKStarLow(KStarLo),
    fKStarHigh(KStarHi),
    fP1Mass(0.1),
    fP2Mass(0.1),
    fNumerator(nullptr),
    fDenominator(nullptr),
    fUseVelGate(0),
    fP1EarlierP2Num(nullptr),
    fP1EarlierP2Dum(nullptr),
    fP2EarlierP1Num(nullptr),
    fP2EarlierP1Dum(nullptr),
    fNumHigh3F(nullptr),
    fDenHigh3F(nullptr),
    fHighCF(false),
    fSideBand(false),
    p2Up3Sigma(nullptr),
    p2Low3Sigma(nullptr),
    SideBandUp(nullptr),
    SideBandLow(nullptr),
    //
    A1SideBandNum(nullptr),
    S1SideBandNum(nullptr),
    A2SideBandNum(nullptr),
    A1A2SideBandNum(nullptr),
    SignalAndSideCFNum(nullptr),
    A1SideBandDum(nullptr),
    S1SideBandDum(nullptr),
    A2SideBandDum(nullptr),
    A1A2SideBandDum(nullptr),
    SignalAndSideCFDum(nullptr),
    fUsePt(0),
    fUseDPhiDEtaQA(0),
    //fNumDPhiDEtaQAFailCut(nullptr),
    //fDumDPhiDEtaQAFailCut(nullptr),
    fNumDPhiDEtaAvgQA(nullptr),
    fDumDPhiDEtaAvgQA(nullptr),
    fUseStavinskyMethod(0),
    fStaSkyBkg(nullptr),
    fUsePairCutEtaPhi(0),
    fPairCut_eta(0.017),
    fPairCut_phi(0.017),
    fPassAllPair(0),
    fUseGobalVelGate(0),
    fUse2DpTvsKStar(0),
    IsSameParticlePair(0),
    KStarVspT_P1Hist(nullptr),
    KStarVspT_P2Hist(nullptr),
    fUse3DkTvsKStarvsmT(0),
    fNum3DkTvsKStarvsmT(nullptr),
    fDum3DkTvsKStarvsmT(nullptr), 
    fUse2DkStarVsmT(0),
    fNum2DkStarVsmT(nullptr),
    fDum2DkStarVsmT(nullptr),
    fUseBumpC(0),
    f2DkSVspT(nullptr),
    f2DkSVsMass(nullptr),
    fUsemTCheck(0),
    f3DmTDepkSVspT(nullptr)
{
    
    fNumerator      = new TH1D(TString::Format("Num%s", fTitle.Data()), "fNumerator", nbins, KStarLo, KStarHi);
    fDenominator    = new TH1D(TString::Format("Dum%s", fTitle.Data()), "fDenominator", nbins, KStarLo, KStarHi);
    
    fNumerator->Sumw2();
    fDenominator->Sumw2();

	
}
AliFemtoCorrFctnpdtHe3::AliFemtoCorrFctnpdtHe3(const AliFemtoCorrFctnpdtHe3& aCorrFctn):
    AliFemtoCorrFctn(aCorrFctn),
    isHe3Pair(aCorrFctn.isHe3Pair),
    fTitle(aCorrFctn.fTitle),
    fNbinsKStar(aCorrFctn.fNbinsKStar),
    fKStarLow(aCorrFctn.fKStarLow),
    fKStarHigh(aCorrFctn.fKStarHigh),    
    fP1Mass(aCorrFctn.fP1Mass),
    fP2Mass(aCorrFctn.fP2Mass),
    fNumerator(aCorrFctn.fNumerator ? new TH1D(*aCorrFctn.fNumerator) : nullptr),
    fDenominator(aCorrFctn.fDenominator ? new TH1D(*aCorrFctn.fDenominator) : nullptr),
    fUseVelGate(aCorrFctn.fUseVelGate),
    fP1EarlierP2Num(aCorrFctn.fP1EarlierP2Num),
    fP1EarlierP2Dum(aCorrFctn.fP1EarlierP2Dum),
    fP2EarlierP1Num(aCorrFctn.fP2EarlierP1Num),
    fP2EarlierP1Dum(aCorrFctn.fP2EarlierP1Dum),
    fNumHigh3F(aCorrFctn.fNumHigh3F),
    fDenHigh3F(aCorrFctn.fDenHigh3F),
    fHighCF(aCorrFctn.fHighCF),
    fSideBand(aCorrFctn.fSideBand),
    p2Up3Sigma(aCorrFctn.p2Up3Sigma),
    p2Low3Sigma(aCorrFctn.p2Low3Sigma),
    SideBandUp(aCorrFctn.SideBandUp),
    SideBandLow(aCorrFctn.SideBandLow),
    A1SideBandNum(aCorrFctn.A1SideBandNum),
    S1SideBandNum(aCorrFctn.S1SideBandNum),
    A2SideBandNum(aCorrFctn.A2SideBandNum),
    A1A2SideBandNum(aCorrFctn.A1A2SideBandNum),
    SignalAndSideCFNum(aCorrFctn.SignalAndSideCFNum),
    A1SideBandDum(aCorrFctn.A1SideBandDum),
    S1SideBandDum(aCorrFctn.S1SideBandDum),
    A2SideBandDum(aCorrFctn.A2SideBandDum),
    A1A2SideBandDum(aCorrFctn.A1A2SideBandDum),
    SignalAndSideCFDum(aCorrFctn.SignalAndSideCFDum),
    fUsePt(aCorrFctn.fUsePt),
    fUseDPhiDEtaQA(aCorrFctn.fUseDPhiDEtaQA),
    //fNumDPhiDEtaQAFailCut(aCorrFctn.fNumDPhiDEtaQAFailCut),
    //fDumDPhiDEtaQAFailCut(aCorrFctn.fDumDPhiDEtaQAFailCut),
    fNumDPhiDEtaAvgQA(aCorrFctn.fNumDPhiDEtaAvgQA),
    fDumDPhiDEtaAvgQA(aCorrFctn.fDumDPhiDEtaAvgQA),
    fUseStavinskyMethod(aCorrFctn.fUseStavinskyMethod),
    fStaSkyBkg(aCorrFctn.fStaSkyBkg),
    fUsePairCutEtaPhi(aCorrFctn.fUsePairCutEtaPhi),
    fPairCut_eta(aCorrFctn.fPairCut_eta),
    fPairCut_phi(aCorrFctn.fPairCut_phi),
    fPassAllPair(aCorrFctn.fPassAllPair),
    fUseGobalVelGate(aCorrFctn.fUseGobalVelGate),
    fUse2DpTvsKStar(aCorrFctn.fUse2DpTvsKStar),
    IsSameParticlePair(aCorrFctn.IsSameParticlePair),
    KStarVspT_P1Hist(aCorrFctn.KStarVspT_P1Hist),
    KStarVspT_P2Hist(aCorrFctn.KStarVspT_P2Hist),
    fUse3DkTvsKStarvsmT(aCorrFctn.fUse3DkTvsKStarvsmT),
    fNum3DkTvsKStarvsmT(aCorrFctn.fNum3DkTvsKStarvsmT),
    fDum3DkTvsKStarvsmT(aCorrFctn.fDum3DkTvsKStarvsmT),
    fUse2DkStarVsmT(aCorrFctn.fUse2DkStarVsmT),
    fNum2DkStarVsmT(aCorrFctn.fNum2DkStarVsmT),
    fDum2DkStarVsmT(aCorrFctn.fDum2DkStarVsmT),
    fUseBumpC(aCorrFctn.fUseBumpC),
    f2DkSVspT(aCorrFctn.f2DkSVspT),
    f2DkSVsMass(aCorrFctn.f2DkSVsMass),
    fUsemTCheck(aCorrFctn.fUsemTCheck),
    f3DmTDepkSVspT(aCorrFctn.f3DmTDepkSVspT)
{
    


}
AliFemtoCorrFctnpdtHe3::~AliFemtoCorrFctnpdtHe3()
{
    // destructor
    delete fNumerator;
    delete fDenominator;
    
    delete fP1EarlierP2Num;
    delete fP1EarlierP2Dum;
    delete fP2EarlierP1Num;
    delete fP2EarlierP1Dum;
    
    delete fNumHigh3F;
    delete fDenHigh3F;

    delete p2Up3Sigma;
    delete p2Low3Sigma;
    delete SideBandUp;
    delete SideBandLow;

    delete A1SideBandNum;
    delete S1SideBandNum;
    delete A2SideBandNum;
    delete A1A2SideBandNum;
    delete SignalAndSideCFNum;

    delete A1SideBandDum;
    delete S1SideBandDum;
    delete A2SideBandDum;
    delete A1A2SideBandDum;
    delete SignalAndSideCFDum;
    
    //delete fNumDPhiDEtaQAFailCut;
    //delete fDumDPhiDEtaQAFailCut;    
    delete fNumDPhiDEtaAvgQA;
    delete fDumDPhiDEtaAvgQA;

    delete fStaSkyBkg;

    delete KStarVspT_P1Hist;
    delete KStarVspT_P2Hist;

    delete fNum3DkTvsKStarvsmT;
    delete fDum3DkTvsKStarvsmT;

    delete fNum2DkStarVsmT;
    delete fDum2DkStarVsmT;

    delete f2DkSVspT;
    delete f2DkSVsMass;


    delete f3DmTDepkSVspT;
}
AliFemtoCorrFctnpdtHe3& AliFemtoCorrFctnpdtHe3::operator=(const AliFemtoCorrFctnpdtHe3& aCorrFctn)
{
    // assignment operator
    if (this == &aCorrFctn) {
        return *this;
    }

    AliFemtoCorrFctnpdtHe3::operator=(aCorrFctn);
    isHe3Pair = aCorrFctn.isHe3Pair;
    fTitle = aCorrFctn.fTitle;
    fNbinsKStar = aCorrFctn.fNbinsKStar;
    fKStarLow = aCorrFctn.fKStarLow;
    fKStarHigh = aCorrFctn.fKStarHigh;
    fP1Mass = aCorrFctn.fP1Mass;
    fP2Mass = aCorrFctn.fP2Mass;

    if(fNumerator) delete fNumerator;
        fNumerator = new TH1D(*aCorrFctn.fNumerator);
    if(fDenominator) delete fDenominator;
        fDenominator = new TH1D(*aCorrFctn.fDenominator);
 
    fUseVelGate = aCorrFctn.fUseVelGate;

    if(fP1EarlierP2Num) delete fP1EarlierP2Num;
        fP1EarlierP2Num = new TH1D(*aCorrFctn.fP1EarlierP2Num);
    if(fP1EarlierP2Dum) delete fP1EarlierP2Dum;
        fP1EarlierP2Dum = new TH1D(*aCorrFctn.fP1EarlierP2Dum);
    if(fP2EarlierP1Num) delete fP2EarlierP1Num;
        fP2EarlierP1Num = new TH1D(*aCorrFctn.fP2EarlierP1Num);
    if(fP2EarlierP1Dum) delete fP2EarlierP1Dum;
        fP2EarlierP1Dum = new TH1D(*aCorrFctn.fP2EarlierP1Dum);

    fHighCF 	= aCorrFctn.fHighCF;
    fSideBand 	= aCorrFctn.fSideBand;

    if(fNumHigh3F) delete fNumHigh3F;
    	fNumHigh3F = new TH3F(*aCorrFctn.fNumHigh3F);
    if(fDenHigh3F) delete fDenHigh3F;
    	fDenHigh3F = new TH3F(*aCorrFctn.fDenHigh3F);

    if(p2Up3Sigma) delete p2Up3Sigma;
	p2Up3Sigma = new TF1(*aCorrFctn.p2Up3Sigma);
    if(p2Low3Sigma) delete p2Low3Sigma;
	p2Low3Sigma = new TF1(*aCorrFctn.p2Low3Sigma);
    if(SideBandUp) delete SideBandUp;
	SideBandUp = new TF1(*aCorrFctn.SideBandUp);
    if(SideBandLow) delete SideBandLow;
	SideBandLow = new TF1(*aCorrFctn.SideBandLow);

	if(A1SideBandNum) delete A1SideBandNum;
        		A1SideBandNum = new TH1D(*aCorrFctn.A1SideBandNum);
	if(S1SideBandNum) delete S1SideBandNum;
			S1SideBandNum = new TH1D(*aCorrFctn.S1SideBandNum);
	if(A2SideBandNum) delete A2SideBandNum;
			A2SideBandNum = new TH1D(*aCorrFctn.A2SideBandNum);
	if(A1A2SideBandNum) delete A1A2SideBandNum;
			A1A2SideBandNum = new TH1D(*aCorrFctn.A1A2SideBandNum);
	if(SignalAndSideCFNum) delete SignalAndSideCFNum;
			SignalAndSideCFNum = new TH1D(*aCorrFctn.SignalAndSideCFNum);
	if(A1SideBandDum) delete A1SideBandDum;
			A1SideBandDum 	= new TH1D(*aCorrFctn.A1SideBandDum);
	if(S1SideBandDum) delete S1SideBandDum;
			S1SideBandDum 	= new TH1D(*aCorrFctn.S1SideBandDum);
	if(A2SideBandDum) delete A2SideBandDum;
			A2SideBandDum 	= new TH1D(*aCorrFctn.A2SideBandDum);
	if(A1A2SideBandDum) delete A1A2SideBandDum;
			A1A2SideBandDum = new TH1D(*aCorrFctn.A1A2SideBandDum);
	if(SignalAndSideCFDum) delete SignalAndSideCFDum;
			SignalAndSideCFDum = new TH1D(*aCorrFctn.SignalAndSideCFDum);
    fUsePt = aCorrFctn.fUsePt;
    fUseDPhiDEtaQA = aCorrFctn.fUseDPhiDEtaQA;
/*   
if(fNumDPhiDEtaQAFailCut) delete fNumDPhiDEtaQAFailCut;
        		fNumDPhiDEtaQAFailCut = new TH2F(*aCorrFctn.fNumDPhiDEtaQAFailCut);
    if(fDumDPhiDEtaQAFailCut) delete fDumDPhiDEtaQAFailCut;
        		fDumDPhiDEtaQAFailCut = new TH2F(*aCorrFctn.fDumDPhiDEtaQAFailCut);
*/
	if(fNumDPhiDEtaAvgQA) delete fNumDPhiDEtaAvgQA;
        		fNumDPhiDEtaAvgQA = new TH2F(*aCorrFctn.fNumDPhiDEtaAvgQA);
    	if(fDumDPhiDEtaAvgQA) delete fDumDPhiDEtaAvgQA;
        		fDumDPhiDEtaAvgQA = new TH2F(*aCorrFctn.fDumDPhiDEtaAvgQA);

    if(fUseStavinskyMethod) delete fStaSkyBkg;
		fStaSkyBkg = new TH1F(*aCorrFctn.fStaSkyBkg);


     if(KStarVspT_P1Hist) delete KStarVspT_P1Hist;
		KStarVspT_P1Hist = new TH2F(*aCorrFctn.KStarVspT_P1Hist);

	if(KStarVspT_P2Hist) delete KStarVspT_P2Hist;
                KStarVspT_P2Hist = new TH2F(*aCorrFctn.KStarVspT_P2Hist);

	if(fNum3DkTvsKStarvsmT) delete fNum3DkTvsKStarvsmT;
		fNum3DkTvsKStarvsmT = new TH3F(*aCorrFctn.fNum3DkTvsKStarvsmT);
	
	if(fDum3DkTvsKStarvsmT) delete fDum3DkTvsKStarvsmT;
                fDum3DkTvsKStarvsmT = new TH3F(*aCorrFctn.fDum3DkTvsKStarvsmT);

if(fNum2DkStarVsmT) delete fNum2DkStarVsmT;
        fNum2DkStarVsmT = new TH2F(*aCorrFctn.fNum2DkStarVsmT);

if(fDum2DkStarVsmT) delete fDum2DkStarVsmT;
        fDum2DkStarVsmT = new TH2F(*aCorrFctn.fDum2DkStarVsmT);

if(fUseBumpC){
	delete f2DkSVspT;
	f2DkSVspT = new TH2F(*aCorrFctn.f2DkSVspT);
	delete f2DkSVsMass;
	f2DkSVsMass = new TH2F(*aCorrFctn.f2DkSVsMass); 
}

if(fUsemTCheck){

	delete f3DmTDepkSVspT;
	f3DmTDepkSVspT = new TH3F(*aCorrFctn.f3DmTDepkSVspT);
}
    return *this;

}

AliFemtoString AliFemtoCorrFctnpdtHe3::Report()
{
    TString report = "AliFemtoCorrFctnKStar:\n";
    report += TString::Format("Number of entries in numerator:\t%E\n", fNumerator->GetEntries());
    report += TString::Format("Number of entries in denominator:\t%E\n", fDenominator->GetEntries());
    return AliFemtoString((const char *)report);

}
TList* AliFemtoCorrFctnpdtHe3::GetOutputList()
{
    // Prepare the list of objects to be written to the output
    TList *tOutputList = new TList();

    tOutputList->Add(fNumerator);
    tOutputList->Add(fDenominator);
	if(fUseVelGate){
	    tOutputList->Add(fP1EarlierP2Num);
	    tOutputList->Add(fP1EarlierP2Dum);
	    tOutputList->Add(fP2EarlierP1Num);
	    tOutputList->Add(fP2EarlierP1Dum);
	}
    if(fHighCF){
    	tOutputList->Add(fNumHigh3F);
    	tOutputList->Add(fDenHigh3F);
    }
    if(fSideBand){
	    tOutputList->Add(A1SideBandNum);
	    tOutputList->Add(S1SideBandNum);
	    tOutputList->Add(A2SideBandNum);
	    tOutputList->Add(A1A2SideBandNum);
	    tOutputList->Add(SignalAndSideCFNum);
	    
	    tOutputList->Add(A1SideBandDum);
	    tOutputList->Add(S1SideBandDum);
	    tOutputList->Add(A2SideBandDum);
	    tOutputList->Add(A1A2SideBandDum);
	    tOutputList->Add(SignalAndSideCFDum);
    }
    if(fUseDPhiDEtaQA){
            //tOutputList->Add(fNumDPhiDEtaQAFailCut);
	    //tOutputList->Add(fDumDPhiDEtaQAFailCut);
  	    tOutputList->Add(fNumDPhiDEtaAvgQA);
	    tOutputList->Add(fDumDPhiDEtaAvgQA);
    }
    if(fUseStavinskyMethod) tOutputList->Add(fStaSkyBkg);
    
   if(fUse2DpTvsKStar)
	{
	tOutputList->Add(KStarVspT_P1Hist);
	if(IsSameParticlePair!=1) tOutputList->Add(KStarVspT_P2Hist);
	}

	if(fUse3DkTvsKStarvsmT){
		tOutputList->Add(fNum3DkTvsKStarvsmT);
		tOutputList->Add(fDum3DkTvsKStarvsmT);
	}

if(fUse2DkStarVsmT){
	tOutputList->Add(fNum2DkStarVsmT);
tOutputList->Add(fDum2DkStarVsmT);
}
 
if(fUseBumpC){

tOutputList->Add(f2DkSVspT);
tOutputList->Add(f2DkSVsMass);
}


if(fUsemTCheck){
tOutputList->Add(f3DmTDepkSVspT);
}
    return tOutputList;
}
void AliFemtoCorrFctnpdtHe3::Finish()
{
    //cout<<"AliFemtoCorrFctnpdtHe3::Finish()"<<endl;
}
void AliFemtoCorrFctnpdtHe3::Write()
{
    // Write out neccessary objects
    fNumerator->Write();
    fDenominator->Write();
	if(fUseVelGate){
	    fP1EarlierP2Num->Write();
	    fP1EarlierP2Dum->Write();
	    fP2EarlierP1Num->Write();
	    fP2EarlierP1Dum->Write();
	}
 	if(fHighCF){
	    fNumHigh3F->Write();
	    fDenHigh3F->Write();
	}
	if(fSideBand){
	    A1SideBandNum->Write();
	    S1SideBandNum->Write();
	    A2SideBandNum->Write();
	    A1A2SideBandNum->Write();
	    SignalAndSideCFNum->Write();

	    A1SideBandDum->Write();
	    S1SideBandDum->Write();
	    A2SideBandDum->Write();
	    A1A2SideBandDum->Write();
	    SignalAndSideCFDum->Write();
	    
        }
        if(fUseDPhiDEtaQA){
		  //fNumDPhiDEtaQAFailCut->Write();
		  //fDumDPhiDEtaQAFailCut->Write();
		  fNumDPhiDEtaAvgQA->Write();
		  fDumDPhiDEtaAvgQA->Write();
	}
        if(fUseStavinskyMethod) fStaSkyBkg->Write();
	if(fUse2DpTvsKStar){
	  KStarVspT_P1Hist->Write();
          if(IsSameParticlePair!=1) KStarVspT_P2Hist->Write();

	}

	 if(fUse3DkTvsKStarvsmT){
                fNum3DkTvsKStarvsmT->Write();
                fDum3DkTvsKStarvsmT->Write();
        }
if(fUse2DkStarVsmT){
        fNum2DkStarVsmT->Write();
	fDum2DkStarVsmT->Write();
}

if(fUseBumpC){
	f2DkSVspT->Write();
	f2DkSVsMass->Write();
}


if(fUsemTCheck){

f3DmTDepkSVspT->Write();
}
}
void AliFemtoCorrFctnpdtHe3::AddRealPair(AliFemtoPair* aPair){

     // change momentum of p2!
          AliFemtoPair* fPair = new AliFemtoPair;
     //
	if(isHe3Pair){
    // change momentum of p2!
    AliFemtoPair* fPair = new AliFemtoPair;
        fPair = ChangeP2Mom(aPair);
    }
    else{
        fPair = aPair;
    }
    
    // add true pair
 double tKStar = fabs(fPair->KStar());
 if (fPairCut && !fPairCut->Pass(fPair)) {
       // failed pair QA
        if(fUseDPhiDEtaQA>1 && tKStar<0.2){

          double eta1 = fPair->Track1()->FourMomentum().PseudoRapidity();
          double eta2 = fPair->Track2()->FourMomentum().PseudoRapidity();
          float AvgDPhi = ReAvgDphi(fPair);
          double deta = eta1 - eta2;

          fDumDPhiDEtaAvgQA->Fill(deta,AvgDPhi);

        }

 
        return;
    }
   if(fUsePairCutEtaPhi){
		    if(!PairEtaPhiSelect(fPair)) return;
	}
/*
	if(fPassAllPair==0){
	    if (fPairCut && !fPairCut->Pass(fPair)) {
		if(fUseDPhiDEtaQA==2 && tKStar<0.2){

          double eta1 = fPair->Track1()->FourMomentum().PseudoRapidity();
          double eta2 = fPair->Track2()->FourMomentum().PseudoRapidity();
          float AvgDPhi = ReAvgDphi(fPair);
          double deta = eta1 - eta2;

          fNumDPhiDEtaAvgQA->Fill(deta,AvgDPhi);

        }

		return;
	    }
	    
	    if(fUsePairCutEtaPhi){
		if(!PairEtaPhiSelect(fPair)) return;
	    }
	}

*/
    if(fUseGobalVelGate){
	int VelLabel = ReVelocityGate(fPair);
	if(fUseGobalVelGate == VelLabel){
	}
	else{
		return;
	}

    }
    fNumerator->Fill(tKStar);

    if(fHighCF){
	float FillMom1 = 0.;
	float FillMom2 = 0.;
	if(fUsePt){
		FillMom1 = fPair->Track1()->Track()->Pt();
		FillMom2 = fPair->Track2()->Track()->Pt();
	}
	else{
		FillMom1 = fPair->Track1()->Track()->P().Mag();
		FillMom2 = fPair->Track2()->Track()->P().Mag();
	}
	fNumHigh3F->Fill(FillMom1,FillMom2,tKStar); 
	}
    if(fSideBand){
	FillSideBandNum(fPair);
	}
 
   if(fUseDPhiDEtaQA>0 && tKStar<0.2){
	  double eta1 = fPair->Track1()->FourMomentum().PseudoRapidity();
    	  double eta2 = fPair->Track2()->FourMomentum().PseudoRapidity();
    	  float AvgDPhi = ReAvgDphi(fPair);
	  double deta = eta1 - eta2;

 	  //fNumDPhiDEtaQA->Fill(deta,AvgDPhi);
	  fNumDPhiDEtaAvgQA->Fill(deta,AvgDPhi);
	
	}

	if(fUseVelGate){
	    int VelLabel = ReVelocityGate(fPair);
	    if(VelLabel == 1){
		fP1EarlierP2Num->Fill(tKStar);
	    }
	    else if(VelLabel == 2){
		fP2EarlierP1Num->Fill(tKStar);
	    }
	    else if(VelLabel == 3){
		return;
	    }
	}

      if(fUseStavinskyMethod){
		bool passpair = false;
		AliFemtoPair* SSPair = new AliFemtoPair;
		SSPair = InversePair(fPair);
		if(fPairCut->Pass(SSPair)){
			passpair = true;
			if(fUsePairCutEtaPhi && !PairEtaPhiSelect(fPair)){
				passpair = false;	
			}
		}
		if(passpair){
			float InverseKStar = fabs(SSPair->KStar());
			fStaSkyBkg->Fill(InverseKStar);
		}
	}
     	if(fUse3DkTvsKStarvsmT){
                fNum3DkTvsKStarvsmT->Fill(fPair->KT(),tKStar,CalcMt(fPair));
        }
if(fUse2DkStarVsmT){
        fNum2DkStarVsmT->Fill(tKStar,CalcMt(fPair));
}


	if(fUseBumpC && 
	tKStar>0.18 && tKStar<0.25 &&
	CalcMt(fPair) > 2.0 && CalcMt(fPair)<3.5
	){
		f2DkSVspT->Fill(tKStar,fPair->Track2()->Track()->Pt());
 	 float c = 1.;
    	float beta = fPair->Track2()->Track()->VTOF();
    	if(beta!=0){

    		float tMom = fPair->Track2()->Track()->P().Mag();
    		float massTOF = tMom*tMom/c/c*(1./(beta*beta)-1);
        	f2DkSVsMass->Fill(tKStar,massTOF);
	}
 }

	return;
    
 
}

void AliFemtoCorrFctnpdtHe3::AddMixedPair(AliFemtoPair* aPair)
{
    // change momentum of p2!
    AliFemtoPair* fPair = new AliFemtoPair;

    if(isHe3Pair){
        fPair = ChangeP2Mom(aPair);
    }
    else{
        fPair = aPair;
    }

    // add true pair
double tKStar = fabs(fPair->KStar());

   if (fPairCut && !fPairCut->Pass(fPair)) {

        return;
    }
    if(fUsePairCutEtaPhi){
		    if(!PairEtaPhiSelect(fPair)) return;
	}
 

if(fUseGobalVelGate){
        int VelLabel = ReVelocityGate(fPair);
        if(fUseGobalVelGate == VelLabel){
        }
        else{
                return;
        }

    }

    fDenominator->Fill(tKStar);
    if(fHighCF){
	float FillMom1 = 0.;
	float FillMom2 = 0.;
	if(fUsePt){
		FillMom1 = fPair->Track1()->Track()->Pt();
		FillMom2 = fPair->Track2()->Track()->Pt();
	}
	else{
		FillMom1 = fPair->Track1()->Track()->P().Mag();
		FillMom2 = fPair->Track2()->Track()->P().Mag();
	}
	fDenHigh3F->Fill(FillMom1,FillMom2,tKStar);
    }
    if(fSideBand){
	FillSideBandDum(fPair);
	}
	if(fUseVelGate){
	    int VelLabel = ReVelocityGate(fPair);
	    if(VelLabel == 1){
		fP1EarlierP2Dum->Fill(tKStar);
	    }
	    else if(VelLabel == 2){
		fP2EarlierP1Dum->Fill(tKStar);
	    }
	    else if(VelLabel == 3){
		return;
	    }
	}
	if(fUse2DpTvsKStar){
			KStarVspT_P1Hist->Fill(tKStar,fPair->Track1()->Track()->Pt());
          		if(IsSameParticlePair!=1) KStarVspT_P2Hist->Fill(tKStar,fPair->Track2()->Track()->Pt());
        }

	if(fUse3DkTvsKStarvsmT){
                fDum3DkTvsKStarvsmT->Fill(fPair->KT(),tKStar,CalcMt(fPair));
	 }
if(fUse2DkStarVsmT){
fDum2DkStarVsmT->Fill(tKStar,CalcMt(fPair));
}

	if(fUsemTCheck){
		if(fUsemTCheck==1){
			 f3DmTDepkSVspT->Fill(tKStar,fPair->Track2()->Track()->Pt(),CalcMt(fPair));
			f3DmTDepkSVspT->Fill(tKStar,fPair->Track1()->Track()->Pt(),CalcMt(fPair));
		}
		if(fUsemTCheck==2){
			f3DmTDepkSVspT->Fill(tKStar,fPair->Track2()->Track()->Pt(),fPair->KT());
			f3DmTDepkSVspT->Fill(tKStar,fPair->Track1()->Track()->Pt(),fPair->KT());
		}
	}
	return;
    

}
void AliFemtoCorrFctnpdtHe3::SetIsHe3Pair(int isOrNot)
{
    // 0 is not He3 pair, 1 is 
    isHe3Pair = isOrNot;
}
AliFemtoPair * AliFemtoCorrFctnpdtHe3::ChangeP2Mom(AliFemtoPair* aPair)
{
    AliFemtoPair* fPair = new AliFemtoPair;

    AliFemtoParticle *tPart1 = new AliFemtoParticle(*aPair->Track1());
    fPair->SetTrack1(tPart1);

    // modify momentum!
    AliFemtoParticle *tPart2 = new AliFemtoParticle(*aPair->Track2());
    AliFemtoLorentzVector tFourMom2 = AliFemtoLorentzVector(tPart2->FourMomentum());
    tFourMom2.SetPx(2.*tFourMom2.px());
    tFourMom2.SetPy(2.*tFourMom2.py());
    tFourMom2.SetPz(2.*tFourMom2.pz());

    float TotalP = tFourMom2.px() * tFourMom2.px() + tFourMom2.py() * tFourMom2.py() + tFourMom2.pz() * tFourMom2.pz();
    float TotalE = TMath::Sqrt( He3Mass * He3Mass + TotalP);
    tFourMom2.SetE(TotalE);
    tPart2->ResetFourMomentum(tFourMom2);

    fPair->SetTrack2(tPart2);

    return fPair;
}
int AliFemtoCorrFctnpdtHe3::ReVelocityGate(AliFemtoPair* aPair){
    //\ 1: p1 faster than p2
    //\ 2: p2 faster than p1
    //\ 3: bad
    AliFemtoParticle *tPart1 = new AliFemtoParticle(*aPair->Track1());
    AliFemtoLorentzVector tFourMom1 = AliFemtoLorentzVector(tPart1->FourMomentum());
    float TotalP = tFourMom1.px() * tFourMom1.px() + tFourMom1.py() * tFourMom1.py() + tFourMom1.pz() * tFourMom1.pz();
    float TotalE = TMath::Sqrt( fP1Mass * fP1Mass + TotalP);
    TotalP = TMath::Sqrt( TotalP);
    float P1velocity = TotalP/TotalE;

    //
   TotalP = 0.;
    TotalE = 0.;   
 
    AliFemtoParticle *tPart2 = new AliFemtoParticle(*aPair->Track2());
    AliFemtoLorentzVector tFourMom2 = AliFemtoLorentzVector(tPart2->FourMomentum());
    TotalP = tFourMom2.px() * tFourMom2.px() + tFourMom2.py() * tFourMom2.py() + tFourMom2.pz() * tFourMom2.pz();
    TotalE = TMath::Sqrt( fP2Mass * fP2Mass + TotalP);
    TotalP = TMath::Sqrt( TotalP);
    float P2velocity = TotalP/TotalE;

    int ReLabel = 3;
    if(P1velocity > P2velocity){
        ReLabel = 1;
    }
    if(P1velocity < P2velocity){
        ReLabel = 2;
    }
    return ReLabel;
}
float AliFemtoCorrFctnpdtHe3::ReAvgDphi(AliFemtoPair* aPair){
	
  AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  Double_t magsign = 0.0;
  if (!aodH) {
   return -999;
  }
  else {
    AliAODEvent *fAOD;
    fAOD = aodH->GetEvent();
    magsign = fAOD->GetMagneticField();
  }
  float fMagSign = 1.;
  if (magsign > 1)
    fMagSign = 1.;
  else if ( magsign < 1)
    fMagSign = -1.;
  else
    fMagSign = magsign;
  
  float magval = 0.5;
  float dphiAvg = 0.;
  double phi1 = aPair->Track1()->Track()->P().Phi();
  double phi2 = aPair->Track2()->Track()->P().Phi();
  double chg1 = aPair->Track1()->Track()->Charge();
  double chg2 = aPair->Track2()->Track()->Charge();
  double pt1 = aPair->Track1()->Track()->Pt();
  double pt2 = aPair->Track2()->Track()->Pt();
  for(int i=0;i<9;i++){
     Double_t rad = TPCradii[i];

    // Calculate dPhiStar:
    //double afsi0b = -0.07510020733*chg1*fMagSign*rad/pt1;
    //double afsi1b = -0.07510020733*chg2*fMagSign*rad/pt2;
    double afsi0b = -0.15*magval*chg1*fMagSign*rad/pt1;
    double afsi1b = -0.15*magval*chg2*fMagSign*rad/pt2;
    Double_t dphistar =  phi2 - phi1 + TMath::ASin(afsi1b) - TMath::ASin(afsi0b);
    dphistar = TVector2::Phi_mpi_pi(dphistar); // returns phi angle in the interval [-PI,PI)
    dphiAvg += dphistar;
   }
   dphiAvg = dphiAvg/9.;
   return dphiAvg; 
}
void AliFemtoCorrFctnpdtHe3::SetP1AndP2Mass(float p1Mass,float p2Mass){

    fP1Mass = p1Mass;
    fP2Mass = p2Mass;
}
void AliFemtoCorrFctnpdtHe3::SetfUseVelGate(int aUseVelGate){
	fUseVelGate = aUseVelGate;
}
void AliFemtoCorrFctnpdtHe3::SetVelGateInit(bool aUseVelGate){
	if(aUseVelGate){
	    fP1EarlierP2Num = new TH1D(TString::Format("p1Ep2Num_%s", fTitle.Data()), " ", fNbinsKStar,fKStarLow,fKStarHigh);
	    fP1EarlierP2Dum = new TH1D(TString::Format("p1Ep2Dum_%s", fTitle.Data()), " ", fNbinsKStar,fKStarLow,fKStarHigh);
	    fP2EarlierP1Num = new TH1D(TString::Format("p2Ep1Num_%s", fTitle.Data()), " ", fNbinsKStar,fKStarLow,fKStarHigh);
	    fP2EarlierP1Dum = new TH1D(TString::Format("p2Ep1Dum_%s", fTitle.Data()), " ", fNbinsKStar,fKStarLow,fKStarHigh);

	    fP1EarlierP2Num->Sumw2();
	    fP1EarlierP2Dum->Sumw2();
	    fP2EarlierP1Num->Sumw2();
	    fP2EarlierP1Dum->Sumw2();
	
	    
	}
}

void AliFemtoCorrFctnpdtHe3::SetfHighCF(bool aHighCF){
    fHighCF = aHighCF;
}
void AliFemtoCorrFctnpdtHe3::SetHighCFInit(bool aHighCF){
    if(aHighCF){
	    int xybin = 19;
	    float testxy[20];
	    for(int i=0;i<20;i++){
		testxy[i] = 0.2*float(i+1);
	    }
	    float testz[401];
	    for(int i=0;i<401;i++){
		testz[i] = 0.005*float(i);
	    }
	    fNumHigh3F = new TH3F(TString::Format("fNumHigh3F_%s", fTitle.Data())," ",xybin,testxy,xybin,testxy,400,testz);
	    fDenHigh3F = new TH3F(TString::Format("fDenHigh3F_%s", fTitle.Data())," ",xybin,testxy,xybin,testxy,400,testz);
	    fNumHigh3F->Sumw2();
	    fDenHigh3F->Sumw2();
    }
}


void AliFemtoCorrFctnpdtHe3::SetfSideBand(bool aSideBand){
    fSideBand = aSideBand;
}
void AliFemtoCorrFctnpdtHe3::SetSideBandTF1Init(bool aSideBand){
    if(aSideBand){
	    //\ dowang 2022.2.15 side band
	    // signal/side band region:
	    // mTOF2 = mu + p7 + p8 * f
	    // mu = p0 + p1*(1-p2/pT)^p3
	    // f  = p4*pT + p5*pT^2 + p6*pT
	    // 
	    p2Up3Sigma	= new TF1(TString::Format("TF1_p2Up3Sigma_%s", fTitle.Data()),MassBandFunc,0,5,9);
	    p2Up3Sigma->SetParameters(1,1,1,1, 1,1,1,1,1);

	    p2Low3Sigma = new TF1(TString::Format("TF1_p2Low3Sigma_%s", fTitle.Data()),MassBandFunc,0,5,9);
	    p2Low3Sigma->SetParameters(1,1,1,1, 1,1,1,1,1);
	   

	    SideBandUp = new TF1(TString::Format("TF1_SideBandUp_%s", fTitle.Data()),MassBandFunc,0,5,9);
	    SideBandUp->SetParameters(0,0,0,6, 1,1,1,1,1);
	

	    SideBandLow = new TF1(TString::Format("TF1_SideBandLow_%s", fTitle.Data()),MassBandFunc,0,5,9);
	    SideBandLow->SetParameters(0,0,0,2, 1,1,1,1,1);
	 
    }

}
void AliFemtoCorrFctnpdtHe3::SetSideBandHistInit(bool aSideBand){
    if(aSideBand){
	    A1SideBandNum	= new TH1D(TString::Format("A1SideBandNum%s", fTitle.Data()), " ",fNbinsKStar,fKStarLow,fKStarHigh);A1SideBandNum->Sumw2();
	    S1SideBandNum	= new TH1D(TString::Format("S1SideBandNum%s", fTitle.Data()), " ",fNbinsKStar,fKStarLow,fKStarHigh);S1SideBandNum->Sumw2();
	    A2SideBandNum	= new TH1D(TString::Format("A2SideBandNum%s", fTitle.Data()), " ",fNbinsKStar,fKStarLow,fKStarHigh);A2SideBandNum->Sumw2();
	    A1A2SideBandNum	= new TH1D(TString::Format("A1A2SideBandNum%s", fTitle.Data()), " ",fNbinsKStar,fKStarLow,fKStarHigh);A1A2SideBandNum->Sumw2();
	    SignalAndSideCFNum  = new TH1D(TString::Format("SignalAndSideCFNum%s", fTitle.Data()), " ",fNbinsKStar,fKStarLow,fKStarHigh);SignalAndSideCFNum->Sumw2();

	    A1SideBandDum	= new TH1D(TString::Format("A1SideBandDum%s", fTitle.Data()), " ",fNbinsKStar,fKStarLow,fKStarHigh);A1SideBandDum->Sumw2();
	    S1SideBandDum	= new TH1D(TString::Format("S1SideBandDum%s", fTitle.Data()), " ",fNbinsKStar,fKStarLow,fKStarHigh);S1SideBandDum->Sumw2();
	    A2SideBandDum	= new TH1D(TString::Format("A2SideBandDum%s", fTitle.Data()), " ",fNbinsKStar,fKStarLow,fKStarHigh);A2SideBandDum->Sumw2();
	    A1A2SideBandDum	= new TH1D(TString::Format("A1A2SideBandDum%s", fTitle.Data()), " ",fNbinsKStar,fKStarLow,fKStarHigh);A1A2SideBandDum->Sumw2();
	    SignalAndSideCFDum  = new TH1D(TString::Format("SignalAndSideCFDum%s", fTitle.Data()), " ",fNbinsKStar,fKStarLow,fKStarHigh);SignalAndSideCFDum->Sumw2();
    }
}
void AliFemtoCorrFctnpdtHe3::SetDPhiDEtaQAInit(bool aDPhiDEtaQA){
	if(!aDPhiDEtaQA) return;
//fNumDPhiDEtaQA= new TH2F(TString::Format("fNumDPhiDEtaQA%s", fTitle.Data()), " ",200,-2.0,2.0, 100, -TMath::Pi(), TMath::Pi());
//fDumDPhiDEtaQA= new TH2F(TString::Format("fDumDPhiDEtaQA%s", fTitle.Data()), " ",200,-2.0,2.0, 100, -TMath::Pi(), TMath::Pi());

//fNumDPhiDEtaQAFailCut = new TH2F(TString::Format("fNumDPhiDEtaAvgQAFailCut%s", fTitle.Data()), " ",150, -0.15, 0.15, 150, -0.15, 0.15);
//fDumDPhiDEtaQAFailCut = new TH2F(TString::Format("fDumDPhiDEtaAvgQAFailCut%s", fTitle.Data()), " ",150, -0.15, 0.15, 150, -0.15, 0.15);
fNumDPhiDEtaAvgQA = new TH2F(TString::Format("fNumDPhiDEtaAvgQA%s", fTitle.Data()), " ",150, -0.15, 0.15, 150, -0.15, 0.15);
fDumDPhiDEtaAvgQA = new TH2F(TString::Format("fDumDPhiDEtaAvgQA%s", fTitle.Data()), " ",150, -0.15, 0.15, 150, -0.15, 0.15);
		
}
void AliFemtoCorrFctnpdtHe3::SetTF1ParaUp3Sigma(float *para){
    p2Up3Sigma->SetParameters(para[0],para[1],para[2],para[3],para[4],para[5],para[6],para[7],para[8]);
}
void AliFemtoCorrFctnpdtHe3::SetTF1ParaLow3Sigma(float *para){
    p2Low3Sigma->SetParameters(para[0],para[1],para[2],para[3],para[4],para[5],para[6],para[7],para[8]);
}
void AliFemtoCorrFctnpdtHe3::SetTF1PareSideBandUp(float *para){
    SideBandUp->SetParameters(para[0],para[1],para[2],para[3],para[4],para[5],para[6],para[7],para[8]);
}
void AliFemtoCorrFctnpdtHe3::SetTF1PareSideBandLow(float *para){
    SideBandLow->SetParameters(para[0],para[1],para[2],para[3],para[4],para[5],para[6],para[7],para[8]);
}
void AliFemtoCorrFctnpdtHe3::FillSideBandNum(AliFemtoPair* aPair){

    float c = 1.;
    float beta = aPair->Track2()->Track()->VTOF();
    if(beta==0) return;
	
    float tMom = aPair->Track2()->Track()->P().Mag();
    float massTOF = tMom*tMom/c/c*(1./(beta*beta)-1);

    float EvalMassSideBandUp 	= SideBandUp->Eval(tMom);
    float EvalMassSideBandLow 	= SideBandLow->Eval(tMom);
    // return if out of range
    if(massTOF > EvalMassSideBandUp) return;
    if(massTOF < EvalMassSideBandLow) return;
   
    double dphi = 0.;
    double deta = 0.;

	
    float EvalMassp2Up3Sigma 	= p2Up3Sigma->Eval(tMom);
    float EvalMassp2Low3Sigma 	= p2Low3Sigma->Eval(tMom);

    float tKStar = fabs(aPair->KStar());
    // total signal
    SignalAndSideCFNum->Fill(tKStar);
	//\ 3 sigma
    if(EvalMassp2Up3Sigma <= massTOF){ 
	A1SideBandNum->Fill(tKStar);
	A1A2SideBandNum->Fill(tKStar);
	
    }
    if(EvalMassp2Low3Sigma <= massTOF && massTOF < EvalMassp2Up3Sigma){
	S1SideBandNum->Fill(tKStar);
	}
    if(massTOF < EvalMassp2Low3Sigma){
	A2SideBandNum->Fill(tKStar);
	A1A2SideBandNum->Fill(tKStar);
    }
    



}
void AliFemtoCorrFctnpdtHe3::FillSideBandDum(AliFemtoPair* aPair){

    float c = 1.;
    float beta = aPair->Track2()->Track()->VTOF();
    if(beta==0.) return;

    float tMom = aPair->Track2()->Track()->P().Mag();
    float massTOF = tMom*tMom/c/c*(1./(beta*beta)-1);

    float EvalMassSideBandUp 	= SideBandUp->Eval(tMom);
    float EvalMassSideBandLow 	= SideBandLow->Eval(tMom);
    // return if out of range
    if(massTOF > EvalMassSideBandUp) return;
    if(massTOF < EvalMassSideBandLow) return;

    float EvalMassp2Up3Sigma 	= p2Up3Sigma->Eval(tMom);
    float EvalMassp2Low3Sigma 	= p2Low3Sigma->Eval(tMom);
	

    float tKStar = fabs(aPair->KStar());
    // total signal
    SignalAndSideCFDum->Fill(tKStar);
    
    if(EvalMassp2Up3Sigma <= massTOF){
 	A1SideBandDum->Fill(tKStar);
	A1A2SideBandDum->Fill(tKStar);
    }
    if(EvalMassp2Low3Sigma <= massTOF && massTOF <EvalMassp2Up3Sigma){
 	S1SideBandDum->Fill(tKStar);
	}
    if(massTOF < EvalMassp2Low3Sigma){
 	A2SideBandDum->Fill(tKStar);
	A1A2SideBandDum->Fill(tKStar);
    }
    


}
void AliFemtoCorrFctnpdtHe3::SetfUsePt(int aUsePt){
	fUsePt = aUsePt;
}
void AliFemtoCorrFctnpdtHe3::SetfUseDPhiDEtaQA(int aUseDPhiDEtaQA){
	fUseDPhiDEtaQA = aUseDPhiDEtaQA;
}

Double_t MassBandFunc(Double_t *x, Double_t *par){
    Float_t xx          = x[0];
    float p0            = par[0];
    float p1            = par[1];
    float p2            = par[2];
    float p3            = par[3];
	
    float p4            = par[4];
    float p5            = par[5];
    float p6            = par[6];
    float constIn11     = par[7];
    float Beforef       = par[8];

    float mu = p0 + p1 * TMath::Power((1. - p2/xx),p3);

    Double_t f   = 0.;
    f = p4 * xx + p5 * TMath::Power(xx,2) + p6 * TMath::Power(xx,3);

    Double_t re = 0.;
    re = mu + abs(Beforef) * constIn11 + Beforef * f;
    return re;

}
void AliFemtoCorrFctnpdtHe3::SetUseStavinskyMethod(int aUse){
	fUseStavinskyMethod = aUse;
}
void AliFemtoCorrFctnpdtHe3::SetStaSkyBkgInit(bool aInit){
	fStaSkyBkg = new TH1F(TString::Format("fStaSkyBkg%s", fTitle.Data()), "fStaSkyBkg", fNbinsKStar,fKStarLow,fKStarHigh);
	fStaSkyBkg->Sumw2();
}

AliFemtoPair * AliFemtoCorrFctnpdtHe3::InversePair(AliFemtoPair* aPair)
{
    AliFemtoPair* fPair = new AliFemtoPair;

    AliFemtoParticle *tPart1 = new AliFemtoParticle(*aPair->Track1());
    AliFemtoLorentzVector tFourMom1 = AliFemtoLorentzVector(tPart1->FourMomentum());
    tFourMom1.SetPx(-1.*tFourMom1.px());
    tFourMom1.SetPy(-1.*tFourMom1.py());
    tFourMom1.SetPz(-1.*tFourMom1.pz());
    tPart1->ResetFourMomentum(tFourMom1);
    fPair->SetTrack1(tPart1);

 AliFemtoParticle *tPart2 = new AliFemtoParticle(*aPair->Track2());
    fPair->SetTrack2(tPart2);


    return fPair;
}
void AliFemtoCorrFctnpdtHe3::SetUse2DpTvsKStar(int aUse){
	fUse2DpTvsKStar = aUse;
}
void AliFemtoCorrFctnpdtHe3::Set2DpTvsKStarInit(bool aInit,
 int nbinsks,float lowks,float upks,
int nbinspT,float lowpT,float uppT
){
	KStarVspT_P1Hist = new TH2F(TString::Format("KStarVspT_P1Hist%s", fTitle.Data())," ",nbinsks,lowks,upks, nbinspT,lowpT,uppT);	
	if(IsSameParticlePair!=1) KStarVspT_P2Hist = new TH2F(TString::Format("KStarVspT_P2Hist%s", fTitle.Data())," ",nbinsks,lowks,upks, nbinspT,lowpT,uppT); 
}

void AliFemtoCorrFctnpdtHe3::SetUsePairCutEtaPhi(int aUsePairCutEtaPhi){
	fUsePairCutEtaPhi = aUsePairCutEtaPhi;
}
void AliFemtoCorrFctnpdtHe3::SetPairCutEtaPhi(float aEtaCut,float aPhiCut){
	fPairCut_eta = aEtaCut;
	fPairCut_phi = aPhiCut;
}

bool AliFemtoCorrFctnpdtHe3::PairEtaPhiSelect(AliFemtoPair* aPair){
	
  AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  Double_t magsign = 0.0;
  if (!aodH) {
   return false;
  }
  else {
    AliAODEvent *fAOD;
    fAOD = aodH->GetEvent();
    magsign = fAOD->GetMagneticField();
  }
  float fMagSign = 1.;
  if (magsign > 1)
    fMagSign = 1.;
  else if ( magsign < 1)
    fMagSign = -1.;
  else
    fMagSign = magsign;
  
  float magval = 0.5;
  double phi1 = aPair->Track1()->Track()->P().Phi();
  double phi2 = aPair->Track2()->Track()->P().Phi();
  double chg1 = aPair->Track1()->Track()->Charge();
  double chg2 = aPair->Track2()->Track()->Charge();
  double pt1 = aPair->Track1()->Track()->Pt();
  double pt2 = aPair->Track2()->Track()->Pt();

  double eta1 = aPair->Track1()->Track()->P().PseudoRapidity();
  double eta2 = aPair->Track2()->Track()->P().PseudoRapidity();
  float deta2 = TMath::Power(eta1-eta2,2);
  float tmpCut = TMath::Power(fPairCut_eta,2) + TMath::Power(fPairCut_phi,2);

  for(int i=0;i<9;i++){
     Double_t rad = TPCradii[i];
    // Calculate dPhiStar:
    //double afsi0b = -0.07510020733*chg1*fMagSign*rad/pt1;
    //double afsi1b = -0.07510020733*chg2*fMagSign*rad/pt2;
    double afsi0b = -0.15*magval*chg1*fMagSign*rad/pt1;
    double afsi1b = -0.15*magval*chg2*fMagSign*rad/pt2;
    Double_t dphistar =  phi2 - phi1 + TMath::ASin(afsi1b) - TMath::ASin(afsi0b);
    dphistar = TVector2::Phi_mpi_pi(dphistar); // returns phi angle in the interval [-PI,PI)
   
    float tmp = deta2 + TMath::Power(dphistar,2);
    if(tmp > tmpCut) return false;
   }
   return true;
  
}
void AliFemtoCorrFctnpdtHe3::SetPassAllPair(int aUse){
	fPassAllPair = aUse;
}
void AliFemtoCorrFctnpdtHe3::SetGobalVelGate(int aUse){
	fUseGobalVelGate = aUse;
}
void AliFemtoCorrFctnpdtHe3::SetIsSameParticlePair(int aUse){
	IsSameParticlePair = aUse;
}
void AliFemtoCorrFctnpdtHe3::SetUse3DkTvsKStarvsmT(int aUse){
	fUse3DkTvsKStarvsmT = aUse;
}
float AliFemtoCorrFctnpdtHe3::CalcMt(const AliFemtoPair* aPair)
{
  return 0.5*aPair->FourMomentumSum().mt();
}
void AliFemtoCorrFctnpdtHe3::Set3DkTvsKStarvsmTInit(bool aInit,
int nbinskT,float lowkT,float upkT,
int nbinsks,float lowks,float upks,
int nbinsmT,float lowmT,float upmT

){

fNum3DkTvsKStarvsmT = new TH3F(TString::Format("fNum3DkTvsKStarvsmT%s", fTitle.Data())," ",nbinskT,lowkT,upkT,nbinsks,lowks,upks,nbinsmT,lowmT,upmT);

fDum3DkTvsKStarvsmT = new TH3F(TString::Format("fDum3DkTvsKStarvsmT%s", fTitle.Data())," ",nbinskT,lowkT,upkT,nbinsks,lowks,upks,nbinsmT,lowmT,upmT);



}
void AliFemtoCorrFctnpdtHe3::Set2DKstarVsmT(int aUse){

fUse2DkStarVsmT = aUse;
}
void AliFemtoCorrFctnpdtHe3::Set2DkStarVsmTInit(bool aInit,
                int nbinsks,float lowks,float upks,
                int nbinsmT,float lowmT,float upmT){


fNum2DkStarVsmT = new TH2F(TString::Format("fNum2DkStarVsmT%s", fTitle.Data())," ",nbinsks,lowks,upks,nbinsmT,lowmT,upmT);
fDum2DkStarVsmT = new TH2F(TString::Format("fDum2DkStarVsmT%s", fTitle.Data())," ",nbinsks,lowks,upks,nbinsmT,lowmT,upmT);

}

void AliFemtoCorrFctnpdtHe3::SetdBumpCheck(int aUse){
	fUseBumpC = aUse;
}

void AliFemtoCorrFctnpdtHe3::SetdBumpCheckInit(bool aInit,
 int nbinsks,float lowks,float upks,
int nbinspT,float lowpT,float uppT,
 int nbinsMass,float lowMass,float upMass

){

f2DkSVspT = new TH2F(TString::Format("f2DkSVspT%s", fTitle.Data())," ",nbinsks,lowks,upks,nbinspT,lowpT,uppT);
f2DkSVsMass = new TH2F(TString::Format("f2DkSVsMass%s", fTitle.Data())," ",nbinsks,lowks,upks,nbinsMass,lowMass,upMass);



}
void AliFemtoCorrFctnpdtHe3::SetmTLimitCheck(int aUse){
	fUsemTCheck = aUse;
}
void AliFemtoCorrFctnpdtHe3::SetmTCheckInit(bool aInit,
 int nbinsks,float lowks,float upks,
int nbinspT,float lowpT,float uppT,
int nbinsmT,float lowmT,float upmT
){

// in mix

//fUsemTCheck = 1: mT vs pT
//fUsemTCheck = 2 : kT vs pT 
f3DmTDepkSVspT = new TH3F(TString::Format("f3DmTDepkSVspT%s", fTitle.Data())," ",nbinsks,lowks,upks,nbinspT,lowpT,uppT,nbinsmT,lowmT,upmT);

}

