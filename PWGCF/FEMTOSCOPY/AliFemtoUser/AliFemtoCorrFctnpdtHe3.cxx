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
    fNumHigh3FDPhiDEtaQA(nullptr),
    fDumHigh3FDPhiDEtaQA(nullptr),
    fNumHigh3FDPhiDEtaQAvsKStar(nullptr),
    fDumHigh3FDPhiDEtaQAvsKStar(nullptr),
    fDPhiDEtaQA_A1SideBandNum(nullptr),
    fDPhiDEtaQA_A1SideBandDum(nullptr),
    fDPhiDEtaQA_S1SideBandNum(nullptr),
    fDPhiDEtaQA_S1SideBandDum(nullptr),
    fDPhiDEtaQA_A2SideBandNum(nullptr),
    fDPhiDEtaQA_A2SideBandDum(nullptr)
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
    fNumHigh3FDPhiDEtaQA(fNumHigh3FDPhiDEtaQA),
    fDumHigh3FDPhiDEtaQA(fDumHigh3FDPhiDEtaQA),
    fNumHigh3FDPhiDEtaQAvsKStar(fNumHigh3FDPhiDEtaQAvsKStar),
    fDumHigh3FDPhiDEtaQAvsKStar(fDumHigh3FDPhiDEtaQAvsKStar),
    fDPhiDEtaQA_A1SideBandNum(fDPhiDEtaQA_A1SideBandNum),
    fDPhiDEtaQA_A1SideBandDum(fDPhiDEtaQA_A1SideBandDum),
    fDPhiDEtaQA_S1SideBandNum(fDPhiDEtaQA_S1SideBandNum),
    fDPhiDEtaQA_S1SideBandDum(fDPhiDEtaQA_S1SideBandDum),
    fDPhiDEtaQA_A2SideBandNum(fDPhiDEtaQA_A2SideBandNum),
    fDPhiDEtaQA_A2SideBandDum(fDPhiDEtaQA_A2SideBandDum)
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
    
    delete fNumHigh3FDPhiDEtaQA;
    delete fDumHigh3FDPhiDEtaQA;    
    delete fNumHigh3FDPhiDEtaQAvsKStar;
    delete fDumHigh3FDPhiDEtaQAvsKStar;

    delete fDPhiDEtaQA_A1SideBandNum;
    delete fDPhiDEtaQA_A1SideBandDum;
    delete fDPhiDEtaQA_S1SideBandNum;
    delete fDPhiDEtaQA_S1SideBandDum;
    delete fDPhiDEtaQA_A2SideBandNum;
    delete fDPhiDEtaQA_A2SideBandDum;
 

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
    if(fNumHigh3FDPhiDEtaQA) delete fNumHigh3FDPhiDEtaQA;
        		fNumHigh3FDPhiDEtaQA = new TH3F(*aCorrFctn.fNumHigh3FDPhiDEtaQA);
    if(fDumHigh3FDPhiDEtaQA) delete fDumHigh3FDPhiDEtaQA;
        		fDumHigh3FDPhiDEtaQA = new TH3F(*aCorrFctn.fDumHigh3FDPhiDEtaQA);
	if(fNumHigh3FDPhiDEtaQAvsKStar) delete fNumHigh3FDPhiDEtaQAvsKStar;
        		fNumHigh3FDPhiDEtaQAvsKStar = new TH3F(*aCorrFctn.fNumHigh3FDPhiDEtaQAvsKStar);
    	if(fDumHigh3FDPhiDEtaQAvsKStar) delete fDumHigh3FDPhiDEtaQAvsKStar;
        		fDumHigh3FDPhiDEtaQAvsKStar = new TH3F(*aCorrFctn.fDumHigh3FDPhiDEtaQAvsKStar);

    if(fDPhiDEtaQA_A1SideBandNum) delete fDPhiDEtaQA_A1SideBandNum;
        		fDPhiDEtaQA_A1SideBandNum = new TH2F(*aCorrFctn.fDPhiDEtaQA_A1SideBandNum);
    if(fDPhiDEtaQA_A1SideBandDum) delete fDPhiDEtaQA_A1SideBandDum;
        		fDPhiDEtaQA_A1SideBandDum = new TH2F(*aCorrFctn.fDPhiDEtaQA_A1SideBandDum);

    if(fDPhiDEtaQA_S1SideBandNum) delete fDPhiDEtaQA_S1SideBandNum;
        		fDPhiDEtaQA_S1SideBandNum = new TH2F(*aCorrFctn.fDPhiDEtaQA_S1SideBandNum);
    if(fDPhiDEtaQA_S1SideBandDum) delete fDPhiDEtaQA_S1SideBandDum;
        		fDPhiDEtaQA_S1SideBandDum = new TH2F(*aCorrFctn.fDPhiDEtaQA_S1SideBandDum);

    if(fDPhiDEtaQA_A2SideBandNum) delete fDPhiDEtaQA_A2SideBandNum;
        		fDPhiDEtaQA_A2SideBandNum = new TH2F(*aCorrFctn.fDPhiDEtaQA_A2SideBandNum);
    if(fDPhiDEtaQA_A2SideBandDum) delete fDPhiDEtaQA_A2SideBandDum;
        		fDPhiDEtaQA_A2SideBandDum = new TH2F(*aCorrFctn.fDPhiDEtaQA_A2SideBandDum);

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
            tOutputList->Add(fNumHigh3FDPhiDEtaQA);
	    tOutputList->Add(fDumHigh3FDPhiDEtaQA);
	    tOutputList->Add(fNumHigh3FDPhiDEtaQAvsKStar);
	    tOutputList->Add(fDumHigh3FDPhiDEtaQAvsKStar);


		if(fSideBand){
		    tOutputList->Add(fDPhiDEtaQA_A1SideBandNum);
		    tOutputList->Add(fDPhiDEtaQA_A1SideBandDum);
		    tOutputList->Add(fDPhiDEtaQA_S1SideBandNum);
		    tOutputList->Add(fDPhiDEtaQA_S1SideBandDum);
		    tOutputList->Add(fDPhiDEtaQA_A2SideBandNum);
		    tOutputList->Add(fDPhiDEtaQA_A2SideBandDum);
		   
		}
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
		  fNumHigh3FDPhiDEtaQA->Write();
		  fDumHigh3FDPhiDEtaQA->Write();
		  fNumHigh3FDPhiDEtaQAvsKStar->Write();
		  fDumHigh3FDPhiDEtaQAvsKStar->Write();
		if(fSideBand){
		     fDPhiDEtaQA_A1SideBandNum->Write();
		     fDPhiDEtaQA_A1SideBandDum->Write();
		     fDPhiDEtaQA_S1SideBandNum->Write();
		     fDPhiDEtaQA_S1SideBandDum->Write();
		     fDPhiDEtaQA_A2SideBandNum->Write();
		     fDPhiDEtaQA_A2SideBandDum->Write();
		    
		}
	}
	

}
void AliFemtoCorrFctnpdtHe3::AddRealPair(AliFemtoPair* aPair)
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
    if (fPairCut && !fPairCut->Pass(fPair)) {
        return;
    }

    double tKStar = fabs(fPair->KStar());
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
    if(fUseDPhiDEtaQA){
	  double phi1 = fPair->Track1()->FourMomentum().Phi();
	  double phi2 = fPair->Track2()->FourMomentum().Phi();
	  double eta1 = fPair->Track1()->FourMomentum().PseudoRapidity();
	  double eta2 = fPair->Track2()->FourMomentum().PseudoRapidity();

	  double dphi = phi1 - phi2;
	  while (dphi<fphiL) dphi+=PIT;
	  while (dphi>fphiT) dphi-=PIT;

	  double deta = eta1 - eta2;

 	  fNumHigh3FDPhiDEtaQA->Fill(dphi,deta,fPair->Track2()->Track()->Pt());
	  fNumHigh3FDPhiDEtaQAvsKStar->Fill(dphi,deta,tKStar);
	
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
    if (fPairCut && !fPairCut->Pass(fPair)) {
        return;
    } 

    double tKStar = fabs(fPair->KStar());
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
	if(fUseDPhiDEtaQA){

	  double phi1 = fPair->Track1()->FourMomentum().Phi();
	  double phi2 = fPair->Track2()->FourMomentum().Phi();
	  double eta1 = fPair->Track1()->FourMomentum().PseudoRapidity();
	  double eta2 = fPair->Track2()->FourMomentum().PseudoRapidity();

	  double dphi = phi1 - phi2;
	  while (dphi<fphiL) dphi+=PIT;
	  while (dphi>fphiT) dphi-=PIT;

	  double deta = eta1 - eta2;

	  fDumHigh3FDPhiDEtaQA->Fill(dphi,deta,fPair->Track2()->Track()->Pt());
	  fDumHigh3FDPhiDEtaQAvsKStar->Fill(dphi,deta,tKStar);
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
    float P1velocity = TotalP/TotalE;

    //
    
    AliFemtoParticle *tPart2 = new AliFemtoParticle(*aPair->Track2());
    AliFemtoLorentzVector tFourMom2 = AliFemtoLorentzVector(tPart2->FourMomentum());
    TotalP = tFourMom2.px() * tFourMom2.px() + tFourMom2.py() * tFourMom2.py() + tFourMom2.pz() * tFourMom2.pz();
    TotalE = TMath::Sqrt( fP2Mass * fP2Mass + TotalP);
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
fNumHigh3FDPhiDEtaQA= new TH3F(TString::Format("fNumHigh3FDPhiDEtaQA%s", fTitle.Data()), " ",20,fphiL,fphiT,20,-2.0,2.0,100,0,5);
fDumHigh3FDPhiDEtaQA= new TH3F(TString::Format("fDumHigh3FDPhiDEtaQA%s", fTitle.Data()), " ",20,fphiL,fphiT,20,-2.0,2.0,100,0,5);
fNumHigh3FDPhiDEtaQAvsKStar = new TH3F(TString::Format("fNumHigh3FDPhiDEtaQAvsKStar%s", fTitle.Data()), " ",20,fphiL,fphiT,20,-2.0,2.0,200,0,1);
fDumHigh3FDPhiDEtaQAvsKStar = new TH3F(TString::Format("fDumHigh3FDPhiDEtaQAvsKStar%s", fTitle.Data()), " ",20,fphiL,fphiT,20,-2.0,2.0,200,0,1);
		if(fSideBand){
fDPhiDEtaQA_A1SideBandNum = new TH2F(TString::Format("fDPhiDEtaQA_A1SideBandNum%s", fTitle.Data()), " ",20,fphiL,fphiT,20,-2.0,2.0);
fDPhiDEtaQA_A1SideBandDum = new TH2F(TString::Format("fDPhiDEtaQA_A1SideBandDum%s", fTitle.Data()), " ",20,fphiL,fphiT,20,-2.0,2.0);
fDPhiDEtaQA_S1SideBandNum = new TH2F(TString::Format("fDPhiDEtaQA_S1SideBandNum%s", fTitle.Data()), " ",20,fphiL,fphiT,20,-2.0,2.0);
fDPhiDEtaQA_S1SideBandDum = new TH2F(TString::Format("fDPhiDEtaQA_S1SideBandDum%s", fTitle.Data()), " ",20,fphiL,fphiT,20,-2.0,2.0);
fDPhiDEtaQA_A2SideBandNum = new TH2F(TString::Format("fDPhiDEtaQA_A2SideBandNum%s", fTitle.Data()), " ",20,fphiL,fphiT,20,-2.0,2.0);
fDPhiDEtaQA_A2SideBandDum = new TH2F(TString::Format("fDPhiDEtaQA_A2SideBandDum%s", fTitle.Data()), " ",20,fphiL,fphiT,20,-2.0,2.0);

		}
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

	if(fUseDPhiDEtaQA){

	  double phi1 = aPair->Track1()->FourMomentum().Phi();
	  double phi2 = aPair->Track2()->FourMomentum().Phi();
	  double eta1 = aPair->Track1()->FourMomentum().PseudoRapidity();
	  double eta2 = aPair->Track2()->FourMomentum().PseudoRapidity();

	  dphi = phi1 - phi2;
	  while (dphi<fphiL) dphi+=PIT;
	  while (dphi>fphiT) dphi-=PIT;

	  deta = eta1 - eta2;
		    
	}
    float EvalMassp2Up3Sigma 	= p2Up3Sigma->Eval(tMom);
    float EvalMassp2Low3Sigma 	= p2Low3Sigma->Eval(tMom);

    float tKStar = fabs(aPair->KStar());
    // total signal
    SignalAndSideCFNum->Fill(tKStar);
	//\ 3 sigma
    if(EvalMassp2Up3Sigma <= massTOF){ 
	A1SideBandNum->Fill(tKStar);
	A1A2SideBandNum->Fill(tKStar);
	if(fUseDPhiDEtaQA){
		fDPhiDEtaQA_A1SideBandNum->Fill(dphi,deta);
	}
    }
    if(EvalMassp2Low3Sigma <= massTOF && massTOF < EvalMassp2Up3Sigma){
	S1SideBandNum->Fill(tKStar);
	if(fUseDPhiDEtaQA) fDPhiDEtaQA_S1SideBandNum->Fill(dphi,deta);
	}
    if(massTOF < EvalMassp2Low3Sigma){
	A2SideBandNum->Fill(tKStar);
	A1A2SideBandNum->Fill(tKStar);
	if(fUseDPhiDEtaQA){
		fDPhiDEtaQA_A2SideBandNum->Fill(dphi,deta);
		
		}
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

	double dphi = 0.;
    	double deta = 0.;

	if(fUseDPhiDEtaQA){

	  double phi1 = aPair->Track1()->FourMomentum().Phi();
	  double phi2 = aPair->Track2()->FourMomentum().Phi();
	  double eta1 = aPair->Track1()->FourMomentum().PseudoRapidity();
	  double eta2 = aPair->Track2()->FourMomentum().PseudoRapidity();

	  dphi = phi1 - phi2;
	  while (dphi<fphiL) dphi+=PIT;
	  while (dphi>fphiT) dphi-=PIT;

	  deta = eta1 - eta2;
		    
	}

    float EvalMassp2Up3Sigma 	= p2Up3Sigma->Eval(tMom);
    float EvalMassp2Low3Sigma 	= p2Low3Sigma->Eval(tMom);
	

    float tKStar = fabs(aPair->KStar());
    // total signal
    SignalAndSideCFDum->Fill(tKStar);
    
    if(EvalMassp2Up3Sigma <= massTOF){
 	A1SideBandDum->Fill(tKStar);
	A1A2SideBandDum->Fill(tKStar);
	if(fUseDPhiDEtaQA){
		fDPhiDEtaQA_A1SideBandDum->Fill(dphi,deta);
	}
    }
    if(EvalMassp2Low3Sigma <= massTOF && massTOF <EvalMassp2Up3Sigma){
 	S1SideBandDum->Fill(tKStar);
	if(fUseDPhiDEtaQA) fDPhiDEtaQA_S1SideBandDum->Fill(dphi,deta);
	}
    if(massTOF < EvalMassp2Low3Sigma){
 	A2SideBandDum->Fill(tKStar);
	A1A2SideBandDum->Fill(tKStar);
	if(fUseDPhiDEtaQA){
		fDPhiDEtaQA_A2SideBandDum->Fill(dphi,deta);
		}
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

