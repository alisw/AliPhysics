///
/// \file AliFemtoCorrFctnpdtHe3.cxx
///

#include "AliFemtoCorrFctnpdtHe3.h"
const float He3Mass = 2.8084;
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
    fNumerator(nullptr),
    fDenominator(nullptr),
    fQANumerator(3, nullptr),
    fQADenominator(3, nullptr)
{
    
    fNumerator      = new TH1D(TString::Format("Num%s", fTitle.Data()), "fNumerator", nbins, KStarLo, KStarHi);
    fDenominator    = new TH1D(TString::Format("Dum%s", fTitle.Data()), "fDenominator", nbins, KStarLo, KStarHi);
    
    fNumerator->Sumw2();
    fDenominator->Sumw2();

    for (int ihist = 0; ihist < 3; ihist++) {

        TString suffix = TString::Format("fQANumerator_%d_%s",ihist, fTitle.Data());
        fQANumerator[ihist] = new TH1D(suffix, "fQANumerator", nbins, KStarLo, KStarHi);
        suffix = TString::Format("fQADenominator_%d_%s",ihist, fTitle.Data());
        fQADenominator[ihist] = new TH1D(suffix, "fQADenominator", nbins, KStarLo, KStarHi);
        

        fQANumerator[ihist]->Sumw2();
        fQADenominator[ihist]->Sumw2();

    }



}
AliFemtoCorrFctnpdtHe3::AliFemtoCorrFctnpdtHe3(const AliFemtoCorrFctnpdtHe3& aCorrFctn):
    AliFemtoCorrFctn(aCorrFctn),
    isHe3Pair(aCorrFctn.isHe3Pair),
    fTitle(aCorrFctn.fTitle),
    fNbinsKStar(aCorrFctn.fNbinsKStar),
    fKStarLow(aCorrFctn.fKStarLow),
    fKStarHigh(aCorrFctn.fKStarHigh),
    fNumerator(aCorrFctn.fNumerator ? new TH1D(*aCorrFctn.fNumerator) : nullptr),
    fDenominator(aCorrFctn.fDenominator ? new TH1D(*aCorrFctn.fDenominator) : nullptr),
    fQANumerator(3, nullptr),
    fQADenominator(3, nullptr)
{
    
    // copy constructor
    for (int ihist = 0; ihist < 3; ihist++) {
        fQANumerator[ihist]     = new TH1D(*aCorrFctn.fQANumerator[ihist]);
        fQADenominator[ihist]   = new TH1D(*aCorrFctn.fQADenominator[ihist]);
    }


}
AliFemtoCorrFctnpdtHe3::~AliFemtoCorrFctnpdtHe3()
{
    // destructor
    delete fNumerator;
    delete fDenominator;
    for (int ihist = 0; ihist < 3; ihist++) {
        delete fQANumerator[ihist];
        delete fQADenominator[ihist];
       
    }
    
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

    if(fNumerator) delete fNumerator;
    fNumerator = new TH1D(*aCorrFctn.fNumerator);
    if(fDenominator) delete fDenominator;
    fDenominator = new TH1D(*aCorrFctn.fDenominator);

    if (fQANumerator.size() == aCorrFctn.fQANumerator.size()) {
        // simple assignments if same number of histograms

        for (int ihist = 0; ihist < 3; ihist++) {
            *fQANumerator[ihist] = *aCorrFctn.fQANumerator[ihist];
            *fQADenominator[ihist] = *aCorrFctn.fQADenominator[ihist];
        }

    } 
    else {
        // delete and resize to match
        for (size_t i = 0; i < 3; ++i) {
            delete fQANumerator[i];
            delete fQADenominator[i];
        }

        fQANumerator.resize(3, nullptr);
        fQADenominator.resize(3, nullptr);

        for (int ihist = 0; ihist < 3; ihist++) {
            fQANumerator[ihist]     = new TH1D(*aCorrFctn.fQANumerator[ihist]);
            fQADenominator[ihist]   = new TH1D(*aCorrFctn.fQADenominator[ihist]);
        }
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

    for (int ihist = 0; ihist < 3; ihist++) {

        tOutputList->Add(fQANumerator[ihist]);
        tOutputList->Add(fQADenominator[ihist]);

    }
    return tOutputList;
}
void AliFemtoCorrFctnpdtHe3::Finish()
{
    cout<<"AliFemtoCorrFctnpdtHe3::Finish()"<<endl;
}
void AliFemtoCorrFctnpdtHe3::Write()
{
    // Write out neccessary objects
    fNumerator->Write();
    fDenominator->Write();
    for (int ihist = 0; ihist < 3; ihist++) {
        fQANumerator[ihist]->Write();
        fQADenominator[ihist]->Write();
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

    float p1pT = fPair->Track1()->Track()->Pt();
    if(p1pT > 0.2) fQANumerator[0]->Fill(tKStar);
    if(p1pT > 0.3) fQANumerator[1]->Fill(tKStar);
    if(p1pT > 0.4) fQANumerator[2]->Fill(tKStar);
  
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

    float p1pT = fPair->Track1()->Track()->Pt();
    if(p1pT > 0.2) fQADenominator[0]->Fill(tKStar);
    if(p1pT > 0.3) fQADenominator[1]->Fill(tKStar);
    if(p1pT > 0.4) fQADenominator[2]->Fill(tKStar);



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


