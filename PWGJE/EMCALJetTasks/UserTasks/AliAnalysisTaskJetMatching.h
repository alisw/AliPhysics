#ifndef AliAnalysisTaskJetMatching_H
#define AliAnalysisTaskJetMatching_H

#include <AliAnalysisTaskEmcalJetDev.h>
#include <AliEmcalJet.h>
#include <AliVTrack.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TRandom3.h>

class TF1;
class THF1;
class THF2;
class TProfile;
class AliRhoParameter;
class AliLocalRhoParameter;
class TClonesArray;

class AliAnalysisTaskJetMatching : public AliAnalysisTaskEmcalJetDev
{
    public:
        // enumerators
        enum matchingSceme      {kGeoEtaPhi, kGeoR, kGeoEtaPhiArea, kGeoRArea, kDeepMatching};
        enum duplicateRecovery  {kDoNothing, kTraceDuplicates, kRemoveDuplicates}; 
        enum sourceBKG          {kNoSourceBKG, kSourceRho, kSourceLocalRho};
        enum targetBKG          {kNoTargetBKG, kTargetRho, kTargetLocalRho};
        // constructors, destructor
                                AliAnalysisTaskJetMatching();
                                AliAnalysisTaskJetMatching(const char *name);
        virtual                 ~AliAnalysisTaskJetMatching();
        // setting up the task and technical aspects
        void                    ExecOnce();
        virtual void            UserCreateOutputObjects();
        TH1F*                   BookTH1F(const char* name, const char* x, Int_t bins, Double_t min, Double_t max, Bool_t append = kTRUE);
        TH2F*                   BookTH2F(const char* name, const char* x, const char* y, Int_t binsx, Double_t minx, Double_t maxx, Int_t binsy, Double_t miny, Double_t maxy, Bool_t append = kTRUE);
        virtual Bool_t          Run();
        /* inline */    Double_t PhaseShift(Double_t x) const {  
            while (x>=TMath::TwoPi())x-=TMath::TwoPi();
            while (x<0.)x+=TMath::TwoPi();
            return x; }
        /* inline */    Double_t PhaseShift(Double_t x, Double_t n) const {
            x = PhaseShift(x);
            if(TMath::Nint(n)==2) while (x>TMath::Pi()) x-=TMath::Pi();
            if(TMath::Nint(n)==3) {
                if(x>2.*TMath::TwoPi()/n) x = TMath::TwoPi() - x;
                if(x>TMath::TwoPi()/n) x = TMath::TwoPi()-(x+TMath::TwoPi()/n);
            }
            return x; }
        /* inline */    Bool_t CompareTracks(AliVParticle* a, AliVParticle* b) const {
            return (TMath::AreEqualAbs(a->Pt(), b->Pt(), 1e-2) && TMath::AreEqualAbs(a->Eta(), b->Eta(), 1e-2) && TMath::AreEqualAbs(a->Phi(), b->Phi(), 1e-2)) ? kTRUE : kFALSE; }
        /* inline */    Double_t GetR(AliVParticle* a, AliVParticle* b) const {
               if(!(a&&b)) return 999;
               Double_t phiA(a->Phi()), phiB(b->Phi()), etaA(a->Eta()), etaB(b->Eta());
               if(TMath::Abs(phiA-phiB) > TMath::Abs(phiA-phiB + TMath::TwoPi())) phiA+=TMath::TwoPi();
               if(TMath::Abs(phiA-phiB) > TMath::Abs(phiA-phiB - TMath::TwoPi())) phiA-=TMath::TwoPi();
               return TMath::Sqrt(TMath::Abs((etaA-etaB)*(etaA-etaB)+(phiA-phiB)*(phiA-phiB))); }

        // setters - setup how to run
        void                    SetDebugMode(Int_t d)                           {fDebug = d;}
        void                    SetMatchingScheme(matchingSceme m)              {fMatchingScheme = m;}
        void                    SetMatchConstituents(Bool_t m)                  {fMatchConstituents = m;}
        void                    SetMinFracRecoveredConstituents(Float_t f)      {fMinFracRecoveredConstituents = f;}
        void                    SetDuplicateRecoveryScheme(duplicateRecovery d) {fDuplicateJetRecoveryMode = d;}
        void                    SetUseEmcalBaseJetCuts(Bool_t b)                {fUseEmcalBaseJetCuts = b;}
        void                    SetSourceBKG(sourceBKG b)                       {fSourceBKG = b;}
        void                    SetTargetBKG(targetBKG b)                       {fTargetBKG = b;}
        void                    SetSourceJetsName(const char* n)                {fSourceJetsName = n;}
        void                    SetTargetJetsName(const char* n)                {fTargetJetsName = n; }
        void                    SetMatchedJetsName(const char* n)               {fMatchedJetsName = n;}
        void                    SetSourceLocalRhoName(const char* n)            {fSourceRhoName = n;}
        void                    SetTargetLocalRhoName(const char* n)            {fTargetRhoName = n;}
        void                    SetSourceRadius(Float_t r)                      {fSourceRadius = r;}
        void                    SetTargetRadius(Float_t r)                      {fTargetRadius = r;}
        void                    SetMatchEta(Float_t f)                          {fMatchEta = f;}
        void                    SetMatchPhi(Float_t f)                          {fMatchPhi = f;}
        void                    SetMatchR(Float_t f)                            {fMatchR = f;}
        void                    SetMatchRelArea(Float_t f)                      {fMatchArea = f;}
        void                    SetMaxRelEnergyDiff(Float_t f)                  {fMaxRelEnergyDiff = f;}
        void                    SetMaxAbsEnergyDiff(Float_t f)                  {fMaxAbsEnergyDiff = f;}
        // methods
        void                    DoGeometricMatchingEtaPhi(Bool_t pairCuts = kFALSE);
        void                    DoGeometricMatchingR(Bool_t pairCuts = kFALSE);
        void                    DoDeepMatching();
        void                    DuplicateJetRecovery();
        void                    PostMatchedJets();
        // fill histogrmas
        void                    FillAnalysisSummaryHistogram() const;
        void                    FillMatchedJetHistograms() const;
        // setters - analysis details
        /* inline */    Bool_t PassesCuts(const AliVTrack* track) const {
            return (!track || TMath::Abs(track->Eta()) > 0.9 || track->Phi() < 0 || track->Phi() > TMath::TwoPi()) ? kFALSE : kTRUE; }
        /* inline */    Bool_t PassesCuts(AliEmcalJet* jet) const { return (jet) ? kTRUE : kFALSE; }
        /* inline */    Bool_t PassesCuts(AliEmcalJet* a, AliEmcalJet* b) const {
            if (fMatchArea > 0) { return (TMath::Abs(a->Area()/b->Area()) < fMatchArea) ? kTRUE : kFALSE; }
            if (fMaxRelEnergyDiff > 0) { return (TMath::Abs(a->E()/b->E()) > fMaxRelEnergyDiff) ? kTRUE : kFALSE; }
            if (fMaxAbsEnergyDiff > 0) { return (TMath::Abs(a->E()-b->E()) > fMaxAbsEnergyDiff) ? kTRUE : kFALSE; }
            return kTRUE; }
        /* inline */    void ClearMatchedJetsCache() {
            for (Int_t i(0); i < fNoMatchedJets; i++) {
                fMatchedJetContainer[i][0] = 0x0; fMatchedJetContainer[i][1] = 0x0; }
            fNoMatchedJets = 0;
        }
        void                    PrintInfo() const;
        virtual void            Terminate(Option_t* option);

    private: 
        Int_t                   fDebug;                 // debug level (0 none, 1 fcn calls, 2 verbose)
        TClonesArray*           fSourceJets;            //! array with source jets
        TString                 fSourceJetsName;        // name of array with source jets
        TClonesArray*           fTargetJets;            //! array with target jets
        TString                 fTargetJetsName;        // name of array with target jets
        TClonesArray*           fMatchedJets;           //! final list of matched jets which is added to event
        TString                 fMatchedJetsName;       // name of list of matched jets
        AliLocalRhoParameter*   fSourceRho;             //! source rho
        TString                 fSourceRhoName;         // source rho  name
        AliLocalRhoParameter*   fTargetRho;             //! target rho
        TString                 fTargetRhoName;         // target rho name
        Bool_t                  fUseScaledRho;          // use scaled rho
        Float_t                 fSourceRadius;          // source radius 
        Float_t                 fTargetRadius;          // target radius
        matchingSceme           fMatchingScheme;        // select your favorite matching algorithm
        duplicateRecovery       fDuplicateJetRecoveryMode;      // what to do with duplicate matches
        Bool_t                  fUseEmcalBaseJetCuts;   // use the emcal jet base class for jet cuts
        sourceBKG               fSourceBKG;             // subtracted background of source jets
        targetBKG               fTargetBKG;             // subtracted background of target jets
        // additional jet cuts (most are inherited)
        Float_t                 fLocalJetMinEta;        // local eta cut for jets
        Float_t                 fLocalJetMaxEta;        // local eta cut for jets
        Float_t                 fLocalJetMinPhi;        // local phi cut for jets
        Float_t                 fLocalJetMaxPhi;        // local phi cut for jets
        // transient object pointers
        TList*                  fOutputList;            //! output list
        TH1F*                   fHistUnsortedCorrelation;       //! dphi correlation of unsorted jets
        TH1F*                   fHistMatchedCorrelation;        //! dphi correlation of matched jets
        TH1F*                   fHistSourceJetPt;       //! pt of source jets
        TH1F*                   fHistTargetJetPt;       //! pt of target jets
        TH1F*                   fHistMatchedJetPt;      //! pt of matched jets
        TH1F*                   fHistNoConstSourceJet;  //! number of constituents of source jet
        TH1F*                   fHistNoConstTargetJet;  //! number of constituents of target jet
        TH1F*                   fHistNoConstMatchJet;   //! number of constituents of matched jets
        TProfile*               fProfFracPtMatched;     //! sum pt fraction for matched tracks / jet
        TProfile*               fProfFracPtJets;        //! sum pt fraction for matched jets
        TProfile*               fProfFracNoMatched;     //! no of constituents fraction found / jet
        TProfile*               fProfFracNoJets;        //! no of consstituents fraction jet / jet
        TH1F*                   fHistAnalysisSummary;   //! flags
        TProfile*               fProfQAMatched;         //! QA spreads of matched jets
        TProfile*               fProfQA;                //! QA spreads of source and target jets
        Int_t                   fNoMatchedJets;         //! number of matched jets
        AliEmcalJet*            fMatchedJetContainer[100][2];   //! pointers to matched jets
        // geometric matching parameters
        Float_t                 fMatchEta;              // max eta distance between centers of matched jets
        Float_t                 fMatchPhi;              // max phi distance between centers of matched jets
        Float_t                 fMatchR;                // max distance between matched jets
        Float_t                 fMatchArea;             // max relative area mismatch between matched jets
        Float_t                 fMaxRelEnergyDiff;      // max relative energy difference between matched jets
        Float_t                 fMaxAbsEnergyDiff;      // max absolute energy difference between matched jets
        Bool_t                  fMatchConstituents;     // match constituents
        Float_t                 fMinFracRecoveredConstituents;  // minimium fraction of constituents that needs to be found

        AliAnalysisTaskJetMatching(const AliAnalysisTaskJetMatching&);                  // not implemented
        AliAnalysisTaskJetMatching& operator=(const AliAnalysisTaskJetMatching&);       // not implemented

        ClassDef(AliAnalysisTaskJetMatching, 4);
};

#endif
