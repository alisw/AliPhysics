#ifndef ALIFEMTOCORRFCTNPDTHE3_H
#define ALIFEMTOCORRFCTNPDTHE3_H

#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH3F.h"
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
        virtual AliFemtoCorrFctnpdtHe3* Clone() const  { return new AliFemtoCorrFctnpdtHe3(*this); }
    protected:
        int isHe3Pair;
        TString fTitle;
        int fNbinsKStar;
        double fKStarLow, fKStarHigh;
        TH1D* fNumerator;          // numerator - real pairs
        TH1D* fDenominator;        // denominator - mixed pairs

        std::vector<TH1D*> fQANumerator;    // QA numerator [0]: for p1 pt > 0.2GeV/c; [1]: for p1 pt > 0.3GeV/c; [2]: for p1 pt > 0.4GeV/c
        std::vector<TH1D*> fQADenominator;  // QA denominator
      
};
#endif


