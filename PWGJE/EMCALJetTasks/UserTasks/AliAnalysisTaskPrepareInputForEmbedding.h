/// \class AliAnalysisTaskPrepareInputForEmbedding
/// \brief Save in a TTree the 4-momentum vector of the jets at reco and particle level
///
/// Two jet container are read, the area- and constituent-based background subtracted one
/// The correspondent particle level jet is also stored. The task requires the JetTagger for PYTHIA matching to be run first
/// The output tree is thought as an input to the single track embedding
///
/// \author Chiara Bianchin
/// \date Aug, 2015

#ifndef ALIANALYSISTASKPREPAREINPUTFOREMBEDDING_H
#define ALIANALYSISTASKPREPAREINPUTFOREMBEDDING_H

class TLorentzVector;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskPrepareInputForEmbedding : public AliAnalysisTaskEmcalJet {

public:
   AliAnalysisTaskPrepareInputForEmbedding();
   AliAnalysisTaskPrepareInputForEmbedding(const char* name);
   virtual ~AliAnalysisTaskPrepareInputForEmbedding();

   void        UserCreateOutputObjects();
   void        Terminate(Option_t *option);
   
   //setters
   void        SetContainerArea(Int_t number)          {fContainerArea = number;}
   void        SetContainerCons(Int_t number)          {fContainerConst = number;}
   void        SetMinFracShared(Double_t minfrac)      {fMinFractionShared = minfrac;}
   void        SetLeadingJetOnly(Bool_t b = kTRUE)     {fLeadingJetOnly = b;}
   
protected:
   Bool_t      Run();
   Bool_t      FillHistograms();

private:
   Int_t        fContainerArea;               ///< container number for area sub jet
   Int_t        fContainerConst;              ///< container number for const sub jet
   Double_t     fMinFractionShared;           ///< minimum pT fraction shared between reco and part level jets
   Bool_t       fLeadingJetOnly;              ///< fill TTrees with leading jet only
   TTree       *fTreeJetsA;                   ///!<! tree with the TLorentzVector of the jet detector and particle level
   TTree       *fTreeJetsC;                   ///!<! tree with the TLorentzVector of the jet detector and particle level
   TLorentzVector *fJetGenSub;                ///< reconstucted jet, area based method
   TLorentzVector *fJetConSub;                ///< reconstucted jet, constituent subtraction method
   TLorentzVector *fJetPart1;                 ///< particle level jet
   TLorentzVector *fJetPart2;                 ///< particle level jet
   TLorentzVector *fJetGenSubL;               ///< reconstucted leading jet, area based method
   TLorentzVector *fJetConSubL;               ///< reconstucted leading jet, constituent subtraction method
   TLorentzVector *fJetPart1L;                ///< particle level leading jet
   TLorentzVector *fJetPart2L;                ///< particle level leading jet
   TH1F         *fNumberOfJets;               ///!<! histograms with number of jets selected
   TH2F         *fhFractionSharedpTA;         ///!<! histogram pT jet (area based), fraction shared pT with Pythia jet
   
   AliAnalysisTaskPrepareInputForEmbedding(const AliAnalysisTaskPrepareInputForEmbedding&);
   AliAnalysisTaskPrepareInputForEmbedding &operator=(const AliAnalysisTaskPrepareInputForEmbedding&);
   
   ClassDef(AliAnalysisTaskPrepareInputForEmbedding, 1)
};
#endif


