/// \class AliAnalysisTaskPrepareInputForEmbedding
/// \brief Save in a TTree the 4-momentum vector of the jets at reco and particle level
///
/// \ingroup PWGJEUSER
/// Reads the reconstructed jet container (no background subtraction since it's thought to run on PYTHIA)
/// The correspondent particle level jet is also stored. 
/// The task requires the AliAnalysisTaskEmcalJetTagger for PYTHIA matching to be run first
/// The output tree is thought as an input to the single track embedding
///
/// \author Chiara Bianchin
/// \date Aug, 2015

#ifndef ALIANALYSISTASKPREPAREINPUTFOREMBEDDING_H
#define ALIANALYSISTASKPREPAREINPUTFOREMBEDDING_H

class TLorentzVector;
class THnSparse;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskPrepareInputForEmbedding : public AliAnalysisTaskEmcalJet {

public:
   AliAnalysisTaskPrepareInputForEmbedding();
   AliAnalysisTaskPrepareInputForEmbedding(const char* name);
   virtual ~AliAnalysisTaskPrepareInputForEmbedding();

   void        UserCreateOutputObjects();
   void        Terminate(Option_t *option);
   
   //setters
   void        SetContainerNumber(Int_t number)        {fContainer = number;}
   void        SetMinFracShared(Double_t minfrac)      {fMinFractionShared = minfrac;}
   void        SetLeadingJetOnly(Bool_t b = kTRUE)     {fLeadingJetOnly = b;}
   void        SetDetHardCoreTagging(Bool_t b = kTRUE) {fHardCoreTag = b;}
   
protected:
   Bool_t      Run();
   Bool_t      FillHistograms();

private:
   Int_t        fContainer;                   ///< reco jet container number
   Double_t     fMinFractionShared;           ///< minimum pT fraction shared between reco and part level jets
   Bool_t       fLeadingJetOnly;              ///< fill TTrees with leading jet only
   Bool_t       fHardCoreTag;                 ///< Request also hard-core tagging of the det level jet (requires tagging task to be run). Default is off
   TTree        *fTreeJets;                   //!<! tree with the TLorentzVector of the jet detector and particle level
   TLorentzVector *fJetDet;                   //!<! reconstucted jets
   TLorentzVector *fJetPart;                  //!<! particle level jet
   TLorentzVector *fJetDetL;                  //!<! reconstucted leading jet
   TLorentzVector *fJetPartL;                 //!<! particle level leading jet
   TH1F         *fNumberOfJets;               //!<! histograms with number of jets selected
   TH2F         *fhFractionSharedpT;         //!<! histogram pT jet (area based), fraction shared pT with Pythia jet
   TH1F         *fNAccJets;                  //!<! number of jets per event
   THnSparse    *fhResponse;                 //!<! response matrix
   
   AliAnalysisTaskPrepareInputForEmbedding(const AliAnalysisTaskPrepareInputForEmbedding&);
   AliAnalysisTaskPrepareInputForEmbedding &operator=(const AliAnalysisTaskPrepareInputForEmbedding&);
   
   ClassDef(AliAnalysisTaskPrepareInputForEmbedding, 4)
};
#endif


