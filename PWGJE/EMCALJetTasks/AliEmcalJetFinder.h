#ifndef ALIEMCALJETFINDER_H
#define ALIEMCALJETFINDER_H

// $Id$

namespace fastjet {
  class PseudoJet;
}

class AliEmcalJet;
class AliFJWrapper;
class TNamed;
class TH1;

class AliEmcalJetFinder : public TNamed
{
  public:
    AliEmcalJetFinder();
    AliEmcalJetFinder(const char* name);
    ~AliEmcalJetFinder();
    
    Bool_t                        FindJets();
    void                          AddInputVector(Double_t px, Double_t py, Double_t pz);
    void                          FillPtHistogram(TH1* histogram);
    void                          FillEtaHistogram(TH1* histogram);
    void                          FillPhiHistogram(TH1* histogram);

    Int_t                         GetJets(AliEmcalJet*& jets)     {jets = fJetArray[0]; return fJetCount;}
    AliEmcalJet*                  GetJet(Int_t index)             {return fJetArray[index];}
    void                          SetGhostArea(Double_t val)      {fGhostArea = val;}
    void                          SetRadius(Double_t val)         {fRadius = val;}
    void                          SetJetAlgorithm(Int_t val)      {fJetAlgorithm = val;}
    void                          SetTrackMaxEta(Double_t val)    {fTrackMaxEta = val;}
    void                          SetJetMaxEta(Double_t val)      {fJetMaxEta = val;}
    void                          SetJetMinPt(Double_t val)       {fJetMinPt = val;}
    void                          SetJetMinArea(Double_t val)     {fJetMinArea = val;}

    void                          SetManualIndex(Int_t val)       {fInputVectorIndex = val;}
  private:
    // General properties
    AliFJWrapper*                 fFastjetWrapper;
    Int_t                         fInputVectorIndex;
    Int_t                         fJetCount;
    std::vector<AliEmcalJet*>     fJetArray;
    // Settings for fastjet
    Double_t                      fGhostArea;
    Double_t                      fRadius;
    Int_t                         fJetAlgorithm;
    Double_t                      fTrackMaxEta;
    // Jet cuts
    Double_t                      fJetMaxEta;
    Double_t                      fJetMinPt;
    Double_t                      fJetMinArea;

    ClassDef(AliEmcalJetFinder, 1); // Lightweight fastjet implementation outside analysis tasks
};

#endif
