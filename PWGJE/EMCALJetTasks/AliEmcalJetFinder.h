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
    void                          AddInputVector(Double_t px, Double_t py, Double_t pz, Double_t E);
    void                          FillPtHistogram(TH1* histogram);
    void                          FillEtaHistogram(TH1* histogram);
    void                          FillPhiHistogram(TH1* histogram);

    Int_t                         GetJets(AliEmcalJet*& jets)     {jets = fJetArray[0]; return fJetCount;}
    AliEmcalJet*                  GetJet(Int_t index)             {return fJetArray[index];}
    Int_t                         GetNumberOfJets()               {return fJetCount;}
    void                          SetGhostArea(Double_t val)      {fGhostArea = val;}
    void                          SetRadius(Double_t val)         {fRadius = val;}
    void                          SetJetAlgorithm(Int_t val)      {fJetAlgorithm = val;}
    void                          SetRecombSheme(Int_t val)       {fRecombScheme = val;}
    void                          SetTrackMaxEta(Double_t val)    {fTrackMaxEta = val;}
    void                          SetJetMaxEta(Double_t val)      {fJetMaxEta = val;}
    void                          SetJetMinPt(Double_t val)       {fJetMinPt = val;}
    void                          SetJetMinArea(Double_t val)     {fJetMinArea = val;}

    void                          SetManualIndex(Int_t val)       {fInputVectorIndex = val;}
  private:
    // General properties
    AliFJWrapper*                 fFastjetWrapper;                // Interface object to fastjet
    Int_t                         fInputVectorIndex;              // Current index of input vectors (by default: count of vectors)
    Int_t                         fJetCount;                      // Found jets within the given acceptances
    std::vector<AliEmcalJet*>     fJetArray;                      // Internal array for the jets
    // Settings for fastjet
    Double_t                      fGhostArea;                     // setting for ghost area in FJ
    Double_t                      fRadius;                        // Radius parameter
    Int_t                         fJetAlgorithm;                  // var for algorithm (0=antikt, 1=kt)
    Int_t                         fRecombScheme;                  // Recombination scheme for Fastjet
    Double_t                      fTrackMaxEta;                   // cut for |track-eta| < fTrackMaxEta
    // Jet cuts
    Double_t                      fJetMaxEta;                     // cut for |jet-eta| < fJetMaxEta
    Double_t                      fJetMinPt;                      // cut for  jet-pT > fJetMinPt
    Double_t                      fJetMinArea;                    // cut for  jet-area > fJetMinArea 

    ClassDef(AliEmcalJetFinder, 1); // Lightweight fastjet implementation outside analysis tasks
};

#endif
