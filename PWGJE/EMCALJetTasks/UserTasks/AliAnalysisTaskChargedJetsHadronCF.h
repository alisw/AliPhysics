#ifndef ALIANALYSISTASKCHARGEDJETSHADRONCF_H
#define ALIANALYSISTASKCHARGEDJETSHADRONCF_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class THn;
class AliBasicParticle;
class AliAODPid;

//###############################################################################################################################################3
/**
 * \class AliChargedJetsHadronCFCuts
 * \brief Container class of cuts for AliAnalysisTaskChargedJetsHadronCF
 *
 *
 *
 * \author Ruediger Haake <ruediger.haake@cern.ch>, CERN
 * \date Nov 20, 2016
 */
// 
class AliChargedJetsHadronCFCuts
{
 public:
  AliChargedJetsHadronCFCuts() {;}
  AliChargedJetsHadronCFCuts(const char* cutName, const char* outputName, Double_t minPt, Double_t maxPt, Double_t minCent, Double_t maxCent, Double_t lowerOffset, Double_t offset2D, Double_t slope2D, Double_t jetVetoPt)
  {
    fCutName = cutName;
    fOutputName = outputName;
    fPtMin = minPt;
    fPtMax = maxPt;
    fCentMin = minCent;
    fCentMax = maxCent;

    fCutLowerPercentageOffset = lowerOffset;
    fCut2DOffset = offset2D;
    fCut2DSlope = slope2D;

    fJetVetoPt = jetVetoPt;

    fAcceptedJets = 0;
    fArrayIndex = -1;
  }
  Bool_t IsCutFulfilled(Double_t pt, Double_t mcPt, Double_t cent, Double_t ptRatio, Double_t vetoPt)
  {
    // Simple kinematic cuts
    if(pt < fPtMin || pt >= fPtMax)
      return kFALSE;
    if(cent < fCentMin || cent >= fCentMax)
      return kFALSE;

    // jet veto
    if(vetoPt >= fJetVetoPt)
      return kFALSE;

    // Lower MC percentage cut
    if(ptRatio < fCutLowerPercentageOffset)
      return kFALSE;

    // Value of the 2D cut line. We cut everything above that
    Double_t cutLineValue = fCut2DSlope*mcPt + fCut2DOffset;
    if(ptRatio > cutLineValue)
      return kFALSE;

    return kTRUE;
  }

  ~AliChargedJetsHadronCFCuts();
  TString     fCutName;                              ///< name of this cut
  TString     fOutputName;                           ///< array name that is used to output objects associated with this cuts
  Double_t    fPtMin;                                ///< valid for jets above this pT
  Double_t    fPtMax;                                ///< valid for jets below this pT
  Double_t    fCentMin;                              ///< valid for centralities above
  Double_t    fCentMax;                              ///< valid for centralities below
  Double_t    fJetVetoPt;                            ///< jet veto pt

  Double_t    fCutLowerPercentageOffset;             ///< cut value
  Double_t    fCut2DOffset;                          ///< cut value
  Double_t    fCut2DSlope;                           ///< cut value


  Int_t       fAcceptedJets;                         ///< temporary var that holds how many jets passed
  Int_t       fArrayIndex;                           ///< array index that holds the output array index
};

//###############################################################################################################################################3

/**
 * \class AliAnalysisTaskChargedJetsHadronCF
 * \brief Support task for (charged) jet-hadron correlations
 *
 * This task selects the jets and optionally the particles which are saved as a TClonesArray to the event
 * and which are read by a correlation task. The tasks does a special event selection (e.g. dijet events)
 * and jet selection (e.g. fakefactor rejection). In addition, it fills auxiliary histograms
 *
 * \author Ruediger Haake <ruediger.haake@cern.ch>, CERN
 * \date Apr 21, 2016
 */
// 
class AliAnalysisTaskChargedJetsHadronCF : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskChargedJetsHadronCF();
  AliAnalysisTaskChargedJetsHadronCF(const char *name);
  virtual ~AliAnalysisTaskChargedJetsHadronCF();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

  // ######### SETTERS/GETTERS
  void                        ActivateJetEmbedding(const char* embArray, const char* embTrackArray, Double_t maxDistance, Int_t numMatchedJets)
  {
    fJetEmbeddingMaxDistance      = maxDistance;
    fJetEmbeddingArrayName        = embArray;
    fJetEmbeddingTrackArrayName   = embTrackArray;
    fJetEmbeddingNumMatchedJets   = numMatchedJets;
  }
  void                        AddJetEmbeddingCut(const char* cutName, const char* outputName, Double_t minPt, Double_t maxPt, Double_t minCent, Double_t maxCent, Double_t lowerOffset, Double_t offset2D, Double_t slope2D, Double_t jetVetoPt)
  {
    AliChargedJetsHadronCFCuts tmpCut(cutName, outputName, minPt, maxPt, minCent, maxCent, lowerOffset, offset2D, slope2D, jetVetoPt);
    fJetEmbeddingCuts.push_back(tmpCut);
  }

  void                        SetJetVetoArrayName(const char* name)   { fJetVetoArrayName = name; }
  void                        SetJetVetoJetByJet(Bool_t val)          { fJetVetoJetByJet = val; }
  void                        SetUsePerTrackMCPercentage(Bool_t val)   { fJetEmbeddingUsePerTrackMCPercentage = val; }
  void                        SetUseBgrdForMCPercentage(Bool_t val)   { fJetEmbeddingUseBgrdForMCPercentage = val; }
  void                        SetCreateEmbeddingPtPlotPerCut(Bool_t val)   { fJetEmbeddingCreatePtPlotPerCut = val; }
  void                        SetConstPtFilterBit(Int_t val)   { fConstPtFilterBit = val; }
  void                        SetNumberOfCentralityBins(Int_t val)   { fNumberOfCentralityBins = val; }
  void                        SetJetParticleArrayName(const char* name)   { fJetParticleArrayName = name; }
  void                        SetTrackParticleArrayName(const char* name) { fTrackParticleArrayName = name; }

  void                        SetJetOutputMode(Int_t mode) {fJetOutputMode = mode;}
  void                        ActivateEventExtraction(Double_t percentage, Double_t minJetPt, Double_t maxJetPt) {fEventExtractionPercentage = percentage; fEventExtractionMinJetPt = minJetPt; fEventExtractionMaxJetPt = maxJetPt;}
  void                        SetTrackExtractionPercentagePower(Double_t val)   { fTrackExtractionPercentagePower = val; }

  void                        SetNumRandomConesPerEvent(Int_t val)   { fNumRandomConesPerEvent = val; }

  void                        SetUseConstituents(Bool_t data, Bool_t mc)   { fUseDataConstituents = data; fUseMCConstituents = mc;}

 protected:
  void                        ExecOnce();
  Bool_t                      Run();

  // ######### META FUNCTIONS
  void                        BinLogAxis(const THn *h, Int_t axisNumber);
  void                        AddEventToTree();
  void                        AddJetToOutputArray(AliEmcalJet* jet, Int_t arrayIndex, Int_t& jetsAlreadyInArray);
  void                        AddTrackToOutputArray(AliVTrack* track);
  void                        AddTrackToTree(AliVTrack* track);
  void                        FillHistogramsTracks(AliVTrack* track);
  void                        FillHistogramsJets(AliEmcalJet* jet, const char* cutName);
  Bool_t                      IsJetSelected(AliEmcalJet* jet);

  AliJetContainer            *fJetsCont;                                //!<! Jets
  AliTrackContainer          *fTracksCont;                              //!<! Tracks
  Double_t                    fEventExtractionPercentage;               ///< percentage of events that is recorded
  Double_t                    fEventExtractionMinJetPt;                 ///< minimum jet pt of recorded events
  Double_t                    fEventExtractionMaxJetPt;                 ///< maximum jet pt of recorded events
  
  Int_t                       fConstPtFilterBit;                        ///< For const pt plot, filter bit
  Int_t                       fNumberOfCentralityBins;                  ///< Number of centrality bins
  std::vector<TClonesArray*>  fJetsOutput;                              //!<! vector of arrays of basic correlation particles attached to the event (jets)
  TClonesArray               *fTracksOutput;                            //!<! Array of basic correlation particles attached to the event (tracks)
  TString                     fJetParticleArrayName;                    ///< Name of fJetsOutput array (if one uses only one)
  TString                     fTrackParticleArrayName;                  ///< Name of fTracksOutput array
  TClonesArray               *fJetEmbeddingArray;                       //!<! Array of generated jets imported into task (for embedding)
  TString                     fJetEmbeddingArrayName;                   ///< Name of array used to match jets
  TString                     fJetEmbeddingTrackArrayName;              ///< Name of array used to match tracks of jets
  Double_t                    fJetEmbeddingMaxDistance;                 ///< Max distance allowed to accept an embedded jet
  Int_t                       fJetEmbeddingNumMatchedJets;              ///< Number of matched leading jets that will be used
  Bool_t                      fJetEmbeddingUsePerTrackMCPercentage;     ///< When cutting on MC percentage, calculate it per track and not for all MC tracks
  Bool_t                      fJetEmbeddingUseBgrdForMCPercentage;      ///< When cutting on MC percentage, use bgrd. corr to calculate MC percentage
  Bool_t                      fJetEmbeddingCreatePtPlotPerCut;          ///< create TH3 per cut or only once
  std::vector<AliChargedJetsHadronCFCuts> fJetEmbeddingCuts;            ///< Cuts used in jet embedding

  TClonesArray               *fJetVetoArray;                            //!<! Array of jets imported into task used for veto a matching/embedding
  TString                     fJetVetoArrayName;                        ///< Name of array used for veto jets
  Bool_t                      fJetVetoJetByJet;                         ///< If true, the jet veto will be applied on a jet-by-jet basis
  std::vector<AliEmcalJet*>   fMatchedJets;                             ///< Jets matched in an event (embedded)
  std::vector<AliEmcalJet*>   fMatchedJetsReference;                    ///< Jets matched in an event (reference)
  TRandom3*                   fRandom;                                  ///< random number generator
  TTree*                      fTracksTree;                              //!<! Tree of extracted jets
  AliBasicParticle*           fTreeBufferTrack;                         //!<! Tree of extracted jets (buffer)
  AliAODPid*                  fTreeBufferPID;                           //!<! Tree of extracted jets (buffer)
  Int_t                       fTreeBufferPDG;                           //!<! Tree of extracted jets (buffer)
  Double_t                    fTrackExtractionPercentagePower;          ///< Extraction percentage for tracks
  Int_t                       fNumRandomConesPerEvent;                  ///< Number of random cones thrown in one event

  Bool_t                      fUseDataConstituents;                     ///< If true, tracks with labels <  10000 will be processed
  Bool_t                      fUseMCConstituents;                       ///< If true, tracks with labels >= 10000 will be processed

  // Criteria for the selection of jets that are passed to the correlation task
  Int_t                       fJetOutputMode;                           ///< mode which jets are written to array (0: all accepted, 1: leading,  2: subleading, 3: leading+subleading)

  // Event properties
  AliEmcalJet*                fLeadingJet;                              //!<!  leading jet (calculated event-by-event)
  AliEmcalJet*                fSubleadingJet;                           //!<!  subleading jet (calculated event-by-event)
  AliEmcalJet*                fInitialPartonMatchedJet1;                //!<!  On PYTHIA data and fJetOutputMode=6, this holds the PDG code of the initial collisions that was matched to this jet
  AliEmcalJet*                fInitialPartonMatchedJet2;                //!<!  On PYTHIA data and fJetOutputMode=6, this holds the PDG code of the initial collisions that was matched to this jet
  Int_t                       fAcceptedJets;                            //!<!  number accepted jets (calculated event-by-event)
  Int_t                       fAcceptedTracks;                          //!<!  number accepted tracks (calculated event-by-event)

  // ######### HISTOGRAM FUNCTIONS
  void                        FillHistogram(const char * key, Double_t x);
  void                        FillHistogram(const char * key, Double_t x, Double_t y);
  void                        FillHistogram(const char * key, Double_t x, Double_t y, Double_t add);
  void                        FillHistogram3D(const char * key, Double_t x, Double_t y, Double_t z, Double_t add = 0);

  template <class T> T*       AddHistogram1D(const char* name = "CustomHistogram", const char* title = "NO_TITLE", const char* options = "", Int_t xBins = 100, Double_t xMin = 0.0, Double_t xMax = 20.0, const char* xTitle = "x axis", const char* yTitle = "y axis");
  template <class T> T*       AddHistogram2D(const char* name = "CustomHistogram", const char* title = "NO_TITLE", const char* options = "", Int_t xBins = 100, Double_t xMin = 0.0, Double_t xMax = 20.0, Int_t yBins = 100, Double_t yMin = 0.0, Double_t yMax = 20.0,  const char* xTitle = "x axis", const char* yTitle = "y axis", const char* zTitle = "z axis");
  template <class T> T*       AddHistogram3D(const char* name = "CustomHistogram", const char* title = "NO_TITLE", const char* options = "", Int_t xBins = 100, Double_t xMin = 0.0, Double_t xMax = 20.0, Int_t yBins = 100, Double_t yMin = 0.0, Double_t yMax = 20.0, Int_t zBins = 100, Double_t zMin = 0.0, Double_t zMax = 20.0, const char* xTitle = "x axis", const char* yTitle = "y axis", const char* zTitle = "z axis");

  // ######### HELPER FUNCTIONS
  void                        CalculateJetType(AliEmcalJet* jet, Int_t& typeIC, Int_t& typeHM);
  void                        GetTrackMCRatios(AliEmcalJet* jet, AliEmcalJet* mcJet, Double_t& trackRatio, Double_t& ptRatio);
  AliEmcalJet*                GetVetoJet(AliEmcalJet* jet);
  AliEmcalJet*                GetLeadingVetoJet();

  AliEmcalJet*                GetReferenceJet(AliEmcalJet* jet);
  Bool_t                      IsTrackInCone(AliVParticle* track, Double_t eta, Double_t phi, Double_t radius);
  void                        GetInitialCollisionJets();
  void                        GetMatchingJets();

  void                        GetLeadingJets(const char* opt, AliEmcalJet*& jetLeading, AliEmcalJet*& jetSubLeading);
  void                        GetLeadingJetsInArray(TClonesArray* arr, const char* opt, AliEmcalJet*& jetLeading, AliEmcalJet*& jetSubLeading);
  void                        CalculateEventProperties();


 private:
  AliAnalysisTaskChargedJetsHadronCF(const AliAnalysisTaskChargedJetsHadronCF&);            // not implemented
  AliAnalysisTaskChargedJetsHadronCF &operator=(const AliAnalysisTaskChargedJetsHadronCF&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskChargedJetsHadronCF, 12) // Charged jet+h analysis support task
  /// \endcond
};



#endif
