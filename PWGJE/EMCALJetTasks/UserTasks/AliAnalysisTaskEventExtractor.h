#ifndef ALIANALYSISTASKEVENTEXTRACTOR_H
#define ALIANALYSISTASKEVENTEXTRACTOR_H

/* Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//###############################################################################################################################################3
class AliEventTree;

/**
 * \class AliAnalysisTaskEventExtractor
 * \brief Analysis task that implements AliEventTree to extract events containing MC particles + associated rec. particles to a tree
 *
 * Usage: Task will by default extract particles + tracks.
 * To adjust event extraction percentage, use task->SetEventPercentage()
 *
 * \author Ruediger Haake <ruediger.haake@cern.ch>, Yale
 * \date Apr 16, 2019
 */
// 

class AliAnalysisTaskEventExtractor : public AliAnalysisTaskEmcal {
 public:
  AliAnalysisTaskEventExtractor();
  AliAnalysisTaskEventExtractor(const char *name);
  virtual ~AliAnalysisTaskEventExtractor();
  static AliAnalysisTaskEventExtractor* AddTaskEventExtractor(TString trackArray, TString particleArray, const char* taskNameSuffix);
  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);
  AliEventTree*               GetTree() {return fTree;}

  void                        SetSaveTracks(Bool_t val) {fSaveTracks = val;}
  void                        SetCustomStartupScript(const char* path)            { fCustomStartupScript = path; }
  void                        SetEventPercentage(Double_t val)                    { fEventPercentage  = val; }
  void                        SetRandomSeed(ULong_t val)                          { fRandomSeed  = val; }

 protected:
  void                        ExecOnce();
  Bool_t                      Run();

  // ################## CUTS AND SETTINGS: Should be set during task initialization

  Bool_t                      fSaveTracks;                              ///< save arrays of tracks for each particle
  Double_t                    fEventPercentage;                         ///< percentage (0, 1] which will be extracted
  ULong_t                     fRandomSeed;                              ///< random seed
  TString                     fCustomStartupScript;                     ///< startup script
  // ################## BUFFERS: Variables set during execution time
  AliEventTree*               fTree;
  Int_t                       fMultiplicity;                            ///< Multiplicity (number tracks, also for multiple containers)
  TRandom3*                   fRandomGenerator;                         //!<! Random number generator, used for event + jet efficiency

  // ################## HELPER FUNCTIONS
  Double_t                    GetDistance(Double_t eta1, Double_t eta2, Double_t phi1, Double_t phi2)
  {
    Double_t deltaPhi = TMath::Min(TMath::Abs(phi1-phi2),TMath::TwoPi() - TMath::Abs(phi1-phi2));
    return TMath::Sqrt((eta1-eta2)*(eta1-eta2) + deltaPhi*deltaPhi);
  }
  void                        FillHistogram(const char * key, Double_t x);
  void                        FillHistogram(const char * key, Double_t x, Double_t y);
  void                        FillHistogram(const char * key, Double_t x, Double_t y, Double_t add);
  void                        FillHistogram3D(const char * key, Double_t x, Double_t y, Double_t z, Double_t add = 0);
  template <class T> T*       AddHistogram1D(const char* name = "CustomHistogram", const char* title = "NO_TITLE", const char* options = "", Int_t xBins = 100, Double_t xMin = 0.0, Double_t xMax = 20.0, const char* xTitle = "x axis", const char* yTitle = "y axis");
  template <class T> T*       AddHistogram2D(const char* name = "CustomHistogram", const char* title = "NO_TITLE", const char* options = "", Int_t xBins = 100, Double_t xMin = 0.0, Double_t xMax = 20.0, Int_t yBins = 100, Double_t yMin = 0.0, Double_t yMax = 20.0,  const char* xTitle = "x axis", const char* yTitle = "y axis", const char* zTitle = "z axis");
  template <class T> T*       AddHistogram3D(const char* name = "CustomHistogram", const char* title = "NO_TITLE", const char* options = "", Int_t xBins = 100, Double_t xMin = 0.0, Double_t xMax = 20.0, Int_t yBins = 100, Double_t yMin = 0.0, Double_t yMax = 20.0, Int_t zBins = 100, Double_t zMin = 0.0, Double_t zMax = 20.0, const char* xTitle = "x axis", const char* yTitle = "y axis", const char* zTitle = "z axis");

 private:
  AliAnalysisTaskEventExtractor(const AliAnalysisTaskEventExtractor&);            // not implemented
  AliAnalysisTaskEventExtractor &operator=(const AliAnalysisTaskEventExtractor&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEventExtractor, 1) // Track extraction task
  /// \endcond
};

//###############################################################################################################################################3
/**
 * \class AliEventTree
 * \brief Class managing creation of a tree containing certain parts of events
 *        The tree is intended to be used for ML tasks
 *
 *
 * \author Ruediger Haake <ruediger.haake@cern.ch>, Yale
 * \date Apr 16, 2019
 */
// 
class AliEventTree : public TNamed
{
  public:
    AliEventTree();
    AliEventTree(const char* name);
    void            InitializeTree(Bool_t saveTracks);
    Bool_t          AddEventToTree(Long64_t eventID, std::vector<Double_t> vertex, std::vector<AliVParticle*> particles, std::vector<AliVTrack*> tracks);
    TTree*          GetTreePointer() {return fTree;}
  protected:
    TTree*          fTree;                                //!<! tree structure
    Bool_t          fInitialized;                         ///< init state of tree

    // Buffers that will be added to the tree
    Long64_t        fBuffer_Event_ID;                     //!<! array buffer
    Int_t           fBuffer_Event_Number_Particles;       //!<! array buffer
    Int_t           fBuffer_Event_Number_Tracks;          //!<! array buffer
    Float_t         fBuffer_Event_Vertex_X;               //!<! array buffer
    Float_t         fBuffer_Event_Vertex_Y;               //!<! array buffer
    Float_t         fBuffer_Event_Vertex_Z;               //!<! array buffer

    Float_t*        fBuffer_Particles_Pt;                 //!<! array buffer
    Float_t*        fBuffer_Particles_Eta;                //!<! array buffer
    Float_t*        fBuffer_Particles_Phi;                //!<! array buffer
    Int_t*          fBuffer_Particles_Charge;             //!<! array buffer
    Int_t*          fBuffer_Particles_Label;              //!<! array buffer
    Int_t*          fBuffer_Particles_PdgCode;            //!<! array buffer
    Float_t*        fBuffer_Particles_E;                  //!<! array buffer
    Bool_t*         fBuffer_Particles_IsPrimary;          //!<! array buffer

    Float_t*        fBuffer_Tracks_Pt;                    //!<! array buffer
    Float_t*        fBuffer_Tracks_Eta;                   //!<! array buffer
    Float_t*        fBuffer_Tracks_Phi;                   //!<! array buffer
    Int_t*          fBuffer_Tracks_Charge;                //!<! array buffer
    Int_t*          fBuffer_Tracks_Label;                 //!<! array buffer
    Bool_t*         fBuffer_Tracks_IsGlobal;              //!<! array buffer

    /// \cond CLASSIMP
    ClassDef(AliEventTree, 1) // Event tree class
    /// \endcond
};

#endif
