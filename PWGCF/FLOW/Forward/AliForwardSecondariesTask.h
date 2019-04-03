//
// Calculate the flow in the forward regions using the Q cumulants method
//
#ifndef AliForwardSecondariesTask_H
#define AliForwardSecondariesTask_H
/**
 * @file AliForwardSecondariesTask.h
 * @author Freja Thoresen <freja.thoresen@cern.ch>
 *
 * @brief
 *
 * @ingroup pwgcf_forward_flow
 */
#include "AliAnalysisTaskSE.h"
#include "TString.h"
#include <TH2D.h>
#include <TH3F.h>
#include "TRandom.h"
#include "AliForwardSettings.h"
#include "AliForwardFlowUtil.h"
#include "AliFMDMCTrackDensity.h"
//#include "AliEventCuts.h"
#include <TF1.h>
class AliMCParticle;
class THn;
class AliTrackReference;
class TParticle;
class TCutG;
class AliAODForwardMult;
class TH2D;
class AliESDEvent;
/**
 * @defgroup pwglf_forward_tasks_flow Flow tasks
 *
 * Code to do flow
 *
 * @ingroup pwglf_forward_tasks
 */
/**
 * Calculate the flow in the forward regions using the Q cumulants method
 *
 * @par Inputs:
 *   - AliAODEvent
 *
 * Outputs:
 *   - forward_flow.root
 *
 * @ingroup pwglf_forward_flow
 *
 */
class AliForwardSecondariesTask : public AliAnalysisTaskSE
{
public:
  /**
   * Constructor
   */
  AliForwardSecondariesTask();
  /**
   * Constructor
   *
   * @param name Name of task
   */
  AliForwardSecondariesTask(const char* name);
  /**
   * Destructor
   */
  virtual ~AliForwardSecondariesTask() {}

  /**
   * Copy constructor
   *
   * @param o Object to copy from
   */
  AliForwardSecondariesTask(const AliForwardSecondariesTask& o);

  /**
   * @{
   * @name Task interface methods
   */

  /**
   * Create output objects
   */
  virtual void UserCreateOutputObjects();


  /**
   * Process each event
   *
   * @param option Not used
   */
  virtual void UserExec(Option_t *option);

  /**
   * End of job
   *
   * @param option Not used
   */
  virtual void Terminate(Option_t *option);

  // Enums for possible interaction/decay chains
  enum {kDIRECTMATERIAL,
  kACCUMULATEDMATERIAL,
  kDIRECTWEAK,
  kACCUMULATEDWEAK,
  kCOMPAREDTOPHYSICALPRIMARY,
  kCOMPAREDTOCHARGEDPHYSICALPRIMARY,
  kCOMPAREDTONEUTRALPHYSICALPRIMARY,
  kCOMPAREDTOPI0PRIMARY,
  kEXTENDEDISPHYSICALPRIMARY,
  kNCHAINS};
  // Species of the primary particle. Make sure that they are orthogonal!
  struct cPrimaryType {
    enum {kPI0,
    kPICHARGED,
    kOTHERS,
    kNSPECIES};
  };
  enum {kMCTRUTH,
  kMCRECON,
  kNMCMODES};
  struct cOriginType {
    enum {kPRIMARY,
    kEARLYDECAY,
    kPIPE,
    kITS,
    kFMD,
    kOTHER,
    kNORIGINTYPES
    };
  };
Double_t WrapPi(Double_t phi);
  Double_t GetTrackReferenceEta(AliTrackReference* tr);

  Bool_t AddMotherIfFirstTimeSeen(AliMCParticle* p, std::vector<Int_t> v);


  // Check if a given particle itself hit the FMD. If so, return the
  // (first) track reference of such a hit
  AliTrackReference* IsHitITS(AliMCParticle* p);

  // Get an iterable container of all the daughters of a given particle
  //std::vector< AliMCParticle* > GetDaughters(AliMCParticle* p);
  // Get the number of hits which p's chain causes on the FMD
  //Int_t ParticleProducedNHitsOnFMD(AliMCParticle* p);
  // Modified IsPhysicalPrimary check to regard pi0s as stable. Nont that this implementation
  // is not "orthogonal" to AliStack::IsFromWeakDecay and AliStack::IsSecondaryFromMaterial
  Bool_t IsRedefinedPhysicalPrimary(AliMCParticle* p);

  // check if a particle has a history of material interaction
  Bool_t hasParticleMaterialInteractionInAncestors(AliMCParticle* p);

  // Returns the region where the given particle's vertex is
  // located. Region is encoded in cOriginType enum
  Int_t GetOriginType(AliMCParticle *p);

  //private:
  TList*                  fOutputList;    //! output list
  //TList* fDiffList; //!
  TList* fEventList; //!
  TList* fDeltaList; //!
  TRandom fRandom;
  AliFMDMCTrackDensity* fTrackDensity; //!

  // A class combining all the settings for this analysis
  AliForwardSettings fSettings;
  AliForwardFlowUtil fUtil;

TF1 *fMultTOFLowCut; //!
TF1 *fMultTOFHighCut; //!
TF1 *fMultCentLowCut; //!


  // Simple dN/deta of particles hitting the ITS or FMD
  TH1F *fdNdeta;//!

  // Check to see the abundance of pi0 and pich in a sample
  TH1F *fPiCheck;//!

  // dN/deta distribution based on origin of particles
  TH2F *fdNdetaOrigin;//!
  // X-ray plot showing the origin of secondary particles
  TH2F *fxray;//!
  // Distribution of observed particles relative to their primary particle
  THn *fNsecondaries; //!
  // Efficiency of various particle species to produces hits on the FMD
  THn *fNprimaries; //!


  // Cuts for various detector regions
  TCutG *fITS;  //!
  TCutG *fFMD1;  //!
  TCutG *fFMD2;  //!
  TCutG *fFMD3;  //!
  TCutG *fPipe;  //!
  TCutG *fEarlyDecay;  //!
  TH1D* phihist; //!
  AliTrackReference* fStored; //! Last stored
  enum {
    kTPCOnly = 128, // TPC only tracks
    kHybrid = 768, // TPC only tracks
    kphiAcceptanceBin = 21 // phi acceptance bin in the FMD histogram (dNdetadphi)
  };

  static Double_t Wrap02pi(Double_t angle);
  static Double_t Mod(Double_t x, Double_t y);

protected:

  /**
   * Structure holding the state of the `tracker'
   *
  */
  mutable struct State
  {
    Double_t angle;            // Angle
    UShort_t oldDetector;      // Last detector
    Char_t   oldRing;          // Last ring
    UShort_t oldSector;        // Last sector
    UShort_t oldStrip;         // Last strip
    UShort_t startStrip;       // First strip
    UShort_t nRefs;            // Number of references
    UShort_t nStrips;          // Number of strips
    UShort_t count;            // Count of hit strips
    AliTrackReference* longest; //! Longest track through
    /**
     * Clear this state
     *
     * @param alsoCount If true, also clear count
    */
    void Clear(Bool_t alsoCount=false);
    /**
     * Assignment operator
     *
     * @param o Object to assign from
     *
     * @return Reference to this object
     */
    State& operator=(const State& o);
  } fState; //! State
  UShort_t fMaxConsequtiveStrips;
  Double_t fLowCutvalue;
  Bool_t            fTrackGammaToPi0;
  AliTrackReference*  ProcessRef(AliMCParticle*       particle,
    AliMCParticle* mother,
    AliTrackReference*   ref,std::vector< Int_t > listOfMothers, Double_t randomInt, Float_t event_vtx_z, Double_t v0cent);

  void BeginTrackRefs();
  void EndTrackRefs();

  void StoreParticle(AliMCParticle*       particle,
	  AliMCParticle* mother,
	  AliTrackReference*   ref,std::vector< Int_t > listOfMothers, Double_t randomInt, Float_t event_vtx_z, Double_t v0cent) ;


  Bool_t ProcessTrack(AliMCParticle* particle, AliMCParticle* mother, std::vector< Int_t > listOfMothers, Double_t randomInt, Float_t event_vtx_z, Double_t v0cent);

  Double_t GetTrackRefTheta(const AliTrackReference* ref) const;

  // Find the primary particle of a decay chain. If `p` is alreay the primary return p.
  // If it was not possible to find the mother, return NULL.
  AliMCParticle* GetMother(AliMCParticle* p);
  AliMCParticle* GetMother(Int_t iTr, const AliMCEvent* event) const;
  // Find the primary particle of a decay chain if it is charged.
  // If `p` is alreay the primary return p. If it was not possible
  // to find the mother or if the mother was not charged, return NULL
  AliMCParticle* GetChargedMother(AliMCParticle*);

  // Complimentary to `GetChargedMother`
  AliMCParticle* GetNeutralMother(AliMCParticle*);

  // Return pi0 if decay chain terminated with pi0 -> gamma where gamma IsPhysicalPrimary.
  // Else, return NULL
  AliMCParticle* GetPi0Mother(AliMCParticle*);

  // Re-define the pi0 a physical primary. This returns the IsPhysicalPrimary particle
  // or a pi0 if the IsPhysicalPrimary was a gamma from pi0 -> gamma + gamma
  // AliMCParticle* GetMotherExtendedPrimaryDef(AliMCParticle* p);

  // Get the eta and phi coordinates where the FMD track reference was created
  // etaPhi is a 2-element array where the values will be written in.
  // If no FMD-reference was created, etaPhi will be NULL
  void GetTrackRefEtaPhi(AliMCParticle* p, Double_t* etaPhi);
  void GetTrackRefEtaPhi(AliTrackReference* ref, Double_t* etaPhi);

  // Get the phi coordinate where the track ref was created
  // Double_t GetTrackRefPhi(AliTrackReference* ref);
  // Get the eta coordinate where the track ref was created
  // Double_t GetTrackRefEta(AliTrackReference* ref);

  // Find the first particle that interacted with material. If there was no material
  // interaction return NULL.
  AliMCParticle* GetIncidentParticleFromFirstMaterialInteraction(AliMCParticle* p);

  // Find the first non-primary particle mother. Returns NULL if no
  // such mother was found (includes the case where the given particle
  // is already a primary
  AliMCParticle* GetFirstNonPrimaryMother(AliMCParticle* p);


  ClassDef(AliForwardSecondariesTask, 1); // Analysis task for flow analysis
};

#endif
