#ifndef AliAnalysisTaskMCSmearing_cxx
#define AliAnalysisTaskMCSmearing_cxx

#include "AliAnalysisTaskSE.h"
class AliMCParticle;
class THn;
class AliTrackReference;
class TParticle;
class TCutG;
class AliAnalysisTaskMCSmearing : public AliAnalysisTaskSE {

 public:
  AliAnalysisTaskMCSmearing();
  AliAnalysisTaskMCSmearing(const char *name);
  virtual ~AliAnalysisTaskMCSmearing() {};

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);

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
  // Simple dN/deta of particles hitting the ITS or FMD
  TH1F *fdNdeta;

  // Event counter with z-pos
  TH1F *fEventCounter;
  
  // Check to see the abundance of pi0 and pich in a sample
  TH1F *fPiCheck;

  // dN/deta distribution based on origin of particles
  TH2F *fdNdetaOrigin;
  // X-ray plot showing the origin of secondary particles
  TH2F *fxray;
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
  // Check if a given particle itself hit the FMD. If so, return the
  // (first) track reference of such a hit
  AliTrackReference* IsHitFMD(AliMCParticle* p);

  // Check if a given particle itself hit the FMD. If so, return the
  // (first) track reference of such a hit
  AliTrackReference* IsHitITS(AliMCParticle* p);
  
  // Get an iterable container of all the daughters of a given particle
  std::vector< AliMCParticle* > GetDaughters(AliMCParticle* p);
  // Get the number of hits which p's chain causes on the FMD
  Int_t ParticleProducedNHitsOnFMD(AliMCParticle* p);
  // Modified IsPhysicalPrimary check to regard pi0s as stable. Nont that this implementation
  // is not "orthogonal" to AliStack::IsFromWeakDecay and AliStack::IsSecondaryFromMaterial
  Bool_t IsRedefinedPhysicalPrimary(AliMCParticle* p);

  // Returns the region where the given particle's vertex is
  // located. Region is encoded in cOriginType enum
  Int_t GetOriginType(AliMCParticle *p);

protected:
  TList *fOutputList;  //!
  // Find the primary particle of a decay chain. If `p` is alreay the primary return p.
  // If it was not possible to find the mother, return NULL.
  AliMCParticle* GetMother(AliMCParticle* p);

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

  // Setup detector region cuts
  void SetupCuts();
  // Declaring these shuts up warnings from Weffc++
  AliAnalysisTaskMCSmearing(const AliAnalysisTaskMCSmearing&); // not implemented
  AliAnalysisTaskMCSmearing& operator=(const AliAnalysisTaskMCSmearing&); // not implemented

  ClassDef(AliAnalysisTaskMCSmearing, 1);
};

#endif
