#ifndef AliAnalysisTaskMCSmearing_cxx
#define AliAnalysisTaskMCSmearing_cxx

#include "AliAnalysisTaskSE.h"
class AliMCParticle;
class THn;
class AliTrackReference;

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
  // Check to see the abundance of pi0 and pich in a sample
  TH1F *fPiCheck;
  
  // Distribution of observed particles relative to their primary particle
  THn *fNsecondaries; //!
  // Efficiency of various particle species to produces hits on the FMD
  THn *fNprimaries; //!
  // Find the Track reference pointing to the hit on the FMD; NULL if there was no hit
  Bool_t IsHitFMD(AliMCParticle* p);
  // Get an iterable container of all the daughters of a given particle
  std::vector< AliMCParticle* > GetDaughters(AliMCParticle* p);
  // Get the number of hits which p's chain causes on the FMD
  Int_t ParticleProducedNHitsOnFMD(AliMCParticle* p);
  // Modified IsPhysicalPrimary check to regard pi0s as stable. Nont that this implementation
  // is not "orthogonal" to AliStack::IsFromWeakDecay and AliStack::IsSecondaryFromMaterial
  Bool_t IsRedefinedPhysicalPrimary(AliMCParticle* p);
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
  AliMCParticle* GetIncidentParticleFromFirstMaterialInteraction(AliMCParticle*);
  // Declaring these shuts up warnings from Weffc++
  AliAnalysisTaskMCSmearing(const AliAnalysisTaskMCSmearing&); // not implemented
  AliAnalysisTaskMCSmearing& operator=(const AliAnalysisTaskMCSmearing&); // not implemented

  ClassDef(AliAnalysisTaskMCSmearing, 1);
};

#endif
