#ifndef ALIDIELECTRONSIGNALMC_H
#define ALIDIELECTRONSIGNALMC_H

#include <TNamed.h>
#include <TMCProcess.h>

/*
   Ionut Cristian Arsene, iarsene@cern.ch
 */

/*
   Monte Carlo signal definition:
      Leg #1  <-- Mother #1  <--  Grandmother #1
                      |
      Leg #2  <-- Mother #2  <--  Grandmother #2

   For every leg, mother or grand-mother, a PDG code and a source can be specified.

   1.) For the PDG codes, the PYTHIA standard is used.
   A few non-existent PYTHIA codes are used to select more than one PYTHIA code. All these are described below
   and implemented in AliDielectronMC::ComparePDG() function:
      0 - default, accepts all PYTHIA codes
    100 - light unflavoured mesons in the code range 100-199
    200 -        --"--                               200-299
    300 - strange mesons in the code range           300-399
    400 - charmed mesons in the code range           400-499
    401 - open charm mesons (all D and D* mesons)    400-439
    402 - open charm mesons and baryons together     400-439, 4000-4399
    403 - all charm hadrons (mesons and baryons)     400-499, 4000-4999
    500 - beauty mesons in the code range            500-599
    501 - open beauty mesons                         500-549
    502 - open beauty mesons and baryons             500-549, 5000-5499
    503 - all beauty hadrons                         500-599, 5000-5999
    902 - all open charm open beauty mesons+baryons  400-439, 500-549, 4000-4399, 5000-5499
    903 - all hadrons in the code range              100-599, 1000-5999
   1000 - light unflavoured baryons in the code range 1000-1999
   2000 -        --"--                                2000-2999
   3000 - strange baryons in the code range           3000-3999
   4000 - charmed baryons in the code range           4000-4999
   4001 - open charm baryons                          4000-4399
   5000 - beauty baryons in the code range            5000-5999
   5001 - open beauty baryons                         5000-5499

   2.) If the exclusion flags are turned ON then the PDG codes required and the conventional codes described above
       are used to exclude the selected particles.

   3.) If the selection of both charges is switched ON then the PDG codes act on both particles and anti-particles.

   4.) Particles sources implemented:
     (Incomplete list, see AliDielectronMC::CheckParticleSource() for more details.)
     1. Primary   - particle originating in the physics event
     2. FinalState- stable(final state) particles which reach the detector -> according to AliStack::IsPhysicalPrimary()
     3. Direct    - primary particle which has no mother (e.g. J/psi's added to pythia MC events via generator cocktails,
                    particles generated in a sudden freeze-out in thermal models, initial state particles)
     4. Secondary - particle created during the GEANT propagation due to interaction of final state primaries with the material

   5.) The 2 legs can originate from the same or different mother particles. This can be specified via the SetMotherRelation()
       method call.

   6.) The filling of the pure MC step can be switched on using SetFillPureMCStep() method call. This should be used
       with care since at the pure MC information level there is no cut applied and for abundant particles the combinatorics
       can be very high.
*/


//__________________________________________________________________
class AliDielectronSignalMC : public TNamed {

 public:
  enum EBranchRelation {kUndefined=0, kSame, kDifferent};
  enum ESource {kDontCare=0, kPrimary, kFinalState, kDirect, kSecondary, kNoCocktail, kSecondaryFromWeakDecay, kSecondaryFromMaterial, kFromBGEvent, kFinalStateFromBGEvent};
  enum EJpsiRadiativ {kAll=0, kIsRadiative, kIsNotRadiative};
  enum EJpsiSignals {
    kBegin = 0,
    kInclusiveJpsi,
    kBeautyJpsi,
    kPromptJpsi,
    kPromptRadJpsi,
    kPromptNonRadJpsi,
    kDirectJpsi,
    kGammaConv,
    kGammaConvDiffMother,
    kElectrons,
    kDirectElectrons,
    kPrimaryElectrons,
    kFromBgTest,
    kEnd
  };


  AliDielectronSignalMC();
  AliDielectronSignalMC(const Char_t* name, const Char_t* title);
  virtual ~AliDielectronSignalMC();

  void SetLegPDGs(Int_t pdg1, Int_t pdg2, Bool_t exclude1=kFALSE, Bool_t exclude2=kFALSE)
    {fLeg1 = pdg1; fLeg2 = pdg2; fLeg1Exclude=exclude1; fLeg2Exclude=exclude2;}
  void SetMotherPDGs(Int_t pdg1, Int_t pdg2, Bool_t exclude1=kFALSE, Bool_t exclude2=kFALSE)
    {fMother1 = pdg1; fMother2 = pdg2; fMother1Exclude=exclude1; fMother2Exclude=exclude2;}
  void SetGrandMotherPDGs(Int_t pdg1, Int_t pdg2, Bool_t exclude1=kFALSE, Bool_t exclude2=kFALSE)
    {fGrandMother1 = pdg1; fGrandMother2 = pdg2; fGrandMother1Exclude=exclude1; fGrandMother2Exclude=exclude2;}
  //both used for niece discrimination, nieces are daughters of the sisters
  void SetLegSources(ESource s1, ESource s2)                       {fLeg1Source = s1;                      fLeg2Source = s2;}
  void SetMotherSources(ESource s1, ESource s2)                    {fMother1Source = s1;                   fMother2Source = s2;}
  void SetGrandMotherSources(ESource s1, ESource s2)               {fGrandMother1Source = s1;              fGrandMother2Source = s2;}
  void SetCheckBothChargesLegs(Bool_t flag1, Bool_t flag2)         {fCheckBothChargesLeg1 = flag1;         fCheckBothChargesLeg2 = flag2;}
  void SetCheckBothChargesMothers(Bool_t flag1, Bool_t flag2)      {fCheckBothChargesMother1 = flag1;      fCheckBothChargesMother2 = flag2;}
  void SetCheckBothChargesGrandMothers(Bool_t flag1, Bool_t flag2) {fCheckBothChargesGrandMother1 = flag1; fCheckBothChargesGrandMother2 = flag2;}
  void SetMothersRelation(EBranchRelation relation)                {fMothersRelation = relation;}
  void SetGrandMothersRelation(EBranchRelation relation)           {fGrandMothersRelation = relation;}
  void SetGEANTProcess(TMCProcess processID)                       {fGEANTProcess = processID; fCheckGEANTProcess=kTRUE;}
  void SetFillPureMCStep(Bool_t fill=kTRUE)                        {fFillPureMCStep = fill;}
  void SetCheckMotherGrandmotherRelation(Bool_t CheckMotherIsGrandmother=kTRUE, Bool_t MotherIsGrandmother=kFALSE)
    {fCheckMotherGrandmother = CheckMotherIsGrandmother; fMotherIsGrandmother = MotherIsGrandmother;}
  void SetCheckMotherGrandmotherDiffPairRelation(Bool_t CheckMotherIsGrandmother=kTRUE, Bool_t MotherIsGrandmother=kFALSE)    // is assuming that particles from pair one have same mother and particles from pair two have same mother (SetMothersRelation(AliDielectronSignalMC::kSame))
    {fCheckMotherGrandmotherDiffPair = CheckMotherIsGrandmother; fMotherIsGrandmotherDiffPair = MotherIsGrandmother;}
  void SetCheckStackForPDG(Bool_t checkStack=kTRUE)                {fCheckStackForPDG = checkStack;}
  void SetPDGforStack(Int_t stackPDG)                              {fStackPDG = stackPDG;}
  void SetCheckUnlikeSign(Bool_t checkULS)                         {fCheckUnlikeSign = checkULS;}
  void SetCheckCorrelatedHF(Bool_t checkHFcorr)                    {fCheckCorrelatedHF = checkHFcorr;}
  void SetCheckLikeSign(Bool_t checkLSpp,Bool_t checkLSmm)         {fCheckLikeSignPP = checkLSpp; fCheckLikeSignMM = checkLSmm;}

  Int_t GetLegPDG(Int_t branch)                                const {return (branch==1 ? fLeg1 : fLeg2);}
  Int_t GetMotherPDG(Int_t branch)                             const {return (branch==1 ? fMother1 : fMother2);}
  Int_t GetGrandMotherPDG(Int_t branch)                        const {return (branch==1 ? fGrandMother1 : fGrandMother2);}
  Bool_t GetLegPDGexclude(Int_t branch)                        const {return (branch==1 ? fLeg1Exclude : fLeg2Exclude);}
  Bool_t GetMotherPDGexclude(Int_t branch)                     const {return (branch==1 ? fMother1Exclude : fMother2Exclude);}
  Bool_t GetGrandMotherPDGexclude(Int_t branch)                const {return (branch==1 ? fGrandMother1Exclude : fGrandMother2Exclude);}
  ESource GetLegSource(Int_t branch)                           const {return (branch==1 ? fLeg1Source : fLeg2Source);}
  ESource GetMotherSource(Int_t branch)                        const {return (branch==1 ? fMother1Source : fMother2Source);}
  ESource GetGrandMotherSource(Int_t branch)                   const {return (branch==1 ? fGrandMother1Source : fGrandMother2Source);}
  Bool_t GetCheckBothChargesLegs(Int_t branch)                 const {return (branch==1 ? fCheckBothChargesLeg1 : fCheckBothChargesLeg2);}
  Bool_t GetCheckBothChargesMothers(Int_t branch)              const {return (branch==1 ? fCheckBothChargesMother1 : fCheckBothChargesMother2);}
  Bool_t GetCheckBothChargesGrandMothers(Int_t branch)         const {return (branch==1 ? fCheckBothChargesGrandMother1 : fCheckBothChargesGrandMother2);}
  EBranchRelation GetMothersRelation()                         const {return fMothersRelation;}
  EBranchRelation GetGrandMothersRelation()                    const {return fGrandMothersRelation;}
  TMCProcess GetGEANTProcess()                                 const {return fGEANTProcess;}
  Bool_t GetCheckGEANTProcess()                                const {return fCheckGEANTProcess;}
  Bool_t GetFillPureMCStep()                                   const {return fFillPureMCStep;}
  Bool_t GetCheckMotherGrandmotherRelation()                   const {return fCheckMotherGrandmother;}
  Bool_t GetCheckMotherGrandmotherDiffPairRelation()           const {return fCheckMotherGrandmotherDiffPair;}
  Bool_t GetMotherIsGrandmother()                              const {return fMotherIsGrandmother;}
  Bool_t GetMotherIsGrandmotherDiffPair()                      const {return fMotherIsGrandmotherDiffPair;}
  Bool_t GetCheckStackForPDG()                                 const {return fCheckStackForPDG;}
  Int_t GetStackPDG()                                          const {return fStackPDG;}
  Bool_t GetCheckUnlikeSign()                                  const {return fCheckUnlikeSign;}
  Bool_t GetCheckLikeSignPP()                                  const {return fCheckLikeSignPP;}
  Bool_t GetCheckLikeSignMM()                                  const {return fCheckLikeSignMM;}
  Bool_t GetCheckCorrelatedHF()                                const {return fCheckCorrelatedHF;}

  static AliDielectronSignalMC* GetJpsiMCsignalDef(EJpsiSignals kSignal);
  static const char* GetJpsiMCsignalDefName(EJpsiSignals kSignal) {return fgkJpsiSignals[kSignal];}

  void SetJpsiRadiative(EJpsiRadiativ rad) { fJpsiRadiative=rad;    }
  EJpsiRadiativ GetJpsiRadiative() const   { return fJpsiRadiative; }
 private:
  // Switches to compare also like sign pairs with the SignalMC definition and to deactivate checking of unlike sign pairs.
  Bool_t fCheckUnlikeSign; // default usage, true by default.
  Bool_t fCheckLikeSignPP; // ++ pairs
  Bool_t fCheckLikeSignMM; // -- pairs

  // PDG codes for legs, mothers and grand-mothers
  Int_t fLeg1;                        // leg 1 PDG
  Int_t fLeg2;                        // leg 2 PDG
  Int_t fMother1;                     // mother 1 PDG
  Int_t fMother2;                     // mother 2 PDG
  Int_t fGrandMother1;                // grandmother 1 PDG
  Int_t fGrandMother2;                // grandmother 2 PDG
  Int_t fStackPDG;                    // PDG to exclude from stack


  // Toggle on/off the use of the PDG codes as inclusion or exclusion
  // Example: if fLeg1=211 and fLeg1Exclude=kTRUE than all codes will be accepted for leg 1 with
  //          the exception of 211 (pions)
  Bool_t fLeg1Exclude;                // leg 1
  Bool_t fLeg2Exclude;                // leg 2
  Bool_t fMother1Exclude;             // mother 1
  Bool_t fMother2Exclude;             // mother 2
  Bool_t fGrandMother1Exclude;        // grandmother 1
  Bool_t fGrandMother2Exclude;        // grandmother 2

  // Particle sources
  ESource fLeg1Source;                // leg 1 source
  ESource fLeg2Source;                // leg 2 source
  ESource fMother1Source;             // mother 1 source
  ESource fMother2Source;             // mother 2 source
  ESource fGrandMother1Source;        // grandmother 1 source
  ESource fGrandMother2Source;        // grandmother 2 source

  // Flaggs whether to check both charges of a given PDG code
  Bool_t fCheckBothChargesLeg1;         // check both charges of the legs pdg
  Bool_t fCheckBothChargesLeg2;         //                leg2
  Bool_t fCheckBothChargesMother1;      //                mother 1
  Bool_t fCheckBothChargesMother2;      //                mother 2
  Bool_t fCheckBothChargesGrandMother1; //              grand mother 1
  Bool_t fCheckBothChargesGrandMother2; //              grand mother 2
  Bool_t fCheckGEANTProcess;            //              GEANT process

  Bool_t fCheckMotherGrandmother;               // check if a mother is also a grandmother to select B -> e D X -> ee X
  Bool_t fMotherIsGrandmother;                  // check if a mother is also a grandmother to select B -> e D X -> ee X
  Bool_t fCheckMotherGrandmotherDiffPair;       // check if a mother is also a grandmother to select Dalitz pairs
  Bool_t fMotherIsGrandmotherDiffPair;          // check if a mother is also a grandmother to select Dalitz pairs

  EBranchRelation fMothersRelation;   // mother 1&2 relation (same, different or whatever)
  EBranchRelation fGrandMothersRelation;   // mother 1&2 relation (same, different or whatever)

  TMCProcess fGEANTProcess;           // GEANT process ID (see roots TMCProcess)
  EJpsiRadiativ fJpsiRadiative;       // check for J/psi radiative decay

  Bool_t fCheckCorrelatedHF;            // check if tracks are closed first mother for correlated charm/beauty in Hijing

  Bool_t fCheckStackForPDG;             // check whole stack to exclude a pdg code to get rid of B feeddown in D sample
  Bool_t fFillPureMCStep;             // check and fill the pure MC step

  const static char *fgkJpsiSignals[EJpsiSignals::kEnd];

  ClassDef(AliDielectronSignalMC,5);
};

#endif
