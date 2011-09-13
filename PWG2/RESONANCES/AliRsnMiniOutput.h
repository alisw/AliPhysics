#ifndef ALIRSNMINIOUTPUT_H
#define ALIRSNMINIOUTPUT_H

//
// Mini-Output
// All the definitions needed for building a RSN histogram
// including:
// -- properties of resonance (mass, PDG code if needed)
// -- properties of daughters (assigned mass, charges)
// -- definition of output histogram
// 

#include "AliRsnDaughter.h"
#include "AliRsnMiniParticle.h"

class THnSparse;
class TList;
class TH1;

class TList;
class TClonesArray;
class AliRsnMiniAxis;
class AliRsnMiniPair;
class AliRsnMiniEvent;

typedef AliRsnDaughter::ESpecies RSNPID;

class AliRsnMiniOutput : public TNamed {
public:

   enum EOutputType {
      kHistogram,
      kHistogramSparse,
      kTypes
   };
   
   enum EComputation {
      kEventOnly,
      kTrackPair,
      kTrackPairMix,
      kTrackPairRotated1,
      kTrackPairRotated2,
      kTruePair,
      kMother,
      kComputations
   };
   
   AliRsnMiniOutput();
   AliRsnMiniOutput(const char *name, EOutputType type, EComputation src = kTrackPair);
   AliRsnMiniOutput(const char *name, const char *outType, const char *compType);
   AliRsnMiniOutput(const AliRsnMiniOutput &copy);
   AliRsnMiniOutput& operator=(const AliRsnMiniOutput &copy);
   
   Bool_t          IsEventOnly()        const {return (fComputation == kEventOnly);}
   Bool_t          IsTrackPair()        const {return (fComputation == kTrackPair);}
   Bool_t          IsTrackPairMix()     const {return (fComputation == kTrackPairMix);}
   Bool_t          IsTruePair()         const {return (fComputation == kTruePair);}
   Bool_t          IsMother()           const {return (fComputation == kMother);}
   Bool_t          IsDefined()          const {return (IsEventOnly() || IsTrackPair() || IsTrackPairMix() || IsTruePair() || IsMother());}
   Bool_t          IsLikeSign()         const {return (fCharge[0] == fCharge[1]);}
   Bool_t          IsSameCut()          const {return (fCutID[0] == fCutID[1]);}
   Bool_t          IsSameDaughter()     const {return (fDaughter[0] == fDaughter[1]);}
   //Bool_t          IsSymmetric()        const {return (IsLikeSign() && IsSameCut());}
   Bool_t          IsSymmetric()        const {return (IsLikeSign() && IsSameDaughter());}
                                        
   EOutputType     GetOutputType()      const {return fOutputType;}
   EComputation    GetComputation()     const {return fComputation;}
   Int_t           GetCutID(Int_t i)    const {if (i <= 0) return fCutID [0]; else return fCutID [1];}
   RSNPID          GetDaughter(Int_t i) const {if (i <= 0) return fDaughter[0]; else return fDaughter[1];}
   Double_t        GetMass(Int_t i)     const {return AliRsnDaughter::SpeciesMass(GetDaughter(i));}
   Int_t           GetPDG(Int_t i)      const {return AliRsnDaughter::SpeciesPDG(GetDaughter(i));}
   Int_t           GetCharge(Int_t i)   const {if (i <= 0) return fCharge[0]; else return fCharge[1];}
   Int_t           GetMotherPDG()       const {return fMotherPDG;}
   Double_t        GetMotherMass()      const {return fMotherMass;}
                   
   void            SetOutputType(EOutputType type)    {fOutputType = type;}
   void            SetComputation(EComputation src)   {fComputation = src;}
   void            SetCutID(Int_t i, Int_t   value)   {if (i <= 0) fCutID [0] = value; else fCutID [1] = value;}
   void            SetDaughter(Int_t i, RSNPID value) {if (i <= 0) fDaughter[0] = value; else fDaughter[1] = value;}
   void            SetCharge(Int_t i, Char_t  value)  {if (i <= 0) fCharge[0] = value; else fCharge[1] = value;}
   void            SetMotherPDG(Int_t pdg)            {fMotherPDG = pdg;}
   void            SetMotherMass(Double_t mass)       {fMotherMass = mass;}
   void            SetPairCuts(AliRsnCutSet *set)     {fPairCuts = set;}
                   
   void            AddAxis(Int_t id, Int_t nbins, Double_t min, Double_t max);
   void            AddAxis(Int_t id, Double_t min, Double_t max, Double_t step);
   void            AddAxis(Int_t id, Int_t nbins, Double_t *values);
   AliRsnMiniAxis* GetAxis(Int_t i)  {if (i >= 0 && i < fAxes.GetEntries()) return (AliRsnMiniAxis*)fAxes[i]; return 0x0;}
   Double_t*       GetAllComputed()  {return fComputed.GetArray();}
   
   AliRsnMiniPair& Pair() {return fPair;}
   Bool_t          Init(const char *prefix, TList *list);
   Bool_t          FillMother(const AliRsnMiniPair *pair, AliRsnMiniEvent *event, TClonesArray *valueList);
   Bool_t          FillEvent(AliRsnMiniEvent *event, TClonesArray *valueList);
   Int_t           FillPair(AliRsnMiniEvent *event1, AliRsnMiniEvent *event2, TClonesArray *valueList, Bool_t refFirst = kTRUE);
                  
private:

   void   CreateHistogram(const char *name);
   void   CreateHistogramSparse(const char *name);
   void   ComputeValues(AliRsnMiniEvent *event, TClonesArray *valueList);
   void   FillHistogram();

   EOutputType      fOutputType;       //  type of output
   EComputation     fComputation;      //  type of computation
   Int_t            fCutID[2];         //  ID of cut set used to select tracks
   RSNPID           fDaughter[2];      //  species of daughters
   Char_t           fCharge[2];        //  required track charge
   Int_t            fMotherPDG;        //  PDG code of resonance
   Double_t         fMotherMass;       //  nominal resonance mass
   AliRsnCutSet    *fPairCuts;         //  cuts on the pair
                           
   Int_t            fOutputID;         //  index of output object in container list
   TClonesArray     fAxes;             //  definitions for the axes of each value
   TArrayD          fComputed;         //! temporary container for all computed values
   AliRsnMiniPair   fPair;             //! minipair for computations
   TList           *fList;             //! pointer to the TList containing the output
   TArrayI          fSel1;             //! list of selected particles for definition 1
   TArrayI          fSel2;             //! list of selected particles for definition 2
   
   ClassDef(AliRsnMiniOutput,1)  // AliRsnMiniOutput class
};

#endif
