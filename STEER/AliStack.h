#ifndef ALI_STACK_H
#define ALI_STACK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TArrayI.h>
#include <TVirtualMCStack.h>

class AliHeader;
class TClonesArray;
class TFile;
class TObjArray;
class TParticle;
class TTree;

class AliStack : public TVirtualMCStack
{
  public:
    // creators, destructors
    AliStack(Int_t size);
    AliStack();
    AliStack(const AliStack& st);
    virtual ~AliStack();
    AliStack& operator=(const AliStack& st)
      {st.Copy(*this); return(*this);}

    // methods
    virtual void  SetTrack(Int_t done, Int_t parent, Int_t pdg, 
  	              Float_t *pmom, Float_t *vpos, Float_t *polar, 
                      Float_t tof, TMCProcess mech, Int_t &ntr,
                      Float_t weight, Int_t is);
    virtual void  SetTrack(Int_t done, Int_t parent, Int_t pdg,
  	              Double_t px, Double_t py, Double_t pz, Double_t e,
  		      Double_t vx, Double_t vy, Double_t vz, Double_t tof,
		      Double_t polx, Double_t poly, Double_t polz,
		      TMCProcess mech, Int_t &ntr, Double_t weight,
		      Int_t is);
    virtual TParticle* GetNextTrack(Int_t& track);
    virtual TParticle* GetCurrentTrack() {return fCurrentTrack;}
    virtual TParticle* GetPrimaryForTracking(Int_t i);    

    void  MakeTree(Int_t event, const char *file);
    void  BeginEvent(Int_t event);
    void  FinishRun();
    Bool_t GetEvent(Int_t nevent);
    void  PurifyKine();
    void  FinishEvent();
    void  FlagTrack(Int_t track);
    void  KeepTrack(Int_t itrack); 
    void  Reset(Int_t size = 0);
    void  DumpPart(Int_t i) const;
    void  DumpPStack ();
    void  DumpLoadedStack () const;

    // set methods
    void  SetNtrack(Int_t ntrack);
    virtual void  SetCurrentTrack(Int_t track);                           
    void  SetHighWaterMark(Int_t hgwmk);    
    // get methods
    virtual Int_t GetNtrack() const;
    Int_t       GetNprimary() const;
    virtual Int_t CurrentTrack() const;
    virtual Int_t CurrentTrackParent() const;
    TObjArray*  Particles() const;
    TParticle*  Particle(Int_t id);
    Int_t       GetPrimary(Int_t id);
    TTree*      TreeK() const {return fTreeK;}
    TParticle*  ParticleFromTreeK(Int_t id) const;
    Int_t       TreeKEntry(Int_t id) const;
    
  protected:
    // methods
    void  CleanParents();
    void  ResetArrays(Int_t size);
    TParticle* GetNextParticle();
    Bool_t KeepPhysics(TParticle* part);
    
  private:
    void Copy(AliStack &st) const;

    // data members
    TClonesArray  *fParticles;         //! Pointer to list of particles
    TObjArray     *fParticleMap;       //! Map of particles in the supporting TClonesArray
    TArrayI        fParticleFileMap;   //  Map for particle ids 
    TParticle     *fParticleBuffer;    //! Pointer to current particle for writing
    TParticle     *fCurrentTrack;      //! Pointer to particle currently transported
    TTree         *fTreeK;             //! Particle stack  
    Int_t          fNtrack;            //  Number of tracks
    Int_t          fNprimary;          //  Number of primaries
    Int_t          fCurrent;           //! Last track returned from the stack
    Int_t          fCurrentPrimary;    //! Last primary track returned from the stack
    Int_t          fHgwmk;             //! Last track purified
    Int_t          fLoadPoint;         //! Next free position in the particle buffer
    
    ClassDef(AliStack,3) //Particles stack
};

// inline

inline void  AliStack::SetNtrack(Int_t ntrack)
{ fNtrack = ntrack; }

inline void  AliStack::SetCurrentTrack(Int_t track)
{ fCurrent = track; }

inline Int_t AliStack::GetNtrack() const
{ return fNtrack; }

inline Int_t AliStack::GetNprimary() const
{ return fNprimary; }

inline Int_t AliStack::CurrentTrack() const 
{ return fCurrent; }

inline TObjArray* AliStack::Particles() const
{ return fParticleMap; }

#endif //ALI_STACK_H
