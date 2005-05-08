#ifndef ALIJFMIXEVENTH
#define ALIJFMIXEVENTH

#include <TMath.h>

class TTree;
class TBranch;
class TBranch;
class TFile;
class TClonesArray;

class AliJFMixEvent 
{
 public:
  AliJFMixEvent(Char_t *file1, Char_t *file2=0);
  AliJFMixEvent(Char_t *files1, Char_t *tname1, Char_t *files2, Char_t *tname2);
  virtual ~AliJFMixEvent();

  void Init(Char_t *file1, Char_t *file2=0);
  void InitChains(Char_t *files1, Char_t *tname1, Char_t *files2=0, Char_t *tname2=0);
  void Clean();

  Int_t CreateNextMixedEvent();
  Int_t CreateMixedEvent(Int_t i=0,Int_t j=0);

  TClonesArray* const GetMixedParticles() {return fMixedParticles;}
  TClonesArray* const GetParticles () {return fParticles1;}
  TClonesArray* const GetParticles2() {return fParticles2;}

  Int_t const GetCurNEvent1() const {return fEvent1;}
  Int_t const GetCurNEvent2() const {return fEvent2;}
  Int_t const GetMaxNEvent1() const {return fMaxEvent1;}
  Int_t const GetMaxNEvent2() const {return fMaxEvent2;}
  Int_t const GetStatus()     const {return fStatus;}

  void Debug();
  void SetMarkPythia(Bool_t b){fMarkPythia=b;}

 protected:
  Int_t MixEvent();

  Int_t fStatus; //positive if initialized

  Int_t fEvent1;
  Int_t fEvent2;
  Int_t fMaxEvent1;
  Int_t fMaxEvent2;

  Bool_t fMarkPythia; //true if pythia particles 
                      //will be marked

  TTree *fTree1; //!
  TTree *fTree2; //!

  TFile *fFile1; //!
  TFile *fFile2; //!

  TBranch *fBranch1; //!
  TBranch *fBranch2; //!

  TClonesArray *fParticles1; //!
  TClonesArray *fParticles2; //!
  TClonesArray *fMixedParticles; //!

  ClassDef(AliJFMixEvent,1) //AliJFMixEvent class
};
 
#endif /*ALIJFMIXEVENTH*/
