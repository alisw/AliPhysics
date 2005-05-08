// $Id$

#ifndef ALIJFJETH
#define ALIJFJETH

#ifndef ROOT_TObject 
#include <TObject.h>
#endif

#include <TMath.h>
#include <TParticle.h>
#include <TClonesArray.h>

class AliJFJet: public TObject
{
 public:
  AliJFJet(Int_t n=250);
  virtual ~AliJFJet();

  inline Double_t GetPhi()   {if(!fIsUpdated) Update(); return fPhi;}  
  inline Double_t GetEta()   {if(!fIsUpdated) Update(); return fEta;}  
  inline Double_t GetY()     {if(!fIsUpdated) Update(); return fY;}  
  inline Double_t GetPt()    {if(!fIsUpdated) Update(); return fPt;}  
  inline Double_t GetPx()    {if(!fIsUpdated) Update(); return fPx;}  
  inline Double_t GetPy()    {if(!fIsUpdated) Update(); return fPy;}  
  inline Double_t GetPz()    {if(!fIsUpdated) Update(); return fPz;}
  inline Double_t GetE()     {if(!fIsUpdated) Update(); return fE;}  
  inline Double_t GetE_()    {if(!fIsUpdated) Update(); return fE_;}  
  inline Double_t GetPtSum() {if(!fIsUpdated) Update(); return fPtSum;}  
  inline Double_t GetPhiSum(){if(!fIsUpdated) Update(); return fPhiSum;}  
  inline Double_t GetEtaSum(){if(!fIsUpdated) Update(); return fEtaSum;}  

  inline Double_t GetPhiC()   {if(!fIsUpdated) Update(); return fPhiC;}  
  inline Double_t GetEtaC()   {if(!fIsUpdated) Update(); return fEtaC;}  
  inline Double_t GetYC()     {if(!fIsUpdated) Update(); return fYC;}  
  inline Double_t GetPtC()    {if(!fIsUpdated) Update(); return fPtC;}  
  inline Double_t GetPxC()    {if(!fIsUpdated) Update(); return fPxC;}  
  inline Double_t GetPyC()    {if(!fIsUpdated) Update(); return fPyC;}  
  inline Double_t GetPzC()    {if(!fIsUpdated) Update(); return fPzC;}
  inline Double_t GetEC()     {if(!fIsUpdated) Update(); return fEC;}  
  inline Double_t GetE_C()    {if(!fIsUpdated) Update(); return fE_C;}  
  inline Double_t GetPtSumC() {if(!fIsUpdated) Update(); return fPtSumC;}  
  inline Double_t GetPhiSumC(){if(!fIsUpdated) Update(); return fPhiSumC;}  
  inline Double_t GetEtaSumC(){if(!fIsUpdated) Update(); return fEtaSumC;}  

  inline Double_t GetPhiN()   {if(!fIsUpdated) Update(); return fPhiN;}  
  inline Double_t GetEtaN()   {if(!fIsUpdated) Update(); return fEtaN;}  
  inline Double_t GetYN()     {if(!fIsUpdated) Update(); return fYN;}  
  inline Double_t GetPtN()    {if(!fIsUpdated) Update(); return fPtN;}  
  inline Double_t GetPxN()    {if(!fIsUpdated) Update(); return fPxN;}  
  inline Double_t GetPyN()    {if(!fIsUpdated) Update(); return fPyN;}  
  inline Double_t GetPzN()    {if(!fIsUpdated) Update(); return fPzN;}
  inline Double_t GetEN()     {if(!fIsUpdated) Update(); return fEN;}  
  inline Double_t GetE_N()    {if(!fIsUpdated) Update(); return fE_N;}  
  inline Double_t GetPtSumN() {if(!fIsUpdated) Update(); return fPtSumN;}  
  inline Double_t GetPhiSumN(){if(!fIsUpdated) Update(); return fPhiSumN;}  
  inline Double_t GetEtaSumN(){if(!fIsUpdated) Update(); return fEtaSumN;}  

  inline Double_t GetPhiEM()   {if(!fIsUpdated) Update(); return fPhiEM;}  
  inline Double_t GetEtaEM()   {if(!fIsUpdated) Update(); return fEtaEM;}  
  inline Double_t GetYEM()     {if(!fIsUpdated) Update(); return fYEM;}  
  inline Double_t GetPtEM()    {if(!fIsUpdated) Update(); return fPtEM;}  
  inline Double_t GetPxEM()    {if(!fIsUpdated) Update(); return fPxEM;}  
  inline Double_t GetPyEM()    {if(!fIsUpdated) Update(); return fPyEM;}  
  inline Double_t GetPzEM()    {if(!fIsUpdated) Update(); return fPzEM;}
  inline Double_t GetEEM()     {if(!fIsUpdated) Update(); return fEEM;}  
  inline Double_t GetE_EM()    {if(!fIsUpdated) Update(); return fE_EM;}  
  inline Double_t GetPtSumEM() {if(!fIsUpdated) Update(); return fPtSumEM;}  
  inline Double_t GetPhiSumEM(){if(!fIsUpdated) Update(); return fPhiSumEM;}  
  inline Double_t GetEtaSumEM(){if(!fIsUpdated) Update(); return fEtaSumEM;}  

  inline Int_t GetNCharged() {if(!fIsUpdated) Update(); return fNCharged;}
  inline Int_t GetNNeutral() {if(!fIsUpdated) Update(); return fNNeutral;}
  inline Int_t GetNEM()      {if(!fIsUpdated) Update(); return fNEM;}

  inline  const Int_t GetNPart() const {return fN;}
  inline  const Int_t GetNJet () const {return fNJet;} 

  TParticle* const GetMaxParticle()  {return &fMaxParticle;}
  TClonesArray* const GetParticles() {return &fParticles;}

  inline  void SetNJet(Int_t n) {fNJet=n;}
  virtual void Update(){};

  virtual void Debug();
  virtual void Clean();

  void Print(Option_t *) const {cout << "Jet " << (int)fNJet << ": " << fPtSum << " " << (int)fN << endl;}
  ULong_t Hash() const {return fNJet;}
  Bool_t IsEqual(const TObject *obj) const {return fNJet == ((AliJFJet*)obj)->GetNJet();}
  Bool_t IsSortable() const {return kTRUE;}
  Int_t  Compare(const TObject *obj) const;

 protected:
  Int_t fNJet;
  Int_t fN;
  Int_t fNCharged;
  Int_t fNNeutral;
  Int_t fNEM;

  Double_t fPhi;
  Double_t fEta;
  Double_t fY;
  Double_t fPt;
  Double_t fPx;
  Double_t fPy;
  Double_t fPz;
  Double_t fE;  
  Double_t fE_; //energy without mass
  Double_t fPtSum;
  Double_t fPhiSum;
  Double_t fEtaSum;

  Double_t fPhiC;
  Double_t fEtaC;
  Double_t fYC;
  Double_t fPtC;
  Double_t fPxC;
  Double_t fPyC;
  Double_t fPzC;
  Double_t fEC;  
  Double_t fE_C; //energy without mass
  Double_t fPtSumC;
  Double_t fPhiSumC;
  Double_t fEtaSumC;

  Double_t fPhiN;
  Double_t fEtaN;
  Double_t fYN;
  Double_t fPtN;
  Double_t fPxN;
  Double_t fPyN;
  Double_t fPzN;
  Double_t fEN;  
  Double_t fE_N; //energy without mass
  Double_t fPtSumN;
  Double_t fPhiSumN;
  Double_t fEtaSumN;

  Double_t fPhiEM;
  Double_t fEtaEM;
  Double_t fYEM;
  Double_t fPtEM;
  Double_t fPxEM;
  Double_t fPyEM;
  Double_t fPzEM;
  Double_t fEEM;  
  Double_t fE_EM; //energy without mass
  Double_t fPtSumEM;
  Double_t fPhiSumEM;
  Double_t fEtaSumEM;

  TParticle fMaxParticle;
  TParticle fMaxParticleC;
  TParticle fMaxParticleN;
  TParticle fMaxParticleEM;
  TClonesArray fParticles;

  Bool_t fIsUpdated;

  ClassDef(AliJFJet,1) //AliJFJet class
};

#endif /*ALIJFJETH*/
