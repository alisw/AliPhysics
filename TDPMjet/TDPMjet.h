#ifndef ROOT_TDPMjet
#define ROOT_TDPMjet

//+SEQ,CopyRight,T=NOINCLUDE.

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// 				  TDPMjet                               //
//                                                                      //
// This class implements an interface to the DPMJET 3.0 event generator.//
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TGenerator
//*KEEP,TGenerator.
#include "TGenerator.h"
//*KEND.
#endif
typedef enum {kDpmMb, kDpmMbNonDiffr} DpmProcess_t;

class TDPMjet : public TGenerator {

public:
   
   TDPMjet();
   TDPMjet(DpmProcess_t ip, Int_t Ip, Int_t Ipz, Int_t It, Int_t Itz, Double_t Epn, Double_t CMEn);

   virtual       ~TDPMjet() {;}

   virtual void  Initialize();

   virtual void  GenerateEvent();

   virtual Int_t ImportParticles(TClonesArray *particles, Option_t *option="");
   virtual TObjArray      *ImportParticles(Option_t * /*option*/) {return 0;}
   
   // Parameters for the generation:
   virtual void  SetNEvent(Int_t iev)    {fNEvent = iev;}
   virtual Int_t GetNEvent()             {return fNEvent;}
   
   virtual void  SetfIp(Int_t Ip)        {fIp = Ip;}
   virtual Int_t GetfIp() const	         {return fIp;}

   virtual void  SetfIpz(Int_t Ipz)      {fIpz = Ipz;}
   virtual Int_t GetfIpz() const         {return fIpz;}

   virtual void  SetfIt(Int_t It)        {fIt = It;}
   virtual Int_t GetfIt() const	         {return fIt;}

   virtual void  SetfItz(Int_t Itz)      {fItz = Itz;}
   virtual Int_t GetfItz() const         {return fItz;}

   virtual void  SetfEpn(Double_t Epn)   {fEpn = Epn;}
   virtual Double_t GetfEpn() const      {return fEpn;}

   virtual void  SetfCMEn(Double_t CMEn) {fCMEn = CMEn;}
   virtual Double_t GetfCMEn() const     {return fCMEn;}

   virtual void  SetfIdp(Int_t idp)      {fIdp = idp;}
   virtual Int_t GetfIdp() const         {return fIdp;}

   virtual void SetbRange(Double_t bmin, Double_t bmax) 
   		{fBmin = bmin; fBmax = bmax;} 
   virtual Double_t GetMinImpactParameter() const {return fBmin;}  
   virtual Double_t GetMaxImpactParameter() const {return fBmax;}

   virtual void  SetfFCentr(Int_t icentr)  {fFCentr = icentr;}
   virtual Int_t GetfFCentr() const        {return fFCentr;}

   virtual void  SetPi0Decay(Int_t iPi0)  {fPi0Decay = iPi0;}

   
   // Access to DPMJET common blocks:
   virtual Int_t    GetEvNum() const;	    	  
   virtual Int_t    GetEntriesNum() const;	    	  
   virtual Int_t    GetNumStablePc() const;	    	  
   virtual Float_t  GetTotEnergy() const;
   virtual Int_t    GetStatusCode(Int_t evnum) const; 
   virtual Int_t    GetPDGCode(Int_t evnum) const; 
   virtual Double_t Getpx(Int_t evnum) const;  
   virtual Double_t Getpy(Int_t evnum) const;  
   virtual Double_t Getpz(Int_t evnum) const;  
   virtual Double_t GetEnergy(Int_t evnum) const;  
   virtual Double_t GetMass(Int_t evnum) const;
   
   virtual Int_t    GetFragmentA(Int_t evnum) const;	
   virtual Int_t    GetFragmentZ(Int_t evnum) const;	
   
   virtual Double_t GetXSFrac() const;
   virtual Double_t GetBImpac() const;
   virtual Double_t GetProjRadius() const;
   virtual Double_t GetTargRadius() const;
   virtual Int_t GetProjWounded() const;
   virtual Int_t GetTargWounded() const;
   virtual Int_t GetProjSpectators() const;
   virtual Int_t GetTargSpectators() const;
   
  
   // Access to DPMJET routines:
   virtual void Dt_Dtuini(int nevts, double epn, int npmass, int npchar, 
   			  int ntmass, int ntchar, int idp, int iemu);
	
   virtual void Dt_Kkinc(int npmass, int npchar, int ntmass, int ntchar, 
   			 int idp, double elab, int kkmat, int irej);

   virtual void Pho_Phist(int imode, double weight);

   virtual void Dt_Dtuout();

   virtual void Dt_Rndm(int idummy);   
   virtual void Dt_Rndmst(int na1, int na2, int na3, int nb1);   
   virtual void Dt_Rndmin(int u, int c, int cd, int cm, int i, int j);   
   virtual void Dt_Rndmou(int u, int c, int cd, int cm, int i, int j);   

protected:

   Int_t        fNEvent;  // Event number to be generated 
   Int_t        fIp;	  // Projectile mass
   Int_t        fIpz;	  // Projectile charge
   Int_t        fIt;	  // Target mass
   Int_t        fItz;	  // Target charge
   Float_t      fEpn;	  // Beam energy
   Float_t      fPpn;	  // Beam momentum
   Float_t      fCMEn;	  // Energy in CM 
   Int_t        fIdp;	  // Internal particle code
   Float_t      fBmin;	  // Minimum impact parameter
   Float_t      fBmax;	  // Maximum impact parameter
   Int_t        fFCentr;  // Flag to force central collisions
   Int_t        fPi0Decay;// Flag for pi0 decays
   DpmProcess_t fProcess; // Process type
   
   ClassDef(TDPMjet,2)  //Interface to DPMJET Event Generator
};

#endif







