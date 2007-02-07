#ifndef ALIHBTMONITORFUNCTION_H
#define ALIHBTMONITORFUNCTION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//__________________________________________________________________
////////////////////////////////////////////////////////////////////
//
// class AliHBTMonitorFunction
//
// class AliHBTMonOneParticleFctn
// class AliHBTMonTwoParticleFctn
//
// class AliHBTMonOneParticleFctn1D
// class AliHBTMonOneParticleFctn2D
// class AliHBTMonOneParticleFctn3D
//
// class AliHBTMonTwoParticleFctn1D
// class AliHBTMonTwoParticleFctn2D
// class AliHBTMonTwoParticleFctn3D
//
// Base Classes for monitoring functions
// author: chajecki@if.pw.edu.pl
//
/******************************************************************/
/*
Base classes for monitor functions

          monitor function
               /    \
              /      \
             /        \
            /          \
           /            \
          /              \
         /                \
   one particle     two particle  
     /  |  \            /  |  \
    /   |   \          /   |   \
   1D  2D   3D        1D  2D   3D

Zbigniew.Chajecki@cern.ch

*/
///////////////////////////////////////////////////////////////////////

#include "AliAODParticleCut.h"

#include <TMath.h>
#include <TH2.h>
#include <TH3.h>

class AliVAODParticle;

class AliHBTMonitorFunction: public TNamed
//Abstract base class for HBT functions
{
  public:
    AliHBTMonitorFunction();
    AliHBTMonitorFunction(const char* name,const char* title);
    AliHBTMonitorFunction(const AliHBTMonitorFunction& /*in*/);
    virtual ~AliHBTMonitorFunction();
    
    AliHBTMonitorFunction& operator=(const AliHBTMonitorFunction& /*in*/);
    
    
    virtual TH1* GetResult() = 0;

    Int_t Write(const char* /*x1*/ = "",Int_t /*x2*/ = 0, Int_t /*x3*/ = 0);
    Int_t Write(const char* x1 = "",Int_t x2 = 0, Int_t x3 = 0) const {return TObject::Write(x1,x2,x3);}
    virtual void Init();
    virtual const char* Name(){return GetName();}
    void Rename(const Char_t * name); 
    void Rename(const Char_t * name, const Char_t * title); 
    
    void SetParticleCut(AliAODParticleCut* cut);

    virtual AliVAODParticle* CheckParticle(AliVAODParticle* particle) const;

  protected:
    AliAODParticleCut*      fParticleCut;//Particle cut
    
  private:  
   ClassDef(AliHBTMonitorFunction,1)
};
/******************************************************************/
/******************************************************************/
inline AliVAODParticle* AliHBTMonitorFunction::CheckParticle(AliVAODParticle* particle) const
{
  //check if particle meets the cut criteria
  if(fParticleCut->Rejected(particle)) //if the particle is BAD
   { 
     return 0x0;//it is BAD as well - so return
   }
  return particle; 
}

/******************************************************************/
/******************************************************************/
/******************************************************************/

class AliHBTMonOneParticleFctn: public AliHBTMonitorFunction
{
  public:
    AliHBTMonOneParticleFctn(){}
    AliHBTMonOneParticleFctn(const Char_t *name, const Char_t *title):AliHBTMonitorFunction(name,title){}
    AliHBTMonOneParticleFctn(const AliHBTMonOneParticleFctn& in):AliHBTMonitorFunction(in){MayNotUse("Cpy Ctor");}
    virtual ~AliHBTMonOneParticleFctn(){}
    
    AliHBTMonOneParticleFctn& operator=(const AliHBTMonOneParticleFctn& /*in*/){MayNotUse("operator=");return *this;} 
    
    virtual void Process(AliVAODParticle* particle) = 0;
    
  protected:
  private:  
   ClassDef(AliHBTMonOneParticleFctn,1)
  
};
/******************************************************************/
class AliHBTMonOneParticleFctn1D: public AliHBTMonOneParticleFctn
{
 public:
  AliHBTMonOneParticleFctn1D();
  AliHBTMonOneParticleFctn1D(Int_t nbins, Double_t maxXval, Double_t minXval);
  AliHBTMonOneParticleFctn1D(const Char_t *name, const Char_t *title,
                      Int_t nbins = 100, Double_t maxXval = 1.4, Double_t minXval = 0.0);
  AliHBTMonOneParticleFctn1D(const AliHBTMonOneParticleFctn1D& in):
               AliHBTMonOneParticleFctn(in),fResult(0x0){MayNotUse("Cpy Ctor");}
  
  virtual ~AliHBTMonOneParticleFctn1D();
  
  AliHBTMonOneParticleFctn1D& operator=(const AliHBTMonOneParticleFctn1D& /*in*/){MayNotUse("operator=");return *this;}   
  TH1* GetResult(){return this->fResult;}

  void Process(AliVAODParticle* particle);

 protected:
  virtual Double_t GetValue(AliVAODParticle* particle) const = 0; 
  TH1D* fResult;//histogram to be filled
 private:
  ClassDef(AliHBTMonOneParticleFctn1D,2)
};
/******************************************************************/
 
class AliHBTMonOneParticleFctn2D: public AliHBTMonOneParticleFctn
{
 public:
  AliHBTMonOneParticleFctn2D(Int_t nXbins = 200, Double_t maxXval = 1.5, Double_t minXval = 0.0, 
                      Int_t nYbins = 200, Double_t maxYval = 1.5, Double_t minYval =-0.1);
  AliHBTMonOneParticleFctn2D(const AliHBTMonOneParticleFctn2D& in):
               AliHBTMonOneParticleFctn(in),fResult(0x0){MayNotUse("Cpy Ctor");}
  virtual ~AliHBTMonOneParticleFctn2D();
  
  AliHBTMonOneParticleFctn2D& operator=(const AliHBTMonOneParticleFctn2D& /*in*/){MayNotUse("operator=");return *this;}   
  TH1* GetResult(){return this->fResult;}
  
  void Process(AliVAODParticle* particle);

 protected:
  virtual void GetValues(AliVAODParticle* particle, Double_t&, Double_t&) const = 0;

  TH2D* fResult;//histogram to be filled
  
 private:
  ClassDef(AliHBTMonOneParticleFctn2D,1)
};
/******************************************************************/
/******************************************************************/
/******************************************************************/

class AliHBTMonOneParticleFctn3D: public AliHBTMonOneParticleFctn
{
 public:
  AliHBTMonOneParticleFctn3D(Int_t nXbins = 200, Double_t maxXval = 1.5, Double_t minXval = 0.0, 
                      Int_t nYbins = 200, Double_t maxYval = .15, Double_t minYval =-0.15, 
                      Int_t nZbins = 200, Double_t maxZval = .15, Double_t minZval =-0.15);
  AliHBTMonOneParticleFctn3D(const AliHBTMonOneParticleFctn3D& in):
               AliHBTMonOneParticleFctn(in),fResult(0x0){MayNotUse("Cpy Ctor");}
	    
  virtual ~AliHBTMonOneParticleFctn3D();

  AliHBTMonOneParticleFctn3D& operator=(const AliHBTMonOneParticleFctn3D& /*in*/){MayNotUse("operator=");return *this;}   
  TH1* GetResult(){return this->fResult;}

 protected:
  TH3D* fResult;//histogram to be filled

 private:
  ClassDef(AliHBTMonOneParticleFctn3D,1)
};
/******************************************************************/
/******************************************************************/
class AliHBTMonTwoParticleFctn: public AliHBTMonitorFunction
{
  public:
    AliHBTMonTwoParticleFctn(){};
    AliHBTMonTwoParticleFctn(const Char_t *name, const Char_t *title):AliHBTMonitorFunction(name,title){}
    AliHBTMonTwoParticleFctn(const AliHBTMonTwoParticleFctn& in):AliHBTMonitorFunction(in){MayNotUse("Cpy Ctor");}
    virtual ~AliHBTMonTwoParticleFctn(){};
    AliHBTMonTwoParticleFctn& operator=(const AliHBTMonTwoParticleFctn& /*in*/){MayNotUse("operator=");return *this;} 
    
    virtual void 
    Process(AliVAODParticle* trackparticle, AliVAODParticle* partparticle) = 0;
	     
  protected:
  private:  
   ClassDef(AliHBTMonTwoParticleFctn,1)
  
};
/******************************************************************/

class AliHBTMonTwoParticleFctn1D: public AliHBTMonTwoParticleFctn
{
 public:
  AliHBTMonTwoParticleFctn1D(Int_t nbins = 200, Double_t maxval = 1.5, Double_t minval = 0.0);
  AliHBTMonTwoParticleFctn1D(const char* name,const char* title,
                      Int_t nbins = 200, Double_t maxval = 1.5, Double_t minval = 0.0);
  AliHBTMonTwoParticleFctn1D(const AliHBTMonTwoParticleFctn1D& in):
               AliHBTMonTwoParticleFctn(in),fResult(0x0){MayNotUse("Cpy Ctor");}
  virtual ~AliHBTMonTwoParticleFctn1D();
  
  AliHBTMonTwoParticleFctn1D& operator=(const AliHBTMonTwoParticleFctn1D& /*in*/){MayNotUse("operator=");return *this;}   
  TH1* GetResult(){return this->fResult;}
  
  void Process(AliVAODParticle* trackparticle, AliVAODParticle* partparticle);
  
 protected:
  virtual Double_t GetValue(AliVAODParticle* trackparticle, AliVAODParticle* partparticle) const = 0;

  TH1D* fResult;//histogram to be filled

 private:
  ClassDef(AliHBTMonTwoParticleFctn1D,1)
};
/******************************************************************/
class AliHBTMonTwoParticleFctn2D: public AliHBTMonTwoParticleFctn
{
 public:
  AliHBTMonTwoParticleFctn2D(Int_t nXbins = 200, Double_t maxXval = 1.5, Double_t minXval = 0.0, 
                       Int_t nYbins = 200, Double_t maxYval = .15, Double_t minYval =-0.15);
  AliHBTMonTwoParticleFctn2D(const AliHBTMonTwoParticleFctn2D& in):
               AliHBTMonTwoParticleFctn(in),fResult(0x0){MayNotUse("Cpy Ctor");}
  virtual ~AliHBTMonTwoParticleFctn2D();
  
  AliHBTMonTwoParticleFctn2D& operator=(const AliHBTMonTwoParticleFctn2D& /*in*/){MayNotUse("operator=");return *this;}   
  TH1* GetResult(){return this->fResult;}
  
  void Process(AliVAODParticle* trackparticle, AliVAODParticle* partparticle);
  
 protected:
  virtual void GetValues(AliVAODParticle*,AliVAODParticle*, Double_t&, Double_t&) const = 0;

  TH2D* fResult;//histogram to be filled
  
 private:
  ClassDef(AliHBTMonTwoParticleFctn2D,1)
};


/******************************************************************/
class AliHBTMonTwoParticleFctn3D: public AliHBTMonTwoParticleFctn
{
 public:
  AliHBTMonTwoParticleFctn3D(Int_t /*nXbins = 200*/, Double_t /*maxXval = 1.5*/, Double_t /*minXval = 0.0*/, 
			     Int_t /*nYbins = 200*/, Double_t /*maxYval = .15*/, Double_t /*minYval =-0.15*/, 
			     Int_t /*nZbins = 200*/, Double_t /*maxZval = .15*/, Double_t /*minZval =-0.15*/){}
  AliHBTMonTwoParticleFctn3D(const AliHBTMonTwoParticleFctn3D& in):
               AliHBTMonTwoParticleFctn(in),fResult(0x0){MayNotUse("Cpy Ctor");}
  virtual ~AliHBTMonTwoParticleFctn3D(){}
  
  AliHBTMonTwoParticleFctn3D& operator=(const AliHBTMonTwoParticleFctn3D& /*in*/){MayNotUse("operator=");return *this;}   
  TH1* GetResult(){return this->fResult;}
  
  void Process(AliVAODParticle* trackparticle, AliVAODParticle* partparticle);

 protected:
  virtual void GetValues(AliVAODParticle*,AliVAODParticle*, Double_t&, Double_t&,Double_t&) const = 0;

  TH3D* fResult; //histogram to be filled
  
 private:
  ClassDef(AliHBTMonTwoParticleFctn3D,1)
};

/******************************************************************/
/******************************************************************/
/******************************************************************/
/******************************************************************/

#endif
