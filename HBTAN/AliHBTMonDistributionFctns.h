#ifndef ALIHBTMONDISTRIBUTIONFCTNS_H
#define ALIHBTMONDISTRIBUTIONFCTNS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//______________________________________________________________
////////////////////////////////////////////////////////////////
//
// class AliHBTMonPxDistributionFctn;
// class AliHBTMonPxDistributionVsPtFctn;
// class AliHBTMonPyDistributionFctn;
// class AliHBTMonPyDistributionVsPtFctn;
// class AliHBTMonPzDistributionFctn;
// class AliHBTMonPzDistributionVsPtFctn;
// class AliHBTMonPDistributionFctn;
// class AliHBTMonPDistributionVsPtFctn;
// class AliHBTMonPtDistributionFctn;
// class AliHBTMonVxDistributionFctn;
// class AliHBTMonVyDistributionFctn;
// class AliHBTMonVzDistributionFctn;
// class AliHBTMonRDistributionFctn;
// class AliHBTMonVyDistributionVsVxFctn;
// class AliHBTMonRtDistributionVsVzFctn;
//
// added by Zbigniew.Chajecki@cern.ch
// this classes create distribution functions of particle momentum
//
/////////////////////////////////////////////////////////////////

class AliHBTMonPxDistributionFctn;
class AliHBTMonPxDistributionVsPtFctn;
class AliHBTMonPyDistributionFctn;
class AliHBTMonPyDistributionVsPtFctn;
class AliHBTMonPzDistributionFctn;
class AliHBTMonPzDistributionVsPtFctn;
class AliHBTMonPDistributionFctn;
class AliHBTMonPDistributionVsPtFctn;
class AliHBTMonPtDistributionFctn;

class AliHBTMonVxDistributionFctn;
class AliHBTMonVyDistributionFctn;
class AliHBTMonVzDistributionFctn;
class AliHBTMonRDistributionFctn;

class AliHBTMonVyDistributionVsVxFctn;
class AliHBTMonRtDistributionVsVzFctn;

#include <TMath.h>

#include "AliHBTMonitorFunction.h"
/***********************************************************************/
/***********************************************************************/
/*************************************************************************************/ 

class AliHBTMonPxDistributionFctn: public AliHBTMonOneParticleFctn1D
{
 public:
   AliHBTMonPxDistributionFctn(Int_t nbins = 200, Double_t maxXval = 1.4, Double_t minXval = -1.4);
   virtual ~AliHBTMonPxDistributionFctn(){};
 protected:
   Double_t GetValue(AliVAODParticle * particle) const;
   ClassDef(AliHBTMonPxDistributionFctn,1)
};
/*************************************************************************************/ 

class AliHBTMonPyDistributionFctn: public AliHBTMonOneParticleFctn1D
{
 public:
   AliHBTMonPyDistributionFctn(Int_t nbins = 200, Double_t maxXval = 1.4, Double_t minXval = -1.4);
   virtual ~AliHBTMonPyDistributionFctn(){};
 protected:
   Double_t GetValue(AliVAODParticle * particle) const { return particle->Py();}
   ClassDef(AliHBTMonPyDistributionFctn,1)
};
/*************************************************************************************/ 

class AliHBTMonPzDistributionFctn: public AliHBTMonOneParticleFctn1D
{
 public:
   AliHBTMonPzDistributionFctn(Int_t nbins = 200, Double_t maxXval = 1.4, Double_t minXval = -1.4);
   virtual ~AliHBTMonPzDistributionFctn(){};
 protected:
   Double_t GetValue(AliVAODParticle * particle) const { return particle->Pz();}
   ClassDef(AliHBTMonPzDistributionFctn,1)
 
};
/*************************************************************************************/ 

class AliHBTMonPDistributionFctn: public AliHBTMonOneParticleFctn1D
{
 public:
   AliHBTMonPDistributionFctn(Int_t nbins = 200, Double_t maxXval = 1.4, Double_t minXval = 0.0);
   virtual ~AliHBTMonPDistributionFctn(){};
 protected:
   Double_t GetValue(AliVAODParticle * particle) const { return particle->P();}
   ClassDef(AliHBTMonPDistributionFctn,1)
 
};
/*************************************************************************************/ 

class AliHBTMonPtDistributionFctn: public AliHBTMonOneParticleFctn1D
{
 public:
   AliHBTMonPtDistributionFctn(Int_t nbins = 200, Double_t maxXval = 1.4, Double_t minXval = 0.0);
   virtual ~AliHBTMonPtDistributionFctn(){};
 protected:
   Double_t GetValue(AliVAODParticle * particle) const { return particle->Pt();}
   ClassDef(AliHBTMonPtDistributionFctn,1)
};

/***********************************************************************/
class AliHBTMonPxDistributionVsPtFctn: public AliHBTMonOneParticleFctn2D
 {
  public: 
   AliHBTMonPxDistributionVsPtFctn(Int_t nXbins = 200, Double_t maxXval = 1.4, Double_t minXval = -0.1, 
                             Int_t nYbins = 200, Double_t maxYval = 1.4, Double_t minYval =-1.4);
   virtual ~AliHBTMonPxDistributionVsPtFctn(){}

   void GetValues(AliVAODParticle* partparticle,  Double_t& x, Double_t& y) const
    {
      x = partparticle->Pt();
      y = partparticle->Px();
    }
   ClassDef(AliHBTMonPxDistributionVsPtFctn,1)
 };

/***********************************************************************/
class AliHBTMonPyDistributionVsPtFctn: public AliHBTMonOneParticleFctn2D
 {
  public: 
   AliHBTMonPyDistributionVsPtFctn(Int_t nXbins = 200, Double_t maxXval = 1.4, Double_t minXval = -0.1, 
                             Int_t nYbins = 200, Double_t maxYval = 1.4, Double_t minYval =-1.4);
   virtual ~AliHBTMonPyDistributionVsPtFctn(){}

   void GetValues(AliVAODParticle* partparticle,  Double_t& x, Double_t& y) const
    {
     x = partparticle->Pt();
     y = partparticle->Py();
    }
  ClassDef(AliHBTMonPyDistributionVsPtFctn,1)
 };
/***********************************************************************/
class AliHBTMonPzDistributionVsPtFctn: public AliHBTMonOneParticleFctn2D
 {
  public: 
   AliHBTMonPzDistributionVsPtFctn(Int_t nXbins = 200, Double_t maxXval = 1.4, Double_t minXval = -0.1, 
                             Int_t nYbins = 200, Double_t maxYval = 1.4, Double_t minYval =-1.4);
   virtual ~AliHBTMonPzDistributionVsPtFctn(){}

   void GetValues(AliVAODParticle* partparticle,  Double_t& x, Double_t& y) const
    {
     x = partparticle->Pt();
     y = partparticle->Pz();
    }
   ClassDef(AliHBTMonPzDistributionVsPtFctn,1)
 };

/***********************************************************************/
class AliHBTMonPDistributionVsPtFctn: public AliHBTMonOneParticleFctn2D
 {
  public: 
   AliHBTMonPDistributionVsPtFctn(Int_t nXbins = 200, Double_t maxXval = 1.4, Double_t minXval = -0.1, 
                             Int_t nYbins = 200, Double_t maxYval = 1.4, Double_t minYval =-1.4);
   virtual ~AliHBTMonPDistributionVsPtFctn(){}

   void GetValues(AliVAODParticle* partparticle,  Double_t& x, Double_t& y) const
    {
     x = partparticle->Pt();
     y = partparticle->P();
    }
   ClassDef(AliHBTMonPDistributionVsPtFctn,1)
 };

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

class AliHBTMonPhiDistributionFctn: public AliHBTMonOneParticleFctn1D
{
 public:
   AliHBTMonPhiDistributionFctn(Int_t nbins = 200, Double_t maxXval = 3.14, Double_t minXval = 0.0);
   virtual ~AliHBTMonPhiDistributionFctn(){};
 protected:
   Double_t GetValue(AliVAODParticle * particle) const { return particle->Phi();}
   ClassDef(AliHBTMonPhiDistributionFctn,1)
};

/***********************************************************************/
class AliHBTMonThetaDistributionFctn: public AliHBTMonOneParticleFctn1D
{
 public:
   AliHBTMonThetaDistributionFctn(Int_t nbins = 200, Double_t maxXval = 3.14, Double_t minXval = 0.0);
   virtual ~AliHBTMonThetaDistributionFctn(){};
 protected:
   Double_t GetValue(AliVAODParticle * particle) const { return particle->Theta();}
   ClassDef(AliHBTMonThetaDistributionFctn,1)
 
};
/***********************************************************************/
/***********************************************************************/
class AliHBTMonPhiDistributionVsPtFctn: public AliHBTMonOneParticleFctn2D
 {
  public: 
   AliHBTMonPhiDistributionVsPtFctn(Int_t nXbins = 200, Double_t maxXval = 1.4, Double_t minXval = -0.1, 
                             Int_t nYbins = 200, Double_t maxYval = 3.14, Double_t minYval =0.0);
   virtual ~AliHBTMonPhiDistributionVsPtFctn(){}

   void GetValues(AliVAODParticle* partparticle,  Double_t& x, Double_t& y) const
    {
     x = partparticle->Pt();
     y = partparticle->Phi();
    }
  ClassDef(AliHBTMonPhiDistributionVsPtFctn,1)
 };

/***********************************************************************/
class AliHBTMonThetaDistributionVsPtFctn: public AliHBTMonOneParticleFctn2D
{
  public: 
   AliHBTMonThetaDistributionVsPtFctn(Int_t nXbins = 200, Double_t maxXval = 1.4, Double_t minXval = -0.1, 
                             Int_t nYbins = 200, Double_t maxYval = 3.14, Double_t minYval =0.0);
   virtual ~AliHBTMonThetaDistributionVsPtFctn(){}

   void GetValues(AliVAODParticle* partparticle,  Double_t& x, Double_t& y) const
    {
      x = partparticle->Pt();
      y = partparticle->Theta();
    }
   ClassDef(AliHBTMonThetaDistributionVsPtFctn,1)
 };

/***********************************************************************/
class AliHBTMonVxDistributionFctn: public AliHBTMonOneParticleFctn1D
{
 public:
   AliHBTMonVxDistributionFctn(Int_t nbins = 200, Double_t maxXval = 500, Double_t minXval = -500);
   virtual ~AliHBTMonVxDistributionFctn(){};
 protected:
   Double_t GetValue(AliVAODParticle * particle)  const{ return particle->Vx();}
   ClassDef(AliHBTMonVxDistributionFctn,1)
};
/***********************************************************************/
class AliHBTMonVyDistributionFctn: public AliHBTMonOneParticleFctn1D
{
 public:
   AliHBTMonVyDistributionFctn(Int_t nbins = 200, Double_t maxXval = 500, Double_t minXval = -500);
   virtual ~AliHBTMonVyDistributionFctn(){};
 protected:
   Double_t GetValue(AliVAODParticle * particle) const { return particle->Vy();}
   ClassDef(AliHBTMonVyDistributionFctn,1)
};
/***********************************************************************/
class AliHBTMonVzDistributionFctn: public AliHBTMonOneParticleFctn1D
{
 public:
   AliHBTMonVzDistributionFctn(Int_t nbins = 200, Double_t maxXval = 300, Double_t minXval = -300);
   virtual ~AliHBTMonVzDistributionFctn(){};
 protected:
   Double_t GetValue(AliVAODParticle * particle) const { return particle->Vz();}
   ClassDef(AliHBTMonVzDistributionFctn,1)
};
/***********************************************************************/
class AliHBTMonRDistributionFctn: public AliHBTMonOneParticleFctn1D
{
 public:
   AliHBTMonRDistributionFctn(Int_t nbins = 200, Double_t maxXval = 500, Double_t minXval = -500);
   virtual ~AliHBTMonRDistributionFctn(){};
 protected:
   Double_t GetValue(AliVAODParticle * p) const { return TMath::Sqrt(p->Vx()*p->Vx() + p->Vy()*p->Vy() + p->Vz()*p->Vz());}
   ClassDef(AliHBTMonRDistributionFctn,1)
};

/***********************************************************************/
class AliHBTMonVyDistributionVsVxFctn: public AliHBTMonOneParticleFctn2D
{
  public: 
   AliHBTMonVyDistributionVsVxFctn(Int_t nXbins = 200, Double_t maxXval = 10.0, Double_t minXval = -10.0, 
                                   Int_t nYbins = 200, Double_t maxYval = 10.0, Double_t minYval =-10.0);
   virtual ~AliHBTMonVyDistributionVsVxFctn(){}

   void GetValues(AliVAODParticle* partparticle,  Double_t& x, Double_t& y) const
    {
      x = partparticle->Vx();
      y = partparticle->Vy();
    }
   ClassDef(AliHBTMonVyDistributionVsVxFctn,1)
 };


class AliHBTMonRtDistributionVsVzFctn: public AliHBTMonOneParticleFctn2D
{
  public: 
   AliHBTMonRtDistributionVsVzFctn(Int_t nXbins = 200, Double_t maxXval = 10.0, Double_t minXval = -10.0, 
                                   Int_t nYbins = 100, Double_t maxYval = 10.0, Double_t minYval = 0.0);
   virtual ~AliHBTMonRtDistributionVsVzFctn(){}
 
   void GetValues(AliVAODParticle* partparticle,  Double_t& x, Double_t& y) const
    {
      x = partparticle->Vz();
      y = TMath::Hypot(partparticle->Vx(),partparticle->Vy());
    }
   ClassDef(AliHBTMonRtDistributionVsVzFctn,1)
 };

/***********************************************************************/
/***********************************************************************/

#endif
