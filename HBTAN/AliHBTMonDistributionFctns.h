#ifndef ALIHBTMONDISTRIBUTIONFCTNS_H
#define ALIHBTMONDISTRIBUTIONFCTNS_H
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

#include "AliHBTMonitorFunction.h"
/***********************************************************************/
/***********************************************************************/
/*************************************************************************************/ 

class AliHBTMonPxDistributionFctn: public AliHBTMonOneParticleFctn1D
{
 public:
   AliHBTMonPxDistributionFctn(Int_t nbins = 200, Double_t maxXval = 1.4, Double_t minXval = -1.4);
   virtual ~AliHBTMonPxDistributionFctn(){};
   TH1* GetResult(){return fResult;}
 protected:
   Double_t GetValue(AliHBTParticle * particle) { return particle->Px();}
   ClassDef(AliHBTMonPxDistributionFctn,1)
};
/*************************************************************************************/ 

class AliHBTMonPyDistributionFctn: public AliHBTMonOneParticleFctn1D
{
 public:
   AliHBTMonPyDistributionFctn(Int_t nbins = 200, Double_t maxXval = 1.4, Double_t minXval = -1.4);
   virtual ~AliHBTMonPyDistributionFctn(){};
   TH1* GetResult(){return fResult;} 
 protected:
   Double_t GetValue(AliHBTParticle * particle) { return particle->Py();}
   ClassDef(AliHBTMonPyDistributionFctn,1)
};
/*************************************************************************************/ 

class AliHBTMonPzDistributionFctn: public AliHBTMonOneParticleFctn1D
{
 public:
   AliHBTMonPzDistributionFctn(Int_t nbins = 200, Double_t maxXval = 1.4, Double_t minXval = -1.4);
   virtual ~AliHBTMonPzDistributionFctn(){};
   TH1* GetResult(){return fResult;} 
 protected:
   Double_t GetValue(AliHBTParticle * particle) { return particle->Pz();}
   ClassDef(AliHBTMonPzDistributionFctn,1)
 
};
/*************************************************************************************/ 

class AliHBTMonPDistributionFctn: public AliHBTMonOneParticleFctn1D
{
 public:
   AliHBTMonPDistributionFctn(Int_t nbins = 200, Double_t maxXval = 1.4, Double_t minXval = 0.0);
   virtual ~AliHBTMonPDistributionFctn(){};
   TH1* GetResult(){return fResult;} 
 protected:
   Double_t GetValue(AliHBTParticle * particle) { return particle->P();}
   ClassDef(AliHBTMonPDistributionFctn,1)
 
};
/*************************************************************************************/ 

class AliHBTMonPtDistributionFctn: public AliHBTMonOneParticleFctn1D
{
 public:
   AliHBTMonPtDistributionFctn(Int_t nbins = 200, Double_t maxXval = 1.4, Double_t minXval = 0.0);
   virtual ~AliHBTMonPtDistributionFctn(){};
   TH1* GetResult(){return fResult;} 
 protected:
   Double_t GetValue(AliHBTParticle * particle) { return particle->Pt();}
   ClassDef(AliHBTMonPtDistributionFctn,1)
};

/***********************************************************************/
class AliHBTMonPxDistributionVsPtFctn: public AliHBTMonOneParticleFctn2D
 {
  public: 
   AliHBTMonPxDistributionVsPtFctn(Int_t nXbins = 200, Double_t maxXval = 1.4, Double_t minXval = -0.1, 
                             Int_t nYbins = 200, Double_t maxYval = 1.4, Double_t minYval =-1.4);
   virtual ~AliHBTMonPxDistributionVsPtFctn(){}

   void GetValues(AliHBTParticle* partparticle,  Double_t& x, Double_t& y)
    {
      x = partparticle->Pt();
      y = partparticle->Px();
    }
   TH1* GetResult(){return fResult;} 
   ClassDef(AliHBTMonPxDistributionVsPtFctn,1)
 };

/***********************************************************************/
class AliHBTMonPyDistributionVsPtFctn: public AliHBTMonOneParticleFctn2D
 {
  public: 
   AliHBTMonPyDistributionVsPtFctn(Int_t nXbins = 200, Double_t maxXval = 1.4, Double_t minXval = -0.1, 
                             Int_t nYbins = 200, Double_t maxYval = 1.4, Double_t minYval =-1.4);
   virtual ~AliHBTMonPyDistributionVsPtFctn(){}

   void GetValues(AliHBTParticle* partparticle,  Double_t& x, Double_t& y)
    {
     x = partparticle->Pt();
     y = partparticle->Py();
    }
  TH1* GetResult(){return fResult;} 
  ClassDef(AliHBTMonPyDistributionVsPtFctn,1)
 };
/***********************************************************************/
class AliHBTMonPzDistributionVsPtFctn: public AliHBTMonOneParticleFctn2D
 {
  public: 
   AliHBTMonPzDistributionVsPtFctn(Int_t nXbins = 200, Double_t maxXval = 1.4, Double_t minXval = -0.1, 
                             Int_t nYbins = 200, Double_t maxYval = 1.4, Double_t minYval =-1.4);
   virtual ~AliHBTMonPzDistributionVsPtFctn(){}

   void GetValues(AliHBTParticle* partparticle,  Double_t& x, Double_t& y)
    {
     x = partparticle->Pt();
     y = partparticle->Pz();
    }
   TH1* GetResult(){return fResult;} 
   ClassDef(AliHBTMonPzDistributionVsPtFctn,1)
 };

/***********************************************************************/
class AliHBTMonPDistributionVsPtFctn: public AliHBTMonOneParticleFctn2D
 {
  public: 
   AliHBTMonPDistributionVsPtFctn(Int_t nXbins = 200, Double_t maxXval = 1.4, Double_t minXval = -0.1, 
                             Int_t nYbins = 200, Double_t maxYval = 1.4, Double_t minYval =-1.4);
   virtual ~AliHBTMonPDistributionVsPtFctn(){}

   void GetValues(AliHBTParticle* partparticle,  Double_t& x, Double_t& y)
    {
     x = partparticle->Pt();
     y = partparticle->P();
    }
   TH1* GetResult(){return fResult;} 
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
   TH1* GetResult(){return fResult;} 
 protected:
   Double_t GetValue(AliHBTParticle * particle) { return particle->Phi();}
   ClassDef(AliHBTMonPhiDistributionFctn,1)
};

/***********************************************************************/
class AliHBTMonThetaDistributionFctn: public AliHBTMonOneParticleFctn1D
{
 public:
   AliHBTMonThetaDistributionFctn(Int_t nbins = 200, Double_t maxXval = 3.14, Double_t minXval = 0.0);
   virtual ~AliHBTMonThetaDistributionFctn(){};
   TH1* GetResult(){return fResult;} 
 protected:
   Double_t GetValue(AliHBTParticle * particle) { return particle->Theta();}
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

   void GetValues(AliHBTParticle* partparticle,  Double_t& x, Double_t& y)
    {
     x = partparticle->Pt();
     y = partparticle->Phi();
    }
  TH1* GetResult(){return fResult;}
  ClassDef(AliHBTMonPhiDistributionVsPtFctn,1)
 };

/***********************************************************************/
class AliHBTMonThetaDistributionVsPtFctn: public AliHBTMonOneParticleFctn2D
{
  public: 
   AliHBTMonThetaDistributionVsPtFctn(Int_t nXbins = 200, Double_t maxXval = 1.4, Double_t minXval = -0.1, 
                             Int_t nYbins = 200, Double_t maxYval = 3.14, Double_t minYval =0.0);
   virtual ~AliHBTMonThetaDistributionVsPtFctn(){}

   void GetValues(AliHBTParticle* partparticle,  Double_t& x, Double_t& y)
    {
      x = partparticle->Pt();
      y = partparticle->Theta();
    }
   TH1* GetResult(){return fResult;} 
   ClassDef(AliHBTMonThetaDistributionVsPtFctn,1)
 };

/***********************************************************************/
class AliHBTMonVxDistributionFctn: public AliHBTMonOneParticleFctn1D
{
 public:
   AliHBTMonVxDistributionFctn(Int_t nbins = 200, Double_t maxXval = 500, Double_t minXval = -500);
   virtual ~AliHBTMonVxDistributionFctn(){};
   TH1* GetResult(){return fResult;}
 protected:
   Double_t GetValue(AliHBTParticle * particle) { return particle->Vx();}
   ClassDef(AliHBTMonVxDistributionFctn,1)
};
/***********************************************************************/
class AliHBTMonVyDistributionFctn: public AliHBTMonOneParticleFctn1D
{
 public:
   AliHBTMonVyDistributionFctn(Int_t nbins = 200, Double_t maxXval = 500, Double_t minXval = -500);
   virtual ~AliHBTMonVyDistributionFctn(){};
   TH1* GetResult(){return fResult;}
 protected:
   Double_t GetValue(AliHBTParticle * particle) { return particle->Vy();}
   ClassDef(AliHBTMonVyDistributionFctn,1)
};
/***********************************************************************/
class AliHBTMonVzDistributionFctn: public AliHBTMonOneParticleFctn1D
{
 public:
   AliHBTMonVzDistributionFctn(Int_t nbins = 200, Double_t maxXval = 300, Double_t minXval = -300);
   virtual ~AliHBTMonVzDistributionFctn(){};
   TH1* GetResult(){return fResult;}
 protected:
   Double_t GetValue(AliHBTParticle * particle) { return particle->Vz();}
   ClassDef(AliHBTMonVzDistributionFctn,1)
};
/***********************************************************************/
class AliHBTMonRDistributionFctn: public AliHBTMonOneParticleFctn1D
{
 public:
   AliHBTMonRDistributionFctn(Int_t nbins = 200, Double_t maxXval = 500, Double_t minXval = -500);
   virtual ~AliHBTMonRDistributionFctn(){};
   TH1* GetResult(){return fResult;}
 protected:
   Double_t GetValue(AliHBTParticle * p) { return TMath::Sqrt(p->Vx()*p->Vx() + p->Vy()*p->Vy() + p->Vz()*p->Vz());}
   ClassDef(AliHBTMonRDistributionFctn,1)
};

/***********************************************************************/
class AliHBTMonVyDistributionVsVxFctn: public AliHBTMonOneParticleFctn2D
{
  public: 
   AliHBTMonVyDistributionVsVxFctn(Int_t nXbins = 200, Double_t maxXval = 10.0, Double_t minXval = -10.0, 
                                   Int_t nYbins = 200, Double_t maxYval = 10.0, Double_t minYval =-10.0);
   virtual ~AliHBTMonVyDistributionVsVxFctn(){}

   void GetValues(AliHBTParticle* partparticle,  Double_t& x, Double_t& y)
    {
      x = partparticle->Vx();
      y = partparticle->Vy();
    }
   TH1* GetResult(){return fResult;} 
   ClassDef(AliHBTMonVyDistributionVsVxFctn,1)
 };


class AliHBTMonRtDistributionVsVzFctn: public AliHBTMonOneParticleFctn2D
{
  public: 
   AliHBTMonRtDistributionVsVzFctn(Int_t nXbins = 200, Double_t maxXval = 10.0, Double_t minXval = -10.0, 
                                   Int_t nYbins = 100, Double_t maxYval = 10.0, Double_t minYval = 0.0);
   virtual ~AliHBTMonRtDistributionVsVzFctn(){}
 
   void GetValues(AliHBTParticle* partparticle,  Double_t& x, Double_t& y)
    {
      x = partparticle->Vz();
      y = TMath::Hypot(partparticle->Vx(),partparticle->Vy());
    }
   TH1* GetResult(){return fResult;}
   ClassDef(AliHBTMonRtDistributionVsVzFctn,1)
 };

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

#endif
