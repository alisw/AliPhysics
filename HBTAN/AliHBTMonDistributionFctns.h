#ifndef ALIHBTMONPXDISTRIBUTIONVSPTFCTN_H
#define ALIHBTMONPXDISTRIBUTIONVSPTFCTN_H

// added by Zbigniew.Chajecki@cern.ch
// this classes create distribution functions of particle momentum 

class AliHBTMonPxDistributionFctn;
class AliHBTMonPxDistributionVsPtFctn;
class AliHBTMonPyDistributionFctn;
class AliHBTMonPyDistributionVsPtFctn;
class AliHBTMonPzDistributionFctn;
class AliHBTMonPzDistributionVsPtFctn;
class AliHBTMonPDistributionFctn;
class AliHBTMonPDistributionVsPtFctn;
class AliHBTMonPtDistributionFctn;

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
  public:
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
  public:
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
  public:
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
  public:
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
  public:
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
  protected:
  private:
  public:
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
  protected:
  private:
  public:
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
  protected:
  private:
  public:
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
  protected:
  private:
  public:
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
  public:
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
  public:
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
  protected:
  private:
  public:
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
  protected:
  private:
  public:
   ClassDef(AliHBTMonThetaDistributionVsPtFctn,1)
 };

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

#endif
