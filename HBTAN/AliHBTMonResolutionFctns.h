#ifndef ALIMONPXRESOLUTIONVSPTFCTN_H
#define ALIMONPXRESOLUTIONVSPTFCTN_H

// added by Zbigniew.Chajecki@cern.ch
// this classes create resolution functions of particle momentum 

class AliHBTMonPxResolutionFctn;
class AliHBTMonPyResolutionFctn;
class AliHBTMonPzResolutionFctn;
class AliHBTMonPResolutionFctn;
class AliHBTMonPtResolutionFctn;
class AliHBTMonPhiResolutionFctn;
class AliHBTMonThetaResolutionFctn;

class AliHBTMonPxResolutionVsPtFctn;
class AliHBTMonPyResolutionVsPtFctn;
class AliHBTMonPzResolutionVsPtFctn;
class AliHBTMonPResolutionVsPtFctn;
class AliHBTMonPtResolutionVsPtFctn;
class AliHBTMonPhiResolutionVsPtFctn;
class AliHBTMonThetaResolutionVsPtFctn;


#include "AliHBTMonitorFunction.h"
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
class AliHBTMonPxResolutionFctn: public AliHBTMonTwoParticleFctn1D
 {
  public: 
   AliHBTMonPxResolutionFctn(Int_t nbins = 200, Double_t maxXval = 0.05, Double_t minXval = -0.05);
   virtual ~AliHBTMonPxResolutionFctn(){}

   Double_t GetValue(AliHBTParticle * trackparticle,AliHBTParticle * partparticle) 
     { 
        return (partparticle->Px()-trackparticle->Px()) ;
     } 
   TH1* GetResult(){return fResult;} 
  protected:
  private:
  public:
    ClassDef(AliHBTMonPxResolutionFctn,1)
 };
/***********************************************************************/
class AliHBTMonPyResolutionFctn: public AliHBTMonTwoParticleFctn1D
 {
  public: 
   AliHBTMonPyResolutionFctn(Int_t nbins = 200, Double_t maxXval = 0.05, Double_t minXval = -0.05);
   virtual ~AliHBTMonPyResolutionFctn(){}

   Double_t GetValue(AliHBTParticle * trackparticle,AliHBTParticle * partparticle) 
     { 
        return (partparticle->Py()-trackparticle->Py()) ;
     } 
   TH1* GetResult(){return fResult;} 
  protected:
  private:
  public:
    ClassDef(AliHBTMonPyResolutionFctn,1)
 };
/***********************************************************************/
class AliHBTMonPzResolutionFctn: public AliHBTMonTwoParticleFctn1D
 {
  public: 
   AliHBTMonPzResolutionFctn(Int_t nbins = 200, Double_t maxXval = 0.05, Double_t minXval = -0.05);
   virtual ~AliHBTMonPzResolutionFctn(){}

   Double_t GetValue(AliHBTParticle * trackparticle,AliHBTParticle * partparticle) 
     { 
        return (partparticle->Pz()-trackparticle->Pz()) ;
     } 
   TH1* GetResult(){return fResult;} 
  protected:
  private:
  public:
    ClassDef(AliHBTMonPzResolutionFctn,1)
 };
/***********************************************************************/
class AliHBTMonPResolutionFctn: public AliHBTMonTwoParticleFctn1D
 {
  public: 
   AliHBTMonPResolutionFctn(Int_t nbins = 200, Double_t maxXval = 0.05, Double_t minXval = -0.05);
   virtual ~AliHBTMonPResolutionFctn(){}

   Double_t GetValue(AliHBTParticle * trackparticle,AliHBTParticle * partparticle) 
     { 
        return (partparticle->P()-trackparticle->P()) ;
     } 
   TH1* GetResult(){return fResult;} 
  protected:
  private:
  public:
    ClassDef(AliHBTMonPResolutionFctn,1)
 };
/***********************************************************************/
class AliHBTMonPtResolutionFctn: public AliHBTMonTwoParticleFctn1D
 {
  public: 
   AliHBTMonPtResolutionFctn(Int_t nbins = 200, Double_t maxXval = 0.05, Double_t minXval = -0.05);
   virtual ~AliHBTMonPtResolutionFctn(){}

   Double_t GetValue(AliHBTParticle * trackparticle,AliHBTParticle * partparticle) 
     { 
        return (partparticle->Pt()-trackparticle->Pt()) ;
     } 
   TH1* GetResult(){return fResult;} 
  protected:
  private:
  public:
    ClassDef(AliHBTMonPtResolutionFctn,1)
 };
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
class AliHBTMonPxResolutionVsPtFctn: public AliHBTMonTwoParticleFctn2D
 {
  public: 
   AliHBTMonPxResolutionVsPtFctn(Int_t nXbins = 200, Double_t maxXval = 1.4, Double_t minXval = -0.1, 
                             Int_t nYbins = 200, Double_t maxYval = 0.05, Double_t minYval =-0.05);
   virtual ~AliHBTMonPxResolutionVsPtFctn(){}
   TH1* GetResult(){return fResult;}
   void GetValues(AliHBTParticle* trackparticle, AliHBTParticle* partparticle, Double_t& x, Double_t& y)
    {
     x = partparticle->Pt();
     y = partparticle->Px()-trackparticle->Px();
    }
  protected:
  private:
  public:
    ClassDef(AliHBTMonPxResolutionVsPtFctn,1)
 };
/***********************************************************************/
class AliHBTMonPyResolutionVsPtFctn: public AliHBTMonTwoParticleFctn2D
 {
  public: 
   AliHBTMonPyResolutionVsPtFctn(Int_t nXbins = 200, Double_t maxXval = 1.4, Double_t minXval = -0.1, 
                             Int_t nYbins = 200, Double_t maxYval = 0.05, Double_t minYval =-0.05);
   virtual ~AliHBTMonPyResolutionVsPtFctn(){}
   TH1* GetResult(){return fResult;}
   void GetValues(AliHBTParticle* trackparticle, AliHBTParticle* partparticle, Double_t& x, Double_t& y)
    {
     x = partparticle->Pt();
     y = partparticle->Py()-trackparticle->Py();
    }
  protected:
  private:
  public:
    ClassDef(AliHBTMonPyResolutionVsPtFctn,1)
 };
/***********************************************************************/
class AliHBTMonPzResolutionVsPtFctn: public AliHBTMonTwoParticleFctn2D
 {
  public: 
   AliHBTMonPzResolutionVsPtFctn(Int_t nXbins = 200, Double_t maxXval = 1.4, Double_t minXval = -0.1, 
                             Int_t nYbins = 200, Double_t maxYval = 0.05, Double_t minYval =-0.05);
   virtual ~AliHBTMonPzResolutionVsPtFctn(){}
   TH1* GetResult(){return fResult;}
   void GetValues(AliHBTParticle* trackparticle, AliHBTParticle* partparticle, Double_t& x, Double_t& y)
    {
     x = partparticle->Pt();
     y = partparticle->Pz()-trackparticle->Pz();
    }
  protected:
  private:
  public:
    ClassDef(AliHBTMonPzResolutionVsPtFctn,1)
 };
/***********************************************************************/
class AliHBTMonPResolutionVsPtFctn: public AliHBTMonTwoParticleFctn2D
 {
  public: 
   AliHBTMonPResolutionVsPtFctn(Int_t nXbins = 200, Double_t maxXval = 1.4, Double_t minXval = -0.1, 
                             Int_t nYbins = 200, Double_t maxYval = 0.05, Double_t minYval =-0.05);
   virtual ~AliHBTMonPResolutionVsPtFctn(){}
   TH1* GetResult(){return fResult;}
   void GetValues(AliHBTParticle* trackparticle, AliHBTParticle* partparticle, Double_t& x, Double_t& y)
    {
     x = partparticle->Pt();
     y = partparticle->P()-trackparticle->P();
    }
  protected:
  private:
  public:
    ClassDef(AliHBTMonPResolutionVsPtFctn,1)
 };
/***********************************************************************/
class AliHBTMonPtResolutionVsPtFctn: public AliHBTMonTwoParticleFctn2D
 {
  public: 
   AliHBTMonPtResolutionVsPtFctn(Int_t nXbins = 200, Double_t maxXval = 1.4, Double_t minXval = -0.1, 
                             Int_t nYbins = 200, Double_t maxYval = 0.05, Double_t minYval =-0.05);
   virtual ~AliHBTMonPtResolutionVsPtFctn(){}
   TH1* GetResult(){return fResult;}
   void GetValues(AliHBTParticle* trackparticle, AliHBTParticle* partparticle, Double_t& x, Double_t& y)
    {
     x = partparticle->Pt();
     y = partparticle->Pt()-trackparticle->Pt();
    }
  protected:
  private:
  public:
    ClassDef(AliHBTMonPtResolutionVsPtFctn,1)
 };
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
class AliHBTMonPhiResolutionFctn: public AliHBTMonTwoParticleFctn1D
 {
  public: 
   AliHBTMonPhiResolutionFctn(Int_t nbins = 200, Double_t maxXval = 0.05, Double_t minXval = -0.05);
   virtual ~AliHBTMonPhiResolutionFctn(){}

   Double_t GetValue(AliHBTParticle * trackparticle,AliHBTParticle * partparticle) 
     { 
        return (partparticle->Phi()-trackparticle->Phi()) ;
     } 
   TH1* GetResult(){return fResult;} 
  protected:
  private:
  public:
    ClassDef(AliHBTMonPhiResolutionFctn,1)
 };
/***********************************************************************/
class AliHBTMonThetaResolutionFctn: public AliHBTMonTwoParticleFctn1D
 {
  public: 
   AliHBTMonThetaResolutionFctn(Int_t nbins = 200, Double_t maxXval = 0.05, Double_t minXval = -0.05);
   virtual ~AliHBTMonThetaResolutionFctn(){}

   Double_t GetValue(AliHBTParticle * trackparticle,AliHBTParticle * partparticle) 
     { 
        return (partparticle->Theta()-trackparticle->Theta()) ;
     } 
   TH1* GetResult(){return fResult;} 
  protected:
  private:
  public:
    ClassDef(AliHBTMonThetaResolutionFctn,1)
 };
/***********************************************************************/
/***********************************************************************/
class AliHBTMonPhiResolutionVsPtFctn: public AliHBTMonTwoParticleFctn2D
 {
  public: 
   AliHBTMonPhiResolutionVsPtFctn(Int_t nXbins = 200, Double_t maxXval = 1.4, Double_t minXval = -0.1, 
                             Int_t nYbins = 200, Double_t maxYval = 0.05, Double_t minYval =-0.05);
   virtual ~AliHBTMonPhiResolutionVsPtFctn(){}
   TH1* GetResult(){return fResult;}
   void GetValues(AliHBTParticle* trackparticle, AliHBTParticle* partparticle, Double_t& x, Double_t& y)
    {
     x = partparticle->Pt();
     y = partparticle->Phi()-trackparticle->Phi();
    }
  protected:
  private:
  public:
    ClassDef(AliHBTMonPhiResolutionVsPtFctn,1)
 };
/***********************************************************************/
class AliHBTMonThetaResolutionVsPtFctn: public AliHBTMonTwoParticleFctn2D
 {
  public: 
   AliHBTMonThetaResolutionVsPtFctn(Int_t nXbins = 200, Double_t maxXval = 1.4, Double_t minXval = -0.1, 
                             Int_t nYbins = 200, Double_t maxYval = 0.05, Double_t minYval =-0.05);
   virtual ~AliHBTMonThetaResolutionVsPtFctn(){}
   TH1* GetResult(){return fResult;}
   void GetValues(AliHBTParticle* trackparticle, AliHBTParticle* partparticle, Double_t& x, Double_t& y)
    {
     x = partparticle->Pt();
     y = partparticle->Theta()-trackparticle->Theta();
    }
  protected:
  private:
  public:
    ClassDef(AliHBTMonThetaResolutionVsPtFctn,1)
 };
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

#endif
