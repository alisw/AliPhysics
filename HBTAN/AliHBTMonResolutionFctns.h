#ifndef ALIHBTMONRESOLUTIONFCTNS_H
#define ALIHBTMONRESOLUTIONFCTNS_H
//_______________________________________________________________________________
/////////////////////////////////////////////////////////////////////////////////
//
// class AliHBTMonPxResolutionFctn;
// class AliHBTMonPyResolutionFctn;
// class AliHBTMonPzResolutionFctn;
// class AliHBTMonPResolutionFctn;
// class AliHBTMonPtResolutionFctn;
// class AliHBTMonPhiResolutionFctn;
// class AliHBTMonThetaResolutionFctn;
// class AliHBTMonPxResolutionVsPtFctn;
// class AliHBTMonPyResolutionVsPtFctn;
// class AliHBTMonPzResolutionVsPtFctn;
// class AliHBTMonPResolutionVsPtFctn;
// class AliHBTMonPtResolutionVsPtFctn;
// class AliHBTMonPhiResolutionVsPtFctn;
// class AliHBTMonThetaResolutionVsPtFctn;
//
// Caution: On 2D plots on X axis in simulated values
// That is contrary to two-particle resolutions where it is reconstructed one
//
// added by Zbigniew.Chajecki@cern.ch
// this classes create resolution functions of particle momentum 
//
//////////////////////////////////////////////////////////////////////////////////

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
   ClassDef(AliHBTMonThetaResolutionFctn,1)
 };
/***********************************************************************/
/***********************************************************************/
class AliHBTMonPhiResolutionVsPtFctn: public AliHBTMonTwoParticleFctn2D
 {
  public: 
   AliHBTMonPhiResolutionVsPtFctn(Int_t nXbins = 200, Double_t maxXval = 1.4, Double_t minXval = 0.0, 
                             Int_t nYbins = 200, Double_t maxYval = 0.05, Double_t minYval =-0.05);
   virtual ~AliHBTMonPhiResolutionVsPtFctn(){}
   TH1* GetResult(){return fResult;}
   void GetValues(AliHBTParticle* trackparticle, AliHBTParticle* partparticle, Double_t& x, Double_t& y)
    {
     x = partparticle->Pt();
     y = partparticle->Phi()-trackparticle->Phi();
    }
   ClassDef(AliHBTMonPhiResolutionVsPtFctn,1)
 };
/***********************************************************************/
class AliHBTMonPhiResolutionVsPhiFctn: public AliHBTMonTwoParticleFctn2D
 {
  public: 
   AliHBTMonPhiResolutionVsPhiFctn(Int_t nXbins = 200, Double_t maxXval = TMath::TwoPi(), Double_t minXval = 0.0,
                             Int_t nYbins = 200, Double_t maxYval = 0.05, Double_t minYval =-0.05);
   virtual ~AliHBTMonPhiResolutionVsPhiFctn(){}
   TH1* GetResult(){return fResult;}
   void GetValues(AliHBTParticle* trackparticle, AliHBTParticle* partparticle, Double_t& x, Double_t& y)
    {
     x = partparticle->Phi();
     y = partparticle->Phi()-trackparticle->Phi();
    }
   ClassDef(AliHBTMonPhiResolutionVsPhiFctn,1)
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
   ClassDef(AliHBTMonThetaResolutionVsPtFctn,1)
 };

/***********************************************************************/
class AliHBTMonThetaResolutionVsThetaFctn: public AliHBTMonTwoParticleFctn2D
 {
  public: 
   AliHBTMonThetaResolutionVsThetaFctn(Int_t nXbins = 200, Double_t maxXval = TMath::PiOver2(), Double_t minXval = -TMath::PiOver2(),
                             Int_t nYbins = 200, Double_t maxYval = 0.05, Double_t minYval =-0.05);
   virtual ~AliHBTMonThetaResolutionVsThetaFctn(){}
   TH1* GetResult(){return fResult;}
   void GetValues(AliHBTParticle* trackparticle, AliHBTParticle* partparticle, Double_t& x, Double_t& y)
    {
     y = partparticle->Theta()-trackparticle->Theta();
     x = partparticle->Theta();
    }
   ClassDef(AliHBTMonThetaResolutionVsThetaFctn,1)
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
/***********************************************************************/

#endif
