#ifndef ALIHBTPIDPURITYFCTNS_H
#define ALIHBTPIDPURITYFCTNS_H
//_______________________________________________________________________________
/////////////////////////////////////////////////////////////////////////////////
//
// class AliHBTMonPIDPurityVsPtFctn;
// class AliHBTMonThetaResolutionVsPtFctn;
//
// file: AliHBTPIDPurityFctns.cxx AliHBTPIDPurityFctns.h
//
// Caution: On 2D plots on X axis in simulated values
// That is contrary to two-particle resolutions where it is reconstructed one
//
// added by Piotr.Skowronski@cern.ch
// 
//
//////////////////////////////////////////////////////////////////////////////////

#include "AliHBTFunction.h"
#include "AliHBTMonitorFunction.h"

class AliHBTMonPIDPurityVsPtFctn: public AliHBTMonTwoParticleFctn1D, public AliHBTCorrelFunction
{
  public: 
    AliHBTMonPIDPurityVsPtFctn(Int_t nbins = 20, Double_t maxXval = 2.0, Double_t minXval = 0.0);
    virtual ~AliHBTMonPIDPurityVsPtFctn();
    void Init();
    void Write();
    void Rename(const Char_t * name);
    void Rename(const Char_t * name, const Char_t * title);
    TH1* GetResult();
    Double_t GetValue(AliHBTParticle * /*track*/,AliHBTParticle * /*part*/) { return 0.0; }
    void Process(AliHBTParticle * track,AliHBTParticle * part);
  protected:
    TH1D* fGood;
    ClassDef(AliHBTMonPIDPurityVsPtFctn,1)
};
/***********************************************************************/

class AliHBTMonPIDContaminationVsPtFctn: public AliHBTMonTwoParticleFctn1D, public AliHBTCorrelFunction
{
  public: 
    AliHBTMonPIDContaminationVsPtFctn(Int_t nbins = 20, Double_t maxXval = 2.0, Double_t minXval = 0.0);
    virtual ~AliHBTMonPIDContaminationVsPtFctn();
    void Init();
    void Write();
    void Rename(const Char_t * name);
    void Rename(const Char_t * name, const Char_t * title);
    TH1* GetResult();
    Double_t GetValue(AliHBTParticle * /*track*/,AliHBTParticle * /*part*/) { return 0.0; }
    void Process(AliHBTParticle * track,AliHBTParticle * part);
  protected:
    TH1D* fWrong;  
    ClassDef(AliHBTMonPIDContaminationVsPtFctn,1)
};
/*************************************************************************************/

class AliHBTPairPIDProbVsQInvFctn: public AliHBTOnePairFctn1D, public AliHBTCorrelFunction
{
//Q Invaraint Correlation Function
//1D two particle function
 public:
   AliHBTPairPIDProbVsQInvFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0);
   virtual ~AliHBTPairPIDProbVsQInvFctn(){};
   void ProcessSameEventParticles(AliHBTPair* pair);
   void ProcessDiffEventParticles(AliHBTPair* pair);
   TH1* GetResult();
 protected:
   Double_t GetValue(AliHBTPair * pair){return pair->GetQInv();}
 private:
   ClassDef(AliHBTPairPIDProbVsQInvFctn,1)
};
/*************************************************************************************/
class AliHBTPairPIDProbVsQOutSQideQLongFctn: public AliHBTOnePairFctn3D, public AliHBTCorrelFunction
{

  public:
    AliHBTPairPIDProbVsQOutSQideQLongFctn(Int_t nXbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0,
                             Int_t nYbins = 100, Double_t maxYval = 0.15, Double_t minYval = 0.0,
                 Int_t nZbins = 100, Double_t maxZval = 0.15, Double_t minZval = 0.0);
    virtual  ~AliHBTPairPIDProbVsQOutSQideQLongFctn(){}

    TH1* GetResult();
    void ProcessSameEventParticles(AliHBTPair* part);
    void ProcessDiffEventParticles(AliHBTPair* pair);

  protected:
    void GetValues(AliHBTPair* /*pair*/, Double_t& /*x*/, Double_t& /*y*/, Double_t& /*z*/){}

    ClassDef(AliHBTPairPIDProbVsQOutSQideQLongFctn,1)
};


#endif
