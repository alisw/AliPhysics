#ifndef ALIEVENTBASECUT_H
#define ALIEVENTBASECUT_H
//________________________________
///////////////////////////////////////////////////////////
//
// class AliEventBaseCut
//
// Base class for cauts that checks only one event property
//
// Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////

#include "TObject.h"

class AliAOD;

enum AliEventCutProperty
 {
   kPrimVertexXCut,
   kPrimVertexYCut,
   kPrimVertexZCut,
   kNChargedCut
 };

class AliEventBaseCut: public TObject
{
  public: 
    AliEventBaseCut();
    AliEventBaseCut(Double_t min,Double_t max);
    virtual ~AliEventBaseCut(){}
    
    virtual Bool_t Pass(AliAOD* aod) const;//returns kTRUE if rejected
  protected:
    virtual Double_t GetValue(AliAOD* aod) const = 0;
    
    Double_t fMin;//Minimum value
    Double_t fMax;//Maximum value
  private:
    ClassDef(AliEventBaseCut,1)
};

/************************************************************/

class AliPrimVertexXCut: public AliEventBaseCut
{
 public: 
   AliPrimVertexXCut(){}
   AliPrimVertexXCut(Double_t min,Double_t max):AliEventBaseCut(min,max){}
   virtual ~AliPrimVertexXCut(){}
 protected:
   Double_t GetValue(AliAOD* aod) const;
   
 private:
   ClassDef(AliPrimVertexXCut,1)
};
/************************************************************/

class AliPrimVertexYCut: public AliEventBaseCut
{
 public: 
   AliPrimVertexYCut(){}
   AliPrimVertexYCut(Double_t min,Double_t max):AliEventBaseCut(min,max){}
   virtual ~AliPrimVertexYCut(){}
   
 protected:
   Double_t GetValue(AliAOD* aod) const;
   
 private:
   ClassDef(AliPrimVertexYCut,1)
};
/************************************************************/

class AliPrimVertexZCut: public AliEventBaseCut
{
 public: 
   AliPrimVertexZCut(){}
   AliPrimVertexZCut(Double_t min,Double_t max):AliEventBaseCut(min,max){}
   virtual ~AliPrimVertexZCut(){}
 protected:
   Double_t GetValue(AliAOD* aod) const;
   
 private:
   ClassDef(AliPrimVertexZCut,1)
};


/************************************************************/

class AliNChargedCut: public AliEventBaseCut
{
 public: 
   AliNChargedCut(){}
   AliNChargedCut(Double_t min, Double_t max, Double_t etamin = -10.0, Double_t etamax = 10.0):
       AliEventBaseCut(min,max),fEtaMin(etamin),fEtaMax(etamax){}
   virtual ~AliNChargedCut(){}
 protected:
   Double_t GetValue(AliAOD* aod) const;
   Double_t fEtaMin;//Defines max of eta range where mult is caclulated
   Double_t fEtaMax;//Defines min of eta range where mult is caclulated
   
 private:
   ClassDef(AliNChargedCut,1)
};


#endif
