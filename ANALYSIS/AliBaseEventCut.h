#ifndef ALIBASEEVENTCUT_H
#define ALIBASEEVENTCUT_H
//________________________________
///////////////////////////////////////////////////////////
//
// class AliBaseEventCut
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

class AliBaseEventCut: public TObject
{
  public: 
    AliBaseEventCut();
    AliBaseEventCut(Double_t min,Double_t max);
    virtual ~AliBaseEventCut(){}
    
    virtual Bool_t Pass(AliAOD* aod) const;//returns kTRUE if rejected
  protected:
    virtual Double_t GetValue(AliAOD* aod) const = 0;
    
    Double_t fMin;//Minimum value
    Double_t fMax;//Maximum value
  private:
    ClassDef(AliBaseEventCut,1)
};

/************************************************************/

class AliPrimVertexXCut: public AliBaseEventCut
{
 public: 
   AliPrimVertexXCut(){}
   AliPrimVertexXCut(Double_t min,Double_t max):AliBaseEventCut(min,max){}
   virtual ~AliPrimVertexXCut(){}
 protected:
   Double_t GetValue(AliAOD* aod) const;
   
 private:
   ClassDef(AliPrimVertexXCut,1)
};
/************************************************************/

class AliPrimVertexYCut: public AliBaseEventCut
{
 public: 
   AliPrimVertexYCut(){}
   AliPrimVertexYCut(Double_t min,Double_t max):AliBaseEventCut(min,max){}
   virtual ~AliPrimVertexYCut(){}
   
 protected:
   Double_t GetValue(AliAOD* aod) const;
   
 private:
   ClassDef(AliPrimVertexYCut,1)
};
/************************************************************/

class AliPrimVertexZCut: public AliBaseEventCut
{
 public: 
   AliPrimVertexZCut(){}
   AliPrimVertexZCut(Double_t min,Double_t max):AliBaseEventCut(min,max){}
   virtual ~AliPrimVertexZCut(){}
 protected:
   Double_t GetValue(AliAOD* aod) const;
   
 private:
   ClassDef(AliPrimVertexZCut,1)
};


/************************************************************/

class AliNChargedCut: public AliBaseEventCut
{
 public: 
   AliNChargedCut(){}
   AliNChargedCut(Double_t min, Double_t max, Double_t etamin = -10.0, Double_t etamax = 10.0):
       AliBaseEventCut(min,max),fEtaMin(etamin),fEtaMax(etamax){}
   virtual ~AliNChargedCut(){}
 protected:
   Double_t GetValue(AliAOD* aod) const;
   Double_t fEtaMin;//Defines max of eta range where mult is caclulated
   Double_t fEtaMax;//Defines min of eta range where mult is caclulated
   
 private:
   ClassDef(AliNChargedCut,1)
};


#endif
