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

class AliEventBaseCut: public TObject
{
  public: 
    enum EEventCutProperty {
       kPrimVertexXCut,kPrimVertexYCut,kPrimVertexZCut,
       kNChargedCut,kNone
     };

    AliEventBaseCut();
    AliEventBaseCut(Double_t min,Double_t max, EEventCutProperty prop = kNone);
    virtual ~AliEventBaseCut(){}
    virtual Bool_t Rejected(AliAOD* aod) const;//returns kTRUE if rejected
    virtual void   SetRange(Double_t min, Double_t max){fMin = min; fMax = max;}

    virtual EEventCutProperty GetProperty()const{return fProperty;}

  protected:
    virtual Double_t GetValue(AliAOD* aod) const = 0;
    
    Double_t fMin;//Minimum value
    Double_t fMax;//Maximum value
    EEventCutProperty fProperty;//Defines the type of the cut - used by the setters cut
    
  private:
    ClassDef(AliEventBaseCut,1)
};

/************************************************************/

class AliPrimVertexXCut: public AliEventBaseCut
{
 public: 
   AliPrimVertexXCut():AliEventBaseCut(0,0,kPrimVertexXCut){}
   AliPrimVertexXCut(Double_t min,Double_t max):AliEventBaseCut(min,max,kPrimVertexXCut){}
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
   AliPrimVertexYCut():AliEventBaseCut(0,0,kPrimVertexYCut){}
   AliPrimVertexYCut(Double_t min,Double_t max):AliEventBaseCut(min,max,kPrimVertexYCut){}
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
   AliPrimVertexZCut():AliEventBaseCut(0,0,kPrimVertexZCut){}
   AliPrimVertexZCut(Double_t min,Double_t max):AliEventBaseCut(min,max,kPrimVertexZCut){}
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
   AliNChargedCut():AliEventBaseCut(0,0,kNChargedCut),fEtaMin(-10.0),fEtaMax(10.0){}
   AliNChargedCut(Int_t min, Int_t max, Double_t etamin = -10.0, Double_t etamax = 10.0):
       AliEventBaseCut(min,max,kNChargedCut),fEtaMin(etamin),fEtaMax(etamax){}
   virtual ~AliNChargedCut(){}
   
   void     SetEtaRange(Double_t min,Double_t max){fEtaMin = min;fEtaMax = max;}
 protected:
   Double_t GetValue(AliAOD* aod) const;
   Double_t fEtaMin;//Defines max of eta range where mult is caclulated
   Double_t fEtaMax;//Defines min of eta range where mult is caclulated
   
 private:
   ClassDef(AliNChargedCut,1)
};


#endif
