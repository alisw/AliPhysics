#ifndef ALIHBTPARTICLECUT_H
#define ALIHBTPARTICLECUT_H
//__________________________________________________________________________
////////////////////////////////////////////////////////////////////////////
//                                                                        //
// class AliHBTParticleCut                                                //
//                                                                        //
// Classes for single particle cuts                                       //
// User should use only AliHBTParticleCut, eventually                     //
// EmptyCut which passes all particles                                    //
// There is all interface for setting cuts on all particle properties     //
// The main method is Pass - which returns                                //
//         True to reject particle                                        //
//         False in case it meets all the criteria of the given cut       //
//                                                                        //
// User should create (and also destroy) cuts himself                     // 
// and then pass them to the Analysis And Function by a proper method     //
//                                                                        //
//                                                                        //
// more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html   //
// resonsible: Piotr Skowronski@cern.ch                                   //
//                                                                        //
////////////////////////////////////////////////////////////////////////////


#include <TObject.h>
#include "AliHBTParticle.h"


class AliHBTEmptyParticleCut;
class AliHBTParticleCut;
class AliHBTPairCut;
class AliHBTPair;
class AliHbtBaseCut;


/******************************************************************/
/******************************************************************/
/******************************************************************/

enum AliHBTCutProperty
 {
//codes particle property
  kHbtP,  //Momentum
  kHbtPt, //Transverse momentum
  kHbtE,  //Energy
  kHbtRapidity, //
  kHbtPseudoRapidity,
  kHbtPx, //X coAnddinate of the momentum
  kHbtPy, //Y coAnddinate of the momentum
  kHbtPz, //Z coAnddinate of the momentum
  kHbtPhi,//angle
  kHbtTheta,//angle
  kHbtVx,  // vertex X coAnddinate
  kHbtVy,  // vertex Y coAnddinate
  kHbtVz,  // vertex Z coAnddinate
  kHbtPid, // vertex Z coAnddinate
//_____________________________
  kHbtNone
 };

/******************************************************************/
/******************************************************************/
/******************************************************************/

class AliHBTParticleCut: public TObject
{
//Class describing cut on pairs of particles
  public:
    AliHBTParticleCut();
    AliHBTParticleCut(const AliHBTParticleCut& in);
    virtual ~AliHBTParticleCut();
    AliHBTParticleCut& operator = (const AliHBTParticleCut& in);
    
    virtual Bool_t Pass(AliHBTParticle* p) const;
    Bool_t IsEmpty() const {return kFALSE;}
    
    void AddBasePartCut(AliHbtBaseCut* basecut);
    
    Int_t GetPID() const { return fPID;}
    void SetPID(Int_t pid){fPID=pid;}
    void SetMomentumRange(Double_t min, Double_t max);
    void SetPRange(Double_t min, Double_t max){SetMomentumRange(min,max);}
    void SetPtRange(Double_t min, Double_t max);
    void SetEnergyRange(Double_t min, Double_t max);
    void SetRapidityRange(Double_t min, Double_t max);
    void SetYRange(Double_t min, Double_t max){SetRapidityRange(min,max);}
    void SetPseudoRapidityRange(Double_t min, Double_t max);
    void SetPxRange(Double_t min, Double_t max);
    void SetPyRange(Double_t min, Double_t max);
    void SetPzRange(Double_t min, Double_t max);
    void SetPhiRange(Double_t min, Double_t max);
    void SetThetaRange(Double_t min, Double_t max);
    void SetVxRange(Double_t min, Double_t max);
    void SetVyRange(Double_t min, Double_t max);
    void SetVzRange(Double_t min, Double_t max);
    
    void Print(void) const;
  protected:
     
    AliHbtBaseCut* FindCut(AliHBTCutProperty property);

    AliHbtBaseCut ** fCuts;//! Array with cuts
    Int_t fNCuts; //number of base cuts stored in fCuts

    Int_t fPID; //particle PID  - if=0 (rootino) all pids are accepted
          
  private:
    static const Int_t fgkMaxCuts; //Size of the fCuts array

    ClassDef(AliHBTParticleCut,1)
};
/******************************************************************/
/******************************************************************/
/******************************************************************/

class AliHBTEmptyParticleCut:  public AliHBTParticleCut
{
//Empty - it passes possitively all particles - it means returns always False
//Class describing cut on pairs of particles
  public:
    AliHBTEmptyParticleCut(){};
    virtual ~AliHBTEmptyParticleCut(){};
    
    Bool_t Pass(AliHBTParticle*) const {return kFALSE;} //accept everything <<CAN NOT BE const!!!!>>
    Bool_t IsEmpty() const {return kTRUE;}

    ClassDef(AliHBTEmptyParticleCut,1)
 
};

/******************************************************************/
/******************************************************************/
/******************************************************************/

class AliHbtBaseCut: public TObject
 {
   //This class defines the range of some property - pure virtual
   //Property is coded by AliHBTCutTypes type
   
   public:
     
     AliHbtBaseCut(Double_t min = 0.0, Double_t max = 0.0,AliHBTCutProperty prop = kHbtNone):
                   fProperty(prop),fMin(min),fMax(max){}

     virtual           ~AliHbtBaseCut(){}
     
     virtual Bool_t    Pass(AliHBTParticle *p) const;
     
     void              SetRange(Double_t min, Double_t max){fMin = min; fMax = max;}

     void              SetMinimum(Double_t min){fMin = min;}
     void              SetMaximum(Double_t max){fMax = max;}
     
     Double_t          GetMinimum() const {return fMin;}
     Double_t          GetMaximum() const {return fMax;}
     
     AliHBTCutProperty GetProperty() const {return fProperty;}
     virtual void Print(void) const;
     
   protected:
     virtual Double_t  GetValue(AliHBTParticle *) const = 0;

     AliHBTCutProperty fProperty; //property that this cut describes
     Double_t fMin;//minimum value
     Double_t fMax;//maximum value
     
   private:
     void PrintProperty(void) const;
     ClassDef(AliHbtBaseCut,1)
   
 };

inline Bool_t
AliHbtBaseCut::Pass(AliHBTParticle *p) const
{
  //cjecks if particle property fits in range
  if ( (GetValue(p) < fMin) || (GetValue(p) > fMax ) ) return kTRUE; //rejected
  else return kFALSE; //accepted
}
/******************************************************************/
/******************************************************************/
/******************************************************************/

 
class AliHBTMomentumCut: public AliHbtBaseCut
 {
  public: 
    AliHBTMomentumCut(Double_t min = 0.0, Double_t max = 0.0):AliHbtBaseCut(min,max,kHbtP){}
    virtual ~AliHBTMomentumCut(){}
  protected:
    Double_t  GetValue(AliHBTParticle * p)const{return p->P();}
    ClassDef(AliHBTMomentumCut,1)
 };

class AliHBTPtCut: public AliHbtBaseCut
 {
  public: 
    AliHBTPtCut(Double_t min = 0.0, Double_t max = 0.0):AliHbtBaseCut(min,max,kHbtPt){}
    virtual ~AliHBTPtCut(){}
  protected:
    Double_t  GetValue(AliHBTParticle * p)const{return p->Pt();}
    ClassDef(AliHBTPtCut,1)
 };


class AliHBTEnergyCut: public AliHbtBaseCut
 {
  public: 
    AliHBTEnergyCut(Double_t min = 0.0, Double_t max = 0.0):AliHbtBaseCut(min,max,kHbtE){}
    virtual ~AliHBTEnergyCut(){}
  protected:
    Double_t  GetValue(AliHBTParticle * p)const {return p->Energy();}
    ClassDef(AliHBTEnergyCut,1)
 };

class AliHBTRapidityCut: public AliHbtBaseCut
 {
  public: 
    AliHBTRapidityCut(Double_t min = 0.0, Double_t max = 0.0):AliHbtBaseCut(min,max,kHbtRapidity){}
    virtual ~AliHBTRapidityCut(){}
  protected:
    Double_t  GetValue(AliHBTParticle * p)const{return p->Y();}
    ClassDef(AliHBTRapidityCut,1)
 };

class AliHBTPseudoRapidityCut: public AliHbtBaseCut
 {
  public: 
    AliHBTPseudoRapidityCut(Double_t min = 0.0, Double_t max = 0.0):AliHbtBaseCut(min,max,kHbtPseudoRapidity){}
    virtual ~AliHBTPseudoRapidityCut(){}
  protected:
    Double_t  GetValue(AliHBTParticle * p)const{return p->Eta();}
    ClassDef(AliHBTPseudoRapidityCut,1)
 };

class AliHBTPxCut: public AliHbtBaseCut
 {
  public: 
    AliHBTPxCut(Double_t min = 0.0, Double_t max = 0.0):AliHbtBaseCut(min,max,kHbtPx){}
    virtual ~AliHBTPxCut(){}
  protected:
    Double_t  GetValue(AliHBTParticle * p)const{return p->Px();}
    ClassDef(AliHBTPxCut,1)
 };

class AliHBTPyCut: public AliHbtBaseCut
 {
  public: 
    AliHBTPyCut(Double_t min = 0.0, Double_t max = 0.0):AliHbtBaseCut(min,max,kHbtPy){}
    virtual ~AliHBTPyCut(){}
  protected:
    Double_t  GetValue(AliHBTParticle * p)const{return p->Py();}
    ClassDef(AliHBTPyCut,1)
 };


class AliHBTPzCut: public AliHbtBaseCut
 {
  public: 
    AliHBTPzCut(Double_t min = 0.0, Double_t max = 0.0):AliHbtBaseCut(min,max,kHbtPz){}
    virtual ~AliHBTPzCut(){}
  protected:
    Double_t  GetValue(AliHBTParticle * p)const{return p->Pz();}
    ClassDef(AliHBTPzCut,1)
 };

class AliHBTPhiCut: public AliHbtBaseCut
 {
  public: 
    AliHBTPhiCut(Double_t min = 0.0, Double_t max = 0.0):AliHbtBaseCut(min,max,kHbtPhi){}
    virtual ~AliHBTPhiCut(){}
  protected:
    Double_t  GetValue(AliHBTParticle * p)const{return p->Phi();}
    ClassDef(AliHBTPhiCut,1)
  
 };

class AliHBTThetaCut: public AliHbtBaseCut
 {
  public: 
    AliHBTThetaCut(Double_t min = 0.0, Double_t max = 0.0):AliHbtBaseCut(min,max,kHbtTheta){}
    virtual ~AliHBTThetaCut(){}
  protected:
    Double_t  GetValue(AliHBTParticle * p)const{return p->Theta();}
    ClassDef(AliHBTThetaCut,1)
  
 };

class AliHBTVxCut: public AliHbtBaseCut
 {
 //Cut of the X coAnddinate of the vertex position
  public: 
    AliHBTVxCut(Double_t min = 0.0, Double_t max = 0.0):AliHbtBaseCut(min,max,kHbtVx){}
    virtual ~AliHBTVxCut(){}
  protected:
    Double_t  GetValue(AliHBTParticle * p)const{return p->Vx();} //retruns value of the vertex
    ClassDef(AliHBTVxCut,1)
  
 };


class AliHBTVyCut: public AliHbtBaseCut
 {
 //Cut of the X coAnddinate of the vertex position
  public: 
    AliHBTVyCut(Double_t min = 0.0, Double_t max = 0.0):AliHbtBaseCut(min,max,kHbtVy){}
    virtual ~AliHBTVyCut(){}
  protected:
    Double_t  GetValue(AliHBTParticle * p)const{return p->Vy();} //retruns value of the vertex
    ClassDef(AliHBTVyCut,1)
  
 };

class AliHBTVzCut: public AliHbtBaseCut
 {
 //Cut of the X coAnddinate of the vertex position
  public: 
    AliHBTVzCut(Double_t min = 0.0, Double_t max = 0.0):AliHbtBaseCut(min,max,kHbtVz){}
    virtual ~AliHBTVzCut(){}
  protected:
    Double_t  GetValue(AliHBTParticle * p)const{return p->Vz();} //retruns value of the vertex
    
    ClassDef(AliHBTVzCut,1)
  
 };

class AliHBTPIDCut:  public AliHbtBaseCut
 {
   public:
     AliHBTPIDCut():AliHbtBaseCut(0.0,0.0,kHbtPid),fPID(0){}
     AliHBTPIDCut(Int_t pid, Double_t min = 0.0, Double_t max = 1.0):AliHbtBaseCut(min,max,kHbtPid),fPID(pid){}
     virtual ~AliHBTPIDCut(){}
   protected:
     Double_t  GetValue(AliHBTParticle * p)const{return p->GetPIDprobability(fPID);}
     Int_t     fPID; //pid of particle that the pid is set 
     ClassDef(AliHBTPIDCut,1)
 };
//___________________________________________________
/////////////////////////////////////////////////////
//                                                 //
// class AliHBTLogicalOperCut                      //
//                                                 //
// This cut is base class fAnd class that perfAndms  //
// logical operations on cuts                      //
//                                                 //
/////////////////////////////////////////////////////
class AliHBTLogicalOperCut:  public AliHbtBaseCut
 {
   public:
     AliHBTLogicalOperCut();
     AliHBTLogicalOperCut(AliHbtBaseCut* first, AliHbtBaseCut* second);
     virtual   ~AliHBTLogicalOperCut();
   protected:
     Double_t  GetValue(AliHBTParticle * /*part*/) const {MayNotUse("GetValue");return 0.0;}
     
     AliHbtBaseCut* fFirst;   //second cut
     AliHbtBaseCut* fSecond;  //first cut
   private:  
    class  AliHBTDummyBaseCut: public AliHbtBaseCut 
     {
       Double_t  GetValue(AliHBTParticle * /*part*/) const {return 0.0;}
       Bool_t    Pass(AliHBTParticle* /*part*/) const;
     };
     
    ClassDef(AliHBTLogicalOperCut,1)
 };

class AliHBTOrCut: public AliHBTLogicalOperCut
{
   public:
     AliHBTOrCut(){}
     AliHBTOrCut(AliHbtBaseCut* first, AliHbtBaseCut* second):AliHBTLogicalOperCut(first,second){}
     virtual   ~AliHBTOrCut(){}
     Bool_t    Pass(AliHBTParticle *p) const;
     ClassDef(AliHBTOrCut,1)
};

class AliHBTAndCut: public AliHBTLogicalOperCut
{
   public:
     AliHBTAndCut(){}
     AliHBTAndCut(AliHbtBaseCut* first, AliHbtBaseCut* second):AliHBTLogicalOperCut(first,second){}
     virtual   ~AliHBTAndCut(){}
     Bool_t    Pass(AliHBTParticle *p) const;
     ClassDef(AliHBTAndCut,1)
};

#endif
