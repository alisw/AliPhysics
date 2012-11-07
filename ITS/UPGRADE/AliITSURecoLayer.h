#ifndef ALIITSURECOLAYER
#define ALIITSURECOLAYER

#include <TNamed.h>
#include <TObjArray.h>
class AliITSUGeomTGeo;
class AliITSsegmentation;
class AliITSURecoSens;

///////////////////////////////////////////////////////////////////////
//                                                                   //
//  Class AliITSURecoLayer                                           //
//  Interface between the framework and reconstruction for ITS layer //
//                                                                   //
///////////////////////////////////////////////////////////////////////


class AliITSURecoLayer : public TNamed
{
 public:
  //
  enum {kPassive=BIT(14)};
  AliITSURecoLayer(const char* name);
  AliITSURecoLayer(const char* name, Int_t activeID,Int_t nsens,AliITSUGeomTGeo* gm);
  virtual ~AliITSURecoLayer();
  //
  Bool_t             Build();
  //
  Int_t              GetID()                       const {return (int)GetUniqueID();}
  Int_t              GetActiveID()                 const {return fActiveID;}
  Int_t              GetNSensors()                 const {return fNSensors;}
  Double_t           GetRMin()                     const {return fRMin;}
  Double_t           GetRMax()                     const {return fRMax;}
  Double_t           GetR()                        const {return fR;}
  Bool_t             IsActive()                    const {return !TestBit(kPassive);}
  Bool_t             IsPassive()                   const {return TestBit(kPassive);}
  //
  void               SetID(Int_t i)                      {SetUniqueID(i);} 
  void               SetActiveID(Int_t i)                {fActiveID = i;} 
  void               SetRMin(Double_t r)                 {fRMin = r;}
  void               SetRMax(Double_t r)                 {fRMax = r;}
  void               SetR(Double_t r)                    {fR = r;}
  void               SetPassive(Bool_t v=kTRUE)          {SetBit(kPassive,v);}

  //
  AliITSURecoSens*   GetSensor(Int_t i)            const {return (AliITSURecoSens*)fSensors.UncheckedAt(i);}
  AliITSURecoSens*   GetSensor(Int_t ld,Int_t is)  const {return GetSensor(ld*fNSensInLadder+is);}
  void               AddSensor(const AliITSURecoSens* mod); 
  //
  Int_t              FindSensors(const double* impPar, AliITSURecoSens **sensors);
  //
  virtual void       Print(Option_t* option = "")  const;

 protected:
  Int_t              fActiveID;  // ID within active layers
  Int_t              fNSensors;  // N of modules
  Int_t              fNSensInLadder; // N sensors in the ladder
  Int_t              fNLadders;  // N ladder
  Double_t           fR;         // mean R
  Double_t           fRMax;      // max  R
  Double_t           fRMin;      // min  R
  Double_t           fZMax;      // max  Z
  Double_t           fZMin;      // min  Z
  Double_t*          fPhiLadMax; // max lab phi of the ladder
  Double_t*          fPhiLadMin; // min lab phi of the ladder
  Double_t           fPhiOffs;   // offset in phi for 1st ladder
  Double_t           fSensDZInv; // inverse mean sensor Z span
  Double_t           fDPhiLadInv;// inverse mean ladder dphi
  TObjArray          fSensors;   // sensors
  AliITSUGeomTGeo*   fITSGeom;   // geometry interface
  //
 private:
  AliITSURecoLayer(const AliITSURecoLayer &source); 
  AliITSURecoLayer& operator=(const AliITSURecoLayer &source); 
  //

  ClassDef(AliITSURecoLayer,1); // helper for layer data used in reco
};

#endif
