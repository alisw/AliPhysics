#ifndef ALIITSURECOLAYER
#define ALIITSURECOLAYER

#include <TNamed.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include "AliITSURecoSens.h"

class AliITSUGeomTGeo;
class AliITSsegmentation;
class AliCluster;


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
  enum {kPassive=BIT(14)                 // layer is passive
	,kOwnsClusterArray=BIT(15)       // owner of cluster array, delete in destructor
  };
  AliITSURecoLayer(const char* name);
  AliITSURecoLayer(const char* name, Int_t activeID,AliITSUGeomTGeo* gm);
  virtual ~AliITSURecoLayer();
  //
  void               ProcessClusters(Int_t mode=0);
  //
  Int_t              GetID()                       const {return (int)GetUniqueID();}
  Int_t              GetActiveID()                 const {return fActiveID;}
  Int_t              GetNSensors()                 const {return fNSensors;}
  Int_t              GetNLadders()                 const {return fNLadders;}
  Int_t              GetNSensInLadder()            const {return fNSensInLadder;}
  Double_t           GetRMin()                     const {return fRMin;}
  Double_t           GetRMax()                     const {return fRMax;}
  Double_t           GetDR()                       const {return fRMax-fRMin;}
  Double_t           GetR()                        const {return fR;}
  Double_t           GetR(Int_t dir)               const {return dir==0 ? fR : (dir<0 ? fRMin : fRMax);}
  Double_t           GetZMin()                     const {return fZMin;}
  Double_t           GetZMax()                     const {return fZMax;}
  Double_t           GetMaxStep()                  const {return fMaxStep;}
  Bool_t             IsActive()                    const {return !TestBit(kPassive);}
  Bool_t             IsPassive()                   const {return TestBit(kPassive);}
  Bool_t             GetOwnsClusterArray()         const {return TestBit(kOwnsClusterArray);}
  //
  void               SetID(Int_t i)                      {SetUniqueID(i);} 
  void               SetActiveID(Int_t i)                {fActiveID = i;} 
  void               SetRMin(Double_t r)                 {fRMin = r;}
  void               SetRMax(Double_t r)                 {fRMax = r;}
  void               SetR(Double_t r)                    {fR = r;}
  void               SetZMin(Double_t z)                 {fZMin = z;}
  void               SetZMax(Double_t z)                 {fZMax = z;}
  void               SetPassive(Bool_t v=kTRUE)          {SetBit(kPassive,v);}
  void               SetOwnsClusterArray(Bool_t v=kTRUE) {SetBit(kOwnsClusterArray,v);}
  void               SetMaxStep(Double_t st)             {fMaxStep = st>0 ? st : 0.1;}
  //
  AliITSURecoSens*   GetSensor(Int_t i)            const {return i<0 ? 0:(AliITSURecoSens*)fSensors[i];}
  AliITSURecoSens*   GetSensor(Int_t ld,Int_t is)  const {return GetSensor(ld*fNSensInLadder+is);}
  AliITSURecoSens*   GetSensorFromID(Int_t i)      const;
  TClonesArray*      GetClusters()                 const {return (TClonesArray*)fClusters;}
  TClonesArray**     GetClustersAddress()                {return (TClonesArray**)&fClusters;}  
  Int_t              GetNClusters()                const {return fClusters ? fClusters->GetEntriesFast() : 0;}
  AliCluster*        GetCluster(Int_t icl)         const {return (AliCluster*)fClusters->UncheckedAt(icl);}
  void               SetClusters(TClonesArray* cl)       {fClusters = cl;}
  //
  Int_t              FindSensors(const double* impPar, AliITSURecoSens *sensors[AliITSURecoSens::kNNeighbors]);
  //
  virtual void       Print(Option_t* option = "")  const;
  virtual Bool_t     IsSortable()                  const {return kTRUE;}
  virtual Bool_t     IsEqual(const TObject* obj)   const;
  virtual Int_t      Compare(const TObject* obj)   const;
  //
 protected:
  void               Build();
  //
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
  Double_t           fMaxStep;   // max step in tracking X allowed within layer
  AliITSURecoSens**  fSensors;   // sensors
  AliITSUGeomTGeo*   fITSGeom;   // geometry interface
  TClonesArray*      fClusters;  // clusters of the layer
  //
 private:
  AliITSURecoLayer(const AliITSURecoLayer &source); 
  AliITSURecoLayer& operator=(const AliITSURecoLayer &source); 
  //

  ClassDef(AliITSURecoLayer,1); // helper for layer data used in reco
};

#endif
