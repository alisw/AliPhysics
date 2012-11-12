#ifndef ALIITSURECOSENS
#define ALIITSURECOSENS

#include <TObject.h>

///////////////////////////////////////////////////////////////////////
//                                                                   //
//  Class AliITSURecoSens                                            //
//  Interface between the framework and reconstruction for           //
//  single ITS sensor                                                //
//                                                                   //
///////////////////////////////////////////////////////////////////////

class AliITSURecoSens : public TObject
{
 public:
  enum {kNghbR,kNghbTR,kNghbT,kNghbTL,kNghbL,kNghbBL,kNghbB,kNghbBR,kNNeighbors}; // neighbors: Top,Left etc
  //
  AliITSURecoSens(Int_t id);
  AliITSURecoSens(const AliITSURecoSens &source); 
  virtual ~AliITSURecoSens() {}
  AliITSURecoSens& operator=(const AliITSURecoSens &source); 
  //
  Int_t              GetID()                       const {return (int)GetUniqueID();}
  Double_t           GetXTF()                      const {return fXTF;}
  Double_t           GetPhiTF()                    const {return fPhiTF;}
  Double_t           GetPhiMin()                   const {return fPhiMin;}
  Double_t           GetPhiMax()                   const {return fPhiMax;}
  Double_t           GetZMin()                     const {return fZMin;}
  Double_t           GetZMax()                     const {return fZMax;}
  //
  Int_t              GetNeighborID(int i)          const {return fNeighbors[i];}
  //
  void               SetID(Int_t i)                      {SetUniqueID(i);}
  void               SetXTF(double v)                    {fXTF = v;}
  void               SetPhiTF(double v)                  {fPhiTF = v;}
  void               SetNeighborID(int i, int id)        {fNeighbors[i] = id;}
  void               SetBoundaries(double phiMn,double phiMx, double zMn, double zMx);
  //
  virtual void       Print(Option_t* option = "")  const;

 protected:
  Int_t              fNeighbors[kNNeighbors];      // id of neighbors  
  Double_t           fXTF;                         // X in tracking frame
  Double_t           fPhiTF;                       // phi of tracking frame
  Double_t           fPhiMin;                      // lab phi min
  Double_t           fPhiMax;                      // lab phi max
  Double_t           fZMin;                        // lab & trk Z min
  Double_t           fZMax;                        // lab & trk Z max
  //
  ClassDef(AliITSURecoSens,1); // helper for sensor data used in reco
};


#endif
