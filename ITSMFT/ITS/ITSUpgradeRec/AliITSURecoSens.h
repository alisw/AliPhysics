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
  //
  enum {kLeft=BIT(1),kRight=BIT(2),kUp=BIT(3),kDown=BIT(4)};
  //
  AliITSURecoSens(Int_t id=0);
  AliITSURecoSens(const AliITSURecoSens &source); 
  virtual ~AliITSURecoSens() {}
  AliITSURecoSens& operator=(const AliITSURecoSens &source); 
  //
  Int_t              GetID()                       const {return (int)GetUniqueID();}
  Int_t              CheckCoverage(double phi, double z) const;
  Double_t           GetXTF()                      const {return fXTF;}
  Double_t           GetPhiTF()                    const {return fPhiTF;}
  Double_t           GetPhiMin()                   const {return fPhiMin;}
  Double_t           GetPhiMax()                   const {return fPhiMax;}
  Double_t           GetZMin()                     const {return fZMin;}
  Double_t           GetZMax()                     const {return fZMax;}
  //
  Int_t              GetNClusters()                const {return fNClusters;}
  Int_t              GetFirstClusterId()           const {return fFirstClusterId;}
  //
  void               SetID(Int_t i)                      {SetUniqueID(i);}
  void               SetXTF(double v)                    {fXTF = v;}
  void               SetPhiTF(double v)                  {fPhiTF = v;}
  void               SetBoundaries(double phiMn,double phiMx, double zMn, double zMx);
  //
  void               SetNClusters(Int_t ncl)             {fNClusters = ncl;}
  void               IncNClusters()                      {fNClusters++;}
  void               SetFirstClusterId(Int_t id)         {fFirstClusterId = id;}
  void               ResetClusters(); 
  void               ProcessClusters(Int_t mode=0);
  //
  virtual void       Print(Option_t* option = "")  const;
  //
  virtual Bool_t     IsSortable()                 const {return kTRUE;}
  virtual Int_t	     Compare(const TObject* obj)  const;
  virtual Bool_t     IsEqual(const TObject* obj)  const {return Compare(obj)==0;}
  //
 protected:
  Int_t              fNClusters;                   // number of clusters
  Int_t              fFirstClusterId;              // index of the 1st cluster in the layer's clusters array
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
