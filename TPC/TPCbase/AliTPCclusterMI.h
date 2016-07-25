#ifndef ALITPCCLUSTERMI_H
#define ALITPCCLUSTERMI_H

/// \class AliTPCclusterMI
/// \brief TPC Cluster Class
///
/// Parallel tracking
///
/// \author Marian Ivanov

/* $Id$ */


#include "AliCluster.h"
#include "TMath.h"
//#include "AliTPCclusterInfo.h"
#include <AliTrackPointArray.h>

//_____________________________________________________________________________
class AliTPCclusterMI : public AliCluster {
  enum Status{ kDisabled = 0x7F};
  enum {
    kSectorChanged=BIT(14)          // to flag sector change due to the distortions
  };
public:
  AliTPCclusterMI();
  AliTPCclusterMI(const AliTPCclusterMI & cluster);
  AliTPCclusterMI &operator = (const AliTPCclusterMI & cluster); //assignment operator
  AliTPCclusterMI(Int_t *lab, Float_t *hit);
  virtual ~AliTPCclusterMI();
  virtual void	Clear(const Option_t*) { ResetBit(0xffffffff);}
  virtual Bool_t IsSortable() const;
  virtual Int_t Compare(const TObject* obj) const;
  inline  void Use(Int_t inc=10);
  inline  void Disable(){fUsed=kDisabled;}
  inline  Bool_t IsDisabled() const {return (fUsed==kDisabled);}

  Bool_t  IsSectorChanged()                const {return TestBit(kSectorChanged);}
  void    SetSectorChanged(Bool_t v=kTRUE)       {SetBit(kSectorChanged,v);}

  virtual Int_t GetDetector() const {return fDetector;}
  virtual Int_t GetRow() const {return fRow;}
  virtual void SetDetector(Int_t detector);
  virtual void SetRow(Int_t row){fRow = (UChar_t)(row%256);}
  virtual void SetTimeBin(Float_t timeBin){ fTimeBin= timeBin;}
  virtual void SetPad(Float_t pad){ fPad = pad;}
  //
  void SetQ(Float_t q) {fQ=(UShort_t)q;}
  void SetType(Char_t type) {fType=type;}
  void SetMax(UShort_t max) {fMax=max;}
  Int_t IsUsed(Int_t th=10) const {return (fUsed>=th) ? 1 : 0;}
  Float_t GetQ() const {return TMath::Abs(fQ);}
  Float_t GetMax() const {return fMax;}
  Char_t  GetType()const {return fType;}
  Float_t GetTimeBin() const { return fTimeBin;}
  Float_t GetPad() const { return fPad;}
  //
  void    SetDistortions(float dx, float dy, float dz);
  void    GetDistortions(float& dx,float& dy, float& dz)  const;
  void    SetDistortionDispersion(float d);
  Float_t GetDistortionX() const;
  Float_t GetDistortionY() const;
  Float_t GetDistortionZ() const;
  Float_t GetDistortionDispersion() const;

  Bool_t  GetGlobalCov(Float_t cov[6]) const;
  //  AliTPCclusterInfo * GetInfo() const { return fInfo;}
  //  void SetInfo(AliTPCclusterInfo * info);
  //
  AliTPCclusterMI*  MakeCluster(AliTrackPoint* point);
  AliTrackPoint*    MakePoint();
  static void     SetGlobalTrackPoint(const AliCluster &cl, AliTrackPoint &point);

 protected:
  enum{ // constants for storing x,y,z distortion in AliCluster::fSigmaYZ 
    kScaleDX=50,kScaleDY=100,kScaleDZ=100,kScaleDisp=85, // 1./kScale gives rounding in cm
    kNBitsDX=10, kNBitsDY=11,kNBitsDZ=11, 
    kMaxDX = (0x1<<(kNBitsDX-1))-1, kMaxDY = (0x1<<(kNBitsDY-1))-1,kMaxDZ = (0x1<<(kNBitsDZ-1))-1,
    kMaxDisp = 0xff,
    kMaskDX = (0x1<<kNBitsDX)-1, kMaskDY = (0x1<<kNBitsDY)-1, kMaskDZ = (0x1<<kNBitsDZ)-1
  };

private:
  //  AliTPCclusterInfo * fInfo;  ///< pointer to the cluster debug info
  Float_t   fTimeBin;  ///< time bin coordinate
  Float_t   fPad;  ///< pad coordinate
  Short_t   fQ ;       ///< Q of cluster (in ADC counts)
  Short_t   fMax;      ///< maximal amplitude in cluster
  Char_t    fType;     ///< different meaning depending on whether hlt of offline cluster finder was used:
                       //   offline: type of the cluster 0 means golden
                       //   hlt: 0 = Not Split, 1 (bit 1) = Split in Pad Direction, 2 (bit 2) = Split in Time Direction, 3 = Split in both directions
                       //        bit 3 and 8: Edge clusters. Bit 8 is set to make the value of fType negative, to match the offline behavior
  Char_t    fUsed;     ///< counter of usage
  UChar_t   fDisp;     ///< dispersion of applied correction
  UChar_t   fDetector; ///< detector  number
  UChar_t   fRow;      ///< row number number
  /// \cond CLASSIMP
  ClassDef(AliTPCclusterMI,7)  // Time Projection Chamber clusters
  /// \endcond
};

void AliTPCclusterMI::Use(Int_t inc)
{
  if (inc>0)  fUsed+=inc;
  else
    fUsed=0;
}



#endif


