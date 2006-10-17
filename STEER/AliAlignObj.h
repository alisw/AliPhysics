#ifndef ALIALIGNOBJ_H
#define ALIALIGNOBJ_H

//************************************************************************
// AliAlignObj: alignment base class for the storage of the alignment    *
//   constants for a single volume:                                      *
//   -  a displacement (a shift and a rotation) either as                *
//      - the 6 doubles which identify it or as                          *
//      - the matrix which identifies it                                 *
//   -  the identity of the volume itself in form of a symbolic volume   *
//      name for alignable volumes, in form of a TGeo path otherwise,    *
//      and as a unique integer identifier                               *
//************************************************************************
#include "TObject.h"
#include "TString.h"
#include "TGeoMatrix.h"

class AliTrackPoint;
class AliTrackPointArray;

class AliAlignObj : public TObject {

 public:

  enum ELayerID{kInvalidLayer=0,
		kFirstLayer=1,
		kSPD1=1, kSPD2=2,
		kSDD1=3, kSDD2=4,
		kSSD1=5, kSSD2=6,
		kTPC1=7, kTPC2=8,
		kTRD1=9, kTRD2=10, kTRD3=11, kTRD4=12, kTRD5=13, kTRD6=14,
		kTOF=15,
		kPHOS1=16, kPHOS2=17,
		kRICH=18,
		kMUON=19,
		kLastLayer=20};
  AliAlignObj();
  AliAlignObj(const char* symname, UShort_t voluid);
  AliAlignObj(const char* symname, ELayerID detId, Int_t modId);
  AliAlignObj(const AliAlignObj& theAlignObj);
  AliAlignObj& operator= (const AliAlignObj& theAlignObj);
  AliAlignObj& operator*=(const AliAlignObj& theAlignObj);
  virtual ~AliAlignObj();

  //Setters
  virtual void SetTranslation(Double_t x, Double_t y, Double_t z) = 0;
  virtual void SetTranslation(const TGeoMatrix& m) = 0;
  virtual void SetRotation(Double_t psi, Double_t theta, Double_t phi) = 0;
  virtual Bool_t SetRotation(const TGeoMatrix& m) = 0;
  virtual void SetPars(Double_t x, Double_t y, Double_t z, Double_t psi,
               Double_t theta, Double_t phi);
  virtual Bool_t SetLocalPars(Double_t x, Double_t y, Double_t z,
			      Double_t psi, Double_t theta, Double_t phi);
  virtual Bool_t SetMatrix(const TGeoMatrix& m);
  virtual Bool_t SetLocalMatrix(const TGeoMatrix& m);
  void  SetSymName(const TString& symname) {fVolPath=symname;}
  void  SetVolUID(UShort_t voluid) {fVolUID=voluid;}
  void  SetVolUID(ELayerID layerId, Int_t modId);

  //Getters
  const char  *GetSymName()    const {return fVolPath.Data();}
  UShort_t     GetVolUID()     const {return fVolUID;}
  void         GetVolUID(ELayerID &layerId, Int_t &modId) const;
  virtual void GetTranslation(Double_t* tr)  const=0;
  virtual Bool_t GetAngles(Double_t* angles) const=0;
  virtual Bool_t GetPars(Double_t transl[], Double_t rot[]) const;
  virtual void GetMatrix(TGeoHMatrix& m) const=0;

  Bool_t   IsSortable() const {return kTRUE;}
  Int_t         GetLevel() const;
  virtual Int_t Compare(const TObject* obj) const;

  virtual AliAlignObj& Inverse() const=0;

  void  Transform(AliTrackPoint &p) const;
  void  Transform(AliTrackPointArray &array) const;

  void  Print(Option_t *) const;

  static Int_t       LayerSize(Int_t layerId);
  static const char* LayerName(Int_t layerId);

  static UShort_t LayerToVolUID(ELayerID layerId, Int_t modId);
  static UShort_t LayerToVolUID(Int_t    layerId, Int_t modId);
  static ELayerID VolUIDToLayer(UShort_t voluid, Int_t &modId);
  static ELayerID VolUIDToLayer(UShort_t voluid);

  static const char* SymName(UShort_t voluid);
  static const char* SymName(ELayerID layerId, Int_t modId);

  Bool_t ApplyToGeometry();
  static Bool_t   GetFromGeometry(const char *symname, AliAlignObj &alobj);

  static AliAlignObj* GetAlignObj(UShort_t voluid);
  static AliAlignObj* GetAlignObj(ELayerID layerId, Int_t modId);

 protected:

  void AnglesToMatrix(const Double_t *angles, Double_t *rot) const;
  Bool_t MatrixToAngles(const Double_t *rot, Double_t *angles) const;

  static void InitSymNames();
  static void InitAlignObjFromGeometry();

  //Volume identifiers
  TString  fVolPath; // Symbolic volume name; in case could coincide with
      // the volume path inside TGeo geometry (for non-alignable volumes)
  UShort_t fVolUID;  // Unique volume ID

  static Int_t       fgLayerSize[kLastLayer - kFirstLayer]; // Size of layers
  static const char* fgLayerName[kLastLayer - kFirstLayer]; // Name of layers

  static TString*    fgVolPath[kLastLayer - kFirstLayer]; // Symbolic volume names
  static AliAlignObj** fgAlignObjs[kLastLayer - kFirstLayer]; // Alignment objects

  ClassDef(AliAlignObj, 2)
};

#endif
