#ifndef ALIEVEEMCALSMODULEDATA_H
#define ALIEVEEMCALSMODULEDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \class AliEveEMCALSModuleData
/// \brief EMCal event display handling of super-modules data
///
/// Store the data related to each Super Module (SM)
/// Possible storage of hits, digits and clusters per SM
///
/// \author Magali Estienne <magali.estienne@cern.ch>, SUBATECH. EMCal implementation, June 2008
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS. DCal implementation + doxygen, May 2015.
///

#include <vector>

#include <AliEveEMCALData.h>

class TClonesArray;
class TGeoNode;
class TGeoMatrix;
class TStyle;
class TGedFrame;
class TBuffer3DTypes;
class TBuffer3D;
class TVirtualPad;
class TVirtualViewer3D;

class TEvePointSet;
class TEveQuadSet;
class TEveBoxSet;
class TEveFrameBox;
class TEvePointSet;
class TEveTrans;
class TTree;

class AliEMCALGeometry;

class AliEveEMCALSModuleData : public TObject
{
public:
  
  AliEveEMCALSModuleData(Int_t chamber, AliEMCALGeometry* geom, TGeoNode* node); //, TGeoHMatrix* m);
  
  virtual ~AliEveEMCALSModuleData();

  void        DropData();
  
  void        Init(Int_t sm);
  
  void        RegisterDigit  (Int_t AbsId, Int_t isupMod, Double_t iamp, Double_t ix, Double_t iy, Double_t iz);
  void        RegisterCluster(Int_t isupMod, Double_t iamp, Double_t ix, Double_t iy, Double_t iz); 
  void        RegisterHit    (Int_t AbsId, Int_t isupMod, Double_t iamp, Double_t ix, Double_t iy, Double_t iz); 
  
  Int_t       GetNDigits()   const { return fNDigits   ; } 
  Int_t       GetNClusters() const { return fNClusters ; } 
  Int_t       GetNHits()     const { return fNHits     ; } 
  
  Int_t       GetSmId()      const { return fSmId      ; }
  Int_t       GetNsm()       const { return fNsm       ; }

  std::vector< std::vector<Double_t> > GetDigitBuffer()   const { return fDigitArray;   };  
  std::vector< std::vector<Double_t> > GetClusterBuffer() const { return fClusterArray; };  
  std::vector< std::vector<Float_t> >  GetHitBuffer()     const { return fHitArray;     };  

  void        GetSModuleBigBox(Float_t& bbox0, Float_t& bbox1, Float_t& bbox2) 
  const { bbox0 = fgSModuleBigBox0; bbox1 = fgSModuleBigBox1; bbox2 = fgSModuleBigBox2;}
  
  void        GetSModuleSmallBox(Float_t& bbox0, Float_t& bbox1, Float_t& bbox2) 
  const { bbox0 = fgSModuleSmallBox0; bbox1 = fgSModuleSmallBox1; bbox2 = fgSModuleSmallBox2;}

  void        GetSModuleDCalBox(Float_t& bbox0, Float_t& bbox1, Float_t& bbox2) 
  const { bbox0 = fgSModuleDCalBox0; bbox1 = fgSModuleDCalBox1; bbox2 = fgSModuleDCalBox2;}
  
  void        GetSModuleSmallDBox(Float_t& bbox0, Float_t& bbox1, Float_t& bbox2) 
  const { bbox0 = fgSModuleSmallDBox0; bbox1 = fgSModuleSmallDBox1; bbox2 = fgSModuleSmallDBox2;}
  
  Float_t     GetPhiTileSize()   const { return fPhiTileSize ; }
  Float_t     GetEtaTileSize()   const { return fEtaTileSize ; }
  
  //TGeoMatrix* GetSModuleMatrix() const { return fMatrix      ; }
  TGeoMatrix* GetSModuleMatrix(Int_t sm) const { return (TGeoMatrix*) fNode->GetDaughter(sm)->GetMatrix(); }
  
 protected:
  
  AliEMCALGeometry* fGeom;                 ///< Data member to set/call EMCAL geometry
  TGeoNode*         fNode;                 ///< Node for bbox definition
  Int_t             fSmId;                 ///< number of the chamber, 0 to 19 
  Int_t             fNsm;                  ///< Total number of super modules
  Int_t             fNDigits;              ///< number of found digits 
  Int_t             fNClusters;            ///< number of found rec points 
  Int_t             fNHits;                ///< number of simulation hits 
  
  Float_t           fPhiTileSize;          ///< Typical phi size of a QuadSet (digit)
  Float_t           fEtaTileSize;          ///< Typical eta size of a QuadSet (digit)
  
  std::vector< std::vector<Float_t> >  fHitArray;     //|| Hit coordinates, etc.
  std::vector< std::vector<Double_t> > fDigitArray;   //|| Digit coordinates, etc.
  std::vector< std::vector<Double_t> > fClusterArray; //|| Rec point coordinates, etc.
  
  static Float_t    fgSModuleBigBox0;      ///< SM envelope box, full EMCAL
  static Float_t    fgSModuleBigBox1;      ///< SM envelope box, full EMCAL
  static Float_t    fgSModuleBigBox2;      ///< SM envelope box, full EMCAL
  
  static Float_t    fgSModuleSmallBox0;    ///< SM envelope box, 1/3 EMCAL
  static Float_t    fgSModuleSmallBox1;    ///< SM envelope box, 1/3 EMCAL
  static Float_t    fgSModuleSmallBox2;    ///< SM envelope box, 1/3 EMCAL

  static Float_t    fgSModuleDCalBox0;     ///< SM envelope box, full DCAL
  static Float_t    fgSModuleDCalBox1;     ///< SM envelope box, full DCAL
  static Float_t    fgSModuleDCalBox2;     ///< SM envelope box, full DCAL
  
  static Float_t    fgSModuleSmallDBox0;   ///< SM envelope box, 1/3 DCAL
  static Float_t    fgSModuleSmallDBox1;   ///< SM envelope box, 1/3 DCAL
  static Float_t    fgSModuleSmallDBox2;   ///< SM envelope box, 1/3 DCAL
  
  //  TGeoMatrix*       fMatrix;               ///< Matrix for local to global transformation (needed?)
  //  TGeoHMatrix*      fHMatrix;              ///< Matrix for local to global transformation (needed?)

 private:
  
  AliEveEMCALSModuleData           (const AliEveEMCALSModuleData& esmdata); 
  
  /// Assignment operator not implemented.
  AliEveEMCALSModuleData& operator=(const AliEveEMCALSModuleData& esmdata); 
  
  /// \cond CLASSIMP
  ClassDef(AliEveEMCALSModuleData, 2); 
  /// \endcond

};

#endif //ALIEVEEMCALSMODULEDATA_H
