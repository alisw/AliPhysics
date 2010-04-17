//*************************************************************************
// EMCAL event display
// Store the data related to each Super Module
// For visualization
//
//  Author: Magali Estienne (magali.estienne@cern.ch)
//  June 30 2008
//*************************************************************************

#ifndef ALIEVEEMCALSMODULEDATA_H
#define ALIEVEEMCALSMODULEDATA_H

#include <vector>

#include <TObject.h>
#include <TClonesArray.h>
#include <TGeoNode.h>
#include <TGeoMatrix.h>
#include <TEvePointSet.h>

class AliEMCALGeometry;

class AliEveEMCALSModuleData : public TObject
{
public:
  AliEveEMCALSModuleData(Int_t chamber,AliEMCALGeometry* geom,TGeoNode* node, TGeoHMatrix* m);
  virtual ~AliEveEMCALSModuleData();

  void        DropData();
  void        Init(Int_t sm);
  void        RegisterDigit(Int_t AbsId, Int_t isupMod, Double_t iamp, Double_t ix, Double_t iy, Double_t iz);
  void        RegisterCluster(Int_t isupMod, Double_t iamp, Double_t ix, Double_t iy, Double_t iz); 
  void        RegisterHit(Int_t AbsId, Int_t isupMod, Double_t iamp, Double_t ix, Double_t iy, Double_t iz); 
  Int_t       GetNDigits()   const { return fNDigits;   }; 
  Int_t       GetNClusters() const { return fNClusters; }; 
  Int_t       GetNHits()     const { return fNHits;     }; 
  Int_t       GetSmId()      const { return fSmId; };
  Int_t       GetNsm()       const {return fNsm;};
  Int_t       GetNsmf()      const {return fNsmfull;};
  Int_t       GetNsmh()      const {return fNsmhalf;};
  std::vector< std::vector<Double_t> > GetDigitBuffer() const  { return fDigitArray;   };  
  std::vector< std::vector<Double_t> > GetClusterBuffer() const { return fClusterArray; };  
  std::vector< std::vector<Float_t> > GetHitBuffer() const { return fHitArray;     };  

  void        GetSModuleBigBox(Float_t& bbox0, Float_t& bbox1, Float_t& bbox2) 
  const { bbox0 = fgSModuleBigBox0; bbox1 = fgSModuleBigBox1; bbox2 = fgSModuleBigBox2;};
  void        GetSModuleSmallBox(Float_t& bbox0, Float_t& bbox1, Float_t& bbox2) 
  const { bbox0 = fgSModuleSmallBox0; bbox1 = fgSModuleSmallBox1; bbox2 = fgSModuleSmallBox2;};
  void        GetSModuleCenter(Float_t& bboxCenter0, Float_t& bboxCenter1, Float_t& bboxCenter2) 
  const { bboxCenter0 = fgSModuleCenter0; bboxCenter1 = fgSModuleCenter1; bboxCenter2 = fgSModuleCenter2;};
  Float_t     GetPhiTileSize() const  {return fPhiTileSize;};
  Float_t     GetEtaTileSize() const  {return fEtaTileSize;};
  TGeoMatrix* GetSModuleMatrix() const {return fMatrix;};
  
 protected:
  AliEMCALGeometry* fGeom;                 // Data member to set/call EMCAL geometry
  TGeoNode*         fNode;                 // Node for bbox definition
  Int_t             fSmId;                 // number of the chamber, 0 to 13 
  Int_t             fNsm;                  // Total number of super modules
  Int_t             fNsmfull;              // Number of full super modules
  Int_t             fNsmhalf;              // Number of half super modules
  Int_t             fNDigits;              // number of found digits 
  Int_t             fNClusters;            // number of found rec points 
  Int_t             fNHits;                // number of simulation hits 
  
  Float_t           fPhiTileSize;          // Typical phi size of a QuadSet (digit)
  Float_t           fEtaTileSize;          // Typical eta size of a QuadSet (digit)
  
  std::vector< std::vector<Float_t> > fHitArray;     //|| Hit coordinates, etc.
  std::vector< std::vector<Double_t> > fDigitArray;   //|| Digit coordinates, etc.
  std::vector< std::vector<Double_t> > fClusterArray; //|| Rec point coordinates, etc.
  
  static Float_t    fgSModuleBigBox0;       // SM envelope box
  static Float_t    fgSModuleBigBox1;       // SM envelope box
  static Float_t    fgSModuleBigBox2;       // SM envelope box
  static Float_t    fgSModuleSmallBox0;     // SM envelope box
  static Float_t    fgSModuleSmallBox1;     // SM envelope box
  static Float_t    fgSModuleSmallBox2;     // SM envelope box
  static Float_t    fgSModuleCenter0;       // SM envelope box
  static Float_t    fgSModuleCenter1;       // SM envelope box
  static Float_t    fgSModuleCenter2;       // SM envelope box

  TGeoMatrix*   fMatrix;                   // Matrix for local to global transformation
  TGeoHMatrix*  fHMatrix;                  // Matrix for local to global transformation

 private:
  AliEveEMCALSModuleData(const AliEveEMCALSModuleData& esmdata);            
  AliEveEMCALSModuleData& operator=(const AliEveEMCALSModuleData&); // Not implemented
  
  ClassDef(AliEveEMCALSModuleData, 0);     // class with data for one chamber
};

#endif
