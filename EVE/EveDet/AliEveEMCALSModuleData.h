//*************************************************************************
// EMCAL event display
// Store the data related to each Super Module
//
//  Author: Magali Estienne (magali.estienne@cern.ch)
//  June 30 2008
//*************************************************************************

#ifndef AliEveEMCALSModuleData_H
#define AliEveEMCALSModuleData_H

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

  void     DropData();

  void     Init(Int_t sm);

  void     RegisterDigit(Int_t AbsId, Int_t isupMod, Float_t iamp, Float_t ix, Float_t iy, Float_t iz);
  void     RegisterCluster(Int_t isupMod, Float_t iamp, Float_t ix, Float_t iy, Float_t iz); 
  void     RegisterHit(Int_t AbsId, Int_t isupMod, Float_t iamp, Float_t ix, Float_t iy, Float_t iz); 

  Int_t    GetNDigits()   const { return fNDigits;   }; 
  Int_t    GetNClusters() const { return fNClusters; }; 
  Int_t    GetNHits()     const { return fNHits;     }; 
  Int_t    GetSmId()      const { return fSmId; };
  Int_t    GetNsm()       const {return fNsm;};
  Int_t    GetNsmf()      const {return fNsmfull;};
  Int_t    GetNsmh()      const {return fNsmhalf;};
  vector< vector<Float_t> > GetDigitBuffer()   { return fDigitArray;   };  
  vector< vector<Float_t> > GetClusterBuffer() { return fClusterArray; };  
  vector< vector<Float_t> > GetHitBuffer()     { return fHitArray;     };  

  void     GetSModuleBigBox(Float_t& bbox0, Float_t& bbox1, Float_t& bbox2) 
  { bbox0 = fSModuleBigBox0; bbox1 = fSModuleBigBox1; bbox2 = fSModuleBigBox2;};
  void     GetSModuleSmallBox(Float_t& bbox0, Float_t& bbox1, Float_t& bbox2) 
  { bbox0 = fSModuleSmallBox0; bbox1 = fSModuleSmallBox1; bbox2 = fSModuleSmallBox2;};
  void     GetSModuleCenter(Float_t& bboxCenter0, Float_t& bboxCenter1, Float_t& bboxCenter2) 
  { bboxCenter0 = fSModuleCenter0; bboxCenter1 = fSModuleCenter1; bboxCenter2 = fSModuleCenter2;};
  Float_t  GetPhiTileSize()   {return fPhiTileSize;};
  Float_t  GetEtaTileSize()   {return fEtaTileSize;};
  TGeoMatrix* GetSModuleMatrix() {return fMatrix;};
  
protected:

   AliEMCALGeometry* fGeom;
   TGeoNode*         fNode;
   Int_t             fSmId;                 // number of the chamber, 0 to 13 
   Int_t             fNsm;
   Int_t             fNsmfull;
   Int_t             fNsmhalf;
   Int_t             fNDigits;                   // number of found digits 
   Int_t             fNClusters;                 // number of found rec points 
   Int_t             fNHits;                     // number of simulation hits 

   Float_t           fPhiTileSize;
   Float_t           fEtaTileSize;

   vector< vector<Float_t> > fHitArray;        //|| hits coordinates, etc.
   vector< vector<Float_t> > fDigitArray;      //|| digits coordinates, etc.
   vector< vector<Float_t> > fClusterArray;    //|| cluster coordinates, etc.

   static Float_t    fSModuleBigBox0;   // sm envelope box
   static Float_t    fSModuleBigBox1;   // sm envelope box
   static Float_t    fSModuleBigBox2;   // sm envelope box
   static Float_t    fSModuleSmallBox0;   // sm envelope box
   static Float_t    fSModuleSmallBox1;   // sm envelope box
   static Float_t    fSModuleSmallBox2;   // sm envelope box
   static Float_t    fSModuleCenter0;   // sm envelope box
   static Float_t    fSModuleCenter1;   // sm envelope box
   static Float_t    fSModuleCenter2;   // sm envelope box

   TGeoMatrix*   fMatrix;
   TGeoHMatrix*  fHMatrix;

private:

  AliEveEMCALSModuleData(const AliEveEMCALSModuleData&);            // Not implemented
  AliEveEMCALSModuleData& operator=(const AliEveEMCALSModuleData&); // Not implemented

  ClassDef(AliEveEMCALSModuleData, 0);     // class with data for one chamber
};

#endif
