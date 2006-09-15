/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
//
/// \ingroup sim
/// \class AliMUONTest
/// \brief Class with functions for testing new segmentation
///
/// Author: Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_TEST_H
#define ALI_MUON_TEST_H

#include <TObject.h>

class AliMUONGeometryTransformer;
class AliMUONSegmentation;
class AliMUONGeometrySegmentation;

class TCanvas;
class TString;

class AliMUONTest : public  TObject 
{
  public:
    AliMUONTest(const TString& option);
    AliMUONTest();
    virtual ~AliMUONTest();
   

    // methods for printing pads
    //							  
    void PrintPadsForAll() const;
    void PrintPadsForSegmentation(Int_t moduleId, Int_t cath) const;
    void PrintPadsForDetElement(Int_t detElemId, Int_t cath) const;
    void PrintPad(Int_t& counter, 
                  Int_t detElemId, Int_t ix, Int_t iy,
                  AliMUONGeometrySegmentation* segmentation) const;

    // methods for drawing pads
    //							  
    void DrawPadsForAll() const;
    void DrawPadsForSegmentation(Int_t moduleId, Int_t cath) const;
    void DrawPadsForDetElement(Int_t detElemId, Int_t cath) const;
    void DrawPad(Int_t& counter, 
                 Int_t detElemId, Int_t ix, Int_t iy,
                 AliMUONGeometrySegmentation* segmentation) const;

    // other tests
    //
    void DetElemTransforms() const;

  private:  
    AliMUONTest(const AliMUONTest& rhs);
    AliMUONTest& operator = (const AliMUONTest& rhs);

    // methods
    void BuildWithMUON(const TString& configMacro);
    void BuildWithoutMUON(const TString& option);

    void PrintPad(Int_t& counter, 
                  Int_t detElemId, Int_t ix, Int_t iy,
                  AliMUONGeometrySegmentation* segmentation);
    void DrawPad(Int_t& counter, 
                  Int_t detElemId, Int_t ix, Int_t iy,
                  AliMUONGeometrySegmentation* segmentation);
    // data members
    const AliMUONGeometryTransformer* fkTransformer; ///< Geometry parametrisation
    AliMUONSegmentation*  fSegmentation;  ///< Segmentation
    TCanvas*              fCanvas;        ///< The canvas for drawing				       

    ClassDef(AliMUONTest,0)  // MUON class for tests
};

#endif //ALI_MUON_TEST_H

