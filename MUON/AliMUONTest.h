/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
//
// Class AliMUONTest
// -----------------
// Class with functions for testing
//
// Author: Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_TEST_H
#define ALI_MUON_TEST_H

#include <TObject.h>

class TCanvas;
class AliMUONGeometrySegmentation;

enum AliMUONTests {
  kPrintPads,
  kDrawPads
};  

class AliMUONTest : public  TObject 
{
  public:
    AliMUONTest(const TString& configMacro);
    AliMUONTest();
    virtual ~AliMUONTest();
   
    // create segmentation
    AliMUONGeometrySegmentation* CreateSegmentation(
                                       Int_t chamberId, Int_t cath);
                                                          
    // other tests
    //
    void DetElemTransforms();
    void PrintPadPositionsOld();


    // selected tests
    //							  
    void ForWhole(AliMUONTests test);
    void ForSegmentation(
                  AliMUONTests test,
                  AliMUONGeometrySegmentation* segmentation);
    void ForDetElement(
                  AliMUONTests test,
                  Int_t detElemId,
                  AliMUONGeometrySegmentation* segmentation);
    void Before(AliMUONTests test);
    void After(AliMUONTests test);
 
    // tests per pad
    //
    void PrintPad(Int_t& counter, 
                  Int_t detElemId, Int_t ix, Int_t iy,
                  AliMUONGeometrySegmentation* segmentation);
    void DrawPad(Int_t& counter, 
                  Int_t detElemId, Int_t ix, Int_t iy,
                  AliMUONGeometrySegmentation* segmentation);
 

    void DrawSegmentation(AliMUONGeometrySegmentation *seg);
             // TBR			  
			  

  protected:
    AliMUONTest(const AliMUONTest& rhs);
    AliMUONTest& operator = (const AliMUONTest& rhs);

  private:
    AliMUONGeometrySegmentation* CreateSt1Segmentation(
                                       Int_t chamberId, Int_t cathod);
    AliMUONGeometrySegmentation* CreateSlatSegmentation(
                                       Int_t chamberId, Int_t cathod);
				       
    // data members
    TCanvas* fCanvas; // The canvas for drawing				       

    ClassDef(AliMUONTest,0)  // MUON class for tests
};

#endif //ALI_MUON_TEST_H

