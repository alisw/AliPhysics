#ifndef ALICALMODULE_H
#define ALICALMODULE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

#include "AliSignal.h"
 
class AliCalmodule : public AliSignal
{
 public:
  AliCalmodule();                                          // Default constructor
  virtual ~AliCalmodule();                                 // Default destructor
  AliCalmodule(AliCalmodule& m);                           // Copy constructor
  AliCalmodule(Int_t row,Int_t col,Float_t sig);           // Create a module and initialise data
  virtual void SetSignal(Int_t row,Int_t col,Float_t sig); // Set or change data for certain module
  virtual void AddSignal(Int_t row,Int_t col,Float_t sig); // Add signal to a certain module
  void SetRow(Int_t i);                                    // Set the row number of the module
  void SetColumn(Int_t i);                                 // Set the column number of the module
  Int_t GetRow();                                          // Return the row number of the module
  Int_t GetColumn();                                       // Return the column number of the module
  void SetClusteredSignal(Float_t val);                    // Set the signal of the module after clustering
  Float_t GetClusteredSignal();                            // Return module signal after clustering
  virtual AliCalmodule* MakeCopy(AliCalmodule& m);         // Make a deep copy and provide its pointer
 
 protected:
  Int_t fRow;        // The current row number
  Int_t fCol;        // The current column number
  Float_t fSigc;     // The signal after clustering
 
 ClassDef(AliCalmodule,4) // Description of a module in a calorimeter system.
};
#endif
