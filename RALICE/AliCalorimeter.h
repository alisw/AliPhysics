#ifndef ALICALORIMETER_H
#define ALICALORIMETER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

#include <math.h>
 
#include "TNamed.h"
#include "TObjArray.h"
#include "TH2.h"
#include "TString.h"

#include "AliDevice.h"
#include "AliObjMatrix.h"
#include "AliCalmodule.h"
#include "AliCalcluster.h"
#include "AliPositionObj.h"
#include "AliAttribObj.h"
 
class AliCalorimeter : public AliDevice
{
 public:
  AliCalorimeter();                                      // Default constructor
  AliCalorimeter(Int_t nrow,Int_t ncol);                 // Create a calorimeter matrix
  virtual ~AliCalorimeter();                             // Destructor
  AliCalorimeter(const AliCalorimeter& c);               // Copy constructor
  virtual TObject* Clone(const char* name="") const;     // Make a deep copy and provide pointer of the copy
  Int_t GetNrows();                                      // Return number of rows of the matrix
  Int_t GetNcolumns();                                   // Return number of columns of the matrix
  using AliDevice::SetSignal;
  void SetSignal(Int_t row,Int_t col,Float_t s);         // Set signal for a certain module
  using AliDevice::AddSignal;
  void AddSignal(Int_t row,Int_t col,Float_t s);         // Add signal to a certain module
  void AddSignal(AliCalmodule* m);                       // Add module signal to current calorimeter
  void Reset(Int_t row,Int_t col);                       // Reset signal for a certain module
  virtual void Reset(Int_t mode=0);                      // Reset the complete calorimeter
  using AliDevice::GetSignal;
  virtual Float_t GetSignal(Int_t row,Int_t col=0) { return GetSignal(row,col,0); }
  Float_t GetSignal(Int_t row,Int_t col,Int_t mode);     // Provide signal of a certain module
  Int_t GetNsignals() const;                             // Return number of modules with a signal
  void Group(Int_t n=1,Int_t mode=1);                    // Group modules into clusters (n rings)
  Int_t GetNclusters() const;                            // Return number of clusters
  Float_t GetClusteredSignal(Int_t row,Int_t col);       // Provide module signal after clustering
  AliCalcluster* GetCluster(Int_t j) const;              // Access to cluster number j
  AliCalmodule* GetModule(Int_t j) const;                // Access to 'fired' module number j
  AliCalmodule* GetModule(Int_t row,Int_t col);          // Access to module at (row,col)
  using AliDevice::SetEdgeOn;
  void SetEdgeOn(Int_t row,Int_t col);                   // Indicate module as 'edge module'
  using AliDevice::SetEdgeOff;
  void SetEdgeOff(Int_t row,Int_t col);                  // Indicate module as 'non-edge module'
  using AliDevice::GetEdgeValue;
  Int_t GetEdgeValue(Int_t row,Int_t col);               // Provide the edge flag of a module
  using AliDevice::SetDead;
  void SetDead(Int_t row,Int_t col);                     // Indicate module as 'dead module'
  using AliDevice::SetAlive;
  void SetAlive(Int_t row,Int_t col);                    // Indicate module as 'active module'
  using AliDevice::GetDeadValue;
  Int_t GetDeadValue(Int_t row,Int_t col);               // Provide the dead flag of a module
  using AliDevice::SetGain;
  void SetGain(Int_t row,Int_t col,Float_t g);           // Set the gain value for a module
  using AliDevice::SetOffset;
  void SetOffset(Int_t row,Int_t col,Float_t o);         // Set the offset value for a module
  using AliDevice::GetGain;
  Float_t GetGain(Int_t row,Int_t col);                  // Provide the gain value of a module
  using AliDevice::GetGainFlag;
  Int_t GetGainFlag(Int_t row,Int_t col);                // Provide the gain flag value of a module
  using AliDevice::GetOffset;
  Float_t GetOffset(Int_t row,Int_t col);                // Provide the offset value of a module
  using AliDevice::GetOffsetFlag;
  Int_t GetOffsetFlag(Int_t row,Int_t col);              // Provide the offset flag value of a module
  using AliDevice::SetPosition;
  void SetPosition(Int_t row,Int_t col,Float_t* r,TString f); // Set module position
  void SetPosition(Int_t row,Int_t col,Ali3Vector& r);   // Set module position
  using AliDevice::GetPosition;
  void GetPosition(Int_t row,Int_t col,Float_t* r,TString f); // Return module position
  AliPosition* GetPosition(Int_t row,Int_t col);         // Access to module position
  TH2F* DrawModules(Float_t thresh=0.,Int_t mode=0);     // Lego plot of module (corr.) signals above threshold
  TH2F* DrawClusters(Float_t thresh=0.);                 // Lego plot of cluster signals above threshold
  void AddVetoSignal(AliSignal& s);                      // Associate (extrapolated) signal
  void AddVetoSignal(AliSignal* s) { AddVetoSignal(*s); }
  AliSignal* GetVetoSignal(Int_t j) const;               // Access to veto signal number j
  Int_t GetNvetos() const;                               // Provide the number of veto signals
  void SetMatrixSwapMode(Int_t swap=1);                  // Set the swapmode for the storage of the matrices
  Int_t GetMatrixSwapMode() const;                       // Provide the swapmode for the storage of the matrices
 
 protected:
  Int_t fNrows;                              // The number of rows
  Int_t fNcolumns;                           // The number of columns
  AliObjMatrix* fMatrix;                     //! Matrix lookup table of module pointers
  Int_t fSwap;                               // The swapmode for the module and position matrices
  void SortM();                              // Order the modules with decreasing signal (matrix search)
  void SortA();                              // Order the modules with decreasing signal (fired array search)
  TObjArray* fClusters;                      // The array of clusters
  void AddRing(Int_t row,Int_t col,Int_t n); // add signals of n rings around cluster center
  void Ungroup();                            // Restore module matrix as before clustering
  TH2F* fHmodules;                           //! The module 2-D histogram for event display
  TH2F* fHclusters;                          //! The cluster 2-D histogram for event display
  TObjArray* fVetos;                         // The array of associated (extrapolated) veto signals
  TObjArray* fAttributes;                    //! Matrix dbase with module attributes (e.g. gain, offset etc...)
  AliObjMatrix* fPositions;                  //! Matrix dbase of module position pointers
  void LoadMatrix();                         // Loading of matrix lookup table from the linear hit array
 
 ClassDef(AliCalorimeter,11) // Description of a modular calorimeter system.
};
#endif
