#ifndef ALICALORIMETER_H
#define ALICALORIMETER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <iostream.h>
#include <math.h>
 
#include "TObject.h"
#include "TObjArray.h"
#include "TH2.h"
#include "TString.h"
#include "TMatrix.h"

#include "AliCalmodule.h"
#include "AliCalcluster.h"
#include "AliSignal.h"
 
class AliCalorimeter : public TObject
{
 public:
  AliCalorimeter();                                // Default constructor
  AliCalorimeter(Int_t nrow,Int_t ncol);           // Create a calorimeter matrix
  ~AliCalorimeter();                               // Destructor
  Int_t GetNrows();                                // Return number of rows of the matrix
  Int_t GetNcolumns();                             // Return number of columns of the matrix
  void SetSignal(Int_t row,Int_t col,Float_t s);   // Set signal for a certain module
  void AddSignal(Int_t row,Int_t col,Float_t s);   // Add signal to a certain module
  void AddSignal(AliCalmodule* m);                 // Add module signal to current calorimeter
  void Reset(Int_t row,Int_t col);                 // Reset signal for a certain module
  void Reset();                                    // Reset the complete calorimeter
  Float_t GetSignal(Int_t row,Int_t col);          // Provide signal of a certain module
  Int_t GetNsignals();                             // Return number of modules with a signal
  void Group(Int_t n);                             // Group modules into clusters (n rings)
  Int_t GetNclusters();                            // Return number of clusters
  Float_t GetClusteredSignal(Int_t row,Int_t col); // Provide module signal after clustering
  AliCalcluster* GetCluster(Int_t j);              // Access to cluster number j
  AliCalmodule* GetModule(Int_t j);                // Access to 'fired' module number j
  AliCalmodule* GetModule(Int_t row,Int_t col);    // Access to module at (row,col)
  void SetEdgeOn(Int_t row,Int_t col);             // Indicate module as 'edge module'
  void SetEdgeOff(Int_t row,Int_t col);            // Indicate module as 'non-edge module'
  Int_t GetEdgeValue(Int_t row,Int_t col);         // Provide the edge flag of a module
  void SetDead(Int_t row,Int_t col);               // Indicate module as 'dead module'
  void SetAlive(Int_t row,Int_t col);              // Indicate module as 'active module'
  Int_t GetDeadValue(Int_t row,Int_t col);         // Provide the dead flag of a module
  void SetGain(Int_t row,Int_t col,Float_t g);     // Set the gain value for a module
  Float_t GetGain(Int_t row,Int_t col);            // Provide the gain value of a module
  void SetPosition(Int_t row,Int_t col,Float_t* r,TString f); // Set module position
  void SetPosition(Int_t row,Int_t col,Ali3Vector& r); // Set module position
  void GetPosition(Int_t row,Int_t col,Float_t* r,TString f); // Return module position
  AliPosition* GetPosition(Int_t row,Int_t col);   // Access to module position
  TH2F* DrawModules();                             // Draw lego plot of module signals
  TH2F* DrawClusters();                            // Draw lego plot of cluster signals
  void AddVetoSignal(Float_t* r,TString f,Float_t s=0); // Associate (extrapolated) signal
  AliSignal* GetVetoSignal(Int_t j);               // Access to veto signal number j
  Int_t GetNvetos();                               // Provide the number of veto signals
 
 protected:
  Int_t fNrows;                              // The number of rows
  Int_t fNcolumns;                           // The number of columns
  Int_t fNsignals;                           // The number of modules with a signal
  Int_t fNclusters;                          // The number of clusters
  AliCalmodule ***fMatrix;                   //! The matrix of module pointers for internal use
  void Sortm(AliCalmodule** a,Int_t& n);     // Order the modules with decreasing signal
  TObjArray* fClusters;                      // The array of clusters
  void AddRing(Int_t row,Int_t col,Int_t n); // add signals of n rings around cluster center
  TObjArray* fModules;                       // The array of modules for output
  void LoadMatrix();                         // Load calorimeter matrix data from input
  void Ungroup();                            // Restore module matrix as before clustering
  TH2F* fHmodules;                           //! The module 2-D histogram for event display
  TH2F* fHclusters;                          //! The cluster 2-D histogram for event display
  Int_t fNvetos;                             // The number of associated veto signals
  TObjArray* fVetos;                         // The array of associated (extrapolated) veto signals
  TMatrix* fAttributes;                      // Matrix with module attributes (dead+10*edge)
  TMatrix* fGains;                           // Matrix with module gains
  AliPosition ***fPositions;                 //! Matrix of module position pointers for internal use
 
 ClassDef(AliCalorimeter,1) // Description of a modular calorimeter system.
};
#endif
