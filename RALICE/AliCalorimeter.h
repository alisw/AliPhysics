#ifndef ALICALORIMETER_H
#define ALICALORIMETER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////
// Class AliCalorimeter
// Description of a modular calorimeter system.
// A matrix geometry is used in which a module is identified by (row,col).
// Note : First module is identified as (1,1).
//
// This is the way to define and enter signals into a calorimeter :
//
//   AliCalorimeter cal(10,15);  // Calorimeter of 10x15 modules
//                               // All module signals set to 0.
//   cal.AddSignal(5,7,85.4);
//   cal.AddSignal(5,7,25.9);
//   cal.AddSignal(3,5,1000);
//   cal.SetSignal(5,7,10.3);
//   cal.Reset(3,5);             // Reset module (3,5) as being 'not fired'
//                               // All module data are re-initialised.
//   cal.SetEdgeOn(1,1);         // Declare module (1,1) as an 'edge module'
//   cal.SetDead(8,3);
//   cal.SetGain(2,8,3.2);
//
//   Float_t vec[3]={6,1,20};
//   cal.SetPosition(2,8,vec,"car");
//
//   Float_t loc[3]={-1,12,3};
//   cal.AddVetoSignal(loc,"car"); // Associate (extrapolated) position as a veto
//
//   cal.Group(2);      // Group 'fired' modules into clusters
//                      // Perform grouping over 2 rings around the center
//   cal.Reset();       // Reset the complete calorimeter
//                      // Normally to prepare for the next event data
//                      // Note : Module gain, edge and dead flags remain
//
//--- NvE 13-jun-1997 UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////
 
#include <iostream.h>
#include <math.h>
 
#include "TObject.h"
#include "TObjArray.h"
#include "TH2.h"
#include "TString.h"

#include "AliDetector.h"

#include "AliCalmodule.h"
#include "AliCalcluster.h"
#include "AliSignal.h"
 
class AliCalorimeter : public AliDetector
{
 public:
  AliCalorimeter();                                // Default constructor
  AliCalorimeter(Int_t nrow,Int_t ncol);           // Create a calorimeter matrix
  ~AliCalorimeter();                               // Destructor
  Int_t GetNrows();                                // Return number of rows of the matrix
  Int_t GetNcolumns();                             // Return number of columns of the matrix
  void SetSignal(Int_t row,Int_t col,Float_t s);   // Set signal for a certain module
  void AddSignal(Int_t row,Int_t col,Float_t s);   // Add signal to a certain module
  void Reset(Int_t row,Int_t col);                 // Reset signal for a certain module
  void Reset();                                    // Reset the complete calorimeter
  Float_t GetSignal(Int_t row,Int_t col);          // Provide signal of a certain module
  Int_t GetNsignals();                             // Return number of modules with a signal
  void Group(Int_t n);                             // Group modules into clusters (n rings)
  Int_t GetNclusters();                            // Return number of clusters
  Float_t GetClusteredSignal(Int_t row,Int_t col); // Provide module signal after clustering
  AliCalcluster* GetCluster(Int_t j);              // Access to cluster number j
  AliCalmodule* GetModule(Int_t j);                // Access to 'fired' module number j
  void SetEdgeOn(Int_t row,Int_t col);             // Indicate module as 'edge module'
  void SetEdgeOff(Int_t row,Int_t col);            // Indicate module as 'non-edge module'
  Int_t GetEdgeValue(Int_t row,Int_t col);         // Provide the edge flag of a module
  void SetDead(Int_t row,Int_t col);               // Indicate module as 'dead module'
  void SetAlive(Int_t row,Int_t col);              // Indicate module as 'active module'
  Int_t GetDeadValue(Int_t row,Int_t col);         // Provide the dead flag of a module
  void SetGain(Int_t row,Int_t col,Float_t g);     // Set the gain value for a module
  Float_t GetGain(Int_t row,Int_t col);            // Provide the gain value of a module
  void SetPosition(Int_t row,Int_t col,Float_t* r,TString f); // Set module position
  void GetPosition(Int_t row,Int_t col,Float_t* r,TString f); // Return module position
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
  AliCalmodule** fMatrix;                    //! The matrix of modules for internal use
  void Sortm(AliCalmodule*);                 // Order the modules with decreasing signal
  TObjArray* fClusters;                      // The array of clusters
  void AddRing(Int_t row,Int_t col,Int_t n); // add signals of n rings around cluster center
  TObjArray* fModules;                       // The array of modules for output
  void LoadMatrix();                         // Load calorimeter matrix data from input
  void Ungroup();                            // Restore module matrix as before clustering
  TH2F* fHmodules;                           //! The module 2-D histogram
  TH2F* fHclusters;                          //! The cluster 2-D histogram
  Int_t fNvetos;                             // The number of associated veto signals
  TObjArray* fVetos;                         // The array of associated (extrapolated) veto signals
 
 ClassDef(AliCalorimeter,1) // Class definition to enable ROOT I/O
};
#endif
