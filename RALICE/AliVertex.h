#ifndef ALIVERTEX_H
#define ALIVERTEX_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <iostream.h>
#include <math.h>
 
#include "TObject.h"
#include "TObjArray.h"
 
#include "AliJet.h"
#include "AliPosition.h"
 
class AliVertex : public AliJet,public AliPosition
{
 public:
  AliVertex();                            // Default constructor
  AliVertex(Int_t n);                     // Create a vertex to hold initially n tracks
  ~AliVertex();                           // Default destructor
  void Reset();                           // Reset all values
  void Add(AliJet& j);                    // Add a jet of tracks to the vertex
  void Add(AliVertex& v,Int_t connect=1); // Add (and connect) a (sec.) vertex to the current vertex
  void Add(AliJet* j)    { Add(*j); }
  void Add(AliVertex* v,Int_t connect=1) { Add(*v,connect); }
  void Info(TString f="car");             // Print the vertex info within coordinate frame f
  void List(TString f="car");             // Print vertex prim. track information for coord. frame f
  void ListAll(TString f="car");          // Print prim. + sec. vertex full track info for coord. frame f
  Int_t GetNvertices();                   // Return the number of (secondary) vertices
  AliVertex* GetVertex(Int_t i);          // Provide i-th (secondary) vertex
  void SetNvmax(Int_t n=2);               // Set the initial max. number of (secondary) vertices

 protected:
  Int_t fNvmax;         // The maximum number of (secondary) vertices
  Int_t fNvtx;          // The number of (secondary) vertices
  TObjArray* fVertices; // Array to hold the pointers to the (secondary) vertices
  TObjArray* fConnects; // Array to hold the pointers to the auto-generated connecting tracks

 private:
  void Dump(AliVertex* v,Int_t n,TString f); // Recursively print all sec. vertices
 
 ClassDef(AliVertex,1) // Creation and investigation of an AliVertex.
};
#endif
