#ifndef ALIVERTEX_H
#define ALIVERTEX_H
///////////////////////////////////////////////////////////////////////////
// Class AliVertex
// Creation and investigation of an AliVertex.
// An AliVertex can be constructed by adding AliTracks and/or AliJets.
//
// Note : Also (secondary) vertices can be added to a vertex.
//
// Coding example to make 3 vertices v1, v2 and v3.
// ------------------------------------------------
// v1 contains the tracks 1,2,3 and 4
// v2 contains the tracks 5,6 and 7
// v3 contains the jets 1 and 2
//
//        AliTrack t1,t2,t3,t4,t5,t6,t7;
//         ...
//         ... // code to fill the track data
//         ...
//
//        AliJet j1,j2;
//         ...
//         ... // code to fill the jet data
//         ...
//
//        AliVertex v1(5);
//
//        v1.Add(t1);
//        v1.Add(t2);
//        v1.Add(t3);
//        v1.Add(t4);
//
//        Float_t r1[3]={2.4,0.1,-8.5};
//        v1.SetPosition(r1,"car");
//
//        AliVertex v2(2);
//        v2.Add(t5);
//        v2.Add(t6);
//        v2.Add(t7);
//
//        Float_t r2[3]={1.6,-3.2,5.7};
//        v2.SetPosition(r2,"car");
//
//        AliVertex v3;
//
//        v3.Add(j1);
//        v3.Add(j2);
//
//        Float_t r3[3]={6.2,4.8,1.3};
//        v3.SetPosition(r3,"car");
//
//        v1.Info("sph");
//        v2.ListAll();
//        v3.List("cyl");
//
//        Float_t e1=v1.GetEnergy();
//        Ali3Vector p1=v1.Get3Momentum();
//        Float_t loc[3];
//        v1.GetPosition(loc,"sph");
//        AliPosition r=v2.GetPosition();
//        r.Info(); 
//        Int_t nt=v2.GetNtracks();
//        AliTrack* tv=v2.GetTrack(1); // Access track number 1 of Vertex v2
//
// Specify the vertices v2 and v3 as secondary vertices of v1
//
//        v1.Add(v2);
//        v1.Add(v3);
//
//        v1.List();
//
//        Int_t nv=v1.GetNvtx();
//        AliVertex* vx=v1.GetVertex(1); // Access 1st secondary vertex of v1
//        Float_t e=vx->GetEnergy();
//
//        Float_t M=v1.GetInvmass(); 
//
// Reconstruct Vertex v1 from scratch
//
//        v1.Reset();
//        v1.SetNvmax(25); // Increase initial no. of sec. vertices
//        v1.Add(t3);
//        v1.Add(t7);
//        v1.Add(j2);
//        Float_t pos[3]={7,9,4};
//        v1.SetPosition(pos,"car");
//
// Note : All quantities are in GeV, GeV/c or GeV/c**2
//
//--- NvE 04-apr-1998 UU-SAP Utrecht
//--- Modified : NvE 08-apr-1999 UU-SAP Utrecht to inherit from AliJet
///////////////////////////////////////////////////////////////////////////
 
#include <iostream.h>
#include <math.h>
 
#include "TObject.h"
#include "TObjArray.h"
 
#include "AliJet.h"
#include "AliPosition.h"
 
class AliVertex : public AliJet,public AliPosition
{
 public:
  AliVertex();                        // Default constructor
  AliVertex(Int_t n);                 // Create a vertex to hold initially n tracks
  ~AliVertex();                       // Default destructor
  void Reset();                       // Reset all values
  void Add(AliJet& j);                // Add a jet of tracks to the vertex
  void Add(AliVertex& v);             // Add a (secondary) vertex to the current vertex
  void Add(AliJet* j)    { Add(*j); }
  void Add(AliVertex* v) { Add(*v); }
  void Info(TString f="car");         // Print the vertex info within coordinate frame f
  void List(TString f="car");         // Print vertex prim. track information for coord. frame f
  void ListAll(TString f="car");      // Print prim. + sec. vertex full track info for coord. frame f
  Int_t GetNvertices();               // Return the number of (secondary) vertices
  AliVertex* GetVertex(Int_t i);      // Provide i-th (secondary) vertex
  void SetNvmax(Int_t n=2);           // Set the initial max. number of (secondary) vertices

 protected:
  Int_t fNvmax;         // The maximum number of (secondary) vertices
  Int_t fNvtx;          // The number of (secondary) vertices
  TObjArray* fVertices; // Array to hold the pointers to the (secondary) vertices

 private:
  void Dump(AliVertex* v,Int_t n,TString f); // Recursively print all sec. vertices
 
 ClassDef(AliVertex,1) // Class definition to enable ROOT I/O
};
#endif
