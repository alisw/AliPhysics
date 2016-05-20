/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// Interface for all the histogram container classes in JCORRAN analysis

//===========================================================
// AliJHistogramInterface.h
//
//   J
//===========================================================

#ifndef ALIJHISTOGRAMINTERFACE_H
#define ALIJHISTOGRAMINTERFACE_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include <TLorentzVector.h>

#include "AliJHistManager.h"

class AliJCard;
class AliJBaseTrack;
class AliJPhoton;
class AliJTrack;

using namespace std;

/*
 * Stream operator for nicer histogram creation
 */
inline ostream &operator << (ostream &out_file, const TLorentzVector &Vec)
{
  out_file<<"Px="<<Vec.Px()<<" Py="<<Vec.Py()<<" Pz="<<Vec.Pz()<<" E="<<Vec.E()<<" M="<<Vec.M()<<endl;
  out_file<<"Theta="<<Vec.Theta()<<" Phi="<<Vec.Phi()<<" p="<<Vec.Rho()<<endl;
  return(out_file);
}

class AliJHistogramInterface {
  
public:
  AliJHistogramInterface(); //default constructor
  AliJHistogramInterface(AliJCard* cardP); //constructor
  AliJHistogramInterface(const AliJHistogramInterface& obj); // copy constructor
  virtual ~AliJHistogramInterface();    //destructor
  
  // Inheriting classes need to implement this function
  // It should create the histograms below
  // These histograms can be created using the CreateDataManagerHistograms()
  virtual void CreateEventTrackHistos() = 0; // Creation of basic event histograms
  
  // Histograms filled in AliJDataManager.cxx
  AliJTH1D fhTrackSelection; // checking bit convention
  AliJTH1D fhEvents; // event count in centrality classes
  AliJTH1D fhEventTrigger; // event count per trigger bits
  AliJTH1D fhZVertRaw; // z-vertex distribution
  AliJTH2D fhVertexZTriggVtx; // number of vertices and z-vertex distribution
  
  AliJHistManager *fHMG; // Histogram manager
  
protected:
  void CreateDataManagerHistograms(); // Creation of histograms needed in AliJDataManager.cxx
  
  AliJCard *fCard;        // JCrad containing binning information etc.
  
  AliJBin fCentBin;       // Bin of Centrality
  AliJBin fVtxBin;        // Bin of Z-vertex Bin
  AliJBin fPTtBin;        // Bin of pT Trigged
  AliJBin fPTaBin;        // Bin of pT Associated
  AliJBin fXEBin;         // Bin of XE
  AliJBin fKLongBin;      // Bin of Klong
  AliJBin fRGapBin;       // Bin of R-gap
  AliJBin fEtaGapBin;     // Bin of Eta-gap
  AliJBin fPhiGapBin;     // Bin of Phi-gap
  AliJBin fMassBin;       // Bin of Mass
  AliJBin fTypBin;        // Bin of type ( data, mixed ), but is being used for any 2 dimension
  AliJBin fTypBin3;       // Bin of type3, is being used for any 3 dimension
  AliJBin fPairPtBin;     // Bin of pT pair

};

#endif






















