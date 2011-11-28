//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTEVECALO_H
#define ALIHLTEVECALO_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTEveCalo.h
/// @author Svein Lindal
/// @brief  Base class for the HLT eve calorimeter display elements

#include "AliHLTEveBase.h"
#include "TString.h"

class AliHLTCaloClusterReader;
class AliESDEvent;
class AliESDCaloCluster;
class TRefArray;
class TEveElementList;
class TEveBoxSet;
class AliHLTHOMERBlockDesc;
class TH1F;


class AliHLTEveCalo : public AliHLTEveBase {

public:
  
  /** Constructor  **/
  AliHLTEveCalo(Int_t nm, TString name);

  /** Destructor **/
 ~AliHLTEveCalo();
  
  void ProcessBlock(AliHLTHOMERBlockDesc * block);

  /** inherited from AliHLTEveBase */
  void UpdateElements();
  
  /** inherited from AliHLTEveBase */
  void ResetElements();

  /// Process esd event
  void ProcessEvent(AliESDEvent * event);


protected :

  /* Process ESD event block */
  void ProcessEsdBlock(AliHLTHOMERBlockDesc * block);
  

  /** Create the elementlist */
  virtual void CreateElementList() = 0;
  
  /** Add clusters boxset to eve display */
  virtual void AddClusters(Float_t * pos, Int_t module, Float_t energy) = 0;

  /** Add digits boxset to eve display */
  virtual void AddDigits(UShort_t fX, UShort_t fZ, Int_t module, Float_t energy) = 0;

  
  virtual void ProcessESDCluster(AliESDCaloCluster * cluster)  = 0;


  /** Process clusters block */
  void ProcessClusters(AliHLTHOMERBlockDesc * block);
  
  /** Process digits block */
  void ProcessDigits(AliHLTHOMERBlockDesc * block);

  /** Process histogram block */
  void ProcessHistogram(AliHLTHOMERBlockDesc * block );

  
  /** Process and draw histograms */
  void AddHistogramsToCanvas(AliHLTHOMERBlockDesc * block, TCanvas * canvas, Int_t &cdCount );  
  

  Int_t GetPadNumber(TString name);  

  virtual Int_t GetClusters (AliESDEvent * event, TRefArray * clusters) = 0;

  TEveBoxSet * fBoxSetDigits;            //Boxset for clusters and digist
  TEveBoxSet * fBoxSetClusters;            //Boxset for clusters and digist

  const Int_t fNModules;          //Number of modules in calorimeter


  TRefArray * fClustersArray;


private:
  
  /** default constructor prohibited */
  AliHLTEveCalo();
  /** copy constructor prohibited */
  AliHLTEveCalo(const AliHLTEveCalo&);
  /** assignment operator prohibited */
  AliHLTEveCalo& operator = (const AliHLTEveCalo &);

  void DrawInvMassHistogram(TH1F * histo);

  TString fName;  //PHOS or EMCAL
  TString * fPadTitles;
  TCanvas * fInvMassCanvas;
  AliHLTCaloClusterReader * fClusterReader;


  ClassDef(AliHLTEveCalo, 0);
};

#endif
