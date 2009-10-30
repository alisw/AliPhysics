#ifndef ALIHLTGLOBALVERTEXERCOMPONENT_H
#define ALIHLTGLOBALVERTEXEROMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTGlobalVertexerComponent.h
    @author Sergey Gorbunov
    @brief  Component for monitor V0 physics
*/

#include "AliHLTProcessor.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
class TH1F;
class TH2F;
class AliESDtrack;
class AliESDVertex;
class AliTracker;
class AliESDtrack;
class AliESDEvent;

/**
 * @class AliHLTTPCGlobalVertexerComponent
 * Component for reconstruct primary vertex and V0's
 */
class AliHLTGlobalVertexerComponent : public AliHLTProcessor
{
public:
  /** default constructor */
  AliHLTGlobalVertexerComponent();
  /** destructor */
  virtual ~AliHLTGlobalVertexerComponent();

  // Public functions to implement AliHLTComponent's interface.
  // These functions are required for the registration process

  /** interface function, see AliHLTComponent for description */
  const char* GetComponentID();
  /** interface function, see AliHLTComponent for description */
  void GetInputDataTypes(AliHLTComponentDataTypeList& list);
  /** interface function, see AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();
  /** interface function, see AliHLTComponent for description */
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);

  /** interface function, see AliHLTComponent for description */
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  /** interface function, see AliHLTComponent for description */
  AliHLTComponent* Spawn();

protected:

  // Protected functions to implement AliHLTComponent's interface.
  // These functions provide initialization as well as the actual processing
  // capabilities of the component. 

  /** interface function, see AliHLTComponent for description */
  int DoInit( int argc, const char** argv );
  /** interface function, see AliHLTComponent for description */
  int DoDeinit();
  /** interface function, see AliHLTComponent for description */
  int DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& trigData );

  int Reconfigure(const char* cdbEntry, const char* chainId);

  using AliHLTProcessor::DoEvent;
  
private:
  /** copy constructor prohibited */
  AliHLTGlobalVertexerComponent(const AliHLTGlobalVertexerComponent&);
  /** assignment operator prohibited */
  AliHLTGlobalVertexerComponent& operator=(const AliHLTGlobalVertexerComponent&);
  /**
   * Configure the component.
   * Parse a string for the configuration arguments and set the component
   * properties.
   */ 
  int Configure(const char* arguments);
  
  class AliESDTrackInfo
    {
    public:
      AliESDTrackInfo(): fParticle(),fPrimDeviation(0),fPrimUsedFlag(0),fOK(0){}
      
      AliKFParticle fParticle; //* assigned KFParticle
      Double_t fPrimDeviation; //* deviation from the primary vertex
      Bool_t fPrimUsedFlag;    //* flag shows that the particle was used for primary vertex fit
      Bool_t fOK;              //* is the track good enough
    };
  void SetESD( AliESDEvent *event );
  void FindPrimaryVertex();
  void FindV0s();


  AliESDEvent *fESD; // pointer to esd event
  AliESDTrackInfo *fTrackInfos; // information about esd tracks
  AliKFVertex fPrimaryVtx; // reconstructed KF primary vertex

  TH2F *fHistPrimVertexXY; // primary vertex distribution in XY;
  TH2F *fHistPrimVertexZX; // primary vertex distribution in ZX;
  TH2F *fHistPrimVertexZY; // primary vertex distribution in ZY;
  Int_t fNEvents; // n of processed events

  Bool_t fPlotHistograms;// flag to produce histogramms
  Bool_t fFillVtxConstrainedTracks; // flag to store vertex constrained track parameters
  Double_t fConstrainedTrackDeviation; // deviation of a track from prim.vtx <=cut 
  Double_t fV0DaughterPrimDeviation; // v0: daughters deviation from prim vertex >= cut
  Double_t fV0PrimDeviation; // v0: v0 deviation from prim vertex <= cut
  Double_t fV0Chi; // v0: v0 sqrt(chi^2/NDF) <= cut
  Double_t fV0DecayLengthInSigmas; // v0: v0 decay length/sigma_length >= cut


  ClassDef(AliHLTGlobalVertexerComponent, 0);

};
#endif
