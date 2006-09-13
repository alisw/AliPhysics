#ifndef ALIPHOSONLINEMONITOR_H
#define ALIPHOSONLINEMONITOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Simple class for online mnitoring     
//                  
//*-- Author: D.Peressouno (RRC KI)


// --- ROOT system ---
#include "TDialogCanvas.h" 
class TList ;
class TString ;
class TClonesArray;
// --- Standard library ---

// --- AliRoot header files ---
class AliPHOSGeometry ;
class AliPHOSConTableDB ;

class AliPHOSOnlineMonitor: public TDialogCanvas{

public:
  AliPHOSOnlineMonitor() ;          // ctor
  AliPHOSOnlineMonitor(const char * dataname) ;          // ctor
  virtual ~AliPHOSOnlineMonitor() ; // dtor

  virtual void Apply(const char *action=""); //Overloaded function to 
                                             //handle butons
  void SetInputFile(const char * filename = "run1.dat") ;    //Open new input file
  void DrawPedestals();    //Scans pedestals in current file
  void DrawSpectrum(const char * opt);  //Plot spectrum of all registered particles"
  void DrawMinv() ;        //Plot invariant mass distribution
  void DrawTriggers() ;        //Plot invariant mass distribution
  void DrawEdep(Int_t mod,const char * opt) ;//Plots energy deposition per crystal for module
  void Go() ;
  void Clean() ;           //Cleans all histograms
  void Reset() ;           //Removes all canvas and histograms

  void WriteHistograms(const char * filename = "onlineout.root") ;

  void SetConTableDB(const char * filename = "ConTableDB.root") ;

 private:
  void MakeButtons(void) ; //Function to make menu and buttons

  void ScanDigits() ;
  void ScanTrigger(Int_t trig) ;
  void ScanPedestals(TClonesArray * digits) ;
  void ScanEdep(TClonesArray * digits) ;
  void ScanRecon(TClonesArray * recParticles) ;

private:
  //They call Fatal, but they are private, user will have compile time error instead
  //of run-time error. Fatal - because it's not clear, should I copy canvases, 
  //hists etc.
  AliPHOSOnlineMonitor(const AliPHOSOnlineMonitor &);
  AliPHOSOnlineMonitor & operator = (const AliPHOSOnlineMonitor &);

  Bool_t fScanPed ;     //should we analyse pedestal events
  Bool_t fScanSig;      //should we analyse signal events
  Bool_t fReconstruct ; //should we reconstruct events

  Int_t fNevents ;      //Number of events processed so far 
  Int_t fNUpdate ;      //Number of events to update histo 

  TList * fCanvasList ;  //!
  TList * fHistosList ;  //!
  
  TString fInputFile ;

  AliPHOSGeometry   * fGeom ;         //! 
  AliPHOSConTableDB * fcdb ;          //!

  ClassDef(AliPHOSOnlineMonitor,1)  //PHOS online monitoring  

};

#endif // AliPHOSONLINEMONITOR_H
