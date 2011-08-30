/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#ifndef ALIANALYSISTASKEMCALFASTOR2TRIGGER_H
#define ALIANALYSISTASKEMCALFASTOR2TRIGGER_H

class TList;
class TTree;
class AliEMCALFastORPatch;
class AliCaloCalibPedestal;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEMCALFastOR2Trigger : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskEMCALFastOR2Trigger();
  AliAnalysisTaskEMCALFastOR2Trigger(const char *name);
  virtual ~AliAnalysisTaskEMCALFastOR2Trigger();

public:  
  virtual void     UserCreateOutputObjects();
  virtual void     UserExec(Option_t *option);
  
  Bool_t                  GetCheckDeadClusters()                  const { return fCheckDeadClusters;    }
  AliCaloCalibPedestal*   GetPedestal()                           const { return fPedestal;             }
  Int_t                   GetnCol()                               const { return fNcol;                 }
  Int_t                   GetnRow()                               const { return fNrow;                 }
  Int_t                   GetshiftCol()                           const { return fShiftCol;             }
  Int_t                   GetshiftRow()                           const { return fShiftRow;             }
  Int_t                   GetMinL0Time()                          const { return fMinL0Time;            }
  Int_t                   GetMaxL0Time()                          const { return fMaxL0Time;            }
  Bool_t                  GetTimeCutOn()                          const { return fTimeCutOn;            }
  void                    SetCheckDeadClusters(Bool_t c)                { fCheckDeadClusters = c;       }
  void                    SetPedestal(AliCaloCalibPedestal *pds)        { fPedestal = pds;              }
  void                    SetnCol(Int_t n)                              { fNcol = n;                    }
  void                    SetnRow(Int_t n)                              { fNrow = n;                    }
  void                    SetshiftCol(Int_t n)                          { fShiftCol = n;                }
  void                    SetshiftRow(Int_t n)                          { fShiftRow = n;                }
  void                    SetTriggerClustersName(const TString &name)   { fTriggerClustersName = name;  }
  void                    SetMinL0Time(Int_t t)                         { fMinL0Time = t;               }
  void                    SetMaxL0Time(Int_t t)                         { fMaxL0Time = t;               }
  void                    SetTimeCutOn(Bool_t yes)                      { fTimeCutOn = yes;             }
  
protected:
  Int_t                 fNcol;                    // Fixed window number of cells in phi direction
  Int_t                 fNrow;                    // Fixed window number of cells in eta direction
  Int_t                 fShiftCol;                // Shifting number of cells in phi direction
  Int_t                 fShiftRow;                // Shifting number of cells in eta direction
  TString               fTriggerClustersName;     // Name of the TClonesObject that will contain the computed triggers
  Bool_t                fTimeCutOn;               // Determines whether or not apply time cuts
  Int_t                 fMinL0Time;               // Minimum L0 time cut
  Int_t                 fMaxL0Time;               // Maximum L0 time cut
  Bool_t                fCheckDeadClusters;       // Determines whether or not check for dead clusters
  AliCaloCalibPedestal *fPedestal;                //!Pointer to an object containing information about dead clusters

private:
  AliAnalysisTaskEMCALFastOR2Trigger (const AliAnalysisTaskEMCALFastOR2Trigger&);           // not implemented
  AliAnalysisTaskEMCALFastOR2Trigger operator=(const AliAnalysisTaskEMCALFastOR2Trigger&);  // not implemented
  
  ClassDef(AliAnalysisTaskEMCALFastOR2Trigger, 1); 
};

#endif

