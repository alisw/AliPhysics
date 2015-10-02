#ifndef ALIMULTSELECTIONCUTS_H
#define ALIMULTSELECTIONCUTS_H


#include <iostream>
#include "TNamed.h"

using namespace std;

typedef enum {kMultSelTrigger = -200,
    kMultSelVzCut = -201,
    kMultSelINELgtZERO = -202,
    kMultSelTrackletsVsClusters = -203,
    kMultSelRejectPileupInMultBins = -204,
    kMultSelCheckConsistencySPDandTrackVertices = -205
} CutType_t; // FIXEM Is an enum the best choiche here? we return floats in the end...

class AliESDEvent;

class AliMultSelectionCuts : public TNamed {
    
public:
    AliMultSelectionCuts();
    AliMultSelectionCuts(const char * name, const char * title = "Event selection cuts for multiplicity estimators");
    ~AliMultSelectionCuts();
    Bool_t IsEventSelected(AliESDEvent * );
    
    void    Print(Option_t *option="") const;
    
    Float_t GetVzCut()                  { return fVzCut;      }
    void    SetVzCut (Float_t var)      { fVzCut = var;       }
    
    Bool_t GetTriggerCut() { return fEvSel_Trig_kMB; }
    void   SetTriggerCut(Bool_t lSetting) { fEvSel_Trig_kMB = lSetting; }
    Bool_t GetINELgtZEROCut() { return fEvSel_INELgtZERO; }
    void   SetINELgtZEROCut(Bool_t lSetting) { fEvSel_INELgtZERO = lSetting; }
    Bool_t GetTrackletsVsClustersCut() { return fEvSel_TrackletsVsClusters; }
    void   SetTrackletsVsClustersCut(Bool_t lSetting) { fEvSel_TrackletsVsClusters = lSetting; }
    Bool_t GetRejectPileupInMultBinsCut() { return fEvSel_RejectPileupInMultBins; }
    void   SetRejectPileupInMultBinsCut(Bool_t lSetting) { fEvSel_RejectPileupInMultBins = lSetting; }
    Bool_t GetVertexConsistencyCut() { return fEvSel_CheckConsistencySPDandTrackVertices; }
    void   SetVertexConsistencyCut(Bool_t lSetting) { fEvSel_CheckConsistencySPDandTrackVertices = lSetting; }
    
    Float_t GetErrorCode()              { return fErrorCode; }
    
private:
    AliESDEvent * fESD ; //! Pointer to the ESD event, should be transient
    Bool_t IsVzCutSelected();
    Bool_t IsMinBias();
    Bool_t IsINELgtZERO();
    Bool_t HasNoInconsistentSPDandTrackVertices(AliESDEvent * lESD);
    Float_t fVzCut;     // Cut based on vertex Z position
    Float_t fErrorCode; // Error Code used if the selection is not sucessful
    
    Bool_t fEvSel_Trig_kMB; // Select kMB
    Bool_t fEvSel_INELgtZERO; // Select Inelastic greater than zero with tracklets
    Bool_t fEvSel_TrackletsVsClusters; // Apply standard Tracklets Vs Clusters cut
    Bool_t fEvSel_RejectPileupInMultBins; // Reject IsPileupFromSPDInMultBins()
    Bool_t fEvSel_CheckConsistencySPDandTrackVertices; //Check consistency
    
    ClassDef(AliMultSelectionCuts, 1)
    
    
};

#endif
