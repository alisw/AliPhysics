#ifndef ALIMULTSELECTIONCUTS_H
#define ALIMULTSELECTIONCUTS_H
#include "TNamed.h"

using namespace std;

class AliESDEvent;

class AliMultSelectionCuts : public TNamed {
public:
    //Aliases for error codes
    typedef enum {
        kNoCalib                            = 199,
        kRejTrigger                         = 200,
        kRejINELgtZERO                      = 201,
        kRejVzCut                           = 202,
        kRejPileupInMultBins                = 203,
        kRejConsistencySPDandTrackVertices  = 204,
        kRejTrackletsVsClusters             = 205,
        kRejNonZeroNContribs                = 206,
        kRejAsymmetricInVZERO               = 207,
        kRejIncompleteDAQ                   = 208,
        kRejNotGoodVertex2016               = 209
    } CutType_t; // FIXEM Is an enum the best choice here? we return floats in the end...
    
    AliMultSelectionCuts();
    AliMultSelectionCuts(const char * name, const char * title = "Event selection cuts for multiplicity estimators");
    AliMultSelectionCuts(const AliMultSelectionCuts& o)
      : TNamed(o),
        fESD(0),
        fVzCut(o.fVzCut),
        fErrorCode(o.fErrorCode),
        fEvSel_Trig_kMB(o.fEvSel_Trig_kMB),
        fEvSel_INELgtZERO(o.fEvSel_INELgtZERO),
        fEvSel_TrackletsVsClusters(o.fEvSel_TrackletsVsClusters),
        fEvSel_RejectPileupInMultBins(o.fEvSel_RejectPileupInMultBins),
        fEvSel_CheckConsistencySPDandTrackVertices(o.fEvSel_CheckConsistencySPDandTrackVertices),
        fEvSel_NonZeroNContribs(o.fEvSel_NonZeroNContribs),
        fEvSel_IsNotAsymmetricInVZERO(o.fEvSel_IsNotAsymmetricInVZERO),
        fEvSel_IsNotIncompleteDAQ(o.fEvSel_IsNotIncompleteDAQ),
        fEvSel_HasGoodVertex2016(o.fEvSel_HasGoodVertex2016)
    {}
    AliMultSelectionCuts& operator=(const AliMultSelectionCuts& o);
    ~AliMultSelectionCuts();
    //Bool_t IsEventSelected(AliESDEvent * );
    
    void    Print(Option_t *option="") const;
    
    Float_t GetVzCut()                  { return fVzCut;      }
    void    SetVzCut (Float_t var)      { fVzCut = var;       }
    
    Bool_t GetTriggerCut()                      { return fEvSel_Trig_kMB;       }
    void   SetTriggerCut(Bool_t lSetting)       { fEvSel_Trig_kMB = lSetting;   }
    Bool_t GetINELgtZEROCut()                   { return fEvSel_INELgtZERO; }
    void   SetINELgtZEROCut(Bool_t lSetting)    { fEvSel_INELgtZERO = lSetting; }
    Bool_t GetTrackletsVsClustersCut()          { return fEvSel_TrackletsVsClusters; }
    void   SetTrackletsVsClustersCut(Bool_t lSetting)       { fEvSel_TrackletsVsClusters = lSetting; }
    Bool_t GetRejectPileupInMultBinsCut()                   { return fEvSel_RejectPileupInMultBins; }
    void   SetRejectPileupInMultBinsCut(Bool_t lSetting)    { fEvSel_RejectPileupInMultBins = lSetting; }
    Bool_t GetVertexConsistencyCut()                        { return fEvSel_CheckConsistencySPDandTrackVertices; }
    void   SetVertexConsistencyCut(Bool_t lSetting)         { fEvSel_CheckConsistencySPDandTrackVertices = lSetting; }
    Bool_t GetNonZeroNContribs()                        { return fEvSel_NonZeroNContribs; }
    void   SetNonZeroNContribs(Bool_t lSetting)         { fEvSel_NonZeroNContribs = lSetting; }
    Bool_t GetIsNotAsymmetricInVZERO()                        { return fEvSel_IsNotAsymmetricInVZERO; }
    void   SetIsNotAsymmetricInVZERO(Bool_t lSetting)         { fEvSel_IsNotAsymmetricInVZERO = lSetting; }
    Bool_t GetIsNotIncompleteDAQ()                        { return fEvSel_IsNotIncompleteDAQ; }
    void   SetIsNotIncompleteDAQ(Bool_t lSetting)         { fEvSel_IsNotIncompleteDAQ = lSetting; }
    Bool_t GetHasGoodVertex2016()                        { return fEvSel_HasGoodVertex2016; }
    void   SetHasGoodVertex2016(Bool_t lSetting)         { fEvSel_HasGoodVertex2016 = lSetting; }
    
    void SetErrorCode(Int_t lCode)   { fErrorCode = lCode; }
    Int_t GetErrorCode() const   { return fErrorCode; }
    
private:
    AliESDEvent * fESD ; //! Pointer to the ESD event, should be transient
    //Bool_t IsVzCutSelected();
    //Bool_t IsMinBias();
    //Bool_t IsINELgtZERO();
    //Bool_t HasNoInconsistentSPDandTrackVertices(AliESDEvent * lESD);
    Float_t fVzCut;     // Cut based on vertex Z position
    Int_t fErrorCode;   // Error Code used if the selection is not sucessful
                        // 0 if selected, something else if not
    
    Bool_t fEvSel_Trig_kMB; // Select kMB
    Bool_t fEvSel_INELgtZERO; // Select Inelastic greater than zero with tracklets
    Bool_t fEvSel_TrackletsVsClusters; // Apply standard Tracklets Vs Clusters cut
    Bool_t fEvSel_RejectPileupInMultBins; // Reject IsPileupFromSPDInMultBins()
    Bool_t fEvSel_CheckConsistencySPDandTrackVertices; //Check consistency
    Bool_t fEvSel_NonZeroNContribs; //at least one contributor to PV
    Bool_t fEvSel_IsNotAsymmetricInVZERO; //asymmetic in vzero
    Bool_t fEvSel_IsNotIncompleteDAQ; //incomplete DAQ
    Bool_t fEvSel_HasGoodVertex2016; //incomplete DAQ
    
    ClassDef(AliMultSelectionCuts, 4)
    //1 - original implementation 
    //2 - added NonZeroNContribs
    //3 - added IsNotAsymmetricInVZERO, IsNotIncompleteDAQ
    //4 - added HasGoodVertex2016
};
#endif
