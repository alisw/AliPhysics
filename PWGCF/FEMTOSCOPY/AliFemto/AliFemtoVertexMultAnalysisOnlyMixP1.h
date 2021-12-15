#ifndef AliFemtoVertexMultAnalysisOnlyMixP1_hh
#define AliFemtoVertexMultAnalysisOnlyMixP1_hh
//===============================================================================//
// dongfang.wang@cern.ch
// AliFemtoVertexMultAnalysisOnlyMixP1: Making pool which keep update/save
// as long as the first particle has been found in this event!
// Only difference with AliFemtoVertexMultAnalysis.cxx in line 217.
//===============================================================================//

#include "AliFemtoSimpleAnalysisOnlyMixP1.h"

class AliFemtoVertexMultAnalysisOnlyMixP1 : public AliFemtoSimpleAnalysisOnlyMixP1{

    public:
        AliFemtoVertexMultAnalysisOnlyMixP1(UInt_t binsVertex=10,
                                            Double_t minVertex=-100.0,
                                            Double_t maxVertex=+100.0,
                                            UInt_t binsMult=10,
                                            Double_t minMult=-1.0e9,
                                            Double_t maxMult=+1.0e9);
        AliFemtoVertexMultAnalysisOnlyMixP1(const AliFemtoVertexMultAnalysisOnlyMixP1 &OriAnalysis);
        AliFemtoVertexMultAnalysisOnlyMixP1& operator =(const AliFemtoVertexMultAnalysisOnlyMixP1 &OriAnalysis);
        virtual ~AliFemtoVertexMultAnalysisOnlyMixP1();


        
        virtual void ProcessEvent(const AliFemtoEvent* HbtEventToProcess);
        virtual AliFemtoString Report(); 

        virtual TList* ListSettings();

        virtual UInt_t OverflowVertexZ() const;   ///< Number of events above vertex-z range
        virtual UInt_t UnderflowVertexZ() const;  ///< Number of events below vertex-z range
        virtual UInt_t OverflowMult() const;      ///< Number of events above multiplicity range
        virtual UInt_t UnderflowMult() const;     ///< Number of events below multiplicity range

    protected:

        Double_t fVertexZ[2];     ///< min/max z-vertex position allowed to be processed
        UInt_t fVertexZBins;      ///< number of VERTEX mixing bins in z-vertex in EventMixing Buffer
        UInt_t fOverFlowVertexZ;  ///< number of events encountered which had too large z-vertex
        UInt_t fUnderFlowVertexZ; ///< number of events encountered which had too small z-vertex

        Double_t fMult[2];        ///< min/max multiplicity allowed for event to be processed
        UInt_t fMultBins;         ///< number of MULTIPLICITY mixing bins in z-vertex in EventMixing Buffer
        UInt_t fOverFlowMult;     ///< number of events encountered which had too large multiplicity
        UInt_t fUnderFlowMult;    ///< number of events encountered which had too small multiplicity

    //----------------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------------------------
    #ifdef __ROOT__
    /// \cond CLASSIMP
    ClassDef(AliFemtoVertexMultAnalysisOnlyMixP1, 1);
    /// \endcond
    #endif

};

inline UInt_t AliFemtoVertexMultAnalysisOnlyMixP1::OverflowVertexZ() const
{
  return fOverFlowVertexZ;
}

inline UInt_t AliFemtoVertexMultAnalysisOnlyMixP1::UnderflowVertexZ() const
{
  return fUnderFlowVertexZ;
}

inline UInt_t AliFemtoVertexMultAnalysisOnlyMixP1::OverflowMult() const
{
  return fOverFlowMult;
}

inline UInt_t AliFemtoVertexMultAnalysisOnlyMixP1::UnderflowMult() const
{
  return fUnderFlowMult;
}

#endif
