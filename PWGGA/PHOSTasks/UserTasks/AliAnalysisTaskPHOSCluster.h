#ifndef AliAnalysisTaskPHOSCluster_H
#define AliAnalysisTaskPHOSCluster_H

#include "AliCaloClusterContent.h"  // class for cluster content objects

#include <Rtypes.h>
#include <TAxis.h>
#include <TRefArray.h>
#include <TNtuple.h>
#include <vector>

#include "AliVEvent.h"
#include "AliAnalysisTaskSE.h"
#include "AliESDtrackCuts.h"
#include "AliEventCuts.h"


#include "AliVCluster.h"
#include "AliESDCaloCluster.h"
#include "AliAODCaloCluster.h"
#include "AliVCaloCells.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSCalibData.h"

class TList;
class TTree;
class TParticle;
class TH1F;
class TH2F;

using std::vector;






class AliAnalysisTaskPHOSCluster : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskPHOSCluster();
    AliAnalysisTaskPHOSCluster(const char *name);
    virtual ~AliAnalysisTaskPHOSCluster();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(const Option_t*) {}

    void           SetZVertexCut(Float_t Z)         {fZVertex = Z;}    // Set z vertex cut parameter
    void           SetMinCellsPerCluster(Int_t c)   {fMinCells = c;}   // Set minimum number of cells per cluster cut parameter
    void           SetHitmapFillingCellByCell(Bool_t b) {fFillHitmapCellByCell = b;}


//  protected:
  private:
    TList*         fOutput;                 // Output List
    TTree*         fClusterTree;            // Tree to save the AliVClusters

    Bool_t         fFillHitmapCellByCell;    // if true fill hitmap cell by cell, else cluster by cluster

    Int_t          runNumber;               //
    Int_t          fMinCells;               // Variable for MinCells cut
    Int_t          fEventCounter;           // number of analyzed events

    Float_t        fZVertex;                // Variable for Z-Vertex cut

    TString        fFiredTriggerClasses;

    AliVEvent*             fAnyEv;                    //!pointer to input event
    AliPHOSGeometry*       fPHOSGeo;                  //!PHOS geometry
    AliPHOSCalibData*      fPHOSCalibData;            //neccesary for cell by cell calibration, before filling CellID_vs_E histos.
    AliAnalysisUtils*      fUtils;                    //utils for zvtxcut
    AliVCluster*           fCluster;                  // Cluster of current event TODO remove ! maybe it helps to fill Ttree correctly
    AliVCaloCells*         fClusterCells;             // Cluster cells
    AliCaloClusterContent* fClusterCellsInformation;  // Saves information about clusters and their belonging cells

    TRefArray*         fcaloClusters;         // Cluster array

    TH1F*              fEventCluster;         //! Count event with and without clusters
    TH1F*              fEventCuts;			    //! Histogramm for event cuts
    TH1F*              fClusterCuts;          //! Histogramm for cluster cuts

    TH2F*              fPHOSHitmap;				 //! PHOS hitmap
    TH2F*              fClusterEnergyVsNoC;	 //! Histogramm cluster energy vs number of cells in cluster
    TH2F*              fClusterEnergyVsM02;	 //! Histogramm cluster energy vs M02
    TH2F*              fClusterEnergyVsM20;	 //! Histogramm cluster energy vs M20




    AliAnalysisTaskPHOSCluster(const AliAnalysisTaskPHOSCluster&);            //! Not implemented
    AliAnalysisTaskPHOSCluster& operator=(const AliAnalysisTaskPHOSCluster&); //! Not implemented
    ClassDef(AliAnalysisTaskPHOSCluster,4);
};

#endif






