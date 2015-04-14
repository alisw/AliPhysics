#ifndef ALIEMCALAFTERBURNERUF_H
#define ALIEMCALAFTERBURNERUF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
/// After-burner for the EMCAL cluster unfolding algorithm
///
/// Input: TObjArray  *clusArray -- array of AliVClusters;
////       AliVCaloCells  *cellsEMCAL -- EMCAL cells.
///
/// Output is appended to clusArray, the original (unfolded or not) clusters
/// are deleted or moved to another position in clusArray.
///
/// If you want to use particular geometry, you must initialize it _before_
/// creating AliEMCALAfterBurnerUF instance. Add this or similar line to the
/// initialization section:
///
///    AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1_DCAL_8SM");
///
/// gGeoManager must be initialized for this code to work! Do it yourself or
/// provide geometry.root file in the current directory so that
/// AliEMCALAfterBurnerUF will take it by itself.
/// How to use:
///
///   // add this lines to the initialization section of your analysis
///   AliEMCALAfterBurnerUF *abuf = new AliEMCALAfterBurnerUF();
///   TObjArray *clusArray = new TObjArray(100);
///
///
///   AliVEvent *event = InputEvent();
///   AliVCaloCells *cellsEMCAL = event->GetEMCALCells();
///
///   for (Int_t i = 0; i < event->GetNumberOfCaloClusters(); i++)
///   {
///     AliVCluster *clus = event->GetCaloCluster(i);
///
///     clusArray->Add(clus->Clone());   // NOTE _CLONE_ in this line
///   }
///
///   abuf->UnfoldClusters(clusArray, cellsEMCAL);
///
///   // do an analysis with clusArray
///   // ....
///
///   // prevent memory leak
///   clusArray->Delete();
///
///
///  \author: Olga Driga (SUBATECH)
//-------------------------------------------------------------------------

// --- ROOT system ---
class TObjArray;
class TClonesArray;

// --- AliRoot header files ---
class AliEMCALGeometry;
class AliEMCALUnfolding;
class AliVCaloCells;

class AliEMCALAfterBurnerUF {

  public:
    AliEMCALAfterBurnerUF();
    AliEMCALAfterBurnerUF(Float_t logWeight, Float_t locMaxCut, Float_t minEcut);
    virtual ~AliEMCALAfterBurnerUF();

    virtual void Clear();
    virtual void Init();
    virtual void RecPoints2Clusters(TObjArray *clusArray);
    virtual void UnfoldClusters(TObjArray *clusArray, AliVCaloCells *cellsEMCAL);  // does the job

    // getters and setters
    virtual AliEMCALUnfolding *GetClusterUnfoldingInstance() { return fClusterUnfolding; }

  protected:
    AliEMCALGeometry  *fGeom;          // EMCAL geometry
    Float_t            fLogWeight;     // used in AliEMCALRecPoint::EvalGlobalPosition()
    Float_t            fECALocMaxCut;  // this amount of energy must distinguish a local maximum from its neighbours
    Float_t            fMinECut;       // minimum energy of cell   
    TObjArray         *fRecPoints;     //! cluster <=> recPoint
    TClonesArray      *fDigitsArr;     //->   cell <=> digit

    AliEMCALUnfolding *fClusterUnfolding;  // unfolding class instance

  private:
    AliEMCALAfterBurnerUF(const AliEMCALAfterBurnerUF & uf) ; // cpy ctor not needed, put here to avoid compilation warning 
    AliEMCALAfterBurnerUF & operator = (const AliEMCALAfterBurnerUF & uf) ;//cpy assignment, put here to avoid compilation warning 
  


    ClassDef(AliEMCALAfterBurnerUF,2)
} ;

#endif // AliEMCALAFTERBURNERUF_H
