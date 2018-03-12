#ifndef ALICALOCLUSTERCONTENT_H
#define ALICALOCLUSTERCONTENT_H

#include <vector>
#include <Rtypes.h>
#include "TObject.h"
class AliVCluster;
class AliVCaloCells;
class AliPHOSGeometry;




  class AliCaloClusterContent : public TObject{

  public:
    AliCaloClusterContent();
    AliCaloClusterContent(const AliVCluster* clust, AliVCaloCells* cells, const AliPHOSGeometry* fgeom);
    AliCaloClusterContent(const AliCaloClusterContent &src);
    AliCaloClusterContent& operator=(const AliCaloClusterContent &src);
    virtual ~AliCaloClusterContent();


    void     Reset();                                            // Reset the object (clear all variables)
    void     SetClusterAndCells(const AliVCluster* clust, AliVCaloCells* cells, const AliPHOSGeometry* fgeom);        // Define cluster and cells

    void     SetType(Char_t type)                                   {fType = type;}

    void     SetLabe(Int_t label)                                   {fLabel = label;}
    void     SetNCells(Int_t NCells)                                {fNCells = NCells;}
    void     SetNTracksMatched(Int_t NTracksMatched)                {fNTracksMatched = NTracksMatched;}

    void     SetIsFilled(Bool_t Filled)                             {fIsFilled = Filled;}
    void     SetIsExotic(Bool_t Exotic)                             {fIsExotic = Exotic;}
    void     SetIsEMCAL(Bool_t EMCAL)                               {fIsEMCAL = EMCAL;}
    void     SetIsPHOS(Bool_t PHOS)                                 {fIsPHOS  = PHOS;}

    void     SetPosition(Float_t* Pos)                              {fPosition[0] = Pos[0]; fPosition[1] = Pos[1]; fPosition[2] = Pos[2];}

    void     SetCoreEnergy(Double_t CoreEnergy)                     {fCoreEnergy = CoreEnergy;}
    void     SetDispersion(Double_t Dispersion)                     {fDispersion = Dispersion;}
    void     SetDistanceToBadChannel(Double_t DistanceToBadChannel) {fDistanceToBadChannel = DistanceToBadChannel;}
    void     SetEmcCpvDistance(Double_t EmcCpvDistance)             {fEmcCpvDistance = EmcCpvDistance;}
    void     SetEnergy(Double_t Energy)                             {fEnergy = Energy;}
    void     SetM02(Double_t M02)                                   {fM02 = M02;}
    void     SetM20(Double_t M20)                                   {fM20 = M20;}
    void     SetTOF(Double_t TOF)                                   {fTOF = TOF;}
    void     SetTrackDx(Double_t TrackDx)                           {fTrackDx = TrackDx;}
    void     SetTrackDz(Double_t TrackDz)                           {fTrackDz = TrackDz;}

    void     SetCellsAbsID(std::vector<Int_t> CellsAbsID)                {fCellAbsID = CellsAbsID;}
    void     SetCellsDetector(std::vector<Int_t> CellsDetector)          {fCellDetector = CellsDetector;}
    void     SetCellsModule(std::vector<Int_t> CellsModule)              {fCellMod = CellsModule;}
    void     SetCellsRelIDX(std::vector<Int_t> CellsRelIDX)              {fCellRelIDX = CellsRelIDX;}
    void     SetCellsRelIDZ(std::vector<Int_t> CellsRelIDZ)              {fCellRelIDZ = CellsRelIDZ;}
    void     SetCellsEnergy(std::vector<Double_t> CellsEnergy)           {fCellEnergy = CellsEnergy;}
    void     SetCellsAmpFraction(std::vector<Double_t> CellsAmpFraction) {fCellAmpFrac = CellsAmpFraction;}
    void     SetCellsTime(std::vector<Double_t> CellsTime)               {fCellTime = CellsTime;}




	  Char_t GetType()              {return fType;}

    Int_t  GetLabel()             {return fLabel;}
    Int_t  GetNCells()            {return fNCells;}
    Int_t  GetNTracksMatched()    {return fNTracksMatched;}

    Bool_t IsFilled()             {return fIsFilled;}
    Bool_t IsExotic()             {return fIsExotic;}
    Bool_t IsEMCAL()              {return fIsEMCAL;}
    Bool_t IsPHOS()               {return fIsPHOS;}

    void  GetPosition(Float_t* x)   {x[0]=fPosition[0]; x[1]=fPosition[1]; x[2]=fPosition[2];}

    Double_t GetCoreEnergy()           {return fCoreEnergy;}
    Double_t GetDispersion()           {return fDispersion;}
    Double_t GetDistanceToBadChannel() {return fDistanceToBadChannel;}
    Double_t GetEmcCpvDistance()       {return fEmcCpvDistance;}
    Double_t GetEnergy()               {return fEnergy;}
    Double_t GetM02()                  {return fM02;}
    Double_t GetM20()                  {return fM20;}
    Double_t GetTOF()                  {return fTOF;}
    Double_t GetTrackDx()              {return fTrackDx;}
    Double_t GetTrackDz()              {return fTrackDz;}

    std::vector<Int_t> GetCellsAbsID()          {return fCellAbsID;}
    std::vector<Int_t> GetCellsDetector()       {return fCellDetector;}
    std::vector<Int_t> GetCellsModule()         {return fCellMod;}
    std::vector<Int_t> GetCellsRelIDX()         {return fCellRelIDX;}
    std::vector<Int_t> GetCellsRelIDZ()         {return fCellRelIDZ;}
    std::vector<Double_t> GetCellsEnergy()      {return fCellEnergy;}
    std::vector<Double_t> GetCellsAmpFraction() {return fCellAmpFrac;}
    std::vector<Double_t> GetCellsTime()        {return fCellTime;}

    Int_t GetVecSizeTrue() {return fCellAbsID.size();}


  protected:

	  //Char_t
    Char_t  fType;

    //Int_t
    Int_t  fLabel;          // Cluster lable
    Int_t  fNCells;         // Number of cells in cluster
    Int_t  fNTracksMatched; // Number of matched tracks

    //Bool_t
    Bool_t  fIsFilled;      // is object filled with data
	  Bool_t  fIsExotic;
    Bool_t  fIsEMCAL;
	  Bool_t  fIsPHOS;

	  //Float_t
    Float_t  fPosition[3];              // Cluster center of gravity

	  //Double_t
    Double_t  fCoreEnergy;            // Energy of the core of cluster
    Double_t  fDispersion;            // Dispersion of cluster
    Double_t  fDistanceToBadChannel;
    Double_t  fEmcCpvDistance;
    Double_t  fEnergy;                // Cluster energy
    Double_t  fM02;
    Double_t  fM20;
    Double_t  fTOF;
    Double_t  fTrackDx;
    Double_t  fTrackDz;

    std::vector<Int_t>    fCellAbsID;    // Absolute ID of cell
    std::vector<Int_t>    fCellDetector; // Detector type (0==PHOS, -1==CPV)
    std::vector<Int_t>    fCellMod;      // Module number
    std::vector<Int_t>    fCellRelIDX;   // x cell ID in module
    std::vector<Int_t>    fCellRelIDZ;   // z cell ID in module
    std::vector<Double_t> fCellEnergy;   // Energy deposited in cell
    std::vector<Double_t> fCellTime;     // Cell timing
    std::vector<Double_t> fCellAmpFrac;  // Cell amplitude fraction (== Cell energy AFTER(!) recalibration)


    ClassDef(AliCaloClusterContent, 9);

  };


#endif // ALICALOCLUSTERCONTENT_H
