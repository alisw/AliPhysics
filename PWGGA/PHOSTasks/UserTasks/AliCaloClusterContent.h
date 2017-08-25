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
//    AliCaloClusterContent(const AliCaloClusterInfo &src);
//    AliCaloClusterContent& operator=(const AliCaloClusterContent &src);
    virtual ~AliCaloClusterContent();


    void     Reset();                                            // Reset the object (clear all variables)
    void     SetClusterAndCells(const AliVCluster* clust, AliVCaloCells* cells, const AliPHOSGeometry* fgeom);        // Define cluster and cells

	 Char_t GetType()              {return fType;}

    Int_t  GetLabel()             {return fLabel;}
    Int_t  GetNCells()            {return fNCells;}
    Int_t  GetNTracksMatched()    {return fNTracksMatched;}

	 Bool_t IsFilled()             {return fIsFilled;}
    Bool_t IsExotic()             {return fIsExotic;}
    Bool_t IsEMCAL()              {return fIsEMCAL;}
    Bool_t IsPHOS()               {return fIsPHOS;}

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

    std::vector<Int_t> GetCellsAbsID()       {return fCellAbsID;}
    std::vector<Int_t> GetCellsDetector()    {return fCellDetector;}
    std::vector<Int_t> GetCellsModule()      {return fCellMod;}
    std::vector<Int_t> GetCellsRelIDX()      {return fCellRelIDX;}
    std::vector<Int_t> GetCellsRelIDZ()      {return fCellRelIDZ;}
    std::vector<Double_t> GetCellsEnergy()   {return fCellEnergy;}
    std::vector<Double_t> GetCellsTime()     {return fCellTime;}

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


    ClassDef(AliCaloClusterContent, 4);

  };


#endif // ALICALOCLUSTERCONTENT_H
