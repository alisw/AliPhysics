//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTSPACEPOINTCONTAINER_H
#define ALIHLTSPACEPOINTCONTAINER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTSpacePointContainer.h
/// @author Matthias Richter
/// @date   2011-04-29
/// @brief  Base helper class for handling of HLT space point data blocks
///

#include <vector>
#include <cmath>
#include <TObject.h>
#include "AliHLTLogging.h"
#include "AliHLTDataTypes.h"
#include "AliHLTStdIncludes.h"

class AliHLTDataDeflater;
class TArrayC;
class TH1;
class TTree;

/**
 * @class AliHLTSpacePointContainer
 * Base class of helper classes for space point data blocks.
 * The class implements a couple of interface methods to be commonly used
 * for space point data blocks.
 *
 * @ingroup alihlt_base
 */
class AliHLTSpacePointContainer : public TObject, public AliHLTLogging
{
 public:
  /// standard constructor
  AliHLTSpacePointContainer();
  /// copy constructor
  AliHLTSpacePointContainer(const AliHLTSpacePointContainer&);
  /// assignment operator
  AliHLTSpacePointContainer& operator=(const AliHLTSpacePointContainer&);

  /// destructor
  ~AliHLTSpacePointContainer();

  //////////////////////////////////////////////////////////////////////////
  //
  // helper classes for access of space points within a grid
  //
  struct AliHLTSpacePointCell {
    int fCount;
    int fFilled;
    int fStartIndex;
  };

  class AliHLTSpacePointGrid {
  public:
    AliHLTSpacePointGrid(float maxX, float stepX,
			 float maxY, float stepY,
			 float maxZ, float stepZ,
			 int initialDataSize=-1);
    virtual ~AliHLTSpacePointGrid();

    // for now array of spacepoint ids
    typedef AliHLTUInt32_t ValueType;

    int GetDimensionX() const {return fDimX;}
    int GetDimensionY() const {return fDimY;}
    int GetDimensionZ() const {return fDimZ;}
    int GetXIndex(float x) const {
      if (x>fMaxX) return fDimX-1;
      if (x<0) return 0;
      return (int)floor(x/fStepX);
    }
    int GetYIndex(float y) const {
      if (y>fMaxY) return fDimY-1;
      if (y<0) return 0;
      return (int)floor(y/fStepY);
    }
    int GetZIndex(float z) const {
      if (z>fMaxZ) return fDimZ-1;
      if (z<0) return 0;
      return (int)floor(z/fStepZ);
    }
    float GetLowerBoundX(int cell) const {
      if (fDimX==0 || fDimY==0 ||fDimZ==0) return 0.;
      int index=cell/(fDimY*fDimZ);
      return index*fStepX;
    }
    float GetLowerBoundY(int cell) const {
      if (fDimX==0 || fDimY==0 ||fDimZ==0) return 0.;
      int index=cell%(fDimY*fDimZ); index/=fDimZ;
      return index*fStepY;
    }
    float GetLowerBoundZ(int cell) const {
      if (fDimX==0 || fDimY==0 ||fDimZ==0) return 0.;
      int index=cell%(fDimY*fDimZ); index%=fDimZ;
      return index*fStepZ;
    }
    int GetCellIndex(float x, float y, float z) const {
      return GetXIndex(x)*fDimY*fDimZ + (y<0?0:GetYIndex(y))*fDimZ + (z<0?0:GetZIndex(z));
    }
    int GetNumberOfSpacePoints(int index, int endIndex) const {
      if (!fCells) return 0;
      int count=0;
      for (int cell=index; cell<endIndex && cell<fCellDimension && count<fCount; cell++) if (fCells[cell].fCount>0) count+=fCells[cell].fCount;
      return count;
    }

    // increment counter of the cell where the spacepoint is
    int CountSpacePoint(float x, float y, float z);

    // add spacepoint, all spacepoints must have been counted before
    int AddSpacePoint(ValueType t, float x, float y, float z);

    void Clear(const char* option="");
    void Print(const char* option="");

    class iterator {
    public:
      iterator()
	: fData(NULL) {}
      iterator(const ValueType* pData)
	: fData(pData) {}
      iterator(const iterator& i)
	: fData(i.fData) {}
      iterator& operator=(const iterator& i)
      { fData=i.fData; return *this;}
      ~iterator() {fData=NULL;}

      bool operator==(const iterator& i) const  {return (fData!=NULL) && (fData==i.fData);}
      bool operator!=(const iterator& i) const  {return (fData!=NULL) && (fData!=i.fData);}
      // prefix operators
      iterator& operator++() {fData++; return *this;}
      iterator& operator--() {fData--; return *this;}
      // postfix operators
      iterator operator++(int) {iterator i(*this); fData++; return i;}
      iterator operator--(int) {iterator i(*this); fData--; return i;}

      iterator& operator+=(int step) {fData+=step; return *this;}

      const ValueType& Data() const {return *fData;}

    protected:
    private:
      const ValueType* fData; //! data
    };

    // prepare iterator and end marker
    iterator& begin(float x=-1., float y=-1., float z=-1.) {
      fIterator.~iterator();
      fIteratorEnd.~iterator();

      int startIndex=0;
      if (x<0) {
	// get all data
	if (fData) {
	  new (&fIterator) iterator(fData);
	  fIteratorEnd=fIterator;
	  fIteratorEnd+=fCount;
	}
	return fIterator;
      }

      // only search for the start index if specific x selected
      int cell=GetCellIndex(x, y, z);
      if (cell<0 || !fCells || cell>=fCellDimension) return fIterator;
      // get the index of the cell
      startIndex=fCells[cell].fStartIndex;
      if (startIndex<0 || !fData || startIndex>=fDataDimension) return fIterator;

      // get the range end position
      int endCell=cell+1;
      if (x<0) endCell=fCellDimension;
      else if (y<0) endCell=GetCellIndex(x+fStepX, -1., -1.); // all entries for fixed x
      else if (z<0) endCell=GetCellIndex(x, y+fStepY, -1.); // all entries for fixed x and y
      if (endCell<=cell) {
	// cell index returned is never outside the array
	// so this is a special case where we get to the bounds of the array
	endCell=fCellDimension;
      }

      new (&fIterator) iterator(fData+startIndex);
      fIteratorEnd=fIterator;
      fIteratorEnd+=GetNumberOfSpacePoints(cell, endCell);
      return fIterator;
    }

    // get loop end marker
    iterator& end() {
      return fIteratorEnd;
    }

  protected:
  private:
    // standard constructor prohibited
    AliHLTSpacePointGrid();
    // copy constructor prohibited
    AliHLTSpacePointGrid(const AliHLTSpacePointGrid&);
    // assignment operator prohibited
    AliHLTSpacePointGrid& operator=(const AliHLTSpacePointGrid&);

    int IndexCells();

    float fMaxX;
    float fStepX;
    float fMaxY;
    float fStepY;
    float fMaxZ;
    float fStepZ;

    int fDimX;
    int fDimY;
    int fDimZ;

    AliHLTSpacePointCell* fCells; //! cell array
    int fCellDimension; //! size of cell array
    ValueType* fData; //! spacepoint data
    int fDataDimension; //! size of spacepoint data
    int fCount;

    iterator fIterator; //! iterator
    iterator fIteratorEnd; //! end marker iterator

    static const int fgkDefaultDataSize; //! the default data size

    ClassDef(AliHLTSpacePointGrid, 0)
  };

  //////////////////////////////////////////////////////////////////////////
  //
  // interface functions
  //

  /// add input block to the collection
  virtual int AddInputBlock(const AliHLTComponentBlockData* pDesc)=0;
  virtual int PopulateAccessGrid(AliHLTSpacePointGrid* /*pGrid*/, AliHLTUInt32_t /*mask*/) const {return -ENOSYS;}
  virtual const AliHLTSpacePointGrid* GetAccessGrid(AliHLTUInt32_t /*mask*/) const {return NULL;}

  virtual int GetNumberOfSpacePoints() const;
  virtual bool Check(AliHLTUInt32_t clusterID) const;
  virtual int GetClusterIDs(vector<AliHLTUInt32_t>& tgt) const = 0;
  virtual const vector<AliHLTUInt32_t>* GetClusterIDs(AliHLTUInt32_t /*mask*/) {return NULL;}
  virtual float GetX(AliHLTUInt32_t clusterID) const = 0;
  virtual float GetXWidth(AliHLTUInt32_t clusterID) const = 0;
  virtual float GetY(AliHLTUInt32_t clusterID) const = 0;
  virtual float GetYWidth(AliHLTUInt32_t clusterID) const = 0;
  virtual float GetZ(AliHLTUInt32_t clusterID) const = 0;
  virtual float GetZWidth(AliHLTUInt32_t clusterID) const = 0;
  virtual float GetCharge(AliHLTUInt32_t clusterID) const = 0;
  virtual float GetMaxSignal(AliHLTUInt32_t /*clusterID*/) const {return 0.0;}
  virtual float GetPhi(AliHLTUInt32_t /*clusterID*/) const {return 0.0;}

  /// create a collection of clusters for a space point mask
  virtual AliHLTSpacePointContainer* SelectByMask(AliHLTUInt32_t mask, bool bAlloc=false) const;

  /// create a collection of clusters for a specific track
  virtual AliHLTSpacePointContainer* SelectByTrack(int trackId, bool bAlloc=false) const;

  /// create a collection of clusters for a specific MC track
  virtual AliHLTSpacePointContainer* SelectByMC(int mcId, bool bAlloc=false) const;

  /// create a collection of all used clusters
  virtual AliHLTSpacePointContainer* UsedClusters(bool bAlloc=false) const;

  /// create a collection of all unused clusters
  virtual AliHLTSpacePointContainer* UnusedClusters(bool bAlloc=false) const;

  int MarkUsed(AliHLTUInt32_t clusterID) {return MarkUsed(&clusterID, sizeof(clusterID));}
  virtual int MarkUsed(const AliHLTUInt32_t* clusterIDs, int arraySize);

  int SetTrackID(int trackID, AliHLTUInt32_t clusterID) {
    return SetTrackID(trackID, &clusterID, sizeof(clusterID));
  }
  virtual int SetTrackID(int trackID, const AliHLTUInt32_t* clusterIDs, int arraySize);

  virtual int GetTrackID(AliHLTUInt32_t /*clusterID*/) const {return -1;}

  int SetMCID(int mcID, AliHLTUInt32_t clusterID) {
    return SetMCID(mcID, &clusterID, sizeof(clusterID));
  }
  virtual int SetMCID(int clusterID, const AliHLTUInt32_t* clusterIDs, int arraySize);

  /// write blocks to HLT component output
  virtual int Write(AliHLTUInt8_t* outputPtr, AliHLTUInt32_t size,
		    vector<AliHLTComponentBlockData>& outputBlocks,
		    AliHLTDataDeflater* /*pDeflater*/,
		    const char* option="") const {
    return Write(outputPtr, size, outputBlocks, option);
  }

  /// write blocks to HLT component output: old function definition for backward compatibility
  virtual int Write(AliHLTUInt8_t* /*outputPtr*/, AliHLTUInt32_t /*size*/,
		    vector<AliHLTComponentBlockData>& /*outputBlocks*/,
		    const char* /*option*/="") const {return 0;}

  /// add input block from file to collection
  int AddInputBlock(const char* filename, AliHLTComponentDataType dt, unsigned specification);

  /// add input block from list of blank separated files to collection
  int AddInputBlocks(const char* filenames, AliHLTComponentDataType dt);

  /// alloc memory for a space point data block
  AliHLTUInt8_t* Alloc(int size);

  /// inherited from TObject: clear the object and reset pointer references
  virtual void Clear(Option_t * /*option*/ ="");

  /// inherited from TObject
  virtual void Print(Option_t *option="") const;

  virtual void Print(ostream& out, Option_t *option="") const;

  void Draw(Option_t *option);

  TH1* DrawProjection(const char* plane) const {
    vector<AliHLTUInt32_t> selection; // empty list -> no selection
    return DrawProjection(plane, selection);
  }

  TH1* DrawProjection(const char* plane, AliHLTUInt32_t specification) const {
    vector<AliHLTUInt32_t> selection; selection.push_back(specification);
    return DrawProjection(plane, selection);
  }

  TH1* DrawProjection(const char* plane, const vector<AliHLTUInt32_t>& selection) const;

  TTree* FillTree(const char* name, const char* title="");

 protected:

 private:
  vector<TArrayC*> fBuffers; //! buffers of loaded files

  ClassDef(AliHLTSpacePointContainer, 0)
};

ostream& operator<<(ostream &out, const AliHLTSpacePointContainer& c);

#endif
