#ifndef ALIITSUCACELL_H
#define ALIITSUCACELL_H

#include <vector>

class AliITSUCACell {
public:
	AliITSUCACell(int xx = 0u,int yy = 0u, int zz = 0u, int dd0 = 0,
                int dd1 = 0, float curv = 0.f, float n[3] = 0x0)
	  : f1OverR(curv), fd0(dd0), fd1(dd1), fN(), fVector()
	{
		fVector.reserve(4);
		fVector.push_back(xx);
		fVector.push_back(yy);
    fVector.push_back(zz);
		fVector.push_back(1u);
    if(n) {
      fN[0] = n[0];
      fN[1] = n[1];
      fN[2] = n[2];
    }
	}

	int x() const { return fVector[0]; }
	int y() const { return fVector[1]; }
  int z() const { return fVector[2]; }
	int d0() const { return fd0; }
  int d1() const { return fd1; }
  int GetLevel() const { return fVector[3]; }
  float GetCurvature() const { return f1OverR; }
  float* GetN() { return fN; }
  
	void SetLevel(int lev) { fVector[3] = lev; }

	int operator()(const int i) { return fVector[4 + i]; }

	size_t NumberOfNeighbours() { return (fVector.size() - 4u); }

	bool Combine(AliITSUCACell &neigh, int idd)
	{
    // From outside inward
		if (this->y() == neigh.z() && this->x() == neigh.y()) // Cells sharing two points
		{
      fVector.push_back(idd);
      if (neigh.GetLevel() + 1 > GetLevel())
      {
        SetLevel(neigh.GetLevel() + 1u);
      }
      return true;
    }
    return false;
  }

private:
  float f1OverR;
  int fd0,fd1;
  float fN[3];
	std::vector<int> fVector;
};

class AliITSUCARoad {
public:
  AliITSUCARoad() : Elements(), N(0)
  {
    ResetElements();
  }

  AliITSUCARoad(int layer, int idd) : Elements(), N()
  {
    ResetElements();
    N = 1;
    Elements[layer] = idd;
  }

  AliITSUCARoad(const AliITSUCARoad& copy) : Elements(), N(copy.N)
  
  {
    for ( int i=0; i<5; ++i )
    {
    	Elements[i] = copy.Elements[i];
    }
  }

  int &operator[] (const int &i) {
    return Elements[i];
  }

  void ResetElements()
  {
    for ( int i=0; i<5; ++i )
      Elements[i] = -1;
    N = 0;
  }

  void AddElement(int i, int el)
  {
    ++N;
    Elements[i] = el;
  }

  int Elements[5];
  int N;
};

#endif // ALIITSUCACELL_H
