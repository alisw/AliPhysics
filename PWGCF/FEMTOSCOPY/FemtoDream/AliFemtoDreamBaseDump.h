#ifndef ALIFEMTODREAMBASEDUMP_H_
#define ALIFEMTODREAMBASEDUMP_H_

#include "TTree.h"
#include "AliFemtoDreamBasePart.h"
#include "AliFemtoDreamEvent.h"

// ------------------------------------------------------------------------------------------------
class AliFemtoDreamPairDump {
 public:
  AliFemtoDreamPairDump();
  AliFemtoDreamPairDump(AliFemtoDreamBasePart part1,
                        AliFemtoDreamBasePart part2, float kstar);
  AliFemtoDreamPairDump &operator=(const AliFemtoDreamPairDump &obj);
  AliFemtoDreamPairDump(const AliFemtoDreamPairDump&);
  virtual ~AliFemtoDreamPairDump() {
  }

  float GetkStar() const {
    return fkstar;
  }
  AliFemtoDreamBasePart GetParticle1() const {
    return fPart1;
  }
  AliFemtoDreamBasePart GetParticle2() const {
    return fPart2;
  }

 private:
  float fkstar;
  AliFemtoDreamBasePart fPart1;
  AliFemtoDreamBasePart fPart2;

ClassDef(AliFemtoDreamPairDump, 1)
};

// ------------------------------------------------------------------------------------------------
class AliFemtoDreamEventDump {
 public:
  AliFemtoDreamEventDump();
  virtual ~AliFemtoDreamEventDump() {
  }

  void Reset();
  void SetEventProperties(float mult, float zvertex);
  void AddPair(AliFemtoDreamBasePart &part1, AliFemtoDreamBasePart &part2,
               float kstar) {
    fPairDump.push_back( { part1, part2, kstar });
  }

  float GetMultiplicity() const {
    return fMultiplicity;
  }
  float GetVertex() const {
    return fVertexZ;
  }
  size_t GetNPairs() const {
    return fPairDump.size();
  }
  std::vector<AliFemtoDreamPairDump> GetPairs() const {
    return fPairDump;
  }

 private:
  float fMultiplicity;
  float fVertexZ;
  std::vector<AliFemtoDreamPairDump> fPairDump;

ClassDef(AliFemtoDreamEventDump, 1)
};

// ------------------------------------------------------------------------------------------------
class AliFemtoDreamDump {
 public:
  AliFemtoDreamDump();
  AliFemtoDreamDump(const char *name);
  virtual ~AliFemtoDreamDump();

  // When mixing different particles
  void SetEvent(std::vector<AliFemtoDreamBasePart> &vec1,
                std::vector<AliFemtoDreamBasePart> &vec2,
                AliFemtoDreamEvent *evt, const int pdg1, const int pdg2);

  // When mixing the same particles
  void SetEvent(std::vector<AliFemtoDreamBasePart> &vec,
                AliFemtoDreamEvent *evt, const int pdg);

  TTree *GetOutput() {
    return fOutputTree;
  }
  void SetName(const char *name) {
    fOutputTree->SetName(name);
  }
  void SetkstarThreshold(float thresh) {
    fkstarThreshold = thresh;
  }

 private:
  AliFemtoDreamDump &operator=(const AliFemtoDreamDump &obj);
  AliFemtoDreamDump(const AliFemtoDreamDump&);

  float fkstarThreshold;  //
  TTree* fOutputTree;  //!
  AliFemtoDreamEventDump *fEventDump;  //!

ClassDef(AliFemtoDreamDump, 1)

};

#endif //ALIFEMTODREAMBASEDUMP_H_
