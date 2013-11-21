#ifndef ALIITSUVERTEXERZ_H
#define ALIITSUVERTEXERZ_H

#include<AliVertexer.h>

///////////////////////////////////////////////////////////////////
//                                                               //
// Dummy vertexer plugin                                         //
//                                                               //
///////////////////////////////////////////////////////////////////

class AliESDVertex;

class AliITSUVertexer : public AliVertexer {

 public:
  AliITSUVertexer();
  virtual ~AliITSUVertexer();
  virtual AliESDVertex* FindVertexForCurrentEvent(TTree *itsClusterTree);
  void PrintStatus() const {}

 protected:

 private:
  AliITSUVertexer(const AliITSUVertexer& vtxr);
  AliITSUVertexer& operator=(const AliITSUVertexer& vtxr );

  ClassDef(AliITSUVertexer,1 );
};

#endif
