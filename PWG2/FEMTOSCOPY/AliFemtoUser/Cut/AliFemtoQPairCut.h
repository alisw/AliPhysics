/***************************************************************************
 *
 * $Id $
 *
 ***************************************************************************
 *         
 *
 ***************************************************************************
 *
 * $Log$
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 *
 **************************************************************************/


#ifndef AliFemtoQPairCut_hh
#define AliFemtoQPairCut_hh

// do I need these lines ?
//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "Base/AliFemtoPairCut.h"

class AliFemtoQPairCut : public AliFemtoPairCut{
public:
  AliFemtoQPairCut();
  ~AliFemtoQPairCut();

  virtual bool Pass(const AliFemtoPair*);
  virtual AliFemtoString Report();
  void Setqlong(const float& lo, const float& hi);
  void Setqout(const float& lo, const float& hi);
  void Setqside(const float& lo, const float& hi);
  void Setqinv(const float& lo, const float& hi);
  AliFemtoQPairCut* Clone();


private:
  long fNPairsPassed;
  long fNPairsFailed;
  float fQlong[2];
  float fQout[2];
  float fQside[2];
  float fQinv[2];
  

#ifdef __ROOT__
  ClassDef(AliFemtoQPairCut, 1)
#endif
};


inline AliFemtoQPairCut* AliFemtoQPairCut::Clone() 
{ 
    AliFemtoQPairCut* c = new AliFemtoQPairCut(*this); 
    return c;
}
inline void AliFemtoQPairCut::Setqlong(const float& lo,const float& hi){fQlong[0]=lo; fQlong[1]=hi;}
inline void AliFemtoQPairCut::Setqout(const float& lo,const float& hi) {fQout[0]=lo;  fQout[1]=hi;}
inline void AliFemtoQPairCut::Setqside(const float& lo,const float& hi){fQside[0]=lo; fQside[1]=hi;}
inline void AliFemtoQPairCut::Setqinv(const float& lo,const float& hi) {fQinv[0]=lo;  fQinv[1]=hi;}

#endif
