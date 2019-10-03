/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliFemtoQPairCut - a simple cut which selects pairs based on the values //
// of their respective q components                                        /
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

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
 * Revision 1.2.6.1  2007/11/01 17:10:38  akisiel
 * Fix code rule conformace
 *
 * Revision 1.2  2007/05/22 09:01:42  akisiel
 * Add the possibiloity to save cut settings in the ROOT file
 *
 * Revision 1.1  2007/05/16 10:25:06  akisiel
 * Making the directory structure of AliFemtoUser flat. All files go into one common directory
 *
 * Revision 1.1.1.1  2007/04/25 15:38:41  panos
 * Importing the HBT code dir
 *
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 *
 **************************************************************************/


#ifndef ALIFEMTOQPAIRCUT_H
#define ALIFEMTOQPAIRCUT_H

// do I need these lines ?
//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "AliFemtoPairCut.h"

class AliFemtoQPairCut : public AliFemtoPairCut{
public:
  AliFemtoQPairCut();
  ~AliFemtoQPairCut();

  virtual bool Pass(const AliFemtoPair* pair);
  virtual AliFemtoString Report();
  virtual TList *ListSettings();

  void Setqlong(const float& lo, const float& hi);
  void Setqout(const float& lo, const float& hi);
  void Setqside(const float& lo, const float& hi);
  void Setqinv(const float& lo, const float& hi);
  AliFemtoQPairCut* Clone();


private:
  long fNPairsPassed;  // Number of pairs that passed the cut
  long fNPairsFailed;  // Number of pairs that failed the cut
  float fQlong[2];     // Qlong range
  float fQout[2];      // Qout range
  float fQside[2];     // Qside range
  float fQinv[2];      // Qinv range
  

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
