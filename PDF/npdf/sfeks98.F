#include "pdf/pilot.h"
C***************************************************************************
C
C sfeks98.f
C
C An interface for calculating the SCALE DEPENDENT NUCLEAR PARTON
C DISTRIBUTIONS 
C 		f_A(x,Q) 
C where f_A is the distribution of parton flavour f in a proton of a 
C nucleus A, and f_p is the corresponding parton distribution in the 
C free proton.
C  
C Uses eks98 nuclear parton ratios and PDFLIB parton distributions.
C 
C
C-------------------------------------------------------------------
      subroutine sfeks98
     +     (x,Q,A,uvA,dvA,usA,dsA,stA,chA,botA,topA,glA)
#include "pdf/impdp.inc"
      data zerod/0.d0/
C      external structm,eks98
      dx = x
      dq = q
      call structm(dx,dq,uv,dv,us,ds,st,ch,bot,top,gl)
      dx = x
      dq = q
      dA = A
      call eks98(dx,dq,dA,Ruv,Rdv,Rub,Rdb,Rs,Rc,Rb,Rt,Rg)

      uvA = Ruv*uv
      dvA = Rdv*dv
      usA = Rub*us
      dsA = Rdb*ds
      stA = Rs*st
      chA = Rc*ch
      botA= Rb*bot
      topA= Rt*top
      glA = Rg*gl
            
      return
      end




