      SUBROUTINE trig(y,x,iflag)
C
C *** DIGITISATION FOR THE MUON RAW DATA AFTER EACH EVENT ***
C *** NVE 24-SEP-1990 CERN GENEVA ***
C
C CALLED BY : GUDIGI
C ORIGIN    : NICK VAN EIJNDHOVEN
C
c	
c version 2 (open geom + L.U.T. + option DS on)
c Input  : Hits on the muon trigger chambers (RPCs) from GALICE 
c Output : muon trigger decision (Unlike Sign) L0 low and high Pt
c
c 1- TRIGMAP --> DESCRIBE MUON trigger GEOMETRY
c 2- REMPL   --> FILL bit pattern of MUON trigger CIRCUITS  
c 3- ALGO    --> INDICATE MUON trigger DECISION 
c
c
      
      real x(4,1000),y(4,1000)     ! 1000=nhitmax
      common/debug/idebug

      data nhitmax/1000/           ! max possible number of hits
      data idebug/0/               ! for debuging
      
c      call hropen(99,'hist','hist_paw.hbook',' ',1024,ISTAT)
c      open(UNIT=99,file='hist_paw.hbook',recl=1024,form='unformatted',
c     +   access='direct',status='unknown')
c      call hlimit(500000)
c      call hbook2(1,'X-Y hits plane 1(1600)',100,-300.,300.,
c     +                                       100,-300.,300.,0.)
c      call hbook2(2,'X-Y hits plane 2(1615)',100,-300.,300.,
c     +                                       100,-300.,300.,0.)
c      call hbook2(3,'X-Y hits plane 3(1700)',100,-300.,300.,
c     +                                       100,-300.,300.,0.)
c      call hbook2(4,'X-Y hits plane 4(1715)',100,-300.,300.,
c     +                                       100,-300.,300.,0.)
cc
c      call hbook1(10,'NUM of non-empty circuits',346,0.,346.,0.)  
c      call hbook1(11,'DISTR of non-empty circuits',
c     +                346,-173.,173.,0.)  
c      call hbook1(12,'NUM of non-empty X circuits',346,0.,346.,0.)  
c      call hbook1(13,'NUM of non-empty Y circuits',346,0.,346.,0.)  
cc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Describe trigger geometry        
      call TRIGMAPP
c
c Hits on trigger chambers (RPCs) from GALICE are passed through 
c      x(iplan,ihit) and y(iplan,ihit)
c with iplan=1-4 corresp. to Z=1600,1615,1700,1715 cm respectively
c      ihit =number associated to a hit (ihit<1000, see nhitmax)        
c 
c
c init
c        do iplan=1,4            
c         do ihit=1,nhitmax                       
c           x(iplan,ihit)=0          
c           y(iplan,ihit)=0
c         enddo
c        enddo
c
c for debugging proposes : ex of test sequences
c        call seqtest

c for debugging proposes : test full evt pbpb simulated with geant
c        call FULLEVT

c           
c Associate Hits to circuits (bit pattern)
            call REMPL(x,y,nhitmax)      
c Look for possible triggers (L0 , Low and High Pt cut at present)                        
            call ALGO(itrigR,itrigL0,itrigH0)
 
c  TRIGGER ?

       print *,' '
       print *,'DIMUON TRIGGER - DIMUON TRIGGER'
       print *,'ROAD +/-8   ',itrigR
       print *,'Low  pt L0  ',itrigL0
       print *,'High pt L0  ',itrigH0
       
c      call hrfile(99,'hist','n')
c      call hrout(0,icycle,' ')
c      call hrend('hist') 
c      stop
       iflag=100*itrigH0+10*itrigL0+itrigR
      end 	
c
c
c
      subroutine TRIGMAPP
c                    
c to be called at a begining of a run for dimuon trigger GEOM MAPing
c
c          
c      integer NUM(173),CODEX(173),CODEY(173)           !173=Ncirc
c      real xcmax(4,173),xcmin(4,173),
c     +     ycmax(4,173),ycmin(4,173),
c     +     xycmax(4,173),xycmin(4,173) 
c      real xwi_c(173),xwi_m(173),xwi_p(173),ywi(173)        
c      integer nstrip_c(173),nstrip_p(173),nstrip_m(173)
c 
c      common/TRIGMAP/xcmax,xcmin,ycmax,ycmin,xycmax,xycmin, 
c     +               CODEX,CODEY,
c     +               xwi_c,xwi_m,xwi_p,ywi, 
c     +               nstrip_c,nstrip_m,nstrip_p, 
c     +               NUM,Z1,Z2,Z3,Z4,Ncirc 

      common/TRIGMAP/xcmax(4,173),xcmin(4,173),ycmax(4,173),
     +               ycmin(4,173),xycmax(4,173),xycmin(4,173), 
     +               CODEX(173),CODEY(173),
     +               xwi_c(173),xwi_m(173),xwi_p(173),ywi(173), 
     +               nstrip_c(173),nstrip_m(173),nstrip_p(173), 
     +               NUM(173),Z1,Z2,Z3,Z4,Ncirc 
      integer nstrip_c,nstrip_p,nstrip_m
      integer NUM,CODEX,CODEY,Ncirc
         
      common/debug/idebug

      integer ixwi_c,ixwi_p,ixwi_m 
      integer iywi,iywi_c,iywi_p,iywi_m,ixyco_m,ixyco_p
      
      data Z1,Z2,Z3,Z4 /1600.,1615.,1700.,1715./
      data Ncirc /173/		    !Half L0 circuits Y>0 (TOT=346 )
c
c
c
c METHOD :  coding from down to up for X 
c           coding only Y>0 (if Y<0, NUM-->-NUM) 
c          
c with NUM = CODE of TRIGGER CIRCUIT
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cXXXXXXXXXX
c
c MC1 plane 1          
c X11max higher X limit of the circuit  
c X11min lower  X limit of the circuit
c
c MC1 plane 2          
c X12max=X11max*Z2/Z1      
c X12min=X11min*Z2/Z1     
c          				      	   									   
c MC2 plane 1         
c X21max=(X11max+8*XWI_p)*Z3/Z1     ex : CODEX=221 --> XWI_c = 2.12
c X21min=(X11min-8*XWI_m)*Z3/Z1         	   --> XWI_p = 1.06
c         					   --> XWI_m = 2.12
c          		            ex : CODEX=224 --> XWI_c = 2.12 
c          				           --> XWI_p = 4.24
c          				           --> XWI_m = 2.12
c          			    ex : CODEX=440 --> XWI_c = 4.24 
c          				           --> XWI_p = 0.  
c						   --> XWI_m = 4.24		            				                             
c
c MC2 plane 2         
c X22max=(X11max+8*XWI_p)*Z4/Z1                 
c X22min=(X11min-8*XWI_m)*Z4/Z1                 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cYYYYYYYYYYY
c	  
c MC1 plane 1         
c Y11max  higher Y limit of the circuit  
c Y11min  lower  Y limit of the circuit
c XY11max higher X limit of the circuit (in Y) 
c XY11min lower  X limit of the circuit (in Y)
c	  
c MC1 plane 2         
c Y12max=Y11max*Z2/Z1
c Y12min=Y11min*Z2/Z1
c XY12max=XY11max*Z2/Z1 
c XY12min=XY11min*Z2/Z1 
c	  
c MC2 plane 1         
c Y21max=Y11max*Z3/Z1
c Y21min=Y11min*Z3/Z1
c XY21max(n)=XY11max(n)  *Z3/Z1     if IXYCO_p=0
c            XY11max(n+1)*Z3/Z1     if IXYCO_p=1    
c XY21min(n)=XY11min(n)  *Z3/Z1     if IXYCO_m=0  
c            XY11min(n-1)*Z3/Z1     if IXYCO_m=1    	 						          	   
c	  
c MC2 plane 2         
c Y22max=Y11max*Z4/Z1
c Y22min=Y11min*Z4/Z1
c XY22max(n)=XY11max(n)  *Z4/Z1     if IXYCO_p=0
c            XY11max(n+1)*Z4/Z1     if IXYCO_p=1  
c XY22min(n)=XY11min(n)  *Z4/Z1     if IXYCO_m=0  
c            XY11min(n-1)*Z4/Z1     if IXYCO_m=1   
c	  
c    		
c CODEY=C1*1000+C2-*100+C2c*10+C2+*1
c      with C1=2/4 for strip width 2.12/4.24 on MC1  
c           C2-,C2c,C2+=0/1/2/3 on MC2
c           0--> nothing 1-->4strips 2-->8strips 3-->16strips  
c ex : CODEY=2332  
c      C1=2.12 (MC1 plane 1)  
c      C2-= 16 strips (MC2) 
c      C2c= 16 strips (MC2)
c      C2+= 8  strips (MC2)
c      IXYCO_p, IXYCO_m : 0/1 if C2+,C2- =/# 0
c    	     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   numerotation 1/2 plane (-xxx for other 1/2 plane)
c  
c     X
c     ^
c     !
c     ! 121 221 321 421 521 621 715 815 915 1008
c     ! 120                 620 714 814 914 1007  
c     ! 119                 619 713 813 913 1007
c     ! 118                 618 712 812 912 1006
c     ! 117                 617 711 811 911 1006 
c     ! 116                 616 710 810 910 1005
c     ! 115                 615 710 810 910 1005
c     ! 114 214 314 414     614 709 809 909 1005   
c     !         313 413 513 613 709 809 909 1005
c     !             412 512 612 708 808 908 1004
c     ! Beam        411 511 611 708 808 908 1004   --------> Y    
c     ! Shield      410 510 610 707 807 907 1004 
c     !             409 509 609 707 807 907 1004 
c     !         308 408 508 608 706 806 906 1003                                  
c     ! 107 207 307 407     607 706 906 906 1003
c     ! 106                 606 705 805 905 1003
c     ! 105                 605 705 805 905 1003
c     ! 104                 604 704 804 904 1002
c     ! 103                 603 703 803 903 1002  
c     ! 102                 602 702 802 902 1001
c     ! 101                 601 701 801 901 1001
c     ! 100 200 300 400 500 600 700 800 900 1000
c
c
cccccccccccccccccccccccccc DATAs CIRCUITS 100-107 114-121 ccccccccccccccccccc
       
       data (NUM(NN),NN=1,16) /100,101,102,103,104,105,106,107,114,115, 
     +                 116,117,118,119,120,121/
                       
       data (CODEX(NN),NN=1,16) /042,422,222,222,221,211,111,110, 
     +                   011,111,112,122,222,222,224,240/   
       
       data (CODEY(NN),NN=1,16) /4011,4110,4012,2122,2222,
     +                    2220,2020,2020,
     +                   2020,2020,2022,2222,2221,4210,4011,4110/       

       data (xcmax(1,NN),NN=1,16) /-238,-204,-170,-136,-102,-85,
     + -68,-51,
     +                     68,85,102,136,170,204,238,306/                  
       
       data (xcmin(1,NN),NN=1,16) /-306,-238,-204,-170,-136,
     +    -102,-85,-68,
     +                     51,68,85,102,136,170,204,238/ 

       data (ycmax(1,NN),NN=1,16) /17,17,17,17,17,17,17,17,
     +                     17,17,17,17,17,17,17,17/  

       data (ycmin(1,NN),NN=1,16) /0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
     

      data (xycmax(1,NN),NN=1,16) /-238,-170,-170,-136,-102,
     +           -51,-51,-51,
     +                     102,102,102,136,170,238,238,306/ 
     
      data (xycmin(1,NN),NN=1,16) /-306,-238,-238,-170,-136,
     +                      -102,-102,-102,
     +                     51,51,51,102,136,170,170,238/ 
       

cccccccccccccccccccccccccc DATAs CIRCUITS 200-207 214-221 ccccccccccccccccccc
      
      data (NUM(NN),NN=17,32) /200,201,202,203,204,205,206,207,214,215, 
     +                 216,217,218,219,220,221/
                       
      data (CODEX(NN),NN=17,32) /042,422,222,222,221,211,111,110, 
     +                   011,111,112,122,222,222,224,240/   
       
      data (CODEY(NN),NN=17,32) /4011,4110,4012,2122,2222,
     +      2220,2020,2020,
     +                   2020,2020,2022,2222,2221,4210,4011,4110/       

      data (xcmax(1,NN),NN=17,32) /-238,-204,-170,-136,-102,
     +                            -85,-68,-51,
     +                    68,85,102,136,170,204,238,306/                  
       
      data (xcmin(1,NN),NN=17,32) /-306,-238,-204,-170,-136,
     +            -102,-85,-68,
     +                     51,68,85,102,136,170,204,238/ 

      data (ycmax(1,NN),NN=17,32) /34,34,34,34,34,34,34,
     +                    34,34,34,34,34,34,34,34,34/

      data (ycmin(1,NN),NN=17,32) /17,17,17,17,17,17,17,
     +                    17,17,17,17,17,17,17,17,17/  
    
      data (xycmax(1,NN),NN=17,32) /-238,-170,-170,-136,-102,
     +              -51,-51,-51,
     +                      102,102,102,136,170,238,238,306/ 
     
      data (xycmin(1,NN),NN=17,32) /-306,-238,-238,-170,-136,
     +                     -102,-102,-102,
     +                     51,51,51,102,136,170,170,238/ 
 

cccccccccccccccccccccccccc DATAs CIRCUITS 300-308 313-321 ccccccccccccccccccc
      
       data (NUM(NN),NN=33,50) /300,301,302,303,304,305,306,307,308, 
     +                  313,314,315,316,317,318,319,320,321/
                              
       data (CODEX(NN),NN=33,50) /042,422,222,222,221,211,111,111,110, 
     +                    011,111,111,112,122,222,222,224,240/   
       
       data (CODEY(NN),NN=33,50) /4011,4110,4012,2122,2222,2220,
     + 2022,2220,2020,
     +                    2020,2022,2220,2022,2222,2221,4210,4011,4110/       

       data (xcmax(1,NN),NN=33,50) /-238,-204,-170,-136,-102,-85,
     +                   -68,-51,-34, 
     +                     51,68,85,102,136,170,204,238,306/                  
       
       data (xcmin(1,NN),NN=33,50) /-306,-238,-204,-170,-136,
     +  -102,-85,-68,-51,   
     +                     34,51,68,85,102,136,170,204,238/ 

       data (ycmax(1,NN),NN=33,50) /51,51,51,51,51,51,51,
     +                     51,51,51,51,51,51,51,51,51,51,51/  

       data (ycmin(1,NN),NN=33,50) /34,34,34,34,34,34,34,
     +                     34,34,34,34,34,34,34,34,34,34,34/
    
      data (xycmax(1,NN),NN=33,50) /-238,-170,-170,-136,-102,
     +                   -68,-68,-34,-34,
     +                     68,68,102,102,136,170,238,238,306/ 
     
      data (xycmin(1,NN),NN=33,50) /-306,-238,-238,-170,-136,
     +                        -102,-102,-68,-68,
     +                     34,34,68,68,102,136,170,170,238/ 
     

cccccccccccccccccccccccccc DATAs CIRCUITS 400 - 421 ccccccccccccccccccccccccc

       data (NUM(NN),NN=51,72) /400,401,402,403,404,405,406,
     + 407,408,409,410,
     +                  411,412,413,414,415,416,417,418,419,420,421/
       
       data (CODEX(NN),NN=51,72) /042,422,222,222,221,211,111,111,111,
     + 111,111, 
     +                    111,111,111,111,111,112,122,222,222,224,240/   
       
       data (CODEY(NN),NN=51,72) /4011,4110,4012,2122,2222,2220,
     +  2022,2220,2022,
     +                    2220,2022,2220,2022,
     +                    2220,2022,2220,2022,2222,2221,4210,4011,4110/       

       data (xcmax(1,NN),NN=51,72) /-238,-204,-170,-136,-102,-85,
     +                 -68,-51,-34,
     +                     -17,0,17,34, 
     +                     51,68,85,102,136,170,204,238,306/                  
       
       data (xcmin(1,NN),NN=51,72) /-306,-238,-204,-170,-136,-102,
     +                    -85,-68,-51,
     +                     -34,-17,0,17,   
     +                     34,51,68,85,102,136,170,204,238/ 

       data (ycmax(1,NN),NN=51,72) /68,68,68,68,68,68,68,68,68,68,68,
     +                     68,68,68,68,68,68,68,68,68,68,68/

       data (ycmin(1,NN),NN=51,72) /51,51,51,51,51,51,51,51,51,51,51,
     +                     51,51,51,51,51,51,51,51,51,51,51/  
    
      data (xycmax(1,NN),NN=51,72) /-238,-170,-170,-136,-102,-68,
     +        -68,-34,-34,
     +                     0,0,34,34,
     +                     68,68,102,102,136,170,238,238,306/ 
     
      data (xycmin(1,NN),NN=51,72) /-306,-238,-238,-170,-136,
     +                     -102,-102,-68,-68,
     +                     -34,-34,0,0,  
     +                     34,34,68,68,102,136,170,170,238/ 
     

cccccccccccccccccccccccccc DATAs CIRCUITS 500 - 521 ccccccccccccccccccccccccc

      data (NUM(NN),NN=73,94) /500,501,502,503,504,505,506,
     +                   507,508,509,510,
     +                  511,512,513,514,515,516,517,518,519,520,521/
       
      data (CODEX(NN),NN=73,94) /042,422,222,222,221,211,111,
     +            111,111,111,111, 
     +                    111,111,111,111,111,112,122,222,222,224,240/   
       
      data (CODEY(NN),NN=73,94) /4011,4110,4012,2122,2222,
     +               2220,2022,2220,2022,
     +                    2220,2022,2220,2022,
     +                    2220,2022,2220,2022,2222,2221,4210,4011,4110/       

      data (xcmax(1,NN),NN=73,94) /-238,-204,-170,-136,
     +                     -102,-85,-68,-51,-34,
     +                     -17,0,17,34, 
     +                     51,68,85,102,136,170,204,238,306/                  
       
      data (xcmin(1,NN),NN=73,94) /-306,-238,-204,-170,-136,
     +                -102,-85,-68,-51,
     +                     -34,-17,0,17,   
     +                     34,51,68,85,102,136,170,204,238/ 

      data (ycmax(1,NN),NN=73,94) /85,85,85,85,85,85,85,85,85,85,85,
     +                     85,85,85,85,85,85,85,85,85,85,85/  

      data (ycmin(1,NN),NN=73,94) /68,68,68,68,68,68,68,68,68,68,68,
     +                     68,68,68,68,68,68,68,68,68,68,68/    

      data (xycmax(1,NN),NN=73,94) /-238,-170,-170,-136,-102,
     +            -68,-68,-34,-34,
     +                      0,0,34,34,
     +                      68,68,102,102,136,170,238,238,306/ 
     
      data (xycmin(1,NN),NN=73,94) /-306,-238,-238,-170,
     +               -136,-102,-102,-68,
     +                      -68,-34,-34,0,0,  
     +                      34,34,68,68,102,136,170,170,238/ 


cccccccccccccccccccccccccc DATAs CIRCUITS 600 - 621 ccccccccccccccccccccccccc

      data (NUM(NN),NN=95,116)/600,601,602,603,604,605,606,
     +           607,608,609,610,
     +                  611,612,613,614,615,616,617,618,619,620,621/
       
      data (CODEX(NN),NN=95,116)/042,422,222,222,221,211,111,
     +      111,111,111,111, 
     +                    111,111,111,111,111,112,122,222,222,224,240/   
       
      data (CODEY(NN),NN=95,116)/4011,4110,4012,2122,2222,
     +                    2220,2022,2220,2022,
     +                    2220,2022,2220,2022,
     +                    2220,2022,2220,2022,2222,2221,4210,4011,4110/       

      data (xcmax(1,NN),NN=95,116)/-238,-204,-170,-136,
     +                     -102,-85,-68,-51,-34,
     +                     -17,0,17,34, 
     +                     51,68,85,102,136,170,204,238,306/                  
       
      data (xcmin(1,NN),NN=95,116)/-306,-238,-204,-170,-136,
     +                    -102,-85,-68,-51,
     +                     -34,-17,0,17,   
     +                     34,51,68,85,102,136,170,204,238/ 

      data (ycmax(1,NN),NN=95,116)/102,102,102,102,102,102,
     +               102,102,102,102,
     +               102,102,102,102,102,102,102,102,102,102,102,102/    

      data (ycmin(1,NN),NN=95,116)/85,85,85,85,85,85,85,85,85,85,85,
     +                     85,85,85,85,85,85,85,85,85,85,85/  

      data (xycmax(1,NN),NN=95,116)/-238,-170,-170,-136,-102,-68,
     +                -68,-34,-34,
     +                      0,0,34,34,
     +                      68,68,102,102,136,170,238,238,306/ 
     
      data (xycmin(1,NN),NN=95,116)/-306,-238,-238,-170,-136,
     +              -102,-102,-68,
     +                      -68,-34,-34,0,0,  
     +                      34,34,68,68,102,136,170,170,238/ 
     

cccccccccccccccccccccccccc DATAs CIRCUITS 700 - 715 ccccccccccccccccccccccccc

      data (NUM(NN),NN=117,132)/700,701,702,703,704,705,706,707,
     +                  708,709,710,711,712,713,714,715/
       
      data (CODEX(NN),NN=117,132)/042,422,222,222,222,222,222,222, 
     +                    222,222,222,222,222,222,224,240/   
       
      data (CODEY(NN),NN=117,132)/4022,4220,4023,2233,2333,
     +         2333,2333,2333,
     +                    2333,2333,2333,2333,2332,4320,4022,4220/       

      data (xcmax(1,NN),NN=117,132)/-238,-204,-170,-136,-102,-68,-34,0,
     +                     34,68,102,136,170,204,238,306/                  
       
      data (xcmin(1,NN),NN=117,132)/-306,-238,-204,-170,-136,
     +  -102,-68,-34,
     +                     0,34,68,102,136,170,204,238/   

      data (ycmax(1,NN),NN=117,132)/136,136,136,136,136,136,136,136,
     +                     136,136,136,136,136,136,136,136/  

      data (ycmin(1,NN),NN=117,132)/102,102,102,102,102,102,102,102,
     +                     102,102,102,102,102,102,102,102/    

      data (xycmax(1,NN),NN=117,132)/-238,-170,-170,
     +            -136,-102,-68,-34,0,
     +                      34,68,102,136,170,238,238,306/ 
     
      data (xycmin(1,NN),NN=117,132)/-306,-238,-238,-170,
     +            -136,-102,-68,-34,
     +                      0,34,68,102,136,170,170,238/ 
     
     

cccccccccccccccccccccccccc DATAs CIRCUITS 800 - 815 ccccccccccccccccccccccccc

      data (NUM(NN),NN=133,148)/800,801,802,803,804,805,806,807,
     +                  808,809,810,811,812,813,814,815/
       
      data (CODEX(NN),NN=133,148)/042,422,222,222,222,222,222,222, 
     +                    222,222,222,222,222,222,224,240/   
       
      data (CODEY(NN),NN=133,148)/4022,4220,4023,2233,2333,
     +              2333,2333,2333,
     +                    2333,2333,2333,2333,2332,4320,4022,4220/       

      data (xcmax(1,NN),NN=133,148)/-238,-204,-170,-136,-102,-68,-34,0,
     +                     34,68,102,136,170,204,238,306/                  
       
      data (xcmin(1,NN),NN=133,148)/-306,-238,-204,-170,
     +-136,-102,-68,-34,
     +                     0,34,68,102,136,170,204,238/   

      data (ycmax(1,NN),NN=133,148)/170,170,170,170,170,170,170,170,
     +                     170,170,170,170,170,170,170,170/    

      data (ycmin(1,NN),NN=133,148)/136,136,136,136,136,136,136,136,
     +                     136,136,136,136,136,136,136,136/  

      data (xycmax(1,NN),NN=133,148)/-238,-170,-170,-136,
     +-102,-68,-34,0,
     +                      34,68,102,136,170,238,238,306/ 
     
      data (xycmin(1,NN),NN=133,148)/-306,-238,-238,-170,
     +      -136,-102,-68,-34,
     +                      0,34,68,102,136,170,170,238/ 
     
     
cccccccccccccccccccccccccc DATAs CIRCUITS 900 - 915 ccccccccccccccccccccccccc  

      data (NUM(NN),NN=149,164)/900,901,902,903,904,905,906,907,
     +                  908,909,910,911,912,913,914,915/
       
      data (CODEX(NN),NN=149,164)/042,422,222,222,222,222,222,222, 
     +                    222,222,222,222,222,222,224,240/   
       
      data (CODEY(NN),NN=149,164)/4022,4220,4022,4220,
     +4022,4220,4022,4220,
     +                    4022,4220,4022,4220,4022,4220,4022,4220/       

      data (xcmax(1,NN),NN=149,164)/-238,-204,-170,-136,-102,-68,
     +-34,0,
     +                     34,68,102,136,170,204,238,306/                  
       
      data (xcmin(1,NN),NN=149,164)/-306,-238,-204,-170,-136,-102,
     +              -68,-34,
     +                     0,34,68,102,136,170,204,238/   

      data (ycmax(1,NN),NN=149,164)/204,204,204,204,204,204,204,204,
     +                     204,204,204,204,204,204,204,204/  

      data (ycmin(1,NN),NN=149,164)/170,170,170,170,170,170,170,170,
     +                     170,170,170,170,170,170,170,170/    

      data (xycmax(1,NN),NN=149,164)/-238,-170,-170,-102,-102,
     +-34,-34,34,
     +                      34,102,102,170,170,238,238,306/ 
     
      data (xycmin(1,NN),NN=149,164)/-306,-238,-238,-170,-170,
     +-102,-102,-34,
     +                      -34,34,34,102,102,170,170,238/ 
     
     
cccccccccccccccccccccccccc DATAs CIRCUITS 1000 - 1008 ccccccccccccccccccccccc  

      data (NUM(NN),NN=165,173)/1000,1001,1002,1003,1004,1005, 
     +                  1006,1007,1008/
       
      data (CODEX(NN),NN=165,173)/044,444,444,444,444,444,444,444,440/ 
       
      data (CODEY(NN),NN=165,173)/4033,4333,4333,4333,4333,4333,
     +                    4333,4333,4330/       

      data (xcmax(1,NN),NN=165,173)/-238,-170,-102,-34,
     +                     34,102,170,238,306/                  
       
      data (xcmin(1,NN),NN=165,173)/-306,-238,-170,-102,-34,
     +                     34,102,170,238/   

      data (ycmax(1,NN),NN=165,173)/272,272,272,272,272,272,
     +272,272,272/

      data (ycmin(1,NN),NN=165,173)/204,204,204,204,204,204,
     +204,204,204/  

      data (xycmax(1,NN),NN=165,173)/-238,-170,-102,-34,
     +                      34,102,170,238,306/  
     
      data (xycmin(1,NN),NN=165,173)/-306,-238,-170,-102,-34,
     +                      34,102,170,238/ 
     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

calculate other quantities associated to a circuit on the other planes
c          
     
      do i=1,ncirc                           !loop on all circuits (Y>0)
c        
c DECODE CODEX
c
      ixwi_m=CODEX(i)/100                    !Dividing integers
      ixwi_c=(CODEX(i)-ixwi_m*100)/10
      ixwi_p=CODEX(i)-ixwi_m*100-ixwi_c*10
c
      if(ixwi_m.eq.0)then
       xwi_m(i)=0.
      elseif(ixwi_m.eq.1)then 
       xwi_m(i)=1.0625
      elseif(ixwi_m.eq.2)then 
       xwi_m(i)=2.125
      elseif(ixwi_m.eq.4)then 
       xwi_m(i)=4.25
      else
       print *,'WARNING :: BAD CODEX value' 
      endif    
      
      if(ixwi_c.eq.0)then
       xwi_c(i)=0.
      elseif(ixwi_c.eq.1)then 
       xwi_c(i)=1.0625
      elseif(ixwi_c.eq.2)then 
       xwi_c(i)=2.125
      elseif(ixwi_c.eq.4)then 
       xwi_c(i)=4.25
      else
       print *,'WARNING :: BAD CODEX value' 
      endif    
 
      if(ixwi_p.eq.0)then
       xwi_p(i)=0.
      elseif(ixwi_p.eq.1)then 
       xwi_p(i)=1.0625
      elseif(ixwi_p.eq.2)then 
       xwi_p(i)=2.125
      elseif(ixwi_p.eq.4)then 
       xwi_p(i)=4.25
      else
       print *,'WARNING :: BAD CODEX value' 
      endif    
       
c ready to calculate XGEOM parameters
     
      xcmax(2,i)=Z2/Z1*(xcmax(1,i))     
      xcmin(2,i)=Z2/Z1*(xcmin(1,i))     
c           				        	                				      	   									   
c          
      xcmax(3,i)=Z3/Z1*(xcmax(1,i)+8*XWI_p(i))      
      xcmin(3,i)=Z3/Z1*(xcmin(1,i)-8*XWI_m(i))         	        
c         					
      xcmax(4,i)=Z4/Z1*(xcmax(1,i)+8*XWI_p(i))                 
      xcmin(4,i)=Z4/Z1*(xcmin(1,i)-8*XWI_m(i))                 
c
c        
c DECODE CODEY
c
      iywi=CODEY(i)/1000                 !Dividing integers
      iywi_m=(CODEY(i)-iywi*1000)/100 
      iywi_c=(CODEY(i)-iywi*1000-iywi_m*100)/10
      iywi_p=CODEY(i)-iywi*1000-iywi_m*100-iywi_c*10

      if(iywi.eq.2)then
       ywi(i)=2.125
      elseif(iywi.eq.4)then 
       ywi(i)=4.25
      else
       print *,'WARNING :: BAD CODEY value' 
      endif    

      if(iywi_m.eq.0)then
       ixyco_m=0
       nstrip_m(i)=0
      elseif(iywi_m.eq.1)then 
       ixyco_m=1
       nstrip_m(i)=4
      elseif(iywi_m.eq.2)then 
       ixyco_m=1
       nstrip_m(i)=8
      elseif(iywi_m.eq.3)then 
       ixyco_m=1
       nstrip_m(i)=16
      else
       print *,'WARNING :: BAD CODEY value' 
      endif    

      if(iywi_c.eq.0)then
       nstrip_c(i)=0
      elseif(iywi_c.eq.1)then 
       nstrip_c(i)=4
      elseif(iywi_c.eq.2)then 
       nstrip_c(i)=8
      elseif(iywi_c.eq.3)then 
       nstrip_c(i)=16
      else
       print *,'WARNING :: BAD CODEY value' 
      endif    
      
      if(iywi_p.eq.0)then
       ixyco_p=0
       nstrip_p(i)=0
      elseif(iywi_p.eq.1)then 
       ixyco_p=1
       nstrip_p(i)=4
      elseif(iywi_p.eq.2)then 
       ixyco_p=1
       nstrip_p(i)=8
      elseif(iywi_p.eq.3)then 
       ixyco_p=1
       nstrip_p(i)=16
      else
       print *,'WARNING :: BAD CODEY value' 
      endif    
       
c ready to calculate YGEOM parameters
    		
    		
      ycmax(2,i)=Z2/Z1*(ycmax(1,i))
      ycmin(2,i)=Z2/Z1*(ycmin(1,i))
      xycmax(2,i)=Z2/Z1*(xycmax(1,i)) 
      xycmin(2,i)=Z2/Z1*(xycmin(1,i)) 
	  
      ycmax(3,i)=Z3/Z1*(ycmax(1,i))
      ycmin(3,i)=Z3/Z1*(ycmin(1,i))

      if(ixyco_p.eq.0)then
       xycmax(3,i)=Z3/Z1*(xycmax(1,i))
      else 
       xycmax(3,i)=Z3/Z1*(xycmax(1,i+1))
      endif

      if(ixyco_m.eq.0)then
       xycmin(3,i)=Z3/Z1*(xycmin(1,i))
      else 
       xycmin(3,i)=Z3/Z1*(xycmin(1,i-1))
      endif	 						          	   

      ycmax(4,i)=Z4/Z1*(ycmax(1,i))
      ycmin(4,i)=Z4/Z1*(ycmin(1,i))
	  
      if(ixyco_p.eq.0)then
       xycmax(4,i)=Z4/Z1*(xycmax(1,i))
      else 
       xycmax(4,i)=Z4/Z1*(xycmax(1,i+1))
      endif

      if(ixyco_m.eq.0)then
       xycmin(4,i)=Z4/Z1*(xycmin(1,i))
      else 
       xycmin(4,i)=Z4/Z1*(xycmin(1,i-1))
      endif	 						          	   
c    		
c
      enddo               !loop on trigger boards
c
      return 
      end
c
cc
      subroutine REMPL(x,y,nhitmax)
c Associate hits to circuits (bit pattern) 

      real x(4,1000),y(4,1000)  ! 1000=nhitmax

      real xu(346),yu(346)                 ! something   in x / y  
      integer*4 bitpx(4,346),bitpy(4,346)  ! bit-pattern in x / y code DECIMAL                                             
      common /hitcirc/xu,yu,bitpx,bitpy                 

c      integer NUM(173),CODEX(173),CODEY(173)           !173=Ncirc
c      real xcmax(4,173),xcmin(4,173),
c     +     ycmax(4,173),ycmin(4,173),
c     +     xycmax(4,173),xycmin(4,173) 
c      real xwi_c(173),xwi_m(173),xwi_p(173),ywi(173)        
c      integer nstrip_c(173),nstrip_p(173),nstrip_m(173)
c 
c      common/TRIGMAP/xcmax,xcmin,ycmax,ycmin,xycmax,xycmin, 
c     +               CODEX,CODEY,
c     +               xwi_c,xwi_m,xwi_p,ywi, 
c     +               nstrip_c,nstrip_m,nstrip_p, 
c     +               NUM,Z1,Z2,Z3,Z4,Ncirc 
      common/TRIGMAP/xcmax(4,173),xcmin(4,173),ycmax(4,173),
     +               ycmin(4,173),xycmax(4,173),xycmin(4,173), 
     +               CODEX(173),CODEY(173),
     +               xwi_c(173),xwi_m(173),xwi_p(173),ywi(173), 
     +               nstrip_c(173),nstrip_m(173),nstrip_p(173), 
     +               NUM(173),Z1,Z2,Z3,Z4,Ncirc 
      integer nstrip_c,nstrip_p,nstrip_m
      integer NUM,CODEX,CODEY,Ncirc

     
      common/debug/idebug

      character*32 i2bin      


cinit 
      do icircT=1,(ncirc*2)
       xu(icircT)=0
       yu(icircT)=0
       do iplan=1,4
        bitpx(iplan,icircT)=0
        bitpy(iplan,icircT)=0
       enddo
      enddo
c
c fill the hits in histograms
c
c      do iplan=1,4
c       do ihit=1,nhitmax
c        if (x(iplan,ihit).ne.0.or.y(iplan,ihit).ne.0)
c     +	call hfill(iplan,x(iplan,ihit),y(iplan,ihit),1.)
c       enddo
c      enddo    	
c
c fill the circuits in X and Y
c
      do iplan=1,4
       do ihit=1,nhitmax
        if (x(iplan,ihit).ne.0.or.y(iplan,ihit).ne.0)then !save CPU time
c           print *,iplan,ihit,x(iplan,ihit),y(iplan,ihit)
c
ccccccccccccccccccccc BIT PATTERN in X cccccccccccccccccccccccccccc
c       
        do icirc=1,ncirc
         yabs=abs(y(iplan,ihit))
         if(x(iplan,ihit).ge.xcmin(iplan,icirc).and.
     +      x(iplan,ihit).lt.xcmax(iplan,icirc))then
          if(yabs.ge.ycmin(iplan,icirc).and.
     +       yabs.lt.ycmax(iplan,icirc))then
           if(y(iplan,ihit).ge.0)xu(icirc)=1  
           if(y(iplan,ihit).lt.0)xu(icirc+ncirc)=1                   
c 
           if(iplan.eq.1)then
            X11=x(iplan,ihit)-xcmin(iplan,icirc)        
            nx11=ifix(X11/xwi_c(icirc))+1              !no of hitten strip 

            if(nX11.le.0.or.nX11.ge.17)
     +         print *,'WARNING : NX11 out of range',NX11

            if(y(iplan,ihit).ge.0)
     +         call sbit1(bitpx(iplan,icirc),nx11)     !sets bit pos nx11 =1
            if(y(iplan,ihit).lt.0)
     +         call sbit1(bitpx(iplan,icirc+ncirc),nx11)   

           
            if (idebug.eq.2)then 
            print *,' ' 
            print *,'iplan,ihit,icirc,NUM=',iplan,ihit,icirc,num(icirc)
            print *,'x y xcmin xcmax ycmin ycmax ',x(iplan,ihit),
     +       y(iplan,ihit),xcmin(iplan,icirc),xcmax(iplan,icirc),               
     +       ycmin(iplan,icirc),ycmax(iplan,icirc) 
            print *,'X11 nX11 xwidth=',X11,nX11,xwi_c(icirc)
            print *,'bitpx(iplan,icirc) 	',
     +        i2bin(bitpx(iplan,icirc),16)       
            print *,'bitpx(iplan,icirc+ncirc) 	',
     +        i2bin(bitpx(iplan,icirc+ncirc),16)  
            endif      

           elseif(iplan.eq.2)then
c             
            X12=x(iplan,ihit)-xcmin(iplan,icirc) 
            nx12=ifix(X12/(xwi_c(icirc)*Z2/Z1))+1      !no of hitten strip 
            if(nX12.le.0.or.nX12.ge.17)
     +         print *,'WARNING : NX12 out of range',NX12

             if(y(iplan,ihit).ge.0)
     +          call sbit1(bitpx(iplan,icirc),nx12)    !sets bit pos nx12 =1
             if(y(iplan,ihit).lt.0)
     +          call sbit1(bitpx(iplan,icirc+ncirc),nx12)   


            if (idebug.eq.2)then 
            print *,' ' 
            print *,'iplan,ihit,icirc,NUM=',iplan,ihit,icirc,num(icirc)
            print *,'x y xcmin xcmax ycmin ycmax ',x(iplan,ihit),
     +       y(iplan,ihit),xcmin(iplan,icirc),xcmax(iplan,icirc),               
     +       ycmin(iplan,icirc),ycmax(iplan,icirc) 
            print *,'X12 nX12 xwidth',X12,nX12,xwi_c(icirc)
            print *,'bitpx(iplan,icirc) 	',
     +        i2bin(bitpx(iplan,icirc),16)       
            print *,'bitpx(iplan,icirc+ncirc) 	',
     +        i2bin(bitpx(iplan,icirc+ncirc),16)       
            endif
           
           elseif(iplan.eq.3)then
c
            X21=x(iplan,ihit)-xcmin(iplan,icirc)        
             if(X21.le.(xwi_m(icirc)*8*Z3/Z1))then   
              nX21=ifix(X21/(xwi_m(icirc)*Z3/Z1))+1    !no of hitten strip 
             elseif
     +        (X21.le.((xwi_m(icirc)*8+xwi_c(icirc)*16)*Z3/Z1))then  
              XX21=X21-xwi_m(icirc)*8*Z3/Z1 
              nX21=ifix(XX21/(xwi_c(icirc)*Z3/Z1))+1+8  
             else
              XXX21=X21-(xwi_m(icirc)*8+xwi_c(icirc)*16)*Z3/Z1 
              nX21=ifix(XXX21/(xwi_p(icirc)*Z3/Z1))+1+8+16   
             endif
                               
            if(nX21.le.0.or.nX21.ge.33)
     +         print *,'WARNING : NX21 out of range',NX21

            if(y(iplan,ihit).ge.0)
     +         call sbit1(bitpx(iplan,icirc),nX21)     !sets bit pos nX21 =1
            if(y(iplan,ihit).lt.0)
     +         call sbit1(bitpx(iplan,icirc+ncirc),nX21)   


            if (idebug.eq.2)then 
            print *,' ' 
            print *,'iplan,ihit,icirc,NUM=',iplan,ihit,icirc,num(icirc)
            print *,'x y xcmin xcmax ycmin ycmax ',x(iplan,ihit),
     +       y(iplan,ihit),xcmin(iplan,icirc),xcmax(iplan,icirc),               
     +       ycmin(iplan,icirc),ycmax(iplan,icirc) 
            print *,'X21 nX21 xwidth-c xwidth-m xwidth-p =',
     +       X21,nX21,xwi_c(icirc),xwi_m(icirc),xwi_p(icirc)
            print *,'bitpx(iplan,icirc) 	',
     +        i2bin(bitpx(iplan,icirc),32)       
            print *,'bitpx(iplan,icirc+ncirc) 	',
     +        i2bin(bitpx(iplan,icirc+ncirc),32)       
            endif

           elseif(iplan.eq.4)then
c
            x22=x(iplan,ihit)-xcmin(iplan,icirc)        
             if(x22.le.(xwi_m(icirc)*8*Z4/Z1))then   
              nx22=ifix(x22/(xwi_m(icirc)*Z4/Z1))+1    !no of hitten strip  
             elseif
     +        (x22.le.((xwi_m(icirc)*8+xwi_c(icirc)*16)*Z4/Z1))then  
              Xx22=x22-xwi_m(icirc)*8*Z4/Z1 
              nx22=ifix(Xx22/(xwi_c(icirc)*Z4/Z1))+1+8     
             else
              XXx22=x22-(xwi_m(icirc)*8+xwi_c(icirc)*16)*Z4/Z1 
              nx22=ifix(XXx22/(xwi_p(icirc)*Z4/Z1))+1+8+16   
             endif 

            if(nx22.le.0.or.nx22.ge.33)
     +         print *,'WARNING : Nx22 out of range',Nx22
                
            if(y(iplan,ihit).ge.0)
     +         call sbit1(bitpx(iplan,icirc),nx22)     !sets bit pos nx22 =1
            if(y(iplan,ihit).lt.0)
     +         call sbit1(bitpx(iplan,icirc+ncirc),nx22)   



            if (idebug.eq.2)then 
            print *,' ' 
            print *,'iplan,ihit,icirc,NUM=',iplan,ihit,icirc,num(icirc)
            print *,'x y xcmin xcmax ycmin ycmax ',x(iplan,ihit),
     +       y(iplan,ihit),xcmin(iplan,icirc),xcmax(iplan,icirc),               
     +       ycmin(iplan,icirc),ycmax(iplan,icirc) 
            print *,'x22 nx22 xwidth-c xwidth-m xwidth-p =',
     +       x22,nx22,xwi_c(icirc),xwi_m(icirc),xwi_p(icirc)
            print *,'bitpx(iplan,icirc) 	',
     +        i2bin(bitpx(iplan,icirc),32)       
            print *,'bitpx(iplan,icirc+ncirc) 	',
     +        i2bin(bitpx(iplan,icirc+ncirc),32)       
            endif


           endif              
          endif
         endif
         
ccccccccccccccccccccc BIT PATTERN in Y cccccccccccccccccccccccccccc

         if(x(iplan,ihit).ge.xycmin(iplan,icirc).and.
     +      x(iplan,ihit).lt.xycmax(iplan,icirc))then
          if(yabs.ge.ycmin(iplan,icirc).and.
     +       yabs.lt.ycmax(iplan,icirc))then
           if(y(iplan,ihit).ge.0)yu(icirc)=1  
           if(y(iplan,ihit).lt.0)yu(icirc+ncirc)=1                    
c
           if(iplan.eq.1)then
            Y11=yabs-ycmin(iplan,icirc)        
            nY11=ifix(Y11/ywi(icirc))+1                !no of hitten strip  
            nstrip=nstrip_c(icirc) 
            if(nY11.le.0.or.nY11.ge.(nstrip+1))
     +         print *,'WARNING : NY11 out of range',NY11

            if(y(iplan,ihit).ge.0)
     +         call sbit1(bitpy(iplan,icirc),nY11)     !sets bit pos nY11 =1
            if(y(iplan,ihit).lt.0)
     +         call sbit1(bitpy(iplan,icirc+ncirc),nY11)   

             
            if (idebug.eq.2)then 
            print *,' ' 
            print *,'iplan,ihit,icirc,NUM=',iplan,ihit,icirc,num(icirc)
            print *,'x y xycmin xycmax ycmin ycmax ',x(iplan,ihit),
     +       y(iplan,ihit),xycmin(iplan,icirc),xycmax(iplan,icirc),               
     +       ycmin(iplan,icirc),ycmax(iplan,icirc) 
            print *,'Y11 nY11 nstripy ywidth=',
     +               Y11,nY11,nstrip,ywi(icirc)
            print *,'bitpy(iplan,icirc) 	',
     +        i2bin(bitpy(iplan,icirc),nstrip)       
            print *,'bitpy(iplan,icirc+ncirc) 	',
     +        i2bin(bitpy(iplan,icirc+ncirc),nstrip)       
             endif


           elseif(iplan.eq.2)then
c 
            Y12=yabs-ycmin(iplan,icirc)        
            nY12=ifix(Y12/(ywi(icirc)*Z2/Z1))+1        !no of hitten strip 
            nstrip=nstrip_c(icirc) 
            if(nY12.le.0.or.nY12.ge.(nstrip+1))
     +         print *,'WARNING : NY12 out of range',NY12
     
            if(y(iplan,ihit).ge.0)
     +         call sbit1(bitpy(iplan,icirc),nY12)     !sets bit pos nY12 =1
            if(y(iplan,ihit).lt.0)
     +         call sbit1(bitpy(iplan,icirc+ncirc),nY12)   


            if (idebug.eq.2)then 
            print *,' ' 
            print *,'iplan,ihit,icirc,NUM=',iplan,ihit,icirc,num(icirc)
            print *,'x y xycmin xycmax ycmin ycmax ',x(iplan,ihit),
     +       y(iplan,ihit),xycmin(iplan,icirc),xycmax(iplan,icirc),               
     +       ycmin(iplan,icirc),ycmax(iplan,icirc) 
            print *,'Y12 nY12 nstripy ywidth=',
     +               Y12,nY12,nstrip,ywi(icirc)
            print *,'bitpy(iplan,icirc) 	',
     +        i2bin(bitpy(iplan,icirc),nstrip)       
            print *,'bitpy(iplan,icirc+ncirc) 	',
     +        i2bin(bitpy(iplan,icirc+ncirc),nstrip)       
            endif

           elseif(iplan.eq.3)then
c 
            Y21=yabs-ycmin(iplan,icirc) 
            nstrip=nstrip_c(icirc) 
            nY21=ifix(Y21/(ywi(icirc)*Z3/Z1))+1        !no of hitten strip 
             if(nY21.le.0.or.nY21.ge.(nstrip+1))
     +         print *,'WARNING : NY21 out of range',NY21
     
             if(y(iplan,ihit).ge.0)
     +         call sbit1(bitpy(iplan,icirc),nY21)     !sets bit pos nY21 =1
             if(y(iplan,ihit).lt.0)
     +         call sbit1(bitpy(iplan,icirc+ncirc),nY21)   

c various Y strip width in the same circuit
	    if(nstrip_p(icirc).ne.0.and.nstrip_p(icirc).lt.nstrip)then	
             xyextrap=xycmax(1,icirc)*Z3/Z1
              if(x(iplan,ihit).gt.xyextrap)then    
               if((mod(ny21,2)).eq.1)then             
                if(y(iplan,ihit).ge.0)
     +           call sbit1(bitpy(iplan,icirc),nY21+1) !sets bit pos nY21+1 =1
                if(y(iplan,ihit).lt.0)
     +           call sbit1(bitpy(iplan,icirc+ncirc),nY21+1)   
               endif
               if((mod(ny21,2)).eq.0)then            
                if(y(iplan,ihit).ge.0)
     +           call sbit1(bitpy(iplan,icirc),nY21-1) !sets bit pos nY21-1 =1
                if(y(iplan,ihit).lt.0)
     +           call sbit1(bitpy(iplan,icirc+ncirc),nY21-1)   
               endif
              endif  
	    endif

	    if(nstrip_m(icirc).ne.0.and.nstrip_m(icirc).lt.nstrip)then	
             xyextram=xycmin(1,icirc)*Z3/Z1
              if(x(iplan,ihit).lt.xyextram)then    
               if((mod(ny21,2)).eq.1)then             
                if(y(iplan,ihit).ge.0)
     +           call sbit1(bitpy(iplan,icirc),nY21+1) !sets bit pos nY21+1 =1
                if(y(iplan,ihit).lt.0)
     +           call sbit1(bitpy(iplan,icirc+ncirc),nY21+1)   
               endif
               if((mod(ny21,2)).eq.0)then           
                if(y(iplan,ihit).ge.0)
     +           call sbit1(bitpy(iplan,icirc),nY21-1) !sets bit pos nY21-1 =1
                if(y(iplan,ihit).lt.0)
     +           call sbit1(bitpy(iplan,icirc+ncirc),nY21-1)   
               endif
              endif  
	    endif

            if (idebug.eq.2)then 
            print *,' ' 
            print *,'iplan,ihit,icirc,NUM=',iplan,ihit,icirc,num(icirc)
            print *,'x y xycmin xycmax ycmin ycmax ',x(iplan,ihit),
     +       y(iplan,ihit),xycmin(iplan,icirc),xycmax(iplan,icirc),               
     +       ycmin(iplan,icirc),ycmax(iplan,icirc) 
            print *,'Y21 nY21 nstripc ywidth= nstripm nstripp',
     +      Y21,nY21,nstrip,ywi(icirc),nstrip_m(icirc),nstrip_p(icirc)
            print *,'bitpy(iplan,icirc) 	',
     +        i2bin(bitpy(iplan,icirc),nstrip)       
            print *,'bitpy(iplan,icirc+ncirc) 	',
     +        i2bin(bitpy(iplan,icirc+ncirc),nstrip)       
	    if(nstrip_p(icirc).ne.0.and.nstrip_p(icirc).lt.nstrip)then	
            print *,'special changt de largeur Y : p plus larges'
            print *,'xyextrap, mod(ny21,2)=',xyextrap, mod(ny21,2)
            endif 
	    if(nstrip_m(icirc).ne.0.and.nstrip_m(icirc).lt.nstrip)then	
            print *,'special changt de largeur Y : m plus larges'
            print *,'xyextram, mod(ny21,2)=',xyextrap, mod(ny21,2)
            endif 
            endif

           elseif(iplan.eq.4)then
c 
            Y22=yabs-ycmin(iplan,icirc) 
            nstrip=nstrip_c(icirc) 
             nY22=ifix(Y22/(ywi(icirc)*Z4/Z1))+1       !no of hitten strip        
             if(nY22.le.0.or.nY22.ge.(nstrip+1))
     +         print *,'WARNING : NY22 out of range',NY22
     
             if(y(iplan,ihit).ge.0)
     +         call sbit1(bitpy(iplan,icirc),nY22)     !sets bit pos nY22 =1
             if(y(iplan,ihit).lt.0)
     +         call sbit1(bitpy(iplan,icirc+ncirc),nY22)   
c
c various Y strip width in the same circuit
	    if(nstrip_p(icirc).ne.0.and.nstrip_p(icirc).lt.nstrip)then	
             xyextrap=xycmax(1,icirc)*Z4/Z1
              if(x(iplan,ihit).gt.xyextrap)then    
               if((mod(ny22,2)).eq.1)then           
                if(y(iplan,ihit).ge.0)
     +           call sbit1(bitpy(iplan,icirc),nY22+1) !sets bit pos nY22+1 =1
                if(y(iplan,ihit).lt.0)
     +           call sbit1(bitpy(iplan,icirc+ncirc),nY22+1)   
               endif
               if((mod(ny22,2)).eq.0)then             
                if(y(iplan,ihit).ge.0)
     +           call sbit1(bitpy(iplan,icirc),nY22-1) !sets bit pos nY22-1 =1
                if(y(iplan,ihit).lt.0)
     +           call sbit1(bitpy(iplan,icirc+ncirc),nY22-1)   
               endif
              endif  
	    endif

	    if(nstrip_m(icirc).ne.0.and.nstrip_m(icirc).lt.nstrip)then	
             xyextram=xycmin(1,icirc)*Z4/Z1
              if(x(iplan,ihit).lt.xyextram)then   
               if((mod(ny22,2)).eq.1)then            
                if(y(iplan,ihit).ge.0)
     +           call sbit1(bitpy(iplan,icirc),nY22+1) !sets bit pos nY22+1 =1
                if(y(iplan,ihit).lt.0)
     +           call sbit1(bitpy(iplan,icirc+ncirc),nY22+1)   
               endif
               if((mod(ny22,2)).eq.0)then             
                if(y(iplan,ihit).ge.0)
     +           call sbit1(bitpy(iplan,icirc),nY22-1) !sets bit pos nY22-1 =1
                if(y(iplan,ihit).lt.0)
     +           call sbit1(bitpy(iplan,icirc+ncirc),nY22-1)   
               endif
              endif  
	    endif
	    

            if (idebug.eq.2)then 
            print *,' ' 
            print *,'iplan,ihit,icirc,NUM=',iplan,ihit,icirc,num(icirc)
            print *,'x y xycmin xycmax ycmin ycmax ',x(iplan,ihit),
     +       y(iplan,ihit),xycmin(iplan,icirc),xycmax(iplan,icirc),               
     +       ycmin(iplan,icirc),ycmax(iplan,icirc) 
            print *,'Y22 nY22 nstripc ywidth= nstripm nstripp',
     +      Y22,nY22,nstrip,ywi(icirc),nstrip_m(icirc),nstrip_p(icirc)
            print *,'bitpy(iplan,icirc) 	',
     +        i2bin(bitpy(iplan,icirc),nstrip)       
            print *,'bitpy(iplan,icirc+ncirc) 	',
     +        i2bin(bitpy(iplan,icirc+ncirc),nstrip)       
	    if(nstrip_p(icirc).ne.0.and.nstrip_p(icirc).lt.nstrip)then	
            print *,'special changt de largeur Y : p plus larges'
            print *,'xyextrap, mod(ny22,2)=',xyextrap,mod(ny22,2)
            endif 
	    if(nstrip_m(icirc).ne.0.and.nstrip_m(icirc).lt.nstrip)then	
            print *,'special changt de largeur Y : m plus larges'
            print *,'xyextram, mod(ny22,2)=',xyextrap, mod(ny22,2)
            endif 
            endif


           endif                                       !iplan=1-4
c
          endif
         endif                  
c       
        enddo                                          !circuit
        endif                                          !x-y ne 0  
       enddo                                           !hits     
      enddo                                            !plan 

         
      return
      end
c
cc
      subroutine ALGO(itrigR,itrigL0,itrigH0)

     
cc     
c      integer NUM(173),CODEX(173),CODEY(173)           !173=Ncirc
c      real xcmax(4,173),xcmin(4,173),
c     +     ycmax(4,173),ycmin(4,173),
c     +     xycmax(4,173),xycmin(4,173) 
c      real xwi_c(173),xwi_m(173),xwi_p(173),ywi(173)        
c      integer nstrip_c(173),nstrip_p(173),nstrip_m(173)
c 
c      common/TRIGMAP/xcmax,xcmin,ycmax,ycmin,xycmax,xycmin, 
c     +               CODEX,CODEY,
c     +               xwi_c,xwi_m,xwi_p,ywi, 
c     +               nstrip_c,nstrip_m,nstrip_p, 
c     +               NUM,Z1,Z2,Z3,Z4,Ncirc 
      common/TRIGMAP/xcmax(4,173),xcmin(4,173),ycmax(4,173),
     +               ycmin(4,173),xycmax(4,173),xycmin(4,173), 
     +               CODEX(173),CODEY(173),
     +               xwi_c(173),xwi_m(173),xwi_p(173),ywi(173), 
     +               nstrip_c(173),nstrip_m(173),nstrip_p(173), 
     +               NUM(173),Z1,Z2,Z3,Z4,Ncirc 
      integer nstrip_c,nstrip_p,nstrip_m
      integer NUM,CODEX,CODEY,Ncirc
     
      common/debug/idebug
     
cc
      real xu(346),yu(346)                 ! something   in x / y  
      integer*4 bitpx(4,346),bitpy(4,346)  ! bit-pattern in x / y code DECIMAL
                                            
      common /hitcirc/xu,yu,bitpx,bitpy                 

      character*32 i2bin      

      integer*4 dble1_0,dble1_1_14,dble1_15
      integer*4 dble2_0,dble2_1_30,dble2_31
      integer*4 sgle1_0,sgle1_1_14,sgle1_15
      integer*4 sgle2_0,sgle2_1_30,sgle2_31
      integer*4 sgle1A,sgle1B,sgle2A,sgle2B 
      integer*4 x1_1,x1_2a,x1_2b,x1_2c,x2_1,x2_2a,x2_2b,x2_2c
      integer*4 dblex1(346),dblex2(346) 
      integer*4 sglex1(346),sglex2(346) 
      integer*4 dbley1(346),dbley2(346) 
      integer*4 sgley1(346),sgley2(346) 
      integer*4 thrl(346)
      integer*4 co_l(16,346),co_y(16,346)
      integer*4 ib,jb,iabit,ibbit,icbit,ibit_l,ibit_y  
      integer*4 sign_l(346),val_y(346) 
      integer*4 sign_lv  
c L0 L.U.T. (et L2)
      integer*4 dnp,dnm,dsup,dinf 
      integer*4 dev_2(346),num_x2(346),num_y2(346),signdev_2(346)         
      integer*4 devmin(16),stripnum(16),idevmin(16)         
      real yL2(346),x1L2(346),x2L2(346)
c          
      JBIT (IZW,IZP)     = IBITS (IZW,IZP-1,1)
      JBYT (IZW,IZP,NZB) = IBITS (IZW,IZP-1,NZB)
 

      
c Threshold Road MAX +/-8      
      data thrl /346*131071/  ! +/- 8 strips 
c      data thrl /346*4064/  ! +/- 3 strips 
c
c Datas for L.U.T. calculations
      data zF /975./ 
      data ptcalLow,ptcalHigh /.60,1.6/
c
c init 
c
      do k=1,346
         sign_l(k)=0
         val_y(k)=0
         yL2(k)=0
         x1L2(k)=0
         x2L2(k)=0
      enddo
c
c histogram number and distribution of hitten circuits
c

      nHxcirc=0
      nHycirc=0
      nHcirc=0
      do icircT=1,ncirc*2
       if(xu(icircT).eq.1)nHxcirc=nHxcirc+1   
       if(yu(icircT).eq.1)nHycirc=nHycirc+1   
       if(xu(icircT).eq.1.or.yu(icircT).eq.1)nHcirc=nHcirc+1  !circuit non-vide
       if(xu(icircT).eq.1.or.yu(icircT).eq.1)then
        if(icircT.le.ncirc)then
c         call hf1(11,float(icircT),1.)
        else
c         call hf1(11,-float(icircT-ncirc),1.)
        endif
       endif   
      enddo
c      call hf1(10,float(nHcirc),1.)    
c      call hf1(12,float(nHxcirc),1.)    
c      call hf1(13,float(nHycirc),1.)  
        
c
c Calculate Singles and Doubles (including MINI-ROADS)
c
      do icircT=1,ncirc*2
      if(xu(icircT).eq.1)then  !loop on non-empty X circuit 
c       
c memo : integer=ibits(integer,start_pos,lenght), start_pos_min=0
c        
c MC1
       x1_1=ibits(bitpx(1,icircT),0,1)     
       x1_2a=ibits(bitpx(2,icircT),0,1)
       x1_2b=ibits(bitpx(2,icircT),1,1)
       dble1_0=iand(x1_1,ior(x1_2a,x1_2b))
       
       x1_1 =ibits(bitpx(1,icircT),1,14)
       x1_2a=ibits(bitpx(2,icircT),0,14)          
       x1_2b=ibits(bitpx(2,icircT),1,14)          
       x1_2c=ibits(bitpx(2,icircT),2,14)          
       dble1_1_14=iand(x1_1,ior(x1_2a,ior(x1_2b,x1_2c)))
       
       x1_1=ibits(bitpx(1,icircT),15,1)
       x1_2a=ibits(bitpx(2,icircT),14,1)
       x1_2b=ibits(bitpx(2,icircT),15,1)
       dble1_15=iand(x1_1,ior(x1_2a,x1_2b))
       
       dblex1(icircT)=dble1_0+2*dble1_1_14+(2**15)*dble1_15  !dble X on MC1

       x1_1=ibits(bitpx(1,icircT),0,16)
       sgle1A=ieor(x1_1,dblex1(icircT))

       x1_2a=ibits(dblex1(icircT),0,1)
       x1_2b=ibits(dblex1(icircT),1,1)
       sgle1_0=ior(x1_2a,x1_2b)

       x1_2a=ibits(dblex1(icircT),0,14)          
       x1_2b=ibits(dblex1(icircT),1,14)          
       x1_2c=ibits(dblex1(icircT),2,14)          
       sgle1_1_14=ior(x1_2a,ior(x1_2b,x1_2c))

       x1_2a=ibits(dblex1(icircT),14,1)
       x1_2b=ibits(dblex1(icircT),15,1)
       sgle1_15=ior(x1_2a,x1_2b)
              
       sgle1B=sgle1_0+2*sgle1_1_14+(2**15)*sgle1_15
       sgle1B=not(sgle1B)
       sgle1B=iand(sgle1B,bitpx(2,icircT))
       
       sglex1(icircT)=ior(sgle1A,sgle1B)		     !sgle X on MC1
      


       dbley1(icircT)=iand(bitpy(1,icircT),bitpy(2,icircT))     !dble Y on MC1
       sgley1(icircT)=
     + ieor(ior(bitpy(1,icircT),bitpy(2,icircT)),dbley1(icircT))!sgle Y on MC1


       if(idebug.eq.3)then         
         if(icircT.le.ncirc)then
          nocirc=NUM(icircT)
         else
          nocirc=-NUM(icircT-ncirc)
         endif  
         print *,' '
         print *,'no circuit ',icircT,nocirc
         print *,'bitpx(1,icircT)         ',i2bin(bitpx(1,icircT),16)
         print *,'bitpx(2,icircT)         ',i2bin(bitpx(2,icircT),16)
         print *,'dblex1(icircT)          ',i2bin(dblex1(icircT),16)            
         print *,'sglex1(icircT)          ',i2bin(sglex1(icircT),16)            
         print *,' '
         print *,'bitpy(1,icircT)         ',i2bin(bitpy(1,icircT),16)
         print *,'bitpy(2,icircT)         ',i2bin(bitpy(2,icircT),16)
         print *,'dbley1(icircT)          ',i2bin(dbley1(icircT),16)            
         print *,'sgley1(icircT)          ',i2bin(sgley1(icircT),16)            
          
       endif  


c MC2       
       x2_1=ibits(bitpx(3,icircT),0,1)     
       x2_2a=ibits(bitpx(4,icircT),0,1)
       x2_2b=ibits(bitpx(4,icircT),1,1)
       dble2_0=iand(x2_1,ior(x2_2a,x2_2b))
       
       x2_1 =ibits(bitpx(3,icircT),1,30)
       x2_2a=ibits(bitpx(4,icircT),0,30)          
       x2_2b=ibits(bitpx(4,icircT),1,30)          
       x2_2c=ibits(bitpx(4,icircT),2,30)          
       dble2_1_30=iand(x2_1,ior(x2_2a,ior(x2_2b,x2_2c)))
       
       x2_1=ibits(bitpx(3,icircT),31,1)
       x2_2a=ibits(bitpx(4,icircT),30,1)
       x2_2b=ibits(bitpx(4,icircT),31,1)
       dble2_31=iand(x2_1,ior(x2_2a,x2_2b))
       
       dblex2(icircT)=dble2_0+2*dble2_1_30+(2**31)*dble2_31  !dble X on MC2

       x2_1=ibits(bitpx(3,icircT),0,32)
       sgle2A=ieor(x2_1,dblex2(icircT))

       x2_2a=ibits(dblex2(icircT),0,1)
       x2_2b=ibits(dblex2(icircT),1,1)
       sgle2_0=ior(x2_2a,x2_2b)

       x2_2a=ibits(dblex2(icircT),0,30)          
       x2_2b=ibits(dblex2(icircT),1,30)          
       x2_2c=ibits(dblex2(icircT),2,30)          
       sgle2_1_30=ior(x2_2a,ior(x2_2b,x2_2c))

       x2_2a=ibits(dblex2(icircT),30,1)
       x2_2b=ibits(dblex2(icircT),31,1)
       sgle2_31=ior(x2_2a,x2_2b)
              
       sgle2B=sgle2_0+2*sgle2_1_30+(2**31)*sgle2_31
       sgle2B=not(sgle2B)
       sgle2B=iand(sgle2B,bitpx(4,icircT))
       
       sglex2(icircT)=ior(sgle2A,sgle2B)		      !sgle X on MC2
       

       dbley2(icircT)=iand(bitpy(3,icircT),bitpy(4,icircT))     !dble Y on MC2
       sgley2(icircT)=
     + ieor(ior(bitpy(3,icircT),bitpy(4,icircT)),dbley2(icircT))!sgle Y on MC2


       if(idebug.eq.3)then         
         print*,' '
         print *,'bitpx(3,icircT) ',i2bin(bitpx(3,icircT),32)
         print *,'bitpx(4,icircT) ',i2bin(bitpx(4,icircT),32)
         print *,'dblex2(icircT)  ',i2bin(dblex2(icircT),32)            
         print *,'sglex2(icircT)  ',i2bin(sglex2(icircT),32)            
         print *,' '
         print *,'bitpy(3,icircT) ',i2bin(bitpy(3,icircT),16)
         print *,'bitpy(4,icircT) ',i2bin(bitpy(4,icircT),16)
         print *,'dbley2(icircT)  ',i2bin(dbley2(icircT),16)            
         print *,'sgley2(icircT)  ',i2bin(sgley2(icircT),16)            
       endif  

       endif                    !non-empty X circuit
      enddo                     !circuits
c
c
c method DS (X only) : reduce the sensibility of the trigger to single hits  
c                      without loosing signal  
c
       do icircT=1,ncirc*2
        if(xu(icircT).eq.1)then  !loop on non-empty X circuit to save CPU time
         if(dblex1(icircT).ne.0)sglex1(icircT)=0
         if(dblex2(icircT).ne.0)sglex2(icircT)=0
        endif
       enddo 
c         
c
c
c coincidence 3/4 in +/- 8 strips road	
c	
      do icircT=1,2*ncirc
      if(xu(icircT).eq.1)then  !loop on non-empty X circuit   
       do ib=1,16
        do jb=ib,ib+16 
         iabit=iand(jbit(sglex1(icircT),ib),jbit(dblex2(icircT),jb)) 
         ibbit=ior(jbit(sglex2(icircT),jb),jbit(dblex2(icircT),jb))
         icbit=ior(iabit,iand(ibbit,jbit(dblex1(icircT),ib)))
         ibit_l=iand(icbit,jbit(thrl(icircT),jb-ib+1))
	 call sbit(ibit_l,co_l(ib,icircT),jb-ib+1)         
        enddo
       enddo
        
       if(idebug.eq.4)then         
        if (xu(icircT).ne.0)then
         if(icircT.le.ncirc)then
          nocirc=NUM(icircT)
         else
          nocirc=-NUM(icircT-ncirc)
         endif  
         print *,' '
         print *,'no circuit ',icircT,nocirc
         print *,'sglex1(icircT)          ',i2bin(sglex1(icircT),16)            
         print *,'dblex1(icircT)          ',i2bin(dblex1(icircT),16)            
         print *,'sglex2(icircT)  ',i2bin(sglex2(icircT),32)            
         print *,'dblex2(icircT)  ',i2bin(dblex2(icircT),32)            
         print *,' '
         print *,'thrl(icircT)             ',i2bin(thrl(icircT),17)            
          do ib=1,16
         print *,'co_l(ib,icircT)          ',i2bin(co_l(ib,icircT),17)            
          enddo

        endif
       endif
c
c signe info (X) before L.U.T.  
c
       ibit_l=0
       do ib=1,16
        ibit_l=ior(ibit_l,co_l(ib,icircT))
       enddo
        ibgauche=ibits(ibit_l,9,8)
        ibcentre=ibits(ibit_l,8,1)
        ibdroite=ibits(ibit_l,0,8)
       if(ibgauche.eq.0.and.ibcentre.eq.0.and.ibdroite.eq.0)then
        sign_l(icircT)=0        			!00
       elseif(ibgauche.eq.0.and.ibcentre.eq.0)then
        sign_l(icircT)=1        			!01
       elseif(ibdroite.eq.0.and.ibcentre.eq.0)then
        sign_l(icircT)=2        			!10
       else
        sign_l(icircT)=3     			        !11
       endif 


       if(idebug.eq.4)then         
        if (xu(icircT).ne.0)then
         print*,' '
         print *,'sign_l(icircT)          ',i2bin(sign_l(icircT),2)            
        endif
       endif
      endif                       !non-empty X circuits 
      enddo                       !circuits 
c 
c coincidences 3/4 Y (constant ROAD +/- 1)
c
      do icircT=1,ncirc*2
      if(xu(icircT).eq.1)then  !loop on non-empty X circuit  
      
c      
         if(icircT.le.ncirc)then
          nstripy=nstrip_c(icircT)
         else
          nstripy=nstrip_c(icircT-ncirc)
         endif  
c
       do ib=1,nstripy
        do jb=ib,ib+2                      !+/- 1 strip  
            
         iabit=iand(jbit(sgley1(icircT),ib),jbit(2*dbley2(icircT),jb)) 
         ibbit=ior(jbit(2*sgley2(icircT),jb),jbit(2*dbley2(icircT),jb))
         ibit_y=ior(iabit,iand(ibbit,jbit(dbley1(icircT),ib)))
	 call sbit(ibit_y,co_y(ib,icircT),jb-ib+1)         
        enddo
       enddo
c validation Y (OUTPUT L0 Y) 
       ibit_y=0
       do ib=1,nstripy
        ibit_y=ior(ibit_y,co_y(ib,icircT))  
       enddo
       

       if(ibit_y.eq.0)then 
        val_y(icircT)=0	          !nothing in Y 
       else
        val_y(icircT)=1           !something in Y
       endif
         
       if(idebug.eq.5)then         
        if (xu(icircT).ne.0)then
         if(icircT.le.ncirc)then
          nocirc=NUM(icircT)
         else
          nocirc=-NUM(icircT-ncirc)
         endif  
         print *,' '
         print *,'no circuit ',icircT,nocirc
         print *,'sgley1(icircT)  ',i2bin(sgley1(icircT),nstripy)            
         print *,'dbley1(icircT)  ',i2bin(dbley1(icircT),nstripy)            
         print *,'sgley2(icircT)  ',i2bin(sgley2(icircT),nstripy)            
         print *,'dbley2(icircT)  ',i2bin(dbley2(icircT),nstripy)            
         print *,' '
          do ib=1,nstripy
         print *,'co_y(ib,icircT)          ',i2bin(co_y(ib,icircT),3)            
          enddo
         print *,' '
         print *,'val_y(icircT)          ',i2bin(val_y(icircT),1)            

        endif
       endif
      endif                     !non-empty X circuits 
      enddo                     !circuits
c       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c infos in +/- 8 road are available : sign_l(346) and val_y(346)           
c where sign_l = 2 bits integer (00 nothing, 01 right dev.,              
c	        	         10 left dev., 11 zero dev. or ambiguity)
c and val_y=0 (nothing) or 1 (validation of X) 				    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Calculate trigger response in +/-8 strip road (for debugging)   
c
      igl_0=0
      idl_0=0 
      igh_0=0
      idh_0=0 
      il_0=0
      ih_0=0
      itrigR=0 
       
      

      do icircT=1,2*ncirc
      if(xu(icircT).eq.1)then  !loop on non-empty X circuit 

c with Y validation
       sign_lv=sign_l(icircT)*val_y(icircT) 
       if (sign_lv.eq.1.or.sign_lv.eq.3)idl_0=idl_0+1        
       if (sign_lv.eq.2.or.sign_lv.eq.3)igl_0=igl_0+1        
       if (sign_lv.ge.1)il_0=il_0+1        

      endif                   !non-empty X circuits 
      enddo                   !circuits


      if(il_0.ge.2.and.
     + idl_0.ne.0.and.igl_0.ne.0)itrigR=1       !U.S. trigger Road +/-8 
*       if(il_0.ge.1)itrigR=1 
c
c
c
c
c INTERFACE L0-L.U.T. (and L2) in X : 
c       nO CIRCUIT + Hitten STRIP (MC1) + DEVMIN (MC2) + SIGNDEV
c
c  L.U.T. (and L2) calculation ARE done only if a LOCAL trigger 
c                in +/-8 strips (sign_lv.ne.0) is found
c
c
      do icircT=1,ncirc*2
      sign_lv=sign_l(icircT)*val_y(icircT) 
       if(sign_lv.ne.0)then                        !local trigger +/- 8 strips
c
        do ib=1,16
         if(co_l(ib,icircT).eq.0)then
          devmin(ib)=100                            !dummy    
         elseif(iand(256,co_l(ib,icircT)).ne.0)then  
          devmin(ib)=0                              !central strip
         elseif(iand(128,co_l(ib,icircT)).ne.0)then    
          devmin(ib)=-1                             ! DOWN priority and <0 
         elseif(iand(512,co_l(ib,icircT)).ne.0)then    
          devmin(ib)=1                              ! UP >0 (convention x>0)  
         elseif(iand(64,co_l(ib,icircT)).ne.0)then    
          devmin(ib)=-2                                
         elseif(iand(1024,co_l(ib,icircT)).ne.0)then    
          devmin(ib)=2                                
         elseif(iand(32,co_l(ib,icircT)).ne.0)then    
          devmin(ib)=-3                                
         elseif(iand(2048,co_l(ib,icircT)).ne.0)then    
          devmin(ib)=3                                
         elseif(iand(16,co_l(ib,icircT)).ne.0)then    
          devmin(ib)=-4                                
         elseif(iand(4096,co_l(ib,icircT)).ne.0)then    
          devmin(ib)=4                                
         elseif(iand(8,co_l(ib,icircT)).ne.0)then    
          devmin(ib)=-5                                
         elseif(iand(8192,co_l(ib,icircT)).ne.0)then    
          devmin(ib)=5                                
         elseif(iand(4,co_l(ib,icircT)).ne.0)then    
          devmin(ib)=-6                                
         elseif(iand(16384,co_l(ib,icircT)).ne.0)then    
          devmin(ib)=6                                
         elseif(iand(2,co_l(ib,icircT)).ne.0)then    
          devmin(ib)=-7                                
         elseif(iand(32768,co_l(ib,icircT)).ne.0)then    
          devmin(ib)=7                                
         elseif(iand(1,co_l(ib,icircT)).ne.0)then    
          devmin(ib)=-8                                
         elseif(iand(65536,co_l(ib,icircT)).ne.0)then    
          devmin(ib)=8
         endif 
        enddo                                 
c order deviation (cf CERNLIB)
        do ib =1,16
          idevmin(ib)=iabs(devmin(ib))
        enddo
        call SORTZV(idevmin,stripnum,16,-1,0,0)
c
        num_x2(icircT)=stripnum(1)
        dev_2(icircT)=devmin(stripnum(1))
         if(dev_2(icirct).lt.0)then
           signdev_2(icircT)=1
         elseif(dev_2(icirct).eq.0)then
           signdev_2(icircT)=3   
         else
           signdev_2(icircT)=2
         endif   


c        
c        
c INTERFACE L0-L.U.T. (and L2) in Y : nO CIRCUIT + Hitten STRIP (MC1)  
c
         if(icircT.le.ncirc)then
          nstripy=nstrip_c(icircT)
         else
          nstripy=nstrip_c(icircT-ncirc)
         endif  
c
       do ib=nstripy,1,-1       !keep lower weight (convention)
        if(co_y(ib,icircT).ne.0)num_y2(icircT)=ib
       enddo 
        
         

        if(idebug.eq.6)then 
          if(icircT.le.ncirc)then
           nocirc=NUM(icircT)
          else
           nocirc=-NUM(icircT-ncirc)
          endif  
          print *,' '
          print *,'no circuit ',icircT,nocirc
         do ib=1,16
         print *,'strip no ',ib,':	',i2bin(co_l(ib,icircT),17),
     +   '	dev_min=',devmin(ib)
         enddo 
         print *,'num_x2(icircT)= ',num_x2(icircT)
         print *,' dev_2(icircT)= ',dev_2(icircT)
c
         print *,' ' 
         do ib=1,nstripy
         print *,'strip no ',ib,':	',
     +   i2bin(co_y(ib,icircT),3)
         enddo 
         print *,'num_y2(icircT)= ',num_y2(icircT)
        endif  


       endif                               !end local trigger +/- 8 strips
      enddo 

c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Datas for L0===>L.U.T. (and L2) :  
c                       icircT
c			num_x2(icircT)===>X MC1
c		        dev_2(icircT) ===>X MC2
c		     signdev_2(icircT)===>sign dev_2(icircT)
c 	                num_y2(icircT)===>Y MC1                                             
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Now calculate X MC1 and X MC2 and Y MC1 
c ===> Back to trigger GEOMETRY
c 
      do icircT=1,2*ncirc
      sign_lv=sign_l(icircT)*val_y(icircT) 
       if(sign_lv.ne.0)then                !local trigger +/- 8 strips
calculation of YL2 (middle of the strip)
        if(icircT.le.ncirc)then            !Y>0
         yL2(icircT)=ycmin(1,icircT)+ywi(icircT)*(num_y2(icircT)-.5) 
        else 				   ! Y<0
         yL2(icircT)=-ycmin(1,icircT-ncirc)
     +               -ywi(icircT-ncirc)*(num_y2(icircT)-.5) 
        endif

calculation of X1L2 (MC1) (middle of the strip)           
        if(icircT.le.ncirc)then            ! Y>0
         x1L2(icircT)=xcmin(1,icircT)+xwi_c(icircT)*(num_x2(icircT)-.5) 
        else  		                   ! Y<0
         x1L2(icircT)=xcmin(1,icircT-ncirc)
     +               +xwi_c(icircT-ncirc)*(num_x2(icircT)-.5) 
        endif

calculation of X2L2 (MC2) (middle of the strip)           

          if(icircT.le.ncirc)then
           xxwi_c=xwi_c(icircT)
           xxwi_p=xwi_p(icircT)
           xxwi_m=xwi_m(icircT)
          else
           xxwi_c=xwi_c(icircT-ncirc)
           xxwi_p=xwi_p(icircT-ncirc)
           xxwi_m=xwi_m(icircT-ncirc)
          endif  

         if(dev_2(icircT).ge.0)then                      
          Dnp=16-num_x2(icircT)
          Dnm=0						!RAZ, not used
          dsup=dev_2(icircT)-Dnp
          dinf=0                                        !RAZ, not used
           if(dsup.le.0)then                              
            x2L2(icircT)=(x1L2(icircT)+xxwi_c*dev_2(icircT))*Z3/Z1
           else                                               
            x2L2(icircT)=(x1L2(icircT)+xxwi_c*Dnp+xxwi_p*dsup)*Z3/Z1 
           endif 
         else                                           
          Dnp=0						!RAZ, not used
          Dnm=num_x2(icircT)-1
          dinf=iabs(dev_2(icircT))-Dnm
          dsup=0                                        !RAZ, not used
           if(dinf.le.0)then                              
            x2L2(icircT)=(x1L2(icircT)-xxwi_c*iabs(dev_2(icircT)))*Z3/Z1
           else                                               
            x2L2(icircT)=(x1L2(icircT)-xxwi_c*Dnm-xxwi_m*dinf)*Z3/Z1 
           endif 
         endif

         if(idebug.eq.7)then 
          if(icircT.le.ncirc)then
           nocirc=NUM(icircT)
           yywi=ywi(icircT)
          else
           nocirc=-NUM(icircT-ncirc)
           yywi=ywi(icircT-ncirc)
          endif  
          print *,' '
          print *,'no circuit ',icircT,nocirc
          print *,'num_y2 ywi yL2  = ',num_y2(icircT),yywi,yL2(icircT)                    
          print *,'num_x2 xwiC x1L2 = ',
     +             num_x2(icircT),xxwi_c,x1L2(icircT)                    
          print *,'dev_2 DnP Dnm dsup dinf ',
     +             dev_2(icircT),Dnp,Dnm,dsup,dinf                    
          print *,'xwiM xwiC xwiP x2L2 = ',
     +             xxwi_m,xxwi_c,xxwi_p,x2L2(icircT)                    
          print *,'sign dev =',i2bin(sign_lv,2)
         endif
         
       endif                                  !end local trigger +/- 8 strips
      enddo				      !circuits

c
c
ccc   NOW pt calculation (for the L.U.T)  
c
cinit
      idlc_2p=0        
      idhc_2p=0        
      idlc_2m=0        
      idhc_2m=0        
      idlc_2=0
      idhc_2=0            
c      
      do icircT=1,2*ncirc
       sign_lv=sign_l(icircT)*val_y(icircT) 
       if(sign_lv.ne.0)then                    !local pt trigger +/- 8 strips 
c calculate tetadev (x-z plane)
          anum=z3*(x2L2(icircT)-x1L2(icircT))/(z3-z1)
          anum=anum-x2L2(icircT)
          tetadev=anum/zF                   
c calculate ptcalc (need rF)
	  xF=(x2L2(icircT)-x1L2(icircT))*(z3-zF)/(z3-z1)                    
          xF=x2L2(icircT)-xF    
          yF=yL2(icircT)*zF/z1
          rF=sqrt(xF*xF+yF*yF)
          if(abs(tetadev).ge.0.00001)then
           ptcalc=(rF/zF)*abs(0.9/tetadev)
          else                                             ! infinity
           ptcalc=10000.
          endif 
 
c  for L0 pt cut using the sign of the deviation signdev_2(icircT)

       if (ptcalc.ge.ptcalLow.and.
     +    (signdev_2(icircT).eq.1.or.signdev_2(icircT).eq.3))
     +     idlc_2p=idlc_2p+1
       if (ptcalc.ge.ptcalLow.and.
     +    (signdev_2(icircT).eq.2.or.signdev_2(icircT).eq.3))
     +     idlc_2m=idlc_2m+1
       if (ptcalc.ge.ptcalLow.and.signdev_2(icircT).ge.1)
     +     idlc_2=idlc_2+1
     
       if (ptcalc.ge.ptcalHigh.and.
     +    (signdev_2(icircT).eq.1.or.signdev_2(icircT).eq.3))
     +     idhc_2p=idhc_2p+1
       if (ptcalc.ge.ptcalHigh.and.
     +    (signdev_2(icircT).eq.2.or.signdev_2(icircT).eq.3))
     +     idhc_2m=idhc_2m+1
       if (ptcalc.ge.ptcalHigh.and.signdev_2(icircT).ge.1)
     +     idhc_2=idhc_2+1
       
         if(idebug.eq.8)then 
          if(icircT.le.ncirc)then
           nocirc=NUM(icircT)
          else
           nocirc=-NUM(icircT-ncirc)
          endif  
          print *,' '
          print *,'no circuit ',icircT,nocirc
          print *,'compteurs trigger L0 (boucle circuit ) '
          print *,'signdev_2 =',i2bin(signdev_2(icircT),2)
          print *,'ptcalc,ptcalLow',ptcalc,ptcalLow 
          print *,'LOW  pt calc + et -, mult:',idlc_2p,idlc_2m,idlc_2
          print *,'ptcalc,ptcalHigh',ptcalc,ptcalHigh 
          print *,'HIGH pt calc + et -, mult:',idhc_2p,idhc_2m,idhc_2
         endif
       
                 
       endif			         !end local trigger +/- 8 strips
      enddo				 !circuits

c
c for L0 L.U.T. trigger response
c
       itrigL0=0
       itrigH0=0

      if(idlc_2.ge.2.and.
     + idlc_2p.ne.0.and.idlc_2m.ne.0)itrigL0=1
      if(idhc_2.ge.2.and.
     + idhc_2p.ne.0.and.idhc_2m.ne.0)itrigH0=1
*       if(idlc_2.ge.1)itrigL0=1
*       if(idhc_2.ge.1)itrigH0=1
c
      return
      end
c
c
cc
      subroutine seqtest 

      real x(4,1000),y(4,1000)     ! 1000=nhitmax
      common /hitRPC/ x,y 
c
c generation of test sequence of hits for debuging proposes
c       
       x(1,1)= 58.8
       y(1,1)=-52.3
       x(2,1)= 60.6
       y(2,1)=-52.7
       x(3,1)= 64.8
       y(3,1)=-55.5
       x(4,1)= 65.5
       y(4,1)=-56.0

       x(1,2)=-105.4 
       y(1,2)=-26.7
       x(2,2)=-106.7 
       y(2,2)=-27.0
       x(3,2)=-114.4 
       y(3,2)=-28.9
       x(4,2)=-115.7
       y(4,2)=-29.2

      return 
      end 
c
c 
c
c function that returns a string corresp. to the integer coded in bits
c
      CHARACTER*(*) FUNCTION i2bin(iw,kbit)    
c      i2bin='' ???????????????

      i2bin=' '
      do ib= 0,kbit-1
       if(ibits(iw,ib,1).eq.0)then 
        i2bin='0'//i2bin  
       else 
        i2bin='1'//i2bin  
       endif 
      enddo
      return 
      end   

      SUBROUTINE SORTZV (A,INDEX,N1,MODE,NWAY,NSORT)
C
C CERN PROGLIB# M101    SORTZV          .VERSION KERNFOR  3.15  820113
C ORIG. 02/10/75
C
      DIMENSION A(N1),INDEX(N1)
C
C
      N = N1
      IF (N.LE.0)            RETURN
      IF (NSORT.NE.0) GO TO 2
      DO 1 I=1,N
    1 INDEX(I)=I
C
    2 IF (N.EQ.1)            RETURN
      IF (MODE)    10,20,30
   10 CALL SORTTI (A,INDEX,N)
      GO TO 40
C
   20 CALL SORTTC(A,INDEX,N)
      GO TO 40
C
   30 CALL SORTTF (A,INDEX,N)
C
   40 IF (NWAY.EQ.0) GO TO 50
      N2 = N/2
      DO 41 I=1,N2
      ISWAP = INDEX(I)
      K = N+1-I
      INDEX(I) = INDEX(K)
   41 INDEX(K) = ISWAP
   50 RETURN
      END
*     ========================================
      SUBROUTINE SORTTF (A,INDEX,N1)
C
      DIMENSION A(N1),INDEX(N1)
C
      N = N1
      DO 3 I1=2,N
      I3 = I1
      I33 = INDEX(I3)
      AI = A(I33)
    1 I2 = I3/2
      IF (I2) 3,3,2
    2 I22 = INDEX(I2)
      IF (AI.LE.A (I22)) GO TO 3
      INDEX (I3) = I22
      I3 = I2
      GO TO 1
    3 INDEX (I3) = I33
    4 I3 = INDEX (N)
      INDEX (N) = INDEX (1)
      AI = A(I3)
      N = N-1
      IF (N-1) 12,12,5
    5 I1 = 1
    6 I2 = I1 + I1
      IF (I2.LE.N) I22= INDEX(I2)
      IF (I2-N) 7,9,11
    7 I222 = INDEX (I2+1)
      IF (A(I22)-A(I222)) 8,9,9
    8 I2 = I2+1
      I22 = I222
    9 IF (AI-A(I22)) 10,11,11
   10 INDEX(I1) = I22
      I1 = I2
      GO TO 6
   11 INDEX (I1) = I3
      GO TO 4
   12 INDEX (1) = I3
      RETURN
      END
*     ========================================
      SUBROUTINE SORTTI (A,INDEX,N1)
C
      INTEGER A,AI
      DIMENSION A(N1),INDEX(N1)
C
      N = N1
      DO 3 I1=2,N
      I3 = I1
      I33 = INDEX(I3)
      AI = A(I33)
    1 I2 = I3/2
      IF (I2) 3,3,2
    2 I22 = INDEX(I2)
      IF (AI.LE.A (I22)) GO TO 3
      INDEX (I3) = I22
      I3 = I2
      GO TO 1
    3 INDEX (I3) = I33
    4 I3 = INDEX (N)
      INDEX (N) = INDEX (1)
      AI = A(I3)
      N = N-1
      IF (N-1) 12,12,5
    5 I1 = 1
    6 I2 = I1 + I1
      IF (I2.LE.N) I22= INDEX(I2)
      IF (I2-N) 7,9,11
    7 I222 = INDEX (I2+1)
      IF (A(I22)-A(I222)) 8,9,9
    8 I2 = I2+1
      I22 = I222
    9 IF (AI-A(I22)) 10,11,11
   10 INDEX(I1) = I22
      I1 = I2
      GO TO 6
   11 INDEX (I1) = I3
      GO TO 4
   12 INDEX (1) = I3
      RETURN
      END
*     ========================================
      SUBROUTINE SORTTC (A,INDEX,N1)
C
      INTEGER A,AI
      DIMENSION A(N1),INDEX(N1)
C
      N = N1
      DO 3 I1=2,N
      I3 = I1
      I33 = INDEX(I3)
      AI = A(I33)
    1 I2 = I3/2
      IF (I2) 3,3,2
    2 I22 = INDEX(I2)
      IF(ICMPCH(AI,A(I22)))3,3,21
   21 INDEX (I3) = I22
      I3 = I2
      GO TO 1
    3 INDEX (I3) = I33
    4 I3 = INDEX (N)
      INDEX (N) = INDEX (1)
      AI = A(I3)
      N = N-1
      IF (N-1) 12,12,5
    5 I1 = 1
    6 I2 = I1 + I1
      IF (I2.LE.N) I22= INDEX(I2)
      IF (I2-N) 7,9,11
    7 I222 = INDEX (I2+1)
      IF (ICMPCH(A(I22),A(I222))) 8,9,9
    8 I2 = I2+1
      I22 = I222
    9 IF (ICMPCH(AI,A(I22))) 10,11,11
   10 INDEX(I1) = I22
      I1 = I2
      GO TO 6
   11 INDEX (I1) = I3
      GO TO 4
   12 INDEX (1) = I3
      RETURN
      END
*     ========================================
      FUNCTION ICMPCH(IC1,IC2)
C     FUNCTION TO COMPARE TWO 4 CHARACTER EBCDIC STRINGS - IC1,IC2
C     ICMPCH=-1 IF HEX VALUE OF IC1 IS LESS THAN IC2
C     ICMPCH=0  IF HEX VALUES OF IC1 AND IC2 ARE THE SAME
C     ICMPCH=+1 IF HEX VALUES OF IC1 IS GREATER THAN IC2
      I1=IC1
      I2=IC2
      IF(I1.GE.0.AND.I2.GE.0)GOTO 40
      IF(I1.GE.0)GOTO 60
      IF(I2.GE.0)GOTO 80
      I1=-I1
      I2=-I2
      IF(I1-I2)80,70,60
 40   IF(I1-I2)60,70,80
 60   ICMPCH=-1
      RETURN
 70   ICMPCH=0
      RETURN
 80   ICMPCH=1
      RETURN
      END
      

      SUBROUTINE SBIT (IT,IZW,IZP)
C
C CERN PROGLIB# M421    SBIT            .VERSION KERNFOR  4.23  891215
C MOD. true default 24/2/89, JZ
C
C     This non-ANSI code is a default which may be slow, if so
C     it should be replaced by a machine-specific fast routine

      IF(IAND(IT,1).EQ.0) THEN
         IZW = IBCLR(IZW,IZP-1)
      ELSE
         IZW = IBSET(IZW,IZP-1)
      ENDIF

      END

      SUBROUTINE SBIT1 (IZW,IZP)
C
C CERN PROGLIB# M421    SBIT1           .VERSION KERNFOR  4.23  891215
C MOD. true default 24/2/89, JZ
C
C     This non-ANSI code is a default which may be slow, if so
C     it should be replaced by a machine-specific fast routine

      IZW = IBSET (IZW, IZP-1)

      END

