;
; $Id$
;
; $Log$
; Revision 1.2  1996/04/26 12:31:53  cernlib
; Correct comment leader and comment/remove cpp ifdef lines
;
; Revision 1.1.1.1  1996/04/01 15:02:54  mclareni
; Mathlib gen
;
;
; Compile this only with:
; #if (defined(CERNLIB_VAX))&&(!defined(CERNLIB_FORTRAN))
;
 .TITLE  NRAN
;
;       SUBROUTINE NRAN (VEC,N)
;       UNIFORM RANDOM NUMBER GENERATOR FOR VAX 11-780
;       REWRITTEN FROM CERN IBM 370 (RNDM) VERSION
;       WITH THE IBM RANDOMNUMBERSEQUENCE.
;       FILLS THE VECTOR VEC WITH N RANDOMNUMBERS
;       ADAPTED AT WUPPERTAL BY H.FORSBACH, JUNE 82
;                       LAST MODIFICATION : JUNE 82
;
MCGN:   .LONG   ^D12345
NRAN::
        .WORD   ^M<R2,R3,R4>    ;SAVE R2,R3,R4
        MOVL    @8(AP),R0       ;GET N
        BGTR    GOOD            ;CHECK N
        RET                     ;RETURN, IF N=0
GOOD:   DECL    R0              ;\
        CLRL    R1              ;_INITIALIZE LOOP
        MOVL    MCGN,R2         ;MOVE MCGN -> R2
LOOP:   MULL2   #^D690069,R2    ;MULTIPLY WITH 690069
        EXTZV   #8,#24,R2,R3    ;MANTISSA INTO R3
        CVTLF   R3,R4           ;MANTISSA TO VAX-FLOATING
        EXTZV   #7,#5,R4,R3     ;GET NORMALIZATION SHIFT
        ADDL2   #^X68,R3        ;ADD 128-EXCESS AND AJUST
        INSV    R3,#7,#8,R4     ;PACK EXPONENT INTO R4
        MOVL    R4,@4(AP)[R1]   ;COPY R4 ONTO ARRAY
        AOBLEQ  R0,R1,LOOP      ;LOOP OVER ARRAY
        MOVL    R2,MCGN         ;FINAL STORE OF NEW MCGN
        RET
NRANUT::
;
;       SUBROUTINE NRANUT (MCGN)
;       MCGN IS THE LAST USED VALUE OF NRAN
;
        .WORD   ^M<>
        MOVL    MCGN,@4(AP)
        RET
NRANIN::
;
;       SUBROUTINE NRANIN (MCGN)
;       MCGN IS THE STARTING VALUE OF NRAN
;
        .WORD   ^M<>
        MOVL    @4(AP),MCGN
        RET
        .END
