;
; $Id$
;
; $Log$
; Revision 1.2  1996/04/26 12:31:54  cernlib
; Correct comment leader and comment/remove cpp ifdef lines
;
; Revision 1.1.1.1  1996/04/01 15:02:54  mclareni
; Mathlib gen
;
;
; Compile this only with:
; #if (defined(CERNLIB_VAX))&&(!defined(CERNLIB_FORTRAN))
;
 .TITLE  RN32
IY:     .LONG   ^X00010003
RN32::
        .WORD   ^M<>
        MOVL    IY,R0
        MULL2   #69069,R0
        BGTR    OK
        ADDL2   #-2147483648, R0
OK:     MOVL    R0,IY
        BICL2   #^XFF,R0
        CVTLF   R0,R0
        MULF2   #^X3100,R0
        RET
RN32IN::
        .WORD   ^M<IV>
        MOVL    @4(AP),IY
        RET
RN32OT::
        .WORD   ^M<IV>
        MOVL    IY, @4(AP)
        RET
        .END
