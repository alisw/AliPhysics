;
; $Id$
;
; $Log$
; Revision 1.2  1996/04/26 12:31:51  cernlib
; Correct comment leader and comment/remove cpp ifdef lines
;
; Revision 1.1.1.1  1996/04/01 15:02:53  mclareni
; Mathlib gen
;
;
; Compile this only with:
;#if (defined(CERNLIB_VAXVMS))&&(!defined(CERNLIB_FORTRAN))
;
 .TITLE  NORRAN
;
;     Normal random number generator for VAX 11-780
;     rewritten from CERN IBM 370 version
;     with the IBM-Randomnumbersequence
;     two parameter version (MCGN,SRGN)
;
;     adapted at Wuppertal by H.Forsbach, March 84
;                     last modification : March 84
;
;     FORTRAN external used : RNORTH
;
;     ----------------------------------------------------------------------
      .PSECT  NORRAN$LOCAL, PIC, CON, REL, LCL, NOSHR, NOEXE, RD, WRT, LONG
MCGN: .LONG   ^D12345
SRGN: .LONG   ^D01073
;
ARGUMENT:       .BLKL   1
;
;     ----------------------------------------------------------------------
      .PSECT  NORRAN$CONST, PIC, CON, REL, LCL, SHR, NOEXE, RD, NOWRT, LONG
;
;     lookup table for NORRAN
;
NTBL: .WORD   ^X0000          ; Floating point : 0.00000
      .WORD   ^X3E80          ; Floating point : 0.06250
      .WORD   ^X3F00          ; Floating point : 0.12500
      .WORD   ^X3F00          ; Floating point : 0.12500
      .WORD   ^X3F40          ; Floating point : 0.18750
      .WORD   ^X3F40          ; Floating point : 0.18750
      .WORD   ^X3F40          ; Floating point : 0.18750
      .WORD   ^X3F40          ; Floating point : 0.18750
      .WORD   ^X3F80          ; Floating point : 0.25000
      .WORD   ^X3F80          ; Floating point : 0.25000
      .WORD   ^X3F80          ; Floating point : 0.25000
      .WORD   ^X3F80          ; Floating point : 0.25000
      .WORD   ^X3F80          ; Floating point : 0.25000
      .WORD   ^X4010          ; Floating point : 0.56250
      .WORD   ^X4020          ; Floating point : 0.62500
      .WORD   ^X4020          ; Floating point : 0.62500
      .WORD   ^X4020          ; Floating point : 0.62500
      .WORD   ^X4020          ; Floating point : 0.62500
      .WORD   ^X4020          ; Floating point : 0.62500
      .WORD   ^X4060          ; Floating point : 0.87500
      .WORD   ^X4060          ; Floating point : 0.87500
      .WORD   ^X4060          ; Floating point : 0.87500
      .WORD   ^X4090          ; Floating point : 1.12500
      .WORD   ^X40B8          ; Floating point : 1.43750
      .WORD   ^X0000          ; Floating point : 0.00000
      .WORD   ^X0000          ; Floating point : 0.00000
      .WORD   ^X0000          ; Floating point : 0.00000
      .WORD   ^X0000          ; Floating point : 0.00000
      .WORD   ^X0000          ; Floating point : 0.00000
      .WORD   ^X3E80          ; Floating point : 0.06250
      .WORD   ^X3E80          ; Floating point : 0.06250
      .WORD   ^X3E80          ; Floating point : 0.06250
      .WORD   ^X3E80          ; Floating point : 0.06250
      .WORD   ^X3E80          ; Floating point : 0.06250
      .WORD   ^X3F00          ; Floating point : 0.12500
      .WORD   ^X3F00          ; Floating point : 0.12500
      .WORD   ^X3F00          ; Floating point : 0.12500
      .WORD   ^X3F00          ; Floating point : 0.12500
      .WORD   ^X3F40          ; Floating point : 0.18750
      .WORD   ^X3F40          ; Floating point : 0.18750
      .WORD   ^X3F80          ; Floating point : 0.25000
      .WORD   ^X3FA0          ; Floating point : 0.31250
      .WORD   ^X3FA0          ; Floating point : 0.31250
      .WORD   ^X3FA0          ; Floating point : 0.31250
      .WORD   ^X3FA0          ; Floating point : 0.31250
      .WORD   ^X3FA0          ; Floating point : 0.31250
      .WORD   ^X3FC0          ; Floating point : 0.37500
      .WORD   ^X3FC0          ; Floating point : 0.37500
      .WORD   ^X3FC0          ; Floating point : 0.37500
      .WORD   ^X3FC0          ; Floating point : 0.37500
      .WORD   ^X3FC0          ; Floating point : 0.37500
      .WORD   ^X3FE0          ; Floating point : 0.43750
      .WORD   ^X3FE0          ; Floating point : 0.43750
      .WORD   ^X3FE0          ; Floating point : 0.43750
      .WORD   ^X3FE0          ; Floating point : 0.43750
      .WORD   ^X3FE0          ; Floating point : 0.43750
      .WORD   ^X4000          ; Floating point : 0.50000
      .WORD   ^X4000          ; Floating point : 0.50000
      .WORD   ^X4000          ; Floating point : 0.50000
      .WORD   ^X4000          ; Floating point : 0.50000
      .WORD   ^X4000          ; Floating point : 0.50000
      .WORD   ^X4010          ; Floating point : 0.56250
      .WORD   ^X4010          ; Floating point : 0.56250
      .WORD   ^X4010          ; Floating point : 0.56250
      .WORD   ^X4010          ; Floating point : 0.56250
      .WORD   ^X4030          ; Floating point : 0.68750
      .WORD   ^X4030          ; Floating point : 0.68750
      .WORD   ^X4030          ; Floating point : 0.68750
      .WORD   ^X4030          ; Floating point : 0.68750
      .WORD   ^X4040          ; Floating point : 0.75000
      .WORD   ^X4040          ; Floating point : 0.75000
      .WORD   ^X4040          ; Floating point : 0.75000
      .WORD   ^X4040          ; Floating point : 0.75000
      .WORD   ^X4050          ; Floating point : 0.81250
      .WORD   ^X4050          ; Floating point : 0.81250
      .WORD   ^X4050          ; Floating point : 0.81250
      .WORD   ^X4050          ; Floating point : 0.81250
      .WORD   ^X4060          ; Floating point : 0.87500
      .WORD   ^X4070          ; Floating point : 0.93750
      .WORD   ^X4070          ; Floating point : 0.93750
      .WORD   ^X4070          ; Floating point : 0.93750
      .WORD   ^X4080          ; Floating point : 1.00000
      .WORD   ^X4080          ; Floating point : 1.00000
      .WORD   ^X4080          ; Floating point : 1.00000
      .WORD   ^X4088          ; Floating point : 1.06250
      .WORD   ^X4088          ; Floating point : 1.06250
      .WORD   ^X4088          ; Floating point : 1.06250
      .WORD   ^X4090          ; Floating point : 1.12500
      .WORD   ^X4090          ; Floating point : 1.12500
      .WORD   ^X4098          ; Floating point : 1.18750
      .WORD   ^X4098          ; Floating point : 1.18750
      .WORD   ^X40A0          ; Floating point : 1.25000
      .WORD   ^X40A0          ; Floating point : 1.25000
      .WORD   ^X40A8          ; Floating point : 1.31250
      .WORD   ^X40A8          ; Floating point : 1.31250
      .WORD   ^X40B0          ; Floating point : 1.37500
      .WORD   ^X40B0          ; Floating point : 1.37500
      .WORD   ^X40B8          ; Floating point : 1.43750
      .WORD   ^X40C0          ; Floating point : 1.50000
      .WORD   ^X40C8          ; Floating point : 1.56250
      .WORD   ^X40D0          ; Floating point : 1.62500
      .WORD   ^X40D8          ; Floating point : 1.68750
      .WORD   ^X40E0          ; Floating point : 1.75000
      .WORD   ^X40E8          ; Floating point : 1.81250
      .WORD   ^X3FA0          ; Floating point : 0.31250
      .WORD   ^X3FA0          ; Floating point : 0.31250
      .WORD   ^X3FA0          ; Floating point : 0.31250
      .WORD   ^X3FA0          ; Floating point : 0.31250
      .WORD   ^X3FA0          ; Floating point : 0.31250
      .WORD   ^X3FA0          ; Floating point : 0.31250
      .WORD   ^X3FA0          ; Floating point : 0.31250
      .WORD   ^X3FA0          ; Floating point : 0.31250
      .WORD   ^X3FA0          ; Floating point : 0.31250
      .WORD   ^X3FA0          ; Floating point : 0.31250
      .WORD   ^X3FC0          ; Floating point : 0.37500
      .WORD   ^X3FC0          ; Floating point : 0.37500
      .WORD   ^X3FC0          ; Floating point : 0.37500
      .WORD   ^X3FC0          ; Floating point : 0.37500
      .WORD   ^X3FC0          ; Floating point : 0.37500
      .WORD   ^X3FC0          ; Floating point : 0.37500
      .WORD   ^X3FC0          ; Floating point : 0.37500
      .WORD   ^X3FE0          ; Floating point : 0.43750
      .WORD   ^X3FE0          ; Floating point : 0.43750
      .WORD   ^X3FE0          ; Floating point : 0.43750
      .WORD   ^X3FE0          ; Floating point : 0.43750
      .WORD   ^X3FE0          ; Floating point : 0.43750
      .WORD   ^X4000          ; Floating point : 0.50000
      .WORD   ^X4000          ; Floating point : 0.50000
      .WORD   ^X4030          ; Floating point : 0.68750
      .WORD   ^X4030          ; Floating point : 0.68750
      .WORD   ^X4030          ; Floating point : 0.68750
      .WORD   ^X4030          ; Floating point : 0.68750
      .WORD   ^X4030          ; Floating point : 0.68750
      .WORD   ^X4030          ; Floating point : 0.68750
      .WORD   ^X4030          ; Floating point : 0.68750
      .WORD   ^X4030          ; Floating point : 0.68750
      .WORD   ^X4030          ; Floating point : 0.68750
      .WORD   ^X4040          ; Floating point : 0.75000
      .WORD   ^X4040          ; Floating point : 0.75000
      .WORD   ^X4040          ; Floating point : 0.75000
      .WORD   ^X4040          ; Floating point : 0.75000
      .WORD   ^X4040          ; Floating point : 0.75000
      .WORD   ^X4050          ; Floating point : 0.81250
      .WORD   ^X4070          ; Floating point : 0.93750
      .WORD   ^X4070          ; Floating point : 0.93750
      .WORD   ^X4070          ; Floating point : 0.93750
      .WORD   ^X4070          ; Floating point : 0.93750
      .WORD   ^X4070          ; Floating point : 0.93750
      .WORD   ^X4070          ; Floating point : 0.93750
      .WORD   ^X4070          ; Floating point : 0.93750
      .WORD   ^X4070          ; Floating point : 0.93750
      .WORD   ^X4070          ; Floating point : 0.93750
      .WORD   ^X4070          ; Floating point : 0.93750
      .WORD   ^X4080          ; Floating point : 1.00000
      .WORD   ^X4080          ; Floating point : 1.00000
      .WORD   ^X4080          ; Floating point : 1.00000
      .WORD   ^X4080          ; Floating point : 1.00000
      .WORD   ^X4080          ; Floating point : 1.00000
      .WORD   ^X4080          ; Floating point : 1.00000
      .WORD   ^X4080          ; Floating point : 1.00000
      .WORD   ^X4088          ; Floating point : 1.06250
      .WORD   ^X4088          ; Floating point : 1.06250
      .WORD   ^X4088          ; Floating point : 1.06250
      .WORD   ^X4098          ; Floating point : 1.18750
      .WORD   ^X4098          ; Floating point : 1.18750
      .WORD   ^X4098          ; Floating point : 1.18750
      .WORD   ^X4098          ; Floating point : 1.18750
      .WORD   ^X4098          ; Floating point : 1.18750
      .WORD   ^X4098          ; Floating point : 1.18750
      .WORD   ^X4098          ; Floating point : 1.18750
      .WORD   ^X4098          ; Floating point : 1.18750
      .WORD   ^X4098          ; Floating point : 1.18750
      .WORD   ^X4098          ; Floating point : 1.18750
      .WORD   ^X4098          ; Floating point : 1.18750
      .WORD   ^X4098          ; Floating point : 1.18750
      .WORD   ^X40A0          ; Floating point : 1.25000
      .WORD   ^X40A0          ; Floating point : 1.25000
      .WORD   ^X40A0          ; Floating point : 1.25000
      .WORD   ^X40A0          ; Floating point : 1.25000
      .WORD   ^X40A0          ; Floating point : 1.25000
      .WORD   ^X40A0          ; Floating point : 1.25000
      .WORD   ^X40A0          ; Floating point : 1.25000
      .WORD   ^X40A0          ; Floating point : 1.25000
      .WORD   ^X40A0          ; Floating point : 1.25000
      .WORD   ^X40A8          ; Floating point : 1.31250
      .WORD   ^X40A8          ; Floating point : 1.31250
      .WORD   ^X40A8          ; Floating point : 1.31250
      .WORD   ^X40A8          ; Floating point : 1.31250
      .WORD   ^X40A8          ; Floating point : 1.31250
      .WORD   ^X40B0          ; Floating point : 1.37500
      .WORD   ^X40B0          ; Floating point : 1.37500
      .WORD   ^X40C0          ; Floating point : 1.50000
      .WORD   ^X40C0          ; Floating point : 1.50000
      .WORD   ^X40C0          ; Floating point : 1.50000
      .WORD   ^X40C0          ; Floating point : 1.50000
      .WORD   ^X40C0          ; Floating point : 1.50000
      .WORD   ^X40C0          ; Floating point : 1.50000
      .WORD   ^X40C0          ; Floating point : 1.50000
      .WORD   ^X40C0          ; Floating point : 1.50000
      .WORD   ^X40C0          ; Floating point : 1.50000
      .WORD   ^X40C0          ; Floating point : 1.50000
      .WORD   ^X40C0          ; Floating point : 1.50000
      .WORD   ^X40C0          ; Floating point : 1.50000
      .WORD   ^X40C0          ; Floating point : 1.50000
      .WORD   ^X40C8          ; Floating point : 1.56250
      .WORD   ^X40C8          ; Floating point : 1.56250
      .WORD   ^X40C8          ; Floating point : 1.56250
      .WORD   ^X40C8          ; Floating point : 1.56250
      .WORD   ^X40C8          ; Floating point : 1.56250
      .WORD   ^X40C8          ; Floating point : 1.56250
      .WORD   ^X40C8          ; Floating point : 1.56250
      .WORD   ^X40C8          ; Floating point : 1.56250
      .WORD   ^X40C8          ; Floating point : 1.56250
      .WORD   ^X40C8          ; Floating point : 1.56250
      .WORD   ^X40D0          ; Floating point : 1.62500
      .WORD   ^X40D0          ; Floating point : 1.62500
      .WORD   ^X40D0          ; Floating point : 1.62500
      .WORD   ^X40D0          ; Floating point : 1.62500
      .WORD   ^X40D0          ; Floating point : 1.62500
      .WORD   ^X40D0          ; Floating point : 1.62500
      .WORD   ^X40D0          ; Floating point : 1.62500
      .WORD   ^X40D8          ; Floating point : 1.68750
      .WORD   ^X40D8          ; Floating point : 1.68750
      .WORD   ^X40D8          ; Floating point : 1.68750
      .WORD   ^X40D8          ; Floating point : 1.68750
      .WORD   ^X40D8          ; Floating point : 1.68750
      .WORD   ^X40E0          ; Floating point : 1.75000
      .WORD   ^X40E0          ; Floating point : 1.75000
      .WORD   ^X40F0          ; Floating point : 1.87500
      .WORD   ^X40F0          ; Floating point : 1.87500
      .WORD   ^X40F0          ; Floating point : 1.87500
      .WORD   ^X40F0          ; Floating point : 1.87500
      .WORD   ^X40F0          ; Floating point : 1.87500
      .WORD   ^X40F0          ; Floating point : 1.87500
      .WORD   ^X40F0          ; Floating point : 1.87500
      .WORD   ^X40F0          ; Floating point : 1.87500
      .WORD   ^X40F0          ; Floating point : 1.87500
      .WORD   ^X40F0          ; Floating point : 1.87500
      .WORD   ^X40F0          ; Floating point : 1.87500
      .WORD   ^X40F0          ; Floating point : 1.87500
      .WORD   ^X40F0          ; Floating point : 1.87500
      .WORD   ^X40F0          ; Floating point : 1.87500
      .WORD   ^X40F0          ; Floating point : 1.87500
      .WORD   ^X40F8          ; Floating point : 1.93750
      .WORD   ^X40F8          ; Floating point : 1.93750
      .WORD   ^X40F8          ; Floating point : 1.93750
      .WORD   ^X40F8          ; Floating point : 1.93750
      .WORD   ^X40F8          ; Floating point : 1.93750
      .WORD   ^X40F8          ; Floating point : 1.93750
      .WORD   ^X40F8          ; Floating point : 1.93750
      .WORD   ^X40F8          ; Floating point : 1.93750
      .WORD   ^X40F8          ; Floating point : 1.93750
      .WORD   ^X40F8          ; Floating point : 1.93750
      .WORD   ^X40F8          ; Floating point : 1.93750
      .WORD   ^X40F8          ; Floating point : 1.93750
      .WORD   ^X40F8          ; Floating point : 1.93750
      .WORD   ^X4100          ; Floating point : 2.00000
      .WORD   ^X4100          ; Floating point : 2.00000
      .WORD   ^X4100          ; Floating point : 2.00000
      .WORD   ^X4100          ; Floating point : 2.00000
      .WORD   ^X4100          ; Floating point : 2.00000
      .WORD   ^X4100          ; Floating point : 2.00000
      .WORD   ^X4100          ; Floating point : 2.00000
      .WORD   ^X4100          ; Floating point : 2.00000
      .WORD   ^X4100          ; Floating point : 2.00000
      .WORD   ^X4100          ; Floating point : 2.00000
      .WORD   ^X4100          ; Floating point : 2.00000
      .WORD   ^X4100          ; Floating point : 2.00000
      .WORD   ^X4104          ; Floating point : 2.06250
      .WORD   ^X4104          ; Floating point : 2.06250
      .WORD   ^X4104          ; Floating point : 2.06250
      .WORD   ^X4104          ; Floating point : 2.06250
      .WORD   ^X4104          ; Floating point : 2.06250
      .WORD   ^X4104          ; Floating point : 2.06250
      .WORD   ^X4104          ; Floating point : 2.06250
      .WORD   ^X4104          ; Floating point : 2.06250
      .WORD   ^X4104          ; Floating point : 2.06250
      .WORD   ^X4104          ; Floating point : 2.06250
      .WORD   ^X4108          ; Floating point : 2.12500
      .WORD   ^X4108          ; Floating point : 2.12500
      .WORD   ^X4108          ; Floating point : 2.12500
      .WORD   ^X4108          ; Floating point : 2.12500
      .WORD   ^X4108          ; Floating point : 2.12500
      .WORD   ^X4108          ; Floating point : 2.12500
      .WORD   ^X4108          ; Floating point : 2.12500
      .WORD   ^X4108          ; Floating point : 2.12500
      .WORD   ^X4108          ; Floating point : 2.12500
      .WORD   ^X410C          ; Floating point : 2.18750
      .WORD   ^X410C          ; Floating point : 2.18750
      .WORD   ^X410C          ; Floating point : 2.18750
      .WORD   ^X410C          ; Floating point : 2.18750
      .WORD   ^X410C          ; Floating point : 2.18750
      .WORD   ^X410C          ; Floating point : 2.18750
      .WORD   ^X410C          ; Floating point : 2.18750
      .WORD   ^X410C          ; Floating point : 2.18750
      .WORD   ^X4110          ; Floating point : 2.25000
      .WORD   ^X4110          ; Floating point : 2.25000
      .WORD   ^X4110          ; Floating point : 2.25000
      .WORD   ^X4110          ; Floating point : 2.25000
      .WORD   ^X4110          ; Floating point : 2.25000
      .WORD   ^X4110          ; Floating point : 2.25000
      .WORD   ^X4110          ; Floating point : 2.25000
      .WORD   ^X4114          ; Floating point : 2.31250
      .WORD   ^X4114          ; Floating point : 2.31250
      .WORD   ^X4114          ; Floating point : 2.31250
      .WORD   ^X4114          ; Floating point : 2.31250
      .WORD   ^X4114          ; Floating point : 2.31250
      .WORD   ^X4114          ; Floating point : 2.31250
      .WORD   ^X4118          ; Floating point : 2.37500
      .WORD   ^X4118          ; Floating point : 2.37500
      .WORD   ^X4118          ; Floating point : 2.37500
      .WORD   ^X4118          ; Floating point : 2.37500
      .WORD   ^X4118          ; Floating point : 2.37500
      .WORD   ^X411C          ; Floating point : 2.43750
      .WORD   ^X411C          ; Floating point : 2.43750
      .WORD   ^X411C          ; Floating point : 2.43750
      .WORD   ^X411C          ; Floating point : 2.43750
      .WORD   ^X4120          ; Floating point : 2.50000
      .WORD   ^X4120          ; Floating point : 2.50000
      .WORD   ^X4120          ; Floating point : 2.50000
      .WORD   ^X4124          ; Floating point : 2.56250
      .WORD   ^X4124          ; Floating point : 2.56250
      .WORD   ^X4124          ; Floating point : 2.56250
      .WORD   ^X4128          ; Floating point : 2.62500
      .WORD   ^X4128          ; Floating point : 2.62500
      .WORD   ^X412C          ; Floating point : 2.68750
      .WORD   ^X412C          ; Floating point : 2.68750
;
;     ----------------------------------------------------------------------
;
      .PSECT  NORRAN$CODE, PIC, CON, REL, LCL, SHR, EXE, RD, NOWRT, LONG
;
      .ENTRY  NORRAN, ^M<R2>
;
;     normal random number generator, FORTRAN callable CALL NORRAN (RANDOM)
;
      MOVL    SRGN,R0         ;move SRGN -> R0
      MOVL    R0,R1           ;move R0 -> R1
      EXTZV   #15,#17,R1,R2   ;\
      MOVL    R2,R1           ;_shift right R1 -> R1 (15 bits)
      XORL2   R1,R0           ;exclusive or of R1,R0 -> R0
      MOVL    R0,R1           ;move R0 -> R1
      EXTZV   #0,#15,R1,R2    ;\
      ROTL    #17,R2,R1       ;_shift left R1 -> R1 (17 bits)
      XORL2   R1,R0           ;exclusive or of R1,R0 -> R0
      MOVL    R0,SRGN         ;save the new SRGN
      MOVL    MCGN,R2         ;get MCGN -> R2
      MULL2   #^D69069,R2     ;69069*R2 -> R2
      MOVL    R2,MCGN         ;save new MCGN
      XORL2   R0,R2           ;exclusive or of R0 [SRGN], R1 [MCGN] -> R2
;
;     ----------------------------------------------------------------------
;
      CLRL    R0              ;\
      EXTZV   #24, #8, R2, R0 ; if R2 is greater then
      CMPL    R0, #^X68       ; 68 00 00 00 00 00 go to ND2
      BGEQ    ND2             ;/
;
      EXTZV   #0,#24,R2,R1    ;mantissa into R1
      CVTLF   R1,R2           ;mantissa to VAX-floating
      EXTZV   #7,#5,R2,R1     ;get normalization shift
      ADDL2   #^X64,R1        ;add 128-excess and ajust
      INSV    R1,#7,#8,R2     ;pack exponent into R2
      CLRL    R1              ;\
      MOVW    NTBL[R0], R1    ;- build floating point value from table
      ADDF3   R1, R2, @4(AP)  ;add two randomnumbers and store as argument
      RET
;     ----------------------------------------------------------------------
ND2:  CMPL    R0, #^XD0       ;  if R2 is greater then
      BGEQ    ND3             ;  D0 00 00 00 00 go to ND3
;
      SUBL2   #^X68, R0       ;subtract to fit into table
      EXTZV   #0,#24,R2,R1    ;mantissa into R1
      CVTLF   R1,R2           ;mantissa to VAX-floating
      EXTZV   #7,#5,R2,R1     ;get normalization shift
      ADDL2   #^X64,R1        ;add 128-excess and ajust
      INSV    R1,#7,#8,R2     ;pack exponent into R2
      CLRL    R1              ;\
      MOVW    NTBL[R0], R1    ;- build floating point from NTBL
      ADDF2   R1, R2          ;add table value and random number
      MNEGF   R2, @4(AP)      ;negate it and write onto argument
      RET
;     ----------------------------------------------------------------------
ND3:  CLRL    R0              ;\
      EXTZV   #20, #12, R2, R0; if R2 is greater then
      CMPL    R0, #^XE2F      ; E2 F0 00 00 00 go to ND4
      BGEQ    ND4             ;/
      SUBL2   #^XCE8, R0      ;subtract to fit into table
      EXTZV   #0,#20,R2,R1    ;mantissa into R1
      ROTL    #4, R1, R1      ;shift 4 bits
      CVTLF   R1,R2           ;mantissa to VAX-floating
      EXTZV   #7,#5,R2,R1     ;get normalization shift
      ADDL2   #^X64,R1        ;add 128-excess and ajust
      INSV    R1,#7,#8,R2     ;pack exponent into R2
      CLRL    R1              ;\
      MOVW    NTBL[R0], R1    ;- build floating point value from table
      ADDF3   R1, R2, @4(AP)  ;add and store as argument
      RET
;     ----------------------------------------------------------------------
ND4:  CMPL    R0, #^XF5E      ; if R2 is greater then
      BGEQ    NTTHTL          ; F5 E0 00 00 go to NTTHTL
;
      SUBL2   #^XE17, R0      ;subtract to fit into table
      EXTZV   #0,#20,R2,R1    ;mantissa into R1
      ROTL    #4, R1, R1      ;shift 4 bits
      CVTLF   R1,R2           ;mantissa to VAX-floating
      EXTZV   #7,#5,R2,R1     ;get normalization shift
      ADDL2   #^X64,R1        ;add 128-excess and ajust
      INSV    R1,#7,#8,R2     ;pack exponent into R2
      CLRL    R1              ;\
      MOVW    NTBL[R0], R1    ;- build floating point value from table
      ADDF2   R1, R2          ;add two random numbers
      MNEGF   R2, @4(AP)      ;negate and store as argument
      RET
;     ----------------------------------------------------------------------
NTTHTL:
      MOVL    R2, ARGUMENT    ;save random digits as argument for RNORTH
      PUSHAL  ARGUMENT        ;push address onto stack
      CALLS   #1, G^RNORTH    ;call RNORTH
      MOVL    R0, @4(AP)      ;move result onto argument of NORRAN
      RET
;
;     ============================================================
;
      .ENTRY  UNI, ^M<R2>     ;save reg R2
;
;     UNI : uniform random number generator [0., 1.] for RNORTH
;
      MOVL    SRGN,R0         ;move SRGN -> R0
      MOVL    R0,R1           ;move R0 -> R1
      EXTZV   #15,#17,R1,R2   ;\
      MOVL    R2,R1           ;_shift right R1 -> R1 (15 bits)
      XORL2   R1,R0           ;exclusive or of R1,R0 -> R0
      MOVL    R0,R1           ;move R0 -> R1
      EXTZV   #0,#15,R1,R2    ;\
      ROTL    #17,R2,R1       ;_shift left R1 -> R1 (17 bits)
      XORL2   R1,R0           ;exclusive or of R1,R0 -> R0
      MOVL    R0,SRGN         ;save the new SRGN
      MOVL    MCGN,R2         ;get MCGN -> R2
      MULL2   #^D69069,R2     ;69069*R2 -> R2
      MOVL    R2,MCGN         ;save new MCGN
      XORL2   R0,R2           ;exclusive or of R0 [SRGN], R1 [MCGN] -> R2
;
      EXTZV   #8,#24,R2,R1    ;mantissa into R1
      CVTLF   R1,R2           ;mantissa to VAX-floating
      EXTZV   #7,#5,R2,R1     ;get normalization shift
      ADDL2   #^X68,R1        ;add 128-excess and ajust
      INSV    R1,#7,#8,R2     ;pack exponent into R2
      MOVL    R2,R0           ;copy onto R0 (RNDM2)
      RET
;
;     ======================================================================
;
      .ENTRY  VNI, ^M<R2, R3> ;save reg R2 and R3
;
;     VNI : uniform random number generator [-1., 1.] for RNORTH
;
      MOVL    SRGN,R0         ;move SRGN -> R0
      MOVL    R0,R1           ;move R0 -> R1
      EXTZV   #15,#17,R1,R2   ;\
      MOVL    R2,R1           ;_shift right R1 -> R1 (15 bits)
      XORL2   R1,R0           ;exclusive or of R1,R0 -> R0
      MOVL    R0,R1           ;move R0 -> R1
      EXTZV   #0,#15,R1,R2    ;\
      ROTL    #17,R2,R1       ;_shift left R1 -> R1 (17 bits)
      XORL2   R1,R0           ;exclusive or of R1,R0 -> R0
      MOVL    R0,SRGN         ;save the new SRGN
      MOVL    MCGN,R2         ;get MCGN -> R2
      MULL2   #^D69069,R2     ;69069*R2 -> R2
      MOVL    R2,MCGN         ;save new MCGN
      XORL2   R0,R2           ;exclusive or of R0 [SRGN], R1 [MCGN] -> R2
;
      EXTZV   #7,#24,R2,R1    ;mantissa into R1
      EXTZV   #31,#1,R2,R3    ;sign into R3
      CVTLF   R1,R2           ;mantissa to VAX-floating
      EXTZV   #7,#5,R2,R1     ;get normalization shift
      ADDL2   #^X68,R1        ;add 128-excess and ajust
      INSV    R1,#7,#8,R2     ;pack exponent into R2
      INSV    R3,#15,#1,R2    ;pack sign into R2
      MOVL    R2,R0           ;copy onto R0 (RNDM2)
      RET
;
;     ======================================================================
;
      .ENTRY  NORRIN, ^M<>
;
;     SUBROUTINE NORRIN (MCGN,SRGN)
;     MCGN, SRGN are the starting integers of NORRAN
;
      MOVL    @4(AP),MCGN
      MOVL    @8(AP),SRGN
      RET
;
;     ======================================================================
;
      .ENTRY  NORRUT, ^M<>
;
;     SUBROUTINE NORRUT (MCGN,SRGN)
;     MCGN, SRGN are the last used interges of NORRAN
;
      MOVL    MCGN,@4(AP)
      MOVL    SRGN,@8(AP)
      RET
      .END
