#ifndef CERNLIB_ISAJET_LISTSS_INC
#define CERNLIB_ISAJET_LISTSS_INC
*
*
* listss.inc
*
C          LISTSS IDENT and JETTYPE codes
C       ISGL  ISUPL -ISUPL  ISDNL -ISDNL  ISSTL -ISSTL  ISCHL -ISCHL
C          1      2      3      4      5      6      7      8      9
C      ISBT1 -ISBT1  ISTP1 -ISTP1  ISUPR -ISUPR  ISDNR -ISDNR  ISSTR
C         10     11     12     13     14     15     16     17     18
C     -ISSTR  ISCHR -ISCHR  ISBT2 -ISBT2  ISTP2 -ISTP2   ISW1  -ISW1
C         19     20     21     22     23     24     25     26     27
C       ISW2  -ISW2   ISZ1   ISZ2   ISZ3   ISZ4  ISNEL -ISNEL   ISEL
C         28     29     30     31     32     33     34     35     36
C      -ISEL  ISNML -ISNML  ISMUL -ISMUL  ISNTL -ISNTL ISTAU1-ISTAU1
C         37     38     39     40     41     42     43     44     45
C       ISER  -ISER  ISMUR -ISMUR ISTAU2-ISTAU2      9      1     -1
C         46     47     48     49     50     51     52     53     54
C          2     -2      3     -3      4     -4      5     -5      6
C         55     56     57     58     59     60     61     62     63
C         -6     11    -11     12    -12     13    -13     14    -14
C         64     65     66     67     68     69     70     71     72
C         15    -15     16    -16     10     80    -80     90   ISHL
C         73     74     75     76     77     78     79     80     81
C       ISHH   ISHA   ISHC  -ISHC
C         82     83     84     85
      COMMON/LISTSS/LISTSS(85)
      INTEGER LISTSS
      SAVE /LISTSS/

#endif
