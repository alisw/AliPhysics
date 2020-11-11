#ifndef CORRELATION_CALCULATOR_H
#define CORRELATION_CALCULATOR_H

#include <iostream>
#include "TComplex.h"

class CorrelationCalculator {

    public:
    
    static const int MaxHarm = 20;
    static const int MaxPow = 20; 

    TComplex Qvector[MaxHarm][MaxPow];
    TComplex Qvector0M[MaxHarm][MaxPow];
    TComplex Qvector0P[MaxHarm][MaxPow];
    TComplex Qvector2M[MaxHarm][MaxPow];
    TComplex Qvector2P[MaxHarm][MaxPow];
    TComplex Qvector4M[MaxHarm][MaxPow];
    TComplex Qvector4P[MaxHarm][MaxPow];
    TComplex Qvector6M[MaxHarm][MaxPow];
    TComplex Qvector6P[MaxHarm][MaxPow];
    TComplex Qvector8M[MaxHarm][MaxPow];
    TComplex Qvector8P[MaxHarm][MaxPow];
    TComplex Qvector10M[MaxHarm][MaxPow];
    TComplex Qvector10P[MaxHarm][MaxPow];
    TComplex Qvector14M[MaxHarm][MaxPow];
    TComplex Qvector14P[MaxHarm][MaxPow];
    TComplex Qvector16M[MaxHarm][MaxPow];
    TComplex Qvector16P[MaxHarm][MaxPow];
    TComplex Qvector18M[MaxHarm][MaxPow];
    TComplex Qvector18P[MaxHarm][MaxPow];
    TComplex QvectorSubLeft[MaxHarm][MaxPow];
    TComplex QvectorSubMiddle[MaxHarm][MaxPow];
    TComplex QvectorSubRight[MaxHarm][MaxPow];
    TComplex QvectorSubGap2Left[MaxHarm][MaxPow];
    TComplex QvectorSubGap2Middle[MaxHarm][MaxPow];
    TComplex QvectorSubGap2Right[MaxHarm][MaxPow];
    TComplex pvector[MaxHarm][MaxPow];
    TComplex pvectorM[MaxHarm][MaxPow];
    TComplex pvectorP[MaxHarm][MaxPow];
    TComplex qvector[MaxHarm][MaxPow];
    TComplex pvector0M[MaxHarm][MaxPow];
    TComplex pvector0P[MaxHarm][MaxPow];
    TComplex pvector4M[MaxHarm][MaxPow];
    TComplex pvector4P[MaxHarm][MaxPow];
    TComplex pvector8M[MaxHarm][MaxPow];
    TComplex pvector8P[MaxHarm][MaxPow];

    TComplex Two(int n1, int n2);
    TComplex TwoGap0(int n1, int n2);
    TComplex TwoGap2(int n1, int n2);
    TComplex TwoGap4(int n1, int n2);
    TComplex TwoGap6(int n1, int n2);
    TComplex TwoGap8(int n1, int n2);
    TComplex TwoGap10(int n1, int n2);
    TComplex TwoGap14(int n1, int n2);
    TComplex TwoGap16(int n1, int n2);
    TComplex TwoGap18(int n1, int n2);
    TComplex Two_3SubLM(int n1, int n2);
    TComplex Two_3SubRM(int n1, int n2);
    TComplex Two_3SubLR(int n1, int n2);
    TComplex Two_3SubGap2LM(int n1, int n2);
    TComplex Two_3SubGap2RM(int n1, int n2);
    TComplex TwoGap0M(int n1, int n2);
    TComplex TwoGap2M(int n1, int n2);
    TComplex TwoGap4M(int n1, int n2);
    TComplex TwoGap6M(int n1, int n2);
    TComplex TwoGap8M(int n1, int n2);
    TComplex TwoGap10M(int n1, int n2);
    TComplex TwoGap0P(int n1, int n2);
    TComplex TwoGap2P(int n1, int n2);
    TComplex TwoGap4P(int n1, int n2);
    TComplex TwoGap6P(int n1, int n2);
    TComplex TwoGap8P(int n1, int n2);
    TComplex TwoGap10P(int n1, int n2);

    TComplex Three(int n1, int n2, int n3);
    TComplex ThreeGap0A(int n1, int n2, int n3);
    TComplex ThreeGap0B(int n1, int n2, int n3);
    TComplex ThreeGap2A(int n1, int n2, int n3);
    TComplex ThreeGap2B(int n1, int n2, int n3);
    TComplex ThreeGap4A(int n1, int n2, int n3);
    TComplex ThreeGap4B(int n1, int n2, int n3);
    TComplex ThreeGap6A(int n1, int n2, int n3);
    TComplex ThreeGap6B(int n1, int n2, int n3);
    TComplex ThreeGap8A(int n1, int n2, int n3);
    TComplex ThreeGap8B(int n1, int n2, int n3);
    TComplex ThreeGap10A(int n1, int n2, int n3);
    TComplex ThreeGap10B(int n1, int n2, int n3);

    TComplex ThreeGap0_subM(int n1, int n2, int n3);
    TComplex ThreeGap0_subP(int n1, int n2, int n3);
    TComplex ThreeGap2_subM(int n1, int n2, int n3);
    TComplex ThreeGap2_subP(int n1, int n2, int n3);
    TComplex ThreeGap4_subM(int n1, int n2, int n3);
    TComplex ThreeGap4_subP(int n1, int n2, int n3);
    TComplex ThreeGap6_subM(int n1, int n2, int n3);
    TComplex ThreeGap6_subP(int n1, int n2, int n3);
    TComplex ThreeGap8_subM(int n1, int n2, int n3);
    TComplex ThreeGap8_subP(int n1, int n2, int n3);
    TComplex ThreeGap10_subM(int n1, int n2, int n3);
    TComplex ThreeGap10_subP(int n1, int n2, int n3);
    TComplex Four(int n1, int n2, int n3, int n4);
    TComplex FourGap0M(int n1, int n2, int n3, int n4);
    TComplex FourGap0P(int n1, int n2, int n3, int n4);
    TComplex FourGap0(int n1, int n2, int n3, int n4);
    TComplex FourGap2(int n1, int n2, int n3, int n4);
    TComplex FourGap4(int n1, int n2, int n3, int n4);
    TComplex FourGap6(int n1, int n2, int n3, int n4);
    TComplex FourGap8(int n1, int n2, int n3, int n4);
    TComplex FourGap10(int n1, int n2, int n3, int n4);
    TComplex Four_3SubMMLR(int n1, int n2, int n3, int n4);
    TComplex Four_3SubLLMR(int n1, int n2, int n3, int n4);
    TComplex Four_3SubRRML(int n1, int n2, int n3, int n4);
    TComplex Four_3SubGap2Evts(int n1, int n2, int n3, int n4);
    TComplex Five(int n1, int n2, int n3, int n4, int n5);
    TComplex FiveGap0A_2(int n1, int n2, int n3, int n4, int n5);
    TComplex FiveGap0A(int n1, int n2, int n3, int n4, int n5);
    TComplex FiveGap0B(int n1, int n2, int n3, int n4, int n5);
    TComplex FiveGap2A(int n1, int n2, int n3, int n4, int n5);
    TComplex FiveGap2B(int n1, int n2, int n3, int n4, int n5);
    TComplex FiveGap4A(int n1, int n2, int n3, int n4, int n5);
    TComplex FiveGap4B(int n1, int n2, int n3, int n4, int n5);
    TComplex FiveGap6A(int n1, int n2, int n3, int n4, int n5);
    TComplex FiveGap6B(int n1, int n2, int n3, int n4, int n5);
    TComplex FiveGap8A(int n1, int n2, int n3, int n4, int n5);
    TComplex FiveGap8B(int n1, int n2, int n3, int n4, int n5);
    TComplex FiveGap10A(int n1, int n2, int n3, int n4, int n5);
    TComplex FiveGap10B(int n1, int n2, int n3, int n4, int n5);
    TComplex Six(int n1, int n2, int n3, int n4, int n5, int n6);
    TComplex SixGap0(int n1, int n2, int n3, int n4, int n5, int n6);
    TComplex SixGap2(int n1, int n2, int n3, int n4, int n5, int n6);
    TComplex SixGap4(int n1, int n2, int n3, int n4, int n5, int n6);
    TComplex SixGap6(int n1, int n2, int n3, int n4, int n5, int n6);
    TComplex SixGap8(int n1, int n2, int n3, int n4, int n5, int n6);
    TComplex SixGap10(int n1, int n2, int n3, int n4, int n5, int n6);
    TComplex Seven(int n1, int n2, int n3, int n4, int n5, int n6, int n7);
    TComplex Eight(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8);
    TComplex EightGap0(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8);

   TComplex Q(int n, int p);
    TComplex QGap0M(int n, int p);
    TComplex QGap0P(int n, int p);
    TComplex QGap2M(int n, int p);
    TComplex QGap2P(int n, int p);
    TComplex QGap4M(int n, int p);
    TComplex QGap4P(int n, int p);
    TComplex QGap6M(int n, int p);
    TComplex QGap6P(int n, int p);
    TComplex QGap8M(int n, int p);
    TComplex QGap8P(int n, int p);
    TComplex QGap10M(int n, int p);
    TComplex QGap10P(int n, int p);
    TComplex QGap14M(int n, int p);
    TComplex QGap14P(int n, int p);
    TComplex QGap16M(int n, int p);
    TComplex QGap16P(int n, int p);
    TComplex QGap18M(int n, int p);
    TComplex QGap18P(int n, int p);
    TComplex QsubLeft(int n, int p);
    TComplex QsubMiddle(int n, int p);
    TComplex QsubRight(int n, int p);
    TComplex QsubGap2Left(int n, int p);
    TComplex QsubGap2Middle(int n, int p);
    TComplex QsubGap2Right(int n, int p);
    TComplex p(int n, int p);
    TComplex pGap10M(int n, int p);
    TComplex pGap10P(int n, int p);
    TComplex q(int n, int p);
    TComplex pGap0M(int n, int p);
    TComplex pGap0P(int n, int p);
    TComplex pGap4M(int n, int p);
    TComplex pGap4P(int n, int p);
    TComplex pGap6M(int n, int p);
    TComplex pGap6P(int n, int p);
    TComplex pGap8M(int n, int p);
    TComplex pGap8P(int n, int p);

    void ResetQ(const int nMaxHarm, const int nMaxPow);
    void FillQVector(TComplex _Qvector[MaxHarm][MaxPow], double _Qcos[MaxHarm][MaxPow], double _Qsin[MaxHarm][MaxPow]);

};
#endif
