/// \ingroup STAT/test
/// \brief  test methods of AliPainter
/// Example usage
/*
\code
.L $AliRoot_SRC/STAT/test/AliPainterTest.C+
AliPainterTest();
root.exe -b -q  $AliRoot_SRC/STAT/test/AliPainterTest.C+ | tee AliPainterTest.log
\endcode
*/

#include "AliPainter.h"

void DivideTPadTest();
void SetMultiGraphTimeAxisTest();
void DrawHistogramTest();

void AliPainterTest(){
  DivideTPadTest();
  SetMultiGraphTimeAxisTest();
  DrawHistogramTest();
}
