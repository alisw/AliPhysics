#include "TF1.h"
TF1* tmp_function_ptr = 0;
double tmp_function(double x=0) { 
    return (tmp_function_ptr) ? tmp_function_ptr->Eval(x) : 0; 
    //return 0;
}
