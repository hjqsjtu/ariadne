add_library(ariadne-function OBJECT
    affine.cpp
    polynomial.cpp
    formula.cpp
    procedure.cpp
    function.cpp
    function_model.cpp
    affine_model.cpp
    taylor_model.cpp
    taylor_function.cpp
    c1_taylor_model.cpp
    chebyshev_model.cpp
)

cotire(ariadne-function)
