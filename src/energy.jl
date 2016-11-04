
# DATA FROM:

# Xia, T., SantaLucia, J., Jr., Burkard, M.E., Kierzek, R., Schroeder, S.J., Jiao, X., Cox, C. and Turner, D.H. (1998) Thermodynamic parameters for an expanded nearest-neighbor model for formation of RNA duplexes with Watson-Crick pairs. Biochemistry, 37, 14719-14735.
#
# Mathews, D.H., Disney, M.D., Childs, J.L., Schroeder, S.J., Zuker, M. and Turner, D.H. (1998) Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of RNA secondary structure. Proc. Natl. Acad. Sci. USA, 101, 7287-7292.

const TURNER_1998_STACK = fill( 0.0, (16,16) )

      TURNER_1998_STACK[index(AU_PAIR),index(AU_PAIR)] = -0.9
      TURNER_1998_STACK[index(AU_PAIR),index(UA_PAIR)] = -1.1
      TURNER_1998_STACK[index(AU_PAIR),index(CG_PAIR)] = -2.2
      TURNER_1998_STACK[index(AU_PAIR),index(GC_PAIR)] = -2.1
      TURNER_1998_STACK[index(AU_PAIR),index(GU_PAIR)] = -0.6
      TURNER_1998_STACK[index(AU_PAIR),index(UG_PAIR)] = -1.4

      TURNER_1998_STACK[index(CG_PAIR),index(AU_PAIR)] = -2.1
      TURNER_1998_STACK[index(CG_PAIR),index(UA_PAIR)] = -2.1
      TURNER_1998_STACK[index(CG_PAIR),index(CG_PAIR)] = -3.3
      TURNER_1998_STACK[index(CG_PAIR),index(GC_PAIR)] = -2.4
      TURNER_1998_STACK[index(CG_PAIR),index(GU_PAIR)] = -1.4
      TURNER_1998_STACK[index(CG_PAIR),index(UG_PAIR)] = -2.1

      TURNER_1998_STACK[index(GC_PAIR),index(AU_PAIR)] = -2.4
      TURNER_1998_STACK[index(GC_PAIR),index(UA_PAIR)] = -2.2
      TURNER_1998_STACK[index(GC_PAIR),index(CG_PAIR)] = -3.4
      TURNER_1998_STACK[index(GC_PAIR),index(GC_PAIR)] = -3.3
      TURNER_1998_STACK[index(GC_PAIR),index(GU_PAIR)] = -1.5
      TURNER_1998_STACK[index(GC_PAIR),index(UG_PAIR)] = -2.5

      TURNER_1998_STACK[index(GU_PAIR),index(AU_PAIR)] = -1.3
      TURNER_1998_STACK[index(GU_PAIR),index(UA_PAIR)] = -1.4
      TURNER_1998_STACK[index(GU_PAIR),index(CG_PAIR)] = -2.5
      TURNER_1998_STACK[index(GU_PAIR),index(GC_PAIR)] = -2.1
      TURNER_1998_STACK[index(GU_PAIR),index(GU_PAIR)] = -0.5 # +0.5
      TURNER_1998_STACK[index(GU_PAIR),index(UG_PAIR)] = +1.3

      TURNER_1998_STACK[index(UA_PAIR),index(AU_PAIR)] = -1.3
      TURNER_1998_STACK[index(UA_PAIR),index(UA_PAIR)] = -0.9
      TURNER_1998_STACK[index(UA_PAIR),index(CG_PAIR)] = -2.4
      TURNER_1998_STACK[index(UA_PAIR),index(GC_PAIR)] = -2.1
      TURNER_1998_STACK[index(UA_PAIR),index(GU_PAIR)] = -1.0
      TURNER_1998_STACK[index(UA_PAIR),index(UG_PAIR)] = -1.3

      TURNER_1998_STACK[index(UG_PAIR),index(AU_PAIR)] = -1.0
      TURNER_1998_STACK[index(UG_PAIR),index(UA_PAIR)] = -0.6
      TURNER_1998_STACK[index(UG_PAIR),index(CG_PAIR)] = -1.5
      TURNER_1998_STACK[index(UG_PAIR),index(GC_PAIR)] = -1.4
      TURNER_1998_STACK[index(UG_PAIR),index(GU_PAIR)] = +0.3
      TURNER_1998_STACK[index(UG_PAIR),index(UG_PAIR)] = -0.5


