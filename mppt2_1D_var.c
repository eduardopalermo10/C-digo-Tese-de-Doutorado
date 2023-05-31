/*----------------------------------------------------------------*/
/* Tese de Doutorado - GEPAE/CTG/UFPE                             */
/* Filename:    GMPCC_PO_6x1.c                                    */
/* Authors:     Eduardo José Barbosa                              */
/* Date:       10.11.2021                                         */
/* Description: MPPT GLOBAL BASEADO EM MODELO DE DIODO ÚNICO      */
/* Status: Final                                                  */
/*----------------------------------------------------------------*/

#define S_FUNCTION_NAME mppt2_1D_var
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include "math.h"
#include "matrix.h"

#ifndef MATLAB_MEX_FILE
#include <brtenv.h>
#endif

#define U(element) (*uPtrs[element])  /* Pointer to Input Port0 */

/* Define input and output widths  */
#define NINPUTS  11
#define NOUTPUTS 9

/* Define number of states  and assign labels*/
#define DSTATES 12

#define IErro_Vout         	                x[0]
#define IErro_Vout0                         x0[0]

#define Gref0                               x[1]
#define Gref00                              x0[1]

#define Gref1                               x[2]
#define Gref10                              x0[2]

#define Gref2                               x[3]
#define Gref20                              x0[3]

#define Gref3                               x[4]
#define Gref30                              x0[4]

#define Gref4                               x[5]
#define Gref40                              x0[5]

#define Gref5                               x[6]
#define Gref50                              x0[6]

#define Vant                                x[7]
#define Vant0                               x0[7]

#define Iant                                x[8]
#define Iant0                               x0[8]

#define T_ant                               x[9]
#define T_ant0                              x0[9]

#define Vref                                x[10]
#define Vref0                               x0[10]

#define contaux                             x[11]
#define contaux0                            x0[11]

/*-------------------------------------------*/

/* Set parameters   */
#define NPARAMS 7
#define DEF_PARAM1(S)   ssGetSFcnParam(S, 0)
#define DEF_PARAM2(S)   ssGetSFcnParam(S, 1)
#define DEF_PARAM3(S)   ssGetSFcnParam(S, 2)
#define DEF_PARAM4(S)   ssGetSFcnParam(S, 3)
#define DEF_PARAM5(S)   ssGetSFcnParam(S, 4)
#define DEF_PARAM6(S)   ssGetSFcnParam(S, 5)
#define DEF_PARAM7(S)   ssGetSFcnParam(S, 6)

/* Constants general	*/
#define q         1.60217657e-19	 // carga elementar
#define K         1.3806488e-23	     // boltzmann
#define A1        1
#define Ns        54
#define ki        0.003251000000000
#define kv        -0.126600000000000
#define Eg_0      1.1557
#define a         7.021e-4
#define b         1108
#define Gr        1000
#define Eg_r      1.1114
#define k_Rs      0.066848496093750
#define k_Rp      1.499999999998725e-05
#define gama_Rs   1.31
#define gama_Rp   0.7
#define a_Voc     0.004842
#define I01_STC   3.4e-10
#define Rp_STC    1.1886e+02
#define Ig_STC    8.21
#define Isc_STC   8.184488654025490
#define Voc_STC   33.096179183135700
#define Rs_T      0.135
#define Rs_S      0.135

/*====================*
 * S-function methods *
 *====================*/

/* Function: mdlInitializeSizes ============================================== */
static void mdlInitializeSizes(SimStruct *S) {
    ssSetNumSFcnParams(S, NPARAMS);  /* Number of expected parameters */
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        return; /* Parameter mismatch will be reported by Simulink */
    }
    
    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, DSTATES);
    
    if (!ssSetNumInputPorts(S, 1)) return;
    ssSetInputPortWidth(S, 0, NINPUTS);
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    
    if (!ssSetNumOutputPorts(S, 1)) return;
    ssSetOutputPortWidth(S, 0, NOUTPUTS);
    
    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);
    
    /* Take care when specifying exception free code - see sfuntmpl_doc.c */
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}

/* Function: mdlInitializeSampleTimes =========================================*/
static void mdlInitializeSampleTimes(SimStruct *S) {
    ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_INITIALIZE_CONDITIONS
/* Function: mdlInitializeConditions ========================================*/
static void mdlInitializeConditions(SimStruct *S) {
    real_T *x0 = ssGetRealDiscStates(S);
    
    IErro_Vout0  = 0;
    Gref00       = 100;
    Gref10       = 100;
    Gref20       = 100;
    Gref30       = 100;
    Gref40       = 100;
    Gref50       = 100;
    Vant0        = 5;
    Iant0        = 1;
    T_ant0       = 25;
    Vref0        = 25;
    contaux0     = 0;
}

/* Function: mdlOutputs =======================================================*/
static void mdlOutputs(SimStruct *S, int_T tid) {
    real_T Tr   		          = *mxGetPr(DEF_PARAM1(S));
    real_T Ts		              = *mxGetPr(DEF_PARAM2(S));
    real_T Iout_satup		      = *mxGetPr(DEF_PARAM3(S));
    real_T Iout_satdown		      = *mxGetPr(DEF_PARAM4(S));
    real_T TdoC		              = *mxGetPr(DEF_PARAM5(S));
    real_T Kp		              = *mxGetPr(DEF_PARAM6(S));
    real_T Kint		              = *mxGetPr(DEF_PARAM7(S));
    
    real_T            *y      = ssGetOutputPortRealSignal(S, 0);
    real_T            *x      = ssGetRealDiscStates(S);
    InputRealPtrsType uPtrs   = ssGetInputPortRealSignalPtrs(S, 0);
    
    //Parâmetros dos Módulos
    float   Vt;
    float   Eg_T;
    float   I01;
    float   Rs[6];
    float   Rp[6];
    float   Ig[6];
    float   Isc[6];
    float   Voc[6];
    
    //Condições de Operação dos Módulos
    float   G[6];
    float   T;

    //Variáveis do MPCC
    float   Gaux;
    float   Vaux;
    float   Iaux;
    float   aux;
    float   Pot_aux;
    float   Pot_aux2;
    float   F;
    float   dFdI;
    float   dI;
    float   dFdV;
    float   dV;
    float   Pot_ref;
    float   Vmpp[6];
    float   Impp[6];
    float   V[6];
    float   Vtot[6];
    float   Pot[6];
    float   Var_T;
    
    //Contadores
    int_T   i;
    int_T   j;
    int_T   k;
    
    //Variáveis do Controlador
    float   V_barr;
    float   Istring;
    float   ativa_INC;
    float   Erro_Vout;
    float   D;
    float   NPS;
    float   enable_control;
    float   Dif_V;
    float   Dif_I;
    float   delta_V = 0.1;
    
//* Entradas */
   
    enable_control  = U(0);
    T               = U(1);
    G[0]            = U(2);
    G[1]            = U(3);
    G[2]            = U(4);
    G[3]            = U(5);
    G[4]            = U(6);
    G[5]            = U(7);
    V_barr          = U(8);
    Istring         = U(9);
    NPS             = U(10);
   
    if (enable_control >= 1) {
        
        T = T+273;
        Var_T = T-T_ant;
        Eg_T = Eg_0 - (a*(pow(T, 2)))/(T+b);
        Vt = NPS*(K*T*Ns)/q;
        I01 = I01_STC*exp(-(q*(Eg_T/T-Eg_r/Tr))/(A1*K))*pow(T/Tr, 3.0/A1);
        //ORGANIZAÇÃO DAS IRRADIÂNCIAS DO STRING EM ORDEM DECRESCENTE
        for (i = 0; i < 6; i++) {
            Gaux = G[i];
            for (j = i+1; j < 6; j++) {
                if (G[j] > G[i]) {
                    G[i] = G[j];
                    G[j] = Gaux;
                    Gaux = G[i];
                }
            }
        }
        //CÁLCULO DOS PARÂMETROS DE CADA MÓDULO DO STRING
        for (i = 0; i < 6; i++) {
            Rs[i] = NPS*(Rs_T*(1 + k_Rs*(T-Tr)) +  Rs_S*pow((G[i]/1000),(-gama_Rs)));
            Rp[i] = NPS*(Rp_STC*(1 - k_Rp*(T-Tr))*pow((G[i]/1000),(-gama_Rp)));
            Voc[i] = NPS*(Voc_STC + kv*(T-Tr) + a_Voc*T*log(G[i]/1000));
            Isc[i] = (Isc_STC+ki*(T-Tr))*(G[i]/1000);
            Ig[i] = Isc[i]*(1+(Rs[i]/Rp[i])) + I01*(exp((Isc[i]*Rs[i])/(A1*Vt))-1);
        }
        
        if ((fabs(G[0] - Gref0) >= 10) || (fabs(G[1] - Gref1) >= 10) || (fabs(G[2] - Gref2) >= 10) || (fabs(G[3] - Gref3) >= 10) || (fabs(G[4] - Gref4) >= 10) || (fabs(G[5] - Gref5) >= 10) || (Var_T >= 2) || (Var_T <= -2)) {
            
            //REALIZAÇÃO DO MPCC CONVENCIONAL PARA CADA MÓDULO DO STRING, VARRENDO ENTRE 0.8*Voc A Voc
            for (i = 0; i < 6; i++) {
                Iaux = Isc[i]; // corrente de curto de cada módulo
                Pot_aux = 0;
                Vaux = 0.8*Voc[i];
                while (Vaux < Voc[i]) { // Varre até Voc de cada módulo.
                    for (j = 0; j < 100; j++) {
                        aux = Vaux + Iaux*Rs[i]; // variavél auxiliar utilizada para várias contas;
                        F = Ig[i] - aux/Rp[i] - I01*(exp(aux/(A1*Vt)) - 1) - Iaux; // função objetivo da equação (21) tese de aguinaldo;
                        dFdI = -Rs[i]/Rp[i] - Rs[i]*I01/(A1*Vt)*exp(aux/(A1*Vt)) - 1;// derivada em relação a I;
                        dI = -F/dFdI;
                        Iaux = Iaux + dI;// resolver newton raphson;
                        if (fabs(dI) < 1e-6) {// critério de parada - erro
                            if (Iaux < 0) {
                                Iaux = 0;
                            }
                            if (Iaux > Isc[i]) {
                                Iaux = Isc[i];
                            }
                            break;
                        }
                    }
                    Pot_aux2 = Vaux*Iaux;
                    if (Pot_aux2 > Pot_aux) {
                        Pot_aux = Pot_aux2;
                        Vmpp[i] = Vaux;// salva v máximo de cada módulo
                        Impp[i] = Iaux;// salva i máximo de cada módulo
                    }
                    else {
                        break; // Novo critério de parada. Com esta parada o esforço computacional da técnica é reduzido em 7 vezes. 
                    }          // O objetivo desta parada é evitar que o algoritmo "percorra" pontos da curva desnecessários, isto é, pontos além do ponto de máxima potência de cada um dos módulos 
                    Vaux = Vaux + 0.5;
                    contaux = contaux + 1;
                }
            }
            //CÁLCULO DA TENSÃO TOTAL NO STRING PARA CADA UM DOS Vmpp OBTIDOS ANTERIORMENTE
            //quando um dos módulos está sombreado a corrente que passa nele é inferior a corrente do arranjo;
            // o objetivo é calcular as tensões de cada módulo quando a corrente de cada um dos módulos (calculadas na etapa anterior) é a corrente de MPP.
            for (i = 0; i < 6; i++) {
                for (j = 0; j < 6; j++) {
                    if (Impp[i] < Isc[j]) {// verifica, a partir da corrente de máxima potência e da corrente de curto circuito calculada pelo modelo, quantos paineis estão sombreados;
                        V[j] = 0.9*Voc[j];// ponto inicial para garantir convergência;
                        if (i != j) {
                            for (k = 0; k < 100; k++) {
                                aux = V[j] + Impp[i]*Rs[j];
                                F = Ig[j] - aux/Rp[j] - I01*(exp(aux/(A1*Vt)) - 1) - Impp[i];//eq (21);
                                dFdV = -1/Rp[j] - I01/(A1*Vt)*exp(aux/(A1*Vt));// derivada em relação a V;
                                dV = -F/dFdV;
                                V[j] = V[j] + dV;// NR
                                if (fabs(dV) < 1e-6) {// critério de parada;
                                    break;
                                }
                            }
                        }
                        else {
                            V[j] = Vmpp[j];
                        }
                    }
                    else {
                        V[j] = -0.7; // tensão no diodo de desvio;
                    }
                    Vtot[i] = Vtot[i] + V[j]; // tensão total
                }
                Pot[i] = Vtot[i]*Impp[i];
            }
            //SELEÇÃO DA TENSÃO TOTAL QUE GERA A MAIOR POTÊNCIA
            Pot_ref = 0;
            for (i = 0; i < 6; i++) {
                if (Pot[i] > Pot_ref) {
                    Pot_ref = Pot[i];
                    Vref = Vtot[i];
                }
            }
            
            Gref0 = G[0];
            Gref1 = G[1];
            Gref2 = G[2];
            Gref3 = G[3];
            Gref4 = G[4];
            Gref5 = G[5];
            T_ant = T;
            ativa_INC = 0;
        }
        
         // Condição para ativação de Indutância incremental
        if (fabs(V_barr - Vant) <= 0.1) { //Se a tensão medida no ciclo atual for igual à medida no ciclo anterior, a tensão chegou ao valor de referência e o algoritmo INC pode começar
            ativa_INC = 1;
        }
        
//      %%%%%%% Final do MPPC %%%%%%%%%
         
        if (ativa_INC == 1) {
//             contaux = 0;
            Dif_V = V_barr - Vant;
            Dif_I = Istring - Iant;
            if (Dif_V == 0) {
                if (Dif_I > 0) {
                    Vref = Vref + delta_V;
                }
                if (Dif_I < 0) {
                    Vref = Vref - delta_V;
                }
            }
            else {
                if (Istring > -Dif_I/Dif_V*V_barr) {
                    Vref = Vref + delta_V;
                }
                if (Istring < -Dif_I/Dif_V*V_barr) {
                    Vref = Vref - delta_V;
                }
            }
            
//             else {
//                 if (Dif_I/Dif_V > -Istring/V_barr) {
//                     Vref = Vref + delta_V;
//                 }
//                 if (Dif_I/Dif_V < -Istring/V_barr) {
//                     Vref = Vref - delta_V;
//                 }
//             }
            ativa_INC = 0;
        }
        Vant = V_barr;
        Iant = Istring;
        
        //COMEÇO DO CONTROLADOR DA TENSÃO DO BARRAMENTO
        Erro_Vout = V_barr - Vref;    /* Negativo do Erro de Tensão */
        IErro_Vout = IErro_Vout + Erro_Vout*Ts;
        D = Kp*Erro_Vout + Kint*IErro_Vout;
        
        if (D > Iout_satup) {
            D = Iout_satup;
            IErro_Vout = IErro_Vout - Erro_Vout*Ts;
        }
        if (D < Iout_satdown) {
            D = Iout_satdown;
            IErro_Vout = IErro_Vout - Erro_Vout*Ts;
        }
        //* Saidas */
        
        y[0] = Pot_ref;
        y[1] = Vref;
        y[2] = D;
        y[3] = contaux;//Pot[0];
        y[4] = Rs[1];//Pot[1];
        y[5] = Rs[2];//Pot[2];
        y[6] = Rs[3];//Pot[3];
        y[7] = Rs[4];//Pot[4];
        y[8] = Rs[5];//Pot[5];
        
     }
}

#define MDL_UPDATE
/* Function: mdlUpdate ====================================================== */
static void mdlUpdate(SimStruct *S, int_T tid) {
    real_T            *x       = ssGetRealDiscStates(S);
    InputRealPtrsType uPtrs    = ssGetInputPortRealSignalPtrs(S, 0);
    
}

/* Function: mdlTerminate ===================================================== */
static void mdlTerminate(SimStruct *S) {
    UNUSED_ARG(S); /* unused input argument */
}

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif