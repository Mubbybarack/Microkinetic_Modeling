/*-------------------------------------------------------------------
 *
 * Copyright (C) 2012 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *-------------------------------------------------------------------
 *
 * $Id$
 */
 /*------------------------------------------------------------------
 * Brief description of this file:
 *
 * This file contains the code for the user defined likelihood data
 * class and the user defined likelihood routine.
 *-----------------------------------------------------------------*/

#include <gravity_likelihood.h>
#include <cmath>
#include <stdio.h>
#include <fstream>
#include <iomanip>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fit.h>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>
#include <math.h>

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */

/* Problem Constants */
#define NEQ    47              /* number of equations  */
#define YY1   RCONST(1.0)      /* initial y components */
#define YY2    RCONST(0.0)
#define YY3    RCONST(0.0)
#define YY4    RCONST(0.0)
#define YY5    RCONST(0.0)
#define YY6    RCONST(0.0)
#define YY7    RCONST(0.0)
#define YY8    RCONST(0.0)
#define YY9    RCONST(0.0)
#define YY10   RCONST(0.0)
#define YY11   RCONST(0.0)
#define YY12   RCONST(0.0)
#define YY13   RCONST(0.0)
#define YY14   RCONST(0.0)
#define YY15   RCONST(0.0)
#define YY16   RCONST(0.0)
#define YY17   RCONST(0.0)
#define YY18   RCONST(0.0)
#define YY19   RCONST(0.0)
#define YY20   RCONST(0.0)
#define YY21   RCONST(0.0)
#define YY22   RCONST(0.0)
#define YY23   RCONST(0.0)
#define YY24   RCONST(0.0)
#define YY25   RCONST(0.0)
#define YY26   RCONST(0.0)
#define YY27   RCONST(0.0)
#define YY28   RCONST(0.0)
#define YY29   RCONST(0.0)
#define YY30   RCONST(0.0)
#define YY31   RCONST(0.0)
#define YY32   RCONST(0.0)
#define YY33   RCONST(0.0)
#define YY34   RCONST(0.0)
#define YY35   RCONST(0.0)
#define YY36   RCONST(0.0)
#define YY37   RCONST(0.0)
#define YY38   RCONST(0.0)
#define YY39   RCONST(0.0)
#define YY40   RCONST(0.0)
#define YY41   RCONST(0.0)
#define YY42   RCONST(0.0)
#define YY43   RCONST(0.0)
#define YY44   RCONST(0.0)
#define YY45   RCONST(0.0)
#define YY46   RCONST(0.0)
#define YY47   RCONST(0.0)  //I have 46 intermediates and a site, so 47 equations

#define RTOL  RCONST(1.0e-12)   /* scalar relative tolerance            */
#define	ATOL1	RCONST(1.0E-12)
#define	ATOL2	RCONST(1.0E-12)
#define	ATOL3	RCONST(1.0E-12)
#define	ATOL4	RCONST(1.0E-12)
#define	ATOL5	RCONST(1.0E-12)
#define	ATOL6	RCONST(1.0E-12)
#define	ATOL7	RCONST(1.0E-12)
#define	ATOL8	RCONST(1.0E-12)
#define	ATOL9	RCONST(1.0E-12)
#define	ATOL10	RCONST(1.0E-12)
#define	ATOL11	RCONST(1.0E-12)
#define	ATOL12	RCONST(1.0E-12)
#define	ATOL13	RCONST(1.0E-12)
#define	ATOL14	RCONST(1.0E-12)
#define	ATOL15	RCONST(1.0E-12)
#define	ATOL16	RCONST(1.0E-12)
#define	ATOL17	RCONST(1.0E-12)
#define	ATOL18	RCONST(1.0E-12)
#define	ATOL19	RCONST(1.0E-12)
#define	ATOL20	RCONST(1.0E-12)
#define	ATOL21	RCONST(1.0E-12)
#define	ATOL22	RCONST(1.0E-12)
#define	ATOL23	RCONST(1.0E-12)
#define	ATOL24	RCONST(1.0E-12)
#define	ATOL25	RCONST(1.0E-12)
#define	ATOL26	RCONST(1.0E-12)
#define	ATOL27	RCONST(1.0E-12)
#define	ATOL28	RCONST(1.0E-12)
#define	ATOL29	RCONST(1.0E-12)
#define	ATOL30	RCONST(1.0E-12)
#define	ATOL31	RCONST(1.0E-12)
#define	ATOL32	RCONST(1.0E-12)
#define	ATOL33	RCONST(1.0E-12)
#define	ATOL34	RCONST(1.0E-12)
#define	ATOL35	RCONST(1.0E-12)
#define	ATOL36	RCONST(1.0E-12)
#define	ATOL37	RCONST(1.0E-12)
#define	ATOL38	RCONST(1.0E-12)
#define	ATOL39	RCONST(1.0E-12)
#define	ATOL40	RCONST(1.0E-12)
#define	ATOL41	RCONST(1.0E-12)
#define	ATOL42	RCONST(1.0E-12)
#define	ATOL43	RCONST(1.0E-12)
#define	ATOL44	RCONST(1.0E-12)
#define	ATOL45	RCONST(1.0E-12)
#define	ATOL46	RCONST(1.0E-12)
#define	ATOL47	RCONST(1.0E-12)

#define TT0    RCONST(0.0)      /* initial time           */
#define TT1    RCONST(1E-40)      /* first output time      */
#define TMULT RCONST(4.0)     /* output time factor     */
#define NOUT  70               /* number of output times */
double mxsteps = 50000;

double p[176];

double PCH3CH2CH3_beta;
double PCH3CHCH2_beta;
double PH2_beta;
double PCH3CCH_beta;
double PCH3CH3_beta;
double PCH2CH2_beta;
double PCHCH_beta;
double PCH4_beta;

double alph, beta, gamm, delta, epsi, zeta, eta, theta; //GAS PHASE uncertainty modifiers. ALPH = Propane, BETA = PROPYLENE, GAMM = H2 DELTA = PROPYNE, ect
double T_beta;
double h=4.135667516E-15;

double ff[138];                                                                                                                        
double bb[138];                                                                                      
double S[138]; 
double Sp[138];
double TS[138];
double KK[138];
double rr[138];


double	CH3CH2CH3_energy;
double	CH3CHCH3_energy	;
double	CH3CH2CH2_energy;
double	CH3CHCH2_energy	;
double	CH3CH2CH_energy	;
double	CH2CH2CH2_energy;
double	CH3CCH3_energy	;

//C5	and
double	CH3CH2C_energy	;
double	CH2CH2CH_energy	;
double	CH2CHCH2_energy	;
double	CH3CHCH_energy	;
double	CH3CCH2_energy	;
double	CH3CHC_energy	;
double	CH2CH2C_energy	;
double	CHCH2CH_energy	;
double	CH2CHCH_energy	;
double	CH2CCH2_energy	;
double	CH3CCH_energy	;
//C3H3to	C3
double	CH3CC_energy	;
double	CH2CHC_energy	;
double	CHCHCH_energy	;
double	CHCH2C_energy	;
double	CH2CCH_energy	;
double	CH2CC_energy	;
double	CHCHC_energy	;
double	CCH2C_energy	;
double	CHCCH_energy	;
double	CCHC_energy	;
double	CHCC_energy	;
double	CCC_energy	;
//CH3CH3	->
double	CH3CH3_energy	;
double	CH3CH2_energy	;
double	CH3CH_energy	;
double	CH3C_energy	;
double	CH2CH2_energy	;
double	CH2CH_energy	;
double	CH2C_energy	;
double	CHCH_energy	;
double	CHC_energy	;
double	CC_energy	;

double	CH4_energy	;
double	CH3_energy	;
double	CH2_energy	;
double	CH_energy	;
double	C_energy	;
double	H_energy	;


double	dEA9	;
double 	dEA10	;
double	dEA11	;
double 	dEA12	;
double	dEA13	;
double 	dEA14	;
double	dEA15	;
double 	dEA16	;
double	dEA17	;
double 	dEA18	;
double	dEA19	;
double 	dEA20	;
double	dEA21	;
double 	dEA22	;
double	dEA23	;
double 	dEA24	;
double	dEA25	;
double 	dEA26	;
double	dEA27	;
double 	dEA28	;
double	dEA29	;
double 	dEA30	;
double	dEA31	;
double 	dEA32	;
double	dEA33	;
double 	dEA34	;
double	dEA35	;
double 	dEA36	;
double	dEA37	;
double 	dEA38	;
double	dEA39	;
double 	dEA40	;
double	dEA41	;
double 	dEA42	;
double	dEA43	;
double 	dEA44	;
double	dEA45	;
double 	dEA46	;
double	dEA47	;
double 	dEA48	;
double	dEA49	;
double 	dEA50	;
double	dEA51	;
double 	dEA52	;
double	dEA53	;
double 	dEA54	;
double	dEA55	;
double 	dEA56	;
double	dEA57	;
double 	dEA58	;
double	dEA59	;
double 	dEA60	;
double	dEA61	;
double 	dEA62	;
double	dEA63	;
double 	dEA64	;
double	dEA65	;
double 	dEA66	;
double	dEA67	;
double 	dEA68	;
double	dEA69	;
double 	dEA70	;
double	dEA71	;
double 	dEA72	;
double	dEA73	;
double 	dEA74	;
double	dEA75	;
double 	dEA76	;
double	dEA77	;
double 	dEA78	;
double	dEA79	;
double 	dEA80	;
double	dEA81	;
double 	dEA82	;
double	dEA83	;
double 	dEA84	;
double	dEA85	;
double 	dEA86	;
double	dEA87	;
double 	dEA88	;
double	dEA89	;
double 	dEA90	;
double	dEA91	;
double 	dEA92	;
double	dEA93	;
double 	dEA94	;
double	dEA95	;
double 	dEA96	;
double	dEA97	;
double 	dEA98	;
double	dEA99	;
double 	dEA100	;
double	dEA101	;
double 	dEA102	;
double	dEA103	;
double 	dEA104	;
double	dEA105	;
double 	dEA106	;
double	dEA107	;
double 	dEA108	;
double	dEA109	;
double 	dEA110	;
double	dEA111	;
double 	dEA112	;
double	dEA113	;
double 	dEA114	;
double	dEA115	;
double 	dEA116	;
double	dEA117	;
double 	dEA118	;
double	dEA119	;
double 	dEA120	;
double	dEA121	;
double 	dEA122	;
double	dEA123	;
double 	dEA124	;
double	dEA125	;
double 	dEA126	;
double	dEA127	;
double 	dEA128	;
double	dEA129	;
double 	dEA130	;
double	dEA131	;
double 	dEA132	;
double	dEA133	;
double 	dEA134	;
double	dEA135	;
double 	dEA136	;
double	dEA137	;
double 	dEA138	;


double b[138];
double r[138];
double K[138];

double CH3CH2CH3_gp 	= 1.7536 ;
double CH3CHCH2_gp 		= 2.5892 ;
double CH3CCH_gp 		= 4.0418 ;
double CH3CH3_gp 		= 0.7522 ;
double CH2CH2_gp 		= 1.7580 ;
double CHCH_gp 			= 3.8069 ;
double CH4_gp 			= -0.3694 ;
double H2_gp 			= -0.7988 ;


/* Functions Called by the Solver */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

/* Private function to check function return values */

static int check_flag(void *flagvalue, char *funcname, int opt);

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
//  realtype y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, y15, y16, y17,y18, y19, y20, y21, y22, y23, y24, y25, y26, y27, y28, y29, y30, y31, y32, y33, y34, y35, y36, y37, y38, y39, y40, y41, y42, y43, y44, y45, y46, y47,yd1, yd2, yd3, yd4, yd5, yd6, yd7, yd8, yd9, yd10, yd11, yd12, yd13, yd14, yd15, yd16, yd17, yd18, yd19, yd20, yd21, yd22, yd23, yd24, yd25, yd26, yd27, yd28, yd29, yd30, yd31, yd32, yd33, yd34, yd35, yd36, yd37, yd38, yd39, yd40, yd41, yd42, yd43, yd44, yd45, yd46, yd47;

 realtype y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, y15, y16, y17,y18, y19, y20, y21, y22, y23, y24, y25, y26, y27, y28, y29, y30, y31, y32, y33, y34, y35, y36, y37, y38, y39, y40, y41, y42, y43, y44, y45, y46, y47, y48, y49, y50 ,y51,y52,y53,y54,y55,yd1, yd2, yd3, yd4, yd5, yd6, yd7, yd8, yd9, yd10, yd11, yd12, yd13, yd14, yd15, yd16, yd17, yd18, yd19, yd20, yd21, yd22, yd23, yd24, yd25, yd26, yd27, yd28, yd29, yd30, yd31, yd32, yd33, yd34, yd35, yd36, yd37, yd38, yd39, yd40, yd41, yd42, yd43, yd44, yd45, yd46, yd47,yd48,yd49,yd50,yd51,yd52,yd53,yd54,yd55;



 realtype f[138];
  realtype b[138];

  realtype r[138];
  realtype PCH3CH2CH3_0, PCH3CHCH2_0 , PH2_0 , PCH3CCH_0, PCH3CH3_0 , PCH2CH2_0, PCHCH_0, PCH4_0,Ptot;
  realtype T, PCH3CH2CH3, PCH3CHCH2, PH2, PCH3CCH, PCH3CH3, PCH2CH2, PCHCH, PCH4, kB;
  realtype K[138];



  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3); y4 = Ith(y,4); y5 = Ith(y,5); y6 = Ith(y,6); y7 = Ith(y,7); y8 = Ith(y,8); y9 = Ith(y,9); y10 = Ith(y,10);
  y11 =Ith(y,11); y12=Ith(y,12); y13=Ith(y,13); y14=Ith(y,14); y15=Ith(y,15); y16=Ith(y,16); y17=Ith(y,17); y18=Ith(y,18); y19=Ith(y,19); y20=Ith(y,20);
  y21 =Ith(y,21); y22=Ith(y,22); y23=Ith(y,23); y24=Ith(y,24); y25=Ith(y,25); y26=Ith(y,26); y27=Ith(y,27); y28=Ith(y,28); y29=Ith(y,29); y30=Ith(y,30);
  y31 =Ith(y,31); y32=Ith(y,32); y33=Ith(y,33); y34=Ith(y,34); y35=Ith(y,35); y36=Ith(y,36); y37=Ith(y,37); y38=Ith(y,38); y39=Ith(y,39); y40=Ith(y,40);
  y41 =Ith(y,41); y42=Ith(y,42); y43=Ith(y,43); y44=Ith(y,44); y45=Ith(y,45); y46=Ith(y,46); y47 = Ith(y,47); //Pt100 adsorbates and site

for (int i = 0; i < 138; i++) 
{f[i] = ff[i];
 }




// As these all have lateral interactions, theorectically there only f1 should be read into the function
  T = T_beta; kB = 8.617385e-5;
  alph = alph; beta = beta; gamm = gamm;
  PCH3CH2CH3_0 = PCH3CH2CH3_beta;  PCH3CHCH2_0 = PCH3CHCH2_beta;  PH2_0 = PH2_beta;  PCH3CCH_0 = PCH3CCH_beta; PCH3CH3_0 = PCH3CH3_beta; PCH2CH2_0 = PCH2CH2_beta; PCHCH_0 = PCHCH_beta; PCH4_0 = PCH4_beta;
  Ptot = PCH3CH2CH3_0 + PCH3CHCH2_0 + PH2_0 + PCH3CCH_0 + PCH3CH3_0 + PCH2CH2_0 + PCHCH_0 + PCH4_0;


//H-Lateral_interactions
double H_H_slope                       	= 0.000;
double CH3CH2CH3_H_slope               	= 0.000;
double CH3CHCH3_H_slope                	= 0.000;
double CH3CH2CH2_H_slope               	= 0.000;
double CH3CHCH2_H_slope                	= 0.000;
double CH3CH2CH_H_slope                	= 0.000;
double CH2CH2CH2_H_slope               	= 0.000;
double CH3CCH3_H_slope                 	= 0.000;
	
double CH3CH2C_H_slope                 	= 0.000;
double CH2CH2CH_H_slope                	= 0.000;
double CH2CHCH2_H_slope                	= 0.000;
double CH3CHCH_H_slope                 	= 0.000;
double CH3CCH2_H_slope                 	= 0.000;
double CH3CHC_H_slope                  	= 0.000;
double CH2CH2C_H_slope                 	= 0.000;
double CH2CHCH_H_slope                 	= 0.000;
double CHCH2CH_H_slope                 	= 0.000;
double CH2CCH2_H_slope                 	= 0.000;
double CH3CCH_H_slope                  	= 0.000;
	
double CH3CC_H_slope                   	= 0.000;
double CH2CHC_H_slope                  	= 0.000;
double CHCH2C_H_slope                  	= 0.000;
double CHCHCH_H_slope                  	= 0.000;
double CH2CCH_H_slope                  	= 0.000;
double CH2CC_H_slope                   	= 0.000;
double CHCHC_H_slope                   	= 0.000;
double CCH2C_H_slope                   	= 0.000;
double CHCCH_H_slope                   	= 0.000;
	
double CCHC_H_slope                    	= 0.000;
double CHCC_H_slope                    	= 0.000;
double CCC_H_slope                     	= 0.000;
	
double CH3CH3_H_slope                  	= 0.000;
double CH3CH2_H_slope                  	= 0.000;
double CH3CH_H_slope                   	= 0.000;
double CH3C_H_slope                    	= 0.000;
double CH2CH2_H_slope                  	= 0.000;
double CH2CH_H_slope                   	= 0.000;
double CH2C_H_slope                    	= 0.000;
double CHCH_H_slope                    	= 0.000;
double CHC_H_slope                     	= 0.000;
double CC_H_slope                      	= 0.000;
double CH4_H_slope                     	= 0.000;
double CH3_H_slope                     	= 0.000;
double CH2_H_slope                     	= 0.000;
double CH_H_slope                      	= 0.000;
double C_H_slope                       	= 0.000;

double CH3CH2CH3_CHCH_slope   			=      0    ;
double CH3CHCH3_CHCH_slope             =       0    ;
double CH3CH2CH2_CHCH_slope            =       0    ;
double CH3CHCH2_CHCH_slope             =       0    ;
double CH3CH2CH_CHCH_slope             =       0    ;
double CH2CH2CH2_CHCH_slope            =       0    ;
double CH3CCH3_CHCH_slope              =       0    ;
double CH3CH2C_CHCH_slope              =       0    ;
double CH2CH2CH_CHCH_slope             =       0    ;
double CH2CHCH2_CHCH_slope             =       0    ;
double CH3CHCH_CHCH_slope              =       0    ;
double CH3CCH2_CHCH_slope              =       0    ;
double CH3CHC_CHCH_slope               =       0    ;
double CH2CH2C_CHCH_slope              =       0    ;
double CH2CHCH_CHCH_slope              =       0    ;
double CHCH2CH_CHCH_slope              =       0    ;
double CH2CCH2_CHCH_slope              =       0    ;
double CH3CCH_CHCH_slope               =       0    ;
double CH3CC_CHCH_slope                =       0    ;
double CH2CHC_CHCH_slope               =       0    ;
double CHCH2C_CHCH_slope               =       0    ;
double CHCHCH_CHCH_slope               =       0    ;
double CH2CCH_CHCH_slope               =       0    ;
double CH2CC_CHCH_slope                =       0    ;
double CHCHC_CHCH_slope                =       0    ;
double CCH2C_CHCH_slope                =       0    ;
double CHCCH_CHCH_slope                =       0    ;
double CCHC_CHCH_slope                 =       0    ;
double CHCC_CHCH_slope                 =       0    ;
double CCC_CHCH_slope                  =       0    ;
double CH3CH3_CHCH_slope               =       0    ;
double CH3CH2_CHCH_slope               =       0    ;
double CH3CH_CHCH_slope                =       0    ;
double CH3C_CHCH_slope                 =       0    ;
double CH2CH2_CHCH_slope               =       0    ;
double CH2CH_CHCH_slope                =       0    ;
double CH2C_CHCH_slope                 =       0    ;
double CHCH_CHCH_slope                 =       0	;
double CHC_CHCH_slope                  =       0    ;
double CC_CHCH_slope                   =       0    ;
double CH4_CHCH_slope                  =       0    ;
double CH3_CHCH_slope                  =       0    ;
double CH2_CHCH_slope                  =       0    ;
double CH_CHCH_slope                   =       0    ;
double C_CHCH_slope                    =       0    ;
double H_CHCH_slope                    =       0    ;

double  CH3CH2CH3_C_slope              =       0  ;
double  CH3CHCH3_C_slope               =       0  ;
double  CH3CH2CH2_C_slope              =       0  ;
double  CH3CHCH2_C_slope               =       0  ;
double  CH3CH2CH_C_slope               =       0  ;
double  CH2CH2CH2_C_slope              =       0  ;
double  CH3CCH3_C_slope                =       0  ;
double  CH3CH2C_C_slope                =       0  ;
double  CH2CH2CH_C_slope               =       0  ;
double  CH2CHCH2_C_slope               =       0  ;
double  CH3CHCH_C_slope                =       0  ;
double  CH3CCH2_C_slope                =       0  ;
double  CH3CHC_C_slope                 =       0  ;
double  CH2CH2C_C_slope                =       0  ;
double  CH2CHCH_C_slope                =       0  ;
double  CHCH2CH_C_slope                =       0  ;
double  CH2CCH2_C_slope                =       0  ;
double  CH3CCH_C_slope                 =       0  ;
double  CH3CC_C_slope                  =       0  ;
double  CH2CHC_C_slope                 =       0  ;
double  CHCH2C_C_slope                 =       0  ;
double  CHCHCH_C_slope                 =       0  ;
double  CH2CCH_C_slope                 =       0  ;
double  CH2CC_C_slope                  =       0  ;
double  CHCHC_C_slope                  =       0  ;
double  CCH2C_C_slope                  =       0  ;
double  CHCCH_C_slope                  =       0  ;
double  CCHC_C_slope                   =       0  ;
double  CHCC_C_slope                   =       0  ;
double  CCC_C_slope                    =       0  ;
double  CH3CH3_C_slope                 =       0  ;
double  CH3CH2_C_slope                 =       0  ;
double  CH3CH_C_slope                  =       0  ;
double  CH3C_C_slope                   =       0  ;
double  CH2CH2_C_slope                 =       0  ;
double  CH2CH_C_slope                  =       0  ;
double  CH2C_C_slope                   =       0  ;
double  CHCH_C_slope                   =       0  ;
double  CHC_C_slope                    =       0  ;
double  CC_C_slope                     =       0  ;
double  CH4_C_slope                    =       0  ;
double  CH3_C_slope                    =       0  ;
double  CH2_C_slope                    =       0  ;
double  CH_C_slope                     =       0  ;
double  C_C_slope                      =       0  ;
double  H_C_slope                      =       0  ;




/*double H =  H_H_slope*y6;
*/
double CH3CH2CH3 	=p[0] + CH3CH2CH3_energy 	+ CH3CH2CH3_H_slope*y2+CH3CH2CH3_CHCH_slope*y40*4 + CH3CH2CH3_C_slope*y47*4;
double CH3CHCH3 	=p[1]+CH3CHCH3_energy 		+ CH3CHCH3_H_slope*y2+CH3CHCH3_CHCH_slope*y40*4 + CH3CHCH3_C_slope*y47*4;
double CH3CH2CH2 	=p[2]+ CH3CH2CH2_energy 	+ CH3CH2CH2_H_slope*y2+CH3CH2CH2_CHCH_slope*y40*4 + CH3CH2CH2_C_slope*y47*4;
double CH3CHCH2	 	=p[3]+CH3CHCH2_energy		+ CH3CHCH2_H_slope*y2+CH3CHCH2_CHCH_slope*y40*4 + CH3CH2CH2_C_slope*y47*4;
double CH3CH2CH 	=p[4]+CH3CH2CH_energy 		+ CH3CH2CH_H_slope*y2+CH3CH2CH_CHCH_slope*y40*4 + CH3CH2CH_C_slope*y47*4;
double CH2CH2CH2	=p[5]+CH2CH2CH2_energy 		+ CH2CH2CH2_H_slope*y2+CH2CH2CH2_CHCH_slope*y40*4 + CH2CH2CH2_C_slope*y47*4;
double CH3CCH3		=p[6]+CH3CCH3_energy 		+ CH3CCH3_H_slope*y2+CH3CCH3_CHCH_slope*y40*4+ CH3CCH3_C_slope*y47*4;

double CH3CH2C		=p[7]+ CH3CH2C_energy		+ CH3CH2C_H_slope*y2+CH3CH2C_CHCH_slope*y40*4 + CH3CH2C_C_slope*y47*4;
double CH2CH2CH		=p[8]+ CH2CH2CH_energy + CH2CH2CH_H_slope*y2+CH2CH2CH_CHCH_slope*y40*4 + CH2CH2CH_C_slope*y47*4;
double CH2CHCH2		=p[9]+ CH2CHCH2_energy + CH2CHCH2_H_slope*y2+CH2CHCH2_CHCH_slope*y40*4 + CH2CHCH2_C_slope*y47*4;
double CH3CHCH 		=p[10]+ CH3CHCH_energy + CH3CHCH_H_slope*y2+CH3CHCH_CHCH_slope*y40*4 + CH3CHCH_C_slope*y47*4;
double CH3CCH2 		=p[11]+ CH3CCH2_energy +CH3CCH2_H_slope*y2+CH3CCH2_CHCH_slope*y40*4 + CH3CCH2_C_slope*y47*4;
double CH3CHC  		=p[12]+CH3CHC_energy + CH3CHC_H_slope*y2+CH3CHC_CHCH_slope*y40*4 + CH3CHC_C_slope*y47*4;
double CH2CH2C 		=p[13]+CH2CH2C_energy +CH2CH2C_H_slope*y2+CH2CH2C_CHCH_slope*y40*4 + CH2CH2C_C_slope*y47*4;
double CHCH2CH 		=p[14]+CHCH2CH_energy +CHCH2CH_H_slope*y2+CHCH2CH_CHCH_slope*y40*4 + CHCH2CH_C_slope*y47*4;
double CH2CHCH 		=p[15]+CH2CHCH_energy +  CH2CHCH_H_slope*y2+CH2CHCH_CHCH_slope*y40*4 + CH2CHCH_C_slope*y47*4;
double CH2CCH2 		=p[16]+CH2CCH2_energy + CH2CCH2_H_slope*y2+CH2CCH2_CHCH_slope*y40*4 + CH2CCH2_C_slope*y47*4;
double CH3CCH  		=p[17]+ CH3CCH_energy + CH3CCH_H_slope*y2+CH3CCH_CHCH_slope*y40*4 + CH3CCH_C_slope*y47*4;

double CH3CC  		=p[18]+ CH3CC_energy + CH3CC_H_slope*y2+CH3CC_CHCH_slope*y40*4 + CH3CC_C_slope*y47*4;
double CH2CHC 		=p[19]+ CH2CHC_energy + CH2CHC_H_slope*y2+CH2CHC_CHCH_slope*y40*4 + CH2CHC_C_slope*y47*4;
double CHCHCH 		=p[20]+ CHCHCH_energy + CHCHCH_H_slope*y2+CHCHCH_CHCH_slope*y40*4 + CHCHCH_C_slope*y47*4;
double CHCH2C 		=p[21]+CHCH2C_energy + CHCH2C_H_slope*y2+CHCH2C_CHCH_slope*y40*4 + CHCH2C_C_slope*y47*4;
double CH2CCH 		=p[22]+CH2CCH_energy +CH2CCH_H_slope*y2+CH2CCH_CHCH_slope*y40*4 + CH2CCH_C_slope*y47*4;
double CH2CC  		=p[23]+CH2CC_energy + CH2CC_H_slope*y2+CH2CC_CHCH_slope*y40*4 + CH2CC_C_slope*y47*4;
double CHCHC  		=p[24]+CHCHC_energy + CHCHC_H_slope*y2+CHCHC_CHCH_slope*y40*4 + CHCHC_C_slope*y47*4;
double CCH2C  		=p[25]+CCH2C_energy + CCH2C_H_slope*y2+CCH2C_CHCH_slope*y40*4 + CCH2C_C_slope*y47*4;
double CHCCH  		=p[26]+CHCCH_energy + CHCCH_H_slope*y2+CHCCH_CHCH_slope*y40*4 + CHCCH_C_slope*y47*4;
double CCHC   		=p[27]+CCHC_energy + CCHC_H_slope*y2+CCHC_CHCH_slope*y40*4 + CCHC_C_slope*y47*4;
double CHCC   		=p[28]+CHCC_energy + CHCC_H_slope*y2+CHCC_CHCH_slope*y40*4 + CHCC_C_slope*y47*4;
double CCC    		=p[29]+CCC_energy + CCC_H_slope*y2+CCC_CHCH_slope*y40*4 + CCC_C_slope*y47*4;

double CH3CH3		=p[30]+ CH3CH3_energy + CH3CH3_H_slope*y2+CH3CH3_CHCH_slope*y40*4 + CH3CH3_C_slope*y47*4;
double CH3CH2		=p[31]+ CH3CH2_energy + CH3CH2_H_slope*y2+CH3CH2_CHCH_slope*y40*4 + CH3CH2_C_slope*y47*4;
double CH3CH 		=p[32]+ CH3CH_energy + CH3CH_H_slope*y2+CH3CH_CHCH_slope*y40*4 + CH3CH_C_slope*y47*4;
double CH3C  		=p[33]+ CH3C_energy + CH3C_H_slope*y2+CH3C_CHCH_slope*y40*4 + CH3C_C_slope*y47;
double CH2CH2		=p[34]+ CH2CH2_energy + CH2CH2_H_slope*y2+CH2CH2_CHCH_slope*y40*4 + CH2CH2_C_slope*y47*4;
double CH2CH 		=p[35]+CH2CH_energy + CH2CH_H_slope*y2+CH2CH_H_slope*y40*4 + CH2CH_C_slope*y47*4;
double CH2C  		=p[36]+CH2C_energy + CH2C_H_slope*y2+CH2C_CHCH_slope*y40*4 + CH2C_C_slope*y47*4;
double CHCH  		=p[37]+CHCH_energy + CHCH_H_slope*y2+CHCH_CHCH_slope*y40*4 + CHCH_C_slope*y47*4;
double CHC   		=p[38]+CHC_energy + CHC_H_slope*y2+CHC_CHCH_slope*y40*4 + CHC_C_slope*y47*4;
double CC    		=p[39]+CC_energy + CC_H_slope*y2+CC_CHCH_slope*y40*4 + CC_C_slope*y47*4;

double CH4   		=p[40]+CH4_energy + CH4_H_slope*y2+CH4_CHCH_slope*y40*4 + CH4_C_slope*y47*4;
double CH3   		=p[41]+CH3_energy + CH3_H_slope*y2+CH3_CHCH_slope*y40*4 + CH3_C_slope*y47*4;
double CH2   		=p[42]+CH2_energy + CH2_H_slope*y2+CH2_CHCH_slope*y40*4 + CH2_C_slope*y47*4;
double CH    		=p[43]+CH_energy + CH_H_slope*y2+CH_CHCH_slope*y40*4 + CH_C_slope*y47*4;
double C     		=p[44]+C_energy + C_H_slope*y2+C_CHCH_slope*y40*4 + C_C_slope*y47*4;
double H     		=p[45]+H_energy + H_H_slope*y2+H_CHCH_slope*y40*4 + H_C_slope*y47*4;
//y2 s H coverage. y40 is CHCH coverage
double CH3CH2CH3_0 = CH3CH2CH3_energy + p[0]; //I am not sure about this, as do we want to propagate uncertainty into the reaction energy slopes? If we don;t, I need to exchange XxXxXx_energy with XxXxXx_0
double CH3CHCH3_0 =p[1]+CH3CHCH3_energy ;
double CH3CH2CH2_0 =p[2]+ CH3CH2CH2_energy ;
double CH3CHCH2_0 =p[3]+CH3CHCH2_energy ;
double CH3CH2CH_0 =p[4]+CH3CH2CH_energy ;
double CH2CH2CH2_0 =p[5]+CH2CH2CH2_energy ;
double CH3CCH3_0 =p[6]+CH3CCH3_energy ;

//C5 and C4
double CH3CH2C_0=p[7]+ CH3CH2C_energy;
double CH2CH2CH_0=p[8]+ CH2CH2CH_energy;
double CH2CHCH2_0=p[9]+ CH2CHCH2_energy;
double CH3CHCH_0 =p[10]+ CH3CHCH_energy;
double CH3CCH2_0 =p[11]+ CH3CCH2_energy;
double CH3CHC_0  =p[12]+CH3CHC_energy ;
double CH2CH2C_0 =p[13]+CH2CH2C_energy ;
double CHCH2CH_0 =p[14]+CHCH2CH_energy ;
double CH2CHCH_0 =p[15]+CH2CHCH_energy ;
double CH2CCH2_0 =p[16]+CH2CCH2_energy ;
double CH3CCH_0  =p[17]+ CH3CCH_energy ;
////C3H3to C3
double CH3CC_0  =p[18]+ CH3CC_energy ;
double CH2CHC_0 =p[19]+ CH2CHC_energy ;
double CHCHCH_0 =p[20]+ CHCHCH_energy ;
double CHCH2C_0 =p[21]+CHCH2C_energy ;
double CH2CCH_0 =p[22]+CH2CCH_energy ;
double CH2CC_0  =p[23]+CH2CC_energy ;
double CHCHC_0  =p[24]+CHCHC_energy ;
double CCH2C_0  =p[25]+CCH2C_energy ;
double CHCCH_0  =p[26]+CHCCH_energy ;
double CCHC_0   =p[27]+CCHC_energy ;
double CHCC_0   =p[28]+CHCC_energy ;
double CCC_0    =p[29]+CCC_energy ;
////CH3CH3 -> C
double CH3CH3_0 =p[30]+ CH3CH3_energy ;
double CH3CH2_0 =p[31]+ CH3CH2_energy ;
double CH3CH_0  =p[32]+ CH3CH_energy ;
double CH3C_0   =p[33]+ CH3C_energy ;
double CH2CH2_0 =p[34]+ CH2CH2_energy ;
double CH2CH_0  =p[35]+CH2CH_energy ;
double CH2C_0  =p[36]+CH2C_energy ;
double CHCH_0  =p[37]+CHCH_energy ;
double CHC_0   =p[38]+CHC_energy ;
double CC_0    =p[39]+CC_energy ;

  double CH4_0 =p[40]+CH4_energy ;
  double CH3_0 =p[41]+CH3_energy ;
  double CH2_0 =p[42]+CH2_energy ;
  double CH_0  =p[43]+CH_energy ;
  double C_0   =p[44]+C_energy ;
  double H_0   =p[45]+H_energy ;




//I know this looks a little strange, but these next 130 reactions basically are 0.5(drxnG(theta1,theta2) - dG(0,0), which is the correction term for activation energy barriers
double CH3CH2CH3_to_CH3CHCH3_H_H_slope = 0.5*((H+CH3CHCH3 - CH3CH2CH3)-(H_0     +CH3CHCH3_0      - CH3CH2CH3_0     ));
double CH3CH2CH3_to_CH3CH2CH2_H_H_slope= 0.5*((H+CH3CH2CH2 - CH3CH2CH3)-(H_0     +CH3CH2CH2_0      - CH3CH2CH3_0     ));
double CH3CHCH3_to_CH3CHCH2_H_H_slope  = 0.5*((H+CH3CHCH2 - CH3CHCH3)-(H_0     +CH3CHCH2_0      - CH3CHCH3_0     ));
double CH3CH2CH2_to_CH3CHCH2_H_H_slope = 0.5*((H+CH3CHCH2 - CH3CH2CH2)-(H_0     +CH3CHCH2_0      - CH3CH2CH2_0     ));
double CH3CH2CH3_to_CH3_CH3CH2_H_slope = 0.5*((CH3+CH3CH2 - CH3CH2CH3)-(CH3_0      + CH3CH2_0      - CH3CH2CH3_0     ));
double CH3CHCH3_to_CH3_CH3CH_H_slope   = 0.5*((CH3+CH3CH - CH3CHCH3) -(CH3_0      + CH3CH_0      - CH3CHCH3_0     ));
double CH3CHCH3_to_CH3CCH3_H_H_slope     = 0.5*((H + CH3CCH3 - CH3CHCH3) - (H_0      + CH3CCH3_0      - CH3CHCH3_0     ));
double CH3CH2CH2_to_CH3CH2_CH2_H_slope   = 0.5*((CH3CH2 + CH2 - CH3CH2CH2) -(CH3CH2_0      + CH2_0      - CH3CH2CH2_0     ));
double CH3CH2CH2_to_CH3_CH2CH2_H_slope   = 0.5*((CH3+CH2CH2 - CH3CH2CH2) - (CH3_0      + CH2CH2_0      - CH3CH2CH2_0     ));
double CH3CH2CH2_to_CH2CH2CH2_H_H_slope  = 0.5*((H + CH2CH2CH2 - CH3CH2CH2) - (H_0     +CH2CH2CH2_0      - CH3CH2CH2_0     ));
double CH3CH2CH2_to_CH3CH2CH_H_H_slope   = 0.5*((H + CH3CH2CH - CH3CH2CH2) - (H_0      + CH3CH2CH_0      - CH3CH2CH2_0     ));

double CH3CH2CH_to_CH3CH2_CH_H_slope   = 0.5*((CH3CH2 + CH - CH3CH2CH) - (CH3CH2_0      + CH_0      - CH3CH2CH_0     ));
double CH3CH2CH_to_CH3_CH2CH_H_slope   = 0.5*((CH2CH + CH3 - CH3CH2CH) - (CH2CH_0      + CH3_0      - CH3CH2CH_0     ));
double CH3CH2CH_to_CH3CH2C_H_H_slope   = 0.5*((H+CH3CH2C - CH3CH2CH) - (H_0      + CH3CH2C_0      - CH3CH2CH_0     ));
double CH3CH2CH_to_CH3CHCH_H_H_slope   = 0.5*((H+CH3CHCH - CH3CH2CH) - (H_0      + CH3CHCH_0      - CH3CH2CH_0     ));
double CH3CH2CH_to_CH2CH2CH_H_H_slope  = 0.5*((H+CH2CH2CH - CH3CH2CH) - (H_0      + CH2CH2CH_0      - CH3CH2CH_0     ));
double CH2CH2CH2_to_CH2_CH2CH2_H_slope = 0.5*((CH2+CH2CH2 - CH2CH2CH2) - (CH2_0      + CH2CH2_0      - CH2CH2CH2_0     ));
double CH2CH2CH2_to_CH2CH2CH_H_H_slope = 0.5*((CH2CH2CH + H - CH2CH2CH2) -(CH2CH2CH_0      + H_0      - CH2CH2CH2_0     ));
double CH2CH2CH2_to_CH2CHCH2_H_H_slope = 0.5*((H + CH2CHCH2 - CH2CH2CH2) - (CH2CHCH2_0      +H_0      - CH2CH2CH2_0     ));
double CH3CHCH2_to_CH3CH_CH2_H_slope = 0.5*((CH3CH+CH2 - CH3CHCH2) - (CH3CH_0      + CH2_0      - CH3CHCH2_0     ));
double CH3CHCH2_to_CH3_CHCH2_H_slope = 0.5*((CH3 + CH2CH - CH3CHCH2) - (CH3_0      + CH2CH_0      - CH3CHCH2_0     ));
double CH3CHCH2_to_CH3CCH2_H_H_slope = 0.5*((CH3CCH2 + H - CH3CHCH2) - (CH3CCH2_0      + H_0      - CH3CHCH2_0     ));
double CH3CHCH2_to_CH3CHCH_H_H_slope = 0.5*((CH3CHCH + H - CH3CHCH2) - (CH3CHCH_0      + H_0      - CH3CHCH2_0     ));
double CH3CHCH2_to_CH2CHCH2_H_H_slope = 0.5*((CH2CHCH2 + H - CH3CHCH2) - (CH2CHCH2_0      + H_0      - CH3CHCH2_0     ));
double CH3CCH3_to_CH3_CCH3_H_slope = 0.5*((CH3 + CH3C - CH3CCH3)-(CH3_0      +CH3C_0      - CH3CCH3_0     ));
double CH3CCH3_to_CH3CCH2_H_H_slope =0.5*((H + CH3CCH2 - CH3CCH3)-(H_0      + CH3CCH2_0      - CH3CCH3_0     ));

double CH3CH2C_to_CH3_CH2C_H_slope = 0.5*((CH3 + CH2C - CH3CH2C) -(CH3_0      + CH2C_0      - CH3CH2C_0     ));
double CH3CH2C_to_CH3CH2_C_H_slope = 0.5*((CH3CH2+C - CH3CH2C) - (CH3CH2_0      + C_0      - CH3CH2C_0     ));
double CH3CH2C_to_CH2CH2C_H_H_slope  = 0.5*((CH2CH2C +H - CH3CH2C) - (CH2CH2C_0      + H_0      - CH3CH2C_0     ));
double CH3CH2C_to_CH3CHC_H_H_slope = 0.5*((CH3CHC +H - CH3CH2C) -(CH3CHC_0      +H_0      - CH3CH2C_0     ));
double CH2CH2CH_to_CH2_CH2CH_H_slope = 0.5*((CH2 + CH2CH - CH2CH2CH) -(CH2_0      + CH2CH_0      - CH2CH2CH_0     ));
double CH2CH2CH_to_CH2CH2_CH_H_slope = 0.5*((CH2CH2 + CH - CH2CH2CH) -(CH2CH2_0      + CH_0      - CH2CH2CH_0     ));
double CH2CH2CH_to_CH2CH2C_H_H_slope = 0.5*((CH2CH2C +H - CH2CH2CH) - (CH2CH2C_0      + H_0      - CH2CH2CH_0     ));
double CH2CH2CH_to_CH2CHCH_H_H_slope = 0.5*((CH2CHCH +H - CH2CH2CH) - (CH2CHCH_0      + H_0      - CH2CH2CH_0     ));
double CH2CH2CH_to_CHCH2CH_H_H_slope = 0.5*((CHCH2CH +H - CH2CH2CH) - (CHCH2CH_0      + H_0      - CH2CH2CH_0     ));
double CH2CHCH2_to_CH2_CHCH2_H_slope = 0.5*((CH2 +CH2CH - CH2CHCH2) - (CH2_0      + CH2CH_0      - CH2CHCH2_0     ));
double CH2CHCH2_to_CH2CHCH_H_H_slope = 0.5*((CH2CHCH +H - CH2CHCH2) - (CH2CHCH_0      + CH2CH_0      - CH2CHCH2_0     ));
double CH2CHCH2_to_CH2CCH2_H_H_slope = 0.5*((CH2CCH2 +H - CH2CHCH2) - (CH2CCH2_0      + H_0      - CH2CHCH2_0     ));
double CH3CHCH_to_CH3_CHCH_H_slope = 0.5*((CH3 + CHCH - CH3CHCH) - (CH3_0      + CHCH_0      - CH3CHCH_0     ));
double CH3CHCH_to_CH3CH_CH_H_slope = 0.5*((CH3CH + CH - CH3CHCH) - (CH3CH_0      + CH_0      - CH3CHCH_0     ));
double CH3CHCH_to_CH3CHC_H_H_slope = 0.5*((CH3CHC +H - CH3CHCH) - (CH3CHC_0      + H_0      - CH3CHCH_0     ));
double CH3CHCH_to_CH3CCH_H_H_slope = 0.5*((CH3CCH +H - CH3CHCH) - (CH3CCH_0      + H_0      - CH3CHCH_0     ));
double CH3CHCH_to_CH2CHCH_H_H_slope = 0.5*((CH2CHCH +H - CH3CHCH) - (CH2CHCH_0      +H_0      - CH3CHCH_0     ));
double CH3CCH2_to_CH3_CCH2_H_slope = 0.5*((CH3 + CH2C - CH3CCH2) - (CH3_0      + CH2C_0      - CH3CCH2_0     ));
double CH3CCH2_to_CH3C_CH2_H_slope = 0.5*((CH3C + CH2 - CH3CCH2) - (CH3C_0      + CH2_0      - CH3CCH2_0     ));
double CH3CCH2_to_CH2CCH2_H_H_slope = 0.5*((CH2CCH2 + H - CH3CCH2) - (CH2CCH2_0      + H_0      - CH3CCH2_0     ));
double CH3CCH2_to_CH3CCH_H_H_slope = 0.5*((CH3CCH +H - CH3CCH2) - (CH3CCH_0      + H_0      - CH3CCH2_0      ));

double CH3CHC_to_CH3_CHC_H_slope = 0.5*((CH3 + CHC - CH3CHC) -(CH3_0      + CHC_0      - CH3CHC_0     ));
double CH3CHC_to_CH3CH_C_H_slope = 0.5*((CH3CH + C - CH3CHC) -(CH3CH_0      + C_0      - CH3CHC_0     ));
double CH3CHC_to_CH3CC_H_H_slope = 0.5*((CH3CC +H - CH3CHC) - (CH3CC_0      + H_0      - CH3CHC_0     ));
double CH3CHC_to_CH2CHC_H_H_slope = 0.5*((CH2CHC +H - CH3CHC) - (CH2CHC_0      +H_0      - CH3CHC_0     ));
double CH2CH2C_to_CH2CH2_C_H_slope = 0.5*((CH2CH2 + C - CH2CH2C) -(CH2CH2_0      +C_0      - CH2CH2C_0     ));
double CH2CH2C_to_CH2_CH2C_H_slope = 0.5*((CH2 + CH2C - CH2CH2C) -(CH2_0      + CH2C_0      - CH2CH2C_0     ));
double CH2CH2C_to_CH2CHC_H_H_slope = 0.5*((CH2CHC +H - CH2CH2C) - (CH2CHC_0      + H_0      - CH2CH2C_0     ));
double CH2CH2C_to_CHCH2C_H_H_slope = 0.5*((CHCH2C +H - CH2CH2C) - (CHCH2C_0      + H_0      - CH2CH2C_0     ));
double CHCH2CH_to_CHCH2_CH_H_slope = 0.5*((CH2CH + CH - CHCH2CH) -(CH2CH_0      + CH_0      - CHCH2CH_0     ));
double CHCH2CH_to_CHCH2C_H_H_slope = 0.5*((CHCH2C +H - CHCH2CH) - (CHCH2C_0      + H_0      - CHCH2CH_0     ));
double CHCH2CH_to_CHCHCH_H_H_slope = 0.5*((CHCHCH + H - CHCH2CH) - (CHCHCH_0      +H_0      - CHCH2CH_0     ));
double CH2CHCH_to_CH2_CHCH_H_slope = 0.5*((CH2 +CHCH - CH2CHCH) - (CH2_0      + CHCHC_0      - CH2CHCH_0     ));
double CH2CHCH_to_CH2CH_CH_H_slope = 0.5*((CH2CH + CH - CH2CHCH) - (CH2CH_0      + CH_0      - CH2CHCH_0     ));
double CH2CHCH_to_CH2CHC_H_H_slope = 0.5*((CH2CHC +H - CH2CHCH)-(CH2CHC_0      + H_0      - CH2CHCH_0     ));
double CH2CHCH_to_CH2CCH_H_H_slope = 0.5*((CH2CCH +H - CH2CHCH) -(CH2CCH_0      +H_0      - CH2CHCH_0     ));
double CH2CHCH_to_CHCHCH_H_H_slope = 0.5*((CHCHCH +H - CH2CHCH) -(CHCHCH_0      + H_0      - CH2CHCH_0     ));
double CH2CCH2_to_CH2C_CH2_H_slope = 0.5*((CH2C + CH2 - CH2CCH2) - (CH2C_0      + CH2_0      - CH2CCH2_0     ));
double CH2CCH2_to_CH2CCH_H_H_slope = 0.5*((CH2CCH + H - CH2CCH2) - (CH2CCH_0      + H_0      - CH2CCH2_0     ));
double CH3CCH_to_CH3C_CH_H_slope = 0.5*((CH3C + CH - CH3CCH) - (CH3C_0      + CH_0      - CH3CCH_0     ));
double CH3CCH_to_CH3_CHC_H_slope = 0.5*((CH3 + CHC - CH3CCH) - (CH3_0      + CHC_0      - CH3CCH_0     ));
double CH3CCH_to_CH3CC_H_H_slope = 0.5*((CH3CC +H - CH3CCH) - ( CH3CC_0      + H_0      - CH3CCH_0     ));
double CH3CCH_to_CH2CCH_H_H_slope =0.5*((CH2CCH +H - CH3CCH) -(CH2CCH_0      + H_0      - CH2CCH_0     ));

double CH3CC_to_CH3_CC_H_slope = 0.5*((CH3+CC - CH3CC) - (CH3_0      +CC_0      -   CH3CC_0     ));
double CH3CC_to_CH3C_C_H_slope = 0.5*((CH3C +C - CH3CC) - (CH3C_0      +C_0      - CH3CC_0     ));
double CH3CC_to_CH2CC_H_H_slope = 0.5*((CH2CC +H - CH3CC) - (CH2CC_0      +H_0      - CH3CC_0     ));
double CH2CHC_to_CH2_CHC_H_slope = 0.5*((CH2 + CHC - CH2CHC) -(CH2_0      +CHC_0      - CH2CHC_0     ));
double CH2CHC_to_CH2CH_C_H_slope = 0.5*((CH2CH + C - CH2CHC)-(CH2CH_0      + C_0      - CH2CHC_0     ));
double CH2CHC_to_CH2CC_H_H_slope = 0.5*((CH2CC +H - CH2CHC) -(CH2CC_0      + H_0      - CH2CHC_0     ));
double CH2CHC_to_CHCHC_H_H_slope = 0.5*((CHCHC +H - CH2CHC) - (CHCHC_0      + H_0      - CH2CHC_0     ));
double CHCH2C_to_CH_CH2C_H_slope = 0.5*((CH2C + CH - CHCH2C) - (CH2C_0      + CH_0      - CHCH2C_0     ));
double CHCH2C_to_CH2CH_C_H_slope = 0.5*((CH2CH + C - CHCH2C) - (CH2CH_0      + C_0      - CHCH2C_0     ));
double CHCH2C_to_CHCHC_H_H_slope = 0.5*((CHCHC + H - CHCH2C) - (CHCHC_0      + H_0      - CHCH2C_0     ));
double CHCH2C_to_CCH2C_H_H_slope = 0.5*((CCH2C +H - CHCH2C) - (CCH2C_0      + H_0      - CHCH2C_0     ));
double CHCHCH_to_CH_CHCH_H_slope = 0.5*((CH + CHCH - CHCHCH) - (CH_0      +CHCH_0      - CHCHCH_0     ));
double CHCHCH_to_CHCHC_H_H_slope = 0.5*((CHCHC +H - CHCHCH) - (CHCHC_0      +H_0      -CHCHCH_0     ));
double CHCHCH_to_CHCCH_H_H_slope = 0.5*((CHCCH +H - CHCHCH) - (CHCCH_0      +H_0      -CHCHCH_0     ));
double CH2CCH_to_CH2_CCH_H_slope = 0.5*((CH2+CHC - CH2CCH) -(CH2_0      + CHC_0      - CH2CCH_0     ));
double CH2CCH_to_CH2C_CH_H_slope = 0.5*((CH2C +CH -CH2CCH) -(CH2C_0      + CH_0      - CH2CCH_0     ));
double CH2CCH_to_CH2CC_H_H_slope = 0.5*((CH2CC + H - CH2CCH) - (CH2CC_0      + H_0      - CH2CCH_0     ));
double CH2CCH_to_CHCCH_H_H_slope = 0.5*((CHCCH + H - CH2CCH) - (CHCCH_0      + H_0      - CH2CCH_0     ));

double CH2CC_to_CH2_CC_H_slope = 0.5*((CH2 + CC - CH2CC) -(CH2_0      +CC_0      -CH2CC_0     ));
double CH2CC_to_CH2C_C_H_slope = 0.5*((CH2C + C - CH2CC) -(CH2C_0      +C_0      - CH2CC_0     ));
double CH2CC_to_CHCC_H_H_slope = 0.5*((CHCC + H -CH2CC) - (CHCC_0      +H_0      - CH2CC_0     ));
double CHCHC_to_CH_CHC_H_slope = 0.5*((CH+CHC - CHCHC) - (CH_0      + CHC_0      - CHCHC_0     ));
double CHCHC_to_CHCH_C_H_slope = 0.5*((CHCH + C - CHCHC) - (CHCH_0      +C_0      - CHCHC_0     ));
double CHCHC_to_CHCC_H_H_slope = 0.5*((CHCC + H - CHCHC) - (CHCC_0      +H_0      - CHCHC_0     ));
double CHCHC_to_CCHC_H_H_slope = 0.5*((CCHC + H - CHCHC) - (CCHC_0      +H_0      - CHCHC_0     ));
double CCH2C_to_C_CH2C_H_slope = 0.5*((CH2C + C - CCH2C) - (CH2C_0      + C_0      - CCH2C_0     ));
double CCH2C_to_CCHC_H_H_slope = 0.5*((CCHC + H - CCH2C) - (CCHC_0      + H_0      - CCH2C_0     ));
double CHCCH_to_CH_CHC_H_slope = 0.5*((CH + CHC - CHCCH) - (CH_0      + CHC_0      - CHCCH_0     ));
double CHCCH_to_CHCC_H_H_slope = 0.5*((CHCC + H - CHCCH) - (CHCC_0      + H_0      - CHCCH_0     ));
double CCHC_to_CHC_C_H_slope = 0.5*((CHC +C - CCHC) - (CHC_0      + C_0      - CCHC_0     ));
double CCHC_to_CCC_H_H_slope = 0.5*((CCC +H - CCHC) - (CCC_0      +H_0      - CCHC_0     ));
double CHCC_to_CH_CC_H_slope = 0.5*((CH +CC - CHCC) - (CH_0      + CC_0      - CHCC_0     ));
double CHCC_to_CHC_C_H_slope = 0.5*((CHC +C - CHCC) - (CHC_0      + C_0      - CHCC_0     ));
double CHCC_to_CCC_H_H_slope = 0.5*((CCC +H - CHCC) - (CCC_0      + H_0      - CHCC_0     ));
double CCC_to_CC_C_H_slope = 0.5*((CC+C - CCC) - (CC_0      + C_0      - CCC_0     ));

double CH3CH3_to_CH3_CH3_H_slope = 0.5*((CH3 +CH3 - CH3CH3) - (CH3_0      + CH3_0      -CH3CH3_0     ));
double CH3CH3_to_CH3CH2_H_H_slope = 0.5*((CH3CH2 + H - CH3CH3) - (CH3CH2_0      +H_0      - CH3CH3_0     ));
double CH3CH2_to_CH3_CH2_H_slope= 0.5*((CH3+CH2 - CH3CH2) -( CH3_0      + CH2_0      - CH3CH2_0     ));
double CH3CH2_to_CH3CH_H_H_slope = 0.5*((CH3CH + H - CH3CH2) - (CH3CH_0      +H_0      - CH3CH2_0     ));
double CH3CH2_to_CH2CH2_H_H_slope = 0.5*((CH2CH2 +H - CH3CH2) - (CH2CH2_0      + H_0      - CH3CH2_0     ));
double CH3CH_to_CH3_CH_H_slope = 0.5*((CH3 + CH - CH3CH) - (CH3_0      + CH_0      - CH3CH_0     ));
double CH3CH_to_CH3C_H_H_slope = 0.5*((CH3C + H - CH3CH) - (CH3C_0      + H_0      - CH3CH_0     ));
double CH3CH_to_CH2CH_H_H_slope = 0.5*((CH2CH +H - CH3CH) - (CH2CH_0      +H_0      - CH3CH_0     ));
double CH3C_to_CH3_C_H_slope = 0.5*((CH3+C-CH3C) - (CH3_0      + C_0      - CH3C_0     ));
double CH3C_to_CH2C_H_H_slope = 0.5 *(( CH2C + H - CH3C) - (CH2C_0    + H_0      - CH3C_0     ));
double CH2CH2_to_CH2_CH2_H_slope =0.5 *((CH2 + CH2 - CH2CH2) - (2*CH2_0      - CH2CH2_0     ));
double CH2CH2_to_CH2CH_H_H_slope =0.5*((CH2CH +H - CH2CH2) - (CH2CH_0      + H_0      - CH2CH2_0     ));
double CH2CH_to_CH2_CH_H_slope = 0.5*((CH2 + CH - CH2CH) - (CH2_0      + CH_0      -CH2CH_0     ));
double CH2CH_to_CH2C_H_H_slope = 0.5*((CH2C + H - CH2CH) - (CH2C_0      + H_0      - CH2CH_0     ));
double CH2CH_to_CHCH_H_H_slope = 0.5*((CHCH + H - CH2CH) - (CHCH_0      + H_0      - CH2CH_0     ));
double CH2C_to_CH2_C_H_slope = 0.5*((CH2 + C - CH2C) - (CH2_0      +C_0      - CH2C_0     ));
double CH2C_to_CHC_H_H_slope = 0.5*((CHC + H - CH2C) - (CHC_0      +H_0      - CH2C_0     ));
double CHCH_to_CH_CH_H_slope = 0.5*((CH+CH-CHCH) - (2*CH_0      - CHCH_0     ));
double CHCH_to_CHC_H_H_slope = 0.5*((CHC+H - CHCH) - (CHC_0      +H_0      - CHCH_0     ));
double CHC_to_CH_C_H_slope = 0.5*((CH +C -CHC) - (CH_0      +C_0      - CHC_0     ));
double CHC_to_CC_H_H_slope = 0.5*((CC +H - CHC) - (CC_0      + H_0      - CHC_0     ));
double CC_to_C_C_H_slope = 0.5*((2*C - CC) - (2*C_0      - CC_0     ));

double CH4_to_CH3_H_H_slope = 0.5*((CH3+H-CH4) - (CH3_0      +H_0      - CH4_0     ));
double CH3_to_CH2_H_H_slope = 0.5*((CH2+H-CH3) - (CH2_0      +H_0      - CH3_0     ));
double CH2_to_CH_H_H_slope = 0.5*((CH+H-CH2) - (CH_0      +H_0      - CH2_0     ));
double CH_to_C_H_H_slope = 0.5*((C+H-CH) - (C_0      +H_0      - CH_0     ));

//The following is an example of Eric's stuff. The 0.77 is the Activation barrier for the reaction.
/*double H2O_to_OH_H = 0.773369890999959 + H2O_to_OH_H_CO_slope*y1 + H2O_to_OH_H_H_slope*y4 - H2O;

if (H2O_to_OH_H < 0)
    H2O_to_OH_H = 0.0;*/

//The below looks like the following
//TST_energy = uncertainty + TST_energy (raw, not subtracted from reactant) + lateral interaction term - Reactant energy (with uncertainty?)



double CH3CH2CH3_to_CH3CHCH3_H = p[46]+ dEA9 + CH3CH2CH3_to_CH3CHCH3_H_H_slope - CH3CH2CH3_0;
if (CH3CH2CH3_to_CH3CHCH3_H < (CH3CHCH3 + H - CH3CH2CH3))
   CH3CH2CH3_to_CH3CHCH3_H = CH3CHCH3 + H - CH3CH2CH3;
if (CH3CH2CH3_to_CH3CHCH3_H <0)
   CH3CH2CH3_to_CH3CHCH3_H = 0.001;

double CH3CH2CH3_to_CH3CH2CH2_H = p[47]+ dEA10 + CH3CH2CH3_to_CH3CH2CH2_H_H_slope - CH3CH2CH3_0;
if ( CH3CH2CH3_to_CH3CH2CH2_H <( CH3CH2CH2 + H - CH3CH2CH3))
   CH3CH2CH3_to_CH3CH2CH2_H = CH3CH2CH2 + H - CH3CH2CH3;
if (CH3CH2CH3_to_CH3CH2CH2_H <0)
        CH3CH2CH3_to_CH3CH2CH2_H = 0.001;

double CH3CHCH3_to_CH3CHCH2_H = p[48]+ dEA11 + CH3CHCH3_to_CH3CHCH2_H_H_slope - CH3CHCH3_0;
if (CH3CHCH3_to_CH3CHCH2_H < (CH3CHCH2 + H - CH3CHCH3))
    CH3CHCH3_to_CH3CHCH2_H = CH3CHCH2 + H - CH3CHCH3;

if (CH3CHCH3_to_CH3CHCH2_H <0)
        CH3CHCH3_to_CH3CHCH2_H = 0.001;

double CH3CH2CH2_to_CH3CHCH2_H = p[49]+ dEA12 + CH3CH2CH2_to_CH3CHCH2_H_H_slope - CH3CH2CH2_0;
if (CH3CH2CH2_to_CH3CHCH2_H < (CH3CHCH2 + H - CH3CH2CH2))
    CH3CH2CH2_to_CH3CHCH2_H = CH3CHCH2 + H - CH3CH2CH2;

if (CH3CH2CH2_to_CH3CHCH2_H <0)
        CH3CH2CH2_to_CH3CHCH2_H = 0.001;

double CH3CH2CH3_to_CH3CH2_CH3 = p[50] +dEA13 + CH3CH2CH3_to_CH3_CH3CH2_H_slope - CH3CH2CH3_0;
if (CH3CH2CH3_to_CH3CH2_CH3 < (CH3 + CH3CH2 - CH3CH2CH3))
    CH3CH2CH3_to_CH3CH2_CH3  = CH3 + CH3CH2 - CH3CH2CH3;

if (CH3CH2CH3_to_CH3CH2_CH3 <0)
        CH3CH2CH3_to_CH3CH2_CH3 =0.001;

double CH3CHCH3_to_CH3_CH3CH = p[51] + dEA14 + CH3CHCH3_to_CH3_CH3CH_H_slope - CH3CHCH3_0;
if (CH3CHCH3_to_CH3_CH3CH <( CH3CH + CH3 - CH3CHCH3))
   CH3CHCH3_to_CH3_CH3CH = CH3CH + CH3 - CH3CHCH3;

if (CH3CHCH3_to_CH3_CH3CH <0 )
        CH3CHCH3_to_CH3_CH3CH =0.001;

double CH3CHCH3_to_CH3CCH3_H = p[52] + dEA15 + CH3CHCH3_to_CH3CCH3_H_H_slope - CH3CHCH3_0;
if (CH3CHCH3_to_CH3CCH3_H < (CH3CCH3 + H - CH3CHCH3 ))
	CH3CHCH3_to_CH3CCH3_H = CH3CCH3 + H - CH3CHCH3;

if (CH3CHCH3_to_CH3CCH3_H <0 )
        CH3CHCH3_to_CH3CCH3_H =0.001;


double CH3CH2CH2_to_CH3CH2_CH2 = p[53]+ dEA16 + CH3CH2CH2_to_CH3CH2_CH2_H_slope - CH3CH2CH2_0;
if (CH3CH2CH2_to_CH3CH2_CH2 < (CH3CH2 + CH2 - CH3CH2CH2))
	CH3CH2CH2_to_CH3CH2_CH2 = CH3CH2 + CH2 - CH3CH2CH2;

if (CH3CH2CH2_to_CH3CH2_CH2 <0)
        CH3CH2CH2_to_CH3CH2_CH2 = 0.001;

double CH3CH2CH2_to_CH3_CH2CH2 = p[54] + dEA17 + CH3CH2CH2_to_CH3_CH2CH2_H_slope- CH3CH2CH2_0;
if (CH3CH2CH2_to_CH3_CH2CH2 < (CH3 + CH2CH2 - CH3CH2CH2))
	CH3CH2CH2_to_CH3_CH2CH2 = CH3 + CH2CH2 - CH3CH2CH2;

if (CH3CH2CH2_to_CH3_CH2CH2 <0)
        CH3CH2CH2_to_CH3_CH2CH2 =0.001;

double CH3CH2CH2_to_CH2CH2CH2_H = p[55] + dEA18 +CH3CH2CH2_to_CH2CH2CH2_H_H_slope - CH3CH2CH2_0;
if (CH3CH2CH2_to_CH2CH2CH2_H < (CH2CH2CH2 + H - CH3CH2CH2))
	CH3CH2CH2_to_CH2CH2CH2_H = CH2CH2CH2 + H - CH3CH2CH2;

if (CH3CH2CH2_to_CH2CH2CH2_H <0)
        CH3CH2CH2_to_CH2CH2CH2_H =0.001;

double CH3CH2CH2_to_CH3CH2CH_H = p[56] + dEA19 + CH3CH2CH2_to_CH3CH2CH_H_H_slope - CH3CH2CH2_0;
if (CH3CH2CH2_to_CH3CH2CH_H < (CH3CH2CH + H - CH3CH2CH2))
	CH3CH2CH2_to_CH3CH2CH_H = CH3CH2CH + H - CH3CH2CH2;

if (CH3CH2CH2_to_CH3CH2CH_H <0)
        CH3CH2CH2_to_CH3CH2CH_H =0.001;

double CH3CH2CH_to_CH3CH2_CH = p[57] +dEA20 + CH3CH2CH_to_CH3CH2_CH_H_slope - CH3CH2CH_0;
if (CH3CH2CH_to_CH3CH2_CH < (CH3CH2 + CH - CH3CH2CH))
    CH3CH2CH_to_CH3CH2_CH = CH3CH2 + CH - CH3CH2CH;

if (CH3CH2CH_to_CH3CH2_CH <0)
        CH3CH2CH_to_CH3CH2_CH =0.001;

double CH3CH2CH_to_CH3_CH2CH = p[58] + dEA21 + CH3CH2CH_to_CH3_CH2CH_H_slope - CH3CH2CH_0;
if (CH3CH2CH_to_CH3_CH2CH < (CH2CH + CH3 - CH3CH2CH))
	CH3CH2CH_to_CH3_CH2CH = CH2CH + CH3 - CH3CH2CH;

if (CH3CH2CH_to_CH3_CH2CH <0)
        CH3CH2CH_to_CH3_CH2CH =0.001;

double CH3CH2CH_to_CH3CH2C_H = p[59] + dEA22 + CH3CH2CH_to_CH3CH2C_H_H_slope - CH3CH2CH_0;
if (CH3CH2CH_to_CH3CH2C_H < (CH3CH2C + H - CH3CH2CH))
	CH3CH2CH_to_CH3CH2C_H = CH3CH2C + H - CH3CH2CH;

if (CH3CH2CH_to_CH3CH2C_H <0)
        CH3CH2CH_to_CH3CH2C_H =0.001;

double CH3CH2CH_to_CH3CHCH_H = p[60] + dEA23 + CH3CH2CH_to_CH3CHCH_H_H_slope - CH3CH2CH_0;
if (CH3CH2CH_to_CH3CHCH_H < (CH3CHCH + H - CH3CH2CH))
    CH3CH2CH_to_CH3CHCH_H = (CH3CHCH + H - CH3CH2CH);

if (CH3CH2CH_to_CH3CHCH_H <0)
        CH3CH2CH_to_CH3CHCH_H =0.001;

double CH3CH2CH_to_CH2CH2CH_H = p[61] + dEA24 + CH3CH2CH_to_CH2CH2CH_H_H_slope - CH3CH2CH_0;
if (CH3CH2CH_to_CH2CH2CH_H < (CH2CH2CH + H - CH3CH2CH))
	CH3CH2CH_to_CH2CH2CH_H = (CH2CH2CH + H - CH3CH2CH);

if (CH3CH2CH_to_CH2CH2CH_H <0)
        CH3CH2CH_to_CH2CH2CH_H =0.001;

double CH2CH2CH2_to_CH2CH2_CH2 = p[62] + dEA25 + CH2CH2CH2_to_CH2_CH2CH2_H_slope - CH2CH2CH2_0;
if (CH2CH2CH2_to_CH2CH2_CH2 < (CH2CH2 + CH2 - CH2CH2CH2))
    CH2CH2CH2_to_CH2CH2_CH2 = (CH2CH2 + CH2 - CH2CH2CH2);

if (CH2CH2CH2_to_CH2CH2_CH2 <0)
        CH2CH2CH2_to_CH2CH2_CH2 =0.001;

double CH2CH2CH2_to_CH2CH2CH_H = p[63] + dEA26 + CH2CH2CH2_to_CH2CH2CH_H_H_slope - CH2CH2CH2_0;
if (CH2CH2CH2_to_CH2CH2CH_H < (CH2CH2CH + H - CH2CH2CH2))
	CH2CH2CH2_to_CH2CH2CH_H < (CH2CH2CH + H - CH2CH2CH2);

if (CH2CH2CH2_to_CH2CH2CH_H <0)
        CH2CH2CH2_to_CH2CH2CH_H =0.001;

double CH2CH2CH2_to_CH2CHCH2_H = p[64] + dEA27 + CH2CH2CH2_to_CH2CHCH2_H_H_slope - CH2CH2CH2_0;
if (CH2CH2CH2_to_CH2CHCH2_H < (CH2CHCH2 + H - CH2CH2CH2))
	CH2CH2CH2_to_CH2CHCH2_H = (CH2CHCH2 + H - CH2CH2CH2);

if (CH2CH2CH2_to_CH2CHCH2_H <0)
        CH2CH2CH2_to_CH2CHCH2_H =0.001;

double CH3CHCH2_to_CH3_CHCH2 = p[65] + dEA28 + CH3CHCH2_to_CH3_CHCH2_H_slope - CH3CHCH2_0;
if (CH3CHCH2_to_CH3_CHCH2 < (CH3 + CH2CH - CH3CHCH2))
	CH3CHCH2_to_CH3_CHCH2 = (CH3 + CH2CH - CH3CHCH2);

if (CH3CHCH2_to_CH3_CHCH2 <0)
        CH3CHCH2_to_CH3_CHCH2 =0.001;
	
double CH3CHCH2_to_CH3CH_CH2 = p[66] + dEA29 + CH3CHCH2_to_CH3CH_CH2_H_slope - CH3CHCH2_0;
if (CH3CHCH2_to_CH3CH_CH2 < (CH3CH + CH2 - CH3CHCH2))
	CH3CHCH2_to_CH3CH_CH2 = (CH3CH + CH2 - CH3CHCH2);

if (CH3CHCH2_to_CH3CH_CH2<0)
        CH3CHCH2_to_CH3CH_CH2 =0.001;

double CH3CHCH2_to_CH3CCH2_H = p[67] + dEA30 + CH3CHCH2_to_CH3CCH2_H_H_slope - CH3CHCH2_0;
if (CH3CHCH2_to_CH3CCH2_H < (CH3CCH2 + H - CH3CHCH2))
	CH3CHCH2_to_CH3CCH2_H = (CH3CCH2 + H - CH3CHCH2);

if (CH3CHCH2_to_CH3CCH2_H <0)
        CH3CHCH2_to_CH3CCH2_H =0.001;

double CH3CHCH2_to_CH3CHCH_H = p[68] + dEA31 + CH3CHCH2_to_CH3CHCH_H_H_slope - CH3CHCH2_0;
if (CH3CHCH2_to_CH3CHCH_H < (CH3CHCH + H - CH3CHCH2))
	CH3CHCH2_to_CH3CHCH_H = (CH3CHCH + H - CH3CHCH2);

if (CH3CHCH2_to_CH3CHCH_H <0)
        CH3CHCH2_to_CH3CHCH_H =0.001;

double CH3CHCH2_to_CH2CHCH2_H = p[69] +dEA32 + CH3CHCH2_to_CH2CHCH2_H_H_slope - CH3CHCH2_0;
if (CH3CHCH2_to_CH2CHCH2_H < (CH2CHCH2 + H - CH3CHCH2))
	CH3CHCH2_to_CH2CHCH2_H = (CH2CHCH2 + H - CH3CHCH2);

if (CH3CHCH2_to_CH2CHCH2_H <0)
        CH3CHCH2_to_CH2CHCH2_H =0.001;

double CH3CCH3_to_CH3_CCH3 = p[70] + dEA33 + CH3CCH3_to_CH3_CCH3_H_slope - CH3CCH3_0;
if (CH3CCH3_to_CH3_CCH3 < (CH3C + CH3 - CH3CCH3))
	CH3CCH3_to_CH3_CCH3 = (CH3C + CH3 - CH3CCH3);

if (CH3CCH3_to_CH3_CCH3 <0)
        CH3CCH3_to_CH3_CCH3 =0.001;

double CH3CCH3_to_CH3CCH2_H = p[71] + dEA34 + CH3CCH3_to_CH3CCH2_H_H_slope - CH3CCH3_0;
if (CH3CCH3_to_CH3CCH2_H < (CH3CCH2 + H - CH3CCH3));
	CH3CCH3_to_CH3CCH2_H = (CH3CCH2 + H - CH3CCH3);

if (CH3CCH3_to_CH3CCH2_H <0 )
        CH3CCH3_to_CH3CCH2_H = 0.001;

double CH3CH2C_to_CH3_CH2C = p[72] + dEA35 + CH3CH2C_to_CH3_CH2C_H_slope - CH3CH2C_0;
if (CH3CH2C_to_CH3_CH2C < (CH3 + CH2C - CH3CH2C))
	CH3CH2C_to_CH3_CH2C = (CH3 + CH2C - CH3CH2C);

if (CH3CH2C_to_CH3_CH2C <0)
        CH3CH2C_to_CH3_CH2C = 0.001;

double CH3CH2C_to_CH3CH2_C = p[73] + dEA36 + CH3CH2C_to_CH3CH2_C_H_slope - CH3CH2C_0;
if (CH3CH2C_to_CH3CH2_C < (CH3CH2 + C - CH3CH2C))
	CH3CH2C_to_CH3CH2_C = (CH3CH2 + C - CH3CH2C);

if (CH3CH2C_to_CH3CH2_C <0)
        CH3CH2C_to_CH3CH2_C =0.001;

double CH3CH2C_to_CH2CH2C_H = p[74] + dEA37 + CH3CH2C_to_CH2CH2C_H_H_slope - CH3CH2C_0;
if (CH3CH2C_to_CH2CH2C_H < (CH2CH2C + H - CH3CH2C))
	CH3CH2C_to_CH2CH2C_H = (CH2CH2C + H - CH3CH2C);

if (CH3CH2C_to_CH2CH2C_H <0)
        CH3CH2C_to_CH2CH2C_H =0.001;
	
double CH3CH2C_to_CH3CHC_H = p[75] + dEA38 + CH3CH2C_to_CH3CHC_H_H_slope - CH3CH2C_0;
if (CH3CH2C_to_CH3CHC_H < (CH3CHC + H - CH3CH2C))
	CH3CH2C_to_CH3CHC_H = (CH3CHC + H - CH3CH2C);

if (CH3CH2C_to_CH3CHC_H <0 )
        CH3CH2C_to_CH3CHC_H =0.001;

double CH2CH2CH_to_CH2_CH2CH = p[76] + dEA39 + CH2CH2CH_to_CH2_CH2CH_H_slope - CH2CH2CH_0;
if (CH2CH2CH_to_CH2_CH2CH < (CH2CH + CH2 - CH2CH2CH))
	CH2CH2CH_to_CH2_CH2CH = (CH2CH + CH2 - CH2CH2CH);

if (CH2CH2CH_to_CH2_CH2CH <0)
        CH2CH2CH_to_CH2_CH2CH =0.001;

double CH2CH2CH_to_CH2CH2_CH = p[77] + dEA40 + CH2CH2CH_to_CH2CH2_CH_H_slope - CH2CH2CH_0;
if (CH2CH2CH_to_CH2CH2_CH < (CH2CH2 + CH - CH2CH2CH))
	CH2CH2CH_to_CH2CH2_CH = (CH2CH2 + CH - CH2CH2CH);

if (CH2CH2CH_to_CH2CH2_CH <0)
        CH2CH2CH_to_CH2CH2_CH =0.001;

double CH2CH2CH_to_CH2CH2C_H = p[78] + dEA41 + CH2CH2CH_to_CH2CH2C_H_H_slope  - CH2CH2CH_0;
if (CH2CH2CH_to_CH2CH2C_H < (CH2CH2C + H - CH2CH2CH))
	CH2CH2CH_to_CH2CH2C_H = (CH2CH2C + H - CH2CH2CH);

if (CH2CH2CH_to_CH2CH2C_H <0)
        CH2CH2CH_to_CH2CH2C_H =0.001;

double CH2CH2CH_to_CH2CHCH_H = p[79] + dEA42 + CH2CH2CH_to_CH2CHCH_H_H_slope  - CH2CH2CH_0;
if (CH2CH2CH_to_CH2CHCH_H < (CH2CHCH +H - CH2CH2CH))
	CH2CH2CH_to_CH2CHCH_H = (CH2CHCH +H - CH2CH2CH);

if (CH2CH2CH_to_CH2CHCH_H <0)
        CH2CH2CH_to_CH2CHCH_H =0.001;

double CH2CH2CH_to_CHCH2CH_H = p[80] + dEA43 + CH2CH2CH_to_CHCH2CH_H_H_slope  - CH2CH2CH_0;
if (CH2CH2CH_to_CHCH2CH_H < (CHCH2CH + H - CH2CH2CH))
	CH2CH2CH_to_CHCH2CH_H = (CHCH2CH + H - CH2CH2CH);

if (CH2CH2CH_to_CHCH2CH_H <0)
        CH2CH2CH_to_CHCH2CH_H =0.001;

double CH2CHCH2_to_CH2_CH2CH = p[81] + dEA44 + CH2CHCH2_to_CH2_CHCH2_H_slope - CH2CHCH2_0;
if (CH2CHCH2_to_CH2_CH2CH < (CH2 + CH2CH - CH2CHCH2))
	CH2CHCH2_to_CH2_CH2CH = (CH2 + CH2CH - CH2CHCH2);

if (CH2CHCH2_to_CH2_CH2CH <0)
        CH2CHCH2_to_CH2_CH2CH = 0.001;

double CH2CHCH2_to_CH2CHCH_H = p[82] + dEA45 + CH2CHCH2_to_CH2CHCH_H_H_slope - CH2CHCH2_0;
if (CH2CHCH2_to_CH2CHCH_H < (CH2CHCH + H - CH2CHCH2))
	CH2CHCH2_to_CH2CHCH_H = (CH2CHCH + H - CH2CHCH2);

if (CH2CHCH2_to_CH2CHCH_H <0)
        CH2CHCH2_to_CH2CHCH_H = 0.001;

double CH2CHCH2_to_CH2CCH2_H = p[83] + dEA46 + CH2CHCH2_to_CH2CCH2_H_H_slope - CH2CHCH2_0;
if (CH2CHCH2_to_CH2CCH2_H < (CH2CCH2 + H - CH2CHCH2));
	CH2CHCH2_to_CH2CCH2_H = (CH2CCH2 + H - CH2CHCH2);

if (CH2CHCH2_to_CH2CCH2_H <0)
        CH2CHCH2_to_CH2CCH2_H =0.001;

double CH3CHCH_to_CH3_CHCH = p[84] + dEA47 + CH3CHCH_to_CH3_CHCH_H_slope - CH3CHCH_0;
if (CH3CHCH_to_CH3_CHCH < ( CH3 + CHCH - CH3CHCH));
	CH3CHCH_to_CH3_CHCH = ( CH3 + CHCH - CH3CHCH);

if (CH3CHCH_to_CH3_CHCH <0)
        CH3CHCH_to_CH3_CHCH = 0.001;

double CH3CHCH_to_CH3CH_CH = p[85] + dEA48 + CH3CHCH_to_CH3CH_CH_H_slope - CH3CHCH_0;
if (CH3CHCH_to_CH3CH_CH < (CH3CH + CH - CH3CHCH))
	CH3CHCH_to_CH3CH_CH = (CH3CH + CH - CH3CHCH);

if (CH3CHCH_to_CH3CH_CH <0)
        CH3CHCH_to_CH3CH_CH =0.001;

double CH3CHCH_to_CH3CHC_H = p[86] + dEA49 + CH3CHCH_to_CH3CHC_H_H_slope - CH3CHCH_0;
if (CH3CHCH_to_CH3CHC_H < (CH3CHC + H - CH3CHCH))
	CH3CHCH_to_CH3CHC_H = (CH3CHC + H - CH3CHCH) ;

if (CH3CHCH_to_CH3CHC_H <0)
        CH3CHCH_to_CH3CHC_H =0.001;

double CH3CHCH_to_CH3CCH_H = p[87] + dEA50 + CH3CHCH_to_CH3CCH_H_H_slope - CH3CHCH_0;
if (CH3CHCH_to_CH3CCH_H < ( CH3CCH + H - CH3CHCH))
	CH3CHCH_to_CH3CCH_H = ( CH3CCH + H - CH3CHCH);

if (CH3CHCH_to_CH3CCH_H <0)
        CH3CHCH_to_CH3CCH_H =0.001;

double CH3CHCH_to_CH2CHCH_H = p[88] + dEA51 + CH3CHCH_to_CH2CHCH_H_H_slope - CH3CHCH_0;
if (CH3CHCH_to_CH2CHCH_H < (CH2CHCH + H - CH3CHCH))
	CH3CHCH_to_CH2CHCH_H < (CH2CHCH + H - CH3CHCH);

if (CH3CHCH_to_CH2CHCH_H <0)
        CH3CHCH_to_CH2CHCH_H=0.001;

double CH3CCH2_to_CH3_CCH2 = p[89] + dEA52 + CH3CCH2_to_CH3_CCH2_H_slope - CH3CCH2_0;
if (CH3CCH2_to_CH3_CCH2 < (CH2C + CH3 - CH3CCH2))
	CH3CCH2_to_CH3_CCH2 = (CH2C + CH3 - CH3CCH2);

if (CH3CCH2_to_CH3_CCH2 <0)
        CH3CCH2_to_CH3_CCH2 = 0.001;

double CH3CCH2_to_CH3C_CH2 = p[90] + dEA53 + CH3CCH2_to_CH3C_CH2_H_slope - CH3CCH2_0;
if (CH3CCH2_to_CH3C_CH2 < (CH3C + CH2 - CH3CCH2))
	CH3CCH2_to_CH3C_CH2 = (CH3C + CH2 - CH3CCH2);

if (CH3CCH2_to_CH3C_CH2 <0)
        CH3CCH2_to_CH3C_CH2 =0.001;

double CH3CCH2_to_CH2CCH2_H = p[91] + dEA54 + CH3CCH2_to_CH2CCH2_H_H_slope - CH3CCH2_0;
if (CH3CCH2_to_CH2CCH2_H < (CH2CCH2 + H - CH3CCH2))
	CH3CCH2_to_CH2CCH2_H = (CH2CCH2 + H - CH3CCH2);

if (CH3CCH2_to_CH2CCH2_H <0)
        CH3CCH2_to_CH2CCH2_H =0.001;
	
double CH3CCH2_to_CH3CCH_H = p[92] + dEA55 + CH3CCH2_to_CH3CCH_H_H_slope - CH3CCH2_0;
if (CH3CCH2_to_CH3CCH_H < (CH3CCH + H - CH3CCH2))
	CH3CCH2_to_CH3CCH_H = (CH3CCH + H - CH3CCH2);

if (CH3CCH2_to_CH3CCH_H <0)
        CH3CCH2_to_CH3CCH_H = 0.001;

double CH3CHC_to_CH3_CHC = p[93] + dEA56 + CH3CHC_to_CH3_CHC_H_slope - CH3CHC_0;
if (CH3CHC_to_CH3_CHC < (CH3 + CHC - CH3CHC ))
	CH3CHC_to_CH3_CHC = (CH3 + CHC - CH3CHC );

if (CH3CHC_to_CH3_CHC <0)
        CH3CHC_to_CH3_CHC =0.001;

double CH3CHC_to_CH3CH_C = p[94] + dEA57 + CH3CHC_to_CH3CH_C_H_slope - CH3CHC_0;
if (CH3CHC_to_CH3CH_C < (CH3CH + C - CH3CHC))
	CH3CHC_to_CH3CH_C = (CH3CH + C - CH3CHC);

if (CH3CHC_to_CH3CH_C <0)
        CH3CHC_to_CH3CH_C =0.001;

double CH3CHC_to_CH3CC_H = p[95] + dEA58 + CH3CHC_to_CH3CC_H_H_slope - CH3CHC_0;
if (CH3CHC_to_CH3CC_H < (CH3CC +H - CH3CHC))
	CH3CHC_to_CH3CC_H = (CH3CC +H - CH3CHC);

if (CH3CHC_to_CH3CC_H <0)
        CH3CHC_to_CH3CC_H =0.001;

double CH3CHC_to_CH2CHC_H = p[96] + dEA59 + CH3CHC_to_CH2CHC_H_H_slope - CH3CHC_0;
if (CH3CHC_to_CH2CHC_H < (CH2CHC + H - CH3CHC))
	CH3CHC_to_CH2CHC_H = (CH2CHC + H - CH3CHC);

if (CH3CHC_to_CH2CHC_H <0)
        CH3CHC_to_CH2CHC_H =0.001;

double CH2CH2C_to_CH2CH2_C = p[97] + dEA60 + CH2CH2C_to_CH2CH2_C_H_slope - CH2CH2C_0;
if (CH2CH2C_to_CH2CH2_C < (CH2CH2 + C - CH2CH2C))
	CH2CH2C_to_CH2CH2_C = (CH2CH2 + C - CH2CH2C);

if (CH2CH2C_to_CH2CH2_C <0)
        CH2CH2C_to_CH2CH2_C =0.001;

double CH2CH2C_to_CH2_CH2C = p[98] + dEA61 + CH2CH2C_to_CH2_CH2C_H_slope - CH2CH2C_0;
if (CH2CH2C_to_CH2_CH2C < (CH2 + CH2C - CH2CH2C))
	CH2CH2C_to_CH2_CH2C = (CH2 + CH2C - CH2CH2C);

if (CH2CH2C_to_CH2_CH2C <0)
        CH2CH2C_to_CH2_CH2C =0.001;

double CH2CH2C_to_CH2CHC_H = p[99] + dEA62 + CH2CH2C_to_CH2CHC_H_H_slope - CH2CH2C_0;
if (CH2CH2C_to_CH2CHC_H < (CH2CHC + H - CH2CH2C))
	CH2CH2C_to_CH2CHC_H = (CH2CHC + H - CH2CH2C);

if (CH2CH2C_to_CH2CHC_H <0)
        CH2CH2C_to_CH2CHC_H =0.001;

double CH2CH2C_to_CHCH2C_H = p[100] + dEA63 + CH2CH2C_to_CHCH2C_H_H_slope - CH2CH2C_0;
if (CH2CH2C_to_CHCH2C_H < (CHCH2C + H - CH2CH2C))
	CH2CH2C_to_CHCH2C_H = (CHCH2C + H - CH2CH2C);

if (CH2CH2C_to_CHCH2C_H <0)
        CH2CH2C_to_CHCH2C_H =0.001;

double CHCH2CH_to_CHCH2_CH = p[101] + dEA64 + CHCH2CH_to_CHCH2_CH_H_slope - CHCH2CH_0;
if (CHCH2CH_to_CHCH2_CH < (CH + CH2CH - CHCH2CH))
	CHCH2CH_to_CHCH2_CH = (CH + CH2CH - CHCH2CH);

if (CHCH2CH_to_CHCH2_CH <0)
        CHCH2CH_to_CHCH2_CH =0.001;

double CHCH2CH_to_CHCH2C_H = p[102] + dEA65 + CHCH2CH_to_CHCH2C_H_H_slope - CHCH2CH_0;
if (CHCH2CH_to_CHCH2C_H < (CHCH2C +H - CHCH2CH))
	CHCH2CH_to_CHCH2C_H = (CHCH2C +H - CHCH2CH);

if (CHCH2CH_to_CHCH2C_H <0)
        CHCH2CH_to_CHCH2C_H =0.001;

double CHCH2CH_to_CHCHCH_H = p[103] + dEA66 + CHCH2CH_to_CHCHCH_H_H_slope - CHCH2CH_0;
if (CHCH2CH_to_CHCHCH_H < (CHCHCH + H - CHCH2CH))
	CHCH2CH_to_CHCHCH_H = (CHCHCH + H - CHCH2CH);

if (CHCH2CH_to_CHCHCH_H <0)
        CHCH2CH_to_CHCHCH_H =0.001;

double CH2CHCH_to_CH2_CHCH = p[104] + dEA67 + CH2CHCH_to_CH2_CHCH_H_slope - CH2CHCH_0;
if (CH2CHCH_to_CH2_CHCH < (CH2+CHCH - CH2CHCH))
	CH2CHCH_to_CH2_CHCH = (CH2+CHCH - CH2CHCH) ;

if (CH2CHCH_to_CH2_CHCH <0)
        CH2CHCH_to_CH2_CHCH =0.001;

double CH2CHCH_to_CH2CH_CH = p[105] + dEA68 + CH2CHCH_to_CH2CH_CH_H_slope - CH2CHCH_0;
if (CH2CHCH_to_CH2CH_CH < (CH2CH + CH - CH2CHCH))
	CH2CHCH_to_CH2CH_CH = (CH2CH + CH - CH2CHCH);

if (CH2CHCH_to_CH2CH_CH <0)
        CH2CHCH_to_CH2CH_CH =0.001;

double CH2CHCH_to_CH2CHC_H = p[106] + dEA69 + CH2CHCH_to_CH2CHC_H_H_slope - CH2CHCH_0;
if (CH2CHCH_to_CH2CHC_H < (CH2CHC + H - CH2CHCH))
	CH2CHCH_to_CH2CHC_H = (CH2CHC + H - CH2CHCH);

if (CH2CHCH_to_CH2CHC_H <0)
        CH2CHCH_to_CH2CHC_H = 0.001;

double CH2CHCH_to_CH2CCH_H = p[107] + dEA70 + CH2CHCH_to_CH2CCH_H_H_slope - CH2CHCH_0;
if (CH2CHCH_to_CH2CCH_H < (CH2CCH + H - CH2CHCH))
	CH2CHCH_to_CH2CCH_H = (CH2CCH + H - CH2CHCH);

if (CH2CHCH_to_CH2CCH_H <0)
        CH2CHCH_to_CH2CCH_H =0.001;

double CH2CHCH_to_CHCHCH_H = p[108] + dEA71 + CH2CHCH_to_CHCHCH_H_H_slope - CH2CHCH_0;
if (CH2CHCH_to_CHCHCH_H < (CHCHCH +H - CH2CHCH))
	CH2CHCH_to_CHCHCH_H = (CHCHCH +H - CH2CHCH);

if (CH2CHCH_to_CHCHCH_H <0)
        CH2CHCH_to_CHCHCH_H =0.001;

double CH2CCH2_to_CH2C_CH2 = p[109] + dEA72 + CH2CCH2_to_CH2C_CH2_H_slope - CH2CCH2_0;
if (CH2CCH2_to_CH2C_CH2 < (CH2C + CH2 - CH2CCH2))
	CH2CCH2_to_CH2C_CH2 = (CH2C + CH2 - CH2CCH2);

if (CH2CCH2_to_CH2C_CH2 <0)
        CH2CCH2_to_CH2C_CH2 =0.001;

double CH2CCH2_to_CH2CCH_H = p[110] + dEA73 + CH2CCH2_to_CH2CCH_H_H_slope - CH2CCH2_0;
if (CH2CCH2_to_CH2CCH_H < (CH2CCH + H -CH2CCH2))
	CH2CCH2_to_CH2CCH_H = (CH2CCH + H -CH2CCH2);

if (CH2CCH2_to_CH2CCH_H <0)
        CH2CCH2_to_CH2CCH_H =0.001;

double CH3CCH_to_CH3C_CH = p[111] + dEA74 + CH3CCH_to_CH3C_CH_H_slope - CH3CCH_0;
if (CH3CCH_to_CH3C_CH < (CH3C + CH - CH3CCH))
	CH3CCH_to_CH3C_CH = (CH3C + CH - CH3CCH);

if (CH3CCH_to_CH3C_CH <0)
        CH3CCH_to_CH3C_CH =0.001;

double CH3CCH_to_CH3_CHC = p[112] + dEA75 + CH3CCH_to_CH3_CHC_H_slope - CH3CCH_0;
if (CH3CCH_to_CH3_CHC < (CHC + CH3 - CH3CCH))
	CH3CCH_to_CH3_CHC = (CHC + CH3 - CH3CCH);

if (CH3CCH_to_CH3_CHC<0)
        CH3CCH_to_CH3_CHC =0.001;

double CH3CCH_to_CH3CC_H = p[113] + dEA76 + CH3CCH_to_CH3CC_H_H_slope - CH3CCH_0;
if (CH3CCH_to_CH3CC_H < (CH3CC + H - CH3CCH))
	CH3CCH_to_CH3CC_H = (CH3CC + H - CH3CCH);

if (CH3CCH_to_CH3CC_H<0)
        CH3CCH_to_CH3CC_H =0.001;

double CH3CCH_to_CH2CCH_H = p[114] + dEA77 + CH3CCH_to_CH2CCH_H_H_slope - CH3CCH_0;
if (CH3CCH_to_CH2CCH_H < (CH2CCH + H - CH3CCH))
	CH3CCH_to_CH2CCH_H = (CH2CCH + H - CH3CCH);

if (CH3CCH_to_CH2CCH_H <0)
        CH3CCH_to_CH2CCH_H =0.001;

double CH3CC_to_CH3_CC = p[115] + dEA78 + CH3CC_to_CH3_CC_H_slope - CH3CC_0;
if (CH3CC_to_CH3_CC < (CH3 + CC - CH3CC))
	CH3CC_to_CH3_CC = (CH3 + CC - CH3CC);

if (CH3CC_to_CH3_CC<0)
        CH3CC_to_CH3_CC =0.001;

double CH3CC_to_CH3C_C = p[116] + dEA79 + CH3CC_to_CH3C_C_H_slope - CH3CC_0;
if (CH3CC_to_CH3C_C < (CH3C + C - CH3CC))
	CH3CC_to_CH3C_C = (CH3C + C - CH3CC);

if (CH3CC_to_CH3C_C <0)
        CH3CC_to_CH3C_C = 0.001;

double CH3CC_to_CH2CC_H = p[117] + dEA80 + CH3CC_to_CH2CC_H_H_slope - CH3CC_0;
if (CH3CC_to_CH2CC_H < (CH2CC + H - CH3CC))
	CH3CC_to_CH2CC_H = (CH2CC + H - CH3CC);

if (CH3CC_to_CH2CC_H <0)
        CH3CC_to_CH2CC_H =0.001;

double CH2CHC_to_CH2_CHC = p[118] + dEA81 + CH2CHC_to_CH2_CHC_H_slope - CH2CHC_0;
if (CH2CHC_to_CH2_CHC < (CHC + CH2 - CH2CHC))
	CH2CHC_to_CH2_CHC = (CHC + CH2 - CH2CHC);

if (CH2CHC_to_CH2_CHC <0)
        CH2CHC_to_CH2_CHC =0.001;

double CH2CHC_to_CH2CH_C = p[119] + dEA82 + CH2CHC_to_CH2CH_C_H_slope - CH2CHC_0;
if (CH2CHC_to_CH2CH_C < (C + CH2CH - CH2CHC))
	CH2CHC_to_CH2CH_C = (C + CH2CH - CH2CHC);

if (CH2CHC_to_CH2CH_C <0)
        CH2CHC_to_CH2CH_C =0.001;

double CH2CHC_to_CH2CC_H = p[120] + dEA83 + CH2CHC_to_CH2CC_H_H_slope - CH2CHC_0;
if (CH2CHC_to_CH2CC_H < (CH2CC + H - CH2CHC))
	CH2CHC_to_CH2CC_H = (CH2CC + H - CH2CHC);

if (CH2CHC_to_CH2CC_H <0)
        CH2CHC_to_CH2CC_H =0.001;

double CH2CHC_to_CHCHC_H = p[121] + dEA84 + CH2CHC_to_CHCHC_H_H_slope - CH2CHC_0;
if (CH2CHC_to_CHCHC_H < (CHCHC + H - CH2CHC))
	CH2CHC_to_CHCHC_H = (CHCHC + H - CH2CHC);

if (CH2CHC_to_CHCHC_H <0)
        CH2CHC_to_CHCHC_H =0.001;

double CHCH2C_to_CH_CH2C = p[122] + dEA85 + CHCH2C_to_CH_CH2C_H_slope - CHCH2C_0;
if (CHCH2C_to_CH_CH2C > (CH + CH2C - CHCH2C))
	CHCH2C_to_CH_CH2C = (CH + CH2C - CHCH2C);

if (CHCH2C_to_CH_CH2C <0)
        CHCH2C_to_CH_CH2C =0.001;

double CHCH2C_to_CH2CH_C = p[123] + dEA86 + CHCH2C_to_CH2CH_C_H_slope - CHCH2C_0;
if (CHCH2C_to_CH2CH_C < (CH2CH + C- CHCH2C))
	CHCH2C_to_CH2CH_C = (CH2CH + C- CHCH2C);

if (CHCH2C_to_CH2CH_C <0)
        CHCH2C_to_CH2CH_C =0.001;

double CHCH2C_to_CHCHC_H = p[124] + dEA87 + CHCH2C_to_CHCHC_H_H_slope - CHCH2C_0;
if (CHCH2C_to_CHCHC_H < (CHCHC + H -CHCH2C))
	CHCH2C_to_CHCHC_H = (CHCHC + H -CHCH2C);

if (CHCH2C_to_CHCHC_H <0)
        CHCH2C_to_CHCHC_H =0.001;

double CHCH2C_to_CCH2C_H = p[125] + dEA88 + CHCH2C_to_CCH2C_H_H_slope - CHCH2C_0;
if (CHCH2C_to_CCH2C_H < (CCH2C +H - CHCH2C))
	CHCH2C_to_CCH2C_H = (CCH2C +H - CHCH2C);

if (CHCH2C_to_CCH2C_H <0)
        CHCH2C_to_CCH2C_H = 0.001;

double CHCHCH_to_CHCH_CH = p[126] + dEA89 + CHCHCH_to_CH_CHCH_H_slope - CHCHCH_0;
if (CHCHCH_to_CHCH_CH < (CHCH+CH - CHCHCH))
	CHCHCH_to_CHCH_CH = CHCH + CH - CHCHCH;

if (CHCHCH_to_CHCH_CH <0)
        CHCHCH_to_CHCH_CH =0.001;

double CHCHCH_to_CHCHC_H = p[127] + dEA90 + CHCHCH_to_CHCHC_H_H_slope - CHCHCH_0;
if (CHCHCH_to_CHCHC_H < (CHCHC + H - CHCHCH))
	CHCHCH_to_CHCHC_H = (CHCHC + H - CHCHCH);

if (CHCHCH_to_CHCHC_H <0)
        CHCHCH_to_CHCHC_H =0.001;

double CHCHCH_to_CHCCH_H = p[128] + dEA91 + CHCHCH_to_CHCCH_H_H_slope - CHCHCH_0;
if (CHCHCH_to_CHCCH_H < (CHCCH + H - CHCHCH))
	CHCHCH_to_CHCCH_H = (CHCCH + H - CHCHCH);

if (CHCHCH_to_CHCCH_H <0)
        CHCHCH_to_CHCCH_H =0.001;

double CH2CCH_to_CH2_CCH = p[129] + dEA92 + CH2CCH_to_CH2_CCH_H_slope - CH2CCH_0;
if (CH2CCH_to_CH2_CCH < (CH2 + CHC - CH2CCH))
	CH2CCH_to_CH2_CCH = (CH2 + CHC - CH2CCH);

if (CH2CCH_to_CH2_CCH <0)
        CH2CCH_to_CH2_CCH =0.001;

double CH2CCH_to_CH2C_CH = p[130] + dEA93 + CH2CCH_to_CH2C_CH_H_slope - CH2CCH_0;
if (CH2CCH_to_CH2C_CH < (CH2C + CH - CH2CCH))
	CH2CCH_to_CH2C_CH = (CH2C + CH - CH2CCH);

if (CH2CCH_to_CH2C_CH <0)
        CH2CCH_to_CH2C_CH =0.001;

double CH2CCH_to_CH2CC_H = p[131] + dEA94 + CH2CCH_to_CH2CC_H_H_slope - CH2CCH_0;
if (CH2CCH_to_CH2CC_H < (CH2CC +H - CH2CCH))
	CH2CCH_to_CH2CC_H = (CH2CC +H - CH2CCH);

if (CH2CCH_to_CH2CC_H <0)
        CH2CCH_to_CH2CC_H =0.001;

double CH2CCH_to_CHCCH_H = p[132] + dEA95 + CH2CCH_to_CHCCH_H_H_slope - CH2CCH_0;
if (CH2CCH_to_CHCCH_H < (CHCCH + H - CH2CCH))
	CH2CCH_to_CHCCH_H = (CHCCH + H - CH2CCH);

if (CH2CCH_to_CHCCH_H <0)
        CH2CCH_to_CHCCH_H =0.001;

double CH2CC_to_CH2_CC = p[133] + dEA96 + CH2CC_to_CH2_CC_H_slope - CH2CC_0;
if (CH2CC_to_CH2_CC < (CH2 + CC - CH2CC))
	CH2CC_to_CH2_CC = (CH2 + CC - CH2CC);

if (CH2CC_to_CH2_CC <0)
        CH2CC_to_CH2_CC =0.001;

double CH2CC_to_CH2C_C = p[134] + dEA97 + CH2CC_to_CH2C_C_H_slope - CH2CC_0;
if (CH2CC_to_CH2C_C < (CH2C + C - CH2CC))
	CH2CC_to_CH2C_C = (CH2C + C - CH2CC);

if (CH2CC_to_CH2C_C <0)
        CH2CC_to_CH2C_C = 0.001;

double CH2CC_to_CHCC_H = p[135] + dEA98 + CH2CC_to_CHCC_H_H_slope - CH2CC_0;
if (CH2CC_to_CHCC_H < (CHCC + H - CH2CC))
	CH2CC_to_CHCC_H = (CHCC + H - CH2CC);

if (CH2CC_to_CHCC_H <0)
        CH2CC_to_CHCC_H =0.001;

double CHCHC_to_CH_CHC = p[136] + dEA99 + CHCHC_to_CH_CHC_H_slope - CHCHC_0;
if (CHCHC_to_CH_CHC < (CHC + CH - CHCHC))
	CHCHC_to_CH_CHC = (CHC + CH - CHCHC);

if (CHCHC_to_CH_CHC <0)
        CHCHC_to_CH_CHC =0.001;

double CHCHC_to_CHCH_C = p[137] + dEA100 + CHCHC_to_CHCH_C_H_slope - CHCHC_0;
if (CHCHC_to_CHCH_C < (CHCH + C - CHCHC))
	CHCHC_to_CHCH_C = (CHCH + C - CHCHC);

if (CHCHC_to_CHCH_C <0)
        CHCHC_to_CHCH_C = 0.001;

double CHCHC_to_CHCC_H = p[138] + dEA101 + CHCHC_to_CHCC_H_H_slope - CHCHC_0;
if (CHCHC_to_CHCC_H < (CHCC + H - CHCHC))
	CHCHC_to_CHCC_H = (CHCC + H - CHCHC);

if (CHCHC_to_CHCC_H <0)
        CHCHC_to_CHCC_H =0.001;

double CHCHC_to_CCHC_H = p[139] + dEA102 + CHCHC_to_CCHC_H_H_slope - CHCHC_0;
if (CHCHC_to_CCHC_H < (CCHC + H - CHCHC))
	CHCHC_to_CCHC_H = (CCHC + H - CHCHC);

if (CHCHC_to_CCHC_H <0)
        CHCHC_to_CCHC_H =0.001;

double CCH2C_to_CH2C_C = p[140] + dEA103 + CCH2C_to_C_CH2C_H_slope - CCH2C_0;
if (CCH2C_to_CH2C_C < (CH2C + C - CCH2C))
	CCH2C_to_CH2C_C = (CH2C + C - CCH2C);

if (CCH2C_to_CH2C_C <0)
        CCH2C_to_CH2C_C =0.001;

double CCH2C_to_CCHC_H = p[141] + dEA104 + CCH2C_to_CCHC_H_H_slope - CCH2C_0;
if (CCH2C_to_CCHC_H < (CCHC + H - CCH2C))
	CCH2C_to_CCHC_H = (CCHC + H - CCH2C);

if (CCH2C_to_CCHC_H <0)
        CCH2C_to_CCHC_H =0.001;

double CHCCH_to_CHC_CH = p[142] + dEA105 + CHCCH_to_CH_CHC_H_slope - CHCCH_0;
if (CHCCH_to_CHC_CH < (CHC + CH - CHCCH))
	CHCCH_to_CHC_CH = (CHC + CH - CHCCH);

if (CHCCH_to_CHC_CH <0)
        CHCCH_to_CHC_CH =0.001;

double CHCCH_to_CHCC_H = p[143] + dEA106 + CHCCH_to_CHCC_H_H_slope - CHCCH_0;
if (CHCCH_to_CHCC_H < (CHCC + H - CHCCH))
	CHCCH_to_CHCC_H = (CHCC + H - CHCCH);

if (CHCCH_to_CHCC_H <0)
        CHCCH_to_CHCC_H =0.001;

double CCHC_to_CHC_C = p[144] + dEA107 + CCHC_to_CHC_C_H_slope - CCHC_0;
if (CCHC_to_CHC_C < (C + CHC - CCHC))
	CCHC_to_CHC_C = (C + CHC - CCHC);

if (CCHC_to_CHC_C <0)
        CCHC_to_CHC_C =0.001;

double CCHC_to_CCC_H = p[145] + dEA108 + CCHC_to_CCC_H_H_slope - CCHC_0;
if (CCHC_to_CCC_H < (CCC +H - CCHC))
	CCHC_to_CCC_H = (CCC +H - CCHC);

if (CCHC_to_CCC_H <0)
        CCHC_to_CCC_H =0.001;

double CHCC_to_CH_CC = p[146] + dEA109 + CHCC_to_CH_CC_H_slope - CHCC_0;
if (CHCC_to_CH_CC< (CC+ CH - CHCC))
	CHCC_to_CH_CC= (CC+ CH - CHCC);

if (CHCC_to_CH_CC <0)
        CHCC_to_CH_CC =0.001;
	
double CHCC_to_CHC_C = p[147] + dEA110 + CHCC_to_CHC_C_H_slope - CHCC_0;
if (CHCC_to_CHC_C < (CHC + C - CHCC))
	CHCC_to_CHC_C = (CHC + C - CHCC);

if (CHCC_to_CHC_C <0)
        CHCC_to_CHC_C =0.001;

double CHCC_to_CCC_H = p[148] + dEA111 + CHCC_to_CCC_H_H_slope - CHCC_0;
if (CHCC_to_CCC_H < (CCC + H - CHCC))
	CHCC_to_CCC_H = (CCC + H - CHCC);

if (CHCC_to_CCC_H <0)
        CHCC_to_CCC_H = 0.001;

double CCC_to_CC_C = p[149] + dEA112 + CCC_to_CC_C_H_slope - CCC_0;
if (CCC_to_CC_C < (C +CC - CCC))
	CCC_to_CC_C = (C +CC - CCC);

if (CCC_to_CC_C <0)
      CCC_to_CC_C = 0.001;

double CH4_to_CH3_H = p[150] + dEA113 + CH4_to_CH3_H_H_slope - CH4_0;
if (CH4_to_CH3_H < (H + CH3 - CH4))
	CH4_to_CH3_H = (H + CH3 - CH4);

if (CH4_to_CH3_H <0)
        CH4_to_CH3_H = 0.001;

double CH3_to_CH2_H = p[151] + dEA114 + CH3_to_CH2_H_H_slope - CH3_0;
if (CH3_to_CH2_H < (H+CH2 - CH3))
	CH3_to_CH2_H = (H+CH2 - CH3);

if (CH3_to_CH2_H <0)
        CH3_to_CH2_H =0.001;

double CH2_to_CH_H = p[152] + dEA115 + CH2_to_CH_H_H_slope - CH2_0;
if (CH2_to_CH_H <( CH + H - CH2))
	CH2_to_CH_H = (CH + H -CH2);

if (CH2_to_CH_H <0)
        CH2_to_CH_H =0.001;

double CH_to_C_H = p[153] + dEA116 + CH_to_C_H_H_slope - CH_0;
if (CH_to_C_H < (C+H - CH))
	CH_to_C_H = (C+H - CH);

if (CH_to_C_H <0)
        CH_to_C_H =0.001;

double CH3CH3_to_CH3_CH3 = p[154] + dEA117 + CH3CH3_to_CH3_CH3_H_slope - CH3CH3_0;
if (CH3CH3_to_CH3_CH3 < (CH3 + CH3 - CH3CH3))
	CH3CH3_to_CH3_CH3 = (CH3 + CH3 - CH3CH3);

if (CH3CH3_to_CH3_CH3 <0)
        CH3CH3_to_CH3_CH3 =0.001;

double CH3CH3_to_CH3CH2_H = p[155] + dEA118 + CH3CH3_to_CH3CH2_H_H_slope - CH3CH3_0;
if (CH3CH3_to_CH3CH2_H < (CH3CH2 + H - CH3CH3))
	CH3CH3_to_CH3CH2_H = (CH3CH2 + H - CH3CH3);

if (CH3CH3_to_CH3CH2_H <0)
        CH3CH3_to_CH3CH2_H =0.001;

double CH3CH2_to_CH3_CH2 = p[156] + dEA119 + CH3CH2_to_CH3_CH2_H_slope - CH3CH2_0;
if (CH3CH2_to_CH3_CH2 < (CH3 + CH2 - CH3CH2))
	CH3CH2_to_CH3_CH2 < (CH3 + CH2 - CH3CH2);

if (CH3CH2_to_CH3_CH2 <0)
        CH3CH2_to_CH3_CH2 =0.001;

double CH3CH2_to_CH3CH_H = p[157] + dEA120 + CH3CH2_to_CH3CH_H_H_slope - CH3CH2_0;
if (CH3CH2_to_CH3CH_H < (CH3CH + H -CH3CH2))
	CH3CH2_to_CH3CH_H = (CH3CH + H -CH3CH2);

if (CH3CH2_to_CH3CH_H <0)
        CH3CH2_to_CH3CH_H =0.001;

double CH3CH2_to_CH2CH2_H = p[158] + dEA121 + CH3CH2_to_CH2CH2_H_H_slope - CH3CH2_0;
if (CH3CH2_to_CH2CH2_H < (CH2CH2 + H - CH3CH2))
	CH3CH2_to_CH2CH2_H = (CH2CH2 + H - CH3CH2);

if (CH3CH2_to_CH2CH2_H <0)
        CH3CH2_to_CH2CH2_H =0.001;

double CH3CH_to_CH3_CH = p[159] + dEA122 + CH3CH_to_CH3_CH_H_slope - CH3CH_0;
if (CH3CH_to_CH3_CH < (CH3 + CH - CH3CH))
	CH3CH_to_CH3_CH = (CH3 + CH - CH3CH);

if (CH3CH_to_CH3_CH <0)
        CH3CH_to_CH3_CH =0.001;

double CH3CH_to_CH3C_H = p[160] + dEA123 + CH3CH_to_CH3C_H_H_slope - CH3CH_0;
if (CH3CH_to_CH3C_H < (CH3C + H - CH3CH))
	CH3CH_to_CH3C_H = (CH3C + H - CH3CH);

if (CH3CH_to_CH3C_H <0)
        CH3CH_to_CH3C_H =0.001;

double CH3CH_to_CH2CH_H = p[161] + dEA124 + CH3CH_to_CH2CH_H_H_slope - CH3CH_0;
if (CH3CH_to_CH2CH_H < (CH2CH +H - CH3CH))
	CH3CH_to_CH2CH_H = (CH2CH +H - CH3CH);

if (CH3CH_to_CH2CH_H<0)
 CH3CH_to_CH2CH_H =0.001;

double CH3C_to_CH3_C = p[162] + dEA125 + CH3C_to_CH3_C_H_slope - CH3C_0;
if (CH3C_to_CH3_C < (C+CH3 +CH3C))
	CH3C_to_CH3_C = (C+CH3 +CH3C);

if (CH3C_to_CH3_C<0)
        CH3C_to_CH3_C =0.001;

double CH3C_to_CH2C_H = p[163] + dEA126 + CH3C_to_CH2C_H_H_slope - CH3C_0;
if (CH3C_to_CH2C_H < (CH2C + H - CH3C))
	CH3C_to_CH2C_H = (CH2C + H - CH3C);

if (CH3C_to_CH2C_H <0)
        CH3C_to_CH2C_H =0.001;

double CH2CH2_to_CH2_CH2 = p[164] + dEA127 + CH2CH2_to_CH2_CH2_H_slope - CH2CH2_0;
if (CH2CH2_to_CH2_CH2 < (CH2 + CH2 - CH2CH2))
	CH2CH2_to_CH2_CH2 = (CH2 + CH2 - CH2CH2);

if (CH2CH2_to_CH2_CH2 <0)
        CH2CH2_to_CH2_CH2 =0.001;

double CH2CH2_to_CH2CH_H = p[165] + dEA128 + CH2CH2_to_CH2CH_H_H_slope - CH2CH2_0;
if (CH2CH2_to_CH2CH_H < (CH2CH+H - CH2CH2))
	CH2CH2_to_CH2CH_H = (CH2CH+H - CH2CH2);

if (CH2CH2_to_CH2CH_H <0)
        CH2CH2_to_CH2CH_H = 0.001;

double CH2CH_to_CH2_CH = p[166] + dEA129 + CH2CH_to_CH2_CH_H_slope - CH2CH_0;
if (CH2CH_to_CH2_CH < (CH2 + CH - CH2CH))
	CH2CH_to_CH2_CH = (CH2 + CH - CH2CH);

if (CH2CH_to_CH2_CH <0)
        CH2CH_to_CH2_CH =0.001;

double CH2CH_to_CH2C_H = p[167] + dEA130 + CH2CH_to_CH2C_H_H_slope - CH2CH_0;
if (CH2CH_to_CH2C_H < (CH2C + H - CH2CH))
	CH2CH_to_CH2C_H = (CH2C + H - CH2CH);

if (CH2CH_to_CH2C_H <0)
        CH2CH_to_CH2C_H =0.001;

double CH2CH_to_CHCH_H = p[168] + dEA131 + CH2CH_to_CHCH_H_H_slope - CH2CH_0;
if (CH2CH_to_CHCH_H < (CHCH +H - CH2CH))
	CH2CH_to_CHCH_H = (CHCH +H - CH2CH);

if (CH2CH_to_CHCH_H <0)
        CH2CH_to_CHCH_H =0.001;

double CH2C_to_CH2_C = p[169] + dEA132 + CH2C_to_CH2_C_H_slope - CH2C_0;
if (CH2C_to_CH2_C < (CH2+C - CH2C))
	CH2C_to_CH2_C = (CH2+C - CH2C);

if (CH2C_to_CH2_C <0)
        CH2C_to_CH2_C =0.001;

double CH2C_to_CHC_H = p[170] + dEA133 + CH2C_to_CHC_H_H_slope - CH2C_0;
if (CH2C_to_CHC_H < (CHC +H - CH2C))
	CH2C_to_CHC_H = (CHC +H - CH2C);

if (CH2C_to_CHC_H <0)
        CH2C_to_CHC_H =0.001;

double CHCH_to_CH_CH = p[171] + dEA134 + CHCH_to_CH_CH_H_slope - CHCH_0;
if (CHCH_to_CH_CH < (CH + CH - CHCH))
	CHCH_to_CH_CH = (CH + CH - CHCH);

if (CHCH_to_CH_CH <0)
        CHCH_to_CH_CH =0.001;

double CHCH_to_CHC_H = p[172] + dEA135 + CHCH_to_CHC_H_H_slope - CHCH_0;
if (CHCH_to_CHC_H < ( CHC + H - CHCH))
	CHCH_to_CHC_H = ( CHC + H - CHCH);

if (CHCH_to_CHC_H <0)
        CHCH_to_CHC_H =0.001;

double CHC_to_CH_C = p[173] +dEA136 + CHC_to_CH_C_H_slope - CHC_0;
if (CHC_to_CH_C < (CH + C - CHC))
	CHC_to_CH_C = (CH + C - CHC);

if (CHC_to_CH_C < 0)
        CHC_to_CH_C =0.001;

double CHC_to_CC_H = p[174] + dEA137 + CHC_to_CC_H_H_slope - CHC_0;
if (CHC_to_CC_H < (CC + H - CHC))
	CHC_to_CC_H = (CC + H - CHC);

if (CHC_to_CC_H <0)
        CHC_to_CC_H = 0.001;

double CC_to_C_C = p[175]  + dEA138 + CC_to_C_C_H_slope - CC_0;
if (CC_to_C_C < ( C+ C - CC))
	CC_to_C_C = ( C+ C - CC);

if (CC_to_C_C <0)
        CC_to_C_C = 0.001;


//All the reaction energies and adsorption energies are ZPE corrected but (!!!) are dE
//K1   = (QVIB_forPROPANE_ADS/(3.847E+01*1.928E+09*6.358E+04))*exp(-(Propane_adsorption_reaction_energy+CH3CH2CH3 - vacancy)/(T*8.6173324E-05));
//The above is orginial, Im changing this guy to be in dG. Yay!

K[0]    = exp(-(CH3CH2CH3 - CH3CH2CH3_gp - alph)/(T*8.6173324E-05));
K[1]    = exp(-(CH3CHCH2 - CH3CHCH2_gp - beta)/(T*8.6173324E-05));
K[2]    = exp(-(2*H - H2_gp -gamm)/(T*8.6173324E-05));
K[3]    = exp(-(CH3CCH - CH3CCH_gp)/(T*8.6173324E-05));
K[4]    = exp(-(CH3CH3 - CH3CH3_gp)/(T*8.6173324E-05));
K[5]    = exp(-(CH2CH2 - CH2CH2_gp)/(T*8.6173324E-05));
K[6]    = exp(-(CHCH - CHCH_gp )/(T*8.6173324E-05));
K[7]    = exp(-(CH4 + CH4_gp)/(T*8.6173324E-05));
K[8]    = exp(-(CH3CHCH3 +H -CH3CH2CH3)/(T*8.6173324E-05));
K[9]   = exp(-(CH3CH2CH2 +H -CH3CH2CH3)/(T*8.6173324E-05));
K[10]   = exp(-(CH3CHCH2+H-CH3CHCH3)/(T*8.6173324E-05));
K[11]   = exp(-(CH3CHCH2+H -CH3CH2CH2)/(T*8.6173324E-05));
K[12]   = exp(-(CH3+CH3CH2 -CH3CH2CH3)/(T*8.6173324E-05));
K[13]   = exp(-(CH3+CH3CH -CH3CHCH3)/(T*8.6173324E-05));
K[14]   = exp(-(CH3CCH3 +H -CH3CHCH3)/(T*8.6173324E-05));
K[15]   = exp(-(CH3CH2+CH2- CH3CH2CH2)/(T*8.6173324E-05));
K[16]   = exp(-(CH3+CH2CH2 -CH3CH2CH2)/(T*8.6173324E-05));
K[17]   = exp(-(CH2CH2CH2+H -CH3CH2CH2)/(T*8.6173324E-05));
K[18]   = exp(-(CH3CH2CH +H -CH3CH2CH2)/(T*8.6173324E-05));
K[19]   = exp(-(CH3CH2+CH-CH3CH2CH)/(T*8.6173324E-05));
K[20]   = exp(-(CH3+CH2CH -CH3CH2CH)/(T*8.6173324E-05));
K[21]   = exp(-(CH3CH2C+H -CH3CH2CH)/(T*8.6173324E-05));
K[22]   = exp(-(CH3CHCH+H -CH3CH2CH)/(T*8.6173324E-05));
K[23]   = exp(-(CH2CH2CH+H -CH3CH2CH)/(T*8.6173324E-05));
K[24]   = exp(-(CH2+CH2CH2 -CH2CH2CH2)/(T*8.6173324E-05));
K[25]   = exp(-(CH2CH2CH+H -CH2CH2CH2)/(T*8.6173324E-05));
K[26]   = exp(-(CH2CHCH2+H -CH2CH2CH2)/(T*8.6173324E-05));
K[27]   = exp(-(CH3+CH2CH  -CH3CHCH2)/(T*8.6173324E-05));
K[28]   = exp(-(CH3CH+CH2  -CH3CHCH2)/(T*8.6173324E-05));
K[29]   = exp(-(CH3CCH2   +H  -CH3CHCH2)/(T*8.6173324E-05));
K[30]   = exp(-(CH3CHCH   +H  -CH3CHCH2)/(T*8.6173324E-05));
K[31]   = exp(-(CH2CHCH2   +H  -CH3CHCH2)/(T*8.6173324E-05));
K[32]   = exp(-(CH3C +CH3  -CH3CCH3)/(T*8.6173324E-05));
K[33]   = exp(-(CH3CCH2+H - CH3CCH3)/(T*8.6173324E-05));
K[34]   = exp(-(CH3 +CH2C - CH3CH2C )/(T*8.6173324E-05));
K[35]   = exp(-(CH3CH2 +C - CH3CH2C )/(T*8.6173324E-05));
K[36]   = exp(-(CH2CH2C + H - CH3CH2C )/(T*8.6173324E-05));
K[37]   = exp(-(CH3CHC +H - CH3CH2C )/(T*8.6173324E-05));
K[38]   = exp(-(CH2+CH2CH - CH2CH2CH )/(T*8.6173324E-05));
K[39]   = exp(-(CH2CH2+CH - CH2CH2CH )/(T*8.6173324E-05));
K[40]   = exp(-(CH2CH2C+H - CH2CH2CH )/(T*8.6173324E-05));
K[41]   = exp(-(CH2CHCH+H - CH2CH2CH )/(T*8.6173324E-05));
K[42]   = exp(-(CHCH2CH+H - CH2CH2CH )/(T*8.6173324E-05));
K[43]   = exp(-(CH2+CH2CH - CH2CHCH2 )/(T*8.6173324E-05));
K[44]   = exp(-(CH2CHCH+H - CH2CHCH2 )/(T*8.6173324E-05));
K[45]   = exp(-(CH2CCH2+H - CH2CHCH2 )/(T*8.6173324E-05));
K[46]   = exp(-(CH3+CHCH - CH3CHCH  )/(T*8.6173324E-05));
K[47]   = exp(-(CH3CH+CH - CH3CHCH  )/(T*8.6173324E-05));
K[48]   = exp(-(CH3CHC+H - CH3CHCH  )/(T*8.6173324E-05));
K[49]   = exp(-(CH3CCH+H - CH3CHCH  )/(T*8.6173324E-05));
K[50]   = exp(-(CH2CHCH+H - CH3CHCH  )/(T*8.6173324E-05));
K[51]   = exp(-(CH3+CH2C - CH3CCH2  )/(T*8.6173324E-05));
K[52]   = exp(-(CH3C+CH2 - CH3CCH2  )/(T*8.6173324E-05));
K[53]   = exp(-(CH2CCH2+H - CH3CCH2  )/(T*8.6173324E-05));
K[54]   = exp(-(CH3CCH+H - CH3CCH2  )/(T*8.6173324E-05));
K[55]   = exp(-(CH3+CHC - CH3CHC)/(T*8.617332E-05));
K[56]   = exp(-(CH3CH+C - CH3CHC)/(T*8.617332E-05));
K[57]   = exp(-(CH3CC+H - CH3CHC)/(T*8.617332E-05));
K[58]   = exp(-(CH2CHC+H - CH3CHC)/(T*8.617332E-05));
K[59]   = exp(-(CH2CH2+C - CH2CH2C)/(T*8.617332E-05));
K[60]   = exp(-(CH2+CH2C - CH2CH2C)/(T*8.617332E-05));
K[61]   = exp(-(CH2CHC+H - CH2CH2C)/(T*8.617332E-05));
K[62]   = exp(-(CHCH2C+H - CH2CH2C)/(T*8.617332E-05));
K[63]   = exp(-(CH2CH+CH - CHCH2CH)/(T*8.617332E-05));
K[64]   = exp(-(CHCH2C+H - CHCH2CH)/(T*8.617332E-05));
K[65]   = exp(-(CHCHCH+H - CHCH2CH)/(T*8.617332E-05));
K[66]   = exp(-(CH2+CHCH- CH2CHCH)/(T*8.617332E-05));
K[67]   = exp(-(CH2CH+CH - CH2CHCH)/(T*8.617332E-05));
K[68]   = exp(-(CH2CHC+H - CH2CHCH)/(T*8.617332E-05));
K[69]   = exp(-(CH2CCH+H - CH2CHCH)/(T*8.617332E-05));
K[70]   = exp(-(CHCHCH+H - CH2CHCH)/(T*8.617332E-05));
K[71]   = exp(-(CH2C+CH2 - CH2CCH2)/(T*8.617332E-05));
K[72]   = exp(-(CH2CCH+H - CH2CCH2)/(T*8.617332E-05));
K[73]   = exp(-(CH3C+CH - CH3CCH)/(T*8.617332E-05));
K[74]   = exp(-(CH3+CHC - CH3CCH)/(T*8.617332E-05));
K[75]   = exp(-(CH3CC+H - CH3CCH)/(T*8.617332E-05));
K[76]   = exp(-(CH2CCH+H - CH3CCH)/(T*8.617332E-05));
K[77]   = exp(-(CH3+CC - CH3CC)/(T*8.617332E-05));
K[78]   = exp(-(CH3C+C - CH3CC)/(T*8.617332E-05));
K[79]   = exp(-(CH2CC+H - CH3CC)/(T*8.617332E-05));
K[80]   = exp(-(CH2+CHC  -CH2CHC)/(T*8.617332E-05));
K[81]   = exp(-(CH2CH+C  -CH2CHC)/(T*8.617332E-05));
K[82]   = exp(-(CH2CC+H -CH2CHC)/(T*8.617332E-05));
K[83]   = exp(-(CHCHC+H  -CH2CHC)/(T*8.617332E-05));
K[84]   = exp(-(CH+CH2C -CHCH2C)/(T*8.617332E-05));
K[85]   = exp(-(CH2CH+C -CHCH2C)/(T*8.617332E-05));
K[86]   = exp(-(CHCHC+H -CHCH2C)/(T*8.617332E-05));
K[87]   = exp(-(CCH2C+H -CHCH2C)/(T*8.617332E-05));
K[88]   = exp(-(CHCH+CH -CHCHCH)/(T*8.617332E-05));
K[89]   = exp(-(CHCHC+H -CHCHCH)/(T*8.617332E-05));
K[90]   = exp(-(CHCCH+H -CHCHCH)/(T*8.617332E-05));
K[91]   = exp(-(CH2+CHC -CH2CCH)/(T*8.617332E-05));
K[92]   = exp(-(CH2C+CH -CH2CCH)/(T*8.617332E-05));
K[93]   = exp(-(CH2CC+H -CH2CCH)/(T*8.617332E-05));
K[94]   = exp(-(CHCCH+H -CH2CCH)/(T*8.617332E-05));
K[95]   = exp(-(CH2+CC -CH2CC)/(T*8.617332E-05));
K[96]   = exp(-(CH2C+C -CH2CC)/(T*8.617332E-05));
K[97]   = exp(-(CHCC+H -CH2CC)/(T*8.617332E-05));
K[98]   = exp(-(CH+CHC -CHCHC)/(T*8.617332E-05));
K[99]  = exp(-(CHCH+C -CHCHC)/(T*8.617332E-05));
K[100]  = exp(-(CHCC+H -CHCHC)/(T*8.617332E-05));
K[101]  = exp(-(CCHC+H -CHCHC)/(T*8.617332E-05));
K[102]  = exp(-(CH2C+C -CCH2C)/(T*8.617332E-05));
K[103]  = exp(-(CCHC+H -CCH2C)/(T*8.617332E-05));
K[104]  = exp(-(CH+CHC -CHCCH)/(T*8.617332E-05));
K[105]  = exp(-(CHCC+H -CHCCH)/(T*8.617332E-05));
K[106]  = exp(-(CHC+C -CCHC)/(T*8.617332E-05));
K[107]  = exp(-(CCC+H-CCHC)/(T*8.617332E-05));
K[108]  = exp(-(CH+CC-CHCC)/(T*8.617332E-05));
K[109]  = exp(-(CHC+C-CHCC)/(T*8.617332E-05));
K[110]  = exp(-(CCC+H-CHCC)/(T*8.617332E-05));
K[111]  = exp(-(C+CC-CCC)/(T*8.617332E-05));
K[112]  = exp(-(CH3+H-CH4)/(T*8.617332E-05));
K[113]  = exp(-(CH2+H-CH3)/(T*8.617332E-05));
K[114]  = exp(-(CH+H-CH2)/(T*8.617332E-05));
K[115]  = exp(-(C+H-CH)/(T*8.617332E-05));
K[116]  = exp(-(CH3+CH3 -CH3CH3)/(T*8.617332E-05));
K[117]  = exp(-(CH3CH2+H -CH3CH3)/(T*8.617332E-05));
K[118]  = exp(-(CH3+CH2 -CH3CH2)/(T*8.617332E-05));
K[119]  = exp(-(CH3CH+H -CH3CH2)/(T*8.617332E-05));
K[120]  = exp(-(CH2CH2+H -CH3CH2)/(T*8.617332E-05));
K[121]  = exp(-(CH3+CH -CH3CH )/(T*8.617332E-05));
K[122]  = exp(-(CH3C+H -CH3CH )/(T*8.617332E-05));
K[123]  = exp(-(CH2CH+H -CH3CH )/(T*8.617332E-05));
K[124]  = exp(-(CH3+C -CH3C  )/(T*8.617332E-05));
K[125]  = exp(-(CH2C+H -CH3C  )/(T*8.617332E-05));
K[126]  = exp(-(CH2+CH2 -CH2CH2  )/(T*8.617332E-05));
K[127]  = exp(-(CH2CH+H -CH2CH2  )/(T*8.617332E-05));
K[128]  = exp(-(CH2+CH -CH2CH  )/(T*8.617332E-05));
K[129]  = exp(-(CH2C+H -CH2CH  )/(T*8.617332E-05));
K[130]  = exp(-(CHCH+H -CH2CH  )/(T*8.617332E-05));
K[131]  = exp(-(CH2+C -CH2C  )/(T*8.617332E-05));
K[132]  = exp(-(CHC+H -CH2C  )/(T*8.617332E-05));
K[133]  = exp(-(CH+CH -CHCH  )/(T*8.617332E-05));
K[134]  = exp(-(CHC+H -CHCH  )/(T*8.617332E-05));
K[135]  = exp(-(CH+C-CHC   )/(T*8.617332E-05));
K[136]  = exp(-(CC+H-CHC   )/(T*8.617332E-05));
K[137]  = exp(-(C+C-CC)/(T*8.617332E-05));


//Old code below
f[0] = f[0];                    
f[1] = f[1];                         
f[2] = f[2];
f[3] = f[3];
f[4] = f[4];
f[5] = f[5];
f[6] = f[6];
f[7] = f[7];

f[8]  = (kB*T/h)*exp(-(CH3CH2CH3_to_CH3CHCH3_H)/(kB*T));
f[9]  = (kB*T/h)*exp(-(CH3CH2CH3_to_CH3CH2CH2_H)/(kB*T));
f[10] = (kB*T/h)*exp(-(CH3CHCH3_to_CH3CHCH2_H/(kB*T)));
f[11] = (kB*T/h)*exp(-(CH3CH2CH2_to_CH3CHCH2_H)/(kB*T));
f[12] = (kB*T/h)*exp(-(CH3CH2CH3_to_CH3CH2_CH3)/(kB*T));
f[13] = (kB*T/h)*exp(-(CH3CHCH3_to_CH3_CH3CH  )/(kB*T));
f[14] = (kB*T/h)*exp(-(CH3CHCH3_to_CH3CCH3_H  )/(kB*T));
f[15] = (kB*T/h)*exp(-(CH3CH2CH2_to_CH3CH2_CH2)/(kB*T));
f[16] = (kB*T/h)*exp(-(CH3CH2CH2_to_CH3_CH2CH2)/(kB*T));
f[17] = (kB*T/h)*exp(-(CH3CH2CH2_to_CH2CH2CH2_H)/(kB*T));
f[18] = (kB*T/h)*exp(-(CH3CH2CH2_to_CH3CH2CH_H)/(kB*T));
f[19] = (kB*T/h)*exp(-(CH3CH2CH_to_CH3CH2_CH)/(kB*T));
f[20] = (kB*T/h)*exp(-(CH3CH2CH_to_CH3_CH2CH)/(kB*T));
f[21] = (kB*T/h)*exp(-(CH3CH2CH_to_CH3CH2C_H)/(kB*T));
f[22] = (kB*T/h)*exp(-(CH3CH2CH_to_CH3CHCH_H)/(kB*T));
f[23] = (kB*T/h)*exp(-(CH3CH2CH_to_CH2CH2CH_H)/(kB*T));
f[24] = (kB*T/h)*exp(-(CH2CH2CH2_to_CH2CH2_CH2)/(kB*T));
f[25] = (kB*T/h)*exp(-(CH2CH2CH2_to_CH2CH2CH_H)/(kB*T));
f[26] = (kB*T/h)*exp(-(CH2CH2CH2_to_CH2CHCH2_H)/(kB*T));
f[27] = (kB*T/h)*exp(-(CH3CHCH2_to_CH3_CHCH2)/(kB*T));
f[28] = (kB*T/h)*exp(-(CH3CHCH2_to_CH3CH_CH2)/(kB*T));
f[29] = (kB*T/h)*exp(-(CH3CHCH2_to_CH3CCH2_H)/(kB*T));
f[30] = (kB*T/h)*exp(-(CH3CHCH2_to_CH3CHCH_H)/(kB*T));
f[31] = (kB*T/h)*exp(-(CH3CHCH2_to_CH2CHCH2_H)/(kB*T));
f[32] = (kB*T/h)*exp(-(CH3CCH3_to_CH3_CCH3   )/(kB*T));
f[33] = (kB*T/h)*exp(-(CH3CCH3_to_CH3CCH2_H  )/(kB*T));
f[34] = (kB*T/h)*exp(-(CH3CH2C_to_CH3_CH2C   )/(kB*T));
f[35] = (kB*T/h)*exp(-(CH3CH2C_to_CH3CH2_C   )/(kB*T));
f[36] = (kB*T/h)*exp(-(CH3CH2C_to_CH2CH2C_H   )/(kB*T));
f[37] = (kB*T/h)*exp(-(CH3CH2C_to_CH3CHC_H   )/(kB*T));
f[38] = (kB*T/h)*exp(-(CH2CH2CH_to_CH2_CH2CH )/(kB*T));
f[39] = (kB*T/h)*exp(-(CH2CH2CH_to_CH2CH2_CH )/(kB*T));
f[40] = (kB*T/h)*exp(-(CH2CH2CH_to_CH2CH2C_H )/(kB*T));
f[41] = (kB*T/h)*exp(-(CH2CH2CH_to_CH2CHCH_H )/(kB*T));
f[42] = (kB*T/h)*exp(-(CH2CH2CH_to_CHCH2CH_H )/(kB*T));
f[43] = (kB*T/h)*exp(-(CH2CHCH2_to_CH2_CH2CH )/(kB*T));
f[44] = (kB*T/h)*exp(-(CH2CHCH2_to_CH2CHCH_H )/(kB*T));
f[45] = (kB*T/h)*exp(-(CH2CHCH2_to_CH2CCH2_H )/(kB*T));
f[46] = (kB*T/h)*exp(-(CH3CHCH_to_CH3_CHCH   )/(kB*T));
f[47] = (kB*T/h)*exp(-(CH3CHCH_to_CH3CH_CH   )/(kB*T));
f[48] = (kB*T/h)*exp(-(CH3CHCH_to_CH3CHC_H   )/(kB*T));
f[49] = (kB*T/h)*exp(-(CH3CHCH_to_CH3CCH_H   )/(kB*T));
f[50] = (kB*T/h)*exp(-(CH3CHCH_to_CH2CHCH_H   )/(kB*T));
f[51] = (kB*T/h)*exp(-(CH3CCH2_to_CH3_CCH2   )/(kB*T));
f[52] = (kB*T/h)*exp(-(CH3CCH2_to_CH3C_CH2   )/(kB*T));
f[53] = (kB*T/h)*exp(-(CH3CCH2_to_CH2CCH2_H  )/(kB*T));
f[54] = (kB*T/h)*exp(-(CH3CCH2_to_CH3CCH_H   )/(kB*T));
f[55] = (kB*T/h)*exp(-(CH3CHC_to_CH3_CHC     )/(kB*T));
f[56] = (kB*T/h)*exp(-(CH3CHC_to_CH3CH_C     )/(kB*T));
f[57] = (kB*T/h)*exp(-(CH3CHC_to_CH3CC_H     )/(kB*T));
f[58] = (kB*T/h)*exp(-(CH3CHC_to_CH2CHC_H    )/(kB*T));
f[59] = (kB*T/h)*exp(-(CH2CH2C_to_CH2CH2_C    )/(kB*T));
f[60] = (kB*T/h)*exp(-(CH2CH2C_to_CH2_CH2C    )/(kB*T));
f[61] = (kB*T/h)*exp(-(CH2CH2C_to_CH2CHC_H    )/(kB*T));
f[62] = (kB*T/h)*exp(-(CH2CH2C_to_CHCH2C_H    )/(kB*T));
f[63] = (kB*T/h)*exp(-(CHCH2CH_to_CHCH2_CH    )/(kB*T));
f[64] = (kB*T/h)*exp(-(CHCH2CH_to_CHCH2C_H    )/(kB*T));
f[65] = (kB*T/h)*exp(-(CHCH2CH_to_CHCHCH_H    )/(kB*T));
f[66] = (kB*T/h)*exp(-(CH2CHCH_to_CH2_CHCH    )/(kB*T));
f[67] = (kB*T/h)*exp(-(CH2CHCH_to_CH2CH_CH    )/(kB*T));
f[68] = (kB*T/h)*exp(-(CH2CHCH_to_CH2CHC_H    )/(kB*T));
f[69] = (kB*T/h)*exp(-(CH2CHCH_to_CH2CCH_H    )/(kB*T));
f[70] = (kB*T/h)*exp(-(CH2CHCH_to_CHCHCH_H    )/(kB*T));
f[71] = (kB*T/h)*exp(-(CH2CCH2_to_CH2C_CH2    )/(kB*T));
f[72] = (kB*T/h)*exp(-(CH2CCH2_to_CH2CCH_H    )/(kB*T));
f[73] = (kB*T/h)*exp(-(CH3CCH_to_CH3C_CH      )/(kB*T));
f[74] = (kB*T/h)*exp(-(CH3CCH_to_CH3_CHC      )/(kB*T));
f[75] = (kB*T/h)*exp(-(CH3CCH_to_CH3CC_H      )/(kB*T));
f[76] = (kB*T/h)*exp(-(CH3CCH_to_CH2CCH_H     )/(kB*T));
f[77] = (kB*T/h)*exp(-(CH3CC_to_CH3_CC        )/(kB*T));
f[78] = (kB*T/h)*exp(-(CH3CC_to_CH3C_C        )/(kB*T));
f[79] = (kB*T/h)*exp(-(CH3CC_to_CH2CC_H        )/(kB*T));
f[80] = (kB*T/h)*exp(-(CH2CHC_to_CH2_CHC       )/(kB*T));
f[81] = (kB*T/h)*exp(-(CH2CHC_to_CH2CH_C       )/(kB*T));
f[82] = (kB*T/h)*exp(-(CH2CHC_to_CH2CC_H       )/(kB*T));
f[83] = (kB*T/h)*exp(-(CH2CHC_to_CHCHC_H       )/(kB*T));
f[84] = (kB*T/h)*exp(-(CHCH2C_to_CH_CH2C       )/(kB*T));
f[85] = (kB*T/h)*exp(-(CHCH2C_to_CH2CH_C       )/(kB*T));
f[86] = (kB*T/h)*exp(-(CHCH2C_to_CHCHC_H       )/(kB*T));
f[87] = (kB*T/h)*exp(-(CHCH2C_to_CCH2C_H       )/(kB*T));
f[88] = (kB*T/h)*exp(-(CHCHCH_to_CHCH_CH       )/(kB*T));
f[89] = (kB*T/h)*exp(-(CHCHCH_to_CHCHC_H       )/(kB*T));
f[90] = (kB*T/h)*exp(-(CHCHCH_to_CHCCH_H       )/(kB*T));
f[91] = (kB*T/h)*exp(-(CH2CCH_to_CH2_CCH       )/(kB*T));
f[92] = (kB*T/h)*exp(-(CH2CCH_to_CH2C_CH       )/(kB*T));
f[93] = (kB*T/h)*exp(-(CH2CCH_to_CH2CC_H       )/(kB*T));
f[94] = (kB*T/h)*exp(-(CH2CCH_to_CHCCH_H       )/(kB*T));
f[95] = (kB*T/h)*exp(-(CH2CC_to_CH2_CC   )/(kB*T));
f[96] = (kB*T/h)*exp(-(CH2CC_to_CH2C_C   )/(kB*T));
f[97] = (kB*T/h)*exp(-(CH2CC_to_CHCC_H   )/(kB*T));
f[98] = (kB*T/h)*exp(-(CHCHC_to_CH_CHC   )/(kB*T));
f[99]  = (kB*T/h)*exp(-(CHCHC_to_CHCH_C   )/(kB*T));
f[100] = (kB*T/h)*exp(-(CHCHC_to_CHCC_H   )/(kB*T));
f[101] = (kB*T/h)*exp(-(CHCHC_to_CCHC_H   )/(kB*T));
f[102] = (kB*T/h)*exp(-(CCH2C_to_CH2C_C   )/(kB*T));
f[103] = (kB*T/h)*exp(-(CCH2C_to_CCHC_H   )/(kB*T));
f[104] = (kB*T/h)*exp(-(CHCCH_to_CHC_CH   )/(kB*T));
f[105] = (kB*T/h)*exp(-(CHCCH_to_CHCC_H   )/(kB*T));
f[106] = (kB*T/h)*exp(-(CCHC_to_CHC_C     )/(kB*T));
f[107] = (kB*T/h)*exp(-(CCHC_to_CCC_H     )/(kB*T));
f[108] = (kB*T/h)*exp(-(CHCC_to_CH_CC     )/(kB*T));
f[109] = (kB*T/h)*exp(-(CHCC_to_CHC_C     )/(kB*T));
f[110] = (kB*T/h)*exp(-(CHCC_to_CCC_H     )/(kB*T));
f[111] = (kB*T/h)*exp(-(CCC_to_CC_C       )/(kB*T));
f[112] = (kB*T/h)*exp(-(CH4_to_CH3_H       )/(kB*T));
f[113] = (kB*T/h)*exp(-(CH3_to_CH2_H       )/(kB*T));
f[114] = (kB*T/h)*exp(-(CH2_to_CH_H       )/(kB*T));
f[115] = (kB*T/h)*exp(-(CH_to_C_H       )/(kB*T));
f[116] = (kB*T/h)*exp(-(CH3CH3_to_CH3_CH3      )/(kB*T));
f[117] = (kB*T/h)*exp(-(CH3CH3_to_CH3CH2_H    )/(kB*T));
f[118] = (kB*T/h)*exp(-(CH3CH2_to_CH3_CH2     )/(kB*T));
f[119] = (kB*T/h)*exp(-(CH3CH2_to_CH3CH_H     )/(kB*T));
f[120] = (kB*T/h)*exp(-(CH3CH2_to_CH2CH2_H     )/(kB*T));
f[121] = (kB*T/h)*exp(-(CH3CH_to_CH3_CH        )/(kB*T));
f[122] = (kB*T/h)*exp(-(CH3CH_to_CH3C_H        )/(kB*T));
f[123] = (kB*T/h)*exp(-(CH3CH_to_CH2CH_H       )/(kB*T));
f[124] = (kB*T/h)*exp(-(CH3C_to_CH3_C          )/(kB*T));
f[125] = (kB*T/h)*exp(-(CH3C_to_CH2C_H         )/(kB*T));
f[126] = (kB*T/h)*exp(-(CH2CH2_to_CH2_CH2      )/(kB*T));
f[127] = (kB*T/h)*exp(-(CH2CH2_to_CH2CH_H      )/(kB*T));
f[128] = (kB*T/h)*exp(-(CH2CH_to_CH2_CH        )/(kB*T));
f[129] = (kB*T/h)*exp(-(CH2CH_to_CH2C_H        )/(kB*T));
f[130] = (kB*T/h)*exp(-(CH2CH_to_CHCH_H        )/(kB*T));
f[131] = (kB*T/h)*exp(-(CH2C_to_CH2_C          )/(kB*T));
f[132] = (kB*T/h)*exp(-(CH2C_to_CHC_H          )/(kB*T));
f[133] = (kB*T/h)*exp(-(CHCH_to_CH_CH          )/(kB*T));
f[134] = (kB*T/h)*exp(-(CHCH_to_CHC_H          )/(kB*T));
f[135] = (kB*T/h)*exp(-(CHC_to_CH_C            )/(kB*T));
f[136] = (kB*T/h)*exp(-(CHC_to_CC_H            )/(kB*T));
f[137] = (kB*T/h)*exp(-(CC_to_C_C              )/(kB*T));

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
for (int i = 0; i < 138; i++) 
{
 b[i] = f[i]/K[i];
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

r[0]   	 = f[0]*PCH3CH2CH3_beta*y1 			- b[0]*y3;      //Propane Adsorption           
r[1]   	 = f[1]*PCH3CHCH2_beta*y1*y1 		- b[1]*y6;  //Propylene Adsorption          
r[2]   	 = f[2]*PH2_beta*y1*y1 				- b[2]*y2*y2; //H2 Adsorption         
r[3]   	 = f[3]*PCH3CCH_beta*y1*y1*y1 		- b[3]*y20;                                                     
r[4]   	 = f[4]*PCH3CH3_beta*y1 			- b[4]*y33;                                 
r[5]   	 = f[5]*PCH2CH2_beta*y1*y1 			- b[5]*y37;                          
r[6]   	 = f[6]*PCHCH_beta*y1*y1*y1 		- b[6]*y40;                         
r[7]   	 = f[7]*PCH4_beta*y1 				- b[7]*y43;                               
r[8]   	 = f[8]*y3*y1 						- b[8]*y4*y2;       // 1 + 1 - 1 - 1                                                                 
r[9]   	 = f[9]*y3*y1 						- b[9]*y5*y2;		// 1 + 1 -1 -1
r[10]  	 = f[10]*y4*y1*y1 				- b[10]*y6*y2;		// 2 + 1 - 1 - 2 
r[11]  	 = f[11]*y5*y1*y1 				- b[11]*y6*y2;		// 2 + 1 -1 -2
r[12]  	 = f[12]*y3*y1 					- b[12]*y34*y44;	// 1 + 1 -1 -1
r[13]  	 = f[13]*y4*y1*y1 				- b[13]*y35*y44;	// 1 + 2 -1 - 2
r[14]  	 = f[14]*y4*y1*y1 				- b[14]*y9*y2;		// 2 + 1 -1 -2
r[15]  	 = f[15]*y5*y1*y1 				- b[15]*y34*y45;	// 2 + 1 -1 -2
r[16]  	 = f[16]*y5*y1*y1 				- b[16]*y37*y44;	// 1 + 2 - 1 -2
r[17]  	 = f[17]*y5*y1*y1 				- b[17]*y8*y2;		// 2 + 1 -1 -2
r[18]  	 = f[18]*y5*y1*y1 				- b[18]*y7*y2;		// 2 + 1 -1 -2
r[19]  	 = f[19]*y7*y1*y1 				- b[19]*y34*y46;			// 3 + 1 -2 -2
r[20]  	 = f[20]*y7*y1*y1 				- b[20]*y38*y44;	// 1 + 3 -2 -2
r[21]  	 = f[21]*y7*y1*y1 				- b[21]*y10*y2;			// 3 + 1 -2 -2
r[22]  	 = f[22]*y7*y1 					- b[22]*y13*y2;			// 2 + 1 -2 -1
r[23]  	 = f[23]*y7*y1*y1 				- b[23]*y11*y2;		// 3 + 1 - 2 -2
r[24]  	 = f[24]*y8*y1*y1 				- b[24]*y37*y45;	// 2 + 2 -2 -2
r[25]  	 = f[25]*y8*y1*y1 				- b[25]*y11*y2;		// 3 + 1 -2 -2
r[26]  	 = f[26]*y8*y1 					- b[26]*y12*y2;		// 2 + 1 -2 -1
r[27]  	 = f[27]*y6*y1*y1    				- b[27]*y38*y44;	// 1 + 3 -2 -2
r[28]  	 = f[28]*y6*y1*y1 				- b[28]*y35*y45;	// 2 + 2 -2 -2
r[29]  	 = f[29]*y6*y1*y1 				- b[29]*y14*y2;		// 3 + 1 -2 -2
r[30]  	 = f[30]*y6*y1 					- b[30]*y13*y2;			// 2 + 1 -2 -1
r[31]  	 = f[31]*y6*y1 					- b[31]*y12*y2;		// 2 + 1 -2 -1
r[32]  	 = f[32]*y9*y1*y1			 		- b[32]*y36*y44;	// 1 + 3 -2 -2
r[33]  	 = f[33]*y9*y1*y1		 			- b[33]*y14*y2;		// 3 + 1 -2 -2
r[34]  	 = f[34]*y10*y1 					- b[34]*y39*y44; 	// 1 + 3 -3 -1
r[35]  	 = f[35]*y10*y1 					- b[35]*y34*y47;	// 3 + 1 - 3
r[36]  	 = f[36]*y10*y1    				- b[36]*y16*y2; 	// 3 + 1 -3 -1
r[37]  	 = f[37]*y10*y1 					- b[37]*y15*y2;		// 3 + 1 -3 -1  
r[38]  	 = f[38]*y11*y1*y1 				- b[38]*y38*y45;	// 2 + 3 -3 -2
r[39]  	 = f[39]*y11*y1*y1 				- b[39]*y37*y46;		// 2 + 3 -3 -2 
r[40]  	 = f[40]*y11*y1 					- b[40]*y16*y2; 	// 3 + 1 -3 -1
r[41]  	 = f[41]*y11*y1   				- b[41]*y18*y2;		// 3 + 1 - 3 -1
r[42]  	 = f[42]*y11*y1   				- b[42]*y17*y2;			// 3 + 1 -3 -1
r[43]  	 = f[43]*y12*y1*y1*y1   			- b[43]*y38*y45;	// 2 + 3 -2 -3
r[44]  	 = f[44]*y12*y1*y1  				- b[44]*y18*y2;			// 3 + 1 -2 -2
r[45]  	 = f[45]*y12*y1*y1  				- b[45]*y19*y2;			// 3 + 1 -2 -2
r[46]  	 = f[46]*y13*y1*y1    			- b[46]*y40*y44;		// 1 + 3 -2 -2
r[47]  	 = f[47]*y13*y1*y1*y1  			- b[47]*y35*y46;		// 2 + 3 -2 -3
r[48]  	 = f[48]*y13*y1*y1    			- b[48]*y15*y2;			// 3 + 1 -2 -2
r[49]  	 = f[49]*y13*y1*y1 				- b[49]*y20*y2;		// 3 + 1 -2 -2
r[50]  	 = f[50]*y13*y1*y1 				- b[50]*y18*y2; 	// 3 + 1 -2 -2
r[51]  	 = f[51]*y14*y1 					- b[51]*y39*y44;	// 1 + 3 -3 -1
r[52]  	 = f[52]*y14*y1*y1		 		- b[52]*y36*y45;	// 2 + 3 -3 -2
r[53]  	 = f[53]*y14*y1		 			- b[53]*y19*y2;		// 3 + 1 -3 -1 
r[54]  	 = f[54]*y14*y1   				- b[54]*y20*y2;		// 3 + 1 -3 -1
r[55]  	 = f[55]*y15*y1		   			- b[55]*y41*y44; 	// 1 + 3 -3 -1
r[56]  	 = f[56]*y15*y1*y1		  		- b[56]*y35*y47; 	// 2 + 3 -3 -2
r[57]  	 = f[57]*y15*y1		 			- b[57]*y21*y2;		// 3 + 1 -3 -1
r[58]  	 = f[58]*y15*y1   				- b[58]*y22*y2; 	// 3 + 1 -3 -1
r[59]  	 = f[59]*y16*y1*y1   				- b[59]*y37*y47;	// 2 + 3 -3 -2 
r[60]  	 = f[60]*y16*y1*y1     			- b[60]*y39*y45;		// 2 + 3 -3 -2
r[61]  	 = f[61]*y16*y1  					- b[61]*y22*y2;			// 3 + 1 -3 -1
r[62]  	 = f[62]*y16*y1*y1   				- b[62]*y24*y2;			// 4 + 1 - 3 -2
r[63]  	 = f[63]*y17*y1*y1*y1   			- b[63]*y38*y46;		// 3 + 3 -3 -3
r[64]  	 = f[64]*y17*y1*y1     			- b[64]*y24*y2; 		// 4 + 1 - 3 -2
r[65]  	 = f[65]*y17*y1     				- b[65]*y23*y2; 		// 3 + 1 -3 -1
r[66]  	 = f[66]*y18*y1*y1		    	- b[66]*y40*y45;		// 2 + 3 -3 -2
r[67]  	 = f[67]*y18*y1*y1*y1    			- b[67]*y38*y46;		// 3 + 3 -3 -3
r[68]  	 = f[68]*y18*y1     				- b[68]*y22*y2;			// 3 + 1 -3 -1 
r[69]  	 = f[69]*y18*y1    				- b[69]*y25*y2;		// 3 + 1 -3 -1
r[70]  	 = f[70]*y18*y1		   			- b[70]*y23*y2;		// 3 + 1 -3 -1
r[71]  	 = f[71]*y19*y1*y1   				- b[71]*y39*y45;	// 3 + 2 -3 -2
r[72]  	 = f[72]*y19*y1     				- b[72]*y25*y2;		// 3 + 1 -3 -1 
r[73]  	 = f[73]*y20*y1*y1*y1	   			- b[73]*y36*y46;		// 3 + 3 -3 -3
r[74]  	 = f[74]*y20*y1    				- b[74]*y41*y44;	// 1 + 3 -3 -1 
r[75]  	 = f[75]*y20*y1    				- b[75]*y21*y2;		// 3 + 1 -3 -1 
r[76]  	 = f[76]*y20*y1        			- b[76]*y25*y2;		// 3 + 1 -3 -1
r[77]  	 = f[77]*y21*y1        			- b[77]*y42*y44;		// 3 + 1 -3 -1
r[78]  	 = f[78]*y21*y1*y1*y1  			- b[78]*y36*y47; 		// 3 + 3 -3 -3
r[79]  	 = f[79]*y21*y1     				- b[79]*y26*y2;		// 3 + 1 -3 -1
r[80]  	 = f[80]*y22*y1*y1		        - b[80]*y41*y45;	// 3 + 2 -3 -2
r[81]  	 = f[81]*y22*y1*y1*y1        		- b[81]*y38*y47;	// 3 + 3 -3 -3
r[82]  	 = f[82]*y22*y1		     		- b[82]*y26*y2;		// 3 + 1 -3 -1
r[83]  	 = f[83]*y22*y1		     		- b[83]*y27*y2;		// 3 +1 -3 -1
r[84]  	 = f[84]*y24*y1*y1		        - b[84]*y39*y46;	// 3 + 3 -4 -2
r[85]  	 = f[85]*y24*y1*y1       			- b[85]*y38*y47;	// 3 + 3 -4 -2
r[86]  	 = f[86]*y24	     				- b[86]*y27*y2; 	// 3 + 1 -4
r[87]  	 = f[87]*y24		   				- b[87]*y28*y2;		// 3 + 1 -4
r[88]  	 = f[88]*y23*y1*y1*y1	   			- b[88]*y40*y46;		// 3 + 3 -3 -3 
r[89]  	 = f[89]*y23*y1     				- b[89]*y27*y2;		// 3 + 1 -3 -1
r[90]  	 = f[90]*y23*y1    				- b[90]*y29*y2;		// 3 + 1 -3 -1
r[91]  	 = f[91]*y25*y1*y1				- b[91]*y41*y45; 	// 3 + 2 -3 -2
r[92]  	 = f[92]*y25*y1*y1*y1	 			- b[92]*y39*y46;	// 3 + 3 -3 -3 
r[93]  	 = f[93]*y25*y1		     		- b[93]*y26*y2;		// 3 + 1 -3 -1
r[94]  	 = f[94]*y25*y1		     		- b[94]*y29*y2;		// 3 + 1 -3 -1
r[95]  	 = f[95]*y26*y1*y1        		- b[95]*y42*y45;	// 2 + 3 -3 -2
r[96]  	 = f[96]*y26*y1*y1*y1      		- b[96]*y39*y47; 	// 3 + 3 -3 -3
r[97]  	 = f[97]*y26*y1   				- b[97]*y31*y2;		// 3 + 1 -3 -1
r[98]  	 = f[98]*y27*y1*y1*y1	   			- b[98]*y41*y46;	// 3 + 3 -3 -3
r[99] 	 = f[99]*y27*y1*y1*y1		 	- b[99]*y40*y47; 	// 3 + 3 -3 -3
r[100] 	 = f[100]*y27*y1   				- b[100]*y31*y2;	// 3 + 1 -3 -1
r[101] 	 = f[101]*y27*y1*y1	 			- b[101]*y30*y2;	// 4 + 1 -3 -2
r[102] 	 = f[102]*y28*y1*y1*y1      		- b[102]*y39*y47;	// 3 + 3 -3 -3
r[103] 	 = f[103]*y28*y1*y1		   		- b[103]*y30*y2;	// 4 + 1 -3 -2 
r[104] 	 = f[104]*y29*y1*y1*y1			- b[104]*y41*y46;	// 3 + 3 -3 -3
r[105] 	 = f[105]*y29*y1   				- b[105]*y31*y2;	// 3 + 1 -3 -1
r[106] 	 = f[106]*y30*y1*y1      			- b[106]*y41*y47;	// 3 + 3 -4 -2
r[107] 	 = f[107]*y30*y1   				- b[107]*y32*y2;	// 4 + 1 -4 -1
r[108] 	 = f[108]*y31*y1*y1*y1			- b[108]*y42*y46;	// 3 + 3 -3 -3
r[109] 	 = f[109]*y31*y1*y1*y1  	   		- b[109]*y41*y47;	// 3 + 3 -3 -3
r[110] 	 = f[110]*y31*y1*y1		        - b[110]*y32*y2;	// 4 + 1 -3 -2
r[111] 	 = f[111]*y32*y1*y1      			- b[111]*y42*y47;	// 3 + 3 - 4 -2
r[112] 	 = f[112]*y43*y1   				- b[112]*y44*y2; 	// 1 + 1 -1 -1  
r[113] 	 = f[113]*y44*y1*y1   			- b[113]*y45*y2;	// 2 + 1 - 1 -2
r[114] 	 = f[114]*y45*y1*y1		   		- b[114]*y46*y2;	// 3 + 1 -2 -2
r[115] 	 = f[115]*y46*y1   				- b[115]*y47*y2; 	// 3 + 1 -3 -1
r[116] 	 = f[116]*y33*y1   				- b[116]*y44*y44;	// 1 + 1 -1 -1
r[117] 	 = f[117]*y33*y1   				- b[117]*y34*y2;  	// 1 + 1 -1 -1
r[118] 	 = f[118]*y34*y1*y1    			- b[118]*y44*y45; 	// 1 + 2 -1 -2
r[119] 	 = f[119]*y34*y1*y1   			- b[119]*y35*y2;  	// 2 + 1 -1 -2
r[120] 	 = f[120]*y34*y1*y1  				- b[120]*y37*y2;	// 2 + 1 -1 -2
r[121] 	 = f[121]*y35*y1*y1		   		- b[121]*y44*y46; 	// 1 + 3 -2 -2
r[122] 	 = f[122]*y35*y1*y1	 			- b[122]*y36*y2;	// 3 + 1 -2 -2 
r[123] 	 = f[123]*y35*y1*y1  				- b[123]*y38*y2;	// 3 + 1 -2 -2
r[124] 	 = f[124]*y36*y1      			- b[124]*y44*y47;	// 1 + 3 - 3 -1
r[125] 	 = f[125]*y36*y1      			- b[125]*y39*y2; 	// 3 + 1 - 3 -1
r[126] 	 = f[126]*y37*y1*y1 				- b[126]*y45*y45;		// 2 + 2 -2 -2
r[127] 	 = f[127]*y37*y1*y1   			- b[127]*y38*y2;		// 3 + 1 -2 -2
r[128] 	 = f[128]*y38*y1*y1			 	- b[128]*y45*y46;		// 2 + 3 -3 -2
r[129] 	 = f[129]*y38*y1 					- b[129]*y39*y2;	// 3 + 1 -3 -1
r[130] 	 = f[130]*y38*y1			 		- b[130]*y40*y2; 	// 3 + 1 -3 -1
r[131] 	 = f[131]*y39*y1*y1		  		- b[131]*y45*y47;	// 2 + 3 - 3 -2
r[132] 	 = f[132]*y39*y1		   			- b[132]*y41*y2;	// 3 +1 -3 -1
r[133] 	 = f[133]*y40*y1*y1*y1		    - b[133]*y46*y46;	// 3 + 3 -3 -3
r[134] 	 = f[134]*y40*y1   				- b[134]*y41*y2;	// 3 + 1 -3 -1
r[135] 	 = f[135]*y41*y1*y1*y1		    - b[135]*y46*y47;	// 3 + 3 -3 -3
r[136] 	 = f[136]*y41*y1   				- b[136]*y42*y2;	// 3 +1 -3 -1
r[137] 	 = f[137]*y42*y1*y1*y1	     	- b[137]*y47*y47;	// 3 + 3 -3 -3


  yd2 = Ith(ydot,2) =  	 2*r[2] +r[8] +r[9] + r[10] + r[11] +r[14] + r[17] +r[18] + r[21] +r[22] +r[23] + r[25] +r[26] + r[29] +r[30] +r[31] +r[33] + r[36] +r[37] +r[40] +r[41] +r[42] +r[44] +r[45] + r[48] + r[49] +r[50]+ r[53] + r[54] +r[57] +r[58] + r[61] + r[62] + r[64] + r[65] + r[68] + r[69] +r[70] + r[72] + r[75] +r[76] + r[79] +r[82] + r[83] + r[86] +r[87] +r[89] +r[90] + r[93] + r[94] + r[97] + r[100] + r[101] + r[103] + r[105] +r[107] + r[110] + r[112] +r[113] +r[114] + r[115] + r[117] + r[119] +r[120] + r[122] + r[123] + r[125] + r[127] + r[129] +r[130] + r[132] + r[134] + r[136]; 
  yd3 = Ith(ydot,3) =  	 r[0] - r[8] -r[9] -r[12]; 
  yd4 = Ith(ydot,4) =  	 r[8] - r[10] -r[13] - r[14];
  yd5 = Ith(ydot,5) =  	 r[9] -r[11] - r[15] -r[16] -r[17] -r[18]; 
  yd6 = Ith(ydot,6) =  	 r[1] + r[10] + r[11] -r[27] -r[28] -r[29] - r[30] -r[31] ; 
  yd7 = Ith(ydot,7) =  	 r[18] - r[19] -r[20] -r[21] -r[22] -r[23]; 
  yd8 = Ith(ydot,8) =  	 r[17] - r[24] -r[25] -r[26]; 
  yd9 = Ith(ydot,9) =  	 r[14] - r[32] -r[33]; 
  yd10 = Ith(ydot,10) =  r[21] - r[34] -r[35] -r[36] -r[37]; 
  yd11 = Ith(ydot,11) =  r[23] +r[25] - r[38] - r[39] -r[40] -r[41] -r[42]; 
  yd12 = Ith(ydot,12) =  r[26] + r[31] - r[43] - r[44] -r[45];
  yd13 = Ith(ydot,13) =  r[22] + r[30] - r[46] -r[47] -r[48] - r[49] -r[50]; 
  yd14 = Ith(ydot,14) =  r[29] + r[33] - r[51] -r[52] -r[53] - r[54]; 
  yd15 = Ith(ydot,15) =  r[37] + r[48] -r[55] -r[56] - r[57] -r[58]; 
  yd16 = Ith(ydot,16) =  r[36] + r[40] - r[59] -r[60] -r[61] - r[62];
  yd17 = Ith(ydot,17) =  r[42] - r[63] - r[64] - r[65] ; 
  yd18 = Ith(ydot,18) =  r[41] +r[44] +r[50] - r[66] -r[67] - r[68] - r[69] - r[70]; 
  yd19 = Ith(ydot,19) =  r[45] +r[53] - r[71] -r[72] ; 
  yd20 = Ith(ydot,20) =  r[3] + r[49] + r[54] - r[73] - r[74] -r[75] - r[76]; 
  yd21 = Ith(ydot,21) =  r[57] + r[75] - r[77] -r[78] - r[79]; 
  yd22 = Ith(ydot,22) =  r[58] + r[61] + r[68] - r[80] -r[81] -r[82] -r[83] ; 
  yd23 = Ith(ydot,23) =  r[65] + r[70] - r[88] - r[89] - r[90]; 
  yd24 = Ith(ydot,24) =  r[62] + r[64] - r[84] -r[85] -r[86] -r[87]; 
  yd25 = Ith(ydot,25) =  r[69] + r[72] + r[76] -r[91] - r[92] - r[93] - r[94]; 
  yd26 = Ith(ydot,26) =  r[79] + r[82] + r[93] - r[95] -r[96] -r[97]; 
  yd27 = Ith(ydot,27) =  r[83] + r[86] +r[89] - r[98] - r[99] -r[100] - r[101];
  yd28 = Ith(ydot,28) =  r[87] - r[102] - r[103]; 
  yd29 = Ith(ydot,29) =  r[90] + r[94] - r[104] - r[105]; 
  yd30 = Ith(ydot,30) =  r[101] + r[103] - r[106] -r[107]; 
  yd31 = Ith(ydot,31) =  r[97] + r[100] + r[105] - r[108] -r[109] -r[110]; 
  yd32 = Ith(ydot,32) =  r[107] + r[110] - r[111]; 
  yd33 = Ith(ydot,33) =  r[4] - r[116] - r[117]; 
  yd34 = Ith(ydot,34) =  r[12] + r[15] + r[19] + r[35] + r[117] -r[118] -r[119] -r[120] ;
  yd35 = Ith(ydot,35) =  r[13] + r[28] +  r[47] + r[56] + r[119] - r[121] -r[122] -r[123];
  yd36 = Ith(ydot,36) =  r[32] + r[52] + r[73] + r[78] +r[122] - r[124]- r[125]; 
  yd37 = Ith(ydot,37) =  r[5] + r[16] +r[24] +r[39] + r[59] +r[120] - r[126] -r[127]; 
  yd38 = Ith(ydot,38) =  r[20] +r[27] + r[38] +r[43] + r[63] +r[67] + r[81] + r[85] +r[123] + r[127] - r[128] - r[129] -r[130]; 
  yd39 = Ith(ydot,39) =  r[34] + r[51] +r[60] + r[71] +r[84] + r[92] + r[96] + r[102] + r[125] + r[129] - r[131] -r[132]; 
  yd40 = Ith(ydot,40) =  r[6] + r[46] +r[66] + r[88] + r[99] + r[130] - r[133] - r[134]; 
  yd41 = Ith(ydot,41) =  r[55] + r[74] +r[80] + r[91] + r[98] + r[104] + r[106] + r[109] + r[132] + r[134] -r[135] -r[136];
  yd42 = Ith(ydot,42) =  r[77] + r[95] + r[108] + r[111] + r[136] - r[137]; 
  yd43 = Ith(ydot,43) =  r[7] - r[112]; 
  yd44 = Ith(ydot,44) =  r[12] + r[13] + r[16] + r[20] +r[27] + r[32] +r[34] +r[46] + r[51] +r[55] + r[74] + r[77] + r[112] -r[113] + 2*r[116] + r[118] + r[121] +r[124];
  yd45 = Ith(ydot,45) =  r[15] +r[24] + r[28] + r[38] +r[43] + r[52] + r[60] + r[66] + r[71] + r[80] + r[91] + r[95] +r[113]  -r[114] + r[118] + 2*r[126] + r[128] + r[131];
  yd46 = Ith(ydot,46) =  r[19] + r[39] + r[47] + r[63] + r[67] + r[73] +r[84] + r[88] + r[92]  + r[98] + r[104] + r[108] + r[114] - r[115] + r[121] + r[128] + 2*r[133] + r[135];
  yd47 = Ith(ydot,47) =  r[35] + r[56] + r[59] + r[78] + r[81] + r[85] + r[96] + r[99] +r[102] + r[106] + r[109] +r[111] + r[115] + r[124] + r[131] + r[135] + 2*r[137];
  yd1 = Ith(ydot,1) = -(yd2 + yd3 +yd4 + yd5 + 2*yd6 + 2*yd7 + 2*yd8 + 2*yd9 + 3*yd10  + 3*yd11  + 2*yd12 + 2*yd13 + 3*yd14 + 3*yd15 +  3*yd16 + 3*yd17 +  3*yd18 + 3*yd19 + 3*yd20 + 3*yd21 + 3*yd22 + 3*yd23 + 4*yd24 + 3*yd25 + 3*yd26 + 3*yd27 + 3*yd28 + 3*yd29 + 4*yd30 + 3*yd31 + 4*yd32 + yd33 + yd34 + 2*yd35 + 3*yd36 + 2*yd37 + 3*yd38 + 3*yd39 + 3*yd40 + 3*yd41 + 3*yd42 + yd43 + yd44 + 2*yd45 + 3*yd46 + 3*yd47) ;
  
for (int i = 0; i< 138; i++)
  {rr[i] = r[i];
   ff[i] = f[i];
   bb[i] = b[i];
   KK[i] = K[i];
 }

/*if (t > 500000)
{
std::cout<< " f0 from code " << f[0] << '\n' ;
std::cout<< std::setprecision (15) << " f0y1PCH3CH2CH3 " << f[0]*y1*PCH3CH2CH3 << '\n' ;
std::cout<< std::setprecision (15) << " b0y3 " << b[0]*y3 << '\n' ;
std::cout<< std::setprecision (15) << " y3 " << y3 << '\n' ;
std::cout<< std::setprecision (15) << " y1 " << y1 << '\n' ;
}*
if (t > 2000) {
std::cout << "Time = " << t << '\n';
std::cout << "Propylene TOF = " << r[1] << '\n';
std::cout << "Propylene TOF using rr = " << rr[1] << '\n';
std::cout << "Propane TOF = " << r[0]  << '\n';
std::cout << "Propane TOF using rr = " << rr[0]  << '\n';
std:: cout << "Site Balance = " << y1 + y2 + y3 +y4 + y5 + 2*y6 + 2*y7 + 2*y8 + y9 + 4*y10  + 3*y11  + 2*y12 + 3*y13 + 3*y14 + 3*y15 +  4*y16 + 4*y17 +  3*y18 + 3*y19 + 4*y20 + 4*y21 + 3*y22 + 4*y23 + 4*y24 + 3*y25 + 4*y26 + 4*y27 + 4*y28 + 4*y29 + 6*y30 + 4*y31 + 6*y32 + y33 + y34 + 2*y35 + 4*y36 + 2*y37 + 2*y38 + 3*y39 + 4*y40 + 4*y41 + 4*y42 + y43 + y44 + 2*y45 + 4*y46 + 4*y47 << '\n' ;
}
*/

  return(0);
}

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errlag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errlag = (int *) flagvalue;
    if (*errlag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errlag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  return(0);
}


// Construtor
likelihoodRoutine_DataClass::likelihoodRoutine_DataClass(const uqBaseEnvironmentClass& env)
  :
  m_heights(0),
  m_times  (0),
  m_stdDevs(0),
  m_env    (&env)
{
  // Data available in /inputData/data02.dat
  double const heights[] = {10,20,30,40,50,60,70,80,90,100,110,120,130,140};
  double const times  [] = {1.41,2.14,2.49,2.87,3.22,3.49,3.81,4.07,4.32,4.47,
      4.75,4.99,5.16,5.26};
  double const stdDevs[] = {0.020,0.120,0.020,0.010,0.030,0.010,0.030,0.030,
      0.030,0.050,0.010,0.040,0.010,0.09};

  std::size_t const n = sizeof(heights)/sizeof(*heights);
  m_heights.assign(heights, heights + n);
  m_times.assign  (times,   times   + n);
  m_stdDevs.assign(stdDevs, stdDevs + n);
}

// Destructor
likelihoodRoutine_DataClass::~likelihoodRoutine_DataClass()
{
}

//------------------------------------------------------
// The user defined likelihood routine
//------------------------------------------------------
double likelihoodRoutine(
  const uqGslVectorClass& paramValues,
  const uqGslVectorClass* paramDirection,
  const void*             functionDataPtr,
  uqGslVectorClass*       gradVector,
  uqGslMatrixClass*       hessianMatrix,
  uqGslVectorClass*       hessianEffect)
{
  const uqBaseEnvironmentClass& env = *(((likelihoodRoutine_DataClass*) functionDataPtr)->m_env);

  if (paramDirection && functionDataPtr && gradVector && hessianMatrix && hessianEffect)
  {
    // Just to eliminate INTEL compiler warnings
  }

  env.subComm().Barrier();

  // The user, at the aplication level, should have set
  // the vector 'paramValues' to have size 1.
  /*UQ_FATAL_TEST_MACRO(paramValues.sizeGlobal() != 1,
                      env.fullRank(),
                      "likelihoodRoutine()",
                      "paramValues vector does not have size 1");*/

  // Compute likelihood

 for (int i = 0; i < 176; i++){
  p[i] = paramValues[i];
   }
 T_beta = 792;
double T = T_beta;

PCH3CH2CH3_beta = 0.29;
PH2_beta = 0.09;
PCH3CHCH2_beta = 0.0;
PCH3CCH_beta = 0.0;
PCH3CH3_beta = 0.0;
PCH2CH2_beta = 0.0;
PCHCH_beta = 0.0;
PCH4_beta = 0.0;

  ff[0] = 1.127655298E+08 ;
  ff[1] = 1.154403922E+08 ;
  ff[2] = 5.295179191E+08 ;
  ff[3] = 1.183093444E+08 ;
  ff[4] = 1.365617097E+08 ;
  ff[5] = 1.413934289E+08 ;
  ff[6] = 1.467490064E+08 ;
  ff[7] = 1.869792775E+08 ;

CH3CH2CH3_energy       	=	2.59050941356 ;
CH3CHCH3_energy        	=	3.06257564721 ;
CH3CH2CH2_energy       	=	3.04149164775 ;
CH3CHCH2_energy        	=	3.21962913029 ;
CH3CH2CH_energy        	=	3.24075713351 ;
CH2CH2CH2_energy       	=	3.60739964562 ;
CH3CCH3_energy 			=	3.49368951336 ;	
CH3CH2C_energy 			=	3.34844550814 ;	
CH2CH2CH_energy        	=	3.83089121638 ;
CH2CHCH2_energy        	=	3.33450679383 ;
CH3CHCH_energy 			=	3.46601689070 ;	
CH3CCH2_energy 			=	3.63743961011 ;	
CH3CHC_energy  			=	3.75835111123 ;	
CH2CH2C_energy 			=	4.18141639387 ;	
CHCH2CH_energy 			=	4.23001097089 ;	
CH2CHCH_energy 			=	3.74221380370 ;	
CH2CCH2_energy 			=	3.87954575181 ;	
CH3CCH_energy  			=	3.40702690846 ;	
CH3CC_energy   			=	4.03865764877 ;	
CH2CHC_energy  			=	4.32050706092 ;	
CHCHCH_energy  			=	4.06952672765 ;	
CHCH2C_energy  			=	4.99394092605 ;	
CH2CCH_energy  			=	4.19263027264 ;	
CH2CC_energy   			=	4.60475210104 ;	
CHCHC_energy   			=	4.65067462672 ;	
CCH2C_energy   			=	5.74308284958 ;	
CHCCH_energy   			=	4.75617276922 ;	
CCHC_energy    			=	5.95807597302 ;	
CHCC_energy    			=	5.46291407737 ;
CCC_energy     			=	6.86211436388 ;	
CH3CH3_energy  			=	1.37999228307 ;
CH3CH2_energy  			=	1.89646618424 ;
CH3CH_energy   			=	2.18709731489 ;
CH3C_energy    			=	2.30188435704 ;
CH2CH2_energy  			=	2.19093425219 ;
CH2CH_energy   			=	2.72140743972 ;
CH2C_energy    			=	2.81804174750 ;
CHCH_energy    			=	2.54746085247 ;
CHC_energy     			=	3.24187638157 ;
CC_energy      			=	4.19389439833 ;
CH4_energy     			=	0.16252677405 ;
CH3_energy     			=	0.77370877934 ;
CH2_energy     			=	1.20089328362 ;
CH_energy      			=	1.55315225512 ;
C_energy       			=	1.96054251724 ;
H_energy       			=	-0.30691026143 ;
			
dEA9	=	3.51052537313 ;
dEA10	=	3.60513565578 ;
dEA11	=	3.95936667611 ;
dEA12	=	3.91122562919 ;
dEA13	=	5.21504744419 ;
dEA14	=	4.95640810023 ;
dEA15	=	3.72862222835 ;
dEA16	=	4.72699907189 ;
dEA17	=	5.42571372965 ;
dEA18	=	4.06938435205 ;
dEA19	=	3.77904801802 ;
dEA20	=	4.98220393588 ;
dEA21	=	5.55539950149 ;
dEA22	=	3.88372898416 ;
dEA23	=	4.09964310857 ;
dEA24	=	4.40902152824 ;
dEA25	=	5.12577699072 ;
dEA26	=	4.35074520989 ;
dEA27	=	4.29782855361 ;
dEA28	=	4.74522419524 ;
dEA29	=	5.15946622269 ;
dEA30	=	3.99859969734 ;
dEA31	=	4.07366316766 ;
dEA32	=	4.10623041730 ;
dEA33	=	5.04384382582 ;
dEA34	=	4.27669526278 ;
dEA35	=	5.48114724246 ;
dEA36	=	5.18983085880 ;
dEA37	=	4.80557054576 ;
dEA38	=	4.33593716552 ;
dEA39	=	5.58100603535 ;
dEA40	=	5.13414423244 ;
dEA41	=	4.89243819260 ;
dEA42	=	4.62498507753 ;
dEA43	=	4.65835118871 ;
dEA44	=	5.82792568722 ;
dEA45	=	4.44346908456 ;
dEA46	=	5.65912738600 ;
dEA47	=	4.99423965625 ;	
dEA48	=	5.51547190885 ;	
dEA49	=	4.38662494673 ;	
dEA50	=	4.71391468257 ;	
dEA51	=	4.82438389880 ;	
dEA52	=	5.15607978429 ;	
dEA53	=	4.62613596535 ;	
dEA54	=	4.63609663406 ;	
dEA55	=	4.59068306799 ;	
dEA56	=	5.87358966889 ;	
dEA57	=	5.71785785007 ;	
dEA58	=	4.78816160643 ;	
dEA59	=	5.14209312040 ;	
dEA60	=	5.42506808202 ;	
dEA61	=	5.54968763877 ;	
dEA62	=	5.41030879182 ;	
dEA63	=	5.45875876344 ;	
dEA64	=	6.15740646496 ;	
dEA65	=	5.66989710230 ;	
dEA66	=	6.24175957935 ;	
dEA67	=	5.52602340740 ;	
dEA68	=	5.88710131624 ;	
dEA69	=	5.02829989336 ;	
dEA70	=	4.86470626452 ;	
dEA71	=	5.02806847170 ;	
dEA72	=	6.10740846594 ;	
dEA73	=	4.57189316860 ;	
dEA74	=	5.27319392889 ;	
dEA75	=	5.54797459974 ;	
dEA76	=	5.21648290014 ;	
dEA77	=	5.01669585444 ;	
dEA78	=	6.51748051292 ;	
dEA79	=	6.01151620022 ;	
dEA80	=	5.47069553812 ;	
dEA81	=	6.07045232631 ;	
dEA82	=	6.29326998784 ;	
dEA83	=	5.51502362674 ;	
dEA84	=	5.30542818293 ;	
dEA85	=	6.65830810223 ;	
dEA86	=	6.12044673613 ;	
dEA87	=	6.26947430170 ;	
dEA88	=	6.20321657277 ;	
dEA89	=	6.19955099260 ;	
dEA90	=	5.55486048656 ;	
dEA91	=	6.27888261723 ;	
dEA92	=	5.95929916350 ;	
dEA93	=	6.37288836472 ;	
dEA94	=	5.66842873437 ;	
dEA95	=	5.88670207940 ;	
dEA96	=	6.51461316990 ;	
dEA97	=	6.26200851762 ;	
dEA98	=	6.66333041598 ;	
dEA99	=	6.43612979715 ;	
dEA100	=	6.35132299667 ;	
dEA101	=	6.71337029567 ;	
dEA102	=	6.78047342508 ;	
dEA103	=	6.71104672185 ;	
dEA104	=	6.89337297783 ;	
dEA105	=	6.25895159503 ;	
dEA106	=	6.19855701647 ;	
dEA107	=	6.55485592476 ;	
dEA108	=	6.34971136488 ;	
dEA109	=	6.92941980690 ;	
dEA110	=	6.32351799777 ;	
dEA111	=	7.00832305621 ;	
dEA112	=	7.26743411838 ;	
dEA113	=	1.28748131620 ;	
dEA114	=	1.55549727823 ;	
dEA115	=	2.13194864417 ;	
dEA116	=	2.64746708498 ;	
dEA117	=	4.03084486646 ;	
dEA118	=	2.50995852843 ;	
dEA119	=	3.54663189783 ;	
dEA120	=	2.68345770838 ;	
dEA121	=	2.72167917220 ;	
dEA122	=	3.99481116212 ;	
dEA123	=	3.13149872676 ;	
dEA124	=	3.20614048867 ;	
dEA125	=	4.04361524981 ;	
dEA126	=	3.56637299856 ;	
dEA127	=	3.97913196393 ;	
dEA128	=	3.05577314317 ;	
dEA129	=	4.59799889613 ;	
dEA130	=	3.47639446709 ;	
dEA131	=	3.55401358557 ;	
dEA132	=	5.03855264325 ;	
dEA133	=	4.01311351286 ;	
dEA134	=	4.49945354151 ;	
dEA135	=	3.89200580568 ;	
dEA136	=	5.30086954234 ;	
dEA137	=	4.66296880060 ;	
dEA138	=	6.20053453537 ;	



  double a = 0.15;
  double ba;
  ba = 0.04- a*(2);


  double alphaa = paramValues[180];
  double betaa  = paramValues[181];
  double gammaa = paramValues[182];
/*  double delta  = paramValues[179]; //These and below are unconstrained. Should  they  be constrained? //
  double epsi   = paramValues[180];
  double zeta   = paramValues[181];
  double eta    = paramValues[182];
  double theta  = paramValues[183];
*/ // I turned this off for math reasons - I'm setting this to 0 until I figure out whats going on
  alph = alphaa  / (alphaa + betaa + gammaa) ;
  beta = betaa    / (alphaa + betaa + gammaa) ;
  gamm = gammaa  / (alphaa + betaa + gammaa) ;

  double delta = 0;
  double epsi = 0;
  double zeta = 0;
  double eta  = 0;
  double theta = 0;
  alph =-1*( alph*(ba-a) + a);
  beta  = (beta*(ba-a) + a);
  gamm= (gamm*(ba-a) + a);


/*  alph = 0;
    beta = 0;
    gamma = 0;
    delta = 0.43;
  */
  realtype reltol, t, tout;
  N_Vector y, abstol;
  void *cvode_mem;
  int flag, flagr, iout;
  int rootsfound[2];

  y = abstol = NULL;
  cvode_mem = NULL;

    /* Create serial vector of length NEQ for I.C. and abstol */
  y = N_VNew_Serial(NEQ);
  abstol = N_VNew_Serial(NEQ);
  /* Initialize y */
  Ith(y,1) = YY1;
  Ith(y,2) = YY2;
  Ith(y,3) = YY3;
  Ith(y,4) = YY4;
  Ith(y,5) = YY5;
  Ith(y,6) = YY6;
  Ith(y,7) = YY7;
  Ith(y,8) = YY8;
  Ith(y,9) = YY9;

  Ith(y,10) = YY10;
  Ith(y,11) = YY11;
  Ith(y,12) = YY12;
  Ith(y,13) = YY13;
  Ith(y,14) = YY14;
  Ith(y,15) = YY15;
  Ith(y,16) = YY16;
  Ith(y,17) = YY17;
  Ith(y,18) = YY18;
  Ith(y,19) = YY19;

  Ith(y,20) = YY20;
  Ith(y,21) = YY21;
  Ith(y,22) = YY22;
  Ith(y,23) = YY23;
  Ith(y,24) = YY24;
  Ith(y,25) = YY25;
  Ith(y,26) = YY26;
  Ith(y,27) = YY27;
  Ith(y,28) = YY28;
  Ith(y,29) = YY29;

  Ith(y,30) = YY30;
  Ith(y,31) = YY31;
  Ith(y,32) = YY32;
  Ith(y,33) = YY33;
  Ith(y,34) = YY34;
  Ith(y,35) = YY35;
  Ith(y,36) = YY36;
  Ith(y,37) = YY37;
  Ith(y,38) = YY38;
  Ith(y,39) = YY39;

  Ith(y,40) = YY40;
  Ith(y,41) = YY41;
  Ith(y,42) = YY42;
  Ith(y,43) = YY43;
  Ith(y,44) = YY44;
  Ith(y,45) = YY45;
  Ith(y,46) = YY46;
  Ith(y,47) = YY47;



  /* Set the scalar relative tolerance */
  reltol = RTOL;
  /* Set the vector absolute tolerance */
  Ith(abstol,1) = ATOL1;
  Ith(abstol,2) = ATOL2;
  Ith(abstol,3) = ATOL3;
  Ith(abstol,4) = ATOL4;
  Ith(abstol,5) = ATOL5;
  Ith(abstol,6) = ATOL6;
  Ith(abstol,7) = ATOL7;
  Ith(abstol,8) = ATOL8;
  Ith(abstol,9) = ATOL9;

  Ith(abstol,10) = ATOL10;
  Ith(abstol,11) = ATOL11;
  Ith(abstol,12) = ATOL12;
  Ith(abstol,13) = ATOL13;
  Ith(abstol,14) = ATOL14;
  Ith(abstol,15) = ATOL15;
  Ith(abstol,16) = ATOL16;
  Ith(abstol,17) = ATOL17;
  Ith(abstol,18) = ATOL18;
  Ith(abstol,19) = ATOL19;

  Ith(abstol,20) = ATOL20;
  Ith(abstol,21) = ATOL21;
  Ith(abstol,22) = ATOL22;
  Ith(abstol,23) = ATOL23;
  Ith(abstol,24) = ATOL24;
  Ith(abstol,25) = ATOL25;
  Ith(abstol,26) = ATOL26;
  Ith(abstol,27) = ATOL27;
  Ith(abstol,28) = ATOL28;
  Ith(abstol,29) = ATOL29;

  Ith(abstol,30) = ATOL30;
  Ith(abstol,31) = ATOL31;
  Ith(abstol,32) = ATOL32;
  Ith(abstol,33) = ATOL33;
  Ith(abstol,34) = ATOL34;
  Ith(abstol,35) = ATOL35;
  Ith(abstol,36) = ATOL36;
  Ith(abstol,37) = ATOL37;
  Ith(abstol,38) = ATOL38;
  Ith(abstol,39) = ATOL39;

  Ith(abstol,40) = ATOL40;
  Ith(abstol,41) = ATOL41;
  Ith(abstol,42) = ATOL42;
  Ith(abstol,43) = ATOL43;
  Ith(abstol,44) = ATOL44;
  Ith(abstol,45) = ATOL45;
  Ith(abstol,46) = ATOL46;
  Ith(abstol,47) = ATOL47;
  /* Set the scalar relative tolerance */
/*reltol = RTOL;
*//* Set the vector absolute tolerance */
/*Ith(abstol,1) = ATOL1;
  Ith(abstol,2) = ATOL2;
  Ith(abstol,3) = ATOL3;
  Ith(abstol,4) = ATOL4;
  Ith(abstol,5) = ATOL5;
  Ith(abstol,6) = ATOL6;
 */
  /* Call CVodeCreate to create the solver memory and specify the
 *  *    * Backward Differentiation Formula and the use of a Newton iteration */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
 //if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  /* Call CVodeInit to initialize the integrator memory and specify the
      *    * user's right hand side function in y'=f(t,y), the inital time T0, and
     *       * the initial dependent variable vector y. */
       flag = CVodeInit(cvode_mem, f, TT0, y);
    // if (check_flag(&flag, "CVodeInit", 1)) return(1);


   /* Call CVodeSVtolerances to specify the scalar relative tolerance
    *    * and vector absolute tolerances */
    flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
      // if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

   /* Call CVDense to specify the CVDENSE dense linear solver */
    flag = CVDense(cvode_mem, NEQ);
    // if (check_flag(&flag, "CVDense", 1)) return(1);
   flag = CVodeSetMaxNumSteps(cvode_mem,50000);

    /* Set the Jacobian routine to Jac (user-suplied) */
    flag = CVDlsSetDenseJacFn(cvode_mem, NULL);
    // if (check_flag(&flag, "CVDlsSetDenseJacFn", 1)) return(1);

     iout = 0;  tout = TT1;
     while(1) {
     flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);

     if (check_flag(&flag, "CVode", 1)) break;
     if (flag == CV_SUCCESS) {
       iout++;
       tout *= TMULT;
     }

     if (iout == NOUT) break;
   }


T = T_beta; double kB = 8.617385e-5;
double  PCH3CH2CH3 = PCH3CH2CH3_beta; double  PCH3CHCH2 = PCH3CHCH2_beta; double  PH2 = PH2_beta; double  PCH3CCH = PCH3CCH_beta; double PCH3CH3 = PCH3CH3_beta;double PCH2CH2 = PCH2CH2_beta;double PCHCH = PCHCH_beta;double PCH4 = PCH4_beta;

double rate_723 = (-rr[1]);
double propane_rate_723 = (rr[0]);
double selectivity_723 = rate_723/propane_rate_723;


//std::cout << "New Set of Data " << '\n';
   ///////////////////////////////////////////////////////////////////////////////////////////
  /*  double    PPCH3CH2CH3[] = {1E-01,4.00,70 };
    double    rates[3];

    for (int i=0; i<3; i++){
            flag = CVodeInit(cvode_mem, f, TT0, y);
            flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
            flag = CVDense(cvode_mem, NEQ);
            flag = CVodeSetMaxNumSteps(cvode_mem,50000);
            PCH3CH2CH3_beta = PPCH3CH2CH3[i];
            iout = 0;  tout = TT1;
            while(1) {
                    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);

                    if (check_flag(&flag, "CVode", 1)) break;
                    if (flag == CV_SUCCESS) {
                            iout++;
                            tout *= TMULT;
                    }

                    if (iout == NOUT) break;
            }

     rates[i] =log(rr[10]+rr[11]-rr[27] -rr[28] -rr[29] - rr[30] -rr[31]);
    };

    int n=3;
    double P[3] = {log(1E-01),log(4.00),log(70)};
    double c0, c1, cov00, cov01, cov11, sumsq;
    gsl_fit_linear(P, 1, rates, 1, n, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
    double CH3CH2CH3r = c1;
    //std::cout << "Propane Reaction Order =" << CH3CH2CH3r << '\n';
    PCH3CH2CH3_beta = 0.04;


 //   double ap_act = -1000*kB*c1;
    /////////////////////////////////////////////////////////////////////////////////////////////


    /////////////////////////////////////////////////////////////////////////////////////////////


    //double PPCO2[] = {0.01, 0.05, 0.10};
    double PPH2[] = {1E-05,1E-01,1000 };

    for (int i=0; i<3; i++){
            flag = CVodeInit(cvode_mem, f, TT0, y);
            flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
            flag = CVodeSetMaxNumSteps(cvode_mem,50000);
            flag = CVDense(cvode_mem, NEQ);
            PH2_beta = PPH2[i];
            iout = 0;  tout = TT1;
            while(1) {
                    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);

                    if (check_flag(&flag, "CVode",1)) break;
                    if (flag == CV_SUCCESS) {
                            iout++;      tout *= TMULT;    }

                    if (iout == NOUT) break;
            }
rates[i] = log(rr[10]+rr[11]-rr[27] -rr[28] -rr[29] - rr[30] -rr[31]); 

    };

    double fizz[3] = {log(1E-05), log(1E-01), log(1000  )};
    gsl_fit_linear(fizz, 1, rates, 1, n, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
    double H2r = c1;

    //std::cout << "H2 Reaction Order =" << H2r << '\n';
    PH2_beta = 2.00;


   //ouble r2 = -1000*kB*c1;

*/
    ////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////////
/*
    T_beta=793;
    T = T_beta;

 ff[0] = 1.081839501E+08;
  ff[1] = 2.215002697E+08;
  ff[2] = 5.080039995E+08;
  ff[3] = 4.540100944E+08;
  ff[4] = 1.310133089E+08;
  ff[5] = 2.712974378E+08;
  ff[6] = 5.631468130E+08;
  ff[7] = 1.793824484E+08;


  ff_Pt111[0] = 1.065504212E+08; //These ar_Pt111e all r_Pt111ef_Pt111er_Pt111enced at 1 b_Pt111ar_Pt111, which is why_Pt111_ PCH3CH2CH3_Pt111 and other_Pt111 pr_Pt111essur_Pt111es ar_Pt111e unitless, as they_Pt111_'ve also b_Pt111een r_Pt111ef_Pt111er_Pt111enced at 1 b_Pt111ar_Pt111.
  ff_Pt111[1]= 2.181557154E+08; //The units of_Pt111 the forwar_Pt111d r_Pt111ates ar_Pt111e 1/s (JFYI)
  ff_Pt111[2] = 5.003333679E+08;
  ff_Pt111[3] = 3.353660599E+08;
  ff_Pt111[4] = 1.290350669E+08;
  ff_Pt111[5] = 2.672009687E+08;
  ff_Pt111[6] = 4.159826624E+08;
  ff_Pt111[7] = 1.766738542E+08;



  ff_Pt211[0] = 8.80581993E+07;
  ff_Pt211[1] = 1.80293980E+08;
  ff_Pt211[2] = 4.13498651E+08;
  ff_Pt211[3] = 2.77162033E+08;
  ff_Pt211[4] = 1.06640551E+08;
  ff_Pt211[5] = 2.20827247E+08;
  ff_Pt211[6] = 3.43787324E+08;
  ff_Pt211[7] = 1.46011450E+08;

CH3CH2CH3_energy=1.36057901687381;dEA9=1.94061221446164;CH3CH2CH3_Pt111_energy=1.23788405671863;dEA_Pt111_9=2.39180974258212;CH3CH2CH3_Pt211_energy =1.46904369227322;dEA_Pt211_9=2.31690222995425;
CH3CHCH3_energy=0.892242587559282;dEA10=2.05496359746526;CH3CHCH3_Pt111_energy=1.03841872106281;dEA_Pt111_10=2.23382372158212;CH3CHCH3_Pt211_energy  =1.18361219580048;dEA_Pt211_10=2.34388222995425;
CH3CH2CH2_energy=0.939423227183711;dEA11=1.50870324514364;CH3CH2CH2_Pt111_energy=1.07912376831216;dEA_Pt111_11=1.92463310568732;CH3CH2CH2_Pt211_energy =1.20484266555771;dEA_Pt211_11=1.70046856953309;
CH3CHCH2_energy=0.311473412822851;dEA12=1.50952468412353;CH3CHCH2_Pt111_energy=0.550413568420995;dEA_Pt111_12=1.9546920016873;CH3CHCH2_Pt211_energy  =0.542995428754539;dEA_Pt211_12=1.66859856953312;
CH3CH2CH_energy=0.278712549431804;dEA13=3.98341484359775;CH3CH2CH_Pt111_energy=0.757655558090741;dEA_Pt111_13=3.93337768058211;CH3CH2CH_Pt211_energy  =0.632835300294667;dEA_Pt211_13=4.30076222995425;
CH2CH2CH2_energy=0.515707089955059;dEA14=2.38294479273634;CH2CH2CH2_Pt111_energy=0.786435096350582;dEA_Pt111_14=3.0365023256873;CH2CH2CH2_Pt211_energy =0.852250742049373;dEA_Pt211_14=3.19126856953314;
CH3CCH3_energy=0.189923473454927;dEA15=1.38869379474325;CH3CCH3_Pt111_energy =0.491694300805948;dEA_Pt111_15=1.8563039476873;CH3CCH3_Pt211_energy   =0.478266532249082;dEA_Pt211_15=2.1088485695331;
CH3CH2C_energy=-0.4043503533515;dEA16=2.39858018761331;CH3CH2C_Pt111_energy =-0.436108733777996;dEA_Pt111_16=3.09832372068731;CH3CH2C_Pt211_energy  =-0.124401259562124;dEA_Pt211_16=3.11840856953311;
CH2CH2CH_energy=0.577001305912265;dEA17=3.66372972092204;CH2CH2CH_Pt111_energy=0.675476751445107;dEA_Pt111_17=3.6356487076873;CH2CH2CH_Pt211_energy  =0.605133311590361;dEA_Pt211_17=3.29610856953314;
CH2CHCH2_energy=-0.312104992734898;dEA18=1.77651010731803;CH2CHCH2_Pt111_energy=0.394093164490631;dEA_Pt111_18=2.14210620768729;CH2CHCH2_Pt211_energy  =-0.0952384779303195;dEA_Pt211_18=2.06105856953311;
CH3CHCH_energy=-0.0989507940644181;dEA19=1.49422292256266;CH3CHCH_Pt111_energy =0.272302199058718;dEA_Pt111_19=1.94414158268731;CH3CHCH_Pt211_energy   =0.448783951798642;dEA_Pt211_19=2.23731856953312;
CH3CCH2_energy=-0.0413756534846037;dEA20=1.94702276931877;CH3CCH2_Pt111_energy =0.146126657081222;dEA_Pt111_20=3.37866384979252;CH3CCH2_Pt211_energy   =0.278863899026693;dEA_Pt211_20=2.58560490911198;
CH3CHC_energy=-0.174962107211496;dEA21=2.30808531527874;CH3CHC_Pt111_energy  =-0.29669008117742;dEA_Pt111_21=3.18871956679251;CH3CHC_Pt211_energy    =-0.221098772009906;dEA_Pt211_21=2.147364909112;
CH2CH2C_energy=-0.432207699465788;dEA22=0.742126737390942;CH2CH2C_Pt111_energy =0.0211202648734989;dEA_Pt111_22=1.36278013557259;CH2CH2C_Pt211_energy   =0.00796196692072826;dEA_Pt211_22=1.240494909112;
CHCH2CH_energy=-0.703823058425694;dEA23=1.16066105880748;CHCH2CH_Pt111_energy =0.0928665245758573;dEA_Pt111_23=1.6251379567925;CHCH2CH_Pt211_energy  =-0.378995639842838;dEA_Pt211_23=1.56126490911197;
CH2CHCH_energy=-0.687456640062162;dEA24=1.41060426616363;CH2CHCH_Pt111_energy =0.9394711014987;dEA_Pt111_24=1.91938362279251;CH2CHCH_Pt211_energy   =0.181401386227076;dEA_Pt211_24=1.64014490911198;
CH2CCH2_energy=-0.509615991511634;dEA25=1.96522808896361;CH2CCH2_Pt111_energy =0.0509048156469696;dEA_Pt111_25=2.27187850879252;CH2CCH2_Pt211_energy   =-0.216556152840749;dEA_Pt211_25=2.07861490911197;
CH3CCH_energy=-1.28807341763281;dEA26=1.23840128079049;CH3CCH_Pt111_energy  =-0.235986823587059;dEA_Pt111_26=1.8146878707925;CH3CCH_Pt211_energy    =-0.294060853146046;dEA_Pt211_26=2.138624909112;
CH3CC_energy=-1.3815850404672;dEA27=1.32162222475571;CH3CC_Pt111_energy  =-0.0924853990271082;dEA_Pt111_27=1.4737100627925;CH3CC_Pt211_energy    =-1.07109963503829;dEA_Pt211_27=1.51821490911197;
CH2CHC_energy=-0.844153512412358;dEA28=1.7104726242973;CH2CHC_Pt111_energy  =-0.0147922907680602;dEA_Pt111_28=2.47318051779251;CH2CHC_Pt211_energy    =-0.569421839967344;dEA_Pt211_28=2.43540490911199;
CHCHCH_energy=0.417686224554647;dEA29=1.79979054198697;CHCHCH_Pt111_energy  =-0.2293714669682;dEA_Pt111_29=2.59439466279252;CHCHCH_Pt211_energy    =-0.743960591187889;dEA_Pt211_29=2.58362490911201;
CHCH2C_energy=-0.818775579094474;dEA30=1.0054827109923;CHCH2C_Pt111_energy  =0.912030295020362;dEA_Pt111_30=1.55861931479251;CHCH2C_Pt211_energy    =-0.397935005322982;dEA_Pt211_30=1.64813490911202;
CH2CCH_energy=-1.00368719709399;dEA31=0.92356611163615;CH2CCH_Pt111_energy  =-0.359301466968201;dEA_Pt111_31=1.50019053979253;CH2CCH_Pt211_energy    =-0.458445018171485;dEA_Pt211_31=1.625734909112;
CH2CC_energy=-1.44988937145482;dEA32=1.04176243316028;CH2CC_Pt111_energy   =-0.19466120396417;dEA_Pt111_32=1.49922515579252;CH2CC_Pt211_energy     =-1.15600699592985;dEA_Pt211_32=1.75201490911202;
*/
/*
*/
  CVodeFree(&cvode_mem);

   /*Free y and abstol vectors */
    N_VDestroy_Serial(y);
    N_VDestroy_Serial(abstol);

  /* Free integrator memory */
    CVodeFree(&cvode_mem);


////////////////////////////////////////////////////////////////////////////////////////////
if (rate_723<0) {
    rate_723 = 1E-10;  //to prevent NaN in the likelihood but to prevent QUESO from sampling    there.
}

if (selectivity_723 <0.001) {
selectivity_723 = 0.001;
}

if (selectivity_723 > 0.999){
selectivity_723 = 0.999;
}


//std::cout << "Selectivity_729 = " << selectivity_723 << '\n';


double Propylene_TOF_term;
double selectivity_term;
double corr = 1;

if (rate_723 != rate_723){
Propylene_TOF_term = 1000000000;
corr = corr+1;
}
else {
Propylene_TOF_term = pow((log10(rate_723)-(log10(0.6)))/paramValues[176],2) + log(2*3.14159*paramValues[176]);// +  pow((H2r_D2-(-0.51))/paramValues[179],2);//
}

if (selectivity_723 !=selectivity_723){
selectivity_term = 1000000000;
corr = corr+1;
}
else {
selectivity_term = log(2*3.14159*paramValues[177])+pow((log(selectivity_723/(1-selectivity_723))-3.47)/paramValues[177],2); //+ +pow((log(selectivity_823/(1-selectivity_823))-1.33)/paramValues[529],2)
}

double result =  Propylene_TOF_term + selectivity_term;
//  if (alph > 0.2 | beta > 0.2 | gamm > 0.2 | alph < -0.2 | beta < -0.2 | gamm < -0.2 ) // | delta < 0.1
 // {
 //    result = 1000000000;
 // };

if (result != result){
result = 200.;
}
if (result > 2000) {
result = 200.;
}
#ifdef QUESO_EXPECTS_LN_LIKELIHOOD_INSTEAD_OF_MINUS_2_LN
  return -0.5*result;
#else

#endif
}
