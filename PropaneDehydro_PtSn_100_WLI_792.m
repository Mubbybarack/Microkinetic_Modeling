function F = PropaneDehydro_PtSn_100_WLI_792(t,y)

PCH3CH2CH3 = 0.29 ;
PH2 = 0.09 ;
PCH3CHCH2 = 0 ;
PCH3CCH = 0 ;
PCH3CH3 = 0 ;
PCH2CH2 = 0 ;
PCHCH = 0 ;
PCH4 = 0 ;
T = 792.0;


% Initializing the forward rate variable
f = zeros(138,1);

% Adsorption energy (dG) for all 46 intermediates and 130 TST

CH3CH2CH3_energyy	=	2.6459 ;
CH3CHCH3_energyy	=	3.3230 ;
CH3CH2CH2_energyy	=	3.2242 ;
CH3CHCH2_energyy	=	3.3979 ;
CH3CH2CH_energyy	=	4.1839 ;
CH2CH2CH2_energyy	=	4.0149 ;
CH3CCH3_energyy		=	3.8911 ;
CH3CH2C_energyy		=	4.7901 ;
CH2CH2CH_energyy	=	4.7750 ;
CH2CHCH2_energyy	=	3.9891 ;
CH3CHCH_energyy		=	3.9235 ;
CH3CCH2_energyy		=	3.8842 ;
CH3CHC_energyy		=	4.7367 ;
CH2CH2C_energyy		=	5.5880 ;
CHCH2CH_energyy		=	5.3348 ;
CH2CHCH_energyy		=	4.5049 ;
CH2CCH2_energyy		=	4.7692 ;
CH3CCH_energyy		=	4.4163 ;
CH3CC_energyy		=	4.6685 ;
CH2CHC_energyy		=	5.6049 ;
CHCHCH_energyy		=	4.7918 ;
CHCH2C_energyy		=	6.5780 ;
CH2CCH_energyy		=	4.8409 ;
CH2CC_energyy		=	5.5176 ;
CHCHC_energyy		=	6.0082 ;
CCH2C_energyy		=	7.3713 ;
CHCCH_energyy		=	5.5169 ;
CCHC_energyy		=	6.9940 ;
CHCC_energyy		=	6.2271 ;
CCC_energyy			=	7.2278 ;
CH3CH3_energyy		=	1.4447 ;
CH3CH2_energyy		=	2.1407 ;
CH3CH_energyy		=	3.1495 ;
CH3C_energyy		=	3.7290 ;
CH2CH2_energyy		=	2.3564 ;
CH2CH_energyy		=	2.8916 ;
CH2C_energyy		=	3.9329 ;
CHCH_energyy		=	3.5453 ;
CHC_energyy			=	3.8692 ;
CC_energyy			=	4.8926 ;
CH4_energyy			=	0.0747 ;
CH3_energyy			=	0.9176 ;
CH2_energyy			=	2.0700 ;
CH_energyy			=	2.7256 ;
C_energyy			=	3.3245 ;
H_energyy			=	-0.2477 ;

ddEA9	=	3.9772 ;
ddEA10	=	3.9882 ;
ddEA11	=	4.7176 ;
ddEA12	=	4.6363 ;
ddEA13	=	4.8293 ;
ddEA14	=	5.6804 ;
ddEA15	=	4.8852 ;
ddEA16	=	5.6120 ;
ddEA17	=	5.9038 ;
ddEA18	=	4.6139 ;
ddEA19	=	4.5918 ;
ddEA20	=	6.8630 ;
ddEA21	=	6.3507 ;
ddEA22	=	5.8801 ;
ddEA23	=	4.7771 ;
ddEA24	=	5.4685 ;
ddEA25	=	5.5874 ;
ddEA26	=	5.1594 ;
ddEA27	=	5.5812 ;
ddEA28	=	6.3935 ;
ddEA29	=	5.9120 ;
ddEA30	=	5.1159 ;
ddEA31	=	4.5127 ;
ddEA32	=	4.1914 ;
ddEA33	=	6.3566 ;
ddEA34	=	5.3338 ;
ddEA35	=	6.9449 ;
ddEA36	=	6.5660 ;
ddEA37	=	6.5069 ;
ddEA38	=	5.6553 ;
ddEA39	=	7.1749 ;
ddEA40	=	6.4121 ;
ddEA41	=	6.1811 ;
ddEA42	=	6.0905 ;
ddEA43	=	5.7353 ;
ddEA44	=	6.1621 ;
ddEA45	=	5.1064 ;
ddEA46	=	6.1379 ;
ddEA47	=	6.8120 ;	
ddEA48	=	7.4115 ;	
ddEA49	=	5.4108 ;	
ddEA50	=	6.0757 ;
ddEA51	=	5.7974 ;	
ddEA52	=	6.6723 ;	
ddEA53	=	6.8880 ;	
ddEA54	=	5.4349 ;	
ddEA55	=	6.0324 ;	
ddEA56	=	7.1353 ;	
ddEA57	=	8.0324 ;	
ddEA58	=	6.2423 ;	
ddEA59	=	6.1865 ;	
ddEA60	=	6.7815 ;	
ddEA61	=	7.1694 ;	
ddEA62	=	7.0484 ;	
ddEA63	=	7.1757 ;	
ddEA64	=	7.4770 ;	
ddEA65	=	7.3061 ;	
ddEA66	=	6.5677 ;	
ddEA67	=	6.7515 ;	
ddEA68	=	6.8114 ;	
ddEA69	=	5.9908 ;	
ddEA70	=	6.9681 ;
ddEA71	=	5.5593 ;	
ddEA72	=	7.1021 ;	
ddEA73	=	5.7217 ;	
ddEA74	=	7.5066 ;	
ddEA75	=	6.4659 ;	
ddEA76	=	5.9048 ;
ddEA77	=	5.9747 ;
ddEA78	=	7.2513 ;	
ddEA79	=	8.2249 ;	
ddEA80	=	6.1992 ;
ddEA81	=	7.5071 ;	
ddEA82	=	7.0936 ;	
ddEA83	=	7.1597 ;	
ddEA84	=	6.7606 ;	
ddEA85	=	8.0578 ;	
ddEA86	=	7.5292 ;	
ddEA87	=	7.8906 ;	
ddEA88	=	8.2521 ;	
ddEA89	=	7.3160 ;	
ddEA90	=	6.7153 ;	
ddEA91	=	7.2585 ;	
ddEA92	=	7.2980 ;	
ddEA93	=	7.9280 ;	
ddEA94	=	6.2394 ;	
ddEA95	=	6.0144 ;	
ddEA96	=	8.0423 ;	
ddEA97	=	9.0031 ;	
ddEA98	=	7.0535 ;	
ddEA99	=	8.6019 ;	
ddEA100	=	8.4709 ;	
ddEA101	=	7.9230 ;	
ddEA102	=	7.9632 ;	
ddEA103	=	9.2791 ;	
ddEA104	=	8.8546 ;	
ddEA105	=	8.6453 ;	
ddEA106	=	6.9084 ;	
ddEA107	=	8.8110 ;	
ddEA108	=	9.1639 ;	
ddEA109	=	9.2121 ;	
ddEA110	=	8.6216 ;	
ddEA111	=	8.2291 ;	
ddEA112	=	10.3712 ;	
ddEA113	=	1.7024 ;	
ddEA114	=	2.4368 ;	
ddEA115	=	3.4330 ;	
ddEA116	=	3.8433 ;	
ddEA117	=	4.4597 ;	
ddEA118	=	2.8748 ;	
ddEA119	=	4.7502 ;	
ddEA120	=	3.5094 ;	
ddEA121	=	3.7792 ;	
ddEA122	=	5.4856 ;	
ddEA123	=	4.7682 ;	
ddEA124	=	4.5106 ;	
ddEA125	=	5.3634 ;	
ddEA126	=	4.6871 ;	
ddEA127	=	5.4143 ;	
ddEA128	=	3.5823 ;	
ddEA129	=	6.0073 ;	
ddEA130	=	4.2890 ;	
ddEA131	=	4.3950 ;
ddEA132	=	6.4878 ;	
ddEA133	=	5.0653 ;
ddEA134	=	7.2965 ;	
ddEA135	=	5.1389 ;
ddEA136	=	7.5007 ;	
ddEA137	=	5.5255 ;	
ddEA138	=	8.3086 ;	

% Gas phase energy
CH3CH2CH3_gp 	= 1.7536 ;
CH3CHCH2_gp 	= 2.5892 ;
CH3CCH_gp 		= 4.0418 ;
CH3CH3_gp 		= 0.7522 ;
CH2CH2_gp 		= 1.7580 ;
CHCH_gp 		= 3.8069 ;
CH4_gp 			= -0.3694 ;
H2_gp 			= -0.7988 ;

% Collision theory derived forward rate 
f(1) = 1.127655298E+08 ;
f(2) = 1.154403922E+08 ;
f(3) = 5.295179191E+08 ;
f(4) = 1.183093444E+08 ;
f(5) = 1.365617097E+08 ;
f(6) = 1.413934289E+08 ;
f(7) = 1.467490064E+08 ;
f(8) = 1.869792775E+08 ;

% P is ads and TST energy uncertainity and it's equal zero since no uncertainity is added
pp = zeros(176); 

% Lateral interaction energy (dG) for H, CHCH and C with all 46 intermediates



% Calculating total energy (dG) with lateral interaction for all 16 intermediates
CH3CH2CH3 	= pp(0+1) 	+ CH3CH2CH3_energyy ; %+ CH3CH2CH3_H_sloppe*y(2)+CH3CH2CH3_CHCH_sloppe*y(40)*4 + CH3CH2CH3_C_sloppe*y(47)*4;
CH3CHCH3 	= pp(1+1)	+ CH3CHCH3_energyy 	; %+ CH3CHCH3_H_sloppe*y(2)+CH3CHCH3_CHCH_sloppe*y(40)*4 + CH3CHCH3_C_sloppe*y(47)*4;
CH3CH2CH2 	= pp(2+1)	+ CH3CH2CH2_energyy ; %+ CH3CH2CH2_H_sloppe*y(2)+CH3CH2CH2_CHCH_sloppe*y(40)*4 + CH3CH2CH2_C_sloppe*y(47)*4;
CH3CHCH2 	= pp(3+1)	+ CH3CHCH2_energyy	; %+ CH3CHCH2_H_sloppe*y(2)+CH3CHCH2_CHCH_sloppe*y(40)*4 + CH3CHCH2_C_sloppe*y(47)*4;
CH3CH2CH 	= pp(4+1)	+ CH3CH2CH_energyy 	; %+ CH3CH2CH_H_sloppe*y(2)+CH3CH2CH_CHCH_sloppe*y(40)*4 + CH3CH2CH_C_sloppe*y(47)*4;
CH2CH2CH2	= pp(5+1)	+ CH2CH2CH2_energyy ; %+ CH2CH2CH2_H_sloppe*y(2)+CH2CH2CH2_CHCH_sloppe*y(40)*4 + CH2CH2CH2_C_sloppe*y(47)*4;
CH3CCH3		= pp(6+1)	+ CH3CCH3_energyy 	; %+ CH3CCH3_H_sloppe*y(2)+CH3CCH3_CHCH_sloppe*y(40)*4+ CH3CCH3_C_sloppe*y(47)*4;

CH3CH2C		= pp(7+1)	+ CH3CH2C_energyy	; %+ CH3CH2C_H_sloppe*y(2)+CH3CH2C_CHCH_sloppe*y(40)*4 + CH3CH2C_C_sloppe*y(47)*4;
CH2CH2CH	= pp(8+1)	+ CH2CH2CH_energyy 	; %+ CH2CH2CH_H_sloppe*y(2)+CH2CH2CH_CHCH_sloppe*y(40)*4 + CH2CH2CH_C_sloppe*y(47)*4;
CH2CHCH2	= pp(9+1)	+ CH2CHCH2_energyy 	; %+ CH2CHCH2_H_sloppe*y(2)+CH2CHCH2_CHCH_sloppe*y(40)*4 + CH2CHCH2_C_sloppe*y(47)*4;
CH3CHCH 	= pp(10+1)	+ CH3CHCH_energyy 	; %+ CH3CHCH_H_sloppe*y(2)+CH3CHCH_CHCH_sloppe*y(40)*4 + CH3CHCH_C_sloppe*y(47)*4;
CH3CCH2 	= pp(11+1)	+ CH3CCH2_energyy 	; %+CH3CCH2_H_sloppe*y(2)+CH3CCH2_CHCH_sloppe*y(40)*4 + CH3CCH2_C_sloppe*y(47)*4;
CH3CHC  	= pp(12+1)	+ CH3CHC_energyy 	; %+ CH3CHC_H_sloppe*y(2)+CH3CHC_CHCH_sloppe*y(40)*4 + CH3CHC_C_sloppe*y(47)*4;
CH2CH2C 	= pp(13+1)	+ CH2CH2C_energyy 	; %+CH2CH2C_H_sloppe*y(2)+CH2CH2C_CHCH_sloppe*y(40)*4 + CH2CH2C_C_sloppe*y(47)*4;
CHCH2CH 	= pp(14+1)	+ CHCH2CH_energyy 	; %+CHCH2CH_H_sloppe*y(2)+CHCH2CH_CHCH_sloppe*y(40)*4 + CHCH2CH_C_sloppe*y(47)*4;
CH2CHCH 	= pp(15+1)	+ CH2CHCH_energyy 	; %+  CH2CHCH_H_sloppe*y(2)+CH2CHCH_CHCH_sloppe*y(40)*4 + CH2CHCH_C_sloppe*y(47)*4;
CH2CCH2 	= pp(16+1)	+ CH2CCH2_energyy 	; %+ CH2CCH2_H_sloppe*y(2)+CH2CCH2_CHCH_sloppe*y(40)*4 + CH2CCH2_C_sloppe*y(47)*4;
CH3CCH  	= pp(17+1)	+ CH3CCH_energyy 	; %+ CH3CCH_H_sloppe*y(2)+CH3CCH_CHCH_sloppe*y(40)*4 + CH3CCH_C_sloppe*y(47)*4;

CH3CC  		= pp(18+1)  + CH3CC_energyy 	; %+ CH3CC_H_sloppe*y(2)+CH3CC_CHCH_sloppe*y(40)*4 + CH3CC_C_sloppe*y(47)*4;
CH2CHC 		= pp(19+1)  + CH2CHC_energyy 	; %+ CH2CHC_H_sloppe*y(2)+CH2CHC_CHCH_sloppe*y(40)*4 + CH2CHC_C_sloppe*y(47)*4;
CHCHCH 		= pp(20+1)  + CHCHCH_energyy 	; %+ CHCHCH_H_sloppe*y(2)+CHCHCH_CHCH_sloppe*y(40)*4 + CHCHCH_C_sloppe*y(47)*4;
CHCH2C 		= pp(21+1)  + CHCH2C_energyy 	; %+ CHCH2C_H_sloppe*y(2)+CHCH2C_CHCH_sloppe*y(40)*4 + CHCH2C_C_sloppe*y(47)*4;
CH2CCH 		= pp(22+1)  + CH2CCH_energyy 	; %+ CH2CCH_H_sloppe*y(2)+CH2CCH_CHCH_sloppe*y(40)*4 + CH2CCH_C_sloppe*y(47)*4;
CH2CC  		= pp(23+1)  + CH2CC_energyy 	; %+ CH2CC_H_sloppe*y(2)+CH2CC_CHCH_sloppe*y(40)*4 + CH2CC_C_sloppe*y(47)*4;
CHCHC  		= pp(24+1)  + CHCHC_energyy 	; %+ CHCHC_H_sloppe*y(2)+CHCHC_CHCH_sloppe*y(40)*4 + CHCHC_C_sloppe*y(47)*4;
CCH2C  		= pp(25+1)  + CCH2C_energyy 	; %+ CCH2C_H_sloppe*y(2)+CCH2C_CHCH_sloppe*y(40)*4 + CCH2C_C_sloppe*y(47)*4;
CHCCH  		= pp(26+1)  + CHCCH_energyy 	; %+ CHCCH_H_sloppe*y(2)+CHCCH_CHCH_sloppe*y(40)*4 + CHCCH_C_sloppe*y(47)*4;
CCHC   		= pp(27+1)  + CCHC_energyy 		; %+ CCHC_H_sloppe*y(2)+CCHC_CHCH_sloppe*y(40)*4 + CCHC_C_sloppe*y(47)*4;
CHCC   		= pp(28+1)  + CHCC_energyy 		; %+ CHCC_H_sloppe*y(2)+CHCC_CHCH_sloppe*y(40)*4 + CHCC_C_sloppe*y(47)*4;
CCC    		= pp(29+1)  + CCC_energyy 		; %+ CCC_H_sloppe*y(2)+CCC_CHCH_sloppe*y(40)*4 + CCC_C_sloppe*y(47)*4;

CH3CH3		= pp(30+1)  + CH3CH3_energyy 	; %+ CH3CH3_H_sloppe*y(2)+CH3CH3_CHCH_sloppe*y(40)*4 + CH3CH3_C_sloppe*y(47)*4;
CH3CH2		= pp(31+1)  + CH3CH2_energyy 	; %+ CH3CH2_H_sloppe*y(2)+CH3CH2_CHCH_sloppe*y(40)*4 + CH3CH2_C_sloppe*y(47)*4;
CH3CH 		= pp(32+1)  + CH3CH_energyy 	; %+ CH3CH_H_sloppe*y(2)+CH3CH_CHCH_sloppe*y(40)*4 + CH3CH_C_sloppe*y(47)*4;
CH3C  		= pp(33+1)  + CH3C_energyy 		; %+ CH3C_H_sloppe*y(2)+CH3C_CHCH_sloppe*y(40)*4 + CH3C_C_sloppe*y(47);
CH2CH2		= pp(34+1)  + CH2CH2_energyy 	; %+ CH2CH2_H_sloppe*y(2)+CH2CH2_CHCH_sloppe*y(40)*4 + CH2CH2_C_sloppe*y(47)*4;
CH2CH 		= pp(35+1)  + CH2CH_energyy 	; %+ CH2CH_H_sloppe*y(2)+CH2CH_CHCH_sloppe*y(40)*4 + CH2CH_C_sloppe*y(47)*4;
CH2C  		= pp(36+1)  + CH2C_energyy 		; %+ CH2C_H_sloppe*y(2)+CH2C_CHCH_sloppe*y(40)*4 + CH2C_C_sloppe*y(47)*4;
CHCH  		= pp(37+1)  + CHCH_energyy 		; %+ CHCH_H_sloppe*y(2)+CHCH_CHCH_sloppe*y(40)*4 + CHCH_C_sloppe*y(47)*4;
CHC   		= pp(38+1)  + CHC_energyy 		; %+ CHC_H_sloppe*y(2)+CHC_CHCH_sloppe*y(40)*4 + CHC_C_sloppe*y(47)*4;
CC    		= pp(39+1)  + CC_energyy 		; %+ CC_H_sloppe*y(2)+CC_CHCH_sloppe*y(40)*4 + CC_C_sloppe*y(47)*4;

CH4   		= pp(40+1)	+ CH4_energyy 		; %+ CH4_H_sloppe*y(2)+CH4_CHCH_sloppe*y(40)*4 + CH4_C_sloppe*y(47)*4;
CH3   		= pp(41+1)	+ CH3_energyy 		; %+ CH3_H_sloppe*y(2)+CH3_CHCH_sloppe*y(40)*4 + CH3_C_sloppe*y(47)*4;
CH2   		= pp(42+1)	+ CH2_energyy 		; %+ CH2_H_sloppe*y(2)+CH2_CHCH_sloppe*y(40)*4 + CH2_C_sloppe*y(47)*4;
CH    		= pp(43+1)	+ CH_energyy 		; %+ CH_H_sloppe*y(2)+CH_CHCH_sloppe*y(40)*4 + CH_C_sloppe*y(47)*4;
C     		= pp(44+1)	+ C_energyy 		; %+ C_H_sloppe*y(2)+C_CHCH_sloppe*y(40)*4 + C_C_sloppe*y(47)*4;
H     		= pp(45+1)	+ H_energyy 		; %+ H_H_sloppe*y(2)+H_CHCH_sloppe*y(40)*4 + H_C_sloppe*y(47)*4;

% calculating total energy (dG) without lateral interaction to be used in calculation of TST lateral interaction

CH3CH2CH3_0 	= pp(0+1) 	+ CH3CH2CH3_energyy  ;
CH3CHCH3_0 		= pp(1+1)	+ CH3CHCH3_energyy ;
CH3CH2CH2_0 	= pp(2+1)	+ CH3CH2CH2_energyy ;
CH3CHCH2_0 		= pp(3+1)	+ CH3CHCH2_energyy ;
CH3CH2CH_0 		= pp(4+1)	+ CH3CH2CH_energyy ;
CH2CH2CH2_0 	= pp(5+1)	+ CH2CH2CH2_energyy ;
CH3CCH3_0 		= pp(6+1)	+ CH3CCH3_energyy ;

CH3CH2C_0		= pp(7+1)	+ CH3CH2C_energyy;
CH2CH2CH_0		= pp(8+1)	+ CH2CH2CH_energyy;
CH2CHCH2_0		= pp(9+1)	+ CH2CHCH2_energyy;
CH3CHCH_0 		= pp(10+1)	+ CH3CHCH_energyy;
CH3CCH2_0 		= pp(11+1)	+ CH3CCH2_energyy;
CH3CHC_0  		= pp(12+1)	+ CH3CHC_energyy ;
CH2CH2C_0 		= pp(13+1)	+ CH2CH2C_energyy ;
CHCH2CH_0 		= pp(14+1)	+ CHCH2CH_energyy ;
CH2CHCH_0 		= pp(15+1)	+ CH2CHCH_energyy ;
CH2CCH2_0 		= pp(16+1)	+ CH2CCH2_energyy ;
CH3CCH_0  		= pp(17+1)	+ CH3CCH_energyy ;

CH3CC_0  		= pp(18+1)	+ CH3CC_energyy ;
CH2CHC_0 		= pp(19+1)	+ CH2CHC_energyy ;
CHCHCH_0 		= pp(20+1)	+ CHCHCH_energyy ;
CHCH2C_0 		= pp(21+1)	+ CHCH2C_energyy ;
CH2CCH_0 		= pp(22+1)	+ CH2CCH_energyy ;
CH2CC_0  		= pp(23+1)	+ CH2CC_energyy ;
CHCHC_0  		= pp(24+1)	+ CHCHC_energyy ;
CCH2C_0  		= pp(25+1)	+ CCH2C_energyy ;
CHCCH_0  		= pp(26+1)	+ CHCCH_energyy ;
CCHC_0   		= pp(27+1)	+ CCHC_energyy ;
CHCC_0   		= pp(28+1)	+ CHCC_energyy ;
CCC_0    		= pp(29+1)	+ CCC_energyy ;

CH3CH3_0 		= pp(30+1)	+ CH3CH3_energyy ;
CH3CH2_0 		= pp(31+1)	+ CH3CH2_energyy ;
CH3CH_0  		= pp(32+1)	+ CH3CH_energyy ;
CH3C_0   		= pp(33+1)	+ CH3C_energyy ;
CH2CH2_0 		= pp(34+1)	+ CH2CH2_energyy ;
CH2CH_0  		= pp(35+1)	+ CH2CH_energyy ;
CH2C_0  		= pp(36+1)	+ CH2C_energyy ;
CHCH_0  		= pp(37+1)	+ CHCH_energyy ;
CHC_0   		= pp(38+1)	+ CHC_energyy ;
CC_0    		= pp(39+1)	+ CC_energyy ;

CH4_0 			= pp(40+1)	+ CH4_energyy ;
CH3_0 			= pp(41+1)	+ CH3_energyy ;
CH2_0 			= pp(42+1)	+ CH2_energyy ;
CH_0  			= pp(43+1)	+ CH_energyy ;
C_0   			= pp(44+1)	+ C_energyy ;
H_0   			= pp(45+1)	+ H_energyy ;


% Calculating TST energy (dG) with lateral interactions for all 130 rxns

CH3CH2CH3_to_CH3CHCH3_H_H_sloppe 	= 0.5*((H+CH3CHCH3 - CH3CH2CH3)-(H_0     +CH3CHCH3_0      - CH3CH2CH3_0     ));
CH3CH2CH3_to_CH3CH2CH2_H_H_sloppe	= 0.5*((H+CH3CH2CH2 - CH3CH2CH3)-(H_0     +CH3CH2CH2_0      - CH3CH2CH3_0     ));
CH3CHCH3_to_CH3CHCH2_H_H_sloppe  	= 0.5*((H+CH3CHCH2 - CH3CHCH3)-(H_0     +CH3CHCH2_0      - CH3CHCH3_0     ));
CH3CH2CH2_to_CH3CHCH2_H_H_sloppe 	= 0.5*((H+CH3CHCH2 - CH3CH2CH2)-(H_0     +CH3CHCH2_0      - CH3CH2CH2_0     ));
CH3CH2CH3_to_CH3_CH3CH2_H_sloppe 	= 0.5*((CH3+CH3CH2 - CH3CH2CH3)-(CH3_0      + CH3CH2_0      - CH3CH2CH3_0     ));
CH3CHCH3_to_CH3_CH3CH_H_sloppe   	= 0.5*((CH3+CH3CH - CH3CHCH3) -(CH3_0      + CH3CH_0      - CH3CHCH3_0     ));
CH3CHCH3_to_CH3CCH3_H_H_sloppe     	= 0.5*((H + CH3CCH3 - CH3CHCH3) - (H_0      + CH3CCH3_0      - CH3CHCH3_0     ));
CH3CH2CH2_to_CH3CH2_CH2_H_sloppe   	= 0.5*((CH3CH2 + CH2 - CH3CH2CH2) -(CH3CH2_0      + CH2_0      - CH3CH2CH2_0     ));
CH3CH2CH2_to_CH3_CH2CH2_H_sloppe   	= 0.5*((CH3+CH2CH2 - CH3CH2CH2) - (CH3_0      + CH2CH2_0      - CH3CH2CH2_0     ));
CH3CH2CH2_to_CH2CH2CH2_H_H_sloppe  	= 0.5*((H + CH2CH2CH2 - CH3CH2CH2) - (H_0     +CH2CH2CH2_0      - CH3CH2CH2_0     ));
CH3CH2CH2_to_CH3CH2CH_H_H_sloppe   	= 0.5*((H + CH3CH2CH - CH3CH2CH2) - (H_0      + CH3CH2CH_0      - CH3CH2CH2_0     ));

CH3CH2CH_to_CH3CH2_CH_H_sloppe   	= 0.5*((CH3CH2 + CH - CH3CH2CH) - (CH3CH2_0      + CH_0      - CH3CH2CH_0     ));
CH3CH2CH_to_CH3_CH2CH_H_sloppe   	= 0.5*((CH2CH + CH3 - CH3CH2CH) - (CH2CH_0      + CH3_0      - CH3CH2CH_0     ));
CH3CH2CH_to_CH3CH2C_H_H_sloppe   	= 0.5*((H+CH3CH2C - CH3CH2CH) - (H_0      + CH3CH2C_0      - CH3CH2CH_0     ));
CH3CH2CH_to_CH3CHCH_H_H_sloppe   	= 0.5*((H+CH3CHCH - CH3CH2CH) - (H_0      + CH3CHCH_0      - CH3CH2CH_0     ));
CH3CH2CH_to_CH2CH2CH_H_H_sloppe  	= 0.5*((H+CH2CH2CH - CH3CH2CH) - (H_0      + CH2CH2CH_0      - CH3CH2CH_0     ));
CH2CH2CH2_to_CH2_CH2CH2_H_sloppe 	= 0.5*((CH2+CH2CH2 - CH2CH2CH2) - (CH2_0      + CH2CH2_0      - CH2CH2CH2_0     ));
CH2CH2CH2_to_CH2CH2CH_H_H_sloppe 	= 0.5*((CH2CH2CH + H - CH2CH2CH2) -(CH2CH2CH_0      + H_0      - CH2CH2CH2_0     ));
CH2CH2CH2_to_CH2CHCH2_H_H_sloppe 	= 0.5*((H + CH2CHCH2 - CH2CH2CH2) - (CH2CHCH2_0      +H_0      - CH2CH2CH2_0     ));
CH3CHCH2_to_CH3CH_CH2_H_sloppe 		= 0.5*((CH3CH+CH2 - CH3CHCH2) - (CH3CH_0      + CH2_0      - CH3CHCH2_0     ));
CH3CHCH2_to_CH3_CHCH2_H_sloppe 		= 0.5*((CH3 + CH2CH - CH3CHCH2) - (CH3_0      + CH2CH_0      - CH3CHCH2_0     ));
CH3CHCH2_to_CH3CCH2_H_H_sloppe 		= 0.5*((CH3CCH2 + H - CH3CHCH2) - (CH3CCH2_0      + H_0      - CH3CHCH2_0     ));
CH3CHCH2_to_CH3CHCH_H_H_sloppe 		= 0.5*((CH3CHCH + H - CH3CHCH2) - (CH3CHCH_0      + H_0      - CH3CHCH2_0     ));
CH3CHCH2_to_CH2CHCH2_H_H_sloppe 	= 0.5*((CH2CHCH2 + H - CH3CHCH2) - (CH2CHCH2_0      + H_0      - CH3CHCH2_0     ));
CH3CCH3_to_CH3_CCH3_H_sloppe 		= 0.5*((CH3 + CH3C - CH3CCH3)-(CH3_0      +CH3C_0      - CH3CCH3_0     ));
CH3CCH3_to_CH3CCH2_H_H_sloppe 		= 0.5*((H + CH3CCH2 - CH3CCH3)-(H_0      + CH3CCH2_0      - CH3CCH3_0     ));

CH3CH2C_to_CH3_CH2C_H_sloppe 		= 0.5*((CH3 + CH2C - CH3CH2C) -(CH3_0      + CH2C_0      - CH3CH2C_0     ));
CH3CH2C_to_CH3CH2_C_H_sloppe 		= 0.5*((CH3CH2+C - CH3CH2C) - (CH3CH2_0      + C_0      - CH3CH2C_0     ));
CH3CH2C_to_CH2CH2C_H_H_sloppe  		= 0.5*((CH2CH2C +H - CH3CH2C) - (CH2CH2C_0      + H_0      - CH3CH2C_0     ));
CH3CH2C_to_CH3CHC_H_H_sloppe 		= 0.5*((CH3CHC +H - CH3CH2C) -(CH3CHC_0      +H_0      - CH3CH2C_0     ));
CH2CH2CH_to_CH2_CH2CH_H_sloppe 		= 0.5*((CH2 + CH2CH - CH2CH2CH) -(CH2_0      + CH2CH_0      - CH2CH2CH_0     ));
CH2CH2CH_to_CH2CH2_CH_H_sloppe 		= 0.5*((CH2CH2 + CH - CH2CH2CH) -(CH2CH2_0      + CH_0      - CH2CH2CH_0     ));
CH2CH2CH_to_CH2CH2C_H_H_sloppe 		= 0.5*((CH2CH2C +H - CH2CH2CH) - (CH2CH2C_0      + H_0      - CH2CH2CH_0     ));
CH2CH2CH_to_CH2CHCH_H_H_sloppe 		= 0.5*((CH2CHCH +H - CH2CH2CH) - (CH2CHCH_0      + H_0      - CH2CH2CH_0     ));
CH2CH2CH_to_CHCH2CH_H_H_sloppe 		= 0.5*((CHCH2CH +H - CH2CH2CH) - (CHCH2CH_0      + H_0      - CH2CH2CH_0     ));
CH2CHCH2_to_CH2_CHCH2_H_sloppe 		= 0.5*((CH2 +CH2CH - CH2CHCH2) - (CH2_0      + CH2CH_0      - CH2CHCH2_0     ));
CH2CHCH2_to_CH2CHCH_H_H_sloppe 		= 0.5*((CH2CHCH +H - CH2CHCH2) - (CH2CHCH_0      + CH2CH_0      - CH2CHCH2_0     ));
CH2CHCH2_to_CH2CCH2_H_H_sloppe 		= 0.5*((CH2CCH2 +H - CH2CHCH2) - (CH2CCH2_0      + H_0      - CH2CHCH2_0     ));
CH3CHCH_to_CH3_CHCH_H_sloppe 		= 0.5*((CH3 + CHCH - CH3CHCH) - (CH3_0      + CHCH_0      - CH3CHCH_0     ));
CH3CHCH_to_CH3CH_CH_H_sloppe 		= 0.5*((CH3CH + CH - CH3CHCH) - (CH3CH_0      + CH_0      - CH3CHCH_0     ));
CH3CHCH_to_CH3CHC_H_H_sloppe 		= 0.5*((CH3CHC +H - CH3CHCH) - (CH3CHC_0      + H_0      - CH3CHCH_0     ));
CH3CHCH_to_CH3CCH_H_H_sloppe 		= 0.5*((CH3CCH +H - CH3CHCH) - (CH3CCH_0      + H_0      - CH3CHCH_0     ));
CH3CHCH_to_CH2CHCH_H_H_sloppe 		= 0.5*((CH2CHCH +H - CH3CHCH) - (CH2CHCH_0      +H_0      - CH3CHCH_0     ));
CH3CCH2_to_CH3_CCH2_H_sloppe 		= 0.5*((CH3 + CH2C - CH3CCH2) - (CH3_0      + CH2C_0      - CH3CCH2_0     ));
CH3CCH2_to_CH3C_CH2_H_sloppe 		= 0.5*((CH3C + CH2 - CH3CCH2) - (CH3C_0      + CH2_0      - CH3CCH2_0     ));
CH3CCH2_to_CH2CCH2_H_H_sloppe 		= 0.5*((CH2CCH2 + H - CH3CCH2) - (CH2CCH2_0      + H_0      - CH3CCH2_0     ));
CH3CCH2_to_CH3CCH_H_H_sloppe 		= 0.5*((CH3CCH +H - CH3CCH2) - (CH3CCH_0      + H_0      - CH3CCH2_0      ));

CH3CHC_to_CH3_CHC_H_sloppe 			= 0.5*((CH3 + CHC - CH3CHC) -(CH3_0      + CHC_0      - CH3CHC_0     ));
CH3CHC_to_CH3CH_C_H_sloppe 			= 0.5*((CH3CH + C - CH3CHC) -(CH3CH_0      + C_0      - CH3CHC_0     ));
CH3CHC_to_CH3CC_H_H_sloppe 			= 0.5*((CH3CC +H - CH3CHC) - (CH3CC_0      + H_0      - CH3CHC_0     ));
CH3CHC_to_CH2CHC_H_H_sloppe 		= 0.5*((CH2CHC +H - CH3CHC) - (CH2CHC_0      +H_0      - CH3CHC_0     ));
CH2CH2C_to_CH2CH2_C_H_sloppe 		= 0.5*((CH2CH2 + C - CH2CH2C) -(CH2CH2_0      +C_0      - CH2CH2C_0     ));
CH2CH2C_to_CH2_CH2C_H_sloppe 		= 0.5*((CH2 + CH2C - CH2CH2C) -(CH2_0      + CH2C_0      - CH2CH2C_0     ));
CH2CH2C_to_CH2CHC_H_H_sloppe 		= 0.5*((CH2CHC +H - CH2CH2C) - (CH2CHC_0      + H_0      - CH2CH2C_0     ));
CH2CH2C_to_CHCH2C_H_H_sloppe 		= 0.5*((CHCH2C +H - CH2CH2C) - (CHCH2C_0      + H_0      - CH2CH2C_0     ));
CHCH2CH_to_CHCH2_CH_H_sloppe 		= 0.5*((CH2CH + CH - CHCH2CH) -(CH2CH_0      + CH_0      - CHCH2CH_0     ));
CHCH2CH_to_CHCH2C_H_H_sloppe 		= 0.5*((CHCH2C +H - CHCH2CH) - (CHCH2C_0      + H_0      - CHCH2CH_0     ));
CHCH2CH_to_CHCHCH_H_H_sloppe 		= 0.5*((CHCHCH + H - CHCH2CH) - (CHCHCH_0      +H_0      - CHCH2CH_0     ));
CH2CHCH_to_CH2_CHCH_H_sloppe 		= 0.5*((CH2 +CHCH - CH2CHCH) - (CH2_0      + CHCHC_0      - CH2CHCH_0     ));
CH2CHCH_to_CH2CH_CH_H_sloppe 		= 0.5*((CH2CH + CH - CH2CHCH) - (CH2CH_0      + CH_0      - CH2CHCH_0     ));
CH2CHCH_to_CH2CHC_H_H_sloppe 		= 0.5*((CH2CHC +H - CH2CHCH)-(CH2CHC_0      + H_0      - CH2CHCH_0     ));
CH2CHCH_to_CH2CCH_H_H_sloppe 		= 0.5*((CH2CCH +H - CH2CHCH) -(CH2CCH_0      +H_0      - CH2CHCH_0     ));
CH2CHCH_to_CHCHCH_H_H_sloppe 		= 0.5*((CHCHCH +H - CH2CHCH) -(CHCHCH_0      + H_0      - CH2CHCH_0     ));
CH2CCH2_to_CH2C_CH2_H_sloppe 		= 0.5*((CH2C + CH2 - CH2CCH2) - (CH2C_0      + CH2_0      - CH2CCH2_0     ));
CH2CCH2_to_CH2CCH_H_H_sloppe 		= 0.5*((CH2CCH + H - CH2CCH2) - (CH2CCH_0      + H_0      - CH2CCH2_0     ));
CH3CCH_to_CH3C_CH_H_sloppe 			= 0.5*((CH3C + CH - CH3CCH) - (CH3C_0      + CH_0      - CH3CCH_0     ));
CH3CCH_to_CH3_CHC_H_sloppe 			= 0.5*((CH3 + CHC - CH3CCH) - (CH3_0      + CHC_0      - CH3CCH_0     ));
CH3CCH_to_CH3CC_H_H_sloppe 			= 0.5*((CH3CC +H - CH3CCH) - ( CH3CC_0      + H_0      - CH3CCH_0     ));
CH3CCH_to_CH2CCH_H_H_sloppe 		= 0.5*((CH2CCH +H - CH3CCH) -(CH2CCH_0      + H_0      - CH2CCH_0     ));

CH3CC_to_CH3_CC_H_sloppe 			= 0.5*((CH3+CC - CH3CC) - (CH3_0      +CC_0      -   CH3CC_0     ));
CH3CC_to_CH3C_C_H_sloppe 			= 0.5*((CH3C +C - CH3CC) - (CH3C_0      +C_0      - CH3CC_0     ));
CH3CC_to_CH2CC_H_H_sloppe 			= 0.5*((CH2CC +H - CH3CC) - (CH2CC_0      +H_0      - CH3CC_0     ));
CH2CHC_to_CH2_CHC_H_sloppe 			= 0.5*((CH2 + CHC - CH2CHC) -(CH2_0      +CHC_0      - CH2CHC_0     ));
CH2CHC_to_CH2CH_C_H_sloppe 			= 0.5*((CH2CH + C - CH2CHC)-(CH2CH_0      + C_0      - CH2CHC_0     ));
CH2CHC_to_CH2CC_H_H_sloppe 			= 0.5*((CH2CC +H - CH2CHC) -(CH2CC_0      + H_0      - CH2CHC_0     ));
CH2CHC_to_CHCHC_H_H_sloppe 			= 0.5*((CHCHC +H - CH2CHC) - (CHCHC_0      + H_0      - CH2CHC_0     ));
CHCH2C_to_CH_CH2C_H_sloppe 			= 0.5*((CH2C + CH - CHCH2C) - (CH2C_0      + CH_0      - CHCH2C_0     ));
CHCH2C_to_CH2CH_C_H_sloppe 			= 0.5*((CH2CH + C - CHCH2C) - (CH2CH_0      + C_0      - CHCH2C_0     ));
CHCH2C_to_CHCHC_H_H_sloppe 			= 0.5*((CHCHC + H - CHCH2C) - (CHCHC_0      + H_0      - CHCH2C_0     ));
CHCH2C_to_CCH2C_H_H_sloppe 			= 0.5*((CCH2C +H - CHCH2C) - (CCH2C_0      + H_0      - CHCH2C_0     ));
CHCHCH_to_CH_CHCH_H_sloppe 			= 0.5*((CH + CHCH - CHCHCH) - (CH_0      +CHCH_0      - CHCHCH_0     ));
CHCHCH_to_CHCHC_H_H_sloppe 			= 0.5*((CHCHC +H - CHCHCH) - (CHCHC_0      +H_0      -CHCHCH_0     ));
CHCHCH_to_CHCCH_H_H_sloppe 			= 0.5*((CHCCH +H - CHCHCH) - (CHCCH_0      +H_0      -CHCHCH_0     ));
CH2CCH_to_CH2_CCH_H_sloppe 			= 0.5*((CH2+CHC - CH2CCH) -(CH2_0      + CHC_0      - CH2CCH_0     ));
CH2CCH_to_CH2C_CH_H_sloppe 			= 0.5*((CH2C +CH -CH2CCH) -(CH2C_0      + CH_0      - CH2CCH_0     ));
CH2CCH_to_CH2CC_H_H_sloppe 			= 0.5*((CH2CC + H - CH2CCH) - (CH2CC_0      + H_0      - CH2CCH_0     ));
CH2CCH_to_CHCCH_H_H_sloppe 			= 0.5*((CHCCH + H - CH2CCH) - (CHCCH_0      + H_0      - CH2CCH_0     ));

CH2CC_to_CH2_CC_H_sloppe 			= 0.5*((CH2 + CC - CH2CC) -(CH2_0      +CC_0      -CH2CC_0     ));
CH2CC_to_CH2C_C_H_sloppe 			= 0.5*((CH2C + C - CH2CC) -(CH2C_0      +C_0      - CH2CC_0     ));
CH2CC_to_CHCC_H_H_sloppe 			= 0.5*((CHCC + H -CH2CC) - (CHCC_0      +H_0      - CH2CC_0     ));
CHCHC_to_CH_CHC_H_sloppe 			= 0.5*((CH+CHC - CHCHC) - (CH_0      + CHC_0      - CHCHC_0     ));
CHCHC_to_CHCH_C_H_sloppe 			= 0.5*((CHCH + C - CHCHC) - (CHCH_0      +C_0      - CHCHC_0     ));
CHCHC_to_CHCC_H_H_sloppe 			= 0.5*((CHCC + H - CHCHC) - (CHCC_0      +H_0      - CHCHC_0     ));
CHCHC_to_CCHC_H_H_sloppe 			= 0.5*((CCHC + H - CHCHC) - (CCHC_0      +H_0      - CHCHC_0     ));
CCH2C_to_C_CH2C_H_sloppe 			= 0.5*((CH2C + C - CCH2C) - (CH2C_0      + C_0      - CCH2C_0     ));
CCH2C_to_CCHC_H_H_sloppe 			= 0.5*((CCHC + H - CCH2C) - (CCHC_0      + H_0      - CCH2C_0     ));
CHCCH_to_CH_CHC_H_sloppe 			= 0.5*((CH + CHC - CHCCH) - (CH_0      + CHC_0      - CHCCH_0     ));
CHCCH_to_CHCC_H_H_sloppe 			= 0.5*((CHCC + H - CHCCH) - (CHCC_0      + H_0      - CHCCH_0     ));
CCHC_to_CHC_C_H_sloppe 				= 0.5*((CHC +C - CCHC) - (CHC_0      + C_0      - CCHC_0     ));
CCHC_to_CCC_H_H_sloppe 				= 0.5*((CCC +H - CCHC) - (CCC_0      +H_0      - CCHC_0     ));
CHCC_to_CH_CC_H_sloppe 				= 0.5*((CH +CC - CHCC) - (CH_0      + CC_0      - CHCC_0     ));
CHCC_to_CHC_C_H_sloppe 				= 0.5*((CHC +C - CHCC) - (CHC_0      + C_0      - CHCC_0     ));
CHCC_to_CCC_H_H_sloppe 				= 0.5*((CCC +H - CHCC) - (CCC_0      + H_0      - CHCC_0     ));
CCC_to_CC_C_H_sloppe 				= 0.5*((CC+C - CCC) - (CC_0      + C_0      - CCC_0     ));

CH3CH3_to_CH3_CH3_H_sloppe 			= 0.5*((CH3 +CH3 - CH3CH3) - (CH3_0      + CH3_0      -CH3CH3_0     ));
CH3CH3_to_CH3CH2_H_H_sloppe 		= 0.5*((CH3CH2 + H - CH3CH3) - (CH3CH2_0      +H_0      - CH3CH3_0     ));
CH3CH2_to_CH3_CH2_H_sloppe			= 0.5*((CH3+CH2 - CH3CH2) -( CH3_0      + CH2_0      - CH3CH2_0     ));
CH3CH2_to_CH3CH_H_H_sloppe 			= 0.5*((CH3CH + H - CH3CH2) - (CH3CH_0      +H_0      - CH3CH2_0     ));
CH3CH2_to_CH2CH2_H_H_sloppe 		= 0.5*((CH2CH2 +H - CH3CH2) - (CH2CH2_0      + H_0      - CH3CH2_0     ));
CH3CH_to_CH3_CH_H_sloppe 			= 0.5*((CH3 + CH - CH3CH) - (CH3_0      + CH_0      - CH3CH_0     ));
CH3CH_to_CH3C_H_H_sloppe 			= 0.5*((CH3C + H - CH3CH) - (CH3C_0      + H_0      - CH3CH_0     ));
CH3CH_to_CH2CH_H_H_sloppe 			= 0.5*((CH2CH +H - CH3CH) - (CH2CH_0      +H_0      - CH3CH_0     ));
CH3C_to_CH3_C_H_sloppe 				= 0.5*((CH3+C-CH3C) - (CH3_0      + C_0      - CH3C_0     ));
CH3C_to_CH2C_H_H_sloppe 			= 0.5 *(( CH2C + H - CH3C) - (CH2C_0    + H_0      - CH3C_0     ));
CH2CH2_to_CH2_CH2_H_sloppe 			= 0.5 *((CH2 + CH2 - CH2CH2) - (2*CH2_0      - CH2CH2_0     ));
CH2CH2_to_CH2CH_H_H_sloppe 			= 0.5*((CH2CH +H - CH2CH2) - (CH2CH_0      + H_0      - CH2CH2_0     ));
CH2CH_to_CH2_CH_H_sloppe 			= 0.5*((CH2 + CH - CH2CH) - (CH2_0      + CH_0      -CH2CH_0     ));
CH2CH_to_CH2C_H_H_sloppe 			= 0.5*((CH2C + H - CH2CH) - (CH2C_0      + H_0      - CH2CH_0     ));
CH2CH_to_CHCH_H_H_sloppe 			= 0.5*((CHCH + H - CH2CH) - (CHCH_0      + H_0      - CH2CH_0     ));
CH2C_to_CH2_C_H_sloppe 				= 0.5*((CH2 + C - CH2C) - (CH2_0      +C_0      - CH2C_0     ));
CH2C_to_CHC_H_H_sloppe 				= 0.5*((CHC + H - CH2C) - (CHC_0      +H_0      - CH2C_0     ));
CHCH_to_CH_CH_H_sloppe 				= 0.5*((CH+CH-CHCH) - (2*CH_0      - CHCH_0     ));
CHCH_to_CHC_H_H_sloppe 				= 0.5*((CHC+H - CHCH) - (CHC_0      +H_0      - CHCH_0     ));
CHC_to_CH_C_H_sloppe 				= 0.5*((CH +C -CHC) - (CH_0      +C_0      - CHC_0     ));
CHC_to_CC_H_H_sloppe 				= 0.5*((CC +H - CHC) - (CC_0      + H_0      - CHC_0     ));
CC_to_C_C_H_sloppe 					= 0.5*((2*C - CC) - (2*C_0      - CC_0     ));

CH4_to_CH3_H_H_sloppe 				= 0.5*((CH3+H-CH4) - (CH3_0      +H_0      - CH4_0     ));
CH3_to_CH2_H_H_sloppe 				= 0.5*((CH2+H-CH3) - (CH2_0      +H_0      - CH3_0     ));
CH2_to_CH_H_H_sloppe 				= 0.5*((CH+H-CH2) - (CH_0      +H_0      - CH2_0     ));
CH_to_C_H_H_sloppe 					= 0.5*((C+H-CH) - (C_0      +H_0      - CH_0     ));


% Correction for TST energy so that TST_energy is >or= ads_energy and TST_energy is > 0

CH3CH2CH3_to_CH3CHCH3_H = pp(47)+ ddEA9 + CH3CH2CH3_to_CH3CHCH3_H_H_sloppe - CH3CH2CH3_0;
if (CH3CH2CH3_to_CH3CHCH3_H < (CH3CHCH3 + H - CH3CH2CH3))
   CH3CH2CH3_to_CH3CHCH3_H = CH3CHCH3 + H - CH3CH2CH3;
end
if (CH3CH2CH3_to_CH3CHCH3_H <0)
   CH3CH2CH3_to_CH3CHCH3_H = 0.001;
end

CH3CH2CH3_to_CH3CH2CH2_H = pp(48)+ ddEA10 + CH3CH2CH3_to_CH3CH2CH2_H_H_sloppe - CH3CH2CH3_0;
if ( CH3CH2CH3_to_CH3CH2CH2_H <( CH3CH2CH2 + H - CH3CH2CH3))
   CH3CH2CH3_to_CH3CH2CH2_H = CH3CH2CH2 + H - CH3CH2CH3;
end
if (CH3CH2CH3_to_CH3CH2CH2_H <0)
        CH3CH2CH3_to_CH3CH2CH2_H = 0.001;
end

CH3CHCH3_to_CH3CHCH2_H = pp(49)+ ddEA11 + CH3CHCH3_to_CH3CHCH2_H_H_sloppe - CH3CHCH3_0;
if (CH3CHCH3_to_CH3CHCH2_H < (CH3CHCH2 + H - CH3CHCH3)) 
    CH3CHCH3_to_CH3CHCH2_H = CH3CHCH2 + H - CH3CHCH3;
end
if (CH3CHCH3_to_CH3CHCH2_H <0)
        CH3CHCH3_to_CH3CHCH2_H = 0.001;
end

CH3CH2CH2_to_CH3CHCH2_H = pp(50)+ ddEA12 + CH3CH2CH2_to_CH3CHCH2_H_H_sloppe - CH3CH2CH2_0;
if (CH3CH2CH2_to_CH3CHCH2_H < (CH3CHCH2 + H - CH3CH2CH2))
    CH3CH2CH2_to_CH3CHCH2_H = CH3CHCH2 + H - CH3CH2CH2;
end
if (CH3CH2CH2_to_CH3CHCH2_H <0)
        CH3CH2CH2_to_CH3CHCH2_H = 0.001;
end

CH3CH2CH3_to_CH3CH2_CH3 = pp(51) +ddEA13 + CH3CH2CH3_to_CH3_CH3CH2_H_sloppe - CH3CH2CH3_0;
if (CH3CH2CH3_to_CH3CH2_CH3 < (CH3 + CH3CH2 - CH3CH2CH3)) 
    CH3CH2CH3_to_CH3CH2_CH3  = CH3 + CH3CH2 - CH3CH2CH3;
end
if (CH3CH2CH3_to_CH3CH2_CH3 <0)
        CH3CH2CH3_to_CH3CH2_CH3 =0.001;
end

CH3CHCH3_to_CH3_CH3CH = pp(52) + ddEA14 + CH3CHCH3_to_CH3_CH3CH_H_sloppe - CH3CHCH3_0;
if (CH3CHCH3_to_CH3_CH3CH <( CH3CH + CH3 - CH3CHCH3))
   CH3CHCH3_to_CH3_CH3CH = CH3CH + CH3 - CH3CHCH3;
end
if (CH3CHCH3_to_CH3_CH3CH <0 )
        CH3CHCH3_to_CH3_CH3CH =0.001;
end

CH3CHCH3_to_CH3CCH3_H = pp(53) + ddEA15 + CH3CHCH3_to_CH3CCH3_H_H_sloppe - CH3CHCH3_0;
if (CH3CHCH3_to_CH3CCH3_H < (CH3CCH3 + H - CH3CHCH3 ))
	CH3CHCH3_to_CH3CCH3_H = CH3CCH3 + H - CH3CHCH3;
end
if (CH3CHCH3_to_CH3CCH3_H <0 )
        CH3CHCH3_to_CH3CCH3_H =0.001;
end

CH3CH2CH2_to_CH3CH2_CH2 = pp(54) + ddEA16 + CH3CH2CH2_to_CH3CH2_CH2_H_sloppe - CH3CH2CH2_0;
if (CH3CH2CH2_to_CH3CH2_CH2 < (CH3CH2 + CH2 - CH3CH2CH2))
	CH3CH2CH2_to_CH3CH2_CH2 = CH3CH2 + CH2 - CH3CH2CH2;
end
if (CH3CH2CH2_to_CH3CH2_CH2 <0)
        CH3CH2CH2_to_CH3CH2_CH2 = 0.001;
end

CH3CH2CH2_to_CH3_CH2CH2 = pp(55) + ddEA17 + CH3CH2CH2_to_CH3_CH2CH2_H_sloppe- CH3CH2CH2_0;
if (CH3CH2CH2_to_CH3_CH2CH2 < (CH3 + CH2CH2 - CH3CH2CH2)) 
	CH3CH2CH2_to_CH3_CH2CH2 = CH3 + CH2CH2 - CH3CH2CH2;
end
if (CH3CH2CH2_to_CH3_CH2CH2 <0)
        CH3CH2CH2_to_CH3_CH2CH2 =0.001;
end

CH3CH2CH2_to_CH2CH2CH2_H = pp(56) + ddEA18 +CH3CH2CH2_to_CH2CH2CH2_H_H_sloppe - CH3CH2CH2_0;
if (CH3CH2CH2_to_CH2CH2CH2_H < (CH2CH2CH2 + H - CH3CH2CH2))
	CH3CH2CH2_to_CH2CH2CH2_H = CH2CH2CH2 + H - CH3CH2CH2;
end
if (CH3CH2CH2_to_CH2CH2CH2_H <0)
        CH3CH2CH2_to_CH2CH2CH2_H =0.001;
end

CH3CH2CH2_to_CH3CH2CH_H = pp(57) + ddEA19 + CH3CH2CH2_to_CH3CH2CH_H_H_sloppe - CH3CH2CH2_0;
if (CH3CH2CH2_to_CH3CH2CH_H < (CH3CH2CH + H - CH3CH2CH2))
	CH3CH2CH2_to_CH3CH2CH_H = CH3CH2CH + H - CH3CH2CH2;
end
if (CH3CH2CH2_to_CH3CH2CH_H <0)
        CH3CH2CH2_to_CH3CH2CH_H =0.001;
end

CH3CH2CH_to_CH3CH2_CH = pp(58) +ddEA20 + CH3CH2CH_to_CH3CH2_CH_H_sloppe - CH3CH2CH_0;
if (CH3CH2CH_to_CH3CH2_CH < (CH3CH2 + CH - CH3CH2CH))
    CH3CH2CH_to_CH3CH2_CH = CH3CH2 + CH - CH3CH2CH;
end
if (CH3CH2CH_to_CH3CH2_CH <0)
        CH3CH2CH_to_CH3CH2_CH =0.001;
end

CH3CH2CH_to_CH3_CH2CH = pp(59) + ddEA21 + CH3CH2CH_to_CH3_CH2CH_H_sloppe - CH3CH2CH_0;
if (CH3CH2CH_to_CH3_CH2CH < (CH2CH + CH3 - CH3CH2CH)) 
	CH3CH2CH_to_CH3_CH2CH = CH2CH + CH3 - CH3CH2CH;
end
if (CH3CH2CH_to_CH3_CH2CH <0)
        CH3CH2CH_to_CH3_CH2CH =0.001;
end

CH3CH2CH_to_CH3CH2C_H = pp(60) + ddEA22 + CH3CH2CH_to_CH3CH2C_H_H_sloppe - CH3CH2CH_0;
if (CH3CH2CH_to_CH3CH2C_H < (CH3CH2C + H - CH3CH2CH))
	CH3CH2CH_to_CH3CH2C_H = CH3CH2C + H - CH3CH2CH;
end
if (CH3CH2CH_to_CH3CH2C_H <0)
        CH3CH2CH_to_CH3CH2C_H =0.001;
end

CH3CH2CH_to_CH3CHCH_H = pp(61) + ddEA23 + CH3CH2CH_to_CH3CHCH_H_H_sloppe - CH3CH2CH_0;
if (CH3CH2CH_to_CH3CHCH_H < (CH3CHCH + H - CH3CH2CH))
    CH3CH2CH_to_CH3CHCH_H = (CH3CHCH + H - CH3CH2CH); 
 end   
if (CH3CH2CH_to_CH3CHCH_H <0)
        CH3CH2CH_to_CH3CHCH_H =0.001;
end

  CH3CH2CH_to_CH2CH2CH_H = pp(62) + ddEA24 + CH3CH2CH_to_CH2CH2CH_H_H_sloppe - CH3CH2CH_0;
if (CH3CH2CH_to_CH2CH2CH_H < (CH2CH2CH + H - CH3CH2CH))
	CH3CH2CH_to_CH2CH2CH_H = (CH2CH2CH + H - CH3CH2CH);
end
if (CH3CH2CH_to_CH2CH2CH_H <0)
        CH3CH2CH_to_CH2CH2CH_H =0.001;
end

  CH2CH2CH2_to_CH2CH2_CH2 = pp(63) + ddEA25 + CH2CH2CH2_to_CH2_CH2CH2_H_sloppe - CH2CH2CH2_0;
if (CH2CH2CH2_to_CH2CH2_CH2 < (CH2CH2 + CH2 - CH2CH2CH2))
    CH2CH2CH2_to_CH2CH2_CH2 = (CH2CH2 + CH2 - CH2CH2CH2); 
end
if (CH2CH2CH2_to_CH2CH2_CH2 <0)
        CH2CH2CH2_to_CH2CH2_CH2 =0.001;
end

  CH2CH2CH2_to_CH2CH2CH_H = pp(64) + ddEA26 + CH2CH2CH2_to_CH2CH2CH_H_H_sloppe - CH2CH2CH2_0;
if (CH2CH2CH2_to_CH2CH2CH_H < (CH2CH2CH + H - CH2CH2CH2))
	CH2CH2CH2_to_CH2CH2CH_H = (CH2CH2CH + H - CH2CH2CH2);
end
if (CH2CH2CH2_to_CH2CH2CH_H <0)
        CH2CH2CH2_to_CH2CH2CH_H =0.001;
end

  CH2CH2CH2_to_CH2CHCH2_H = pp(65) + ddEA27 + CH2CH2CH2_to_CH2CHCH2_H_H_sloppe - CH2CH2CH2_0;
if (CH2CH2CH2_to_CH2CHCH2_H < (CH2CHCH2 + H - CH2CH2CH2))
	CH2CH2CH2_to_CH2CHCH2_H = (CH2CHCH2 + H - CH2CH2CH2);
end
if (CH2CH2CH2_to_CH2CHCH2_H <0)
        CH2CH2CH2_to_CH2CHCH2_H =0.001;
end

CH3CHCH2_to_CH3_CHCH2 = pp(66) + ddEA28 + CH3CHCH2_to_CH3_CHCH2_H_sloppe - CH3CHCH2_0;
if (CH3CHCH2_to_CH3_CHCH2 < (CH3 + CH2CH - CH3CHCH2))
	CH3CHCH2_to_CH3_CHCH2 = (CH3 + CH2CH - CH3CHCH2);
end
if (CH3CHCH2_to_CH3_CHCH2 <0)
        CH3CHCH2_to_CH3_CHCH2 =0.001;
end

  CH3CHCH2_to_CH3CH_CH2 = pp(67) + ddEA29 + CH3CHCH2_to_CH3CH_CH2_H_sloppe - CH3CHCH2_0;
if (CH3CHCH2_to_CH3CH_CH2 < (CH3CH + CH2 - CH3CHCH2))
	CH3CHCH2_to_CH3CH_CH2 = (CH3CH + CH2 - CH3CHCH2);
end
if (CH3CHCH2_to_CH3CH_CH2<0)
        CH3CHCH2_to_CH3CH_CH2 =0.001;
end

  CH3CHCH2_to_CH3CCH2_H = pp(68) + ddEA30 + CH3CHCH2_to_CH3CCH2_H_H_sloppe - CH3CHCH2_0;
if (CH3CHCH2_to_CH3CCH2_H < (CH3CCH2 + H - CH3CHCH2))
	CH3CHCH2_to_CH3CCH2_H = (CH3CCH2 + H - CH3CHCH2);
end
if (CH3CHCH2_to_CH3CCH2_H <0)
        CH3CHCH2_to_CH3CCH2_H =0.001;
end

  CH3CHCH2_to_CH3CHCH_H = pp(69) + ddEA31 + CH3CHCH2_to_CH3CHCH_H_H_sloppe - CH3CHCH2_0;
if (CH3CHCH2_to_CH3CHCH_H < (CH3CHCH + H - CH3CHCH2))
	CH3CHCH2_to_CH3CHCH_H = (CH3CHCH + H - CH3CHCH2);
end
if (CH3CHCH2_to_CH3CHCH_H <0)
        CH3CHCH2_to_CH3CHCH_H =0.001;
end

  CH3CHCH2_to_CH2CHCH2_H = pp(70) +ddEA32 + CH3CHCH2_to_CH2CHCH2_H_H_sloppe - CH3CHCH2_0;
if (CH3CHCH2_to_CH2CHCH2_H < (CH2CHCH2 + H - CH3CHCH2))
	CH3CHCH2_to_CH2CHCH2_H = (CH2CHCH2 + H - CH3CHCH2);
end
if (CH3CHCH2_to_CH2CHCH2_H <0)
        CH3CHCH2_to_CH2CHCH2_H =0.001;
end

  CH3CCH3_to_CH3_CCH3 = pp(71) + ddEA33 + CH3CCH3_to_CH3_CCH3_H_sloppe - CH3CCH3_0;
if (CH3CCH3_to_CH3_CCH3 < (CH3C + CH3 - CH3CCH3))
	CH3CCH3_to_CH3_CCH3 = (CH3C + CH3 - CH3CCH3);
end
if (CH3CCH3_to_CH3_CCH3 <0)
        CH3CCH3_to_CH3_CCH3 =0.001;
end

  CH3CCH3_to_CH3CCH2_H = pp(72) + ddEA34 + CH3CCH3_to_CH3CCH2_H_H_sloppe - CH3CCH3_0;
if (CH3CCH3_to_CH3CCH2_H < (CH3CCH2 + H - CH3CCH3))
	CH3CCH3_to_CH3CCH2_H = (CH3CCH2 + H - CH3CCH3);
end
if (CH3CCH3_to_CH3CCH2_H <0 )
        CH3CCH3_to_CH3CCH2_H = 0.001;
end

  CH3CH2C_to_CH3_CH2C = pp(73) + ddEA35 + CH3CH2C_to_CH3_CH2C_H_sloppe - CH3CH2C_0;
if (CH3CH2C_to_CH3_CH2C < (CH3 + CH2C - CH3CH2C))
	CH3CH2C_to_CH3_CH2C = (CH3 + CH2C - CH3CH2C);
end
if (CH3CH2C_to_CH3_CH2C <0)
        CH3CH2C_to_CH3_CH2C = 0.001;
end

  CH3CH2C_to_CH3CH2_C = pp(74) + ddEA36 + CH3CH2C_to_CH3CH2_C_H_sloppe - CH3CH2C_0;
if (CH3CH2C_to_CH3CH2_C < (CH3CH2 + C - CH3CH2C))
	CH3CH2C_to_CH3CH2_C = (CH3CH2 + C - CH3CH2C);
end
if (CH3CH2C_to_CH3CH2_C <0)
        CH3CH2C_to_CH3CH2_C =0.001;
end

  CH3CH2C_to_CH2CH2C_H = pp(75) + ddEA37 + CH3CH2C_to_CH2CH2C_H_H_sloppe - CH3CH2C_0;
if (CH3CH2C_to_CH2CH2C_H < (CH2CH2C + H - CH3CH2C))
	CH3CH2C_to_CH2CH2C_H = (CH2CH2C + H - CH3CH2C);
end
if (CH3CH2C_to_CH2CH2C_H <0)
        CH3CH2C_to_CH2CH2C_H =0.001;
end

  CH3CH2C_to_CH3CHC_H = pp(76) + ddEA38 + CH3CH2C_to_CH3CHC_H_H_sloppe - CH3CH2C_0;
if (CH3CH2C_to_CH3CHC_H < (CH3CHC + H - CH3CH2C))
	CH3CH2C_to_CH3CHC_H = (CH3CHC + H - CH3CH2C);
end
if (CH3CH2C_to_CH3CHC_H <0 )
        CH3CH2C_to_CH3CHC_H =0.001;
end

  CH2CH2CH_to_CH2_CH2CH = pp(77) + ddEA39 + CH2CH2CH_to_CH2_CH2CH_H_sloppe - CH2CH2CH_0;
if (CH2CH2CH_to_CH2_CH2CH < (CH2CH + CH2 - CH2CH2CH))
	CH2CH2CH_to_CH2_CH2CH = (CH2CH + CH2 - CH2CH2CH);
end
if (CH2CH2CH_to_CH2_CH2CH <0)
        CH2CH2CH_to_CH2_CH2CH =0.001;
end

  CH2CH2CH_to_CH2CH2_CH = pp(78) + ddEA40 + CH2CH2CH_to_CH2CH2_CH_H_sloppe - CH2CH2CH_0;
if (CH2CH2CH_to_CH2CH2_CH < (CH2CH2 + CH - CH2CH2CH))
	CH2CH2CH_to_CH2CH2_CH = (CH2CH2 + CH - CH2CH2CH);
end
if (CH2CH2CH_to_CH2CH2_CH <0)
        CH2CH2CH_to_CH2CH2_CH =0.001;
end

  CH2CH2CH_to_CH2CH2C_H = pp(79) + ddEA41 + CH2CH2CH_to_CH2CH2C_H_H_sloppe  - CH2CH2CH_0;
if (CH2CH2CH_to_CH2CH2C_H < (CH2CH2C + H - CH2CH2CH))
	CH2CH2CH_to_CH2CH2C_H = (CH2CH2C + H - CH2CH2CH);
end
if (CH2CH2CH_to_CH2CH2C_H <0)
        CH2CH2CH_to_CH2CH2C_H =0.001;
end

  CH2CH2CH_to_CH2CHCH_H = pp(80) + ddEA42 + CH2CH2CH_to_CH2CHCH_H_H_sloppe  - CH2CH2CH_0;
if (CH2CH2CH_to_CH2CHCH_H < (CH2CHCH +H - CH2CH2CH))
	CH2CH2CH_to_CH2CHCH_H = (CH2CHCH +H - CH2CH2CH);
end
if (CH2CH2CH_to_CH2CHCH_H <0)
        CH2CH2CH_to_CH2CHCH_H =0.001;
end

  CH2CH2CH_to_CHCH2CH_H = pp(81) + ddEA43 + CH2CH2CH_to_CHCH2CH_H_H_sloppe  - CH2CH2CH_0;
if (CH2CH2CH_to_CHCH2CH_H < (CHCH2CH + H - CH2CH2CH))
	CH2CH2CH_to_CHCH2CH_H = (CHCH2CH + H - CH2CH2CH);
end
if (CH2CH2CH_to_CHCH2CH_H <0)
        CH2CH2CH_to_CHCH2CH_H =0.001;
end

  CH2CHCH2_to_CH2_CH2CH = pp(82) + ddEA44 + CH2CHCH2_to_CH2_CHCH2_H_sloppe - CH2CHCH2_0;
if (CH2CHCH2_to_CH2_CH2CH < (CH2 + CH2CH - CH2CHCH2))
	CH2CHCH2_to_CH2_CH2CH = (CH2 + CH2CH - CH2CHCH2);
end
if (CH2CHCH2_to_CH2_CH2CH <0)
        CH2CHCH2_to_CH2_CH2CH = 0.001;
end

  CH2CHCH2_to_CH2CHCH_H = pp(83) + ddEA45 + CH2CHCH2_to_CH2CHCH_H_H_sloppe - CH2CHCH2_0;
if (CH2CHCH2_to_CH2CHCH_H < (CH2CHCH + H - CH2CHCH2))
	CH2CHCH2_to_CH2CHCH_H = (CH2CHCH + H - CH2CHCH2);
end
if (CH2CHCH2_to_CH2CHCH_H <0)
        CH2CHCH2_to_CH2CHCH_H = 0.001;
end

  CH2CHCH2_to_CH2CCH2_H = pp(84) + ddEA46 + CH2CHCH2_to_CH2CCH2_H_H_sloppe - CH2CHCH2_0;
if (CH2CHCH2_to_CH2CCH2_H < (CH2CCH2 + H - CH2CHCH2))
	CH2CHCH2_to_CH2CCH2_H = (CH2CCH2 + H - CH2CHCH2);
end
if (CH2CHCH2_to_CH2CCH2_H <0)
        CH2CHCH2_to_CH2CCH2_H =0.001;
end

  CH3CHCH_to_CH3_CHCH = pp(85) + ddEA47 + CH3CHCH_to_CH3_CHCH_H_sloppe - CH3CHCH_0;
if (CH3CHCH_to_CH3_CHCH < ( CH3 + CHCH - CH3CHCH))
	CH3CHCH_to_CH3_CHCH = ( CH3 + CHCH - CH3CHCH);
end
if (CH3CHCH_to_CH3_CHCH <0)
        CH3CHCH_to_CH3_CHCH = 0.001;
end

  CH3CHCH_to_CH3CH_CH = pp(86) + ddEA48 + CH3CHCH_to_CH3CH_CH_H_sloppe - CH3CHCH_0;
if (CH3CHCH_to_CH3CH_CH < (CH3CH + CH - CH3CHCH))
	CH3CHCH_to_CH3CH_CH = (CH3CH + CH - CH3CHCH);
end
if (CH3CHCH_to_CH3CH_CH <0)
        CH3CHCH_to_CH3CH_CH =0.001;
end

  CH3CHCH_to_CH3CHC_H = pp(87) + ddEA49 + CH3CHCH_to_CH3CHC_H_H_sloppe - CH3CHCH_0;
if (CH3CHCH_to_CH3CHC_H < (CH3CHC + H - CH3CHCH))
	CH3CHCH_to_CH3CHC_H = (CH3CHC + H - CH3CHCH) ;
end
if (CH3CHCH_to_CH3CHC_H <0)
        CH3CHCH_to_CH3CHC_H =0.001;
end

  CH3CHCH_to_CH3CCH_H = pp(88) + ddEA50 + CH3CHCH_to_CH3CCH_H_H_sloppe - CH3CHCH_0;
if (CH3CHCH_to_CH3CCH_H < ( CH3CCH + H - CH3CHCH))
	CH3CHCH_to_CH3CCH_H = ( CH3CCH + H - CH3CHCH);
end
if (CH3CHCH_to_CH3CCH_H <0)
        CH3CHCH_to_CH3CCH_H =0.001;
end

  CH3CHCH_to_CH2CHCH_H = pp(89) + ddEA51 + CH3CHCH_to_CH2CHCH_H_H_sloppe - CH3CHCH_0;
if (CH3CHCH_to_CH2CHCH_H < (CH2CHCH + H - CH3CHCH))
	CH3CHCH_to_CH2CHCH_H = (CH2CHCH + H - CH3CHCH);
end
if (CH3CHCH_to_CH2CHCH_H <0)
        CH3CHCH_to_CH2CHCH_H=0.001;
end

  CH3CCH2_to_CH3_CCH2 = pp(90) + ddEA52 + CH3CCH2_to_CH3_CCH2_H_sloppe - CH3CCH2_0;
if (CH3CCH2_to_CH3_CCH2 < (CH2C + CH3 - CH3CCH2))
	CH3CCH2_to_CH3_CCH2 = (CH2C + CH3 - CH3CCH2);
end
if (CH3CCH2_to_CH3_CCH2 <0)
        CH3CCH2_to_CH3_CCH2 = 0.001;
end

  CH3CCH2_to_CH3C_CH2 = pp(91) + ddEA53 + CH3CCH2_to_CH3C_CH2_H_sloppe - CH3CCH2_0;
if (CH3CCH2_to_CH3C_CH2 < (CH3C + CH2 - CH3CCH2))
	CH3CCH2_to_CH3C_CH2 = (CH3C + CH2 - CH3CCH2);
end
if (CH3CCH2_to_CH3C_CH2 <0)
        CH3CCH2_to_CH3C_CH2 =0.001;
end

  CH3CCH2_to_CH2CCH2_H = pp(92) + ddEA54 + CH3CCH2_to_CH2CCH2_H_H_sloppe - CH3CCH2_0;
if (CH3CCH2_to_CH2CCH2_H < (CH2CCH2 + H - CH3CCH2))
	CH3CCH2_to_CH2CCH2_H = (CH2CCH2 + H - CH3CCH2);
end
if (CH3CCH2_to_CH2CCH2_H <0)
        CH3CCH2_to_CH2CCH2_H =0.001;
end

  CH3CCH2_to_CH3CCH_H = pp(93) + ddEA55 + CH3CCH2_to_CH3CCH_H_H_sloppe - CH3CCH2_0;
if (CH3CCH2_to_CH3CCH_H < (CH3CCH + H - CH3CCH2))
	CH3CCH2_to_CH3CCH_H = (CH3CCH + H - CH3CCH2);
end
if (CH3CCH2_to_CH3CCH_H <0)
        CH3CCH2_to_CH3CCH_H = 0.001;
end

  CH3CHC_to_CH3_CHC = pp(94) + ddEA56 + CH3CHC_to_CH3_CHC_H_sloppe - CH3CHC_0;
if (CH3CHC_to_CH3_CHC < (CH3 + CHC - CH3CHC ))
	CH3CHC_to_CH3_CHC = (CH3 + CHC - CH3CHC );
end
if (CH3CHC_to_CH3_CHC <0)
        CH3CHC_to_CH3_CHC =0.001;
end

  CH3CHC_to_CH3CH_C = pp(95) + ddEA57 + CH3CHC_to_CH3CH_C_H_sloppe - CH3CHC_0;
if (CH3CHC_to_CH3CH_C < (CH3CH + C - CH3CHC))
	CH3CHC_to_CH3CH_C = (CH3CH + C - CH3CHC);
end
if (CH3CHC_to_CH3CH_C <0)
        CH3CHC_to_CH3CH_C =0.001;
end

  CH3CHC_to_CH3CC_H = pp(96) + ddEA58 + CH3CHC_to_CH3CC_H_H_sloppe - CH3CHC_0;
if (CH3CHC_to_CH3CC_H < (CH3CC +H - CH3CHC))
	CH3CHC_to_CH3CC_H = (CH3CC +H - CH3CHC);
end
if (CH3CHC_to_CH3CC_H <0)
        CH3CHC_to_CH3CC_H =0.001;
end

  CH3CHC_to_CH2CHC_H = pp(97) + ddEA59 + CH3CHC_to_CH2CHC_H_H_sloppe - CH3CHC_0;
if (CH3CHC_to_CH2CHC_H < (CH2CHC + H - CH3CHC))
	CH3CHC_to_CH2CHC_H = (CH2CHC + H - CH3CHC);
end
if (CH3CHC_to_CH2CHC_H <0)
        CH3CHC_to_CH2CHC_H =0.001;
end

  CH2CH2C_to_CH2CH2_C = pp(98) + ddEA60 + CH2CH2C_to_CH2CH2_C_H_sloppe - CH2CH2C_0;
if (CH2CH2C_to_CH2CH2_C < (CH2CH2 + C - CH2CH2C))
	CH2CH2C_to_CH2CH2_C = (CH2CH2 + C - CH2CH2C);
end
if (CH2CH2C_to_CH2CH2_C <0)
        CH2CH2C_to_CH2CH2_C =0.001;
end

  CH2CH2C_to_CH2_CH2C = pp(99) + ddEA61 + CH2CH2C_to_CH2_CH2C_H_sloppe - CH2CH2C_0;
if (CH2CH2C_to_CH2_CH2C < (CH2 + CH2C - CH2CH2C))
	CH2CH2C_to_CH2_CH2C = (CH2 + CH2C - CH2CH2C);
end
if (CH2CH2C_to_CH2_CH2C <0)
        CH2CH2C_to_CH2_CH2C =0.001;
end

  CH2CH2C_to_CH2CHC_H = pp(100) + ddEA62 + CH2CH2C_to_CH2CHC_H_H_sloppe - CH2CH2C_0;
if (CH2CH2C_to_CH2CHC_H < (CH2CHC + H - CH2CH2C))
	CH2CH2C_to_CH2CHC_H = (CH2CHC + H - CH2CH2C);
end
if (CH2CH2C_to_CH2CHC_H <0)
        CH2CH2C_to_CH2CHC_H =0.001;
end

  CH2CH2C_to_CHCH2C_H = pp(101) + ddEA63 + CH2CH2C_to_CHCH2C_H_H_sloppe - CH2CH2C_0;
if (CH2CH2C_to_CHCH2C_H < (CHCH2C + H - CH2CH2C))
	CH2CH2C_to_CHCH2C_H = (CHCH2C + H - CH2CH2C);
end
if (CH2CH2C_to_CHCH2C_H <0)
        CH2CH2C_to_CHCH2C_H =0.001;
end

  CHCH2CH_to_CHCH2_CH = pp(102) + ddEA64 + CHCH2CH_to_CHCH2_CH_H_sloppe - CHCH2CH_0;
if (CHCH2CH_to_CHCH2_CH < (CH + CH2CH - CHCH2CH))
	CHCH2CH_to_CHCH2_CH = (CH + CH2CH - CHCH2CH);
end
if (CHCH2CH_to_CHCH2_CH <0)
        CHCH2CH_to_CHCH2_CH =0.001;
end

  CHCH2CH_to_CHCH2C_H = pp(103) + ddEA65 + CHCH2CH_to_CHCH2C_H_H_sloppe - CHCH2CH_0;
if (CHCH2CH_to_CHCH2C_H < (CHCH2C +H - CHCH2CH))
	CHCH2CH_to_CHCH2C_H = (CHCH2C +H - CHCH2CH);
end
if (CHCH2CH_to_CHCH2C_H <0)
        CHCH2CH_to_CHCH2C_H =0.001;
end

  CHCH2CH_to_CHCHCH_H = pp(104) + ddEA66 + CHCH2CH_to_CHCHCH_H_H_sloppe - CHCH2CH_0;
if (CHCH2CH_to_CHCHCH_H < (CHCHCH + H - CHCH2CH))
	CHCH2CH_to_CHCHCH_H = (CHCHCH + H - CHCH2CH);
end
if (CHCH2CH_to_CHCHCH_H <0)
        CHCH2CH_to_CHCHCH_H =0.001;
end

  CH2CHCH_to_CH2_CHCH = pp(105) + ddEA67 + CH2CHCH_to_CH2_CHCH_H_sloppe - CH2CHCH_0;
if (CH2CHCH_to_CH2_CHCH < (CH2+CHCH - CH2CHCH))
	CH2CHCH_to_CH2_CHCH = (CH2+CHCH - CH2CHCH) ;
end
if (CH2CHCH_to_CH2_CHCH <0)
        CH2CHCH_to_CH2_CHCH =0.001;
end

  CH2CHCH_to_CH2CH_CH = pp(106) + ddEA68 + CH2CHCH_to_CH2CH_CH_H_sloppe - CH2CHCH_0;
if (CH2CHCH_to_CH2CH_CH < (CH2CH + CH - CH2CHCH))
	CH2CHCH_to_CH2CH_CH = (CH2CH + CH - CH2CHCH);
end
if (CH2CHCH_to_CH2CH_CH <0)
        CH2CHCH_to_CH2CH_CH =0.001;
end

  CH2CHCH_to_CH2CHC_H = pp(107) + ddEA69 + CH2CHCH_to_CH2CHC_H_H_sloppe - CH2CHCH_0;
if (CH2CHCH_to_CH2CHC_H < (CH2CHC + H - CH2CHCH))
	CH2CHCH_to_CH2CHC_H = (CH2CHC + H - CH2CHCH);
end
if (CH2CHCH_to_CH2CHC_H <0)
        CH2CHCH_to_CH2CHC_H = 0.001;
end

  CH2CHCH_to_CH2CCH_H = pp(108) + ddEA70 + CH2CHCH_to_CH2CCH_H_H_sloppe - CH2CHCH_0;
if (CH2CHCH_to_CH2CCH_H < (CH2CCH + H - CH2CHCH))
	CH2CHCH_to_CH2CCH_H = (CH2CCH + H - CH2CHCH);
end
if (CH2CHCH_to_CH2CCH_H <0)
        CH2CHCH_to_CH2CCH_H =0.001;
end

  CH2CHCH_to_CHCHCH_H = pp(109) + ddEA71 + CH2CHCH_to_CHCHCH_H_H_sloppe - CH2CHCH_0;
if (CH2CHCH_to_CHCHCH_H < (CHCHCH +H - CH2CHCH))
	CH2CHCH_to_CHCHCH_H = (CHCHCH +H - CH2CHCH);
end
if (CH2CHCH_to_CHCHCH_H <0)
        CH2CHCH_to_CHCHCH_H =0.001;
end

  CH2CCH2_to_CH2C_CH2 = pp(110) + ddEA72 + CH2CCH2_to_CH2C_CH2_H_sloppe - CH2CCH2_0;
if (CH2CCH2_to_CH2C_CH2 < (CH2C + CH2 - CH2CCH2))
	CH2CCH2_to_CH2C_CH2 = (CH2C + CH2 - CH2CCH2);
end
if (CH2CCH2_to_CH2C_CH2 <0)
        CH2CCH2_to_CH2C_CH2 =0.001;
end

  CH2CCH2_to_CH2CCH_H = pp(111) + ddEA73 + CH2CCH2_to_CH2CCH_H_H_sloppe - CH2CCH2_0;
if (CH2CCH2_to_CH2CCH_H < (CH2CCH + H -CH2CCH2))
	CH2CCH2_to_CH2CCH_H = (CH2CCH + H -CH2CCH2);
end
if (CH2CCH2_to_CH2CCH_H <0)
        CH2CCH2_to_CH2CCH_H =0.001;
end

  CH3CCH_to_CH3C_CH = pp(112) + ddEA74 + CH3CCH_to_CH3C_CH_H_sloppe - CH3CCH_0;
if (CH3CCH_to_CH3C_CH < (CH3C + CH - CH3CCH))
	CH3CCH_to_CH3C_CH = (CH3C + CH - CH3CCH);
end
if (CH3CCH_to_CH3C_CH <0)
        CH3CCH_to_CH3C_CH =0.001;
end

  CH3CCH_to_CH3_CHC = pp(113) + ddEA75 + CH3CCH_to_CH3_CHC_H_sloppe - CH3CCH_0;
if (CH3CCH_to_CH3_CHC < (CHC + CH3 - CH3CCH))
	CH3CCH_to_CH3_CHC = (CHC + CH3 - CH3CCH);
end
if (CH3CCH_to_CH3_CHC<0)
        CH3CCH_to_CH3_CHC =0.001;
end

  CH3CCH_to_CH3CC_H = pp(114) + ddEA76 + CH3CCH_to_CH3CC_H_H_sloppe - CH3CCH_0;
if (CH3CCH_to_CH3CC_H < (CH3CC + H - CH3CCH))
	CH3CCH_to_CH3CC_H = (CH3CC + H - CH3CCH);
end
if (CH3CCH_to_CH3CC_H<0)
        CH3CCH_to_CH3CC_H =0.001;
end

  CH3CCH_to_CH2CCH_H = pp(115) + ddEA77 + CH3CCH_to_CH2CCH_H_H_sloppe - CH3CCH_0;
if (CH3CCH_to_CH2CCH_H < (CH2CCH + H - CH3CCH))
	CH3CCH_to_CH2CCH_H = (CH2CCH + H - CH3CCH); 
end
if (CH3CCH_to_CH2CCH_H <0)
        CH3CCH_to_CH2CCH_H =0.001;
end

  CH3CC_to_CH3_CC = pp(116) + ddEA78 + CH3CC_to_CH3_CC_H_sloppe - CH3CC_0;
if (CH3CC_to_CH3_CC < (CH3 + CC - CH3CC))
	CH3CC_to_CH3_CC = (CH3 + CC - CH3CC);
end
if (CH3CC_to_CH3_CC<0)
        CH3CC_to_CH3_CC =0.001;
end

  CH3CC_to_CH3C_C = pp(117) + ddEA79 + CH3CC_to_CH3C_C_H_sloppe - CH3CC_0;
if (CH3CC_to_CH3C_C < (CH3C + C - CH3CC))
	CH3CC_to_CH3C_C = (CH3C + C - CH3CC);
end
if (CH3CC_to_CH3C_C <0)
        CH3CC_to_CH3C_C = 0.001;
end

  CH3CC_to_CH2CC_H = pp(118) + ddEA80 + CH3CC_to_CH2CC_H_H_sloppe - CH3CC_0;
if (CH3CC_to_CH2CC_H < (CH2CC + H - CH3CC))
	CH3CC_to_CH2CC_H = (CH2CC + H - CH3CC);
end
if (CH3CC_to_CH2CC_H <0)
        CH3CC_to_CH2CC_H =0.001;
end

  CH2CHC_to_CH2_CHC = pp(119) + ddEA81 + CH2CHC_to_CH2_CHC_H_sloppe - CH2CHC_0;
if (CH2CHC_to_CH2_CHC < (CHC + CH2 - CH2CHC))
	CH2CHC_to_CH2_CHC = (CHC + CH2 - CH2CHC);
end
if (CH2CHC_to_CH2_CHC <0)
        CH2CHC_to_CH2_CHC =0.001;
end

  CH2CHC_to_CH2CH_C = pp(120) + ddEA82 + CH2CHC_to_CH2CH_C_H_sloppe - CH2CHC_0;
if (CH2CHC_to_CH2CH_C < (C + CH2CH - CH2CHC))
	CH2CHC_to_CH2CH_C = (C + CH2CH - CH2CHC);
end
if (CH2CHC_to_CH2CH_C <0)
        CH2CHC_to_CH2CH_C =0.001;
end

  CH2CHC_to_CH2CC_H = pp(121) + ddEA83 + CH2CHC_to_CH2CC_H_H_sloppe - CH2CHC_0;
if (CH2CHC_to_CH2CC_H < (CH2CC + H - CH2CHC))
	CH2CHC_to_CH2CC_H = (CH2CC + H - CH2CHC);
end
if (CH2CHC_to_CH2CC_H <0)
        CH2CHC_to_CH2CC_H =0.001;
end

  CH2CHC_to_CHCHC_H = pp(122) + ddEA84 + CH2CHC_to_CHCHC_H_H_sloppe - CH2CHC_0;
if (CH2CHC_to_CHCHC_H < (CHCHC + H - CH2CHC))
	CH2CHC_to_CHCHC_H = (CHCHC + H - CH2CHC);	
end
if (CH2CHC_to_CHCHC_H <0)
        CH2CHC_to_CHCHC_H =0.001;
end

 CHCH2C_to_CH_CH2C = pp(123) + ddEA85 + CHCH2C_to_CH_CH2C_H_sloppe - CHCH2C_0;
if (CHCH2C_to_CH_CH2C > (CH + CH2C - CHCH2C))
	CHCH2C_to_CH_CH2C = (CH + CH2C - CHCH2C);
end
if (CHCH2C_to_CH_CH2C <0)
        CHCH2C_to_CH_CH2C =0.001;
end

  CHCH2C_to_CH2CH_C = pp(124) + ddEA86 + CHCH2C_to_CH2CH_C_H_sloppe - CHCH2C_0;
if (CHCH2C_to_CH2CH_C < (CH2CH + C- CHCH2C))
	CHCH2C_to_CH2CH_C = (CH2CH + C- CHCH2C);
end
if (CHCH2C_to_CH2CH_C <0)
        CHCH2C_to_CH2CH_C =0.001;
end

  CHCH2C_to_CHCHC_H = pp(125) + ddEA87 + CHCH2C_to_CHCHC_H_H_sloppe - CHCH2C_0;
if (CHCH2C_to_CHCHC_H < (CHCHC + H -CHCH2C))
	CHCH2C_to_CHCHC_H = (CHCHC + H -CHCH2C);
end
if (CHCH2C_to_CHCHC_H <0)
        CHCH2C_to_CHCHC_H =0.001;
end

  CHCH2C_to_CCH2C_H = pp(126) + ddEA88 + CHCH2C_to_CCH2C_H_H_sloppe - CHCH2C_0;
if (CHCH2C_to_CCH2C_H < (CCH2C +H - CHCH2C))
	CHCH2C_to_CCH2C_H = (CCH2C +H - CHCH2C);
end
if (CHCH2C_to_CCH2C_H <0)
        CHCH2C_to_CCH2C_H = 0.001;
end

  CHCHCH_to_CHCH_CH = pp(127) + ddEA89 + CHCHCH_to_CH_CHCH_H_sloppe - CHCHCH_0;
if (CHCHCH_to_CHCH_CH < (CHCH+CH - CHCHCH))
	CHCHCH_to_CHCH_CH = CHCH + CH - CHCHCH;
end
if (CHCHCH_to_CHCH_CH <0)
        CHCHCH_to_CHCH_CH =0.001;
end

  CHCHCH_to_CHCHC_H = pp(128) + ddEA90 + CHCHCH_to_CHCHC_H_H_sloppe - CHCHCH_0;
if (CHCHCH_to_CHCHC_H < (CHCHC + H - CHCHCH))
	CHCHCH_to_CHCHC_H = (CHCHC + H - CHCHCH);
end
if (CHCHCH_to_CHCHC_H <0)
        CHCHCH_to_CHCHC_H =0.001;
end

  CHCHCH_to_CHCCH_H = pp(129) + ddEA91 + CHCHCH_to_CHCCH_H_H_sloppe - CHCHCH_0;
if (CHCHCH_to_CHCCH_H < (CHCCH + H - CHCHCH))
	CHCHCH_to_CHCCH_H = (CHCCH + H - CHCHCH);
end
if (CHCHCH_to_CHCCH_H <0)
        CHCHCH_to_CHCCH_H =0.001;
end

  CH2CCH_to_CH2_CCH = pp(130) + ddEA92 + CH2CCH_to_CH2_CCH_H_sloppe - CH2CCH_0;
if (CH2CCH_to_CH2_CCH < (CH2 + CHC - CH2CCH))
	CH2CCH_to_CH2_CCH = (CH2 + CHC - CH2CCH);
end
if (CH2CCH_to_CH2_CCH <0)
        CH2CCH_to_CH2_CCH =0.001;
end

  CH2CCH_to_CH2C_CH = pp(131) + ddEA93 + CH2CCH_to_CH2C_CH_H_sloppe - CH2CCH_0;
if (CH2CCH_to_CH2C_CH < (CH2C + CH - CH2CCH))
	CH2CCH_to_CH2C_CH = (CH2C + CH - CH2CCH);
end
if (CH2CCH_to_CH2C_CH <0)
        CH2CCH_to_CH2C_CH =0.001;
end

  CH2CCH_to_CH2CC_H = pp(132) + ddEA94 + CH2CCH_to_CH2CC_H_H_sloppe - CH2CCH_0;
if (CH2CCH_to_CH2CC_H < (CH2CC +H - CH2CCH))
	CH2CCH_to_CH2CC_H = (CH2CC +H - CH2CCH);
end
if (CH2CCH_to_CH2CC_H <0)
        CH2CCH_to_CH2CC_H =0.001;
end

  CH2CCH_to_CHCCH_H = pp(133) + ddEA95 + CH2CCH_to_CHCCH_H_H_sloppe - CH2CCH_0;
if (CH2CCH_to_CHCCH_H < (CHCCH + H - CH2CCH))
	CH2CCH_to_CHCCH_H = (CHCCH + H - CH2CCH);
end
if (CH2CCH_to_CHCCH_H <0)
        CH2CCH_to_CHCCH_H =0.001;
end

  CH2CC_to_CH2_CC = pp(134) + ddEA96 + CH2CC_to_CH2_CC_H_sloppe - CH2CC_0;
if (CH2CC_to_CH2_CC < (CH2 + CC - CH2CC))
	CH2CC_to_CH2_CC = (CH2 + CC - CH2CC);
end
if (CH2CC_to_CH2_CC <0)
        CH2CC_to_CH2_CC =0.001;
end

  CH2CC_to_CH2C_C = pp(135) + ddEA97 + CH2CC_to_CH2C_C_H_sloppe - CH2CC_0;
if (CH2CC_to_CH2C_C < (CH2C + C - CH2CC))
	CH2CC_to_CH2C_C = (CH2C + C - CH2CC);
end
if (CH2CC_to_CH2C_C <0)
        CH2CC_to_CH2C_C = 0.001;
end

  CH2CC_to_CHCC_H = pp(136) + ddEA98 + CH2CC_to_CHCC_H_H_sloppe - CH2CC_0;
if (CH2CC_to_CHCC_H < (CHCC + H - CH2CC))
	CH2CC_to_CHCC_H = (CHCC + H - CH2CC);	
end
if (CH2CC_to_CHCC_H <0)
        CH2CC_to_CHCC_H =0.001;
end

  CHCHC_to_CH_CHC = pp(137) + ddEA99 + CHCHC_to_CH_CHC_H_sloppe - CHCHC_0;
if (CHCHC_to_CH_CHC < (CHC + CH - CHCHC))
	CHCHC_to_CH_CHC = (CHC + CH - CHCHC);
end
if (CHCHC_to_CH_CHC <0)
        CHCHC_to_CH_CHC =0.001;
end

  CHCHC_to_CHCH_C = pp(138) + ddEA100 + CHCHC_to_CHCH_C_H_sloppe - CHCHC_0;
if (CHCHC_to_CHCH_C < (CHCH + C - CHCHC))
	CHCHC_to_CHCH_C = (CHCH + C - CHCHC);
end
if (CHCHC_to_CHCH_C <0)
        CHCHC_to_CHCH_C = 0.001;
end

  CHCHC_to_CHCC_H = pp(139) + ddEA101 + CHCHC_to_CHCC_H_H_sloppe - CHCHC_0;
if (CHCHC_to_CHCC_H < (CHCC + H - CHCHC))
	CHCHC_to_CHCC_H = (CHCC + H - CHCHC);
end
if (CHCHC_to_CHCC_H <0)
        CHCHC_to_CHCC_H =0.001;
end

  CHCHC_to_CCHC_H = pp(140) + ddEA102 + CHCHC_to_CCHC_H_H_sloppe - CHCHC_0;
if (CHCHC_to_CCHC_H < (CCHC + H - CHCHC))
	CHCHC_to_CCHC_H = (CCHC + H - CHCHC);
end
if (CHCHC_to_CCHC_H <0)
        CHCHC_to_CCHC_H =0.001;
end

  CCH2C_to_CH2C_C = pp(141) + ddEA103 + CCH2C_to_C_CH2C_H_sloppe - CCH2C_0;
if (CCH2C_to_CH2C_C < (CH2C + C - CCH2C))
	CCH2C_to_CH2C_C = (CH2C + C - CCH2C);
end
if (CCH2C_to_CH2C_C <0)
        CCH2C_to_CH2C_C =0.001;
end

  CCH2C_to_CCHC_H = pp(142) + ddEA104 + CCH2C_to_CCHC_H_H_sloppe - CCH2C_0;
if (CCH2C_to_CCHC_H < (CCHC + H - CCH2C))
	CCH2C_to_CCHC_H = (CCHC + H - CCH2C);
end
if (CCH2C_to_CCHC_H <0)
        CCH2C_to_CCHC_H =0.001;
end

  CHCCH_to_CHC_CH = pp(143) + ddEA105 + CHCCH_to_CH_CHC_H_sloppe - CHCCH_0;
if (CHCCH_to_CHC_CH < (CHC + CH - CHCCH))
	CHCCH_to_CHC_CH = (CHC + CH - CHCCH);
end
if (CHCCH_to_CHC_CH <0)
        CHCCH_to_CHC_CH =0.001;
end

  CHCCH_to_CHCC_H = pp(144) + ddEA106 + CHCCH_to_CHCC_H_H_sloppe - CHCCH_0;
if (CHCCH_to_CHCC_H < (CHCC + H - CHCCH))
	CHCCH_to_CHCC_H = (CHCC + H - CHCCH);
end
if (CHCCH_to_CHCC_H <0)
        CHCCH_to_CHCC_H =0.001;
end

  CCHC_to_CHC_C = pp(145) + ddEA107 + CCHC_to_CHC_C_H_sloppe - CCHC_0;
if (CCHC_to_CHC_C < (C + CHC - CCHC))
	CCHC_to_CHC_C = (C + CHC - CCHC);
end
if (CCHC_to_CHC_C <0)
        CCHC_to_CHC_C =0.001;
end

  CCHC_to_CCC_H = pp(146) + ddEA108 + CCHC_to_CCC_H_H_sloppe - CCHC_0;
if (CCHC_to_CCC_H < (CCC +H - CCHC))
	CCHC_to_CCC_H = (CCC +H - CCHC);
end
if (CCHC_to_CCC_H <0)
        CCHC_to_CCC_H =0.001;
end

  CHCC_to_CH_CC = pp(147) + ddEA109 + CHCC_to_CH_CC_H_sloppe - CHCC_0;
if (CHCC_to_CH_CC< (CC+ CH - CHCC))
	CHCC_to_CH_CC= (CC+ CH - CHCC);
end
if (CHCC_to_CH_CC <0)
        CHCC_to_CH_CC =0.001;
end

  CHCC_to_CHC_C = pp(148) + ddEA110 + CHCC_to_CHC_C_H_sloppe - CHCC_0;
if (CHCC_to_CHC_C < (CHC + C - CHCC))
	CHCC_to_CHC_C = (CHC + C - CHCC);
end
if (CHCC_to_CHC_C <0)
        CHCC_to_CHC_C =0.001;
end

  CHCC_to_CCC_H = pp(149) + ddEA111 + CHCC_to_CCC_H_H_sloppe - CHCC_0;
if (CHCC_to_CCC_H < (CCC + H - CHCC))
	CHCC_to_CCC_H = (CCC + H - CHCC);
end
if (CHCC_to_CCC_H <0)
        CHCC_to_CCC_H = 0.001;
end

  CCC_to_CC_C = pp(150) + ddEA112 + CCC_to_CC_C_H_sloppe - CCC_0;
if (CCC_to_CC_C < (C +CC - CCC))
	CCC_to_CC_C = (C +CC - CCC);
end
if (CCC_to_CC_C <0)
      CCC_to_CC_C = 0.001;
end

  CH4_to_CH3_H = pp(151) + ddEA113 + CH4_to_CH3_H_H_sloppe - CH4_0;
if (CH4_to_CH3_H < (H + CH3 - CH4))
	CH4_to_CH3_H = (H + CH3 - CH4);
end
if (CH4_to_CH3_H <0)
        CH4_to_CH3_H = 0.001;
end

  CH3_to_CH2_H = pp(152) + ddEA114 + CH3_to_CH2_H_H_sloppe - CH3_0;
if (CH3_to_CH2_H < (H+CH2 - CH3))
	CH3_to_CH2_H = (H+CH2 - CH3);
end
if (CH3_to_CH2_H <0)
        CH3_to_CH2_H =0.001;
end

  CH2_to_CH_H = pp(153) + ddEA115 + CH2_to_CH_H_H_sloppe - CH2_0;
if (CH2_to_CH_H <( CH + H - CH2))
	CH2_to_CH_H = (CH + H -CH2);
end
if (CH2_to_CH_H <0)
        CH2_to_CH_H =0.001;
end

  CH_to_C_H = pp(154) + ddEA116 + CH_to_C_H_H_sloppe - CH_0;
if (CH_to_C_H < (C+H - CH))
	CH_to_C_H = (C+H - CH);
end
if (CH_to_C_H <0)
        CH_to_C_H =0.001;
end

  CH3CH3_to_CH3_CH3 = pp(155) + ddEA117 + CH3CH3_to_CH3_CH3_H_sloppe - CH3CH3_0;
if (CH3CH3_to_CH3_CH3 < (CH3 + CH3 - CH3CH3))
	CH3CH3_to_CH3_CH3 = (CH3 + CH3 - CH3CH3);
end
if (CH3CH3_to_CH3_CH3 <0)
        CH3CH3_to_CH3_CH3 =0.001;
end

  CH3CH3_to_CH3CH2_H = pp(156) + ddEA118 + CH3CH3_to_CH3CH2_H_H_sloppe - CH3CH3_0;
if (CH3CH3_to_CH3CH2_H < (CH3CH2 + H - CH3CH3))
	CH3CH3_to_CH3CH2_H = (CH3CH2 + H - CH3CH3);
end
if (CH3CH3_to_CH3CH2_H <0)
        CH3CH3_to_CH3CH2_H =0.001;
end

  CH3CH2_to_CH3_CH2 = pp(157) + ddEA119 + CH3CH2_to_CH3_CH2_H_sloppe - CH3CH2_0;
if (CH3CH2_to_CH3_CH2 < (CH3 + CH2 - CH3CH2))
	CH3CH2_to_CH3_CH2 = (CH3 + CH2 - CH3CH2);
end
if (CH3CH2_to_CH3_CH2 <0)
        CH3CH2_to_CH3_CH2 =0.001;
end

  CH3CH2_to_CH3CH_H = pp(158) + ddEA120 + CH3CH2_to_CH3CH_H_H_sloppe - CH3CH2_0;
if (CH3CH2_to_CH3CH_H < (CH3CH + H -CH3CH2))
	CH3CH2_to_CH3CH_H = (CH3CH + H -CH3CH2);
end
if (CH3CH2_to_CH3CH_H <0)
        CH3CH2_to_CH3CH_H =0.001;
end

  CH3CH2_to_CH2CH2_H = pp(159) + ddEA121 + CH3CH2_to_CH2CH2_H_H_sloppe - CH3CH2_0;
if (CH3CH2_to_CH2CH2_H < (CH2CH2 + H - CH3CH2))
	CH3CH2_to_CH2CH2_H = (CH2CH2 + H - CH3CH2);
end
if (CH3CH2_to_CH2CH2_H <0)
        CH3CH2_to_CH2CH2_H =0.001;
end

  CH3CH_to_CH3_CH = pp(160) + ddEA122 + CH3CH_to_CH3_CH_H_sloppe - CH3CH_0;
if (CH3CH_to_CH3_CH < (CH3 + CH - CH3CH))
	CH3CH_to_CH3_CH = (CH3 + CH - CH3CH);
end
if (CH3CH_to_CH3_CH <0)
        CH3CH_to_CH3_CH =0.001;
end

  CH3CH_to_CH3C_H = pp(161) + ddEA123 + CH3CH_to_CH3C_H_H_sloppe - CH3CH_0;
if (CH3CH_to_CH3C_H < (CH3C + H - CH3CH))
	CH3CH_to_CH3C_H = (CH3C + H - CH3CH);
end
if (CH3CH_to_CH3C_H <0)
        CH3CH_to_CH3C_H =0.001;
end

  CH3CH_to_CH2CH_H = pp(162) + ddEA124 + CH3CH_to_CH2CH_H_H_sloppe - CH3CH_0;
if (CH3CH_to_CH2CH_H < (CH2CH +H - CH3CH))
	CH3CH_to_CH2CH_H = (CH2CH +H - CH3CH);
end
if (CH3CH_to_CH2CH_H<0)
 CH3CH_to_CH2CH_H =0.001;
end

  CH3C_to_CH3_C = pp(163) + ddEA125 + CH3C_to_CH3_C_H_sloppe - CH3C_0;
if (CH3C_to_CH3_C < (C+CH3 +CH3C))
	CH3C_to_CH3_C = (C+CH3 +CH3C);
end
if (CH3C_to_CH3_C<0)
        CH3C_to_CH3_C =0.001;
end

  CH3C_to_CH2C_H = pp(164) + ddEA126 + CH3C_to_CH2C_H_H_sloppe - CH3C_0;
if (CH3C_to_CH2C_H < (CH2C + H - CH3C))
	CH3C_to_CH2C_H = (CH2C + H - CH3C);
end
if (CH3C_to_CH2C_H <0)
        CH3C_to_CH2C_H =0.001;
end

  CH2CH2_to_CH2_CH2 = pp(165) + ddEA127 + CH2CH2_to_CH2_CH2_H_sloppe - CH2CH2_0;
if (CH2CH2_to_CH2_CH2 < (CH2 + CH2 - CH2CH2))
	CH2CH2_to_CH2_CH2 = (CH2 + CH2 - CH2CH2);
end
if (CH2CH2_to_CH2_CH2 <0)
        CH2CH2_to_CH2_CH2 =0.001;
end

  CH2CH2_to_CH2CH_H = pp(166) + ddEA128 + CH2CH2_to_CH2CH_H_H_sloppe - CH2CH2_0;
if (CH2CH2_to_CH2CH_H < (CH2CH+H - CH2CH2))
	CH2CH2_to_CH2CH_H = (CH2CH+H - CH2CH2);
end
if (CH2CH2_to_CH2CH_H <0)
        CH2CH2_to_CH2CH_H = 0.001;
end

  CH2CH_to_CH2_CH = pp(167) + ddEA129 + CH2CH_to_CH2_CH_H_sloppe - CH2CH_0;
if (CH2CH_to_CH2_CH < (CH2 + CH - CH2CH))
	CH2CH_to_CH2_CH = (CH2 + CH - CH2CH);
end
if (CH2CH_to_CH2_CH <0)
        CH2CH_to_CH2_CH =0.001;
end

  CH2CH_to_CH2C_H = pp(168) + ddEA130 + CH2CH_to_CH2C_H_H_sloppe - CH2CH_0;
if (CH2CH_to_CH2C_H < (CH2C + H - CH2CH))
	CH2CH_to_CH2C_H = (CH2C + H - CH2CH);
end
if (CH2CH_to_CH2C_H <0)
        CH2CH_to_CH2C_H =0.001;
end

  CH2CH_to_CHCH_H = pp(169) + ddEA131 + CH2CH_to_CHCH_H_H_sloppe - CH2CH_0;
if (CH2CH_to_CHCH_H < (CHCH +H - CH2CH))
	CH2CH_to_CHCH_H = (CHCH +H - CH2CH);
end
if (CH2CH_to_CHCH_H <0)
        CH2CH_to_CHCH_H =0.001;
end

  CH2C_to_CH2_C = pp(170) + ddEA132 + CH2C_to_CH2_C_H_sloppe - CH2C_0;
if (CH2C_to_CH2_C < (CH2+C - CH2C))
	CH2C_to_CH2_C = (CH2+C - CH2C);
end
if (CH2C_to_CH2_C <0)
        CH2C_to_CH2_C =0.001;
end

  CH2C_to_CHC_H = pp(171) + ddEA133 + CH2C_to_CHC_H_H_sloppe - CH2C_0;
if (CH2C_to_CHC_H < (CHC +H - CH2C))
	CH2C_to_CHC_H = (CHC +H - CH2C);
end
if (CH2C_to_CHC_H <0)
        CH2C_to_CHC_H =0.001;
end

  CHCH_to_CH_CH = pp(172) + ddEA134 + CHCH_to_CH_CH_H_sloppe - CHCH_0;
if (CHCH_to_CH_CH < (CH + CH - CHCH))
	CHCH_to_CH_CH = (CH + CH - CHCH);
end
if (CHCH_to_CH_CH <0)
        CHCH_to_CH_CH =0.001;
end

  CHCH_to_CHC_H = pp(173) + ddEA135 + CHCH_to_CHC_H_H_sloppe - CHCH_0;
if (CHCH_to_CHC_H < ( CHC + H - CHCH))
	CHCH_to_CHC_H = ( CHC + H - CHCH);
end
if (CHCH_to_CHC_H <0)
        CHCH_to_CHC_H =0.001;
end

  CHC_to_CH_C = pp(174) +ddEA136 + CHC_to_CH_C_H_sloppe - CHC_0;
if (CHC_to_CH_C < (CH + C - CHC))
	CHC_to_CH_C = (CH + C - CHC);
end
if (CHC_to_CH_C < 0)
        CHC_to_CH_C =0.001;
end

  CHC_to_CC_H = pp(175) + ddEA137 + CHC_to_CC_H_H_sloppe - CHC_0;
if (CHC_to_CC_H < (CC + H - CHC))
	CHC_to_CC_H = (CC + H - CHC);
end
if (CHC_to_CC_H <0)
        CHC_to_CC_H = 0.001;
end

  CC_to_C_C = pp(176)  + ddEA138 + CC_to_C_C_H_sloppe - CC_0;
if (CC_to_C_C < ( C+ C - CC))
	CC_to_C_C = ( C+ C - CC);
end
if (CC_to_C_C <0)
        CC_to_C_C = 0.001;
end


% calculating equilibrium rate constant for 5 adsorption rxns and 26 surface rxns

K = zeros(138,1); % Initialization

K(1)    = exp(-(CH3CH2CH3 - CH3CH2CH3_gp)/(T*8.6173324E-05));
K(2)    = exp(-(CH3CHCH2 - CH3CHCH2_gp)/(T*8.6173324E-05));
K(3)    = exp(-(2*H - H2_gp)/(T*8.6173324E-05));
K(4)    = exp(-(CH3CCH - CH3CCH_gp)/(T*8.6173324E-05));
K(5)    = exp(-(CH3CH3 - CH3CH3_gp)/(T*8.6173324E-05));
K(6)    = exp(-(CH2CH2 - CH2CH2_gp)/(T*8.6173324E-05));
K(7)    = exp(-(CHCH - CHCH_gp )/(T*8.6173324E-05));
K(8)    = exp(-(CH4 + CH4_gp)/(T*8.6173324E-05));

K(9)    = exp(-(CH3CHCH3 +H -CH3CH2CH3)/(T*8.6173324E-05));
K(10)   = exp(-(CH3CH2CH2 +H -CH3CH2CH3)/(T*8.6173324E-05));
K(11)   = exp(-(CH3CHCH2+H-CH3CHCH3)/(T*8.6173324E-05));
K(12)   = exp(-(CH3CHCH2+H -CH3CH2CH2)/(T*8.6173324E-05));
K(13)   = exp(-(CH3+CH3CH2 -CH3CH2CH3)/(T*8.6173324E-05));
K(14)   = exp(-(CH3+CH3CH -CH3CHCH3)/(T*8.6173324E-05));
K(15)   = exp(-(CH3CCH3 +H -CH3CHCH3)/(T*8.6173324E-05));
K(16)   = exp(-(CH3CH2+CH2- CH3CH2CH2)/(T*8.6173324E-05));
K(17)   = exp(-(CH3+CH2CH2 -CH3CH2CH2)/(T*8.6173324E-05));
K(18)   = exp(-(CH2CH2CH2+H -CH3CH2CH2)/(T*8.6173324E-05));
K(19)   = exp(-(CH3CH2CH +H -CH3CH2CH2)/(T*8.6173324E-05));
K(20)   = exp(-(CH3CH2+CH-CH3CH2CH)/(T*8.6173324E-05));
K(21)   = exp(-(CH3+CH2CH -CH3CH2CH)/(T*8.6173324E-05));
K(22)   = exp(-(CH3CH2C+H -CH3CH2CH)/(T*8.6173324E-05));
K(23)   = exp(-(CH3CHCH+H -CH3CH2CH)/(T*8.6173324E-05));

K(24)   = exp(-(CH2CH2CH+H -CH3CH2CH)/(T*8.6173324E-05));
K(25)   = exp(-(CH2+CH2CH2 -CH2CH2CH2)/(T*8.6173324E-05));
K(26)   = exp(-(CH2CH2CH+H -CH2CH2CH2)/(T*8.6173324E-05));
K(27)   = exp(-(CH2CHCH2+H -CH2CH2CH2)/(T*8.6173324E-05));
K(28)   = exp(-(CH3+CH2CH  -CH3CHCH2)/(T*8.6173324E-05));
K(29)   = exp(-(CH3CH+CH2  -CH3CHCH2)/(T*8.6173324E-05));
K(30)   = exp(-(CH3CCH2   +H  -CH3CHCH2)/(T*8.6173324E-05));
K(31)   = exp(-(CH3CHCH   +H  -CH3CHCH2)/(T*8.6173324E-05));
K(32)   = exp(-(CH2CHCH2   +H  -CH3CHCH2)/(T*8.6173324E-05));
K(33)   = exp(-(CH3C +CH3  -CH3CCH3)/(T*8.6173324E-05));
K(34)   = exp(-(CH3CCH2+H - CH3CCH3)/(T*8.6173324E-05));
K(35)   = exp(-(CH3 +CH2C - CH3CH2C )/(T*8.6173324E-05));
K(36)   = exp(-(CH3CH2 +C - CH3CH2C )/(T*8.6173324E-05));
K(37)   = exp(-(CH2CH2C + H - CH3CH2C )/(T*8.6173324E-05));
K(38)   = exp(-(CH3CHC +H - CH3CH2C )/(T*8.6173324E-05));
K(39)   = exp(-(CH2+CH2CH - CH2CH2CH )/(T*8.6173324E-05));
K(40)   = exp(-(CH2CH2+CH - CH2CH2CH )/(T*8.6173324E-05));
K(41)   = exp(-(CH2CH2C+H - CH2CH2CH )/(T*8.6173324E-05));

K(42)   = exp(-(CH2CHCH+H - CH2CH2CH )/(T*8.6173324E-05));
K(43)   = exp(-(CHCH2CH+H - CH2CH2CH )/(T*8.6173324E-05));
K(44)   = exp(-(CH2+CH2CH - CH2CHCH2 )/(T*8.6173324E-05));
K(45)   = exp(-(CH2CHCH+H - CH2CHCH2 )/(T*8.6173324E-05));
K(46)   = exp(-(CH2CCH2+H - CH2CHCH2 )/(T*8.6173324E-05));
K(47)   = exp(-(CH3+CHCH - CH3CHCH  )/(T*8.6173324E-05));
K(48)   = exp(-(CH3CH+CH - CH3CHCH  )/(T*8.6173324E-05));
K(49)   = exp(-(CH3CHC+H - CH3CHCH  )/(T*8.6173324E-05));
K(50)   = exp(-(CH3CCH+H - CH3CHCH  )/(T*8.6173324E-05));
K(51)   = exp(-(CH2CHCH+H - CH3CHCH  )/(T*8.6173324E-05));
K(52)   = exp(-(CH3+CH2C - CH3CCH2  )/(T*8.6173324E-05));
K(53)   = exp(-(CH3C+CH2 - CH3CCH2  )/(T*8.6173324E-05));
K(54)   = exp(-(CH2CCH2+H - CH3CCH2  )/(T*8.6173324E-05));
K(55)   = exp(-(CH3CCH+H - CH3CCH2  )/(T*8.6173324E-05));
K(56)   = exp(-(CH3+CHC - CH3CHC)/(T*8.617332E-05));
K(57)   = exp(-(CH3CH+C - CH3CHC)/(T*8.617332E-05));
K(58)   = exp(-(CH3CC+H - CH3CHC)/(T*8.617332E-05));
K(59)   = exp(-(CH2CHC+H - CH3CHC)/(T*8.617332E-05));
K(60)   = exp(-(CH2CH2+C - CH2CH2C)/(T*8.617332E-05));
K(61)   = exp(-(CH2+CH2C - CH2CH2C)/(T*8.617332E-05));

K(62)   = exp(-(CH2CHC+H - CH2CH2C)/(T*8.617332E-05));
K(63)   = exp(-(CHCH2C+H - CH2CH2C)/(T*8.617332E-05));
K(64)   = exp(-(CH2CH+CH - CHCH2CH)/(T*8.617332E-05));
K(65)   = exp(-(CHCH2C+H - CHCH2CH)/(T*8.617332E-05));
K(66)   = exp(-(CHCHCH+H - CHCH2CH)/(T*8.617332E-05));
K(67)   = exp(-(CH2+CHCH- CH2CHCH)/(T*8.617332E-05));
K(68)   = exp(-(CH2CH+CH - CH2CHCH)/(T*8.617332E-05));
K(69)   = exp(-(CH2CHC+H - CH2CHCH)/(T*8.617332E-05));
K(70)   = exp(-(CH2CCH+H - CH2CHCH)/(T*8.617332E-05));
K(71)   = exp(-(CHCHCH+H - CH2CHCH)/(T*8.617332E-05));
K(72)   = exp(-(CH2C+CH2 - CH2CCH2)/(T*8.617332E-05));
K(73)   = exp(-(CH2CCH+H - CH2CCH2)/(T*8.617332E-05));
K(74)   = exp(-(CH3C+CH - CH3CCH)/(T*8.617332E-05));
K(75)   = exp(-(CH3+CHC - CH3CCH)/(T*8.617332E-05));
K(76)   = exp(-(CH3CC+H - CH3CCH)/(T*8.617332E-05));
K(77)   = exp(-(CH2CCH+H - CH3CCH)/(T*8.617332E-05));
K(78)   = exp(-(CH3+CC - CH3CC)/(T*8.617332E-05));
K(79)   = exp(-(CH3C+C - CH3CC)/(T*8.617332E-05));
K(80)   = exp(-(CH2CC+H - CH3CC)/(T*8.617332E-05));
K(81)   = exp(-(CH2+CHC  -CH2CHC)/(T*8.617332E-05));
K(82)   = exp(-(CH2CH+C  -CH2CHC)/(T*8.617332E-05));

K(83)   = exp(-(CH2CC+H -CH2CHC)/(T*8.617332E-05));
K(84)   = exp(-(CHCHC+H  -CH2CHC)/(T*8.617332E-05));
K(85)   = exp(-(CH+CH2C -CHCH2C)/(T*8.617332E-05));
K(86)   = exp(-(CH2CH+C -CHCH2C)/(T*8.617332E-05));
K(87)   = exp(-(CHCHC+H -CHCH2C)/(T*8.617332E-05));
K(88)   = exp(-(CCH2C+H -CHCH2C)/(T*8.617332E-05));
K(89)   = exp(-(CHCH+CH -CHCHCH)/(T*8.617332E-05));
K(90)   = exp(-(CHCHC+H -CHCHCH)/(T*8.617332E-05));
K(91)   = exp(-(CHCCH+H -CHCHCH)/(T*8.617332E-05));
K(92)   = exp(-(CH2+CHC -CH2CCH)/(T*8.617332E-05));
K(93)   = exp(-(CH2C+CH -CH2CCH)/(T*8.617332E-05));
K(94)   = exp(-(CH2CC+H -CH2CCH)/(T*8.617332E-05));
K(95)   = exp(-(CHCCH+H -CH2CCH)/(T*8.617332E-05));
K(96)   = exp(-(CH2+CC -CH2CC)/(T*8.617332E-05));
K(97)   = exp(-(CH2C+C -CH2CC)/(T*8.617332E-05));
K(98)   = exp(-(CHCC+H -CH2CC)/(T*8.617332E-05));
K(99)   = exp(-(CH+CHC -CHCHC)/(T*8.617332E-05));
K(100)  = exp(-(CHCH+C -CHCHC)/(T*8.617332E-05));
K(101)  = exp(-(CHCC+H -CHCHC)/(T*8.617332E-05));
K(102)  = exp(-(CCHC+H -CHCHC)/(T*8.617332E-05));
K(103)  = exp(-(CH2C+C -CCH2C)/(T*8.617332E-05));
K(104)  = exp(-(CCHC+H -CCH2C)/(T*8.617332E-05));
K(105)  = exp(-(CH+CHC -CHCCH)/(T*8.617332E-05));
K(106)  = exp(-(CHCC+H -CHCCH)/(T*8.617332E-05));
K(107)  = exp(-(CHC+C -CCHC)/(T*8.617332E-05));

K(108)  = exp(-(CCC+H-CCHC)/(T*8.617332E-05));
K(109)  = exp(-(CH+CC-CHCC)/(T*8.617332E-05));
K(110)  = exp(-(CHC+C-CHCC)/(T*8.617332E-05));
K(111)  = exp(-(CCC+H-CHCC)/(T*8.617332E-05));
K(112)  = exp(-(C+CC-CCC)/(T*8.617332E-05));
K(113)  = exp(-(CH3+H-CH4)/(T*8.617332E-05));
K(114)  = exp(-(CH2+H-CH3)/(T*8.617332E-05));
K(115)  = exp(-(CH+H-CH2)/(T*8.617332E-05));
K(116)  = exp(-(C+H-CH)/(T*8.617332E-05));
K(117)  = exp(-(CH3+CH3 -CH3CH3)/(T*8.617332E-05));
K(118)  = exp(-(CH3CH2+H -CH3CH3)/(T*8.617332E-05));
K(119)  = exp(-(CH3+CH2 -CH3CH2)/(T*8.617332E-05));
K(120)  = exp(-(CH3CH+H -CH3CH2)/(T*8.617332E-05));
K(121)  = exp(-(CH2CH2+H -CH3CH2)/(T*8.617332E-05));
K(122)  = exp(-(CH3+CH -CH3CH )/(T*8.617332E-05));
K(123)  = exp(-(CH3C+H -CH3CH )/(T*8.617332E-05));
K(124)  = exp(-(CH2CH+H -CH3CH )/(T*8.617332E-05));
K(125)  = exp(-(CH3+C -CH3C  )/(T*8.617332E-05));
K(126)  = exp(-(CH2C+H -CH3C  )/(T*8.617332E-05));
K(127)  = exp(-(CH2+CH2 -CH2CH2  )/(T*8.617332E-05));
K(128)  = exp(-(CH2CH+H -CH2CH2  )/(T*8.617332E-05));
K(129)  = exp(-(CH2+CH -CH2CH  )/(T*8.617332E-05));
K(130)  = exp(-(CH2C+H -CH2CH  )/(T*8.617332E-05));
K(131)  = exp(-(CHCH+H -CH2CH  )/(T*8.617332E-05));
K(132)  = exp(-(CH2+C -CH2C  )/(T*8.617332E-05));
K(133)  = exp(-(CHC+H -CH2C  )/(T*8.617332E-05));
K(134)  = exp(-(CH+CH -CHCH  )/(T*8.617332E-05));
K(135)  = exp(-(CHC+H -CHCH  )/(T*8.617332E-05));
K(136)  = exp(-(CH+C-CHC   )/(T*8.617332E-05));
K(137)  = exp(-(CC+H-CHC   )/(T*8.617332E-05));
K(138)  = exp(-(C+C-CC)/(T*8.617332E-05));

% Calculating forward rates for 26 surface rxns

kB = 8.617385e-5;
hh = 4.135667516E-15;

f(9)  = (kB*T/hh)*exp(-(CH3CH2CH3_to_CH3CHCH3_H)/(kB*T));
f(10) = (kB*T/hh)*exp(-(CH3CH2CH3_to_CH3CH2CH2_H)/(kB*T));
f(11) = (kB*T/hh)*exp(-(CH3CHCH3_to_CH3CHCH2_H/(kB*T)));
f(12) = (kB*T/hh)*exp(-(CH3CH2CH2_to_CH3CHCH2_H)/(kB*T));
f(13) = (kB*T/hh)*exp(-(CH3CH2CH3_to_CH3CH2_CH3)/(kB*T));
f(14) = (kB*T/hh)*exp(-(CH3CHCH3_to_CH3_CH3CH  )/(kB*T));
f(15) = (kB*T/hh)*exp(-(CH3CHCH3_to_CH3CCH3_H  )/(kB*T));
f(16) = (kB*T/hh)*exp(-(CH3CH2CH2_to_CH3CH2_CH2)/(kB*T));
f(17) = (kB*T/hh)*exp(-(CH3CH2CH2_to_CH3_CH2CH2)/(kB*T));
f(18) = (kB*T/hh)*exp(-(CH3CH2CH2_to_CH2CH2CH2_H)/(kB*T));
f(19) = (kB*T/hh)*exp(-(CH3CH2CH2_to_CH3CH2CH_H)/(kB*T));

f(20) = (kB*T/hh)*exp(-(CH3CH2CH_to_CH3CH2_CH)/(kB*T));
f(21) = (kB*T/hh)*exp(-(CH3CH2CH_to_CH3_CH2CH)/(kB*T));
f(22) = (kB*T/hh)*exp(-(CH3CH2CH_to_CH3CH2C_H)/(kB*T));
f(23) = (kB*T/hh)*exp(-(CH3CH2CH_to_CH3CHCH_H)/(kB*T));
f(24) = (kB*T/hh)*exp(-(CH3CH2CH_to_CH2CH2CH_H)/(kB*T));
f(25) = (kB*T/hh)*exp(-(CH2CH2CH2_to_CH2CH2_CH2)/(kB*T));
f(26) = (kB*T/hh)*exp(-(CH2CH2CH2_to_CH2CH2CH_H)/(kB*T));
f(27) = (kB*T/hh)*exp(-(CH2CH2CH2_to_CH2CHCH2_H)/(kB*T));
f(28) = (kB*T/hh)*exp(-(CH3CHCH2_to_CH3_CHCH2)/(kB*T));
f(29) = (kB*T/hh)*exp(-(CH3CHCH2_to_CH3CH_CH2)/(kB*T));


f(30) = (kB*T/hh)*exp(-(CH3CHCH2_to_CH3CCH2_H)/(kB*T));
f(31) = (kB*T/hh)*exp(-(CH3CHCH2_to_CH3CHCH_H)/(kB*T));
f(32) = (kB*T/hh)*exp(-(CH3CHCH2_to_CH2CHCH2_H)/(kB*T));
f(33) = (kB*T/hh)*exp(-(CH3CCH3_to_CH3_CCH3   )/(kB*T));
f(34) = (kB*T/hh)*exp(-(CH3CCH3_to_CH3CCH2_H  )/(kB*T));
f(35) = (kB*T/hh)*exp(-(CH3CH2C_to_CH3_CH2C   )/(kB*T));
f(36) = (kB*T/hh)*exp(-(CH3CH2C_to_CH3CH2_C   )/(kB*T));
f(37) = (kB*T/hh)*exp(-(CH3CH2C_to_CH2CH2C_H   )/(kB*T));
f(38) = (kB*T/hh)*exp(-(CH3CH2C_to_CH3CHC_H   )/(kB*T));
f(39) = (kB*T/hh)*exp(-(CH2CH2CH_to_CH2_CH2CH )/(kB*T));

f(40) = (kB*T/hh)*exp(-(CH2CH2CH_to_CH2CH2_CH )/(kB*T));
f(41) = (kB*T/hh)*exp(-(CH2CH2CH_to_CH2CH2C_H )/(kB*T));
f(42) = (kB*T/hh)*exp(-(CH2CH2CH_to_CH2CHCH_H )/(kB*T));
f(43) = (kB*T/hh)*exp(-(CH2CH2CH_to_CHCH2CH_H )/(kB*T));
f(44) = (kB*T/hh)*exp(-(CH2CHCH2_to_CH2_CH2CH )/(kB*T));
f(45) = (kB*T/hh)*exp(-(CH2CHCH2_to_CH2CHCH_H )/(kB*T));
f(46) = (kB*T/hh)*exp(-(CH2CHCH2_to_CH2CCH2_H )/(kB*T));
f(47) = (kB*T/hh)*exp(-(CH3CHCH_to_CH3_CHCH   )/(kB*T));
f(48) = (kB*T/hh)*exp(-(CH3CHCH_to_CH3CH_CH   )/(kB*T));
f(49) = (kB*T/hh)*exp(-(CH3CHCH_to_CH3CHC_H   )/(kB*T));

f(50) = (kB*T/hh)*exp(-(CH3CHCH_to_CH3CCH_H   )/(kB*T));
f(51) = (kB*T/hh)*exp(-(CH3CHCH_to_CH2CHCH_H   )/(kB*T));
f(52) = (kB*T/hh)*exp(-(CH3CCH2_to_CH3_CCH2   )/(kB*T));
f(53) = (kB*T/hh)*exp(-(CH3CCH2_to_CH3C_CH2   )/(kB*T));
f(54) = (kB*T/hh)*exp(-(CH3CCH2_to_CH2CCH2_H  )/(kB*T));
f(55) = (kB*T/hh)*exp(-(CH3CCH2_to_CH3CCH_H   )/(kB*T));
f(56) = (kB*T/hh)*exp(-(CH3CHC_to_CH3_CHC     )/(kB*T));
f(57) = (kB*T/hh)*exp(-(CH3CHC_to_CH3CH_C     )/(kB*T));
f(58) = (kB*T/hh)*exp(-(CH3CHC_to_CH3CC_H     )/(kB*T));
f(59) = (kB*T/hh)*exp(-(CH3CHC_to_CH2CHC_H    )/(kB*T));

f(60) = (kB*T/hh)*exp(-(CH2CH2C_to_CH2CH2_C    )/(kB*T));
f(61) = (kB*T/hh)*exp(-(CH2CH2C_to_CH2_CH2C    )/(kB*T));
f(62) = (kB*T/hh)*exp(-(CH2CH2C_to_CH2CHC_H    )/(kB*T));
f(63) = (kB*T/hh)*exp(-(CH2CH2C_to_CHCH2C_H    )/(kB*T));
f(64) = (kB*T/hh)*exp(-(CHCH2CH_to_CHCH2_CH    )/(kB*T));
f(65) = (kB*T/hh)*exp(-(CHCH2CH_to_CHCH2C_H    )/(kB*T));
f(66) = (kB*T/hh)*exp(-(CHCH2CH_to_CHCHCH_H    )/(kB*T));
f(67) = (kB*T/hh)*exp(-(CH2CHCH_to_CH2_CHCH    )/(kB*T));
f(68) = (kB*T/hh)*exp(-(CH2CHCH_to_CH2CH_CH    )/(kB*T));
f(69) = (kB*T/hh)*exp(-(CH2CHCH_to_CH2CHC_H    )/(kB*T));

f(70) = (kB*T/hh)*exp(-(CH2CHCH_to_CH2CCH_H    )/(kB*T));
f(71) = (kB*T/hh)*exp(-(CH2CHCH_to_CHCHCH_H    )/(kB*T));
f(72) = (kB*T/hh)*exp(-(CH2CCH2_to_CH2C_CH2    )/(kB*T));
f(73) = (kB*T/hh)*exp(-(CH2CCH2_to_CH2CCH_H    )/(kB*T));
f(74) = (kB*T/hh)*exp(-(CH3CCH_to_CH3C_CH      )/(kB*T));
f(75) = (kB*T/hh)*exp(-(CH3CCH_to_CH3_CHC      )/(kB*T));
f(76) = (kB*T/hh)*exp(-(CH3CCH_to_CH3CC_H      )/(kB*T));
f(77) = (kB*T/hh)*exp(-(CH3CCH_to_CH2CCH_H     )/(kB*T));
f(78) = (kB*T/hh)*exp(-(CH3CC_to_CH3_CC        )/(kB*T));
f(79) = (kB*T/hh)*exp(-(CH3CC_to_CH3C_C        )/(kB*T));

f(80) = (kB*T/hh)*exp(-(CH3CC_to_CH2CC_H        )/(kB*T));
f(81) = (kB*T/hh)*exp(-(CH2CHC_to_CH2_CHC       )/(kB*T));
f(82) = (kB*T/hh)*exp(-(CH2CHC_to_CH2CH_C       )/(kB*T));
f(83) = (kB*T/hh)*exp(-(CH2CHC_to_CH2CC_H       )/(kB*T));
f(84) = (kB*T/hh)*exp(-(CH2CHC_to_CHCHC_H       )/(kB*T));
f(85) = (kB*T/hh)*exp(-(CHCH2C_to_CH_CH2C       )/(kB*T));
f(86) = (kB*T/hh)*exp(-(CHCH2C_to_CH2CH_C       )/(kB*T));
f(87) = (kB*T/hh)*exp(-(CHCH2C_to_CHCHC_H       )/(kB*T));
f(88) = (kB*T/hh)*exp(-(CHCH2C_to_CCH2C_H       )/(kB*T));
f(89) = (kB*T/hh)*exp(-(CHCHCH_to_CHCH_CH       )/(kB*T));

f(90) = (kB*T/hh)*exp(-(CHCHCH_to_CHCHC_H       )/(kB*T));
f(91) = (kB*T/hh)*exp(-(CHCHCH_to_CHCCH_H       )/(kB*T));
f(92) = (kB*T/hh)*exp(-(CH2CCH_to_CH2_CCH       )/(kB*T));
f(93) = (kB*T/hh)*exp(-(CH2CCH_to_CH2C_CH       )/(kB*T));
f(94) = (kB*T/hh)*exp(-(CH2CCH_to_CH2CC_H       )/(kB*T));
f(95) = (kB*T/hh)*exp(-(CH2CCH_to_CHCCH_H       )/(kB*T));
f(96) = (kB*T/hh)*exp(-(CH2CC_to_CH2_CC   )/(kB*T));
f(97) = (kB*T/hh)*exp(-(CH2CC_to_CH2C_C   )/(kB*T));
f(98) = (kB*T/hh)*exp(-(CH2CC_to_CHCC_H   )/(kB*T));
f(99) = (kB*T/hh)*exp(-(CHCHC_to_CH_CHC   )/(kB*T));

f(100) = (kB*T/hh)*exp(-(CHCHC_to_CHCH_C   )/(kB*T));
f(101) = (kB*T/hh)*exp(-(CHCHC_to_CHCC_H   )/(kB*T));
f(102) = (kB*T/hh)*exp(-(CHCHC_to_CCHC_H   )/(kB*T));
f(103) = (kB*T/hh)*exp(-(CCH2C_to_CH2C_C   )/(kB*T));
f(104) = (kB*T/hh)*exp(-(CCH2C_to_CCHC_H   )/(kB*T));
f(105) = (kB*T/hh)*exp(-(CHCCH_to_CHC_CH   )/(kB*T));
f(106) = (kB*T/hh)*exp(-(CHCCH_to_CHCC_H   )/(kB*T));
f(107) = (kB*T/hh)*exp(-(CCHC_to_CHC_C     )/(kB*T));
f(108) = (kB*T/hh)*exp(-(CCHC_to_CCC_H     )/(kB*T));
f(109) = (kB*T/hh)*exp(-(CHCC_to_CH_CC     )/(kB*T));
f(110) = (kB*T/hh)*exp(-(CHCC_to_CHC_C     )/(kB*T));
f(111) = (kB*T/hh)*exp(-(CHCC_to_CCC_H     )/(kB*T));
f(112) = (kB*T/hh)*exp(-(CCC_to_CC_C       )/(kB*T));

f(113) = (kB*T/hh)*exp(-(CH4_to_CH3_H       )/(kB*T));
f(114) = (kB*T/hh)*exp(-(CH3_to_CH2_H       )/(kB*T));
f(115) = (kB*T/hh)*exp(-(CH2_to_CH_H       )/(kB*T));
f(116) = (kB*T/hh)*exp(-(CH_to_C_H       )/(kB*T));
f(117) = (kB*T/hh)*exp(-(CH3CH3_to_CH3_CH3      )/(kB*T));
f(118) = (kB*T/hh)*exp(-(CH3CH3_to_CH3CH2_H    )/(kB*T));
f(119) = (kB*T/hh)*exp(-(CH3CH2_to_CH3_CH2     )/(kB*T));

f(120) = (kB*T/hh)*exp(-(CH3CH2_to_CH3CH_H     )/(kB*T));
f(121) = (kB*T/hh)*exp(-(CH3CH2_to_CH2CH2_H     )/(kB*T));
f(122) = (kB*T/hh)*exp(-(CH3CH_to_CH3_CH        )/(kB*T));
f(123) = (kB*T/hh)*exp(-(CH3CH_to_CH3C_H        )/(kB*T));
f(124) = (kB*T/hh)*exp(-(CH3CH_to_CH2CH_H       )/(kB*T));
f(125) = (kB*T/hh)*exp(-(CH3C_to_CH3_C          )/(kB*T));
f(126) = (kB*T/hh)*exp(-(CH3C_to_CH2C_H         )/(kB*T));
f(127) = (kB*T/hh)*exp(-(CH2CH2_to_CH2_CH2      )/(kB*T));
f(128) = (kB*T/hh)*exp(-(CH2CH2_to_CH2CH_H      )/(kB*T));
f(129) = (kB*T/hh)*exp(-(CH2CH_to_CH2_CH        )/(kB*T));

f(130) = (kB*T/hh)*exp(-(CH2CH_to_CH2C_H        )/(kB*T));
f(131) = (kB*T/hh)*exp(-(CH2CH_to_CHCH_H        )/(kB*T));
f(132) = (kB*T/hh)*exp(-(CH2C_to_CH2_C          )/(kB*T));
f(133) = (kB*T/hh)*exp(-(CH2C_to_CHC_H          )/(kB*T));
f(134) = (kB*T/hh)*exp(-(CHCH_to_CH_CH          )/(kB*T));
f(135) = (kB*T/hh)*exp(-(CHCH_to_CHC_H          )/(kB*T));
f(136) = (kB*T/hh)*exp(-(CHC_to_CH_C            )/(kB*T));
f(137) = (kB*T/hh)*exp(-(CHC_to_CC_H            )/(kB*T));
f(138) = (kB*T/hh)*exp(-(CC_to_C_C              )/(kB*T));

% calculating backward rates for 138 reactions
b = zeros(138,1);

for i = 1:138
	b(i) = f(i)/K(i);
end	

% calculating rates for all 138 rxns

% r = zeros(138,1);

r1   = f(1)*PCH3CH2CH3*y(1) 				- b(1)*y(3);      		%Propane Adsorption           
r2   = f(2)*PCH3CHCH2*y(1) 					- b(2)*y(6);  			%Propylene Adsorption          
r3   = f(3)*PH2*y(1)*y(1) 					- b(3)*y(2)*y(2); 		%H2 Adsorption         
r4   = f(4)*PCH3CCH*y(1)*y(1)				- b(4)*y(20);                                                     
r5   = f(5)*PCH3CH3*y(1) 					- b(5)*y(33);                                 
r6   = f(6)*PCH2CH2*y(1)	 				- b(6)*y(37);                          
r7   = f(7)*PCHCH*y(1)*y(1)			 		- b(7)*y(40);                         
r8   = f(8)*PCH4*y(1) 						- b(8)*y(43);                               
r9 	 = f(9)*y(3)*y(1) 						- b(9)*y(4)*y(2);       % 1 + 1 - 1 - 1                                                                 
r10  = f(10)*y(3)*y(1) 						- b(10)*y(5)*y(2);		% 1 + 1 -1 -1
r11  = f(11)*y(4)*y(1) 						- b(11)*y(6)*y(2);		% 1 + 1 - 1 - 1 
r12  = f(12)*y(5)*y(1) 						- b(12)*y(6)*y(2);		% 1 + 1 -1 -1
r13  = f(13)*y(3)*y(1) 						- b(13)*y(34)*y(44);	% 1 + 1 -1 -1
r14  = f(14)*y(4)*y(1)*y(1) 				- b(14)*y(35)*y(44);	% 1 + 2 -1 - 2
r15  = f(15)*y(4)*y(1) 						- b(15)*y(9)*y(2);		% 1 + 1 -1 -1
r16  = f(16)*y(5)*y(1)*y(1) 				- b(16)*y(34)*y(45);	% 1 + 2 -1 -2
r17  = f(17)*y(5)*y(1)		 				- b(17)*y(37)*y(44);	% 1 + 1 - 1 -1
r18  = f(18)*y(5)*y(1)*y(1) 				- b(18)*y(8)*y(2);		% 2 + 1 -1 -2
r19  = f(19)*y(5)*y(1)*y(1) 				- b(19)*y(7)*y(2);		% 2 + 1 -1 -2
r20  = f(20)*y(7)*y(1)*y(1) 					- b(20)*y(34)*y(46);			% 3 + 1 -2 -2
r21  = f(21)*y(7) 							- b(21)*y(38)*y(44);	% 1 + 1 -2
r22  = f(22)*y(7)*y(1)*y(1) 					- b(22)*y(10)*y(2);			% 3 + 1 -2 -2
r23  = f(23)*y(7)	 						- b(23)*y(13)*y(2);			% 1 + 1 -2 
r24  = f(24)*y(7)*y(1)*y(1) 					- b(24)*y(11)*y(2);		% 3 + 1 - 2 -2
r25  = f(25)*y(8)*y(1) 						- b(25)*y(37)*y(45);	% 1 + 2 -2 -1
r26  = f(26)*y(8)*y(1)*y(1) 					- b(26)*y(11)*y(2);		% 3 + 1 -2 -2
r27  = f(27)*y(8)*y(1) 						- b(27)*y(12)*y(2);		% 2 + 1 -2 -1
r28  = f(28)*y(6)*y(1)	    				- b(28)*y(38)*y(44);	% 1 + 1 -1 -1
r29  = f(29)*y(6)*y(1)*y(1)*y(1) 			- b(29)*y(35)*y(45);	% 2 + 2 -1 -3
r30  = f(30)*y(6)*y(1) 						- b(30)*y(14)*y(2);		% 1 + 1 -1 -1
r31  = f(31)*y(6)*y(1) 						- b(31)*y(13)*y(2);			% 1 + 1 -1 -1
r32  = f(32)*y(6)*y(1)*y(1) 					- b(32)*y(12)*y(2);		% 2 + 1 -1 -2
r33  = f(33)*y(9)*y(1)*y(1)*y(1)			 	- b(33)*y(36)*y(44);	% 1 + 3 -1 -3
r34  = f(34)*y(9)*y(1)			 			- b(34)*y(14)*y(2);		% 1 + 1 -1 -1
r35  = f(35)*y(10)*y(1) 						- b(35)*y(39)*y(44); 	% 1 + 3 -3 -1
r36  = f(36)*y(10)*y(1) 						- b(36)*y(34)*y(47);	% 3 + 1 - 3 - 1
r37  = f(37)*y(10)*y(1)    					- b(37)*y(16)*y(2); 	% 3 + 1 -3 -1
r38  = f(38)*y(10) 							- b(38)*y(15)*y(2);		% 2 + 1 -3  
r39  = f(39)*y(11) 							- b(39)*y(38)*y(45);	% 1 + 2 -3
r40  = f(40)*y(11)*y(1) 						- b(40)*y(37)*y(46);		% 1 + 3 -3 -1 
r41  = f(41)*y(11)*y(1) 						- b(41)*y(16)*y(2); 	% 3 + 1 -3 -1
r42  = f(42)*y(11)   						- b(42)*y(18)*y(2);		% 2 + 1 - 3
r43  = f(43)*y(11)*y(1)*y(1)   				- b(43)*y(17)*y(2);			% 4 + 1 -3 -2
r44  = f(44)*y(12)*y(1)   					- b(44)*y(38)*y(45);	% 1 + 2 -2 -1
r45  = f(45)*y(12)*y(1)  					- b(45)*y(18)*y(2);			% 2 + 1 -2 -1
r46  = f(46)*y(12)*y(1)*y(1)  				- b(46)*y(19)*y(2);			% 3 + 1 -2 -2
r47  = f(47)*y(13)*y(1)*y(1)    				- b(47)*y(40)*y(44);		% 2 + 1 -1 -2
r48  = f(48)*y(13)*y(1)*y(1)*y(1)*y(1)  		- b(48)*y(35)*y(46);		% 2 + 3 -1 -4
r49  = f(49)*y(13)*y(1)*y(1)    				- b(49)*y(15)*y(2);			% 2 + 1 -1 -2
r50  = f(50)*y(13)*y(1)*y(1) 				- b(50)*y(20)*y(2);		% 2 + 1 -1 -2
r51  = f(51)*y(13)*y(1)*y(1) 				- b(51)*y(18)*y(2); 	% 2 + 1 -1 -2
r52  = f(52)*y(14)*y(1)*y(1)*y(1)			- b(52)*y(39)*y(44);	% 1 + 3 -1 -3
r53  = f(53)*y(14)*y(1)*y(1)*y(1)*y(1) 		- b(53)*y(36)*y(45);	% 2 + 3 -1 -4
r54  = f(54)*y(14)*y(1)*y(1)*y(1)			- b(54)*y(19)*y(2);		% 3 + 1 -1 -3 
r55  = f(55)*y(14)*y(1)*y(1)   				- b(55)*y(20)*y(2);		% 2 + 1 -1 -2
r56  = f(56)*y(15)			   				- b(56)*y(41)*y(44); 	% 1 + 1 -2
r57  = f(57)*y(15)*y(1)*y(1)*y(1)  			- b(57)*y(35)*y(47); 	% 2 + 3 -2 -3
r58  = f(58)*y(15)		 					- b(58)*y(21)*y(2);		% 1 + 1 -2
r59  = f(59)*y(15)*y(1)*y(1)   				- b(59)*y(22)*y(2); 	% 3 + 1 -2 -2
r60  = f(60)*y(16)*y(1)   					- b(60)*y(37)*y(47);	% 1 + 3 -3 -1 
r61  = f(61)*y(16)*y(1)*y(1)     			- b(61)*y(39)*y(45);		% 2 + 3 -3 -2
r62  = f(62)*y(16)*y(1)  					- b(62)*y(22)*y(2);			% 3 + 1 -3 -1
r63  = f(63)*y(16)*y(1)*y(1)   				- b(63)*y(24)*y(2);			% 4 + 1 - 3 -2
r64  = f(64)*y(17)			   				- b(64)*y(38)*y(46);		% 1 + 3 -4 
r65  = f(65)*y(17)*y(1)     					- b(65)*y(24)*y(2); 		% 4 + 1 - 4 -1
r66  = f(66)*y(17)     						- b(66)*y(23)*y(2)*y(1); 		% 2 + 1 + 1 -4
r67  = f(67)*y(18)*y(1)*y(1)		    		- b(67)*y(40)*y(45);		% 2 + 2 -2 -2
r68  = f(68)*y(18)*y(1)*y(1)	    			- b(68)*y(38)*y(46);		% 1 + 3 -2 -2
r69  = f(69)*y(18)*y(1)*y(1)     			- b(69)*y(22)*y(2);			% 3 + 1 -2 -2 
r70  = f(70)*y(18)*y(1)    					- b(70)*y(25)*y(2);		% 2 + 1 -2 -1
r71  = f(71)*y(18)*y(1)		   				- b(71)*y(23)*y(2);		% 2 + 1 -2 -1
r72  = f(72)*y(19)*y(1)*y(1)   				- b(72)*y(39)*y(45);	% 3 + 2 -3 -2
r73  = f(73)*y(19)	     					- b(73)*y(25)*y(2);		% 2 + 1 -3  
r74  = f(74)*y(20)*y(1)*y(1)*y(1)*y(1)   	- b(74)*y(36)*y(46);	% 3 + 3 -2 -4
r75  = f(75)*y(20)	    					- b(75)*y(41)*y(44);	% 1 + 1 -2  
r76  = f(76)*y(20)	    					- b(76)*y(21)*y(2);		% 1 + 1 -2  
r77  = f(77)*y(20)*y(1)        				- b(77)*y(25)*y(2);		% 2 + 1 -2 -1
r78  = f(78)*y(21)*y(1)*y(1)        			- b(78)*y(42)*y(44);	% 2 + 1 -1 -2
r79  = f(79)*y(21)*y(1)*y(1)*y(1)*y(1)*y(1)  - b(79)*y(36)*y(47); 		% 3 + 3 -1 -5
r80  = f(80)*y(21)*y(1)*y(1)     			- b(80)*y(26)*y(2);		% 2 + 1 -1 -2
r81  = f(81)*y(22)				        	- b(81)*y(41)*y(45);	% 1 + 2 -3
r82  = f(82)*y(22)*y(1)			        	- b(82)*y(38)*y(47);	% 1 + 3 -3 -1
r83  = f(83)*y(22)			     			- b(83)*y(26)*y(2);		% 2 + 1 -3 
r84  = f(84)*y(22)			     			- b(84)*y(27)*y(2);		% 2 + 1 -3 
r85  = f(85)*y(24)*y(1)*y(1)		        	- b(85)*y(39)*y(46);	% 3 + 3 -4 -2
r86  = f(86)*y(24)			       			- b(86)*y(38)*y(47);	% 1 + 3 -4
r87  = f(87)*y(24)	     					- b(87)*y(27)*y(2)*y(1); 	% 2 + 1 + 1 -4
r88  = f(88)*y(24)*y(1)	   					- b(88)*y(28)*y(2);		% 4 + 1 -4 -1
r89  = f(89)*y(23)*y(1)*y(1)*y(1)	    	- b(89)*y(40)*y(46);		% 2 + 3 -2 -3 
r90  = f(90)*y(23)*y(1)     					- b(90)*y(27)*y(2);		% 2 + 1 -2 -1
r91  = f(91)*y(23)*y(1)*y(1)*y(1) 			- b(91)*y(29)*y(2);		% 4 + 1 -2 -3
r92  = f(92)*y(25)*y(1)						- b(92)*y(41)*y(45); 	% 1 + 2 -2 -1
r93  = f(93)*y(25)*y(1)*y(1)*y(1)*y(1)		- b(93)*y(39)*y(46);	% 3 + 3 -2 -4 
r94  = f(94)*y(25)*y(1)		     			- b(94)*y(26)*y(2);		% 2 + 1 -2 -1
r95  = f(95)*y(25)*y(1)*y(1)*y(1)			- b(95)*y(29)*y(2);		% 4 + 1 -2 -3
r96  = f(96)*y(26)*y(1)*y(1)        			- b(96)*y(42)*y(45);	% 2 + 2 -2 -2
r97  = f(97)*y(26)*y(1)*y(1)*y(1)*y(1)   	- b(97)*y(39)*y(47); 	% 3 + 3 -2 -4
r98  = f(98)*y(26)*y(1)   					- b(98)*y(31)*y(2);		% 2 + 1 -2 -1
r99  = f(99)*y(27)*y(1)*y(1)			    	- b(99)*y(41)*y(46);	% 1 + 3 -2 -2
r100 = f(100)*y(27)*y(1)*y(1)*y(1)			- b(100)*y(40)*y(47); 	% 2 + 3 -2 -3
r101 = f(101)*y(27)*y(1)   					- b(101)*y(31)*y(2);	% 2 + 1 -2 -1
r102 = f(102)*y(27)*y(1)*y(1)	 			- b(102)*y(30)*y(2);	% 2 + 1 -2 -2
r103 = f(103)*y(28)*y(1)*y(1)	        	- b(103)*y(39)*y(47);	% 3 + 3 -4 -2
r104 = f(104)*y(28)*y(1)*y(1)		   		- b(104)*y(30)*y(2);	% 4 + 1 -3 -2 
r105 = f(105)*y(29)*y(1)*y(1)*y(1)			- b(105)*y(41)*y(46);	% 3 + 3 -3 -3
r106 = f(106)*y(29)*y(1)   					- b(106)*y(31)*y(2);	% 3 + 1 -3 -1
r107 = f(107)*y(30)*y(1)*y(1)      			- b(107)*y(41)*y(47);	% 3 + 3 -4 -2
r108 = f(108)*y(30)*y(1)   					- b(108)*y(32)*y(2);	% 4 + 1 -4 -1
r109 = f(109)*y(31)*y(1)*y(1)*y(1)			- b(109)*y(42)*y(46);	% 3 + 3 -3 -3
r110 = f(110)*y(31)*y(1)*y(1)*y(1)  	    - b(110)*y(41)*y(47);	% 3 + 3 -3 -3
r111 = f(111)*y(31)*y(1)			        - b(111)*y(32)*y(2);	% 2 + 1 -2 -1
r112 = f(112)*y(32)*y(1)*y(1)*y(1)      	- b(112)*y(42)*y(47);	% 2 + 3 - 2 -3
r113 = f(113)*y(43)*y(1)   					- b(113)*y(44)*y(2); 	% 1 + 1 -1 -1  
r114 = f(114)*y(44)*y(1)*y(1)   			- b(114)*y(45)*y(2);	% 2 + 1 - 1 -2
r115 = f(115)*y(45)*y(1)*y(1)		   		- b(115)*y(46)*y(2);	% 3 + 1 -2 -2
r116 = f(116)*y(46)*y(1)   					- b(116)*y(47)*y(2); 	% 3 + 1 -3 -1
r117 = f(117)*y(33)*y(1)   					- b(117)*y(44)*y(44);	% 1 + 1 -1 -1
r118 = f(118)*y(33)*y(1)   					- b(118)*y(34)*y(2);  	% 1 + 1 -1 -1
r119 = f(119)*y(34)*y(1)*y(1)    			- b(119)*y(44)*y(45); 	% 1 + 2 -1 -2
r120 = f(120)*y(34)*y(1)*y(1)   			- b(120)*y(35)*y(2);  	% 2 + 1 -1 -2
r121 = f(121)*y(34)*y(1)	  				- b(121)*y(37)*y(2);	% 1 + 1 -1 -1
r122 = f(122)*y(35)*y(1)*y(1)		   		- b(122)*y(44)*y(46); 	% 1 + 3 -2 -2
r123 = f(123)*y(35)*y(1)*y(1)	 			- b(123)*y(36)*y(2);	% 3 + 1 -2 -2 
r124 = f(124)*y(35)  						- b(124)*y(38)*y(2);	% 1 + 1 -2
r125 = f(125)*y(36)*y(1)      				- b(125)*y(44)*y(47);	% 1 + 3 - 3 -1
r126 = f(126)*y(36)*y(1)      				- b(126)*y(39)*y(2); 	% 3 + 1 - 3 -1
r127 = f(127)*y(37)*y(1)*y(1)*y(1) 			- b(127)*y(45)*y(45);	% 2 + 2 -1 -3
r128 = f(128)*y(37)*y(1)  					- b(128)*y(38)*y(2);	% 1 + 1 -1 -1
r129 = f(129)*y(38)*y(1)*y(1)*y(1)*y(1)		- b(129)*y(45)*y(46);	% 2 + 3 -1 -4
r130 = f(130)*y(38)*y(1)*y(1)*y(1)			- b(130)*y(39)*y(2);	% 3 + 1 -1 -3
r131 = f(131)*y(38)*y(1)			 		- b(131)*y(40)*y(2); 	% 2 + 1 -1 -1
r132 = f(132)*y(39)*y(1)*y(1)		  		- b(132)*y(45)*y(47);	% 2 + 3 - 3 -2
r133 = f(133)*y(39)		   					- b(133)*y(41)*y(2)*y(1);	% 1 + 1 + 1 -3
r134 = f(134)*y(40)*y(1)*y(1)*y(1)*y(1)		- b(134)*y(46)*y(46);	% 3 + 3 -2 -4
r135 = f(135)*y(40)		   					- b(135)*y(41)*y(2);	% 1 + 1 -2 
r136 = f(136)*y(41)*y(1)*y(1)*y(1)*y(1)*y(1)     - b(136)*y(46)*y(47);	% 3 + 3 -1 -5
r137 = f(137)*y(41)*y(1)*y(1)   			- b(137)*y(42)*y(2);	% 2 +1 -1 -2
r138 = f(138)*y(42)*y(1)*y(1)*y(1)*y(1)     - b(138)*y(47)*y(47);	% 3 + 3 -2 -4

% calculating rate of coverage for all 16 intermediates and 1 free site

F = zeros(47,1);

F(2)  =  2*r3 +r9 +r10 + r11 + r12 +r15 + r18 +r19 + r22 +r23 +r24 + r26 +r27 + r30 +r31 +r32 +r34 + r37 +r38 +r41 +r42 +r43 +r45 +r46 + r49 + r50 +r51+ r54 + r55 +r58 +r59 + r62 + r63 + r65 + r66 + r69 + r70 +r71 + r73 + r76 +r77 + r80 +r83 + r84 + r87 +r88 +r90 +r91 + r94 + r95 + r98 + r101 + r102 + r104 + r106 +r108 + r111 + r113 +r114 +r115 + r116 + r118 + r120 +r121 + r123 + r124 + r126 + r128 + r130 +r131 + r133 + r135 + r137; 
F(3)  =  r1 - r9 -r10 -r13; 
F(4)  =  r9 - r11 -r14 - r15;
F(5)  =  r10 -r12 - r16 -r17 -r18 -r19; 
F(6)  =  r2 + r11 + r12 -r28 -r29 -r30 - r31 -r32 ; 
F(7)  =  r19 - r20 -r21 -r22 -r23 -r24; 
F(8)  =  r18 - r25 -r26 -r27; 
F(9)  =  r15 - r33 -r34; 
F(10) =  r22 - r35 -r36 -r37 -r38; 
F(11) =  r24 +r26 - r39 - r40 -r41 -r42 -r43; 
F(12) =  r27 + r32 - r44 - r45 -r46;
F(13) =  r23 + r31 - r47 -r48 -r49 - r50 -r51; 
F(14) =  r30 + r34 - r52 -r53 -r54 - r55; 
F(15) =  r38 + r49 -r56 -r57 - r58 -r59; 
F(16) =  r37 + r41 - r60 -r61 -r62 - r63;
F(17) =  r43 - r64 - r65 - r66 ; 
F(18) =  r42 +r45 +r51 - r67 -r68 - r69 - r70 - r71; 
F(19) =  r46 +r54 - r72 -r73 ; 
F(20) =  r4 + r50 + r55 - r74 - r75 -r76 - r77; 
F(21) =  r58 + r76 - r78 -r79 - r80; 
F(22) =  r59 + r62 + r69 - r81 -r82 -r83 -r84 ; 
F(23) =  r66 + r71 - r89 - r90 - r91; 
F(24) =  r63 + r65 - r85 -r86 -r87 -r88; 
F(25) =  r70 + r73 + r77 -r92 - r93 - r94 - r95; 
F(26) =  r80 + r83 + r94 - r96 -r97 -r98; 
F(27) =  r84 + r87 +r90 - r99 - r100 -r101 - r102;
F(28) =  r88 - r103 - r104; 
F(29) =  r91 + r95 - r105 - r106; 
F(30) =  r102 + r104 - r107 -r108; 
F(31) =  r98 + r101 + r106 - r109 -r110 -r111; 
F(32) =  r108 + r111 - r112; 
F(33) =  r5 - r117 - r118; 
F(34) =  r13 + r16 + r20 + r36 + r118 -r119 -r120 -r121 ;
F(35) =  r14 + r29 +  r48 + r57 + r120 - r122 -r123 -r124;
F(36) =  r33 + r53 + r74 + r79 +r123 - r125- r126; 
F(37) =  r6 + r17 +r25 +r40 + r60 +r121 - r127 -r128; 
F(38) =  r21 +r28 + r39 +r44 + r64 +r68 + r82 + r86 +r124 + r128 - r129 - r130 -r131; 
F(39) =  r35 + r52 +r61 + r72 +r85 + r93 + r97 + r103 + r126 + r130 - r132 -r133; 
F(40) =  r7 + r47 +r67 + r89 + r100 + r131 - r134 - r135; 
F(41) =  r56 + r75 +r81 + r92 + r99 + r105 + r107 + r110 + r133 + r135 -r136 -r137;
F(42) =  r78 + r96 + r109 + r112 + r137 - r138; 
F(43) =  r8 - r113; 
F(44) =  r13 + r14 + r17 + r21 +r28 + r33 +r35 +r47 + r52 +r56 + r75 + r78 + r113 -r114 + 2*r117 + r119 + r122 +r125;
F(45) =  r16 +r25 + r29 + r39 +r44 + r53 + r61 + r67 + r72 + r81 + r92 + r96 +r114  -r115 + r119 + 2*r127 + r129 + r132;
F(46) =  r20 + r40 + r48 + r64 + r68 + r74 +r85 + r89 + r93  + r99 + r105 + r109 + r115 - r116 + r122 + r129 + 2*r134 + r136;
F(47) =  r36 + r57 + r60 + r79 + r82 + r86 + r97 + r100 +r103 + r107 + r110 +r112 + r116 + r125 + r132 + r136 + 2*r138;
 
F(1)  = -(F(2) + F(3) +F(4) + F(5) + F(6) + 2*F(7) + 2*F(8) + F(9) + 3*F(10)  + 3*F(11)  + 2*F(12) + F(13) + F(14) + 2*F(15) +  3*F(16) + ...
		4*F(17) +  2*F(18) + 3*F(19) + 2*F(20) + F(21) + 3*F(22) + 2*F(23) + 4*F(24) + 2*F(25) + 2*F(26) + 2*F(27) + 4*F(28) + 4*F(29) + 2*F(30) + ...
		2*F(31) + 2*F(32) + F(33) + F(34) + 2*F(35) + 3*F(36) + F(37) + F(38) + 3*F(39) + 2*F(40) + F(41) + 2*F(42) + F(43) + F(44) + 2*F(45) + 3*F(46) + 3*F(47)) ;

end