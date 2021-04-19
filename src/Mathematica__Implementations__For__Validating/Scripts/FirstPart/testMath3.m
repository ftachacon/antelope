#!/usr/local/bin/MathematicaScript -script

Print["\n******************************************"]
Print["\n******************************************"]
Print["\n       A. Chacon's solver to SBEs\n April 12, 2019"]
Print["\n*********************************\n\n"]

ClearAll["Global`*"] ;
Needs["DifferentialEquations`NDSolveProblems`"];
Needs["DifferentialEquations`NDSolveUtilities`"];


Fehlbergamat = { {1/4}, {3/32, 9/32}, {1932/2197, -7200/2197, 7296/2197}, {439/216, -8,  3680/513, -845/4104}, {-8/27, 2, -3544/2565, 1859/4104, -11/40}}; Fehlbergbvec = {25/216, 0, 1408/2565, 2197/4104, -1/5, 0}; Fehlbergcvec = {1/4, 3/8, 12/13, 1, 1/2}; Fehlbergevec = {-1/360, 0, 128/4275, 2197/75240, -1/50, -2/55}; FehlbergCoefficients[4, p_] :=  N[{Fehlbergamat, Fehlbergbvec, Fehlbergcvec, Fehlbergevec}, p];

DOPRIamat = { {1/5}, {3/40, 9/40}, {44/45, -56/15, 32/9}, {19372/6561, -25360/2187, 64448/6561, -212/729}, {9017/3168, -355/33, 46732/5247, 49/176, -5103/18656}, {35/384, 0, 500/1113, 125/192, -2187/6784, 11/84}}; DOPRIbvec = {35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0}; DOPRIcvec = {1/5, 3/10, 4/5, 8/9, 1, 1}; DOPRIevec = {71/57600, 0, -71/16695, 71/1920, -17253/339200,  22/525, -1/40}; DOPRICoefficients[5, p_] :=  N[{DOPRIamat, DOPRIbvec, DOPRIcvec, DOPRIevec}, p];


rootPathA ="/Users/achacon/Documents/Workplace/MichaelC/TopologyProj/\
CodeDeveloping/March19_2019/";
(*NotebookDirectory[];*)

(*"/Users/achacon/Documents/Workplace/MichaelC/TopologyProj/\
CodeDeveloping/March19_2019/";*)

$PreRead = (# /.
            s_String /;
            StringMatchQ[s, NumberString] &&
            Precision@ToExpression@s == MachinePrecision :>
            s <> "`26" &);

efac        = 27.2;

(************************************************
input = "input.in";
InPutAC0 = Import[rootPathA <> input,"Table" ]
q = Flatten[InPutAC0];
Print["\nInput Data, yindex = ",Flatten[InPutAC0],"\n"]
************************************************)



toFixedWidth[n_Integer, width_Integer] :=
StringJoin[PadLeft[Characters[ToString[n]], width, "0"]];
toNumberedFileName[fname_, n_Integer] :=
ToString@StringForm[fname <> toFixedWidth[n, 5]]



lfile       = "laserParam.dat";
mfile       = "MomGridParam.dat";
bandfile    = "BandsParam.dat";
laserfile   = "outlaserdata.dat";
bandoutput  = "energydiffEg.dat";

lfile       = rootPathA <> lfile;
mfile       = rootPathA <> mfile;
bandfile    = rootPathA <> bandfile;
laserfile   = rootPathA <> laserfile;
bandoutput  = rootPathA <> bandoutput;

Print["\nCurrent directory: ",rootPathA,"\n\n"];

(***********************************************)
(*PARAMETERS FOR: Time Axis and more parameters*)
(* Laser parameters *)
(***********************************************)

e00         = 0.003;            (*laser field amplitude*)
w00         = 0.014;            (*laser freq.*)
T00         = 2.*Pi/w00;        (*laser cycle or period*)
Ncycles     = 4.;                (*No. of cycles*)
alpha0      = Ncycles*T00;      (*time width or integration window*)


(***********************************************)
(******* TIME STEP *******)
dt          = 1.5;           (*T00/160.0/2.50;*)        (*Time step, T00/160.06/2.;*)
tmax        = 3.00*alpha0;  (*2.0*alpha0;*)(*3.50*alpha0;*)
Ntime       = Floor[ 2.*tmax/dt ];
Nsnaper     = 1;

Print["\n\n-------------------------------------\n\nTime-Step dt    =   ", dt, "  a.u."];
Print["No. of Time steps, Ntime     =   ", Ntime];
Print["-------------------------------------\n"]
Print["\n************************\nLASER Parameters"]
Print["E0           =   ",e00, "  a.u."]
Print["w0           =   ",w00, "  a.u."]
Print["Ncycles      =   ",Ncycles, "  a.u."]
Print["*********************************************"]



(***********************************************)
(*Haldane model parameters*)
t1          = 0.075;
t2          = t1/3.0;
M           = 2.54*t2;
phi         = 1.16;
T2          = T00;                  (*Dephasing*)

iregular    = t1/6.;        (*Regularization of Dipole*)
rregular    = t1/50.;        (*t1/30/1000.;*)
creg        = t1/10.;           (*Regularization of connection VERY IMPORTANT*)

(*********************************************)
(********* LATTICE CONSTANT a0 ***********)
a0      = 1.000/0.529;

Print["\n*******************************\n"]
Print["\nLATTICE PARAMETERS:\n"]
Print["Lattice Const.    a0           =   ",a0,  "   a.u."]
Print["Hopping param.NN  t1           =   ",t1, "   a.u."]
Print["Hopping param.NNN t2           =   ",t2, "   a.u."]
Print["Phase or Mag. Flux phi0        =   ",phi,"   a.u."]
Print["Dephasing T2                   =   ",T2,  "   a.u."]



(*********************************************)
(* Momentum grid Parameter *)
(*********************************************)
(**************  NUMBERS OF POINTS AT THE BZ *****************)
dir     = 1; (*direction*)

Nkx     = 501;  (* 301;251,151,101,This controls the momentum grid steps *)
Nky     = 101;

aindex0 = 50;
aindexf = aindex0;(*0*Nky;*)

(***********************************************)
(*Creating momentum axis limits*)

xtemp       = N[ 2.*Pi/Sqrt[3.]/a0 ];
ytemp       = N[ 2.*Pi/3./a0 ];

kxMin       = - xtemp; kxMax = xtemp;
kyMin       = - ytemp; kyMax = ytemp;
kxMax0      = + kxMax; (*Max grid space momentum *)

dky         = 2.*kyMax/(Nky - 1.);
dktest      = 2.*kxMax/(Nkx - 1.);

(***********************************************)
Print["\n\n************************************\nTotal No. of Momen. Steps, Nkx = ", Nkx];
Print["Momen. grid-step dkx           =  ", dktest, "   a.u."];
Print["Nky                            =  ", Nky ];
Print["dky                            =  ", dky,   "    a.u."];
Print["kxMax                          =  ",kxMax,  "    a.u."];
Print["kyMax                          =  ",kyMax,  "    a.u."];
Print["\n/*******************************************************/\n"];
(***********************************************)


(*Next Nearest Neighbor, NNN, vectors of the HoneyComb Lattice*)

b1          = { +Sqrt[3.], 0. }*a0;
b2          = { -Sqrt[3.]/2., +3/2.}*a0;
b3          = { -Sqrt[3.]/2., -3/2.}*a0;
b           = { b1, b2, b3 };

(*Nearest Neighbor, NN, vectors of the HoneyComb Lattice*)
a1          = { 0., 1. }*a0;
a2          = {-Sqrt[3.]/2., -1/2. }*a0;
a3          = { Sqrt[3.]/2., -1/2. }*a0;
a           = { a1, a2, a3 };

ax          = { a1[[1]], a2[[1]], a3[[1]] };
ay          = { a1[[2]], a2[[2]], a3[[2]] };
bx          = { b1[[1]], b2[[1]], b3[[1]] };
by          = { b1[[2]], b2[[2]], b3[[2]] };

(*******
   K' and K points with respect to the gamma points on the BZ
 ***)
 
K1          = {-4.*Pi/(3.0*Sqrt[3.0] a0), 0};
K2          = { 4.*Pi/(3.0*Sqrt[3.0] a0), 0};


MomParameters   = {dktest, Nkx, kxMax, dky, Nky, kyMax,K1[[1]],K1[[2]]};
Export[ mfile, MomParameters, "Data"];



(*Laser field definition*)
Ef := Function[t,
               e00*Cos[w00*t]*Exp[ -2.0*(t/alpha0)^2 ]
               ];

Af := Function[t,
                -e00*Sin[w00*t]*Exp[ -2.0*(t/alpha0)^2 ]/w00
                ];
 

timeaxes := Function[
                     {dt, Nt},
                     
                     tmax = dt*(Nt - 1)/2.;
                     wmax = Pi/dt;
                     
                     dw = 2 Pi/Nt/dt;
                     
                     {Table[-tmax + i dt, {i, 0, Nt - 1}]
                         , Table[-wmax + i dw, {i, 0, Nt - 1}]
                     }
                     
                     ];


 
 (* Momentum axis *)
 momaxis := Function[
                     {NkVar, kminVar},
                     
                     dkVar = Abs[2*kminVar/(NkVar - 1)];
                     Table[kminVar + dkVar i, {i, 0, NkVar - 1}]
                     
                     ];
 
 








(******************************************************)
(*****************************************************)
(*Haldane Model (HM) and topological Chern insulators

Hamiltonian and its "components" B0, B1,...
H(k) = B0*sig0 + B1*sig1 + B2*sig2 + B3*sig3
where, sig1, sig2, sig3 are Pauli's matrices, sig0 identity matrix

*)




(*HALDANE MODEL, DIFINITION OF HAMILTONIAN COHEFICIENTS*)

B0 := Function[{kx, ky}
               , 2.*t2*Cos[phi]*
               Sum[
                   Cos[ kx *bx[[i]] + ky *by[[i]]  ]
                   , {i, 1, 3}
                   ]
               ];


B1 := Function[{kx, ky}
               , t1*Sum [
                         Cos[ kx *ax[[i]] + ky *ay[[i]] ]
                         , {i, 1, 3}
                         ]
               ];


B2 := Function[{kx, ky}
               , t1*Sum[
                        Sin[ kx *ax[[i]] + ky *ay[[i]] ]
                        , {i, 1, 3}
                        ]
               ];


B3 := Function[ {kx, ky}
               , M - 2.*t2*Sin[phi] Sum[
                                          Sin[ kx *bx[[i]] + ky *by[[i]] ]
                                          , {i, 1, 3}
                                          ]
               ];


BNorm := Function[{kx, ky}
                  , Sqrt[
                         B1[kx, ky]^2 + B2[kx, ky]^2 +
                         B3[kx, ky]^2
                         ]
                  ];



(*Definition of quantize indexes*)

nB1 := Function[{kx, ky}
                , B1[kx, ky]/
                BNorm[kx, ky]
                ];


nB2 := Function[{kx, ky}
                , B2[kx, ky]/
                BNorm[kx, ky]
                ];


nB3 := Function[{kx, ky}
                , B3[kx, ky]/
                BNorm[kx, ky]
                ];


Phi := Function[{kx, ky}
                , ArcTan[
                         B1[kx, ky],
                         B2[kx, ky]
                         ]
                ];

Theta := Function[{kx, ky}
                  , ArcCos[nB3[kx, ky]]
                  ];





(*GRADIENTS OF B-SETs HAM. Components*)
B0Grad := Function[{kx, ky}
                   , Grad[ B0[kkx, kky]
                          , {kkx, kky}
                          ] /. kkx -> kx /. kky -> ky
                   ];

 
B1Grad := Function[{kx, ky}
                   , Grad[ B1[kkx, kky]
                          , {kkx, kky}] /. kkx -> kx /. kky -> ky
                   ];


B2Grad := Function[{kx, ky}
                   , Grad[B2[kkx, kky]
                          , {kkx, kky}] /. kkx -> kx /. kky -> ky
                   ];


B3Grad := Function[{kx, ky}
                   , Grad[B3[kkx, kky]
                          , {kkx, kky}] /. kkx -> kx /. kky -> ky
                   ];


BNormGrad := Function[{kx, ky}
                      , Grad[BNorm[kkx, kky]
                             , {kkx, kky}] /. kkx -> kx /. kky -> ky
                      ];




(*GRADIENTS OF Bloch Sphere angles*)

PhiGrad := Function[{kx, ky,qeps},
                    
                    {
                        If[ky < 0.,
                           
                           -(B2Grad[kx, ky][[1]] * B1[kx, ky] -
                             B1Grad[kx, ky][[1]] * B2[kx, ky]) /(B1[kx, ky]^2 + B2[kx, ky]^2 + qeps)
                           ,
                           (B2Grad[kx, ky][[1]] * B1[kx, ky] -
                            B1Grad[kx, ky][[1]] * B2[kx, ky]) /(B1[kx, ky]^2 + B2[kx, ky]^2 + qeps)
                           ]
                        ,
                        
                        
                        (B2Grad[kx, ky][[2]]*B1[kx, ky] - B1Grad[kx, ky][[2]]*B2[kx, ky]) /(B1[kx, ky]^2 + B2[kx, ky]^2 + qeps)
                        
                    }
                    ]  ;



ThetaGrad :=
Function[{kx, ky, qeps},
         
         {
             -(B3Grad[kx, ky][[1]] * BNorm[kx, ky] - BNormGrad[kx, ky][[1]] *B3[kx, ky])
               / (Sqrt[BNorm[kx, ky]^2 - B3[kx, ky]^2 ] + qeps)/(BNorm [kx,ky] + qeps)
             
             
             ,
             
             If[ ky < 0.,
                
                (B3Grad[kx, ky][[2]] * BNorm[kx, ky] - BNormGrad[kx, ky][[2]] * B3[kx, ky])/(Sqrt[ BNorm[kx, ky]^2 - B3[kx, ky]^2 ] + qeps)/(BNorm [kx,ky] + qeps)
                ,
                
                -(B3Grad[kx, ky][[2]] *BNorm[kx, ky] - BNormGrad[kx, ky][[2]] *
                  B3[kx, ky]) /(Sqrt[ BNorm[kx, ky]^2 - B3[kx, ky]^2 ] + qeps)/(BNorm [kx,ky] + qeps)
                
                ]
             
         }
         
         ]  ;



(*Definition of Energy dispersion for conduction and valance band*)

Ec := Function[{kx, ky}
               , B0[kx, ky]  + BNorm[kx, ky]
               ];



Ev := Function[{kx, ky}
               ,B0[kx, ky]  - BNorm[kx, ky]
               ];



(*Group velocities of the conduction and valance band*)

cGroupVel := Function[{kx, ky}
                      , Grad[
                             Ec[kkx, kky], {kkx, kky}] /.
                      kkx -> kx /. kky -> ky
                      ];



vGroupVel := Function[{kx, ky},
                      Grad[Ev[kkx, kky], {kkx,
    kky}] /. kkx -> kx /. kky -> ky
                      ];



(*Berry connection and curvature of the conduction band*)

BerryConnectionReg :=
Function[{kx, ky, qeps}
         ,
         
         1/2.*B3[kx, ky]*PhiGrad[kx, ky, qeps ]/(BNorm[kx, ky] + qeps)
         
         ];




(*Berry curvature of the conduction band*)

cBerryCurvature := Function[{kx, ky}
                            , 1/2.*Cross[
                                        
                                        Grad[nB3[kkx, kky]
                                             , {kkx, kky, kkz} ] /. kkx -> kx /. kky -> ky
                                        
                                        , Grad[Phi[kkx, kky]
                                               , {kkx, kky, kkz}] /. kkx -> kx /. kky -> ky
                                        ]
                            ];


 
 

(*Dipole transition matrix element d_cv = i u^*_c d/dk u_v*)
RDipoleCV := Function[{kx, ky, qeps}
                      ,
                      
                      1/2.*Sin[Theta[kx, ky] ]
                      PhiGrad[kx, ky, qeps]
                      
                      ]; (*Real part*)



IDipoleCV := Function[{kx, ky, qeps}
                      , 1/2 ThetaGrad[kx, ky, qeps]
                      ]; (*Imaginary part*)



Dipole := Function[{kx, ky}
                   ,
                   
                   RDipoleCV[kx, ky, rregular]
                   + I*IDipoleCV[kx, ky, iregular]
                   
                   
                   (*{1.0 + I, 1.0 + I}
                   *)
                   ];





(*Energy Gap, Eg = Ec- Ev and Berry Connection 'Chig' Gap*)

Eg := Function[
               {kx, ky}
               , Ec[kx, ky]
               - Ev[kx, ky]
               ];



(*Chig = Chic - Chiv *)
Chig := Function [
                  {kx, ky, creg},
                  2.*BerryConnectionReg[kx, ky, creg]
                  ];


(*****END OF HALDANE MODEL STRUCTURE*****)
(*******************************************************)



(* Energy Band Gap*)
 
 (***********************)
 (*Max and Min gap*)
 
 
 dkx1 = 0.061;
 dky1 = 0.071;
 
 Tac = Table[
             Ec[kx, ky], {kx, kxMin, kxMax,
                 dkx1}, {ky, kyMin, kyMax, dky1}];
 
 
 Tav = Table[
             Ev[kx, ky], {kx, kxMin, kxMax,
                 dkx1}, {ky, kyMin, kyMax, dky1}];
 
 
 (*******************************************)
 (********************************)
 (*Calculation of the CHERN No.*)
 
 chernNo    = -NIntegrate[cBerryCurvature[kx, ky][[3]], {kx, -kxMax, kxMax}, {ky, -kyMax, kyMax}  ]/(4.*Pi);
 
 
Print["\n**********************************"];
Print["M/t2         =   ", M/t2 ];
Print["t2/t1        =   ", t2/t1];
Print["t1           =   ", t1*efac, "  eV"];
Print["t2           =   ", t2*efac, "  eV"];
Print["M            =   ", M*efac, "  eV"];
Print["phi0         =   ", phi, "  rad."];
 
egmin0 = Min[Tac - Tav];
egmax0 = Max[Tac - Tav];


Print["\n***********************************"];
Print["Energy band gaps and Chern No."];
Print["Eg               =       ", N[egmin0,16]*efac, " eV"];
Print["Max gap          =       ", N[egmax0*efac,16], " eV"];
Print["Chern No.        =       ", chernNo,"\n"];


haldaneparam    = {a0, t1, t2, M, phi, T2, egmax0, egmin0, chernNo};
Export[ bandfile, haldaneparam, "Data"];
 

 
 
 (************************************************)
 (************************************************)
 (*Numerical time axis for time evolution of SBEs*)
 
 tme = timeaxes[dt, Ntime];
 time = tme[[]][[1]];
 freq = tme[[]][[2]];
 

 (*Momentum Axis*)
 kxx = momaxis[Nkx, kxMin]; (*momentum space axis*)
 dkx = kxx[[2]] - kxx[[1]];(*space step => momentum step*)
 
 
 LaserOut = Table[{j - 1, time[[j]], Ef[ time[[j]] ], Af[ time[[j]] ]  }, {j, 1,Ntime}];
 Export[laserfile, LaserOut, "Data"];
 
 BandsOutPut = Table[ { kxx[[i]], Eg[kxx[[i]],0] },{i,1,Nkx} ];
 Export[bandoutput,BandsOutPut,"Data" ];
 
 (**************************)
 (*Rabbi Freq. *)
 ROmega := Function[{kx, ky, t}
                    
                    , Ef[ t ]*Dipole[ kx, ky ][[dir]]
                    
                    ];
 
Print["\n********************************************************"]
Print["CFL condition E0*dt/dk/2.   =    ", e00*dt/dkx/2., "    "];
Print["\n***************************************\n\n"]
 
Print["\n-------------------------------------\nTime-Step dt    =   ", dt, "  a.u."];
Print["No. of Time steps, Ntime     =   ", Ntime];
Print["Nx    =    ",Nkx,";   dkx   =  ",dkx, " a.u."];
Print["Ny    =    ",Nky,";   dky   =  ",dky, " a.u."];
 

(*******************************************)
(****
        Evol. Equations of SBEs       ****)
(*******************************************)
(*********
 Dinamical evolution of populations and SBEs defined by Eqs 1 -- 3,
    according to the constriction
 *****************)
(*[SlideView[*)

solver = Function[ky,
                                  
            NDSolve[
                                          
                    {
                                              
                        (*******************************************)
                            (****
                                        Evol. Equations of SBEs       ****)
                            (*******************************************)
                    
                                              
                        (*Equation1*)
                                              
                        D[pi[kx, t], t]  == -I*( Eg[kx+Af[t], ky]
                                                + Ef[t]*Chig[kx+Af[t], ky, creg][[dir]]
                                                -  I/T2
                                                ) * pi[kx, t]
                                            - I*ROmega[kx+Af[t], ky, t]*( 1. - 2.*nc[kx,t] )
                        
                        
                        (*Equation2*)
                        ,
                                              
                        D[nc[kx, t], t] (*- Ef[t] * D[nc[kx, t], kx]*) == -2.*Im[ Conjugate[ROmega[kx+Af[t], ky, t]] * pi[kx, t] ]
                                
                        
                        (*******************************************)
                        (****         \
                                    Initials1             ****)
                        (*******************************************)
                                
                        , pi[kx, -tmax] == 0.
                        
                                              
                        (*Initials2*)
                        
                        , nc[kx, -tmax] == 0.
                            
                        (*******************************************)
                                              
                        }
                                          
                        , {pi, nc}
                                          
                        , {kx, kxMin, kxMax}
                                          
                        , {t, -tmax, tmax}
                    
                        , MaxStepSize-> {dkx,dt/2.}
                    
                        , AccuracyGoal -> 10
                    
                        , PrecisionGoal -> 10
                    
                        (*, Method -> {"FixedStep",  Method -> {"ImplicitRungeKutta", "DifferenceOrder" -> 4,  "ImplicitSolver" -> {"Newton",  AccuracyGoal -> MachinePrecision,  PrecisionGoal -> MachinePrecision,  "IterationSafetyFactor" -> 1}} }
                         *)
                    
                        ,  Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 5,  "Coefficients" -> DOPRICoefficients, "StiffnessTest" -> False}
                    
                    
                    (*, StartingStepSize -> 1/10*)
                    
                        (*,Method\[Rule]"ImplicitRungeKutta",*)
                        (*,Method\[Rule]
                        RK5*)
                ]
                
];

                

(***********************************************)
(*** Starting Do Loop on ky direction ***)
(***********************************************)
                
Do[
   
   yshift   = kyMin +  dky*qindex; (**** Momentum on ky-direction ****)
   yindex   = qindex;

   Print["\n\n\n\n/*******************************************/\ny-index = ",yindex,";      Ky = ", yshift, "  a.u. \n/************/\n"];
   Print["\nCurrent directory: ",rootPathA,"\n\n"];
   
   
   time01 = AbsoluteTime[];
   Print["\n****************Evaluating************************\ntime1 = ",time01," s"];
   
   
    mdfun = solver[yshift];
                
    pi0 = Table[
                kxtemp=kxx[[i]];
                ttemp=time[[j*Nsnaper]];
                            Evaluate[ pi[ kxtemp, ttemp ] /. mdfun]
                            , {j, 1, Ntime/Nsnaper}, {i, 1, Nkx}
                            ];
   

    nc0 = Table[
                
                xktemp = kxx[[i]];
                ttemp  = time[[j*Nsnaper]];
                            Re[Evaluate[ nc[ xktemp, ttemp ] /. mdfun] ]
                            , {j, 1, Ntime/Nsnaper}, {i, 1, Nkx}
                ];
                
                
    pi1 = Table[
                            pi0[[j]][[i]][[1]],
                            {j, 1, Ntime/Nsnaper}, {i, 1, Nkx}
                            ];
   
   
                
    nc1 = Table[
                            nc0[[j]][[i]][[1]],
                            {j, 1, Ntime/Nsnaper}, {i, 1, Nkx}
                            ];
   
   time02=AbsoluteTime[];
   Print["\nEvaluation time = ",time02-time01," s"];
   Print["\n************\nStarting the output\n"];
                
                
    Repi2       = Re[pi1];
    ImPi2       = Im[pi1];
                
               (* j=9;
                i=9;
                Print["\npi = ", pi0[[j]][[i]] ];
                *)
                
                (* OUTPUTS *)
                
    fname0      = toNumberedFileName["CoherenceMathRealPiNx_"    , Nkx];
    fname0      = toNumberedFileName[fname0 <> "_kyIndex_"       , qindex];
    DataPath0   = rootPathA <> fname0  <>  ".bin";
                
                
    fname01     = toNumberedFileName["CoherenceMathImagPiNx_"   , Nkx];
    fname01     = toNumberedFileName[fname01 <> "_kyIndex_"     , qindex];
    DataPath01  = rootPathA <> fname01  <>  ".dat";
                
    fname1      = toNumberedFileName["OccupationCBMath_ncNx_"    , Nkx];
    fname1      = toNumberedFileName[fname1 <> "_kyIndex_"       , qindex];
    DataPath1   = rootPathA <> fname1  <>  ".bin";
                
                
    BinaryWrite[DataPath0,   Repi2, "Real64"]
                (*BinaryWrite[DataPath01,  Impi2, "Real64"]*)
                
    Export[DataPath01, ImPi2, "Data"];
                
   BinaryWrite[DataPath1,   nc1,   "Real64"]
                
    (*BinaryWrite[..., "Real64"]*)
    (*OutPuts Coherence/Populations*)
                
    Clear[RePi2];
    Clear[ImPi2];
    Clear[nc1];
    Clear[pi0];
    Clear[pi1];
    Clear[nc0];
    Clear[mdfun];
   
    time01 = AbsoluteTime[];
    Print["\nOutput time = ",time01-time02," s"];
   
                
    (**************************************)
    (**DO--Conditions**)
    ,
                
    {qindex,aindex0,aindexf}
    (**************************************)
                
];

Print["\n\n*****************************************\n"]
Print["END OF THE PROGRAM AT, \nky = ",yshift, " a.u.; and\nyindex = ", yindex ]
Print["\n***********************************\n\n\n"]

