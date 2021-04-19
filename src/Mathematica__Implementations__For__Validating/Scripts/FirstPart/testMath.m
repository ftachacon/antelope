#!/usr/local/bin/MathematicaScript -script

Print["\n******************************************"]
Print["\n******************************************"]
Print["\n       A. Chacon's solver to SBEs\n April 12, 2019"]
Print["\n*********************************\n\n"]

ClearAll["Global`*"] ;
Needs["DifferentialEquations`NDSolveProblems`"];
Needs["DifferentialEquations`NDSolveUtilities`"];


rootPathA =
"/Users/achacon/Documents/Workplace/MichaelC/TopologyProj/\
CodeDeveloping/March19_2019/";

$PreRead = (# /.
            s_String /;
            StringMatchQ[s, NumberString] &&
            Precision@ToExpression@s == MachinePrecision :>
            s <> "`26" &);

efac        = 27.2;

(*input = "input.in";
InPutAC0 = Import[rootPathA <> input,"Table" ]
q = Flatten[InPutAC0];
Print["\nInput Data, yindex = ",Flatten[InPutAC0],"\n"]
*)



toFixedWidth[n_Integer, width_Integer] :=
StringJoin[PadLeft[Characters[ToString[n]], width, "0"]];
toNumberedFileName[fname_, n_Integer] :=
ToString@StringForm[fname <> toFixedWidth[n, 5]]



lfile       = "laserParam.dat";
mfile       = "MomGridParam.dat";
bandfile    = "BandsParam.dat";
laserfile   = "outlaserdata.dat";

lfile       = rootPathA <> lfile;
mfile       = rootPathA <> mfile;
bandfile    = rootPathA <> bandfile;
laserfile   = rootPathA <> laserfile;




(***********************************************)
(*PARAMETERS FOR: Time Axis and more parameters*)
(* Laser parameters *)
(***********************************************)

e00         = 0.007;            (*laser field amplitude*)
w00         = 0.014;            (*laser freq.*)
T00         = 2.*Pi/w00;        (*laser cycle or period*)
Ncycles     = 4.;                (*No. of cycles*)
alpha0      = Ncycles*T00;      (*time width or integration window*)


(***********************************************)
(******* TIME STEP *******)
dt          = T00/160.06/2.;        (*Time step*)
tmax        = 2.0*alpha0;
Ntime       = Floor[ 2.*tmax/dt ];
Nsnaper     = 1;



Print["\n\n-------------------------------------\n\nTime-Step dt    =   ", dt, "  a.u."];
Print["No. of Time steps, Ntime     =   ", Ntime];
Print["\n-------------------------------------\n"]
Print["\n\nE0       =   ",e00, " a.u."]
Print["\nw0         =   ",w00, " a.u."]
Print["\nNcycles    =   ",Ncycles, " a.u.\n"]





(***********************************************)
(*Haldane model parameters*)
t01     = 0.075;
t02     = t01/3.0;
M0      = 2.54*t02;
phi0    = 1.16;
T2      = T00;          (*Dephasing*)



regular = 10^-20;       (*Regularization of Dipole*)
creg    = t01/10.;      (*Regularization of connection VERY IMPORTANT*)




(*********************************************)
(********* LATTICE CONSTANT a0 ***********)
a0      = 1.000/0.529;




Print["\n*******************************\n"]
Print["\nLATTICE PARAMETERS:\n"]
Print["Lattice Const.    a0           =   ",a0,  "   a.u."]
Print["\nHopping param.NN  t1         =   ",t01, "   a.u."]
Print["\nHopping param.NNN t2         =   ",t02, "   a.u."]
Print["\nPhase or Mag. Flux phi0      =   ",phi0,"   a.u."]
Print["\nDephasing T2                 =   ",T2,  "   a.u."]






(*********************************************)
(* Momentum grid Parameter *)
(*********************************************)
(**************  NUMBERS OF POINTS AT THE BZ *****************)


dir     = 1; (*direction*)


Nkx     = 251;  (*This controls the momentum grid steps*)
Nky     = 100;




(***********************************************)
(*Creating momentum axis limits*)
xtemp       = N[2.*Pi/Sqrt[3.]/a0];
ytemp       = N[2.*Pi/3./a0];


kxMin       = - xtemp; kxMax = xtemp;
kyMin       = - ytemp; kyMax = ytemp;
kxMax0      = kxMax; (*Max grid space momentum *)


dky         = 2.*kyMax/(Nky - 1.);
dktest      = 2.*kxMax/(Nkx - 1);


(***********************************************)

Print["\n\nTotal No. of Momen. Steps, Nkx = ", Nkx];
Print["\nMomen. grid-step dkx           =  ", dktest, "   a.u."];
Print["\nNky                            =  ", Nky]
Print["\ndky                            =  ", dky,   "    a.u."]
Print["\nkxMax                          =  ",kxMax,  "    a.u."]
Print["\nkyMax                          =  ",kyMax,  "    a.u."]
(***********************************************)






(*Next Nearest Neighbor, NNN, vectors of the HoneyComb Lattice*)

b1      = { +Sqrt[3], 0.}*a0;
b2      = {-Sqrt[3]/2, +3/2}*a0;
b3      = {-Sqrt[3]/2, -3/2}*a0;

b       = {b1, b2, b3};

(*Nearest Neighbor, NN, vectors of the HoneyComb Lattice*)
a1      = {0, 1}*a0;
a2      = {-Sqrt[3]/2, -1/2}*a0;
a3      = {Sqrt[3]/2, -1/2}*a0;

a       = {a1, a2, a3};
aax     = { a1[[1]], a2[[1]], a3[[1]] };
aay     = { a1[[2]], a2[[2]], a3[[2]] };

bbx     = { b1[[1]], b2[[1]], b3[[1]] };
bby     = { b1[[2] ], b2[[2]], b3[[2]] };



(*K' and K points with respect to the gamma points on the BZ*)
 
K1      = {-4 Pi/(3*Sqrt[3] a0), 0};
K2      = {4 Pi/(3*Sqrt[3] a0), 0};
 
 




MomParameters   = {dktest, Nkx, kxMax, dky, Nky, kyMax,K1[[1]],K1[[2]]};
Export[ mfile, MomParameters, "Data"];





(*Laser field definition*)
Ef := Function[{t, w0, e0, alpha}, e0 Cos[w0 t] Exp[-2 (t/alpha)^2]];



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
                     
                     dkVar = Abs[2 kminVar/(NkVar - 1)];
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

B0 := Function[{kx, ky, t2, \[Phi], bx, by}
               , 2 t2 Cos[\[Phi]]*
               Sum[
                   Cos[ kx *bx[[i]] + ky *by[[i]]  ]
                   , {i, 1, 3}
                   ]
               ];


B1 := Function[{kx, ky, t1, ax, ay}
               , t1 Sum [
                         Cos[ kx *ax[[i]] + ky *ay[[i]] ]
                         , {i, 1, 3}
                         ]
               ];


B2 := Function[{kx, ky, t1, ax, ay}
               , t1 Sum[
                        Sin[ kx *ax[[i]] + ky *ay[[i]] ]
                        , {i, 1, 3}
                        ]
               ];


B3 := Function[ {kx, ky, t2, M, \[Phi], bx, by}
               , M - 2 t2*Sin[\[Phi]] Sum[
                                          Sin[ kx *bx[[i]] + ky *by[[i]] ]
                                          , {i, 1, 3}
                                          ]
               ];


BNorm := Function[{kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by}
                  , Sqrt[
                         B1[kx, ky, t1, ax, ay]^2 + B2[kx, ky, t1, ax, ay]^2 +
                         B3[kx, ky, t2, M, \[Phi], bx, by]^2
                         ]
                  ];



(*Definition of quantize indexes*)

nB1 := Function[{kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by}
                , B1[kx, ky, t1, ax, ay]/
                BNorm[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by]
                ];


nB2 := Function[{kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by}
                , B2[kx, ky, t1, ax, ay]/
                BNorm[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by]
                ];


nB3 := Function[{kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by}
                , B3[kx, ky, t2, M, \[Phi], bx, by]/
                BNorm[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by]
                ];


Phi := Function[{kx, ky, t1, t2, M, \[Phi], ax, ay}
                , ArcTan[
                         B1[kx, ky, t1, ax, ay],
                         B2[kx, ky, t1, ax, ay]
                         ]
                ];

Theta := Function[{kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by}
                  , ArcCos[nB3[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by]]
                  ];





(*GRADIENTS OF B-SETs HAM. Components*)
B0Grad := Function[{kx, ky, t2, \[Phi], bx, by}
                   , Grad[ B0[kkx, kky, t2, \[Phi], bx, by]
                          , {kkx, kky}
                          ] /. kkx -> kx /. kky -> ky
                   ];



B1Grad := Function[{kx, ky, t1, ax, ay}
                   , Grad[ B1[kkx, kky, t1, ax, ay]
                          , {kkx, kky}] /. kkx -> kx /. kky -> ky
                   ];


B2Grad := Function[{kx, ky, t1, ax, ay}
                   , Grad[B2[kkx, kky, t1, ax, ay]
                          , {kkx, kky}] /. kkx -> kx /. kky -> ky
                   ];


B3Grad := Function[{kx, ky, t2, M, \[Phi], bx, by}
                   , Grad[B3[kkx, kky, t2, M, \[Phi], bx, by]
                          , {kkx, kky}] /. kkx -> kx /. kky -> ky
                   ];


BNormGrad := Function[{kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by}
                      , Grad[BNorm[kkx, kky, t1, t2, M, \[Phi], ax, ay, bx, by]
                             , {kkx, kky}] /. kkx -> kx /. kky -> ky
                      ];




(*GRADIENTS OF Bloch Sphere angles*)

PhiGrad := Function[{kx, ky, t1, ax, ay, qeps},
                    
                    {
                        If[ky <= 0,
                           
                           -(B2Grad[kx, ky, t1, ax, ay][[1]] B1[kx, ky, t1, ax, ay] -
                             B1Grad[kx, ky, t1, ax, ay][[1]] B2[kx, ky, t1, ax, ay]) /(B1[
                                                                                          kx, ky, t1, ax, ay]^2 + B2[kx, ky, t1, ax, ay]^2 + qeps)
                           ,
                           (B2Grad[kx, ky, t1, ax, ay][[1]] B1[kx, ky, t1, ax, ay] -
                            B1Grad[kx, ky, t1, ax, ay][[1]] B2[kx, ky, t1, ax, ay]) /(B1[
                                                                                         kx, ky, t1, ax, ay]^2 + B2[kx, ky, t1, ax, ay]^2 + qeps)
                           ]
                        ,
                        
                        
                        (B2Grad[kx, ky, t1, ax, ay][[2]] B1[kx, ky, t1, ax, ay] -
                         B1Grad[kx, ky, t1, ax, ay][[2]]  B2[kx, ky, t1, ax, ay]) /(B1[
                                                                                       kx, ky, t1, ax, ay]^2 + B2[kx, ky, t1, ax, ay]^2 + qeps)
                        
                    }
                    ]  ;



ThetaGrad :=
Function[{kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by, qeps},
         
         {
             -(B3Grad[kx, ky, t2, M, \[Phi], bx, by][[1]] *
               BNorm[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by] -
               BNormGrad[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by][[1]] *
               B3[kx, ky, t2, M, \[Phi], bx, by]) /(Sqrt[
                                                         BNorm[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by]^2 -
                                                         B3[kx, ky, t2, M, \[Phi], bx, by]^2 ] + qeps)/(BNorm [kx,
                                                                                                               ky, t1, t2, M, \[Phi], ax, ay, bx, by] + qeps)
             
             
             ,
             
             If[(B3Grad[kx, ky, t2, M, \[Phi], bx, by][[2]] *
                 BNorm[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by] -
                 BNormGrad[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by][[2]] *
                 B3[kx, ky, t2, M, \[Phi], bx, by]) /(Sqrt[
                                                           BNorm[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by]^2 -
                                                           B3[kx, ky, t2, M, \[Phi], bx, by]^2 ] + qeps)/(BNorm [kx,
                                                                                                                 ky, t1, t2, M, \[Phi], ax, ay, bx, by] + qeps)
                ,
                -(B3Grad[kx, ky, t2, M, \[Phi], bx, by][[2]] *
                  BNorm[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by] -
                  BNormGrad[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by][[2]] *
                  B3[kx, ky, t2, M, \[Phi], bx, by]) /(Sqrt[
                                                            BNorm[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by]^2 -
                                                            B3[kx, ky, t2, M, \[Phi], bx, by]^2 ] + qeps)/(BNorm [kx,
                                                                                                                  ky, t1, t2, M, \[Phi], ax, ay, bx, by] + qeps)
                ]
             
         }
         
         ]  ;



(*Definition of Energy dispersion for conduction and valance band*)

Ec := Function[{kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by}
               , B0[kx, ky, t2, \[Phi], bx, by]  +
               BNorm[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by]
               ];



Ev := Function[{kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by}
               , B0[kx, ky, t2, \[Phi], bx, by]  -
               BNorm[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by]
               ];



(*Group velocities of the conduction and valance band*)

cGroupVel := Function[{kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by}
                      , Grad[
                             Ec[kkx, kky, t1, t2, M, \[Phi], ax, ay, bx, by], {kkx, kky}] /.
                      kkx -> kx /. kky -> ky
                      ];



vGroupVel := Function[{kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by},
                      Grad[Ev[kkx, kky, t1, t2, M, \[Phi], ax, ay, bx, by], {kkx,
    kky}] /. kkx -> kx /. kky -> ky
                      ];



(*Berry connection and curvature of the conduction band*)

BerryConnectionReg :=
Function[{kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by, qeps}
         ,
         
         1/2 B3[kx, ky, t2, M, \[Phi], bx, by]
         /(BNorm[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by] + qeps)
         PhiGrad[kx, ky, t1, ax, ay, qeps ]
         
         ];




(*Berry curvature of the conduction band*)

cBerryCurvature := Function[{kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by}
                            , 1/2 Cross[
                                        
                                        Grad[nB3[kkx, kky, t1, t2, M, \[Phi], ax, ay, bx, by]
                                             , {kkx, kky, kkz} ] /. kkx -> kx /. kky -> ky
                                        
                                        , Grad[Phi[kkx, kky, t1, t2, M, \[Phi], ax, ay]
                                               , {kkx, kky, kkz}] /. kkx -> kx /. kky -> ky
                                        ]
                            ];


cBerryCurvatureNew :=
Function[{kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by, qeps}
         ,
         -1/2 Sin[Theta[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by]]
         *Cross[
                
                { ThetaGrad[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by, qeps][[
                                                                              1]]
                    , ThetaGrad[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by, qeps][[2]]
                    , 0}
                
                ,
                
                {PhiGrad[kx, ky, t1, ax, ay, qeps ][[1]]
                    , PhiGrad[kx, ky, t1, ax, ay, qeps ][[2]]
                    , 0}
                
                ]
         ];


(*Dipole transition matrix element d_cv = i u^*_c d/dk u_v*)
RDipoleCV := Function[{kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by, qeps}
                      ,
                      
                      1/2 Sin[Theta[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by] ]
                      PhiGrad[kx, ky, t1, ax, ay, qeps]
                      
                      ]; (*Real part*)



IDipoleCV := Function[{kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by, qeps}
                      , 1/2 ThetaGrad[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by, qeps]
                      ]; (*Imaginary part*)



Dipole := Function[{kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by, qeps}
                   ,
                   
                   RDipoleCV[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by, qeps]
                   + I IDipoleCV[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by, qeps]
                   
                   
                   ];


CrossDipole :=
Function[{kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by, qeps}
         , Abs[
               Conjugate[
                         Dipole[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by, qeps][[2]]
                         ]*
               Dipole[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by, qeps][[1]]
               ]
         ];



(*Energy Gap, Eg = Ec- Ev and Berry Connection 'Chig' Gap*)

Eg := Function[
               {kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by}
               , Ec[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by]
               - Ev[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by]
               ];



(*Chig = Chic - Chiv *)
Chig := Function [
                  {kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by, creg},
                  2 BerryConnectionReg[kx, ky, t1, t2, M, \[Phi], ax, ay, bx, by,
                                       creg]
                  ];


(*****END OF HALDANE MODEL STRUCTURE*****)
(*******************************************************)








(* Energy Band Gap*)
 
 (***********************)
 (*Max and Min gap*)
 
 
 dkx1 = 0.061;
 dky1 = 0.071;
 
 Tac = Table[
             Ec[kx, ky, t01, t02, M0, phi0, aax, aay, bbx, bby], {kx, kxMin, kxMax,
                 dkx1}, {ky, kyMin, kyMax, dky1}];
 
 
 Tav = Table[
             Ev[kx, ky, t01, t02, M0, phi0, aax, aay, bbx, bby], {kx, kxMin, kxMax,
                 dkx1}, {ky, kyMin, kyMax, dky1}];
 
 
 
 
 Print["\n\n\n**********************************\nM/t2      =   ", M0/t02, ";\nt2/t1       =   ", t02/t01, ";\nt1           =   ",
       t01*efac, " eV;\nt2      =   ", t02*efac, "    eV;\nM     =   ",
       M0*efac, " eV; \nphi0        =   ", phi0, "  rad."]
 
 egmin0 = Min[Tac - Tav];
 egmax0 = Max[Tac - Tav];
 Print["\n***********************************\nEnergy band gap, Eg   =     ", N[egmin0,16]*efac, ";\nMax gap      =           ",
       N[egmax0*efac,16], "     eV"]
 
 
 
 
 (*******************************************)
 (********************************)
 (*Calculation of the CHERN No.*)
 
 chernNo    = -NIntegrate[
                            cBerryCurvature[kx, ky, t01, t02, M0, phi0, aax, aay, bbx, bby][[
                                                                                             3]], {kx, -kxMax, kxMax}, {ky, -kyMax, kyMax}  ]/(2 Pi)/2.;


 haldaneparam    = {a0, t01, t02, M0, phi0, T2, egmax0, egmin0, chernNo};
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
 
 
 
 
 LaserOut = Table[{j - 1, time[[j]], Ef[time[[j]], w00, e00, alpha0]}, {j, 1,
    Ntime}];
 Export[laserfile, LaserOut, "Data"];
 
 
 
 
 
 
 (**************************)
 (*Rabbi Freq. *)
 ROmega := Function[{kx, ky, t}
                    
                    , Ef[t, w00, e00, alpha0]*
                    Dipole[ kx, ky, t01, t02, M0, phi0, aax, aay, bbx, bby,
                           regular][[dir]]
                    
                    ];
 

 
 
 
 
(*******************************************)
(*******************************************)
(*******************************************)
(****
        Evol. Equations of SBEs       ****)
(*******************************************)
(***********************************************)
(*** Starting Do Loop on ky direction ***)
(***********************************************)
 
 
Do[
 
    
    
    yshift   = kyMin +  dky qindex; (**** Momentum on ky-direction ****)
    yindex   = qindex;
    
    Print["\n\n\n\n/*******************************************/\ny-index = ",yindex,";      Ky = ", yshift, "  a.u. \n/************/\n"];
    

 
  
 
 
   (*Dinamical evolution of populations and SBEs defined by Eqs 1 -- 3, \
    according to the constriction*)
   (*SlideView[*)
    solver = Function[ky,
                               
            NDSolve[
                                       
                    {
                                           
                                (*******************************************)
                                    (****
                                        Evol. Equations of SBEs       ****)
                                (*******************************************)
                                           \
                                           
                    (*Equation1*)
                                           
                    D[pi[kx, t], t] - Ef[t, w00, e00, alpha0]* D[pi[kx, t], kx] ==
                                    -I ( Eg[kx, ky, t01, t02, M0, phi0, aax, aay, bbx, bby]
                                               +
                                            Ef[t, w00, e00,
                                                  alpha0] Chig[kx, ky, t01, t02, M0, phi0, aax, aay, bbx,
                                                               bby, creg][[dir]]
                                            -  I/T2)* pi[kx, t]
                                        - I ROmega[kx, ky, t] (*pv[kx,t] - pc[kx,t]*)
                                           
                                
                                
                    (*Equation2*)
                                           ,
                    D[nc[kx, t], t] - Ef[t, w00, e00, alpha0] D[nc[kx, t], kx] ==
                                           -2.*Im[
                                                 Conjugate[ROmega[kx, ky, t]]
                                                 * pi[kx, t]
                                                 ]
                                           
                    
                        
                        
                                           
                    (*******************************************)
                            (****         \
                            Initials1             ****)
                    (*******************************************)
                        
                        , pi[kxMin, t] == pi[kxMax, t]
                        
                        , pi[kx, -tmax] == 0.
                                           
                        (*Initials2*)
                        , nc[kxMin, t] == nc[kxMax, t]
                        , nc[kx, -tmax] == 0.
                        (*******************************************)
                                           
                    }
                                       
                    , {pi, nc}
                                       
                    , {kx, kxMin, kxMax}
                                       
                    , {t, -tmax, tmax}
                                       
                    , AccuracyGoal -> 10, PrecisionGoal -> 10
                                       (*,Method\[Rule]"ImplicitRungeKutta",*)
                                       (*,Method\[Rule]
                                        RK5*)
                    ]
                               
    ];
    
    

 
 
    mdfun = solver[yshift];
             
    pi0 = Table[
                         Evaluate[ pi[kxx[[i]], time[[j*Nsnaper]] ] /. mdfun]
                         , {j, 1, Ntime/Nsnaper}, {i, 1, Nkx}
                         ];
             
             
    nc0 = Table[
                         Re[Evaluate[ nc[kxx[[i]], time[[j*Nsnaper]] ] /. mdfun] ]
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

             
             
             
             
    Repi2           = Re[pi1];
    ImPi2           = Im[pi1];
             



(* OUTPUTS *)
             
    fname0 = toNumberedFileName["CoherenceMathRealPiNx_"    , Nkx];
    fname0 = toNumberedFileName[fname0 <> "_kyIndex_"       , qindex];
    fname01 = toNumberedFileName["CoherenceMathImagPiNx_"   , Nkx];
    fname01 = toNumberedFileName[fname01 <> "_kyIndex_"     , qindex];
    fname02 = toNumberedFileName["Integ_Occupation_CMathNx_", Nkx];
    fname02 = toNumberedFileName[fname02 <> "_kyIndex_"     , qindex];
    fname03 = toNumberedFileName["Integ_Coherence_ReMathNx_", Nkx];
    fname03 = toNumberedFileName[fname03 <> "_kyIndex_"     , qindex];
    fname04 = toNumberedFileName["Integ_Coherence_ImMathNx_", Nkx];
    fname04 = toNumberedFileName[fname04 <> "_kyIndex_"     , qindex];
             
    fname1 = toNumberedFileName["OccupationCBMath_ncNx_"    , Nkx];
    fname1 = toNumberedFileName[fname1 <> "_kyIndex_"       , qindex];
             
    DataPath0   = rootPathA <> fname0   <> ".dat";
    DataPath01  = rootPathA <> fname01  <> ".dat";
    DataPath02  = rootPathA <> fname02  <> ".dat";
    DataPath03  = rootPathA <> fname03  <> ".dat";
    DataPath04  = rootPathA <> fname04  <> ".dat";
    DataPath1   = rootPathA <> fname1   <> ".dat";
             
             
(*OutPuts Coherence/Populations*)

    Export[DataPath0,   Repi2, "Data"];
    Export[DataPath01,  ImPi2, "Data"];
    Export[DataPath1,   nc1,   "Data"];
             
             
             
    (*
              
    Export[DataPath02,integOccupC,"Data"];
              Export[DataPath03,integCohRe,"Data"];
              Export[DataPath04,integCohIm,"Data"];
              
    *)
             


             
             
             

    (**************************************)
        (**DO--Conditions**)
    ,
             
    {qindex,0,1}
    (**************************************)
             
];
  
Print["\n\n*****************************************\n"]
Print["END OF THE PROGRAM AT, \nky = ",yshift, " a.u.; and\nyindex = ", yindex ]
Print["\n***********************************\n\n\n"]

