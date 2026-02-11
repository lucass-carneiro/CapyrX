#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Symbols from ADMBaseX*)


(* ::Text:: *)
(*We follow Eqs.(56)-(58) of https://arxiv.org/abs/gr-qc/0507004*)


ClearAll[\[Alpha], \[Beta]]
\[Alpha][][] = alp;
\[Beta][i_][] := {betax, betay, betaz}[[i]];

ClearAll[\[Gamma]Mat, \[Gamma]Imat, \[Gamma]];
\[Gamma]Mat={{gxx, gxy, gxz},{gxy, gyy, gyz}, {gxz, gyz, gzz}};
\[Gamma]iMat=Inverse[\[Gamma]Mat]//FullSimplify;

\[Gamma][][i_,j_] := \[Gamma]Mat[[i,j]];
\[Gamma][i_,j_][] := \[Gamma]iMat[[i,j]];
\[Gamma][][] = Det[\[Gamma]Mat]//FullSimplify;


(* ::Section:: *)
(*Wave equation symbols*)


ClearAll[d];
d[][i_] := {dx, dy, dz}[[i]];

ClearAll[H];
H[i_,j_][] := FullSimplify[\[Gamma][i, j][] - (\[Alpha][][])^(-2) \[Beta][i][]\[Beta][j][]];

ClearAll[F];
F[i_][] := FullSimplify[(\[Gamma][][])^(1/2) \[Beta][i][] \[CapitalPi] + \[Alpha][][](\[Gamma][][])^(1/2) Sum[H[i, j][] d[][j], {j, 1, 3}]];


(* ::Section:: *)
(*Right Hand Sides*)


ClearAll[\[Phi]RHS, \[CapitalPi]RHS, dRHS];
\[Phi]RHS = \[Alpha][][] \[CapitalPi];
\[CapitalPi]RHS = FullSimplify[Sum[\[Beta][i][](\[Alpha][][])^(-1) deriv[FluxA, i] + (\[Gamma][][])^(-1/2) deriv[FluxB[i][], i],{i, 1, 3}]];
dRHS[][i_] = deriv[FluxA, i];


(* ::Section:: *)
(*Code generation*)


(* ::Subsection:: *)
(*Formatting*)


Unprotect[Power];
Format[Power[E, a_], CForm] := exp[a];
Format[Power[a_, 1/2], CForm] := sqrt[a];
Format[Power[a_, b_], CForm] := pow[a,b];
Protect[Power];


(* ::Subsection:: *)
(*Flux B*)


ClearAll[rules];
rules={
  "gxx"->"gxx(pI)",
  "gxy"->"gxy(pI)",
  "gxz"->"gxz(pI)",
  "gyy"->"gyy(pI)",
  "gyz"->"gyz(pI)",
  "gzz"->"gzz(pI)",

  "alp"->"alp(pI)",

  "betax"->"betax(pI)",
  "betay"->"betay(pI)",
  "betaz"->"betaz(pI)",
  
  "\[CapitalPi]" -> "Pi(pI)",
  "dx" -> "Dx(pI)",
  "dy" -> "Dy(pI)",
  "dz" -> "Dz(pI)"
};

Print[
  "const auto flux_b_x = [&](const vect<int, dim> &pI) { return ",
  StringReplace[ToString[F[1][],CForm], rules],
  ";};\n"
]

Print[
  "const auto flux_b_y = [&](const vect<int, dim> &pI) { return ",
  StringReplace[ToString[F[2][],CForm], rules],
  ";};\n"
]

Print[
  "const auto flux_b_z = [&](const vect<int, dim> &pI) { return ",
  StringReplace[ToString[F[3][],CForm], rules],
  ";};"
]

ClearAll[rules];


(* ::Subsection:: *)
(*RHS*)


ClearAll[rules];
rules={
  "gxx"->"gxx(p.I)",
  "gxy"->"gxy(p.I)",
  "gxz"->"gxz(p.I)",
  "gyy"->"gyy(p.I)",
  "gyz"->"gyz(p.I)",
  "gzz"->"gzz(p.I)",

  "alp"->"alp(p.I)",

  "betax"->"betax(p.I)",
  "betay"->"betay(p.I)",
  "betaz"->"betaz(p.I)",
  
  "deriv(FluxA,1)" -> "g_dflux_a.dx",
  "deriv(FluxA,2)" -> "g_dflux_a.dy",
  "deriv(FluxA,3)" -> "g_dflux_a.dz",
  
  "deriv(FluxB(1)(),1)" -> "g_dflux_b_x.dx",
  "deriv(FluxB(2)(),2)" -> "g_dflux_b_y.dy",
  "deriv(FluxB(3)(),3)" -> "g_dflux_b_z.dz"
};

Print["Pi_rhs(p.I) = ", StringReplace[ToString[\[CapitalPi]RHS, CForm], rules],";"]

ClearAll[rules];
