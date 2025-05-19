clc
clear all
close all

function [mX] = AproximacaoJacobi(mA, x0)
  mX(1) = (mA(1,2) * -1 * x0(2)) + (mA(1,3) * -1 * x0(3)) + mA(1,4);
  mX(2) = (mA(2, 1) * -1 * x0(1)) + (mA(2,3) * -1 * x0(3)) + mA(2,4);
  mX(3) = (mA(3,1) * -1 * x0(1)) + (mA(3,2) * -1 * x0(2)) + mA(3,4);

  mX(1) = mX(1)/mA(1,1);
  mX(2) = mX(2)/mA(2,2);
  mX(3) = mX(3)/mA(3,3);
end

function [mX] = AproximacaoSeidel(mA, x0)
  mX(1) = (mA(1,2) * -1 * x0(2)) + (mA(1,3) * -1 * x0(3)) + mA(1,4);
  mX(1) = mX(1)/mA(1,1);
  mX(2) = (mA(2, 1) * -1 * mX(1)) + (mA(2,3) * -1 * x0(3)) + mA(2,4);
  mX(2) = mX(2)/mA(2,2);
  mX(3) = (mA(3,1) * -1 * mX(1)) + (mA(3,2) * -1 * mX(2)) + mA(3,4);

  mX(3) = mX(3)/mA(3,3);
end

function [dr] = Calcula_dr(mX1, mX0)
    dr = max(abs(mX1 .- mX0)) / max(abs(mX1));
end

function [Resp] = DiagDominante(mA, n)
  Resp = 1;
  for i = 1: n
    soma = 0;
    for j = 1: n
      if j != i
        soma = soma + abs(mA(i,j));
      endif
    endfor
    if abs(mA(i,i)) <= soma
      Resp = 0;
    endif
  endfor
end

% IMPLEMENTACAO DO METODO DE GAUSS SEIDEL
function [xi] = GaussSeidel(mA, epsilon)
  dr = 1;
  x0 = [0, 0, 0];
  xi = x0;
  k = 0;
  while epsilon < dr
      aux = xi;
      k = k + 1
      xi = AproximacaoSeidel(mA, xi);
      dr = Calcula_dr(xi, aux);
  endwhile
end

% IMPLEMENTACAO DO METODO DE GAUSS JACOBI
function [xi] = GaussJacobi(mA, epsilon)
  dr = 1;
  x0 = [0, 0, 0];
  xi = x0;
  k = 0;
  converg = DiagDominante(mA, 3);
  if converg = 0
    while epsilon < dr
      aux = xi;
      k = k + 1
      xi = AproximacaoJacobi(mA, xi);
      dr = Calcula_dr(xi, aux);
    endwhile
  endif
end

% MATRIZ DAS EQUACOES LINEARES DO PROBLEMA E
mA = [-132, 22, 0, -1000; 5, -27, 7, -2000; 117, 0, -22, 110];
epsilon = 0.01;
mX = GaussSeidel(mA, epsilon)
mX = GaussJacobi(mA, epsilon)
