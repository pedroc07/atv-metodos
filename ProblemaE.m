clc
clear all
close all

function [mX] = Aproximacao(mA, x0)
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

% IMPLEMENTACAO DO METODO DE GAUSS SEIDEL
function [xi] = GaussSeidel(mA, epsilon)
  dr = 1;
  x0 = [0, 0, 0];
  xi = x0;
  k = 0;
  while epsilon < dr
      aux = xi;
      k = k + 1
      xi = Aproximacao(mA, xi);
      dr = Calcula_dr(xi, aux);
  endwhile
end

% MATRIZ DAS EQUACOES LINEARES DO PROBLEMA E
mA = [-132, 22, 0, -1000; 5, -27, 7, -2000; 117, 0, -22, 110];
epsilon = 0.01;
mX = GaussSeidel(mA, epsilon)
