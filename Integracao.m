clc
clear all
close all

function [y] = fx(x)
  y = 1 - exp(-x);
end

function [I] = TrapezioSimples(a, b)
  h = b - a;

  I = h * (fx(a) + fx(b))/2;
end

function [I] = TrapezioComposto(a, b, n)
  h = (b - a)/n;

  x =[];
  x(1) = a;
  S = fx(x(1));

  for i = 2:n
    x(i) = x(i-1) + h;
    S = S + 2*fx(x(i));
  endfor

  x(n+1) = x(n) + h;
  S = S + fx(x(n+1));
  I = h * S / 2;
end

function [I] = Simpson13(a, b)
  h = (b - a)/2;

  I = h * (fx(1) + 4*fx(2) + fx(3))/3;
end

function [I] = Simpson38(a, b)
  h = (b - a)/3;

  I = 3 * h * (fx(1) + 3*(fx(2) + fx(3)) + fx(4))/8;
end

a = 0;
b = 4;

[Trapezio] = TrapezioSimples(a, b)
[Simpson] = Simpson13(a, b)
[TrapezioComp] = TrapezioComposto(a, b, 4)
[Simpson2] = Simpson38(a, b)
