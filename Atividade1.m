clc
clear all
close all

function [y] = fx(x)
  y = 1 - exp(-x);
end

function [y] = fr(r, r0, n)
  y = r*(1 - r/r0)^(1/n);
end

function [I] = TrapezioSimples(a, b, n)
  h = b - a;

  I = h * (fr(a, b, n) + fr(b, b, n))/2;
end

function [I] = TrapezioComposto(a, b, n, num)
  h = (b - a)/n;

  x =[];
  x(1) = a;
  S = fr(x(1), b, num);

  for i = 2:n
    x(i) = x(i-1) + h;
    S = S + 2*fr(x(i), b, num);
  endfor

  x(n+1) = x(n) + h;
  S = S + fr(x(n+1), b, num);
  I = h * S / 2;
end

function [I] = Simpson13(a, b, n)
  h = (b - a)/2;

  I = h * (fr(1, b, n) + 4*fr(2, b, n) + fr(3, b, n))/3;
end

function [I] = Simpson38(a, b, n)
  h = (b - a)/3;

  I = 3 * h * (fr(1, b, n) + 3*(fr(2, b, n) + fr(3, b, n)) + fr(4, b, n))/8;
end

r0 = 0.5;
n = 8;
a = 0;

[Trapezio] = TrapezioSimples(a, r0, n)
[Simpson] = Simpson13(a, r0, n)
[TrapezioComp] = TrapezioComposto(a, r0, 4, n)
[Simpson2] = Simpson38(a, r0, n)
