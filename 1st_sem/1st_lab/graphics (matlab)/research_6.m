y=@f2;
a=2;
b=3;

options = optimset('TolX', 1e-32, 'Display', 'iter');
f=vpa(fzero(y, (a+b)/2, options),32);
disp(f);
f=vpa(fzero(y, -35, options),32);
disp(f);

g=roots([1,0,-6,12,-8]);
disp(g);