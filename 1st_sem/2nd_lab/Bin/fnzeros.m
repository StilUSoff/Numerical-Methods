g=-0.02:0.001:0.02;
h=2*g.^3-g.^2-60*abs(g)+1
s1=csape(g,h);
j=fnder(s1);
f=fnzeros(j);

m=2*f.^3-f.^2-60*abs(f)+1

i=figure;
plot(g,h, 'Color','#d16608');
ylim([-0.2;1.2])
xlim([-0.5,0.5]);
hold on; grid on;
plot(f,m,'c*','MarkerSize',10, 'Color','#530FAD');
hold off;