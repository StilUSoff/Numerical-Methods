function csaperes
!sin
g=figure;
k=0;
for j=0.01:0.01:0.2
    k=k+1;
    x=0:j:pi;
    y=sin(x);
    s=csape(x,y);
    for i=1:size(x)
        norm(i)=max(abs(s.coefs(1,i)*x(i)^3+s.coefs(2,i)*x(i)^2+...
            s.coefs(3,i)*x(i)+s.coefs(4,i)-y(i)));
    end
    norm_max(k)=max(norm);
end
j=0.01:0.01:0.2;
plot(j,norm_max,'Color', '#d16608', 'LineWidth',2)
title('norm sin(x)', 'FontSize', 12);
xlabel('h', 'FontSize', 12);
ylabel('norm', 'FontSize', 12);
grid on;
saveas(g, 'graphics/sin(x)', 'png')
!cos
c=figure;
k=0;
for j=0.01:0.01:0.2
    k=k+1;
    x_1=0:j:pi;
    y_1=cos(x_1);
    s_1=csape(x_1,y_1);
    for i=1:size(x_1)
        norm(i)=max(abs(s_1.coefs(1,i)*x_1(i)^3+s_1.coefs(2,i)*x_1(i)^2+...
            s_1.coefs(3,i)*x_1(i)+s_1.coefs(4,i)-y_1(i)));
    end
    norm_max(k)=max(norm);
end
j=0.01:0.01:0.2;
plot(j,norm_max,'Color', '#d16608', 'LineWidth',2)
title('norm cos(x)', 'FontSize', 12);
xlabel('h', 'FontSize', 12);
ylabel('norm', 'FontSize', 12);
grid on;
saveas(c, 'graphics/cos(x)', 'png')

end

