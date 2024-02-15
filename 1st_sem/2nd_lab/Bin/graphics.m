function graphics(number)

    for i=5:2:9
        folder = ['bin/func' int2str(number) '/res/'];
        folders = ['bin/graphics/func' int2str(number) '/'];

        A = load([folder int2str(i) 'splres.txt']);
        B = load([folder int2str(i) 'chebres.txt']);
        C = load([folder int2str(i) 'lineres.txt']);

    if number==1
        a=0.5;
        b=2.5;
    elseif number==2
        a=-0.51;
        b=0.49;
    end
    x_nodes=zeros(i,1);
    cheb_nodes=zeros(i,1);
    h=abs(b-a)/(i-1);
    for k=1:1:i
        x_nodes(k)=a+(k-1)*h;
        cheb_nodes(k)=(a+b)/2+(b-a)*cos((2*k-1)/(2*i)*pi)/2;
    end

    sz=40;

    f=figure;   
    x = A(:,1);
    y1 = A(:,2);
    plot(x,y1,'-',...
        'LineWidth',4,...
        'Color',"#7decf0");
    title(['spline' ' ' int2str(i) ' ' 'nodes'], 'FontSize', 12);
    xlabel('x', 'FontSize', 12);
    ylabel('y(x)', 'FontSize', 12);
    hold on; grid on;
    if number==1
        x = 0.5:0.01:2.5;
        y2=atan(x)+0.25*x.^(-3);
        y3=atan(x_nodes)+0.25*x_nodes.^(-3);
        ylim([0,2.465]);
    elseif number==2
        x=-0.51:0.004:0.49;
        y2=2*x.^3-x.^2-60*abs(x)+1;
        y3=2*x_nodes.^3-x_nodes.^2-60*abs(x_nodes)+1;
        ylim([-30.6,1.2]);
    end 
	plot(x,y2, 'Color', '#d16608', 'LineWidth',1);
    scatter(x_nodes,y3,sz,'d',...
        'MarkerEdgeColor',"#530FAD",...
        'MarkerFaceColor',"#530FAD",...
        'LineWidth',1.5);
    legend(['Spline' ' ' int2str(i)], 'Actual', 'nodes');
    hold off;
    saveas(f, [folders int2str(i) 'spl.png']);
    close
    h=figure;
    x = B(:,1);
    y1 = B(:,2);
    plot(x,y1,'-',...
        'LineWidth',4,...
        'Color',"#7decf0");
    title(['chebyshev' ' ' int2str(i) ' ' 'nodes'], 'FontSize', 12);
    xlabel('x', 'FontSize', 12);
    ylabel('y(x)', 'FontSize', 12);
    hold on; grid on;
    if number==1
        x = 0.5:0.01:2.5;
        y2= atan(x)+0.25*x.^(-3);
        y3 = atan(cheb_nodes)+0.25*cheb_nodes.^(-3);
        ylim([0,2.46]);
    elseif number==2
        x=-0.51:0.004:0.49;
        y2=2*x.^3-x.^2-60*abs(x)+1;
        y3=2*cheb_nodes.^3-cheb_nodes.^2-60*abs(cheb_nodes)+1;
        ylim([-30.6,1.2]);
    end 
    plot(x,y2, 'Color', '#d16608', 'LineWidth',1);
    scatter(cheb_nodes,y3,sz,'d',...
        'MarkerEdgeColor',"#530FAD",...
        'MarkerFaceColor',"#530FAD",...
        'LineWidth',1.5);
    legend('LagrangeChebyshev', 'Actual', 'nodes');
    hold off;
    saveas(h,  [folders int2str(i) 'cheb.png']);
    close
    g=figure;
    x = C(:,1);
    y1 = C(:,2);
    plot(x,y1,'-',...
        'LineWidth',4,...
        'Color',"#7decf0");
    title(['linear' ' ' int2str(i) ' ' 'nodes'], 'FontSize', 12);
    xlabel('x', 'FontSize', 12);
    ylabel('y(x)', 'FontSize', 12);
    hold on; grid on;
    if number==1
        x = 0.5:0.01:2.5;
        y2=atan(x)+0.25*x.^(-3);
        y3=atan(x_nodes)+0.25*x_nodes.^(-3);
        ylim([0,2.465]);
    elseif number==2
        x=-0.51:0.004:0.49;
        y2=2*x.^3-x.^2-60*abs(x)+1;
        y3=2*x_nodes.^3-x_nodes.^2-60*abs(x_nodes)+1;
        ylim([-30.6,1.2]);
    end 
    plot(x,y2, 'Color', '#d16608', 'LineWidth',1);
    scatter(x_nodes,y3,sz,'d',...
        'MarkerEdgeColor',"#530FAD",...
        'MarkerFaceColor',"#530FAD",...
        'LineWidth',1.5);
    legend('LagrangeLinear', 'Actual', 'nodes');
    hold off;
    saveas(g,  [folders int2str(i) 'line.png']);
    close
    e=figure;
    x=A(:,1);
    if number==1
        y1=atan(x)+0.25*x.^(-3);
    elseif number==2
        y1=2*x.^3-x.^2-60*abs(x)+1;
    end
    y2=abs(A(:,2)-y1);
    y3=abs(B(:,2)-y1);
    y4=abs(C(:,2)-y1);
    alpha=(a+b)/2;
    if i==5
        if number==1
            y5=(x-x_nodes(1)).*(x-x_nodes(2)).*(x-x_nodes(3)).*(x-x_nodes(4)).*(x-x_nodes(5))*(((atan(alpha)+0.25*alpha.^(-3)).^6)/(factorial(6)));
        elseif number==2
            y5=(x-x_nodes(1)).*(x-x_nodes(2)).*(x-x_nodes(3)).*(x-x_nodes(4)).*(x-x_nodes(5))*(((2*alpha.^3-alpha.^2-60*abs(alpha)+1).^6)/(factorial(6)));
        end            
    elseif i==7
        if number==1
            y5=(x-x_nodes(1)).*(x-x_nodes(2)).*(x-x_nodes(3)).*(x-x_nodes(4)).*(x-x_nodes(5)).*(x-x_nodes(6)).*(x-x_nodes(7))*(((atan(alpha)+0.25*alpha.^(-3)).^8)/(factorial(8)));
        elseif number==2
            y5=(x-x_nodes(1)).*(x-x_nodes(2)).*(x-x_nodes(3)).*(x-x_nodes(4)).*(x-x_nodes(5)).*(x-x_nodes(6)).*(x-x_nodes(7))*(((2*alpha.^3-alpha.^2-60*abs(alpha)+1).^8)/(factorial(8)));
        end            
    else
        if number==1
            y5=(x-x_nodes(1)).*(x-x_nodes(2)).*(x-x_nodes(3)).*(x-x_nodes(4)).*(x-x_nodes(5)).*(x-x_nodes(6)).*(x-x_nodes(7)).*(x-x_nodes(8)).*(x-x_nodes(9))*(((atan(alpha)+0.25*alpha.^(-3)).^10)/(factorial(10)));
        elseif number==2
            y5=(x-x_nodes(1)).*(x-x_nodes(2)).*(x-x_nodes(3)).*(x-x_nodes(4)).*(x-x_nodes(5)).*(x-x_nodes(6)).*(x-x_nodes(7)).*(x-x_nodes(8)).*(x-x_nodes(9))*(((2*alpha.^3-alpha.^2-60*abs(alpha)+1).^10)/(factorial(10)));
        end            
    end
    plot(x,y2, 'Color', '#d16608', 'LineWidth',1);
    title(['error' ' ' int2str(i) ' ' 'nodes'], 'FontSize', 12);
    xlabel('x', 'FontSize', 12);
    ylabel('ln(error)', 'FontSize', 12);
    hold on; grid on;
    plot(x,y3, 'Color', '#7decf0', 'LineWidth',1);
    plot(x,y4, 'Color',  '#530FAD', 'LineWidth',1);
    plot(x,y5, 'Color',  '#FA7080', 'LineWidth',1);
    legend('Spline', 'LagrCheb', 'LagrLine', 'Theor error');
    saveas(e, [folders int2str(i) 'error(x).png']);
    close
    end
    
    
    
    folder1 = ['bin/func' int2str(number) '/norm/'];
    folder2 = ['bin/func' int2str(number) '/fluct/'];
    folder3 = ['bin/func' int2str(number) '/thick/'];
        
    J = load([folder1 'splnorm.txt']);
    K = load([folder1 'chebnorm.txt']);
    L = load([folder1 'linenorm.txt']);
        
    D = load([folder2 'splfluct.txt']);
    E = load([folder2 'chebfluct.txt']);
    F = load([folder2 'linefluct.txt']);
        
    G = load([folder3 'splthick.txt']);
    H = load([folder3 'chebthick.txt']);
    I = load([folder3 'linethick.txt']);
    
    j=figure;
    x = J(:,1);
    y = J(:,2);
    loglog(x,y,'-bo',...
        'LineWidth',1,...
        'Color',"#d16608",...
        'MarkerSize',3,... 
        'MarkerFaceColor',"#d16608");
    title('spline norm', 'fontsize', 12);
    xlabel('n', 'fontsize', 12);
    ylabel('norm', 'fontsize', 12);
    hold off;
    saveas(j, [folders 'splnorm.png']);
    close
    
    
    k=figure;
    x = K(:,1);
    y = K(:,2);
    loglog(x,y,'-bo',...
        'LineWidth',1,...
        'Color',"#d16608",...
        'MarkerSize',3,... 
        'MarkerFaceColor',"#d16608");
    title('cheb norm', 'fontsize', 12);
    xlabel('n', 'fontsize', 12);
    ylabel('norm', 'fontsize', 12);
    hold off;
    saveas(k, [folders 'chebnorm.png']);
    close
    l=figure;
    x = L(:,1);
    y = L(:,2);
    loglog(x,y,'-bo',...
        'LineWidth',1,...
        'Color',"#d16608",...
        'MarkerSize',3,... 
        'MarkerFaceColor',"#d16608");
    title('linear norm', 'fontsize', 12);
    xlabel('n', 'fontsize', 12);
    ylabel('norm', 'fontsize', 12);
    hold off;
    saveas(l, [folders 'linenorm.png']);
    close
    jj=figure;
    x = J(:,1);
    y1 = J(:,3);
    title('spline error', 'fontsize', 12);
    xlabel('n', 'fontsize', 12);
    ylabel('err', 'fontsize', 12);
    hold on;
    loglog(x,y1,'-bo',...
        'LineWidth',1,...
        'Color',"#d16608",...
        'MarkerSize',3,... 
        'MarkerFaceColor',"#d16608");
    y2 = J(:,4);
    hold on;
    loglog(x,y2,'-bo',...
        'LineWidth',1,...
        'Color',"#7decf0",...
        'MarkerSize',3,... 
        'MarkerFaceColor',"#7decf0");
    if number==1
        legend(['in x=0.79985'], ['in x=1.09985']);
    else
        legend(['in x=-0.3001'], ['in x=-0.1001']);
    end
    hold off;
    saveas(jj, [folders 'splerr.png']);
    close
    kk=figure;
    x = K(:,1);
    y1 = K(:,3);
    title('chebyshev error', 'fontsize', 12);
    xlabel('n', 'fontsize', 12);
    ylabel('err', 'fontsize', 12);
    hold on;
    loglog(x,y1,'-bo',...
        'LineWidth',1,...
        'Color',"#d16608",...
        'MarkerSize',3,... 
        'MarkerFaceColor',"#d16608");
    y2 = K(:,4);
    hold on;
    loglog(x,y2,'-bo',...
        'LineWidth',1,...
        'Color',"#7decf0",...
        'MarkerSize',3,... 
        'MarkerFaceColor',"#7decf0");
    if number==1
        legend(['in x=0.79985'], ['in x=1.09985']);
    else
        legend(['in x=-0.3001'], ['in x=-0.1001']);
    end
    hold off;
    saveas(kk, [folders 'cheberr.png']);
    close
    ll=figure;
    x = L(:,1);
    y1 = L(:,3);
    title('linear error', 'fontsize', 12);
    xlabel('n', 'fontsize', 12);
    ylabel('err', 'fontsize', 12);
    hold on;
    loglog(x,y1,'-bo',...
        'LineWidth',1,...
        'Color',"#d16608",...
        'MarkerSize',3,... 
        'MarkerFaceColor',"#d16608");
    y2 = L(:,4);
    hold on;
    plot(x,y2,'-bo',...
        'LineWidth',1,...
        'Color',"#7decf0",...
        'MarkerSize',3,... 
        'MarkerFaceColor',"#7decf0");
    if number==1
        legend(['in x=0.79985'], ['in x=1.09985']);
    else
        legend(['in x=-0.3001'], ['in x=-0.1001']);
    end
    hold off;
    saveas(ll, [folders 'lineerr.png']);
    close
    m=figure;
    x = D(:,1);
    y = D(:,2);
    scatter(x,y, "filled","markerfacecolor","#d16608");
    title('spline fluct', 'fontsize', 12);
    xlabel('fluct', 'fontsize', 12);
    ylabel('%', 'fontsize', 12);
    hold off;
    saveas(m, [folders 'splfluct.png']);
    close
    n=figure;
    x = E(:,1);
    y = E(:,2);
    scatter(x,y, "filled","markerfacecolor","#d16608");
    title('chebyshev fluct', 'fontsize', 12);
    xlabel('fluct', 'fontsize', 12);
    ylabel('%', 'fontsize', 12);
    hold off;
    saveas(n, [folders 'chebfluct.png']);
    close
    o=figure;
    x = F(:,1);
    y = F(:,2);
    scatter(x,y, "filled","markerfacecolor","#d16608");
    title('linear fluct', 'fontsize', 12);
    xlabel('fluct', 'fontsize', 12);
    ylabel('%', 'fontsize', 12);
    hold off;
    saveas(o, [folders 'linefluct.png']);
    close
    p=figure;
    x = G(:,1);
    y = G(:,2);
    loglog(x,y,'-bo',...
        'LineWidth',1,...
        'Color',"#d16608",...
        'MarkerSize',3,... 
        'MarkerFaceColor',"#d16608")
    title('spline thick', 'fontsize', 12);
    xlabel('n', 'fontsize', 12);
    ylabel('norm', 'fontsize', 12);
    hold off;
    saveas(p, [folders 'splthick.png']);
    close
    q=figure;
    x = H(:,1);
    y = H(:,2);
    loglog(x,y,'-bo',...
        'LineWidth',1,...
        'Color',"#d16608",...
        'MarkerSize',3,... 
        'MarkerFaceColor',"#d16608");
    title('chebyshev thick', 'fontsize', 12);
    xlabel('n', 'fontsize', 12);
    ylabel('norm', 'fontsize', 12);
    hold off;
    saveas(q, [folders 'chebthick.png']);
    close
    r=figure;
    x = I(:,1);
    y = I(:,2);
    loglog(x,y,'-bo',...
        'LineWidth',1,...
        'Color',"#d16608",...
        'MarkerSize',3,... 
        'MarkerFaceColor',"#d16608");
    title('linear thick', 'fontsize', 12);
    xlabel('n', 'fontsize', 12);
    ylabel('norm', 'fontsize', 12);
    hold off;
    saveas(r, [folders 'linethick.png']);
    close
end