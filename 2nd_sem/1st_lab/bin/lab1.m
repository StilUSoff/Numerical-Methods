folder1 = ['matlab_res/'];

for i=1:9
    
    g=figure;
    tic;
    R{i}=sprandsym (500, i/10, 1,1); %генерация сим. пол. опред. матрицы
    R_time(i)=toc;
    %spy (R{i}, 'm')
    %saveas(g, [folder1 '(' int2str(i) ')R.png']);
    %close
    
    f=figure;
    tic;
    L{i}=chol(R{i}); %заполнение множителей разложения
    n{i}=nnz(L{i});
    L_time(i)=toc;
    spy(L{i}, 'm')
    saveas(f, [folder1 '(' int2str(i) ')_native_matrix_chol.png']);
    close
    close
    
    h=figure;
    tic;
    p=symrcm(R{i}); %Обычный алгоритм Катхилла-Макки; перестановка строк
    R_1{i}=R{i}(p,p);
    R_1_time(i)=toc;
    %spy(R_1{i}, 'm')
    %saveas(h, [folder1 '(' int2str(i) ')R_1.png']);
    %close
    
    l=figure;
    tic;
    L_1{i} =chol(R_1{i}); %заполнение множителей разложения
    n_1{i}=nnz(L_1{i});
    L_1_time(i)=toc;
    spy(L_1{i}, 'm')
    saveas(l, [folder1 '(' int2str(i) ')_symrcm_matrix_chol.png']);
    close
    close
    
    m=figure;
    tic;
    p=symamd(R{i}); %минимальная степень; перестановка столбцов
    R_2{i}=R{i}(p,p);
    R_2_time(i)=toc;
    %spy(R_2{i}, 'm')
    %saveas(m, [folder1 '(' int2str(i) ')R_2.png']);
    %close
    
    q=figure;
    tic;
    L_2{i} =chol(R_2{i}); %заполнение множителей разложения
    n_2{i}=nnz(L_2{i});
    L_2_time(i)=toc;
    spy(L_2{i}, 'm')
    saveas(q, [folder1 '(' int2str(i) ')_symamd_matrix_chol.png']);
    close
    close
    
end

x=1:i;

% ab=figure;
% plot(x,R_time,'-*',... %вермя на построение исходной матрицы
%     'LineWidth',3,...
%     'Color', '#AFB5FF',...
%     'MarkerSize',1,... 
%     'MarkerFaceColor','#AFB5FF');
% title('R', 'fontsize', 12);
% xlabel('i, м', 'fontsize', 12);
% ylabel('time, K', 'fontsize', 12);
% saveas(ab, [folder1 'native_time.png']);
% close

bc=figure;
plot(x,L_time,'-*',... %вермя на разложения исходной хол
    'LineWidth',3,...
    'Color', '#AFB5FF',...
    'MarkerSize',1,... 
    'MarkerFaceColor','#AFB5FF');
title('native chol time', 'fontsize', 12);
xlabel('i', 'fontsize', 12);
ylabel('time, c', 'fontsize', 12);
saveas(bc, [folder1 'native_chol_time.png']);
close
close

cb=figure;
plot(x,R_1_time,'-*',... %вермя на применение symrcm
    'LineWidth',3,...
    'Color', '#AFB5FF',...
    'MarkerSize',1,... 
    'MarkerFaceColor','#AFB5FF');
title('symrcm time', 'fontsize', 12);
xlabel('i', 'fontsize', 12);
ylabel('time, c', 'fontsize', 12);
saveas(cb, [folder1 'symrcm_time.png']);
close
close

abc=figure;
plot(x,L_1_time,'-*',... %вермя на разложения symrcm хол
    'LineWidth',3,...
    'Color', '#AFB5FF',...
    'MarkerSize',1,... 
    'MarkerFaceColor','#AFB5FF');
title('symrcm chol time', 'fontsize', 12);
xlabel('i', 'fontsize', 12);
ylabel('time, i', 'fontsize', 12);
saveas(abc, [folder1 'symrcm_chol_time.png']);
close
close

abdf=figure;
plot(x,R_2_time,'-*',... %вермя на применение symand
    'LineWidth',3,...
    'Color', '#AFB5FF',...
    'MarkerSize',1,... 
    'MarkerFaceColor','#AFB5FF');
title('symamd time', 'fontsize', 12);
xlabel('i', 'fontsize', 12);
ylabel('time, c', 'fontsize', 12);
saveas(abdf, [folder1 'symamd_time.png']);
close
close

abqr=figure;
plot(x,L_2_time,'-*',... %вермя на разложения symand хол
    'LineWidth',3,...
    'Color', '#AFB5FF',...
    'MarkerSize',1,... 
    'MarkerFaceColor','#AFB5FF');
title('symamd chol time', 'fontsize', 12);
xlabel('i', 'fontsize', 12);
ylabel('time, c', 'fontsize', 12);
saveas(abqr, [folder1 'symamd_chol_time.png']);
close
close

clear