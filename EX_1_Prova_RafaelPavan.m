clc
clear all
close all

% Rafael Pavan
% Exercício 1 - Parte A da Prova
% Disciplina - Otimização e Método de Pontos Interiores

% Método Primal Dual Barreira Logarítmica Previsor-Corretor Mehrotra

%% Inicializa dados do problema

A = [-1 1 1 0 0 0;-1 -1 0 1 0 0;1 0 0 0 1 0; 0 1 0 0 0 1];

b = [1 ; -1; 3; 2];

c = [4;-6;0;0;0;0];

Q = [4 0 0 0 0 0; 0 6 0 0 0 0; 0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0];

e = ones(1,6);

iteracoes=1;

maxi = 4;

Xk = zeros(6,6);
Sk = zeros(6,6);

%% Inicializa Ponto Inicial

xk = [2.75;0.75;1;1;1.5;1];

wk = [0;0;0;0];

sk = [0.5;0.5;0.5;0.5;0.5;0.5;];

xx = [xk];

%% Define Tolerâncias

episolon1 = 9.99e-2;
episolon2 = 9.99e-2;
episolon3 = 9.99e-2;


f_obj =[transpose(xk)*Q*xk+transpose(c)*xk];

%% Teste Inicial

display('Teste Inicial: x0>0?')

x0 = -inv(Q(1:2,1:2))*c(1:2);

if sum(x0>0)<length(x0)
    
    display('Solução do Teste Inicial Não Foi Factível')
else
    display('Solução do Teste Inicial Foi Factível')
end

display('Teste Inicial: A*x0-b = 0 ?')

if sum(A(1:2,1:2)*x0-b(1:2)==0)==length(x0)
    
    display('Solução do Teste Inicial Não Foi Factível')

else

    display('Solução do Teste Inicial Não Foi Factível')

end

display('____________________')


%% Algoritmo

while(iteracoes<maxi)
    
    for i=1:length(xk)
        
        Xk(i,i) = xk(i);
        
    end
    
    
    for i=1:length(sk)
        
        Sk(i,i) = sk(i);
        
    end
    
    
%% Cálculo de mik, uk, tk, vk, pk, thetak
    
    mik = transpose(xk)*sk/length(sk); 

    if iteracoes == 1
         
         mik = 1;
         
    end
        
    tk =  b - A*xk;
    
    uk = Q*xk + c - transpose(A)*wk-sk;
    
    vk = mik*(transpose(e)) - Xk*Sk*transpose(e);
    
    pk = inv(Xk)*vk;
    
    thetak = inv(Q+inv(Xk)*Sk);
    
%% Teste de Otimalidade

   if mik < episolon1 & norm(tk)/(norm(b)+1) < episolon2 & norm(uk)/(norm(Q*xk+c)+1) < episolon3
        
        display('Solução Ótima Encontrada')
        
        break
        
   end
        
%% Passo Previsor    (mik=0)
    
 
    dwk_p = inv(A*thetak*transpose(A))*(tk + A*thetak*(uk));
    
    dxk_p = thetak*(transpose(A)*dwk_p-uk);
    
    dsk_p = -inv(Xk)*Sk*dxk_p;
    
       
    alfa = 0.9995;
   
    
    betap1 = min([1;-alfa*xk(dxk_p<0)./dxk_p(dxk_p<0)]);
    
    betad = min([1;-alfa*sk(dsk_p<0)./dsk_p(dsk_p<0)]);
 
    betaq = -transpose(dxk_p)*(Q*xk+c)/(transpose(dxk_p)*Q*dxk_p);
    
    betap=betap1;
    
    if betaq>0

        betap = min(betap1,betaq);
    
    end
       
   
    xk = xk + betap*dxk_p;
    wk = wk + betad*dwk_p;
    sk = sk + betad*dsk_p; 
      
%% Passo Corretor    

    mik_p = transpose(xk)*sk/length(sk);
    
    Dsk_p = zeros(length(dsk_p),length(dsk_p));
    
    Dxk_p = zeros(length(dxk_p),length(dxk_p));
    
    for i=1:length(dxk_p)
        
        Dxk_p(i,i) = dxk_p(i);
        
    end
    
    
    for i=1:length(dsk_p)
        
        Dsk_p(i,i) = dsk_p(i);
        
    end
    
   
    vk_c = mik_p*transpose(e) - Xk*Sk*transpose(e) - Dxk_p*Dsk_p*transpose(e);
   
    pk=inv(Xk)*vk_c;
    
    dwk = inv(A*thetak*transpose(A))*(tk + A*thetak*(uk-pk));
        
    dxk = thetak*(transpose(A)*dwk-uk+pk);
    
    dsk = inv(Xk)*(vk_c-Sk*dxk);
    
    alfa = 0.9995;
    
    betap1 = min([1;-alfa*xk(dxk<0)./dxk(dxk<0)]);
    
    betad = min([1;-alfa*sk(dsk<0)./dsk(dsk<0)]);
        
    betaq = -transpose(dxk)*(Q*xk+c)/(transpose(dxk)*Q*dxk);
   
    betap=betap1;
    
    if betaq>0
    
        betap = min(betap1,betaq);
    
    end
    
    
    
    % Teste de Otimalidade e Ilimitariedade
    
    if sum(tk==0) == length(tk) & sum(dxk>0)==length(dxk) & transpose(c)*dxk<0
        
        display('Problema Ilimitado')
        break 
    
    end
    
    if uk <= 0.09 & dsk>0 & transpose(b)*dwk>0

        display('Solução Ótima Encontrada')
        break 
    end
    

    
    % Atualização das Variáveis
    
   
    xk = xk + betap*dxk;
    wk = wk + betad*dwk;
    sk = sk + betad*dsk; 
    
    xx = [xx,xk];

    iteracoes = iteracoes + 1;
    
    f_obj = [f_obj,transpose(xk(1:2))*Q(1:2,1:2)*xk(1:2)/2+transpose(c(1:2))*xk(1:2)];
    
    
end

figure

plot(f_obj)
title('Curva de Minimização da Função Objetivo')
xlabel('Iterações')
ylabel('Valor da Função Objetivo')
grid on
grid minor

display(['Iteração Final: ' num2str(iteracoes-1)])

display(['Norma do Resíduo Primal tk: ' num2str(norm(tk))])

display(['Norma do Resíduo Dual uk: ' num2str(norm(uk))])

display(['Norma do Resíduo de Complementaridade vk: ' num2str(max(abs(vk)))])

display(['Valor final do parâmetro de barreira: ' num2str(mik)])

display(['Valor da Função Objetivo: ' num2str(f_obj(iteracoes))])

display('Solução Primal:')
display(xk)
    
display('Solução Dual:')
display(sk)
display(wk)
