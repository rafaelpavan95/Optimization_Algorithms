clc
clear all
close all


% Rafael Pavan
% Exercício 2 - Parte A da Prova
% Disciplina - Otimização e Método de Pontos Interiores

% Método Primal Dual Barreira Logarítmica para PNL

%% Dados Iniciais

x0 = [3;3;3;3;3];  

s0 = [1;1;1;1;1]; 

syms x1 x2 x3 x4 x5

eq = 100*(x1-x2*x2)*(x1-x2*x2) + (1-x1)*(1-x1) + 0*x3 + 0*x4 + 0*x5;

gradiente=gradient(eq,[x1,x2,x3,x4,x5]);

hessiana = hessian(eq,[x1,x2,x3,x4,x5]);

x1 = x0(1);
x2 = x0(2);
x3 = x0(3);
x4 = x0(4);
x5 = x0(5);

gradientevpa = vpa(subs(gradiente));

hessianavpa = vpa(subs(hessiana));

% Inicia Valor de Alfa:

alfa = 0.9995;

% Define a Matriz A:

A = [1, -1, 1, 0, 0; 0, -1, 0, 1, 0; 0, 1, 0, 0, 1];

% Define Vetor b:

b = [15; -2; 15];

% Define valor de tolerancia episolon:

e1 = 9.9999e-2;
e2 = 9.9999e-2;
e3 = 9.9999e-2;

% Define matriz Xk diagonal:

Xk = zeros(length(x0),length(x0));

Sk = zeros(length(s0),length(s0));

% Inicializa vetor de e [1,1,1....,1]

e = ones(1,length(x0));

maxi = 5;

iteracoes = 1;
xs=[x0];
xk = x0;

sk = s0;

wk = [0;0;0];

func = [];


%% Algoritmo

while(iteracoes<maxi)
    
    for i=1:length(xk)
        
        Xk(i,i) = xk(i);
        
    end
    
    
    for i=1:length(sk)
        
        Sk(i,i) = sk(i);
        
    end
    

    %% Cálculo de mik, uk, thetak, vk, tk, pk
    
    if iteracoes ==1
        
        mik=1.5;
    
    else
        
        mik=mik/2; % Mik de Armijjo
    
    end
    
%     
%     
%      mik = transpose(sk)*xk/length(sk);
%      
%      
%      if iteracoes == 1
%          
%          mik=1;
%          
%      end
    
   
    % Cálculo de tk
    
    tk =  b - A*xk;
    
    % Cálculo de uk
    
    uk = gradientevpa - transpose(A)*wk-sk;
    
    % Cálculo de vk
    
    vk = mik*(transpose(e)) - Xk*Sk*transpose(e);
    
    % Cálculo de pk
    
    pk = inv(Xk)*vk;
    
    % Cálculo de thetak
    
    thetak = inv(hessianavpa+inv(Xk)*Sk);
    
    %% Teste de Otimalidade
    
    if mik < e1 & norm(tk)/(norm(b)+1) < e2 & norm(uk)/(norm(gradientevpa)+1) < e3
        
        display('Solução Ótima Encontrada')
        
        break
        
    end
          
  
    %% Cálculo da Direção de Translação 
    
    
    dwk = inv(A*thetak*transpose(A))*(tk + A*thetak*(uk-pk));
       
    dxk = thetak*(transpose(A)*dwk - uk + pk);
    
    dsk = inv(Xk)*(vk-Sk*dxk);
        
    %% Teste de Ilimitariedade e Otimalidade

    if sum(abs(tk)<e1) == length(tk) & sum(abs(dxk)<e2)==length(dxk) & sum(abs(dxk)<e2)==length(dxk) & sum(abs(dsk)<e2)==length(dsk) 
        
        display('Problema Ilimitado')
        break 
    
    end
    
    if uk == 0 & dsk>0 & transpose(b)*dwk>0

        display('Solução Ótima Encontrada')
        break 
    end
    
    %% Cálculo do Passo  
 
    alfa = 0.9995;
    
        
    betap1 = min([1;-alfa*xk(dxk<0)./dxk(dxk<0)]);
    
    betad = min([1;-alfa*sk(dsk<0)./dsk(dsk<0)]);
        
    betaq = -transpose(dxk)*(gradientevpa)/(transpose(dxk)*hessianavpa*dxk);
   
    betap=betap1;

        if betaq>0

            betap = min(betap1,betaq);

        end
        
   
    
    %% Atualização das Variáveis e Novos Pontos
    
    func = [func,vpa(subs(eq))];

    xk = xk + betap*dxk;
    wk = wk + betad*dwk;
    sk = sk + betad*dsk;
    xs=[xs,xk];
    x1 = xk(1);
    x2 = xk(2);
    x3 = xk(3);
    x4 = xk(4);
    x5 = xk(5);

    gradientevpa = vpa(subs(gradiente));

    hessianavpa = vpa(subs(hessiana));
    
    iteracoes = iteracoes + 1;
    
end



func = [func,vpa(subs(eq))];

figure

plot(func)
title('Minimização da Função de Rosenbrock (Exercício 2 da Prova)')
xlabel('Iterações')
ylabel('Valor da Função Objetivo')
grid on
grid minor


display(['Iteração Final: ' num2str(iteracoes-1)])

display('Norma do Resíduo Primal tk: ')
norm(tk)

display('Norma do Resíduo Dual uk: ')
norm(uk)

display('Norma do Resíduo de Complementaridade vk: ')
norm(vk)

display('Valor final do parâmetro de barreira: ')
mik

display('Valor da Função Objetivo: ')

func(iteracoes)
    
display('Solução Dual:')
display(sk)
display(wk)


display('Solução Primal:')
display(xk)