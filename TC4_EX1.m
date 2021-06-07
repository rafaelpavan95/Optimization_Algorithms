clc
clear all
close all

% Rafael Pavan
% Método Primal Dual Barreira Logarítmica

% Exemplo 1


A = [1 -1 1 0; 0 1 0 1];
b = [15 ; 15];
c = [-2; 1; 0; 0];

xk = [2.5;2.5;2.5;2.5];

sk = [2.5;2.5;2.5;2.5];

wk = [1;1];

e = ones(1,length(xk));

episolon1 = 1e-2;
episolon2 = 1e-2;
episolon3 = 1e-2;

iteracoes=1;

maxi = 1000;

Xk = zeros(length(xk),length(xk));
Sk = zeros(length(xk),length(xk));

f_obj = [];

while(iteracoes<maxi)
    
    for i=1:length(xk)
        
        Xk(i,i) = xk(i);
        
    end
    
    
    for i=1:length(sk)
        
        Sk(i,i) = sk(i);
        
    end
    
    
    % Cálculo de mik
    

    mik = 0.9*transpose(xk)*sk/length(sk); % quantidade de variáveis de folga no problema dual
    
    % Cálculo de tk
    
    tk =  b - A*xk;
    
    % Cálculo de uk
    
    uk = c - transpose(A)*wk-sk;
    
    
    % Cálculo de vk
    
    vk = mik*(transpose(e)) - Xk*Sk*transpose(e);
    
    % Cálculo de pk
    
    pk = inv(Xk)*vk;
    
    
    % Cálculo de thetak
    
    thetak = inv(Sk)*Xk;
    
    
    if mik < episolon1 & norm(tk)/(norm(b)+1) < episolon2 & norm(uk)/(norm(c)+1) < episolon3
        
        display('Solução Ótima Encontrada')
        
        break
        
    end
       
        
        
  
    % Calcula Direção de Translação 
    
    
    % Com as equaçeõs da aula o algoritmo não funcionou. Com as equações do
    % Livro da Puthenpura funcionou! (7.89b e 7.89c)
    
    dwk = inv(A*thetak*transpose(A))*(tk + A*thetak*(uk-pk));
        
    dsk = uk - transpose(A)*dwk;
    
    dxk = thetak*(pk-dsk);
    
    % Teste de Otimalidade
    
    if sum(tk==0) == length(tk) & sum(dxk>0)==length(dxk) & transpose(c)*dxk<0
        
        display('Problema Ilimitado')
        break 
    
    end
    
    if uk == 0 & dsk>0 & transpose(b)*dwk>0

        display('Solução Ótima Encontrada')
        break 
    end
    
    % Cálculo do Passo
    
 
    alfa = 0.9995;
   
    lista_x = max(-dxk./alfa.*xk);
    lista_s = max(-dsk./alfa.*sk);
    
    lista_x = [lista_x,1];
    lista_s = [lista_s,1];
    
    betap = (1/max(lista_x));
    betad = (1/max(lista_s));

    
    % Atualização das Variáveis
    
   
    xk = xk + betap*dxk;
    wk = wk + betad*dwk;
    sk = sk + betad*dsk; 
    
    f_obj = [f_obj,transpose(c)*xk];


    iteracoes = iteracoes + 1;
end



display('Solução Ótima: ')

xk

display('Função Objetivo: ')

transpose(c)*xk

plot(f_obj);
grid on
grid minor
xlabel('Iterações')
ylabel('Função Objetivo')
title('Exemplo 1 - Curva')