clc
clear all
close all

% Rafael Pavan
% Método Primal Dual Barreira Logarítmica Previsor Corretor Mehota

% Exemplo 2


A = [-1 1 1 0; 1 1 0 1];
b = [1 ; 5];
c = [-1; -2; 0; 0];

xk = [1;1;3;3];
sk = [2.5;2.5;2.5;2.5];
wk = [0;-0.5];

e = ones(1,length(xk));

episolon1 = 1e-2;
episolon2 = 1e-2;
episolon3 = 1e-2;

iteracoes=1;

maxi = 100;

Xk = zeros(length(xk),length(xk));
Sk = zeros(length(xk),length(xk));

mik0 = transpose(xk)*sk/length(sk); 

f_obj =[];

while(iteracoes<maxi)
    
    for i=1:length(xk)
        
        Xk(i,i) = xk(i);
        
    end
    
    
    for i=1:length(sk)
        
        Sk(i,i) = sk(i);
        
    end
    
    
    % Cálculo de mik
    
    
    if iteracoes>1
        
        mik = transpose(xk)*sk/length(sk); % quantidade de variáveis de folga no problema dual % constante pequena conforme livro puthenpura
    
    end
    
    mik=2;
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
   
    % Passo Previsor    
        
  
    % Calcula Direção de Translação 
    
    pk = zeros(length(pk),1);
    
    dwk_p = inv(A*thetak*transpose(A))*(tk + A*thetak*(uk-pk));
        
    dsk_p = uk - transpose(A)*dwk_p;
    
    dxk_p = thetak*(pk-dsk_p);
    
 
    % Cálculo do Passo
    
    alfa = 0.9995;
   
    lista_x = max(-dxk_p./alfa.*xk);
    lista_s = max(-dsk_p./alfa.*sk);
    
    lista_x = [lista_x,1];
    lista_s = [lista_s,1];
    
   betap = (1/max(lista_x));
   betad = (1/max(lista_s));
    
    
    % Atualização das Variáveis
    
   
    xk = xk + betap*dxk_p;
    wk = wk + betad*dwk_p;
    sk = sk + betad*dsk_p; 
    
        
    mik_p = transpose(xk)*sk/length(sk);
    
    % Passo Corretor
    
    Dsk_p = zeros(length(dsk_p),length(dsk_p));
    
    Dxk_p = zeros(length(dxk_p),length(dxk_p));
    
    for i=1:length(dxk_p)
        
        Dxk_p(i,i) = dxk_p(i);
        
    end
    
    
    for i=1:length(dsk_p)
        
        Dsk_p(i,i) = dsk_p(i);
        
    end
    
   
    vk_c = mik_p*transpose(e) - Xk*Sk*transpose(e) - Dxk_p*Dsk_p*transpose(e);
   
    % Calcula Direção de Translação 
    
    pk=inv(Xk)*vk_c;
    
    dwk = inv(A*thetak*transpose(A))*(tk + A*thetak*(uk-pk));
        
    dsk = uk - transpose(A)*dwk;
    
    dxk = thetak*(pk-dsk);
    
    
    % Cálculo do Passo
    
    alfa = 0.9995;
   
    lista_x = max(-dxk./alfa.*xk);
    lista_s = max(-dsk./alfa.*sk);
    
    lista_x = [lista_x,1];
    lista_s = [lista_s,1];
    
   betap = (1/max(lista_x));
   betad = (1/max(lista_s));
         
    
    % Teste de Otimalidade
    
    if sum(tk==0) == length(tk) & sum(dxk>0)==length(dxk) & transpose(c)*dxk<0
        
        display('Problema Ilimitado')
        break 
    
    end
    
    if uk == 0 & dsk>0 & transpose(b)*dwk>0

        display('Solução Ótima Encontrada')
        break 
    end
    
   if mik < episolon1 & norm(tk)/(norm(b)+1) < episolon2 & norm(uk)/(norm(c)+1) < episolon3
        
        display('Solução Ótima Encontrada')
        
        break
        
    end
    
    % Atualização das Variáveis
    
   
    xk = xk + betap*dxk;
    wk = wk + betad*dwk;
    sk = sk + betad*dsk; 
    
    
    iteracoes = iteracoes + 1;
    
    f_obj = [f_obj,transpose(c)*xk];
    
    xk
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
title('Exemplo 2 - Curva')