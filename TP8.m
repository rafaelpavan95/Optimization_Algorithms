
clear all
clc
close all


% Primal Afim PPQ

% Matriz Q*x = Gradiente
% Matriz Q = Hessiana
% Desconsidera a constante porque ela só translada a função
% Aumentar Q e A com as variáveis de folga X inseridas 
% Inicia Valor de Alfa:

alfa = 0.9995;

f = [];

% Define a Matriz A:

A = [1 1 0 ; 0 1 1];

% Define vetor c:

c = [1; 2; -3];

% Define Vetor b:

b = [5; 10];

% Define Matriz Q

Q = [4 0 0; 0 6 0; 0 0 10];

x0 = -inv(Q)*c;

if (A*x0) == transpose(b)
    
    if length(x0(x0>=0)) == length(x0)

        display("Pare. Solução Ótima Encontrada")

    end
end
 
iteracoes = 1;

x0 = [1.6;3.4;6.6];

Xk = zeros(length(x0),length(x0));

xk = x0;

sk = [1.5;1.5;1.5];

s0 = sk;

wk = [0;0];

e = ones(1,length(xk));

episolon1 = 9.99e-2;
episolon2 = 9.99e-2;
episolon3 = 9.99e-2;

maxi = 4;

miks =[];

Xk = zeros(length(xk),length(xk));
Sk = zeros(length(xk),length(xk));

while(iteracoes<maxi)
    
    for i=1:length(xk)
        
        Xk(i,i) = xk(i);
        
    end
    
    
    for i=1:length(sk)
        
        Sk(i,i) = sk(i);
        
    end
    
    
    % Cálculo de mik

    mik = transpose(sk)*xk/length(sk); % quantidade de variáveis de folga no problema dual
    
    
    if iteracoes == 1
        
        mik=1
        
    end
    miks=[miks,mik];
    
    % Cálculo de tk
    
    tk =  b - A*xk;
    
    % Cálculo de uk
    
    uk = Q*xk + c - transpose(A)*wk-sk;
    
    % Cálculo de vk
    
    vk = mik*(transpose(e)) - Xk*Sk*transpose(e);
    
    % Cálculo de pk
    
    pk = inv(Xk)*vk;
    
    % Cálculo de thetak
    
    thetak = inv(Sk+Xk*Q);
    
    
    if mik < episolon1 & norm(tk)/(norm(b)+1) < episolon2 & norm(uk)/(norm(c+Q*xk)+1) < episolon3
        
        display('Solução Ótima Encontrada')
        
        break
        
    end
       
        
        
  
    % Calcula Direção de Translação 
    
    
    % Com as equaçeõs da aula o algoritmo não funcionou. Com as equações do
    % Livro da Puthenpura funcionou! (7.89b e 7.89c)

    dwk = inv(A*thetak*Xk*transpose(A))*(tk + A*thetak*(Xk*uk-vk));
       
    dxk = thetak*(Xk*(transpose(A)*dwk - uk) + vk);
    
    dsk =  inv(Xk)*(vk-Sk*dxk);
        
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
    
 
    alfa = 0.995;
   
    
    lista_x = max(dxk./(-alfa*xk));
    lista_s = max(dsk./(-alfa*sk));
    

    betap = 1/max([lista_x,1]);
    
    betad = 1/max([lista_s,1]);
    
    
    % Atualização das Variáveis
    
    xa = xk;
    xk = xk + betap*dxk;
    wk = wk + betad*dwk;
    sk = sk + betad*dsk;
    
    xk
 
    iteracoes = iteracoes + 1;
    f = [f; transpose(xk)*Q*xk/2 + transpose(c)*xk];
    
    
end


plot(f)
grid on
grid minor
xlabel('Iterações')
ylabel('Função Objetivo')
title('Primal - Dual para PPQ')