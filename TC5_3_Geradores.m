clear all
clc
close all

% Método Primal Dual para PPQ - Despacho Econômico 3 Geradores
% Rafael Pavan

% Define a Matriz A:

A = [1 1 1 0 0 0 0 0 0; -1 0 0 1 0 0 0 0 0; 1 0 0 0 1 0 0 0 0; 0 -1 0 0 0 1 0 0 0; 0 1 0 0 0 0 1 0 0; 0 0 -1 0 0 0 0 1 0; 0 0 1 0 0 0 0 0 1];

% Define vetor c:

c = [7.92;7.97;7.85;0;0;0;0;0;0];

% Define Vetor b:

b = [850;-100;600;-50;200;-100;400];

% Define Matriz Q

Q = [2*0.001562 0 0 0 0 0 0 0 0; 0 2*0.00482 0 0 0 0 0 0 0; 0 0 2*0.00194 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0];


% Realiza Teste

x0 = -inv(Q)*c;


if (A*x0) == transpose(b)
    
    if length(x0(x0>=0)) == length(x0)

        display("Pare. Solução Ótima Encontrada")

    end
end
 
iteracoes = 1;

% Define Solução Inicial

x0 = [500; 100; 250; 350; 150; 50; 100; 200; 100];

Xk = zeros(length(x0),length(x0));

xk = x0;

e = ones(1,length(xk));

custo = [];

sk = 1./x0;

% Define Tolerâncias

episolon1 = 9.99e-2;
episolon2 = 9.99e-2;
episolon3 = 9.99e-2;

maxi = 10;

miks =[];

wk = zeros(7,1);

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
        
        mik=1;
        
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
    
    % Teste de Otimalidade
    
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
    
    
    % Passos conforme livro da Fang, Puthenpura.
    
    betap = 1/max([lista_x,1]);
    
    betad = 1/max([lista_s,1]);
    
    
    % Atualização das Variáveis
    
    xa = xk;
    xk = xk + betap*dxk;
    wk = wk + betad*dwk;
    sk = sk + betad*dsk;
 
    iteracoes = iteracoes + 1;
    
    custo = [custo,(transpose(xk)*Q*xk/2)+(transpose(c)*xk)+561+78+310];
    
end

display('Solução [MW]')
disp(['Gerador 1: ' num2str(xk(1))])
disp(['Gerador 2: ' num2str(xk(2))])
disp(['Gerador 3: ' num2str(xk(3))])

display(['Função Objetivo [$]: ',num2str((transpose(xk)*Q*xk/2)+(transpose(c)*xk)+561+78+310)])

plot(custo)
grid on
grid minor
title('Custo de Geração')
xlabel('Iterações')
ylabel('Custo ($)')