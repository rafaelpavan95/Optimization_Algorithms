clear all
clc
close all


% Método Primal Dual para PPQ - Despacho Econômico 13 Geradores
% Rafael Pavan

% Define a Matriz A:

A_sh = [1 1 1 1 1 1 1 1 1 1 1 1 1;
        -1 0 0 0 0 0 0 0 0 0 0 0 0;
         1 0 0 0 0 0 0 0 0 0 0 0 0;
         0 -1 0 0 0 0 0 0 0 0 0 0 0;
         0 1 0 0 0 0 0 0 0 0 0 0 0;
         0 0 -1 0 0 0 0 0 0 0 0 0 0;
         0 0 1 0 0 0 0 0 0 0 0 0 0;
         0 0 0 -1 0 0 0 0 0 0 0 0 0;
         0 0 0 1 0 0 0 0 0 0 0 0 0;
         0 0 0 0 -1 0 0 0 0 0 0 0 0;
         0 0 0 0 1 0 0 0 0 0 0 0 0;
         0 0 0 0 0 -1 0 0 0 0 0 0 0;
         0 0 0 0 0 1 0 0 0 0 0 0 0;
         0 0 0 0 0 0 -1 0 0 0 0 0 0;
         0 0 0 0 0 0 1 0 0 0 0 0 0;
         0 0 0 0 0 0 0 -1 0 0 0 0 0;
         0 0 0 0 0 0 0 1 0 0 0 0 0;
         0 0 0 0 0 0 0 0 -1 0 0 0 0;
         0 0 0 0 0 0 0 0 1 0 0 0 0;
         0 0 0 0 0 0 0 0 0 -1 0 0 0;
         0 0 0 0 0 0 0 0 0 1 0 0 0;
         0 0 0 0 0 0 0 0 0 0 -1 0 0;
         0 0 0 0 0 0 0 0 0 0 1 0 0;
         0 0 0 0 0 0 0 0 0 0 0 -1 0;
         0 0 0 0 0 0 0 0 0 0 0 1 0;
         0 0 0 0 0 0 0 0 0 0 0 0 -1;
         0 0 0 0 0 0 0 0 0 0 0 0 1];
     
A_comp = zeros(27,26);

for i=2:27
    
    A_comp(i,i-1)=1;
end

A = [A_sh,A_comp];

% Define vetor c:

csh = [8.1;8.1;8.1;7.74;7.74;7.74;7.74;7.74;7.74;8.6;8.6;8.6;8.6];
ccomp = zeros(26,1);

c = [csh;ccomp];

% Define Vetor b:

b = [2520;0;680;0;360;0;360;-60;180;-60;180;-60;180;-60;180;-60;180;-60;180;-40;120;-40;120;-55;120;-55;120];


% Define Matriz Q

Qcomp=2*[0.00028 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
       0 0.00056 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
       0 0 0.00056 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
       0 0 0 0.00324 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
       0 0 0 0 0.00324 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
       0 0 0 0 0 0.00324 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0.00324 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0.00324 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0 0.00324 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0 0 0.00284 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0 0 0 0.00284 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0 0 0 0 0.00284 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0 0 0 0 0 0.00284 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
   
Qzero = zeros(26,39);

Q = [Qcomp;Qzero];

iteracoes = 1;

% Define Solução Inicial

x0 = [600; 300; 300; 150;150;150;150;150;150;105;105;105;105;600;80;300;60;300;60;150-60;30;150-60;30;150-60;30;150-60;30;150-60;30;150-60;30;105-40;15;105-40;15;105-55;15;105-55;15];
    
Xk = zeros(length(x0),length(x0));

xk = x0;

e = ones(1,length(xk));

sk = 1./xk;


wk = zeros(27,1);

% Define Tolerância

episolon1 = 9.99e-2;
episolon2 = 9.99e-2;
episolon3 = 9.99e-2;

% Máxima Iteração

maxi = 10;

miks = [];

custo = [];

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
    % Livro da Puthenpura funcionou! (7.89b e 7.89c) pag 241

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
  
    custo=[custo,(transpose(xk)*Q*xk/2)+(transpose(c)*xk)+550+309+307+240+240+240+240+240+240+126+126+126+126];


end

display('Solução [MW]')
disp(['Gerador 1: ' num2str(xk(1))])
disp(['Gerador 2: ' num2str(xk(2))])
disp(['Gerador 3: ' num2str(xk(3))])
disp(['Gerador 4: ' num2str(xk(4))])
disp(['Gerador 5: ' num2str(xk(5))])
disp(['Gerador 6: ' num2str(xk(6))])
disp(['Gerador 7: ' num2str(xk(7))])
disp(['Gerador 8: ' num2str(xk(8))])
disp(['Gerador 9: ' num2str(xk(9))])
disp(['Gerador 10: ' num2str(xk(10))])
disp(['Gerador 11: ' num2str(xk(11))])
disp(['Gerador 12: ' num2str(xk(12))])
disp(['Gerador 13: ' num2str(xk(13))])

display(['Função Objetivo [$]: ',num2str((transpose(xk)*Q*xk/2)+(transpose(c)*xk)+550+309+307+240+240+240+240+240+240+126+126+126+126)])

plot(custo)
grid on
grid minor
title('Custo de Geração')
xlabel('Iterações')
ylabel('Custo ($)')