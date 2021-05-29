clc
clear all
close all

% Primal Afim Linear

% Completar desigualdades com X, não esquecer de inserir zeros no vetor C.

A = [1 -1 1 0; 0 1 0 1];
b = [15 ; 15];
c = [-2; 1; 0; 0];
xk = [10;6;11;9];

e = ones(1,length(xk));

tol = 1e-4;

iteracoes=1;
max = 5;

Xk = zeros(length(xk),length(xk));

while(iteracoes<max)
    
    for i=1:length(xk)
        
        Xk(i,i) = xk(i);
        
    end
    
    
    
    % Cálculo de Wk
    
    wk = inv(A*Xk*Xk*transpose(A))*A*Xk*Xk*c;
    
    % Cálculo de rk
    
    rk = c-transpose(A)*wk;
    
    if length(rk(rk>=0)) == length(rk) && e*Xk*rk < tol
        
        display('Solução Ótima Encontrada')
        
        break
        
    end
    
   
    % Calcula Direção de Translação
    
    dyk = -Xk*rk;
    
        
    if length(dyk(dyk==0)) == length(dyk)
        
        display('Solução Ótima Encontrada')
        
        break
        
    end
    
    % Cálculo do Passo
    
    alfak = min(-0.9995./dyk(dyk<0));
    
    yk = transpose(e) + alfak*dyk;
    
    xk = Xk*yk;
    
    iteracoes = iteracoes+1   
    
    % Teste Adicional
    
    [n,m] = size(A);
    
    if length(xk(xk<tol)) == m-n
        
        display('Solução Ótima Encontrada')
        
        break
        
    end
    
    

end