clc
clear all
close all

% Primal Afim Linear

% Completar desigualdades com s, não aumenta A, tem Sk. s > 0, x livre

A = [-1 1; 1 0; 0 1; -1 0; 0 -1];

b = [1;3;2;0;0];

c = [-2;-1];

xk = [2;1];

sk = b - A*xk;

e = ones(1,length(sk));

tol = 1e-4;

iteracoes=1;

max = 5;

Sk = zeros(length(sk),length(sk));

while(iteracoes<max)
    
    for i=1:length(sk)
        
       Sk(i,i) = sk(i);
        
    end

    
    % Cálculo de dxk
    
    dxk = -inv(transpose(A)*inv(Sk*Sk)*A)*c;
    
    % Cálculo de dyk
    
    dsk = -A*dxk;
    
    % Teste de Limitariedade
    
    if length(dsk(dsk>0))==length(dsk)
        
        display('Problema Limitado')
        
        break
        
    end
    
    if length(dsk(dsk<tol))==length(dsk)
        
        display('Solução Ótima Encontrada')
        
        break
        
    end
        
    % Cálculo de Wk
    
    wk = inv(Sk*Sk)*dsk;
    
    % Teste de Otimalidade
    
    [m,n] = size(A);
    
    if length(wk(wk>=0))==length(wk) & (((transpose(c)*xk-transpose(b)*wk)<tol) | e*Sk*wk <= tol | length(sk(sk<tol))==m-n)
        
        display('Solução Ótima Encontrada')
        
        break
        
    end
    
    % Cálculo do Passo
    
    betak = min(-0.9995*sk(dsk<0)./dsk(dsk<0));
    
    xk = xk + betak*dxk;
    
    sk = sk + betak*dsk;
    
    
    iteracoes = iteracoes+1   
    
        
end
