clear all
clc

% Dual Afim PPQ

% Matriz Q*x = Gradiente
% Matriz Q = Hessiana
% Desconsidera a constante porque ela só translada a função
% Não Aumentar Q e A com as variáveis de folga X inseridas 
% Inicia Valor de Alfa:

alfa = 0.9995;

% Define a Matriz A:

A = [-1 1; -1 -1;1 0;0 1; -1 0; 0 -1];

% Define vetor c:

c = [4; -6];

% Define Vetor b:

b = [1;-1;3;2;0;0];

% Define Matriz Q

Q = [4 0; 0 6];

x0 = -inv(Q)*c;

if (A*x0) == transpose(b)
    
    if length(x0(x0>=0)) == length(x0)

        display("Pare. Solução Ótima Encontrada")

    end
end
 
iteracoes = 1;

x0 = [1.5;1];
uk = 1.5;
s0 = b - A*x0;

Sk = zeros(length(x0),length(x0));

xk = x0;
sk = s0;

e = ones(1,length(sk));

maxi = 10;

while(iteracoes<maxi)
    
    for i=1:length(sk)
        
        Sk(i,i) = sk(i);
        
    end
    
     % Cálculo da Direção de Translação 

     dxk = -inv(transpose(A)*inv(Sk*Sk)*A)*(Q*xk+c+uk*transpose(A)*inv(Sk)*transpose(e));
     
     dsk = -A*dxk;

    episolon3 = 0.001;

    % Teste de Ilimitariedade dsk
    
    if length(dsk(dsk>0)) == length(dsk)

        display('Problema Ilimitado')
        condicao_de_parada=1;

    end
     
    episolon3 = 0.01;
    
    % Teste de Otimalidade dsk
    
    if length(dsk(dsk<episolon3)) == length(dsk)

        display('Solucao Otima Encontrada')
        break
              
    end
    
    % Cálculo do Vetor Estimativa Dual
     
     wk = inv(Sk*Sk)*dsk; 

    
    % Teste de Otimalidade

    if length(wk(wk>=0)) == length(wk) & e*Sk*dsk < episolon3

                display('Solução Ótima Encontrada')
                break
                               
    end
    
    % Teste de Otimalidade Adicional
    
    if iteracoes ~= 1

        if max(abs(sk-sant))/max([1],max(abs(sk))) <= episolon3 & max(abs(wk-want))/max([1],max(abs(wk))) <= episolon3

    
                display('Solução Ótima Encontrada')
                break
                            
        end
       
        
        
    end
          
   % Cálculo do Passo Bk
   
   bk1 = min(-alfa*sk(dsk<0)./dsk(dsk<0));
   
   bk2 = -transpose(dxk)*(Q*xk+c)/(transpose(dxk)*Q*dxk);
   
   if bk2<0
   
       bk = bk1;
   
   else
       
       bk = alfa*min(bk1,bk2);   
   
   end
   
   
   sant=sk;
   want=wk;
   
   xk = xk+bk*dxk;
   sk = sk+bk*dsk;
   
   uk = uk/2;
      
   iteracoes=iteracoes+1;

end

