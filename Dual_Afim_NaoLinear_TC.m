
clc
clear all

x0 = [2.5;1.5];

syms x1 x2

eq = (x1-4)*(x1-4)*(x1-4)*(x1-4) + (x1-4*x2)*(x1-4*x2);

gradiente=gradient(eq,[x1,x2]);

hessiana = hessian(eq,[x1,x2]);

x1 = x0(1);
x2 = x0(2);

gradientevpa = vpa(subs(gradiente));

hessianavpa = vpa(subs(hessiana));

% Inicia Valor de Alfa:

alfa = 0.9995;

% Define a Matriz A:

A = [-1,1; -1,-1; 1,0; 0,1; -1,0; 0,-1];

% Define Vetor b:

b = [1;-1;3;2;0;0];

% Define valor de tolerancia episolon:

episolon = 0.1;

iteracoes = 1;

s0 = b - A*x0;

Sk = zeros(length(x0),length(x0));


% Inicializa vetor de e [1,1,1....,1]

e = ones(1,length(s0));

xk = x0;
sk = s0;

uk=1.5;

maxi = 10;



while(iteracoes<maxi)
    
    for i=1:length(sk)
        
        Sk(i,i) = sk(i);
        
    end
    
     % Cálculo da Direção de Translação 

     dxk = -inv(transpose(A)*inv(Sk*Sk)*A)*(gradientevpa+uk*transpose(A)*inv(Sk)*transpose(e));
     
     dsk = -A*dxk;

    episolon3 = 0.001;

    % Teste de Ilimitariedade dsk
    
    if length(dsk(dsk>0)) == length(dsk)

        display('Problema Ilimitado')
        condicao_de_parada=1;

    end
     
    episolon3 = 0.1;
    
    % Teste de Otimalidade dsk
    
    if length(dsk(dsk<episolon3)) == length(dsk)

        display('Solucao Otima Encontrada')
        break
              
    end
    
    % Cálculo do Vetor Estimativa Dual
     
     wk = -inv(Sk*Sk)*dsk; 

    
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
   
   bk2 = -transpose(dxk)*((gradientevpa+uk*transpose(A)*inv(Sk)*transpose(e)))/(transpose(dxk)*(hessianavpa-uk*transpose(A)*inv(Sk*Sk)*transpose(e))*dxk);
   % OU bk2 =
   % -transpose(dxk)*((gradientevpa+uk*transpose(A)*inv(Sk)*transpose(e)))/(transpose(dxk)*(hessianavpa)*dxk);
   % ?????????????????????????????????????
   
   if bk2<0
   
       bk = bk1;
   
   else
       
       bk = alfa*min(bk1,bk2);
   
   end
   
   
   sant=sk;
   want=wk;
   
   xk = xk+bk*dxk;
   sk = sk+bk*dsk;
   
    x1 = xk(1);
    
    x2 = xk(2);
    
    
    gradientevpa = vpa(subs(gradiente));

    hessianavpa = vpa(subs(hessiana));
    
    uk=uk/2;

      
   iteracoes=iteracoes+1;

end
