clear all
clc
close all

% Primal Afim PPQ

% Matriz Q*x = Gradiente
% Matriz Q = Hessiana
% Desconsidera a constante porque ela só translada a função
% Aumentar Q e A com as variáveis de folga X inseridas 
% Inicia Valor de Alfa:

alfa = 0.995;

% Define a Matriz A:

A = [1 1 1 0 0 0 0 0 0; -1 0 0 1 0 0 0 0 0; 1 0 0 0 1 0 0 0 0; 0 -1 0 0 0 1 0 0 0; 0 1 0 0 0 0 1 0 0; 0 0 -1 0 0 0 0 1 0; 0 0 1 0 0 0 0 0 1];

% Define vetor c:

c = [7.92;7.97;7.85;0;0;0;0;0;0];

% Define Vetor b:

b = [850;-100;600;-50;200;-100;400];

% Define Matriz Q

Q = [2*0.001562 0 0 0 0 0 0 0 0; 0 2*0.00482 0 0 0 0 0 0 0; 0 0 2*0.00194 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0];

x0 = -inv(Q)*c;


uk = 1.5;


if (A*x0) == transpose(b)
    
    if length(x0(x0>=0)) == length(x0)

        display("Pare. Solução Ótima Encontrada")

    end
end
 
iteracoes = 1;

x0 = [500; 100; 250; 350; 150; 50; 100; 200; 100];

Xk = zeros(length(x0),length(x0));

xk = x0;

xx = [xk];

max = 100;

e = ones(1,length(xk));

curva = [];

while(iteracoes<max)
    
    for i=1:length(xk)
        
        Xk(i,i) = xk(i);
        
    end
    
     % Cálculo do vetor estimativa dual
     
    Hk = inv(Q+(uk+1)*inv(Xk)^2);    
    wk = inv(A*Hk*transpose(A))*A*Hk*(Q*xk+c-uk*inv(Xk)*transpose(e));

    % Cálculo do vetor custo relativo

    sk = (Q*xk + c - uk*inv(Xk)*transpose(e)) - transpose(A)*wk;
    
    episolon3 = 0.001;

    % Teste de Otimalidade

    if length(xk(xk>0)) == length(x0)

        if length(sk(sk>-episolon3)) == length(sk)

            if transpose(xk)*sk < episolon3

                display('Solução Ótima Encontrada')
                break
                               
                
                % Teste de Factibilidade

                episolon1 = 0.001;

                episolon2 = 0.1;


                if length(xk(xk>=0)) == length(xk)

                            if norm(A*xk - b)/(norm(b)+1) < episolon1

                            display('Factibilidade Primal Atingida')
                            

                    end
                end

                if length(sk(sk>episolon3)) == length(sk)

                    if  norm(sk)/(norm(Q*xk+c-uk*inv(Xk)*transpose(e))+1) < episolon2

                        display('Factibilidade Dual Atingida')
                        

                    end
                end


            end
            
        end
    end


   % Calculo da Direcao de Translacao

    dxk = -Hk*sk;

    % Teste de Ilimitariedade
    
    if length(dxk(dxk>0)) == length(dxk)

        display('Problema Ilimitado')
        condicao_de_parada=1;

    end
     
    % Teste de Otimalidade dxk
    
    if length(dxk(dxk<episolon3)) == length(dxk)

        display('Solucao Otima Encontrada')
        break
              

            % Teste de Factibilidade

            episolon1 = 0.001;

            episolon2 = 0.1;


            if length(xk(xk>=0)) == length(xk)

                        if norm(A*xk - b)/(norm(b)+1) < episolon1

                    display('Factibilidade Primal Atingida')

                end
            end

            if length(sk(sk>=0)) == length(sk)

                if  norm(sk)/(norm(Q*xk+c-uk*inv(Xk)*transpose(e))+1) < episolon2

                    display('Factibilidade Dual Atingida')
                 
                end
            end

    end

     % Cálculo do Comprimento do Passo:

    alfak2 = (-transpose(dxk)*(Q*xk+c)/(transpose(dxk)*(Q)*dxk));

    alfak1 = min(-alfa*xk(dxk<0)./dxk(dxk<0));

    alfak = alfa*min(alfak1,alfak2);

    if length(alfak1) == 0

        alfak = alfak2;

    end

    % Nova Solução:

    xk = xk + alfak*dxk;
    
    uk = uk/2;
      
    iteracoes=iteracoes+1;
    
    curva=[curva,(transpose(xk)*Q*xk/2)+(transpose(c)*xk)+561+78+310];

end


display('Solução [MW]')
disp(['Gerador 1: ' num2str(xk(1))])
disp(['Gerador 2: ' num2str(xk(2))])
disp(['Gerador 3: ' num2str(xk(3))])

display(['Função Objetivo [$]: ',num2str((transpose(xk)*Q*xk/2)+(transpose(c)*xk)+561+78+310)])

plot(curva)
grid on
grid minor
title('Custo de Geração')
xlabel('Iterações')
ylabel('Custo ($)')
