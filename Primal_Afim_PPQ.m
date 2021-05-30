clear all
clc

% Primal Afim PPQ

% Matriz Q*x = Gradiente
% Matriz Q = Hessiana
% Desconsidera a constante porque ela só translada a função
% Aumentar Q e A com as variáveis de folga X inseridas 
% Inicia Valor de Alfa:

alfa = 0.9995;

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

max = 100;

while(iteracoes<max)
    
    for i=1:length(xk)
        
        Xk(i,i) = xk(i);
        
    end
    
     % Cálculo do vetor estimativa dual

    Hk = inv((Q + inv(Xk)*inv(Xk)));

    wk = inv(A*Hk*transpose(A)) * A*Hk*(Q*xk+c);

    % Cálculo do vetor custo relativo

    sk = (Q*xk + c) - transpose(A)*wk;
    
    episolon3 = 0.001;

    % Teste de Otimalidade

    if length(xk(xk>0)) == length(x0)

        if length(sk(sk>0)) == length(sk)

            if x0*sk < episolon3

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

                    if  norm(sk)/(norm(Q*xk+c)+1) < episolon2

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

                if  norm(sk)/(norm(Q*xk+c)+1) < episolon2

                    display('Factibilidade Dual Atingida')
                 
                end
            end

    end

     % Cálculo do Comprimento do Passo:

    alfak2 = -transpose(dxk)*(Q*xk+c)/((transpose(dxk)*Q*dxk));

    alfak1 = min(-alfa*xk(dxk<0)./dxk(dxk<0));

    alfak = alfa*min(alfak1,alfak2);

    if length(alfak1) == 0

        alfak = alfak2;

    end

    % Nova Solução:

    xk = xk + alfak*dxk;
      
    iteracoes=iteracoes+1;

end