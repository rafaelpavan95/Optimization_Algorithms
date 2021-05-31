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

uk = 1.5;

iteracoes = 1;

x0 = [600; 300; 300; 150;150;150;150;150;150;105;105;105;105;600;80;300;60;300;60;150-60;30;150-60;30;150-60;30;150-60;30;150-60;30;150-60;30;105-40;15;105-40;15;105-55;15;105-55;15];
    
Xk = zeros(length(x0),length(x0));

xk = x0;

xx = [xk];

max = 500;

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
    
    alfak = alfa*min([alfak1,alfak2,1]);

    if length(alfak1) == 0

        alfak = alfak2;

    end

    % Nova Solução:

    xk = xk + alfak*dxk;
    
    uk = uk/2;
      
    iteracoes=iteracoes+1;
    
    curva=[curva,(transpose(xk)*Q*xk/2)+(transpose(c)*xk)+550+309+307+240+240+240+240+240+240+126+126+126+126];

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

plot(curva)
grid on
grid minor
title('Custo de Geração')
xlabel('Iterações')
ylabel('Custo ($)')