clc
clear all

x0 = [9;4;10;2;11];

syms x1 x2 x3 x4 x5

eq = 100*(x1-x2*x2)*(x1-x2*x2) + (1-x1)*(1-x1) + 0*x3 + 0*x4 + 0*x5;

gradiente=gradient(eq,[x1,x2,x3,x4,x5]);

hessiana = hessian(eq,[x1,x2,x3,x4,x5]);

x1 = x0(1);
x2 = x0(2);
x3 = x0(3);
x4 = x0(4);
x5 = x0(5);

gradientevpa = vpa(subs(gradiente));

hessianavpa = vpa(subs(hessiana));

% Inicia Valor de Alfa:

alfa = 0.9995;

% Define a Matriz A:

A = [1, -1, 1, 0, 0; 0, -1, 0, 1, 0; 0, 1, 0, 0, 1];

% Define Vetor b:

b = [15; -2; 15];

% Define valor de tolerancia episolon:

episolon = 0.001;

% Define matriz Xk diagonal:

Xk = zeros(length(x0),length(x0));

% Inicializa vetor de e [1,1,1....,1]

e = ones(1,length(x0));

maxi = 100;

iteracoes = 1;

xk = x0;

uk=1.5;

while (iteracoes<maxi)

    for i=1:length(xk)
        
        Xk(i,i) = xk(i);
        
    end
    
     % Cálculo do vetor estimativa dual

    Hk = inv((hessianavpa + uk*inv(Xk*Xk) + inv(Xk)*inv(Xk)));

    wk = inv(A*Hk*transpose(A)) * A*Hk*(gradientevpa - uk*inv(Xk)*transpose(e));

    % Cálculo do vetor custo relativo

    sk = (gradientevpa - uk*inv(Xk)*transpose(e)) - transpose(A)*wk;
    
    episolon3 = 0.001;

    % Teste de Otimalidade

    if length(xk(xk>0)) == length(x0)

        if length(sk(sk>0)) == length(sk)

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

                    if  norm(sk)/(norm((gradientevpa - uk*inv(Xk)*transpose(e)))+1) < episolon2

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

                if  norm(sk)/(norm((gradientevpa - uk*inv(Xk)*transpose(e)))+1) < episolon2

                    display('Factibilidade Dual Atingida')
                 
                end
            end

    end

     % Cálculo do Comprimento do Passo:

    alfak2 = -transpose(dxk)*((gradientevpa - uk*inv(Xk)*transpose(e)))/((transpose(dxk)*(hessianavpa+uk*inv(Xk*Xk))*dxk));

    alfak1 = min(-alfa*xk(dxk<0)./dxk(dxk<0));

    alfak = alfa*min(alfak1,alfak2);

    if length(alfak1) == 0

        alfak = alfak2;

    end

    % Nova Solução:

    xk = xk + alfak*dxk;
      
    x1 = xk(1);
    
    x2 = xk(2);
    
    x3 = xk(3);
    
    x4 = xk(4);
    
    x5 = xk(5);

    iteracoes=iteracoes+1;

    gradientevpa = vpa(subs(gradiente));

    hessianavpa = vpa(subs(hessiana));
    
    uk=uk/2;
    
    iteracoes=iteracoes+1;

end
