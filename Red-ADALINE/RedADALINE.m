clc
clear
iteraciones = input('Ingresa el numero de iteraciones:');
alpha = input('Ingresa el valor de alfa: ');
eit = input('Ingresa el valor de la señal de error: ');
W = rand(1,3); %Asignamos valores aleatorios a la matriz de pesos
%W = [0.84 0.39 0.78]; 
d = [           % Conjunto de entrenamiento
    0 0 0;
    0 0 1; 
    0 1 0; 
    0 1 1; 
    1 0 0; 
    1 0 1; 
    1 1 0; 
    1 1 1
    ];
W_1 = zeros(1,9);   % Vector para guardar los valores de cambio de W[1]
W_2 = zeros(1,9);   % Vector para guardar los valores de cambio de W[2]
W_3 = zeros(1,9);   % Vector para guardar los valores de cambio de W[3]
ErrorGlobal = zeros(1,9);   % Vector para guardar los valores de cambio de la señal de error
t = 0.0;    % Iniciamos t en 0
Eit = 0.0;  % Error global
figure
rango = 0:8;    % Rango para el plot
for i = 1:iteraciones   % Loop sobre el numero maximo de iteraciones
    fprintf('Inicia la iteracion %d\n',i);  
    for j = 1:8     % Loop sobre los elementos del cto. de entrenamiento
        W_1(j) = W(1);  % Asignamos los valores de cambio de las entradas 
        W_2(j) = W(2);  % de W y la señal de error para el ploteo
        W_3(j) = W(3);
        ErrorGlobal(j) = Eit;
        %fprintf('\tSe propaga d%d hacia adelante\n',j);
        a = purelin(W*d(j,:)');
        ej = t - a;
        t = t+1;
        fprintf('\ta = %f\n',a);
        fprintf('\tej = %f\n',ej);
        if ej ~= 0.0 % Si hubo cambios en la señal del error del termino j
            W = W+(2*alpha*ej*d(j,:)); % Se ajustan los valores de la matriz de pesos
        end
        fprintf('\tW =[%f,%f,%f]\n',W(1),W(2),W(3));
        Eit = Eit+(1/8)*ej; % Actualizamos el error global
    end
    W_1(j+1) = W(1);  % Asignamos los valores de cambio de las entradas  
    W_2(j+1) = W(2);  % de W y la señal de error para el ploteo
    W_3(j+1) = W(3);
    ErrorGlobal(j+1) = Eit;
    fprintf('\tEit = %f\n',Eit);
    subplot(2,ceil(iteraciones/2),i);   % Creamos un subplot por iteracion
    plot(rango,W_1);
    text(8,W_1(8),'\leftarrow W[1]');
    titulo = strcat('Iteracion',{' '},num2str(i));
    title(titulo);
    hold on
    plot(rango,W_2);
    text(8,W_2(8),'\leftarrow W[2]');
    plot(rango,W_3);
    text(8,W_3(8),'\leftarrow W[3]');
    plot(rango,ErrorGlobal);
    text(8,ErrorGlobal(8),'\leftarrow S. error');
    axis([0 8 0 inf])
    hold off
    xlabel({'Propagación del elemento i-ésimo','Del conjunto de Entrenamiento'},'FontSize',7,'FontWeight','bold')
    ylabel('Valor','FontWeight','bold')
    if Eit == 0 % Si el error global es cero (caso ideal)
        fprintf('Se ha tenido un aprendizaje exitoso en la iteracion %d\n',i);
        break   % Termina el aprendizaje de la red
    end
    if Eit < eit % Si el error global es menor al valor establecido de la señal del error eit
        fprintf('Se ha tenido un aprendizaje exitoso en la iteracion %d\n',i);
        break    % Termina el aprendizaje de la red 2/3 de probabilidad de asegurar un aprendizaje adecuado
    end
    Eit = 0;    % Se reinician el valor del error global y t
    t = 0;
end
clear