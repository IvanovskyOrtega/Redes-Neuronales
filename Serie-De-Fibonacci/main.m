clc % Limpiamos la ventana de comandos
clear % Eliminamos los items del workspace
resultados = fopen('numerosDeFibonacci.txt','w');
limite = inputdlg('Ingresa el limite para la serie de Fibonacci:'); % Archivo para guardar los terminos
n = str2double(limite{:}); % Convertimos a cadena a numero
A = zeros(1,n+1);   % Arreglo para guardar los valores de la serie
fprintf(resultados,'Primeros %d términos de la Serie de Fibonacci:\n',n);
for i=1:n+1
    A(:,i) = fibonacci(i-1); % Guardamos los valores calculados en un array
    fprintf(resultados,'%d\n',A(:,i));
end
fclose(resultados);
figure % Se plotean los valores
rango = 0:1:n;
stem(rango,A);
xlabel('n')
ylabel('Fibonacci(n)')
titulo = strcat('Serie de Fibonacci de 0 a',{' '},num2str(n));
title(titulo)
