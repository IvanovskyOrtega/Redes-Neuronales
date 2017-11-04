clc 
clear
numProt = input('Ingresa el número de vectores prototipo: ');
r = input('Ingresa el tamaño (r) de los vectores prototipo:  ');
nombreArchivo = input('Ingresa el nombre del archivo de los vectores prototipo (con extension .txt):  ','s');
archivo = fopen(nombreArchivo,'r');
while (archivo == -1)
    disp('El archivo solicitado no existe o no se abrió correctamente.');
    nombreArchivo = input('Ingresa el nombre de un archivo válido:  ','s');
    archivo = fopen(nombreArchivo,'r');
end
disp('Se abrió correctamente el archivo.');
fclose(archivo);
W = dlmread(nombreArchivo);
dimensiones = size(W);
S = numProt;
if(dimensiones(1) ~= numProt || dimensiones(2) ~= r)
    disp('Las dimensiones ingresadas no coinciden.');
    disp('Se utilizaran las siguientes: ');
    cadena = sprintf('\tNumero de Vectores prototipo: %d',dimensiones(1));
    disp(cadena);
    cadena2 = sprintf('\tR: %d',dimensiones(2));
    disp(cadena2);
    numProt = dimensiones(1);
    S = dimensiones(1);
    r = dimensiones(2);
end
disp('W =');
disp(W);
archivoBias = input('Ingresa el nombre del archivo que contiene el vector bias (con extension .txt):  ','s');
archivo = fopen(archivoBias,'r');
while (archivo == -1)
    disp('El archivo solicitado no existe o no se abrió correctamente.');
    archivoBias = input('Ingresa el nombre de un archivo válido:  ','s');
    archivo = fopen(archivoBias,'r');
end
disp('Se abrió correctamente el archivo.');
fclose(archivo);
b = dlmread(archivoBias);
disp('b =');
disp(b);
epsilon = 0;
while epsilon == 0
    epsilon = mod(rand(),(1/S-1));
end
epsilonMatrix = ones(S,S)*epsilon; 
W_2 = tril(epsilonMatrix,-1)+triu(epsilonMatrix,1)+diag(ones(1,S));
disp('Se usará la siguiente matriz en la capa recurrente:');
disp(W_2);
%op = 's';
archivosGrafica = zeros(1,S);
for k=1:S
    nombre = sprintf('valoresDeGraficacion%d.txt',k);
    archivosGrafica(k) = fopen(nombre,'w');
end 
%while strcmp(op,'s');
    disp('Ingrese los valores del vector de entrada a propagar:');
    P = zeros(r,1);
    for i=1:r
        cad = sprintf('P(%d):',i);
        P(i) = input(cad);
    end
    t = 0;
    a1 = purelin(W*P+b);
    a2_1 = a1;
    for k=1:S
        fprintf(archivosGrafica(k),'%f\r\n',a2_1(k));
    end
    a2_2 = ones(S,1)*(-99999);
    while true
        a2_2 = poslin(W_2*a2_1);
        if isequal(a2_1,a2_2)
            break;
        end
        a2_1 = a2_2;
        for k=1:S
            fprintf(archivosGrafica(k),'%f\r\n',a2_1(k));
        end
        t = t+1;
    end
    a2_2 = round(a2_2);
    disp('a2 = ');
    disp(a2_2);
    for k=1:S
        fprintf(archivosGrafica(k),'%f\r\n',a2_2(k));
    end
    msg = sprintf('La red convergió exitosamente en %d iteraciones.',t);
    disp(msg);
    for k=1:S
        if a2_2(k) ~= 0
            msg = sprintf('El vector P converge a la clase %d.',k);
            disp(msg);
        end
    end
    %op = input('¿Desea probar otra entrada?:  ','s');
%end
for k=1:S
    fclose(archivosGrafica(k));
end
rango = 0:t+1;
valores = zeros(S,t+2);
for k=1:S
    nombre = sprintf('valoresDeGraficacion%d.txt',k);
    aux = dlmread(nombre); 
    valores(k,:) = aux';
end
figure
plot(rango,valores(1,:));
grid on
text(0,a1(1),'\leftarrow a2[1]');
title('Red de Hamming');
hold on
for k=2:S
    plot(rango,valores(k,:));
    cad = strcat('\leftarrow a2[',num2str(k),']');
    text(0,a1(k),cad);
end
hold off
xlabel('Iteraciones','FontWeight','bold')
ylabel('Salida de la capa recurrente','FontWeight','bold')
clear