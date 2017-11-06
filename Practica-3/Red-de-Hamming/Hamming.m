clc
clear
% Se pide el número de vectores prototipo para asegurar la lectura correcta
% del archivo
numProt = input('Ingresa el número de vectores prototipo: ');

%De igual forma con el tamaño
r = input('Ingresa el tamaño (r) de los vectores prototipo:  ');

% Se pide el nombre del archivo sin extensión
nombreArchivo = input('Ingresa el nombre del archivo de los vectores prototipo (sin extension .txt):  ','s');
nombreArchivo = sprintf('%s,.txt',nombreArchivo);

% Se abre el archivo, si no lo puede abrir, solicita otro archivo válido
archivo = fopen(nombreArchivo,'r');
while (archivo == -1)
    disp('El archivo solicitado no existe o no se abrió correctamente.');
    nombreArchivo = input('Ingresa el nombre de un archivo válido:  ','s');
    archivo = fopen(nombreArchivo,'r');
end
disp('Se abrió correctamente el archivo.');
fclose(archivo);

% Si se abrió correctamente, se leen los prototipos y se cargan a la matriz
% de pesos.
W = dlmread(nombreArchivo);
dimensiones = size(W);

% El número de prototipos, es el número de neuronas de la red.
S = numProt;

% Si las dimensiones ingresadas no coinciden, se cambian por las correctas
% de acuerdo a lo que se leyó.
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

% Se muestra la matriz de pesos 
disp('W =');
disp(W);

% Se solicita el archivo que contiene los valores del vector bias.
archivoBias = input('Ingresa el nombre del archivo que contiene el vector bias (con extension .txt):  ','s');
archivo = fopen(archivoBias,'r');

% Si no existe, pide otro archivo válido.
while (archivo == -1)
    disp('El archivo solicitado no existe o no se abrió correctamente.');
    archivoBias = input('Ingresa el nombre de un archivo válido:  ','s');
    archivo = fopen(archivoBias,'r');
end
disp('Se abrió correctamente el archivo.');
fclose(archivo);

% Se cargan los valores del archivo al vector y se muestran en pantalla.
b = dlmread(archivoBias);
disp('b =');
disp(b);

% Se establece un valor aleatorio de epsilon para la matriz de pesos de la
% segunda capa.
epsilon = 0;
while epsilon == 0
    epsilon = mod(rand(),(1/S-1));
end
epsilonMatrix = ones(S,S)*epsilon;

% Se calcula la matriz de pesos de la segunda capa sumando las matrices
% triangulares de epsilon (inferior y superior) con una diagonal de 1's y
% se muestra en pantalla.
W_2 = tril(epsilonMatrix,-1)+triu(epsilonMatrix,1)+diag(ones(1,S));
disp('Se usará la siguiente matriz en la capa recurrente:');
disp(W_2);

% Se crea un vector de archivos para guardar los valores de graficación.
archivosGrafica = zeros(1,S);
for k=1:S
    nombre = sprintf('valoresDeGraficacion%d.txt',k);
    archivosGrafica(k) = fopen(nombre,'w');
end

% Se pide un vector de entrada para propagar en la red.
disp('Ingrese los valores del vector de entrada a propagar:');
P = zeros(r,1);
for i=1:r
    cad = sprintf('P(%d):',i);
    P(i) = input(cad);
end

% Inicia el aprendizaje
t = 0;
a1 = purelin(W*P+b);
a2_1 = a1;

% Se guardan los valores de salida de la 1ra capa a aun archivo
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
    
    % Se guardan las salidas de la capa recurrente.
    for k=1:S
        fprintf(archivosGrafica(k),'%f\r\n',a2_1(k));
    end
    
    t = t+1;
end

% Debido a que la salida puede ser un número real, se redondea el valor
% para obtener el entero más cercano y se muestra en pantalla.
a2_2 = round(a2_2);
disp('a2 = ');
disp(a2_2);

% Se guarda el resultado final
for k=1:S
    fprintf(archivosGrafica(k),'%f\r\n',a2_2(k));
end

% Se muestra al usuario en cuantas iteraciones convergió la entrada.
msg = sprintf('La red convergió exitosamente en %d iteraciones.',t);
disp(msg);

% Se busca a que clase pertenece el vector de entrada.
for k=1:S
    if a2_2(k) ~= 0
        msg = sprintf('El vector P converge a la clase %d.',k);
        disp(msg);
    end
end

% Se cierran los archivos de graficación
for k=1:S
    fclose(archivosGrafica(k));
end

% Inicia la graficación
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
