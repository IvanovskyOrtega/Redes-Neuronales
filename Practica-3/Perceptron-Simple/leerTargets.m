% Esta función lee los valores target t de un archivo .txt del conjunto de
% entrenamiento que se encuentren en un formato:
% {[p1 p2 p3 ... pn],[t1 t2 t3 ... tn]}
function A = leerTargets(nombre,targetDim,numProt,protDim)
archivo = fopen(nombre,'r');
while archivo == -1
    nombre = input('Ingrese un nombre válido:  ','s');
    archivo = fopen(nombre,'r');
end
A = zeros(numProt,targetDim);
for i=1:numProt
    fscanf(archivo,'{[');
    fscanf(archivo,' %d');
    fscanf(archivo,'],[');
    A(i,:) = fscanf(archivo,'%d');
    fscanf(archivo,']}\n');
end
fclose(archivo);
end
