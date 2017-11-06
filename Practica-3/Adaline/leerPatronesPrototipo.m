% Esta función lee los valores de los vectores prototipo p de un archivo 
% .txt del conjunto de entrenamiento que se encuentren en un formato:
% {[p1 p2 p3 ... pn],[t1 t2 t3 ... tn]}
function A = leerPatronesPrototipo(nombre,targetDim,numProt,protDim)
    archivo = fopen(nombre,'r');
    A = zeros(numProt,protDim);
    disp([numProt protDim]);
    while archivo == -1
        nombre = input('Ingrese un nombre válido:  ','s');
        archivo = fopen(nombre,'r');
    end
    for i=1:numProt
        fscanf(archivo,'{[');
        A(i,:) = fscanf(archivo,'%d');
        fscanf(archivo,'],[');
        fscanf(archivo,'%d');
        fscanf(archivo,']}\n');
    end
    fclose(archivo);
end
