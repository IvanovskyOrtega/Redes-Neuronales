function A = leerPatronesPrototipo(nombre,targetDim,numProt,protDim)
    archivo = fopen(nombre,'r');
    A = zeros(numProt,protDim);
    while archivo == -1
        nombre = input('Ingrese un nombre válido:  ','s');
    archivo = fopen(nombre,'r');
    end
    for i=1:numProt
        A(i,1) = fscanf(archivo,'{[%d');
        for j=1:protDim-2
            A(i,j) = fscanf(archivo,' %d');
        end
        A(i,protDim) = fscanf(archivo,' %d],');
        fscanf(archivo,'[%d');
        for j=1:targetDim-2
            fscanf(archivo,' %d');
        end
        fscanf(archivo,' %d]},');
        fscanf(archivo,'\n');
    end
    fclose(archivo);
end