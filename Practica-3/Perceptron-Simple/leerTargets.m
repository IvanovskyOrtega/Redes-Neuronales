function A = leerTargets(nombre,targetDim,numProt,protDim)
    archivo = fopen(nombre,'r');
    A = zeros(numProt,targetDim);
    for i=1:numProt
       fscanf(archivo,'{[%d');
        for j=1:protDim-2
            fscanf(archivo,' %d');
        end
        fscanf(archivo,' %d],');
        A(i,1) = fscanf(archivo,'[%d');
        for j=1:targetDim-2
            A(i,j) = fscanf(archivo,' %d');
        end
        A(i,targetDim) = fscanf(archivo,' %d]},');
        fscanf(archivo,'\n');
    end
    fclose(archivo);
end