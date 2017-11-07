function imprimirMatricesEnArchivo(W,b,tipo)
    switch(tipo)
        case 'c'
            imprimirConBias(W,b)
        case 's'
            imprimirSinBias(W)
    end
end

function imprimirConBias(W,b)
t = clock;
nombre = sprintf('resultados_%d_%d_%d_%d_%d.txt',t(4),t(5),t(3),t(2),t(1));
resultados = fopen(nombre,'w');
fprintf(resultados,'W = \r\n');
fclose(resultados);
dlmwrite(nombre,W,'-append','delimiter','\t');
resultados = fopen(nombre,'a');
fprintf(resultados,'\r\n');
fprintf(resultados,'b = \r\n');
fclose(resultados);
dlmwrite(nombre,b,'-append','delimiter','\t');
end

function imprimirSinBias(W)
t = clock;
nombre = sprintf('resultados_%d_%d_%d_%d_%d.txt',t(4),t(5),t(3),t(2),t(1));
resultados = fopen(nombre,'w');
fprintf(resultados,'W = \r\n');
fclose(resultados);
dlmwrite(nombre,W,'-append','delimiter','\t');
end