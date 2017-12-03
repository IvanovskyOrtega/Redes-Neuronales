function M = obtenerConjuntoDeValidacion(entradas,targets,num_datos,num_elem)
M = zeros(num_elem,2);
inc = round(num_datos/(num_elem+1));
j = 1;
for i=inc:inc:num_datos
    M(j,1) = entradas(i);
    M(j,2) = targets(i);
    j = j+1;
    if j == num_elem+1
        break;
    end
end
end