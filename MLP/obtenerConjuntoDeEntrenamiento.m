function M = obtenerConjuntoDeEntrenamiento(entradas,targets,num_datos,num_elem_ent,cto_val,cto_prueba)
M = zeros(num_elem_ent,2);
j = 1;
for i=1:num_datos
    if ismember(entradas(i),cto_val(:,1)) || ismember(entradas(i),cto_prueba(:,1))
    else
        M(j,1) = entradas(i);
        M(j,2) = targets(i);
        j = j+1;
    end
    if j == num_elem_ent+1
        break;
    end
end
end