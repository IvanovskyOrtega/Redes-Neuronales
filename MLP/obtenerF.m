function F = obtenerF(tipo_funcion,num_neuronas,a)
    switch tipo_funcion
        case 1
            F = diag(ones(1,num_neuronas));
        case 2
            F = diag(logsig('dn',a,a));
        case 3
            F = diag(tansig('dn',a,a));
    end
end
