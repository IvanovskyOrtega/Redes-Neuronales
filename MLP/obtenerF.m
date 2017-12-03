function F = obtenerF(tipo_funcion,num_neuronas,a)
    switch tipo_funcion
        case 1
            derivada = ones(1,num_neuronas);
            F = diag(derivada);
        case 2
            unos = ones(1,num_neuronas);
            derivada = diag(unos);
            derivada = derivada - diag(a');
            F = derivada*diag(a');
        case 3
            unos = ones(1,num_neuronas);
            derivada = diag(unos);
            F = derivada - (diag(a')*diag(a'));
    end
end