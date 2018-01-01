function a = funcionDeActivacion(n,tipo)
    switch tipo
        case 1
            a = purelin(n);
        case 2
            a = logsig(n);
        case 3
            a = tansig(n);
    end
end