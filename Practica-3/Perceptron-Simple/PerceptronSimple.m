clc 
clear
disp('PERCEPTR�N SIMPLE');
tipo = input('1)M�todo gr�fico\n2)Regla de aprendizaje\nIngrese su elecci�n:  ','s');
switch(tipo)
    case '1'
        disp('MG');
        numProt = input('Ingresa el n�mero de patrones prototipo:  ');
        res = input('�La dimensi�n es menor o igual a 3? (s/n):  ','s');
        protDim = input('Ingresa la dimensi�n de los vectores prototipo:  ');
        if strcmp(res,'s') || strcmp(res,'S')
            res = input('�La dimensi�n de los targets es mayor a 1?','s');
            if strcmp(res,'s') || strcmp(res,'S')
                targetDim = input('Ingresa la dimensi�n de los targets:  ');
                disp('Ingresa el nombre del archivo que contiene el Cto. de Entrenamiento\n');
                disp('(Sin extensi�n)');
                nombre = input('Nombre del Archivo:  ','s');
                nombre = strcat(nombre,'.txt');
                prototipos = leerPatronesPrototipo(nombre,targetDim,numProt,protDim);
                targets = leerTargets(nombre,targetDim,numProt,protDim);
                disp('Prototipos =');
                disp(prototipos);
                disp('Targets =');
                disp(targets);
            else
                disp('Ingresa el nombre del archivo que contiene el Cto. de Entrenamiento\n');
                disp('(Sin extensi�n)');
                nombre = input('Nombre del Archivo:  ','s');
                nombre = strcar(nombre,'.txt');
                ctoDeEntrenamiento = dlmread(nombre);
            end
        else
            exit(0);
        end
    case '2'
        disp('RDA');
    default
        display('Opci�n inv�lida');
end
