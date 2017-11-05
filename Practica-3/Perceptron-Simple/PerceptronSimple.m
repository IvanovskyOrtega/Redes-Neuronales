clc
clear
disp('PERCEPTRÓN SIMPLE');
tipo = input('1)Método gráfico\n2)Regla de aprendizaje\nIngrese su elección:  ','s');
switch(tipo)
    case '1'
        disp('MG');
        numProt = input('Ingresa el número de patrones prototipo:  ');
        res = input('¿La dimensión es menor o igual a 3? (s/n):  ','s');
        R = input('Ingresa la dimensión de los vectores prototipo:  ');
        while R > 4 && R < 1
            R = input('Ingresa una dimensión adecuada:  ');
        end
        if strcmp(res,'s') || strcmp(res,'S')
            res = input('¿La dimensión de los targets es mayor a 1? (s/n):  ','s');
            if strcmp(res,'s') || strcmp(res,'S')
                targetDim = input('Ingresa la dimensión de los targets:  ');
                disp('Ingresa el nombre del archivo que contiene el Cto. de Entrenamiento\n');
                disp('(Sin extensión)');
                nombre = input('Nombre del Archivo:  ','s');
                nombre = strcat(nombre,'.txt');
                prototipos = leerPatronesPrototipo(nombre,targetDim,numProt,R);
                targets = leerTargets(nombre,targetDim,numProt,R);
                disp('Prototipos =');
                disp(prototipos);
                disp('Targets =');
                disp(targets);
                numeroDeClases = input('Ingresa el número de clases:  ');
                fprintf('Numero de clases: %d\n',numeroDeClases);
            else
                disp('Ingresa el nombre del archivo que contiene el Cto. de Entrenamiento\n');
                disp('(Sin extensión)');
                nombre = input('Nombre del Archivo:  ','s');
                nombre = strcat(nombre,'.txt');
                ctoDeEntrenamiento = dlmread(nombre);
                disp(ctoDeEntrenamiento);
                dimensiones = size(ctoDeEntrenamiento);
                R = dimensiones(2)-1;
                prototipos = ctoDeEntrenamiento(:,1:R);
                targets = ctoDeEntrenamiento(:,R-1:R);
                disp(prototipos);
                disp(targets);
                numeroDeClases = input('Ingresa el número de clases:  ');
                n = 1;
                S = 1;
                while numeroDeClases > pow(2,n)
                    n = n+1;
                    S = S+1;
                end
                
            end
        else
            disp('El método gráfico solo funciona para vecotres prototipo de dimensión menor o igual a 3.');
            exit(0);
        end
    case '2'
        disp('RDA');
        numProt = input('Ingresa el número de patrones prototipo:  ');
        R = input('Ingresa la dimensión de los vectores prototipo:  ');
        res = input('¿La dimensión de los targets es mayor a 1? (s/n):  ','s');
        if strcmp(res,'s') || strcmp(res,'S')
            targetDim = input('Ingresa la dimensión de los targets:  ');
            disp('Ingresa el nombre del archivo que contiene el Cto. de Entrenamiento\n');
            disp('(Sin extensión)');
            nombre = input('Nombre del Archivo:  ','s');
            nombre = strcat(nombre,'.txt');
            prototipos = leerPatronesPrototipo(nombre,targetDim,numProt,R);
            targets = leerTargets(nombre,targetDim,numProt,R);
            disp('Prototipos =');
            disp(prototipos);
            disp('Targets =');
            disp(targets);
            numeroDeClases = input('Ingresa el número de clases:  ');
            n = 1;
            S = 1.0;
            while numeroDeClases > power(2,n)
                n = n+1;
                S = S+1.0;
            end
            itmax = input('Itmax:  ');
            Eit = input('Eit:  ');
            W = ones(S,R)*(mod(rand(),2));
            b = ones(S,1)*(mod(rand(),2));
            disp('Se usarán las siguientes matrices iniciales:');
            disp('W =');
            disp(W);
            disp('b =');
            disp(b);
            for i=1:itmax
                eit = 0;
                t = 0;
                for j=1:numProt
                    a = hardlim(W*prototipos(j,:)'+b);
                    disp(a);
                    ej = targets(j,:)'-a;
                    disp(ej);
                    if ej ~= 0
                        t = t+1;
                        W = W+ej*prototipos(j,:);
                        b = b+ej;
                        eit = eit+(1/numProt)*ej;
                    end
                end
                if isequal(eit,Eit) && t == 0
                    fprintf('Se tuvo un aprendizaje exitoso en la iteración %d.\n',i);
                    break;
                end
            end
            disp('W =');
            disp(W);
            disp('b =');
            disp(b);
        else
            disp('Ingresa el nombre del archivo que contiene el Cto. de Entrenamiento\n');
            disp('(Sin extensión)');
            nombre = input('Nombre del Archivo:  ','s');
            nombre = strcat(nombre,'.txt');
            ctoDeEntrenamiento = dlmread(nombre);
            dimensiones = size(ctoDeEntrenamiento);
            R = dimensiones(2)-1;
            prototipos = ctoDeEntrenamiento(:,1:R);
            targets = ctoDeEntrenamiento(:,R+1);
            disp('Prototipos =');
            disp(prototipos);
            disp('Targets =');
            disp(targets);
            numeroDeClases = input('Ingresa el número de clases:  ');
            n = 1;
            S = 1.0;
            while numeroDeClases > power(2,n)
                n = n+1;
                S = S+1.0;
            end
            itmax = input('Itmax:  ');
            Eit = input('Eit:  ');
            W = ones(S,R)*(mod(rand(),2));
            b = ones(S,1)*(mod(rand(),2));
            disp('Se usarán las siguientes matrices iniciales:');
            disp('W =');
            disp(W);
            disp('b =');
            disp(b);
            for i=1:itmax
                eit = 0;
                t = 0;
                for j=1:numProt
                    a = hardlim(W*prototipos(j,:)'+b);
                    ej = targets(j,:)-a;
                    if ej ~= 0
                        t = t+1;
                        W = W+ej*prototipos(j,:);
                        b = b+ej;
                        eit = eit+(1/numProt)*ej;
                    end
                end
                if eit <= Eit && t == 0
                    fprintf('Se tuvo un aprendizaje exitoso en la iteración %d.\n',i);
                    break;
                end
            end
            disp('W =');
            disp(W);
            disp('b =');
            disp(b);
        end
    otherwise
        display('Opción inválida');
end
