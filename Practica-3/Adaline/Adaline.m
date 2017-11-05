clc
clear
op = input('¿La red usará bias?:  ','s');
switch(op)
    case 's'
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
            S = 1;
            while numeroDeClases > power(2,n)
                n = n+1;
                S = S+1;
            end
            itmax = input('Itmax:  ');
            Eit = input('Eit:  ');
            corM = (1/numProt)*(prototipos*prototipos');
            eigenvalores = eig(corM);
            alfa = 1/max(eigenvalores(:));
            fprintf('Se recomienda un valor de alfa menor a %f\n',alfa);
            alfa = input('alpha =  ');
            W = randn(S,R);
            b = randn(S,1);
            disp('Se usarán las siguientes matrices iniciales:');
            disp('W =');
            disp(W);
            disp('b =');
            disp(b);
            aux = zeros(S,1);
            for i=1:itmax
                eit = zeros(S,1);
                for j=1:numProt
                    a = (W*prototipos(j,:)')+b;
                    ej = targets(j,:)'-a;
                    if ~isequal(ej,aux)
                        W = W+(2*alfa)*(ej*prototipos(j,:));
                        b = b+(2*alfa)*ej;
                        eit = eit+(1/numProt)*ej;
                    end
                end
                if isequal(eit,aux)
                    fprintf('Se tuvo un aprendizaje exitoso en la iteración %d.\n',i);
                    break;
                end
                %if sum(eit) < Eit
                %    fprintf('Se tuvo un aprendizaje con 3/4 de probabilidad de éxito en la iteración %d.\n',i);
                %    break;
                %end
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
            S = 1;
            while numeroDeClases > power(2,n)
                n = n+1;
                S = S+1;
            end
            itmax = input('Itmax:  ');
            Eit = input('Eit:  ');
            try
                corM = (1/numProt)*(prototipos*prototipos');
                eigenvalores = eig(corM);
                alfa = 1/max(eigenvalores(:));
                fprintf('Se recomienda un valor de alfa menor a %f\n',alfa);
            catch
            end
            alfa = input('alfa =  ');
            W = zeros(S,R);
            b = zeros(S,1);
            disp('Se usarán las siguientes matrices iniciales:');
            disp('W =');
            disp(W);
            disp('b =');
            disp(b);
            for i=1:itmax
                eit = 0;
                for j=1:numProt
                    a = W*prototipos(j,:)'+b;
                    ej = targets(j,:)'-a;
                    if ej ~= 0
                        W = W+(2*alfa*(ej*prototipos(j,:)));
                        b = b+(2*alfa)*ej;
                        eit = eit+(1/numProt)*ej;
                    end
                end
                if eit == Eit
                    fprintf('Se tuvo un aprendizaje exitoso en la iteración %d.\n',i);
                    break;
                end
            end
            if i == itmax
                disp('Es probable que los resultados finales no sean los correctos. Prueba con otro alfa.');
            end
            disp('W =');
            disp(W);
            disp('b =');
            disp(b);
            for i=1:numProt
                disp(W*prototipos(i,:)'+b);
            end
        end
    case 'n'
        disp('BCD');
        tam = input('Ingresa el tamaño del codificador:  ');
        max = 2^tam-1;
        targets = (0:max);
        bin = de2bi(targets);
        bin = fliplr(bin);
        disp('El conjunto de entrnamiento es: ');
        disp([bin targets']);
        R = max+1;
        S = 1;
        W = zeros(S,tam);
        disp('Se usará la siguiente matriz inicial de pesos.');
        disp('W =');
        disp(W);
        itmax = input('Itmax =  ');
        Eit = input('Eit =  ');
        alfa = input('alfa =  ');
        for i=1:itmax
            t = 0;
            eit = 0;
            fprintf('Iteracion %d\n',i);
            for j=1:R
                a = W*bin(j,:)';
                disp(a);
                ej = t-a;
                disp(ej);
                if ej ~= 0.0
                    W = W+(2*alfa*ej*bin(j,:));
                    eit = eit+(1/R)*ej;
                end
                t = t+1;
            end
            if eit == Eit
                fprintf('Se ha tenido un aprendizaje exitoso en la iteración %d.\n',i);
                break;
            end
        end
        disp('W =');
        disp(W);
    otherwise
end