clc
clear
disp('PERCEPTRÓN SIMPLE');

% El método gráfico no se pudo implementar, visitar el siguiete repositorio
% para ver una versión (No muy gráfica) de este método.
tipo = input('1)Método gráfico\n2)Regla de aprendizaje\nIngrese su elección:  ','s');
switch(tipo)
    case '1'
        disp('No implementado :(');
    case '2'
        disp('REGLA DE APRENDIZAJE');
        % Se solicita el número de patrones prototipo y su dimensión.
        numProt = input('Ingresa el número de patrones prototipo:  ');
        R = input('Ingresa la dimensión de los vectores prototipo:  ');
        
        % Se pregunta si la dimensión de los targets es mayor a 2, debido
        % al formato especificado en la práctica.
        res = input('¿La dimensión de los targets es mayor a 1? (s/n):  ','s');
        
        % Si es mayor a 1
        if strcmp(res,'s') || strcmp(res,'S')
            % Se pide la dimensión de los targets
            targetDim = input('Ingresa la dimensión de los targets:  ');
            
            % Se solicita el nombre del archivo que contiene el conjunto de
            % entrenamiento sin ingresar la extensión (
            % Se usa .txt por defecto.) 
            disp('Ingresa el nombre del archivo que contiene el Cto. de Entrenamiento\n');
            disp('(Sin extensión)');
            nombre = input('Nombre del Archivo:  ','s');
            nombre = strcat(nombre,'.txt');
            
            % Se leen los patrones prototipo y los targets del archivo
            % indicado y se muestran en pantalla.
            prototipos = leerPatronesPrototipo(nombre,targetDim,numProt,R);
            targets = leerTargets(nombre,targetDim,numProt,R);
            disp('Prototipos =');
            disp(prototipos);
            disp('Targets =');
            disp(targets);
            
            % Se solicita al usuario el número de clases a las que va a
            % clasificar para calcular el número de neuronas. Recordando
            % que S neuronas separan 2^S clases.
            numeroDeClases = input('Ingresa el número de clases:  ');
            n = 1;
            S = 1.0;
            while numeroDeClases > power(2,n)
                n = n+1;
                S = S+1.0;
            end
            
            % Se solicitan los valores de itmax y Eit
            itmax = input('Itmax:  ');
            Eit = input('Eit:  ');
            Eit = ones(S,1)*Eit;
            
            % Se inicializan la matriz de pesos y el vector bias con
            % valores aleatorios y se muestran en pantalla.
            W = randn(S,R);
            b = randn(S,1);
            disp('Se usarán las siguientes matrices iniciales:');
            disp('W =');
            disp(W);
            disp('b =');
            disp(b);
            
            % Se usa un vector auxiliar para comparar con la señal del 
            % error.
            aux = zeros(S,1);
            
            % Se abren los archivos necesarios para guardar los valores de
            % graficación.
            pesos = zeros(1,S*R);
            bias = zeros(1,S);
            errores = zeros(1,S);
            for k=1:S*R
                nombre = sprintf('pesos%d.txt',k);
                pesos(k) = fopen(nombre,'w');
            end
            for k=1:S
                nombre = sprintf('bias%d.txt',k);
                bias(k) = fopen(nombre,'w');
            end
            for k=1:S
                nombre = sprintf('error%d.txt',k);
                errores(k) = fopen(nombre,'w');
            end
            
            % Se imprimen a archivo los valores iniciales.
            k=1;
            for l=1:S
                for m=1:R
                    fprintf(pesos(k),'%f\r\n',W(l,m));
                    k = k+1;
                end
            end
            for l=1:S
                fprintf(bias(l),'%f\r\n',b(l));
            end
            for l=1:S
                fprintf(errores(l),'%f\r\n',0);
            end
            
            % Inicia el aprendizaje
            for i=1:itmax
                eit = zeros(S,1);
                t = 0;
                for j=1:numProt
                    a = hardlim(W*prototipos(j,:)'+b);
                    ej = targets(j,:)'-a;
                    if ~isequal(ej,aux)
                        t = t+1;
                        W = W+ej*prototipos(j,:);
                        b = b+ej;
                        eit = eit+(1/numProt)*ej;
                    end
                end
                
                % Se imprimen valores a archivo
                k=1;
                for l=1:S
                    for m=1:R
                        fprintf(pesos(k),'%f\r\n',W(l,m));
                        k = k+1;
                    end
                end
                for l=1:S
                    fprintf(bias(l),'%f\r\n',b(l));
                end
                for l=1:S
                    fprintf(errores(l),'%f\r\n',eit(l));
                end
                
                % 1er criterio de finalización: Si la señal del error es 0.
                if isequal(eit,Eit) && t == 0
                    fprintf('Se tuvo un aprendizaje exitoso en la iteración %d.\n',i);
                    break;
                end
                
                % 2dor criterio de finalización: Si la señal del error es
                % mayor a 0 y menor a Eit.
                if (1/S)*sum(eit) < (1/S)*sum(Eit) && (1/S)*sum(eit) > 0
                   fprintf('Se tuvo un aprendizaje con 3/4 de probabilidad de éxito en la iteración %d.\n',i);
                    break;
                end
            end
            itmax = i;
            
            % Se muestran los valores finales de la matriz de pesos y el 
            % vector bias.
            disp('W =');
            disp(W);
            disp('b =');
            disp(b);
            
            % Se cierran los archivos
            for k=1:S*R
                fclose(pesos(k));
            end
            for k=1:S
                fclose(bias(k));
            end
            for k=1:S
                fclose(errores(k));
            end
            
            % Inicia la graficación
            numeroDeArchivos = S*R;
            valores = zeros(numeroDeArchivos,itmax+1);
            for k=1:numeroDeArchivos
                nombre = sprintf('pesos%d.txt',k);
                aux = dlmread(nombre);
                valores(k,:) = aux';
            end
            figure
            rango = (0:itmax);
            plot(rango,valores(1,:));
            grid on
            text(itmax,valores(1,itmax+1),'\leftarrow W[1]');
            title('Perceptron Simple');
            hold on
            for k=2:numeroDeArchivos
                plot(rango,valores(k,:));
                cad = strcat('\leftarrow W[',num2str(k),']');
                text(itmax,valores(k,itmax+1),cad);
            end
            hold off
            xlabel('Iteraciones','FontWeight','bold')
            ylabel('Evolución de los pesos','FontWeight','bold')
            numeroDeArchivos = S;
            valores = zeros(numeroDeArchivos,itmax+1);
            for k=1:numeroDeArchivos
                nombre = sprintf('bias%d.txt',k);
                aux = dlmread(nombre);
                valores(k,:) = aux';
            end
            figure
            plot(rango,valores(1,:));
            grid on
            text(itmax,valores(1,itmax+1),'\leftarrow b[1]');
            title('Perceptron Simple');
            hold on
            for k=2:numeroDeArchivos
                plot(rango,valores(k,:));
                cad = strcat('\leftarrow b[',num2str(k),']');
                text(itmax,valores(k,itmax+1),cad);
            end
            hold off
            xlabel('Iteraciones','FontWeight','bold')
            ylabel('Evolución del vector bias','FontWeight','bold')
            numeroDeArchivos = S;
            valores = zeros(numeroDeArchivos,itmax+1);
            for k=1:numeroDeArchivos
                nombre = sprintf('error%d.txt',k);
                aux = dlmread(nombre);
                valores(k,:) = aux';
            end
            figure
            plot(rango,valores(1,:));
            grid on
            text(itmax,valores(1,itmax+1),'\leftarrow eit[1]');
            title('Perceptron Simple');
            hold on
            for k=2:numeroDeArchivos
                plot(rango,valores(k,:));
                cad = strcat('\leftarrow eit[',num2str(k),']');
                text(itmax,valores(k,itmax+1),cad);
            end
            hold off
            xlabel('Iteraciones','FontWeight','bold')
            ylabel('Señal del error','FontWeight','bold')
            imprimirMatricesEnArchivo(W,b,'c');
        
        % Si la dimensión de los targets es 1.
        else 
            % Se solicita el archivo que contiene el conjunto de
            % entrenamiento.
            disp('Ingresa el nombre del archivo que contiene el Cto. de Entrenamiento\n');
            disp('(Sin extensión)');
            nombre = input('Nombre del Archivo:  ','s');
            nombre = strcat(nombre,'.txt');
            ctoDeEntrenamiento = dlmread(nombre);
            dimensiones = size(ctoDeEntrenamiento);
            R = dimensiones(2)-1;
            
            % Se leen los vectores prototipo y los targets y se muestran en
            % pantalla.
            prototipos = ctoDeEntrenamiento(:,1:R);
            targets = ctoDeEntrenamiento(:,R+1);
            disp('Prototipos =');
            disp(prototipos);
            disp('Targets =');
            disp(targets);
            
            % Se solicita el numero de clases a clasificar para calcular el
            % número de neuronas de la red.
            numeroDeClases = input('Ingresa el número de clases:  ');
            n = 1;
            S = 1.0;
            while numeroDeClases > power(2,n)
                n = n+1;
                S = S+1.0;
            end
            
            % Se solicitan los valores itmax y Eit
            itmax = input('Itmax:  ');
            Eit = input('Eit:  ');
            
            % Se inicializan la matriz de pesos y el vector bias con
            % valores aleatorios de la distribución normal y se muestran en
            % pantalla.
            W = randn(S,R);
            b = randn(S,1);
            disp('Se usarán las siguientes matrices iniciales:');
            disp('W =');
            disp(W);
            disp('b =');
            disp(b);
            
            % Se abren los archivos necesarios para guardar los valores de
            % graficación.
            pesos = zeros(1,S*R);
            bias = zeros(1,S);
            errores = zeros(1);
            for k=1:S*R
                nombre = sprintf('pesos%d.txt',k);
                pesos(k) = fopen(nombre,'w');
            end
            for k=1:S
                nombre = sprintf('bias%d.txt',k);
                bias(k) = fopen(nombre,'w');
            end
            for k=1:S
                nombre = sprintf('error%d.txt',k);
                errores(k) = fopen(nombre,'w');
            end
            
            % Se imprimen a archivo los valores iniciales
            k=1;
            for l=1:S
                for m=1:R
                    fprintf(pesos(k),'%f\r\n',W(l,m));
                    k = k+1;
                end
            end
            for l=1:S
                fprintf(bias(l),'%f\r\n',b(l));
            end
            fprintf(errores,'%f\r\n',0);
            
            % Inicia el aprendizaje
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
                
                % Se imprimen valores a archivo
                k=1;
                for l=1:S
                    for m=1:R
                        fprintf(pesos(k),'%f\r\n',W(l,m));
                        k = k+1;
                    end
                end
                for l=1:S
                    fprintf(bias(l),'%f\r\n',b(l));
                end
                for l=1:S
                    fprintf(errores(l),'%f\r\n',eit(l));
                end
                
                % 1er criterio de finalización
                if isequal(eit,Eit) && t == 0
                    fprintf('Se tuvo un aprendizaje exitoso en la iteración %d.\n',i);
                    break;
                end
                
                % Segundo criterio de finalización
                if (1/S)*sum(eit) < (1/S)*sum(Eit) && (1/S)*sum(eit) > 0
                   fprintf('Se tuvo un aprendizaje con 3/4 de probabilidad de éxito en la iteración %d.\n',i);
                    break;
                end
            end
            itmax = i;
            
            % Se muestran la matriz de pesos y el vector bias finales.
            disp('W =');
            disp(W);
            disp('b =');
            disp(b);
            
            % Se cierran los archivos de graficación.
            for k=1:S*R
                fclose(pesos(k));
            end
            for k=1:S
                fclose(bias(k));
            end
            for k=1:S
                fclose(errores(k));
            end
            
            % Inicia la graficación
            numeroDeArchivos = S*R;
            valores = zeros(numeroDeArchivos,itmax+1);
            for k=1:numeroDeArchivos
                nombre = sprintf('pesos%d.txt',k);
                aux = dlmread(nombre);
                valores(k,:) = aux';
            end
            figure
            rango = (0:itmax);
            plot(rango,valores(1,:));
            grid on
            text(itmax,valores(1,itmax+1),'\leftarrow W[1]');
            title('Perceptron Simple');
            hold on
            for k=2:numeroDeArchivos
                plot(rango,valores(k,:));
                cad = strcat('\leftarrow W[',num2str(k),']');
                text(itmax,valores(k,itmax+1),cad);
            end
            hold off
            xlabel('Iteraciones','FontWeight','bold')
            ylabel('Evolución de los pesos','FontWeight','bold')
            numeroDeArchivos = S;
            valores = zeros(numeroDeArchivos,itmax+1);
            for k=1:numeroDeArchivos
                nombre = sprintf('bias%d.txt',k);
                aux = dlmread(nombre);
                valores(k,:) = aux';
            end
            figure
            plot(rango,valores(1,:));
            grid on
            text(itmax,valores(1,itmax+1),'\leftarrow b[1]');
            title('Perceptron Simple');
            hold on
            for k=2:numeroDeArchivos
                plot(rango,valores(k,:));
                cad = strcat('\leftarrow b[',num2str(k),']');
                text(itmax,valores(k,itmax+1),cad);
            end
            hold off
            xlabel('Iteraciones','FontWeight','bold')
            ylabel('Evolución del vector bias','FontWeight','bold')
            numeroDeArchivos = S;
            valores = zeros(numeroDeArchivos,itmax+1);
            for k=1:numeroDeArchivos
                nombre = sprintf('error%d.txt',k);
                aux = dlmread(nombre);
                valores(k,:) = aux';
            end
            figure
            plot(rango,valores(1,:));
            grid on
            text(itmax,valores(1,itmax+1),'\leftarrow eit[1]');
            title('Perceptron Simple');
            hold on
            for k=2:numeroDeArchivos
                plot(rango,valores(k,:));
                cad = strcat('\leftarrow eit[',num2str(k),']');
                text(itmax,valores(k,itmax+1),cad);
            end
            hold off
            xlabel('Iteraciones','FontWeight','bold')
            ylabel('Señal del error','FontWeight','bold')
            imprimirMatricesEnArchivo(W,b,'c');
        end
    otherwise
        display('Opción inválida');
end
