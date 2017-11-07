clc
clear

% Se solicita al usuario si usará una red ADALINE con bias o sin bias.
op = input('¿La red usará bias? (s/n):  ','s');
switch(op)
    
    % Si es con bias
    case 's'
        
        % Se solicita el número de patrones prototipo para una lectura
        % correcta de estos mismos en el archivo, además de su dimensión.
        numProt = input('Ingresa el número de patrones prototipo:  ');
        R = input('Ingresa la dimensión de los vectores prototipo:  ');
        
        % Se solicita indicar si la dimensión de los targets es mayor a 1,
        % para cambiar el modo de lectura de los patrones prototipo y
        % targets según el formato establecido en la práctica.
        res = input('¿La dimensión de los targets es mayor a 1? (s/n):  ','s');
        
        % Si es mayor.
        if strcmp(res,'s') || strcmp(res,'S')
            
            % Se pide la dimensión de los targets.
            targetDim = input('Ingresa la dimensión de los targets:  ');
            
            % Se solicita el archivo del conjunto de entrenamiento para la
            % lectura de los patrones prototipo y los targets sin
            % extensión. 
            % Se usa la extensión .txt por defecto.
            disp('Ingresa el nombre del archivo que contiene el Cto. de Entrenamiento\n');
            disp('(Sin extensión)');
            nombre = input('Nombre del Archivo:  ','s');
            nombre = strcat(nombre,'.txt');
            
            % Se leen los patrones prototipo y los targets y se muestran en
            % pantalla.
            prototipos = leerPatronesPrototipo(nombre,targetDim,numProt,R);
            targets = leerTargets(nombre,targetDim,numProt,R);
            disp('Prototipos =');
            disp(prototipos);
            disp('Targets =');
            disp(targets);
            
            % Se solicita el número de clases a clasificar para calcular el
            % número de neuronas de la red.
            numeroDeClases = input('Ingresa el número de clases:  ');
            n = 1;
            S = 1;
            while numeroDeClases > power(2,n)
                n = n+1;
                S = S+1;
            end
            
            % Se solicitan los valores itmax y Eit.
            itmax = input('Itmax:  ');
            Eit = input('Eit:  ');
            
            % Se intenta hacer un cálculo de un alfa adecuado, según el
            % algoritmo de aprendizaje LMS.
            try
                corM = (1/numProt)*(prototipos*prototipos');
                eigenvalores = eig(corM);
                alfa = 1/max(eigenvalores(:));
                fprintf('Se recomienda un valor de alfa menor a %f\n',alfa);
            catch
            end
            
            % Se solicita el valor de alfa o factor de aprendizaje.
            alfa = input('alpha =  ');
            
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
            
            % Se usa un vector de 0's auxiliar para comparar con la señal
            % del error en el aprendizaje de la red.
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
                for j=1:numProt
                    a = (W*prototipos(j,:)')+b;
                    ej = targets(j,:)'-a;
                    if ~isequal(ej,aux)
                        W = W+(2*alfa)*(ej*prototipos(j,:));
                        b = b+(2*alfa)*ej;
                        eit = eit+(1/numProt)*ej;
                    end
                end
                
                % Se imprimen valores a archivo.
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
                
                % 1er criteri de finalización. Si el error es 0.
                if isequal(eit,aux)
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
            
            % Se cierran los archivos de los valores de graficación.
            for k=1:S*R
                fclose(pesos(k));
            end
            for k=1:S
                fclose(bias(k));
            end
            for k=1:S
                fclose(errores(k));
            end
            
            % Inicia la graficación.
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
            title('ADALINE');
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
            title('ADALINE');
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
            title('ADALINE');
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
            
            %Se solicita el archivo que contiene el conjunto de
            %entrenamiento para leer los patrones prototipo y los targets 
            % sin extensión.
            disp('Ingresa el nombre del archivo que contiene el Cto. de Entrenamiento\n');
            disp('(Sin extensión)');
            nombre = input('Nombre del Archivo:  ','s');
            nombre = strcat(nombre,'.txt');
            ctoDeEntrenamiento = dlmread(nombre);
            dimensiones = size(ctoDeEntrenamiento);
            R = dimensiones(2)-1;
            
            % Se leen los patrones prototipo y los targets y se muestran en
            % pantalla.
            prototipos = ctoDeEntrenamiento(:,1:R);
            targets = ctoDeEntrenamiento(:,R+1);
            disp('Prototipos =');
            disp(prototipos);
            disp('Targets =');
            disp(targets);
            
            % Se solicita el número de clases a clasificar para el cálculo
            % del número de neuronas de la red.
            numeroDeClases = input('Ingresa el número de clases:  ');
            n = 1;
            S = 1;
            while numeroDeClases > power(2,n)
                n = n+1;
                S = S+1;
            end
            
            % Se solicitan los valores itmax y Eit
            itmax = input('Itmax:  ');
            Eit = input('Eit:  ');
            
            % Se intenta hacer el cálculo de un valor de alfa adecuado de
            % acuerdo al algoritmo de aprendizaje LSM.
            try
                corM = (1/numProt)*(prototipos*prototipos');
                eigenvalores = eig(corM);
                alfa = 1/max(eigenvalores(:));
                fprintf('Se recomienda un valor de alfa menor a %f\n',alfa);
            catch
            end
            
            % Se solicita el valor de alfa.
            alfa = input('alfa =  ');
            
            % Se inicializan la matriz de pesos y el vector bias con ceros
            % en sus entradas y se muestran en pantalla.
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
            
            % Se imprimen los valores iniciales.
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
                if eit == Eit
                    fprintf('Se tuvo un aprendizaje exitoso en la iteración %d.\n',i);
                    break;
                end
                % Segundo criterio de finalización
                if eit < Eit && eit > 0
                   fprintf('Se tuvo un aprendizaje con 3/4 de probabilidad de éxito en la iteración %d.\n',i);
                    break;
                end
            end
            
            % Si se llegó al número máximo de iteraciones sin un
            % aprendizaje exitoso, se indica al usuario que posiblemente
            % los valores finales no son los adecuados.
            if i == itmax
                disp('Es probable que los resultados finales no sean los correctos. Prueba con otro alfa.');
            end
            itmax = i;
            
            % Se muestran la matriz de pesos y el vector bias finales.
            disp('W =');
            disp(W);
            disp('b =');
            disp(b);
            
            % Se cierran los archivos de valores de graficación.
            for k=1:S*R
                fclose(pesos(k));
            end
            for k=1:S
                fclose(bias(k));
            end
            for k=1:S
                fclose(errores(k));
            end
            
            % Inicia la graficación.
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
            title('ADALINE');
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
            title('ADALINE');
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
            title('ADALINE');
            hold on
            for k=2:numeroDeArchivos
                plot(rango,valores(k,:));
                cad = strcat('\eit[',num2str(k),']');
                text(itmax,valores(k,itmax+1),cad);
            end
            hold off
            xlabel('Iteraciones','FontWeight','bold')
            ylabel('Señal del error','FontWeight','bold')
            imprimirMatricesEnArchivo(W,b,'c');
        end
    
    % Si es una red sin bias, se resuelve el BCD de tamaño n.
    case 'n'
        disp('BCD');
        
        % Se solicita el tamaño del codificador y se calcula el valor
        % máximo decimal para el codificador.
        tam = input('Ingresa el tamaño del codificador:  ');
        max = 2^tam-1;
        
        % Se crea el vector de targets
        targets = (0:max);
        
        % Se crea la matriz correspondiente a los targets en binario y se
        % invierte la matriz para que los prototipos empiecen desde el bit
        % más significativo.
        bin = de2bi(targets);
        bin = fliplr(bin);
        
        % Se muestra el conjunto de entrenamiento
        disp('El conjunto de entrnamiento es: ');
        disp([bin targets']);
        
        % Se establecen los valores del tamaño de los vectores prototipo,
        % el número de neuronas de la red.
        R = max+1;
        S = 1;
        
        % Se inicializa la matriz de pesos con ceros en sus entradas y se
        % muestra en pantalla
        W = zeros(S,tam);
        disp('Se usará la siguiente matriz inicial de pesos.');
        disp('W =');
        disp(W);
        
        % Se solicitan los valores de itmax, Eit y alfa (factor de
        % aprendizaje).
        itmax = input('Itmax =  ');
        Eit = input('Eit =  ');
        alfa = input('alfa =  ');
        
        % Se abren los archivos necesarios para guardar los valores de
        % graficación.
        pesos = zeros(1,S*tam);
        errores = zeros(1,S);
        for k=1:S*tam
            nombre = sprintf('pesos%d.txt',k);
            pesos(k) = fopen(nombre,'w');
        end
        for k=1:S
            nombre = sprintf('error%d.txt',k);
            errores(k) = fopen(nombre,'w');
        end
        
        % Se imprimen a archivo los valores iniciales.
        k=1;
        for l=1:S
            for m=1:tam
                fprintf(pesos(k),'%f\r\n',W(l,m));
                k = k+1;
            end
        end
        for l=1:S
            fprintf(errores(l),'%f\r\n',0);
        end
        
        % Inicia el aprendizaje de la red.
        for i=1:itmax
            t = 0;
            eit = 0;
            for j=1:R
                a = W*bin(j,:)';
                ej = t-a;
                if ej ~= 0.0
                    W = W+(2*alfa*ej*bin(j,:));
                    eit = eit+(1/R)*ej;
                end
                t = t+1;
            end
            
            % Se imprimen valroes a archivo.
            k=1;
            for l=1:S
                for m=1:tam
                    fprintf(pesos(k),'%f\r\n',W(l,m));
                    k = k+1;
                end
            end
            for l=1:S
                fprintf(errores(l),'%f\r\n',eit(l));
            end
            
            %1er criterio de finalización.
            if eit == Eit
                fprintf('Se ha tenido un aprendizaje exitoso en la iteración %d.\n',i);
                break;
            end
            % Segundo criterio de finalización
            if eit < Eit && eit > 0
               fprintf('Se tuvo un aprendizaje con 3/4 de probabilidad de éxito en la iteración %d.\n',i);
               break;
            end
        end
        itmax = i;
        
        % Se muestra la matriz de pesos final.
        disp('W =');
        disp(W);
        
        % Se cierran los archivos de valores de graficación.
        for k=1:S*tam
            fclose(pesos(k));
        end
        for k=1:S
            fclose(errores(k));
        end
        
        % Inicia la graficación.
        numeroDeArchivos = S*tam;
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
        title('ADALINE');
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
            nombre = sprintf('error%d.txt',k);
            aux = dlmread(nombre);
            valores(k,:) = aux';
        end
        figure
        plot(rango,valores(1,:));
        grid on
        text(itmax,valores(1,itmax+1),'\leftarrow eit[1]');
        title('ADALINE');
        hold on
        for k=2:numeroDeArchivos
            plot(rango,valores(k,:));
            cad = strcat('\leftarrow eit');
            text(itmax,valores(k,itmax+1),cad);
        end
        hold off
        xlabel('Iteraciones','FontWeight','bold')
        ylabel('Señal del error','FontWeight','bold')
        imprimirMatricesEnArchivo(W,0,'s');
    otherwise
end
