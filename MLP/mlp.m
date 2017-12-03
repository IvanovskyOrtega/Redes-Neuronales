clc
clear

% Leemos el archivo de entradas
input1 = input('Ingresa el nombre del archivo que contiene los valores de entrada "*.txt" (sin extension): ','s');
nombreArchivo1 = strcat(input1,'.txt');
p = importdata(nombreArchivo1);

% Leemos el archivo de los valores target
input2 = input('Ingresa el nombre del archivo que contiene los valores target "*.txt" (sin extension): ','s');
nombreArchivo2 = strcat(input2,'.txt');
targets = importdata(nombreArchivo2);
fprintf('\n');

% Se calcula el rango a trabajar
dim_p = size(p);
lim_inf = p(1);
lim_sup = p(dim_p(1));
incremento = abs(p(1)-p(2));
if lim_inf < 0
    if lim_sup <0
        rango = lim_inf:incremento:lim_sup;
    else
        rango = lim_inf:incremento:lim_sup;
    end
else
    rango = lim_inf:incremento:lim_sup;
end

% Se muestra el rango y el numero de datos a trabajar
fprintf('Se trabajara el siguiente rango, de acuerdo al archivo:\n');
fprintf('[%f, %f], con un incremento de %f\n',lim_inf,lim_sup,incremento);
num_datos = dim_p(1);
fprintf('Se trabajara con %d datos, de acuerdo al archivo.\n\n',num_datos);

% Se solicita la arquitectura del M.L.P. 
R = input('Ingresa el tamanio del vector de entrada: ');
num_capas = input('Ingresa el numero de capas del M.L.P.: ');
arq_mlp = zeros(1,num_capas+1);
fun_capa = zeros(1,num_capas);
arq_mlp(1) = R;
fprintf('\n');
fprintf('Se solicitara la arquitectura del M.L.P.\n');
fprintf('Para las funciones de activacion se tienen las siguientes:\n');
fprintf('1) purelin()\n2) logsig()\n3) tansig()\n');
for i=2:num_capas+1
    fprintf('Ingresa el numero de neuronas de la capa %d: ',i-1);
    arq_mlp(i) = input('');
    fprintf('Ingresa la funcion de activacion de la capa %d: ',i-1);
    fun_capa(i-1) = input(''); 
end
disp('La arquitectura del M.L.P. es:');
disp(arq_mlp);
disp(fun_capa);

% Se solicita el valor del factor de aprendizaje
alfa = input('Ingresa el valor del factor de aprendizaje(alfa): ');

% Se solicitan los valores de usuario eit,itmax, itval, numval
itmax = input('Ingresa el numero de iteraciones maximas de la red(itmax): ');
itval = input('¿Cada cuanto se hara una iteracion de validacion? (itval): ');
numval = input('Numero maximo de incrementos consecutivos del error de validacion (numval): ');
eit = input('Ingrese l valor minimo del error en una epoca (eit): ');
fprintf('\n');

% Se solicita la configuracion para dividir en subconjuntos
fprintf('Seleccione una de las siguientes configuraciones a trabajar:\n');
fprintf('1) 80-10-10\n');
fprintf('2) 70-15-15\n');
config = input('Ingrese su seleccion: ');

% Se forman los subconjutnos de acuerdo al tipo de configuracion
switch(config)
    case 1
        num_elem_val = round(num_datos*.2); % Numero de elementos del conjunto de validacion
        num_elem_prueba = num_elem_val; % Numero de elementos del conjunto de prueba
        num_elem_ent = num_datos - 2*num_elem_val; % Numero de elementos del conjunto de entrenamiento
    case 2
        num_elem_val = round(num_datos*.15);
        num_elem_prueba = num_elem_val;
        num_elem_ent = num_datos - 2*num_elem_val;
end
disp('Se usaran los siguientes tamanios en los subconjuntos:');
fprintf('Conjunto de entrenamiento: %d elementos\n',num_elem_ent);
fprintf('Conjunto de validacion: %d elementos\n',num_elem_val);
fprintf('Conjunto de prueba: %d elementos\n',num_elem_prueba);
cto_val = obtenerConjuntoDeValidacion(p,targets,num_datos,num_elem_val);
cto_prueba = obtenerConjuntoDePrueba(p,targets,num_datos,num_elem_prueba);
cto_ent = obtenerConjuntoDeEntrenamiento(p,targets,num_datos,num_elem_ent,cto_val,cto_prueba);
disp('Conjunto de validacion:');
disp(cto_val);
disp('Conjunto de prueba:');
disp(cto_prueba);
disp('Conjunto de entrenamiento:');
disp(cto_ent);

% Se inicializan la matriz de pesos y el vector bias con valores aleatorios
% entre -1 y 1
W = cell(num_capas,1);
b = cell(num_capas,1);
disp('Valores iniciales de las matrices:');
for i=1:num_capas
    temp_W = -1 + (2)*rand(arq_mlp(i+1),arq_mlp(i));
    W{i} = temp_W;
    fprintf('W_%d = \n',i);
    disp(W{i});
    temp_b = -1 + (2)*rand(arq_mlp(i+1),1);
    b{i} = temp_b;
    fprintf('b_%d = \n',i);
    disp(b{i});
end

% Se utiliza una cell para guardar las salidas de cada capa
a = cell(num_capas+1,1);

% Se utiliza una cell para guardas las sensitividades de cada capa
S = cell(num_capas,1);
F_m = cell(num_capas,1);
X = input('Presiona cualquier tecla para comenzar el aprendizaje...');

% Comienza el aprendizaje
Err_val = 0;
count_val = 0;
for it=1:itmax
    Eap = 0; % Error de aprendizaje
    
    % Si no es una iteracion de validacion
    if(mod(it,itval)~=0)
        for dato=1:num_elem_ent
            a{1} = cto_ent(dato,1); % Condicion inicial
            
            % Se propaga hacia adelante el elemento del cto. de
            % entrenamiento
            for k=1:num_capas
                W_temp = cell2mat(W(k));
                a_temp = cell2mat(a(k));
                b_temp = cell2mat(b(k));
                a{k+1} = funcionDeActivacion(W_temp*a_temp+b_temp,fun_capa(k));
            end
            a_temp = cell2mat(a(num_capas+1));
            ej = cto_ent(dato,2)-a_temp;
            Eap = Eap+(ej/num_datos);
            
            % Se calculan las sensitividades y se propagan hacia atras,
            % es decir, inicia el backpropagation.
            F_m{num_capas} = obtenerF(fun_capa(num_capas),arq_mlp(num_capas+1),a_temp);
            F_m_temp = cell2mat(F_m(num_capas));
            S{num_capas} = -2*F_m_temp*(ej);
            for m = num_capas-1:-1:1
                W_temp = cell2mat(W(m+1));
                a_temp = cell2mat(a(m+1));
                S_temp = cell2mat(S(m+1));
                F_m{m} = obtenerF(fun_capa(m),arq_mlp(m+1),a_temp);
                F_m_temp = cell2mat(F_m(m));
                S{m} = F_m_temp*(W_temp')*S_temp;
            end
            
            % Se aplican las reglas de aprendizaje
            for k=num_capas:-1:1
                W_temp = cell2mat(W(k));
                b_temp = cell2mat(b(k));
                a_temp = cell2mat(a(k));
                S_temp = cell2mat(S(k));
                W{k} = W_temp-alfa*S_temp*(a_temp');
                b{k} = b_temp-alfa*S_temp;
            end
            
        end
        
    % Si es una iteracion de validacion    
    else
        E_val = 0;
        for dato=1:num_elem_val
            a{1} = cto_val(dato,1); % Condicion inicial
            % Se propaga hacia adelante el elemento del cto. de
            % validacion.
            for k=1:num_capas
                W_temp = cell2mat(W(k));
                a_temp = cell2mat(a(k));
                b_temp = cell2mat(b(k));
                a{k+1} = funcionDeActivacion(W_temp*a_temp+b_temp,fun_capa(k));
            end
            a_temp = cell2mat(a(num_capas+1));
            e_val = cto_val(dato,2)-a_temp;
            E_val = E_val+(e_val/num_elem_val);
        end
        if count_val == 0
            Err_val = E_val;
            count_val = count_val+1;
            fprintf('Count val = %d\n',count_val);
        else
            if count_val == numval
                fprintf('Early stopping en iteracion %d\n',it);
                break;
            else
                if E_val > Err_val
                    count_val = count_val+1;
                    fprintf('Count val = %d\n',count_val);
                else
                    Err_val = 0;
                    count_val = 0;
                    fprintf('Count val = %d\n',count_val);
                end
            end
        end
    end
    
    % Se comprueban las condiciones de finalizacion
    if Eap < eit && Eap > 0
        fprintf('Aprendizaje exitoso en la iteracion %d\n',it);
        break;
    end
end

% Se propaga el conjunto de prueba
for i=1:num_elem_prueba
    a{1} = cto_prueba(i,1); % Condicion inicial
    % Se propaga hacia adelante el elemento del cto. de
    % entrenamiento
    for k=2:num_capas+1
        W_temp = cell2mat(W(k-1));
        a_temp = cell2mat(a(k-1));
        b_temp = cell2mat(b(k-1));
        a{k} = funcionDeActivacion(W_temp*a_temp+b_temp,fun_capa(k-1));
    end
    dato_entrada = cell2mat(a(1));
    a_temp = cell2mat(a(num_capas+1));
    fprintf('Para %f la salida es %f\n',dato_entrada,a_temp);
end
