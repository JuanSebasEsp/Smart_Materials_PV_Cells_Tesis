% Algoritmo genético. SINGLE DIODE MODEL.
clc; clear; close all

%% Cargar los datos de corriente y voltaje medidos de la celda.
Vcell = 0:0.1:0.5; % Vector de voltaje medida de la celda. (Cambiar si necesario)
Icellm = [0.4266, 0.42, 0.4130, 0.3988, 0.2749, 0]; % Vector de corriente medida de la celda. (Cambiar si necesario)

%% Rangos de búsqueda de los parámetros.
Isc = 0.43; % Corriente en corto circuito de la celda (A).

IphMin = 0.90*Isc; IphMax = 1.1*Isc; % Rango parámetro 1. Photocorriente.
IoMin = 1e-10; IoMax = 1e-7; % Rango parámetro 2. Corriente de saturación inversa.
AMin = 0.05; AMax = 4; % Rango parámetro 3. Factor de idealidad del diodo.
RsMin = 1e-5; RsMax = 2; % Rango parámetro 4. Resistencia en serie.
RshMin = 20; RshMax = 100; % Rango parámetro 5. Resistencia shunt.

%% Inicializa aleatoriamente los parámetros.
r = 5; % Número de parámetros. (Cambiar si necesario)
p = 11;  % Número de individuos. (Cambiar si necesario)

Iph = IphMin + (IphMax - IphMin)*rand(1,p); % Vector aleatorio del parámetro 1.
Io = IoMin + (IoMax - IoMin)*rand(1,p); % Vector aleatorio del parámetro 2.
A = AMin + (AMax - AMin)*rand(1,p); % Vector aleatorio del parámetro 3.
Rs = RsMin + (RsMax - RsMin)*rand(1,p); % Vector aleatorio del parámetro 4.
Rsh = RshMin + (RshMax - RshMin)*rand(1,p); % Vector aleatorio del parámetro 5.

%% Generar población inicial.
iter = 0; % Inicializa las iteraciones.
iter_n = 0; % Inicializa las iteraciones consecutivas del titular.

Ind = zeros(1,r); % Inicializa el individuo.
Pop = zeros(p,r); % Inicializa la población.
for i = 1:p
    Ind = [Iph(i), Io(i), A(i), Rs(i), Rsh(i)]; % El individuo.
    Pop(i,:) = Ind; % La población.
end

%% Evaluar la función objetivo.
% Se quiere minimizar el RMSE entre la corriente medida y el valor estimado.
q = 1.602e-19; % Carga elemental del electrón (C).
k = 1.38e-23; % Constante de Boltzman (J/K).
T = 25 + 273.15; % Temperatura de la celda (K).

VT = (k*T)/q; % Voltaje térmico.

% Calculo del valor estimado de la corriente de la celda.
Icelle = zeros(p,size(Vcell,2)); % Inicializa el vector de corriente estimada.
for i = 1:p
    for j = 1:size(Vcell,2)
        Icelle(i,j) = NRSDM(Iph(i), Io(i), A(i), Rs(i), Rsh(i), VT, Vcell(j)); % Calcula la corriente utilizando el método de Newton Raphson.
    end
end

% Raíz del error cuadrático medio.
OF = rmse(Icelle',Icellm')'; % Se evalua la función objetivo con la corriente estimada y la corriente medida.

%% Seleccionar al titular.
% El titular es el individuo con el menor valor de la función objetivo.
[MinOF, nc] = min(OF); % Encuentra el valor mínimo, y su posición, de la función objetivo evaluada. (El mejor individuo)
Inc = Pop(nc,:); % El titular.

%% Geración descendiente.
max_iter = 876; % Número máximo de iteraciones.
max_iter_n = 206; % Número máximo de iteraciones consecutivas del titular.

while iter < max_iter && iter_n < max_iter_n
    iter = iter + 1; % Aumenta en 1 el número de iteraciones.
    iter_n = iter_n + 1; % Aumenta en 1 el número de iteraciones consecutivas del titular.

    % Seleccionar padres de la población inicial.
    u_p = randperm(p,2); % Dos enteros únicos seleccionados aleatoriamente desde 1 hasta el número de individuos.
    Parent1 = Pop(u_p(1),:); % Primer padre.
    Parent2 = Pop(u_p(2),:); % Segundo padre.

    % Crear hijos al recombinar padres.
    u_i = randi([2 r]); % Punto de recombinación de parámetros.
    Child1 = [Parent1(1:u_i-1), Parent2(u_i:end)]; % Primer hijo.
    Child2 = [Parent2(1:u_i-1), Parent1(u_i:end)]; % Segundo hijo.

    % Mutación.
    mut = randi([0 1]); % Generar un 0 o 1 aleatoriamente.

    if mut == 1 % En este caso hay mutación de los hijos.
        % Seleccionar parámetro a mutar en cada hijo.
        u_i2 = randi(r); % Punto de mutación de parámetro.

        % Seleccionar un valor aleatorio del parámetro dentro del rango del parámetro.
        switch u_i2
            case 1 % Se muta el primer parámetro.
                Child1(u_i2) = IphMin + (IphMax - IphMin)*rand();
                Child2(u_i2) = IphMin + (IphMax - IphMin)*rand();
            case 2 % Se muta el segundo parámetro.
                Child1(u_i2) = IoMin + (IoMax - IoMin)*rand();
                Child2(u_i2) = IoMin + (IoMax - IoMin)*rand();
            case 3 % Se muta el tercer parámetro.
                Child1(u_i2) = AMin + (AMax - AMin)*rand();
                Child2(u_i2) = AMin + (AMax - AMin)*rand();
            case 4 % Se muta el cuarto parámetro.
                Child1(u_i2) = RsMin + (RsMax - RsMin)*rand();
                Child2(u_i2) = RsMin + (RsMax - RsMin)*rand();
            case 5 % Se muta el quinto parámetro.
                Child1(u_i2) = RshMin + (RshMax - RshMin)*rand();
                Child2(u_i2) = RshMin + (RshMax - RshMin)*rand();
        end
    end

    % Seleccionar a un solo hijo de manera aleatoria.
    Children = [Child1; Child2]; % Matriz de hijos.
    u_i3 = randi([1 2]); % Seleccionar un número aleatorio entre 1 y 2.
    Child = Children(u_i3,:); % Hijo seleccionado aleatoriamente.

    % Evaluar la función objetivo del hijo seleccionado.
    Icell_Child = zeros(1,size(Vcell,2)); % Inicializa el vector de corriente estimada del hijo.
    for i = 1:size(Vcell,2)
        Icell_Child(1,i) = NRSDM(Child(1), Child(2), Child(3), Child(4), Child(5), VT, Vcell(i)); % Calcula la corriente del hijo utilizando el método de Newton Raphson.
    end

    OF_Child = rmse(Icell_Child,Icellm); % Se evalua la función objetivo con la corriente estimada del hijo y la corriente medida.
    
    % Actualizar el titular.
    % Criterio de diversificación.
    Check = ~any(ismember(Child, Pop, "rows")); % Verifica que el hijo seleccionado no sea parte de la población.

    if any(Check) % En este caso el hijo seleccionado no hace parte de la población.
        [MaxOF, nd] = max(OF); % Encuentra el valor máximo, y su posición, de la función objetivo evaluada. (El peor individuo)
        Pop(nd,:) = Child; % Reemplaza al peor individuo por el hijo seleccionado en la población.
        OF(nd) = OF_Child; % Reemplaza el peor valor de la función objetivo de la población por el valor de la función objetivo del hijo seleccionado.

        % Verificar si este nuevo individuo tiene el menor valor de la función objetivo.
        if OF_Child < MinOF % En este caso el nuevo individuo es el mejor individuo.
            Inc = Child; % El nuevo individuo se vuelve el nuevo titular.
            [MinOF, nc] = min(OF); % Encuentra el valor mínimo, y su posición, de la función objetivo evaluada. (El mejor individuo)
            iter_n = 0; % Reinicia el conteo de iteraciones consecutivas del titular.
        end
    end
end

%% Imprimir resultados.
disp(['El individuo titular que contiene los parámetros deseados es el individuo ', num2str(nc)])
disp('Los parámetros son: ')
disp(['Iph = ', num2str(Inc(1))])
disp(['Io = ', num2str(Inc(2))])
disp(['A = ', num2str(Inc(3))])
disp(['Rs = ', num2str(Inc(4))])
disp(['Rsh = ', num2str(Inc(5))])

%% Grafica de corriente contra voltaje.
Icell = zeros(1,size(Vcell,2)); % Inicializa el vector de corriente de la celda.
for i = 1:size(Vcell,2)
    Icell(1,i) = NRSDM(Inc(1), Inc(2), Inc(3), Inc(4), Inc(5), VT, Vcell(i)); % Calcula la corriente utilizando el método de Newton Raphson.
end

figure(1)
plot(Vcell,Icellm) % Gráfica real de la corriente en función del voltaje.
hold on
plot(Vcell,Icell) % Grafica estimada de la corriente en función del voltaje.
legend('Curva real', 'Curva estimada')
ylabel('Current (A)')
xlabel('Voltage (V)')