% Algoritmo genético. BISHOP MODEL.
clc; clear; close all

%% Cargar los datos de corriente y voltaje medidos de la celda.
Vcell = [-13.675, -13.674, -13.673, -13.586, -13.584, -13.583, -13.582, -13.5, -13.498, -13.496, -13.494, -13.486, -13.409, -13.407, -13.406, -13.404, -13.318, -13.316,  -13.227, -13.225, -13.134, -13.042, -13.041, -12.854, -12.761, -12.574, -12.572, -12.48, -12.386, -12.292, -12.198, -12.103, -11.820, -11.631, -11.442, -11.253, -10.969, -10.59, -10.117, -9.549, -8.791, -8.129, -7.371, -6.424, -5.288, -3.678, -2.542, -1.027, 0.109, 0.482, 0.487, 0.5452, 0.5519,  0.558, 0.565, 0.569, 0.572]; % Vector de voltaje medida de la celda. (Cambiar si necesario)
Icellm = [9.949, 9.874, 9.748, 9.673, 9.547, 9.447, 9.296, 9.145, 8.844, 8.643, 8.417, 8.291, 8.19, 8.09, 7.989, 7.788, 7.613, 7.386, 7.185, 7.01, 6.859, 6.658, 6.507, 6.331, 6.13, 6.055, 5.879, 5.804, 5.728, 5.678, 5.603, 5.552, 5.477, 5.402, 5.326, 5.276, 5.226, 5.175, 5.1, 5.025, 4.949, 4.887, 4.798, 4.673, 4.572, 4.396, 4.296, 4.12, 4.045, 3.864, 3.417, 2.939, 2.562, 2.16, 1.457, 0.728, 0.251]; % Vector de corriente medida de la celda. (Cambiar si necesario)
 
%% Rangos de búsqueda de los parámetros.
Isc = 4; % Corriente en corto circuito de la celda (A).

IphMin = 0.90*Isc; IphMax = 1.1*Isc; % Rango parámetro 1. Photocorriente.
IoMin = 1e-10; IoMax = 1e-7; % Rango parámetro 2. Corriente de saturación inversa.
AMin = 0.05; AMax = 4; % Rango parámetro 3. Factor de idealidad del diodo.
RsMin = 1e-5; RsMax = 2; % Rango parámetro 4. Resistencia en serie.
RshMin = 20; RshMax = 100; % Rango parámetro 5. Resistencia shunt.
aMin = 1e-3; aMax = 30e-3; % Rango parámetro 6. Fracción ohmica de la corriente relacionada a la ruptura de avalancha.
mMin = 2; mMax = 8; % Rango parámetro 7. Exponente de la ruptura de avalancha.
VbrMin = -10; VbrMax = -50; % Rango parámetro 8. Voltaje de ruptura de unión.

%% Inicializa aleatoriamente los parámetros.
r = 8; % Número de parámetros. (Cambiar si necesario)
p = 14;  % Número de individuos. (Cambiar si necesario)

Iph = IphMin + (IphMax - IphMin)*rand(1,p); % Vector aleatorio del parámetro 1.
Io = IoMin + (IoMax - IoMin)*rand(1,p); % Vector aleatorio del parámetro 2.
A = AMin + (AMax - AMin)*rand(1,p); % Vector aleatorio del parámetro 3.
Rs = RsMin + (RsMax - RsMin)*rand(1,p); % Vector aleatorio del parámetro 4.
Rsh = RshMin + (RshMax - RshMin)*rand(1,p); % Vector aleatorio del parámetro 5.
a = aMin + (aMax - aMin)*rand(1,p); % Vector aleatorio del parámetro 6.
m = mMin + (mMax - mMin)*rand(1,p); % Vector aleatorio del parámetro 7.
Vbr = VbrMin + (VbrMax - VbrMin)*rand(1,p); % Vector aleatorio del parámetro 8.

%% Generar población inicial.
iter = 0; % Inicializa las iteraciones.
iter_n = 0; % Inicializa las iteraciones consecutivas del titular.

Ind = zeros(1,r); % Inicializa el individuo.
Pop = zeros(p,r); % Inicializa la población.
for i = 1:p
    Ind = [Iph(i), Io(i), A(i), Rs(i), Rsh(i), a(i), m(i), Vbr(i)]; % El individuo.
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
        Icelle(i,j) = NRBishop(Iph(i), Io(i), A(i), Rs(i), Rsh(i), a(i), m(i), Vbr(i), VT, Vcell(j)); % Calcula la corriente utilizando el método de Newton Raphson.
    end
end

% Raíz del error cuadrático medio.
OF = rmse(Icelle',Icellm')'; % Se evalua la función objetivo con la corriente estimada y la corriente medida.

%% Seleccionar al titular.
% El titular es el individuo con el menor valor de la función objetivo.
[MinOF, nc] = min(OF); % Encuentra el valor mínimo, y su posición, de la función objetivo evaluada. (El mejor individuo)
Inc = Pop(nc,:); % El titular.

%% Geración descendiente.
max_iter = 709; % Número máximo de iteraciones.
max_iter_n = 161; % Número máximo de iteraciones consecutivas del titular.

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
            case 6 % Se muta el sexto parámetro.
                Child1(u_i2) = aMin + (aMax - aMin)*rand();
                Child2(u_i2) = aMin + (aMax - aMin)*rand();
            case 7 % Se muta el séptimo parámetro.
                Child1(u_i2) = mMin + (mMax - mMin)*rand();
                Child2(u_i2) = mMin + (mMax - mMin)*rand();
            case 8 % Se muta el octavo parámetro.
                Child1(u_i2) = VbrMin + (VbrMax - VbrMin)*rand();
                Child2(u_i2) = VbrMin + (VbrMax - VbrMin)*rand();
        end
    end

    % Seleccionar a un solo hijo de manera aleatoria.
    Children = [Child1; Child2]; % Matriz de hijos.
    u_i3 = randi([1 2]); % Seleccionar un número aleatorio entre 1 y 2.
    Child = Children(u_i3,:); % Hijo seleccionado aleatoriamente.

    % Evaluar la función objetivo del hijo seleccionado.
    Icell_Child = zeros(1,size(Vcell,2)); % Inicializa el vector de corriente estimada del hijo.
    for i = 1:size(Vcell,2)
        Icell_Child(1,i) = NRBishop(Child(1), Child(2), Child(3), Child(4), Child(5), Child(6), Child(7), Child(8), VT, Vcell(i)); % Calcula la corriente del hijo utilizando el método de Newton Raphson.
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
disp(['a = ', num2str(Inc(6))])
disp(['m = ', num2str(Inc(7))])
disp(['Vbr = ', num2str(Inc(8))])

%% Grafica de corriente contra voltaje.
Icell = zeros(1,size(Vcell,2)); % Inicializa el vector de corriente de la celda.
for i = 1:size(Vcell,2)
    Icell(1,i) = NRBishop(Inc(1), Inc(2), Inc(3), Inc(4), Inc(5), Inc(6), Inc(7), Inc(8), VT, Vcell(i)); % Calcula la corriente utilizando el método de Newton Raphson.
end

figure(1)
plot(Vcell,Icellm) % Gráfica real de la corriente en función del voltaje.
hold on
plot(Vcell,Icell) % Grafica estimada de la corriente en función del voltaje.
legend('Curva real', 'Curva estimada')
ylabel('Current (A)')
xlabel('Voltage (V)')