function [I_cell] = NRSDM(I_ph, I_o, A, R_s, R_sh, V_T, V_cell)
    % Se define la ecuación.
    syms I_cell
    
    % SDM
    eq = I_ph - I_o*(exp((V_cell + (I_cell*R_s))/(A*V_T)) - 1) - (V_cell + (I_cell*R_s))/R_sh - I_cell;

    % Se define la derivada de la ecuación.
    deq = diff(eq);

    I_cell_0 = 0;  % Valor inicial de I_cell.
    max_iter = 100; % Número máximo de iteraciones.
    Tol = 1e-6; % Tolerancia.

    % Se inicializan las variables.
    I_cell_prev = I_cell_0;
    iter = 0;

    while iter < max_iter
        iter = iter + 1; % Aumenta en 1 el número de iteraciones.
    
        % Evaluar la función y su derivada en el punto actual.
        f = double(subs(eq, I_cell, I_cell_prev));
        df = double(subs(deq, I_cell, I_cell_prev));
    
        I_cell_next = I_cell_prev - f/df; % Calcula el siguiente valor de I_cell.
    
        % Verificar convergencia.
        if abs(I_cell_next - I_cell_prev) < Tol
            break;
        end

        I_cell_prev = I_cell_next; % Actualiza el valor de I_cell para cada iteración.
    end

    I_cell = I_cell_next; % Devuelve el valor de I_cell estimado.

end