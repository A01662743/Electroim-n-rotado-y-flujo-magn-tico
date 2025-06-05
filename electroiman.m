clear;
clc;   
close all;

%constantes
mu0 = 4 * pi * 1e-7;
I = 10;
L_alambre = 3;

%alambres
ai1 = [0, 0, -L_alambre/2];
af1 = [0, 0, L_alambre/2];

ai2 = af1;
af2 = [L_alambre, 0, L_alambre/2];

ai3 = af2;
af3 = [L_alambre, 0, -L_alambre/2];

ai4 = af3;
af4 = ai1;


ai = [ai1; ai2; ai3; ai4];
af = [af1; af2; af3; af4];
acent = [L_alambre/2, 0, 0];

% rotación

ang1 = pi/4;
ang2 = pi/4;
ang3 = pi/4;

Rx = [1, 0, 0;
       0, cos(ang1), -sin(ang1);
       0, sin(ang1), cos(ang1)];
Rz = [cos(ang2), -sin(ang2), 0;
       sin(ang2), cos(ang2), 0;
       0, 0, 1];
Ry = [cos(ang3), 0, sin(ang3);
       0, 1, 0;
       -sin(ang3), 0, cos(ang3)];

for n = 1:length(ai)
    ai(n,:) = (Ry * (Rz * (Rx*(ai(n,:) - acent).'))).' + acent;
    af(n,:) = (Ry * (Rz * (Rx*(af(n,:) - acent).'))).' + acent;
end

a = ai -af;

%area

circent = [1.5, 1.8, 0];
rad = 1.5;

crange = linspace(0, 2*pi, 100);

circx = circent(1) + rad .* cos(crange);
circz = circent(3) + rad .* sin(crange);
circy = repmat(circent(2), size(crange));

% meshgrid
rango = 5;
N = 20;

x = linspace(-rango, L_alambre + rango, N); % Para incluir el alambre en X
y = linspace(-rango, rango, N); % Para incluir y alrededor de Y=0
z = linspace(-(L_alambre/2 + rango), L_alambre/2 + rango, N); % Para incluir Z alrededor de los alambres
[X, Y, Z] = meshgrid(x, y, z);

% Inicializamos matrices para campo magnético
Bx = zeros(N,N,N);
By = zeros(N,N,N);
Bz = zeros(N,N,N);

%cálculo campo por cada punto del meshgrid
for i = 1:N
    for j = 1:N
        for k = 1:N

            %puntos importantes
            P = [X(i,j,k) Y(i,j,k) Z(i,j,k)];
            B = zeros(length(a), 3);
            for m = 1:length(a)
                b = af(m,:) - P;
                c = ai(m,:) - P;
                
                cxa = cross(c, a(m,:));
                adc = dot(a(m,:), c);
                adb = dot(a(m,:), b);
                cxa2 = norm(cxa)^2;
                
                %campo [x, y, z]
                B(m,:) = ((mu0*I*cxa)/(4*pi*abs(cxa2))) .* ((adc./abs(c)) - (adb./abs(b)));
                
                %evitar errores
                B(isnan(B)) = 0;
                B(isinf(B)) = 0;
            end
            %separar campo por eje
            Bx(i,j,k) = sum(B(:,1));
            By(i,j,k) = sum(B(:,2));
            Bz(i,j,k) = sum(B(:,3));
        end
    end
end

%meshgrid círculo (area para calcular el flujo)
num_circ = 30; 
r_grid = linspace(0, rad, num_circ);
theta_grid = linspace(0, 2*pi, num_circ);
[R_circ, Theta_circ] = meshgrid(r_grid, theta_grid);

% coordenadas polares a cartesianas
X_circ = (R_circ .* cos(Theta_circ)) + circent(1);
Z_circ = (R_circ .* sin(Theta_circ)) + circent(3);
Y_circ = circent(2) * ones(size(X_circ));

%modificar Bx, By y Bz para usar en interp
Bx_perm = permute(Bx, [2 1 3]);
By_perm = permute(By, [2 1 3]);
Bz_perm = permute(Bz, [2 1 3]);

%usar interp para calcular campo en círculo
Bx_circ = interpn(x, y, z, Bx_perm, X_circ, Y_circ, Z_circ, 'linear', 0);
By_circ = interpn(x, y, z, By_perm, X_circ, Y_circ, Z_circ, 'linear', 0);
Bz_circ = interpn(x, y, z, Bz_perm, X_circ, Y_circ, Z_circ, 'linear', 0);

%obtener dr y dtheta
dr = r_grid(2) - r_grid(1);
dtheta = theta_grid(2) - theta_grid(1);

%calcular dA
dA = R_circ .* dr .* dtheta;

% obtener flujo neto
flujo = sum(By_circ(:) .* dA(:));

%imprimir flujo
fprintf('Magnetic Flux through the circle: %.4e Weber\n', ...
        flujo);

% Graficar
figure;
%campo total
quiver3(X, Y, Z, Bx, By, Bz, 2);
%alambres
for n = 1:length(a)
    hold on
    plot3(ai(n,1), ai(n,2), ai(n,3), 'k-o', 'LineWidth', 2, 'MarkerSize', 4); 
    plot3(af(n,1), af(n,2), af(n,3), 'k-o', 'LineWidth', 2, 'MarkerSize', 4);
    line([ai(n, 1) af(n, 1)], [ai(n, 2) af(n, 2)], [ai(n, 3) af(n, 3)], 'Color', 'k', 'LineWidth', 3); %mostrar alambre
end

%círculo
hold on
plot3(circx, circy, circz, 'r--', 'LineWidth', 2, 'DisplayName', 'Círculo en plano XY');
plot3(circent(1), circent(2), circent(3), 'k.');

% nombres
hold on
xlabel('X (m)');
ylabel('Y (m)');
title('Lineas de Campo Magnetico');
axis equal;
grid on;


% graficar flujo
figure('Name', 'Campo Magnético sobre la Superficie del Círculo');

% Graficar área del círculo
surf(X_circ, Y_circ, Z_circ, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', 'Superficie del Círculo');
hold on;

% Graficar campo en el círculo
quiver3(X_circ, Y_circ, Z_circ, Bx_circ, By_circ, Bz_circ, 1, 'Color', 'b', 'LineWidth', 1);

%graficar circunferencia y centro
plot3(circx, circy, circz, 'r--', 'LineWidth', 2, 'DisplayName', 'Contorno del Círculo');
plot3(circent(1), circent(2), circent(3), 'k.', 'MarkerSize', 15, 'DisplayName', 'Centro del Círculo');
