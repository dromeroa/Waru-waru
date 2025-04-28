% Este programa calcula la temperatura del waru-waru de la planta

clear;clc;

function Tc = calculaTemperaturaCultivo(Q,raa,raw,rac,Ra,Ta,Tw)
    sigma = 5.670374419*10^(-8);
    % Esta función calcula la temperatura del cultivo
    Tc=(Q*((Ta/raa)+(Tw/raw)) + ((Ra + 3*sigma*(Ta.^4))*((1/raa) + (1/raw) + (1/rac) )))/( Q*((1/raa) + (1/raw)) + 4*sigma*(Ta.^3)*((1/raa) + (1/raw) + (1/rac) )); 
end


% Definimos las constantes
cd = 0.2; % mean dragg coefficient (Llhomme)
L0 = 4; % indice del area de la hoja (Definir con Jael) (asumimos por ahora 4) 
zh = 1; % altura de la cubierta vegetal (m) 
X = cd*L0;
Wc = 1; % ancho de la plataforma (m) 
Ww = 1; % ancho del canal (m) (asumimos por ahora 1m) Con 2m probar
zr = 2.0; % altura de referencia sobre el nivel del suelo (2m según paper Llhomme 2002)
z0c = 0.01; % Llhomme
z0w = 0.0002; % Llhomme
k = 0.4; % k : Constante de Von Karman (0.4 según Llhomme) 
ua = 2.16; % velocidad del viento a la altura de referencia (m/s) (segun datos de Henry)
alphaW = 2.5; % constante sin dimensión = 2.5 (según Llhomme)
kzh = 0.5; % Coeficiente de distorsión turbulenta a la altura de la cubierta vegetal
w1 = 0.10; % ancho de la hoja (m)  ( escogemos 0.1 m,  segun el libro la quinua y la kañiwa, Mario Tapia mide hasta 0.12 m de ancho por 0.15 de alto)
uzh = 2.0; % rapidez del viento a la altura de la planta (m/s)  (asumimos 2m/s este dato equivale a 7.2 km/h)
alpha0 = 0.005; % constante = 0.005 (ms^-1/2)
rho = 0.832; % densidad del aire (kg/m3) (El valor es 0.832) 
CP = 1004; %(Si ponemos 30000 se ajusta ver las unidades, el valor referido es 1004)  calor específico del aire a presión constante (J/kg K)
c = 0.5; % razón del flujo del calor del suelo a la radiación neta del cultivo (0,50)
epsic = 0.97; % emisividad de la vegetacion (El valor es 0.97) Con 0.1 queda bien
Ra = 150; %Radiación incidente de onda larga (W/m2) 


% Calculamos z0s
z0s = ((Wc*z0c) + (Ww*z0w))/(Wc + Ww);

%Calculamos z0
z0 = z0s + 0.3*zh*(X^(0.5));

%Calculamos d
d = (1.1)*zh*log(1 + X^(0.25) );

%Calaculamos raa
raa =  (1/(k*k*ua))*(log((zr-d)/z0))^(2);

%Calaculamos raw
raw = (zh/(alphaW*kzh))*exp(alphaW)* (exp(-alphaW*z0s/zh) - exp(-alphaW*(d+z0)/zh));

%Calaculamos rac
rac = (alphaW*((w1/uzh)^(0.5)))/(4*L0*alpha0*(1-exp(-alphaW/2)));

%Calaculamos Q
Q = (rho*CP*(Wc + Ww))/(rac*(1-c)*Wc*epsic);




% Llamamos al archivo csv con los datos experimentales
datos = readtable('Waru_replica_del_19al26_enero_2025_modificado', 'VariableNamingRule', 'preserve');


% Suponiendo que ya tienes una tabla 'datos' y quieres llenar la columna 'resultado'
datos.resultado = NaN(height(datos), 1);  % Inicializamos la columna con NaN o ceros

for i = 1:height(datos)  % Recorremos todos los índices de la tabla
    % Realizamos el cálculo con el valor de la columna S3 para el índice i
    datos.resultado(i) = calculaTemperaturaCultivo(Q,raa,raw,rac,Ra,datos.S6(i),datos.S3(i));  
end


% Creamos una tabla con los datos experimentales y con los datos del modelo
nuevaTabla = table(datos.S5, datos.S3, datos.S6, datos.resultado);

Graf1 = plot(datos.Tiempo,datos.resultado,'b','LineWidth', 1);
title('Gráfico combinado');
xlabel('Tiempo');
ylabel('Temperatura (°C)');
hold on;
Graf4 = plot(datos.Tiempo,datos.S5,'m','LineWidth', 1);

grid on;

legend('Temperatura del cultivo teorica', 'Temperatura del cultivo experimental','Location', 'best','FontSize', 14);
   % Modificar los límites de los ejes
ylim([0 40]); % Limitar el eje Y entre -1 y 1