%% initialize the system
numerator=[11.88 4.977 539.6 129.9 5625]
denominator=[1 1.169 50.3 45.94 685.3 391.7 1952]
sys=tf(numerator, denominator);

%% Bode Plot of the system
bode(sys); % change frequency unit to KHz in "View" -> "Property Editor"

%% Step response of the system
plot(step(sys));

%% The relative degree of the system
% relative degree =  6-4 =2

%%
csys= canon(sys,'companion')