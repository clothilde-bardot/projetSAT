function [code_L5I, code_L5Q] = GPS_L5_code_prim(prn)

% Error messages

if nargin~=1
   error('insufficient number of input arguments')
end
if prn<0 || prn > 37
   error('invalid prn: must be between 1 and 37')
end

%XA Coder = cryptage ? ressemble a de la LSFR (linear Feedback Shift
%Register), brouillage ou modulation => mod BPSK
code_XA         =   zeros(1,10230);  %tableau de 0 pour stocker la séquence générée
Init_state_XA   =   [1 1 1 1 1 1 1 1 1 1 1 1 1]; %état initial 
Poly_XA         =   [0 0 0 0 0 0 0 0 1 1 0 1 1]; %polynome générateur/polynôme de rétroaction

State_XA = Init_state_XA; %état initial, state_XA évolue à chaque itération
%for i = 1:6190
for i = 1:8189
    root = Poly_XA.*State_XA; %multiplie par le polynome pour identifier les bits à prendre en compte
    reg = xor(xor(xor(root(9),root(10)),root(12)),root(13));% on applique XOR sur les bits utilis (9,10,12,13  cf polynôme)
    code_XA(i) = State_XA(13);%stocke le dernier bit de l'état dans le code généré
    for j=13:-1:2
        State_XA(j) = State_XA(j-1); %décale tous les bits de l'état vers la droite
    end
    State_XA(1) = reg; %Insère le nouveau bit calculé au début 

end

%code_XA(6191:10230) = code_XA(1:10230-6191+1);
code_XA(8191:10230) = code_XA(1:10230-8191+1); %copie les 1890 bits de code_XA dans les indices 8191 à 10230

 
%XBI Coder
code_XBI = zeros(1,10230);
Poly_XBI = [1 0 1 1 0 1 1 1 0 0 0 1 1]; %polynome de rétroaction
Init_state_XBI = [0 1 0 1 0 1 1 1 0 0 1 0 0  %contient plusieurs états initiaux possibles (plusieurs PRN ???)
                  1 1 0 0 0 0 0 1 1 0 1 0 1
                  0 1 0 0 0 0 0 0 0 1 0 0 0
                  1 0 1 1 0 0 0 1 0 0 1 1 0
                  1 1 1 0 1 1 1 0 1 0 1 1 1
                  0 1 1 0 0 1 1 1 1 1 0 1 0
                  1 0 1 0 0 1 0 0 1 1 1 1 1
                  1 0 1 1 1 1 0 1 0 0 1 0 0
                  1 1 1 1 1 0 0 1 0 1 0 1 1
                  0 1 1 1 1 1 1 0 1 1 1 1 0
                  0 0 0 0 1 0 0 1 1 1 0 1 0
                  1 1 1 0 0 1 1 1 1 1 0 0 1
                  0 0 0 1 1 1 0 0 1 1 1 0 0
                  0 1 0 0 0 0 0 1 0 0 1 1 1
                  0 1 1 0 1 0 1 0 1 1 0 1 0
                  0 0 0 1 1 1 1 0 0 1 0 0 1
                  0 1 0 0 1 1 0 0 0 1 1 1 1 
                  1 1 1 1 0 0 0 0 1 1 1 1 0 
                  1 1 0 0 1 0 0 0 1 1 1 1 1 
                  0 1 1 0 1 0 1 1 0 1 1 0 1
                  0 0 1 0 0 0 0 0 0 1 0 0 0 
                  1 1 1 0 1 1 1 1 0 1 1 1 1
                  1 0 0 0 0 1 1 1 1 1 1 1 0 
                  1 1 0 0 0 1 0 1 1 0 1 0 0 
                  1 1 0 1 0 0 1 1 0 1 1 0 1 
                  1 0 1 0 1 1 0 0 1 0 1 1 0 
                  0 1 0 1 0 1 1 0 1 1 1 1 0 
                  0 1 1 1 1 0 1 0 1 0 1 1 0
                  0 1 0 1 1 1 1 1 0 0 0 0 1
                  1 0 0 0 0 1 0 1 1 0 1 1 1 
                  0 0 0 1 0 1 0 0 1 1 1 1 0
                  0 0 0 0 0 1 0 1 1 1 0 0 1 
                  1 1 0 1 0 1 0 0 0 0 0 0 1 
                  1 1 0 1 1 1 1 1 1 1 0 0 1
                  1 1 1 1 0 1 1 0 1 1 1 0 0 
                  1 0 0 1 0 1 1 0 0 1 0 0 0 
                  0 0 1 1 0 1 0 0 1 0 0 0 0];
              
State_XBI = Init_state_XBI(prn,:);
for i = 1:10230
    root = Poly_XBI.*State_XBI;
    reg = xor(xor(xor(root(1),root(3)),root(4)),root(6));
    reg = xor(xor(xor(xor(reg,root(7)),root(8)),root(12)),root(13));
    code_XBI(i) = State_XBI(13);
    for (j=13:-1:2)
        State_XBI(j) = State_XBI(j-1);
    end
    State_XBI(1) = reg;
end

%XBQ Coder
code_XBQ = zeros(1,10230);
Poly_XBQ = [1 0 1 1 0 1 1 1 0 0 0 1 1]; %polynome de rétrocation
Init_state_XBQ = [1 0 0 1 0 1 1 0 0 1 1 0 0 
                  0 1 0 0 0 1 1 1 1 0 1 1 0 %différents états initiaux pour plusieurs PRN
                  1 1 1 1 0 0 0 1 0 0 0 1 1 
                  0 0 1 1 1 0 1 1 0 1 0 1 0 
                  0 0 1 1 1 1 0 1 1 0 0 1 0 
                  0 1 0 1 0 1 0 1 0 1 0 0 1
                  1 1 1 1 1 1 0 0 0 0 0 0 1 
                  0 1 1 0 1 0 1 1 0 1 0 0 0 
                  1 0 1 1 1 0 1 0 0 0 0 1 1 
                  0 0 1 0 0 1 0 0 0 0 1 1 0 
                  0 0 0 1 0 0 0 0 0 0 1 0 1 
                  0 1 0 1 0 1 1 0 0 0 1 0 1 
                  0 1 0 0 1 1 0 1 0 0 1 0 1
                  1 0 1 0 0 0 0 1 1 1 1 1 1 
                  1 0 1 1 1 1 0 0 0 1 1 1 1 
                  1 1 0 1 0 0 1 0 1 1 1 1 1 
                  1 1 1 0 0 1 1 0 0 1 0 0 0 
                  1 0 1 1 0 1 1 1 0 0 1 0 0 
                  0 0 1 1 0 0 1 0 1 1 0 1 1 
                  1 1 0 0 0 0 1 1 1 0 0 0 1
                  0 1 1 0 1 1 0 0 1 0 0 0 0 
                  0 0 1 0 1 1 0 0 0 1 1 1 0
                  1 0 0 0 1 0 1 1 1 1 1 0 1 
                  0 1 1 0 1 1 1 1 1 0 0 1 1 
                  0 1 0 0 0 1 0 0 1 1 0 1 1 
                  0 1 0 1 0 1 0 1 1 1 1 0 0 
                  1 0 0 0 0 1 1 1 1 1 0 1 0 
                  1 1 1 1 1 0 1 0 0 0 0 1 0 
                  0 1 0 1 0 0 0 1 0 0 1 0 0
                  1 0 0 0 0 0 1 1 1 1 0 0 1
                  0 1 0 1 1 1 1 1 0 0 1 0 1
                  1 0 0 1 0 0 0 1 0 1 0 1 0
                  1 0 1 1 0 0 1 0 0 0 1 0 0
                  1 1 1 1 0 0 1 0 0 0 1 0 0
                  0 1 1 0 0 1 0 1 1 0 0 1 1 
                  0 0 1 1 1 1 0 1 0 1 1 1 1
                  0 0 1 0 0 1 1 0 1 0 0 0 1];
                  
%                   0 1 0 1 0 0 0 1 0 0 1 0 0 
%                   1 0 0 0 0 0 1 1 1 1 0 0 1
%                   0 1 0 1 1 1 1 1 0 0 1 0 1 
%                   0 1 0 1 1 1 1 1 0 0 0 0 1
%                   1 0 1 1 0 0 1 0 0 0 1 0 0 
%                   1 0 1 1 0 0 0 1 0 0 1 1 0 
%                   0 1 1 0 0 1 0 1 1 0 0 1 1 
%                   0 0 1 1 1 1 0 1 0 1 1 1 1 
%                   0 1 0 0 1 1 0 0 0 1 1 1 1];
              
State_XBQ = Init_state_XBQ(prn,:); %état initial séléctionner en fonction du PRN, chaque prn correspond à une ligne de la matrice
for i = 1:10230
    root = Poly_XBQ.*State_XBQ; %état actuel * polynome de rétroaction
    reg = xor(xor(xor(root(1),root(3)),root(4)),root(6));
    reg = xor(xor(xor(xor(reg,root(7)),root(8)),root(12)),root(13)); %comprend pas on recalcule reg ??
    code_XBQ(i) = State_XBQ(13); %dernier bit de state stocké dans code
    for (j=13:-1:2)
        State_XBQ(j) = State_XBQ(j-1); %décalage vers la droite
    end
    State_XBQ(1) = reg;%bit calculé (reg) inséré à gauche
end

%GENERATION OF L5
code_L5I = zeros(1,10230); %voie I du signal
code_L5Q = zeros(1,10230); %voie Q du signal
for i = 1:10230
    code_L5I(i) = (xor(code_XA(i),code_XBI(i))-0.5)*2; %comparaison entre les 2 bits obtenus par XA et XBI (sibits identique = 1, sinon =0)
    code_L5Q(i) = (xor(code_XA(i),code_XBQ(i))-0.5)*2; %et donne 1 ou -1 (modulation BPSK)
end
